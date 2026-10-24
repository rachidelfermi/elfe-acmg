"""Deterministic, offline PM4 decisions from supplied molecular evidence."""

from dataclasses import fields
from typing import Any, Mapping

from .guidelines import GENERAL, RETT, RETT_GENES, RETT_REPEAT_PROTEINS, RETT_REPEAT_REGIONS
from .models import BiologicalFeatures, Evidence, EvidenceState, PM4Result


def _value(raw: Any) -> Any:
    if isinstance(raw, Evidence):
        return raw.value if raw.state == EvidenceState.KNOWN else None
    if isinstance(raw, Mapping) and "state" in raw:
        state = EvidenceState(raw["state"])
        if (state == EvidenceState.KNOWN) != (raw.get("value") is not None):
            raise ValueError("KNOWN evidence needs a value; other states cannot have one")
        return raw.get("value") if state == EvidenceState.KNOWN else None
    return raw


def _state(raw: Any) -> EvidenceState:
    if isinstance(raw, Evidence):
        return raw.state
    if isinstance(raw, Mapping) and "state" in raw:
        return EvidenceState(raw["state"])
    return EvidenceState.KNOWN if raw is not None else EvidenceState.UNKNOWN


def _overlap(start: int, end: int, intervals: tuple[tuple[int, int], ...]) -> bool:
    return any(start <= right and end >= left for left, right in intervals)


def evaluate_pm4(biological_features: BiologicalFeatures | Mapping[str, Any]) -> PM4Result:
    """Evaluate PM4 only; no annotation, database lookup, or overall classification.

    Boolean evidence may be bool/None or {state, value, source}. Conflicting
    essential evidence yields UNDETERMINED. Unknown extra input keys are errors.
    """
    if isinstance(biological_features, Mapping):
        allowed = {f.name for f in fields(BiologicalFeatures)}
        extra = set(biological_features) - allowed
        if extra:
            raise ValueError(f"Unsupported PM4 input fields: {', '.join(sorted(extra))}")
        feature = BiologicalFeatures(**biological_features)
    elif isinstance(biological_features, BiologicalFeatures):
        feature = biological_features
    else:
        raise TypeError("biological_features must be BiologicalFeatures or a mapping")

    used: list[str] = []
    missing: list[str] = []
    exclusions: list[str] = []

    def get(name: str) -> Any:
        value = _value(getattr(feature, name))
        if value is not None and name not in used:
            used.append(name)
        return value

    gene = get("gene_id")
    profile = get("vcep_specification_version")
    guideline = RETT if profile == "rett_angelman_2022" else GENERAL
    if profile is not None and (profile != "rett_angelman_2022" or gene not in RETT_GENES):
        return PM4Result("PM4", "UNDETERMINED", "NONE", None, tuple(used),
                         ("applicable_vcep_specification",), (), GENERAL,
                         "The supplied VCEP profile is unverified or outside its gene scope.")

    def result(status: str, reason: str, path: str | None = None,
               strength: str = "NONE") -> PM4Result:
        return PM4Result("PM4", status, strength, path, tuple(used),
                         tuple(dict.fromkeys(missing)), tuple(dict.fromkeys(exclusions)),
                         guideline, reason)

    consequence = get("variant_consequence")
    if _state(feature.variant_consequence) == EvidenceState.CONFLICTING:
        missing.append("variant_consequence")
        return result("UNDETERMINED", "Conflicting molecular consequence annotations.")
    if consequence is None:
        missing.append("variant_consequence")
        return result("UNDETERMINED", "Molecular consequence is unknown.")
    valid = {"INFRAME_DELETION", "INFRAME_INSERTION", "STOP_LOSS"}
    if consequence not in valid:
        return result("NOT_APPLICABLE", "This consequence is outside PM4 scope.")
    path = consequence

    if gene is None:
        missing.append("gene_id")
    if get("disease_id") is None:
        missing.append("disease_id")
    transcript = get("transcript_id_version")
    relevance = get("transcript_biological_relevance")
    disease_transcript = get("disease_relevant_transcript")
    if relevance is False or (disease_transcript is not None and transcript is not None
                            and disease_transcript != transcript):
        return result("NOT_MET", "The evaluated protein consequence is not on the disease-relevant transcript.", path)
    if transcript is None and relevance is not True:
        missing.append("transcript_id_version or transcript_biological_relevance")
    elif relevance is None and disease_transcript is None:
        missing.append("transcript_biological_relevance or disease_relevant_transcript")
    if feature.relevant_transcript_consequences:
        used.append("relevant_transcript_consequences")
        if any(_value(item) != consequence for item in feature.relevant_transcript_consequences):
            return result("UNDETERMINED", "Relevant transcripts have conflicting protein consequences.", path)
    if feature.evidence_conflicts:
        used.append("evidence_conflicts")
        return result("UNDETERMINED", "Supplied evidence contains unresolved conflicts.", path)

    for name in ("variant_call_quality", "indel_alignment_reliability",
                 "protein_consequence_prediction_quality"):
        value = get(name)
        if value is False or _state(getattr(feature, name)) == EvidenceState.CONFLICTING:
            return result("UNDETERMINED", f"{name} does not support a reliable protein consequence.", path)
    if get("gene_disease_validity") is False or get("protein_length_change_mechanism_compatibility") is False:
        return result("NOT_MET", "The supplied gene-disease or molecular mechanism evidence contradicts PM4 relevance.", path)
    if get("alternative_transcript_rescue") is True:
        exclusions.append("alternative_transcript_rescue")
        return result("NOT_MET", "A disease-relevant alternative transcript rescues the length change.", path)
    if get("pvs1_applied_to_same_consequence") is True:
        exclusions.append("pvs1_applied_to_same_consequence")
        return result("NOT_MET", "The same molecular consequence is already counted under PVS1.", path)
    rna = get("rna_splicing_evidence")
    if rna == "CONTRADICTS_LENGTH_CHANGE":
        return result("NOT_MET", "RNA evidence contradicts the predicted protein length change.", path)
    if rna == "CONFLICTING":
        return result("UNDETERMINED", "RNA and predicted molecular consequences conflict.", path)

    ref = get("reference_protein_length")
    alt = get("alternate_protein_length")
    difference = get("protein_length_difference")
    if ref is not None and alt is not None:
        if ref <= 0 or alt <= 0:
            raise ValueError("Protein lengths must be positive")
        derived = alt - ref
        if difference is not None and difference != derived:
            return result("UNDETERMINED", "Supplied protein lengths and difference disagree.", path)
        difference = derived

    if consequence in {"INFRAME_DELETION", "INFRAME_INSERTION"}:
        frame = get("reading_frame_preserved")
        if frame is False:
            return result("NOT_APPLICABLE", "The reported indel does not preserve the reading frame.", path)
        if frame is None:
            missing.append("reading_frame_preserved")
        count_name = "deleted_amino_acid_count" if consequence == "INFRAME_DELETION" else "inserted_amino_acid_count"
        count = get(count_name)
        if count is not None and count < 0:
            raise ValueError(f"{count_name} cannot be negative")
        if difference is None and count is not None:
            difference = -count if consequence == "INFRAME_DELETION" else count
        if difference is None:
            missing.append("protein_length_difference or protein lengths or amino acid count")
        elif difference == 0:
            return result("NOT_MET", "The protein has no net length change.", path)
        elif (consequence == "INFRAME_DELETION" and difference > 0) or (consequence == "INFRAME_INSERTION" and difference < 0):
            return result("UNDETERMINED", "The reported indel type disagrees with protein length direction.", path)
        if count is not None and difference is not None and abs(difference) != count:
            return result("UNDETERMINED", "Amino acid count disagrees with protein length change.", path)

        repeat = get("repeat_region_status")
        start, end = get("protein_change_start"), get("protein_change_end")
        if (profile == "rett_angelman_2022" and gene in RETT_REPEAT_REGIONS
                and get("protein_accession_version") == RETT_REPEAT_PROTEINS[gene]
                and start is not None and end is not None):
            if _overlap(start, end, RETT_REPEAT_REGIONS[gene]):
                repeat = True
                exclusions.append("vcep_defined_repeat_region")
        if get("known_benign_repeat_length_variation") is True:
            exclusions.append("known_benign_repeat_length_variation")
            return result("NOT_MET", "The change is established benign repeat-length variation.", path)
        if repeat is True:
            exclusions.append("repeat_region_status")
            return result("NOT_MET", "PM4 excludes this in-frame change in a repeat region.", path)
        if repeat is None:
            missing.append("repeat_region_status")
        if get("known_benign_indel_same_region") is True or get("regional_indel_tolerance") is True:
            exclusions.append("regional_benign_length_variation")
            return result("UNDETERMINED", "Regional benign length variation needs case-specific review.", path)
    else:
        native = get("native_stop_codon_disrupted")
        if native is False:
            return result("NOT_MET", "The native termination codon is retained.", path)
        if native is None:
            missing.append("native_stop_codon_disrupted")
        if get("stop_loss_rescue_evidence") is True:
            exclusions.append("stop_loss_rescue_evidence")
            return result("NOT_MET", "The stop-loss effect is rescued.", path)
        extension = get("protein_extension_length")
        predicted = get("predicted_protein_extension")
        if extension is not None and extension < 0:
            raise ValueError("protein_extension_length cannot be negative")
        if predicted is False or extension == 0 or (difference is not None and difference <= 0):
            return result("NOT_MET", "No protein extension is supported.", path)
        if extension is not None and difference is not None and extension != difference:
            return result("UNDETERMINED", "Extension length disagrees with protein length difference.", path)
        if not (extension is not None and extension > 0 or difference is not None and difference > 0 or predicted is True):
            missing.append("predicted_protein_extension or protein_extension_length")
        if get("downstream_inframe_stop_codon") is False and get("protein_product_confirmation") is False:
            return result("UNDETERMINED", "A stable extended protein is not supported.", path)

    if missing:
        return result("UNDETERMINED", "Essential evidence is missing: " + ", ".join(dict.fromkeys(missing)) + ".", path)

    strength = "MODERATE"
    if profile == "rett_angelman_2022":
        if consequence == "STOP_LOSS" and gene in {"MECP2", "UBE3A"}:
            strength = "STRONG"
        elif consequence != "STOP_LOSS" and difference is not None and abs(difference) <= 2:
            # The panel retains Moderate for a small event in a documented
            # functionally important region (its PM1 exception).
            important = get("critical_protein_domain_overlap") is True and get("affected_region_functional_evidence") is True
            strength = "MODERATE" if important else "SUPPORTING"
    return result("ACTIVATED", "Verified PM4 protein length change on a relevant transcript with applicable exclusions addressed.", path, strength)
