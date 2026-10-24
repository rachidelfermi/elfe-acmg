"""Deterministic, conservative BP3 decision procedure."""

from .guidelines import (
    BMPR2_TRANSCRIPT, FOXG1_PROTEIN, FOXG1_REGIONS, FOXG1_TRANSCRIPT,
    select_guideline,
)
from .models import BP3Result, BiologicalFeatures, Evidence, EvidenceState


def _value(value):
    if isinstance(value, Evidence):
        return value.value if value.state == EvidenceState.KNOWN else None
    return value


def _conflicting(value) -> bool:
    return isinstance(value, Evidence) and value.state == EvidenceState.CONFLICTING


def evaluate_bp3(features: BiologicalFeatures) -> BP3Result:
    """Evaluate extracted protein evidence; never infer missing biology from HGVS text."""
    if not isinstance(features, BiologicalFeatures):
        raise TypeError("features must be BiologicalFeatures")
    profile, guideline = select_guideline(features.gene_id, features.disease_id)
    used: list[str] = []
    missing: list[str] = []
    exclusions: list[str] = []

    def seen(name: str):
        if name not in used:
            used.append(name)
        return _value(getattr(features, name))

    def result(status: str, reason: str, path: str | None = None) -> BP3Result:
        return BP3Result("BP3", status, "SUPPORTING" if status == "ACTIVATED" else "NONE",
                         path, tuple(used), tuple(missing), tuple(exclusions), guideline, reason)

    if profile == "SCOPE_UNKNOWN":
        missing.append("disease_id (verified VCEP scope)")
        return result("UNDETERMINED", "Disease scope for the gene-specific BP3 specification is unresolved.")
    if profile in {"RYR1", "MYOC"}:
        return result("NOT_APPLICABLE", "The applicable VCEP specifies BP3 as not applicable.", "VCEP_NOT_APPLICABLE")
    if not features.gene_id or not features.disease_id:
        missing.extend(name for name in ("gene_id", "disease_id") if not getattr(features, name))
        return result("UNDETERMINED", "Gene and disease context is needed to check applicable BP3 specifications.")
    if features.vcep_specification_version and profile != "GENERAL":
        seen("vcep_specification_version")
        if features.vcep_specification_version != guideline["version"]:
            return result("UNDETERMINED", "Requested VCEP version differs from the implemented, verified version.")

    # Explicit conflicts defeat an otherwise complete pathway.
    decisive = (
        "variant_consequence", "reading_frame_preserved", "biologically_relevant_transcript",
        "protein_consequence_confidence", "repeat_region_status", "affected_protein_positions",
        "repeat_region_functional_importance", "repeat_region_disease_relevance",
        "critical_protein_domain_overlap", "active_site_overlap", "binding_site_overlap",
        "protein_interaction_region_overlap", "repeat_region_functional_validation",
        "known_pathogenic_repeat_length_variants", "gene_disease_validity",
        "repeat_length_disease_mechanism", "gene_specific_repeat_pathogenicity",
        "disease_specific_repeat_exceptions", "indel_alignment_reliability",
        "variant_call_quality", "repeat_region_mapping_confidence", "repeat_mapping_ambiguity",
        "protein_consequence_prediction_quality", "alternative_splice_consequence",
        "rna_splicing_evidence", "alternative_transcript_consequence", "complex_allele_evidence",
        "pm4_applied_to_same_consequence",
    )
    for name in decisive:
        if _conflicting(getattr(features, name)):
            seen(name)
            return result("UNDETERMINED", f"Conflicting evidence: {name}.")
    if seen("evidence_conflict_status") is True:
        return result("UNDETERMINED", "Evidence sources conflict; BP3 cannot be resolved.")
    for name in ("evidence_quality", "repeat_region_clinical_evidence_quality"):
        if getattr(features, name) in {"UNRELIABLE", "CONFLICTING"}:
            seen(name)
            return result("UNDETERMINED", f"Evidence quality is insufficient: {name}.")

    consequence = seen("variant_consequence")
    allowed = {"INFRAME_DELETION", "INFRAME_INSERTION"}
    if consequence is not None and consequence not in allowed:
        return result("NOT_MET", "The protein consequence is not a qualifying in-frame insertion or deletion.", "NONQUALIFYING_CONSEQUENCE")
    if seen("reading_frame_preserved") is False:
        if consequence in allowed:
            return result("UNDETERMINED", "The in-frame consequence label conflicts with frame evidence.")
        return result("NOT_MET", "The reading frame is not preserved.", "NONQUALIFYING_CONSEQUENCE")
    length_change = seen("protein_length_change")
    if length_change == 0:
        return result("NOT_MET", "There is no protein length change.", "NONQUALIFYING_CONSEQUENCE")
    if consequence == "INFRAME_DELETION" and length_change is not None and length_change > 0:
        return result("UNDETERMINED", "Deletion label conflicts with the positive protein length change.")
    if consequence == "INFRAME_INSERTION" and length_change is not None and length_change < 0:
        return result("UNDETERMINED", "Insertion label conflicts with the negative protein length change.")
    if seen("biologically_relevant_transcript") is False:
        return result("NOT_MET", "The in-frame consequence is not on a biologically relevant transcript.", "TRANSCRIPT_MISMATCH")
    for name in ("protein_consequence_confidence", "indel_alignment_reliability",
                 "variant_call_quality", "protein_consequence_prediction_quality",
                 "repeat_region_mapping_confidence"):
        if seen(name) is False:
            exclusions.append(name)
            return result("UNDETERMINED", f"Unreliable evidence: {name}.")
    for name in ("repeat_mapping_ambiguity", "alternative_splice_consequence",
                 "rna_splicing_evidence", "alternative_transcript_consequence",
                 "complex_allele_evidence", "pm4_applied_to_same_consequence"):
        if seen(name) is True:
            exclusions.append(name)
            return result("UNDETERMINED", f"{name} conflicts with an unambiguous BP3 protein consequence.")

    if seen("gene_disease_validity") is False:
        return result("NOT_APPLICABLE", "The asserted gene–disease relationship is invalid.")
    for name in ("repeat_length_disease_mechanism", "gene_specific_repeat_pathogenicity",
                 "disease_specific_repeat_exceptions", "known_pathogenic_repeat_length_variants",
                 "repeat_region_disease_relevance"):
        if seen(name) is True:
            exclusions.append(name)
            return result("NOT_MET", f"Known disease-relevant repeat-length evidence excludes BP3: {name}.", "PATHOGENIC_REPEAT_MECHANISM")
    if seen("protein_length_change_mechanism_compatibility") is False:
        return result("NOT_MET", "The protein length change is incompatible with benign repeat tolerance.")
    for name in ("critical_protein_domain_overlap", "active_site_overlap",
                 "binding_site_overlap", "protein_interaction_region_overlap"):
        if seen(name) is True:
            exclusions.append(name)
            return result("NOT_MET", f"The repeat overlaps an established functional region: {name}.", "FUNCTIONAL_REPEAT")
    functional = seen("repeat_region_functional_importance")
    if functional in {"CRITICAL", "ESTABLISHED_FUNCTION"}:
        return result("NOT_MET", "The repeat has established functional importance.", "FUNCTIONAL_REPEAT")
    if seen("repeat_region_functional_validation") == "DAMAGING":
        return result("NOT_MET", "Validated damaging repeat-region function conflicts with BP3.", "FUNCTIONAL_REPEAT")

    if consequence is None:
        missing.append("variant_consequence")
    # A confirmed INFRAME_DELETION/INSERTION consequence establishes frame preservation.
    if seen("biologically_relevant_transcript") is None:
        missing.append("biologically_relevant_transcript")
    if seen("protein_consequence_confidence") is None:
        missing.append("protein_consequence_confidence")
    if profile == "FOXG1":
        if not features.transcript_id_version:
            missing.append("transcript_id_version")
        elif features.transcript_id_version != FOXG1_TRANSCRIPT:
            return result("UNDETERMINED", "FOXG1 transcript does not match the VCEP reference.")
        else:
            seen("transcript_id_version")
        if not features.protein_accession_version:
            missing.append("protein_accession_version")
        elif features.protein_accession_version != FOXG1_PROTEIN:
            return result("UNDETERMINED", "FOXG1 protein accession does not match the VCEP coordinate reference.")
        else:
            seen("protein_accession_version")
    if profile == "BMPR2":
        if not features.transcript_id_version:
            missing.append("transcript_id_version")
        elif features.transcript_id_version != BMPR2_TRANSCRIPT:
            return result("UNDETERMINED", "BMPR2 transcript differs from the VCEP reference.")
        else:
            seen("transcript_id_version")

    repeat = seen("repeat_region_status")
    positions = seen("affected_protein_positions")
    boundaries = features.repeat_region_boundaries
    if positions is not None and (len(positions) != 2 or positions[0] < 1 or positions[0] > positions[1]):
        return result("UNDETERMINED", "Affected protein coordinates are invalid.")
    if repeat is True and positions and boundaries and not (boundaries[0] <= positions[0] <= positions[1] <= boundaries[1]):
        return result("UNDETERMINED", "Repeat annotation and affected protein coordinates conflict.")
    if profile == "FOXG1":
        if positions is None:
            missing.append("affected_protein_positions")
        elif not any(start <= positions[0] <= positions[1] <= end for start, end in FOXG1_REGIONS):
            return result("NOT_MET", "The indel is outside the FOXG1 BP3 repeat regions.", "VCEP_REPEAT_DEFINITION")
        if repeat is False:
            return result("UNDETERMINED", "Repeat annotation conflicts with the FOXG1 VCEP coordinates.")
        # The VCEP explicitly identifies these polymorphic repeats as eligible for BP3.
        if functional not in (None, "UNKNOWN", "NO_KNOWN_FUNCTION"):
            return result("UNDETERMINED", "Functional assessment conflicts with the VCEP-defined repeat.")
    else:
        if repeat is False:
            return result("NOT_MET", "The in-frame indel is outside a repetitive region.", "NONREPEAT_INDEL")
        if repeat is None:
            missing.append("repeat_region_status")
        if functional != "NO_KNOWN_FUNCTION":
            missing.append("repeat_region_functional_importance (affirmative no-known-function assessment)")
        # General ACMG wording does not create an independent tolerance-only pathway.
        if not any(seen(name) is False for name in (
            "repeat_length_disease_mechanism", "disease_specific_repeat_exceptions",
            "gene_specific_repeat_pathogenicity")):
            missing.append("disease-specific pathogenic repeat mechanism assessment")
    if missing:
        return result("UNDETERMINED", "Essential BP3 evidence is missing or unresolved.")
    if profile == "FOXG1":
        return result("ACTIVATED", "Validated in-frame indel lies wholly in a FOXG1 VCEP-defined polymorphic repeat.", "VCEP_FOXG1_POLYMORPHIC_REPEAT")
    return result("ACTIVATED", "Reliable in-frame indel lies in a verified repeat assessed to have no known function, with disease-specific repeat mechanisms excluded.", "STANDARD_NONFUNCTIONAL_REPEAT")
