"""Small, conservative PP2 decision procedure over supplied evidence."""

from __future__ import annotations

from typing import Mapping

from .guidelines import ACMG, find_rule
from .models import Assessment, BiologicalFeatures, Decision, EvidenceState, Status, Strength

_ASSESSMENTS = {
    "relevant_transcript", "consequence_reliable", "alternative_splice_consequence",
    "mapping_reliable", "gene_disease_valid", "missense_mechanism_common",
    "benign_missense_uncommon", "mechanism_compatible",
}


def _assessment(value: object) -> Assessment:
    if isinstance(value, Assessment):
        return value
    if value is None:
        return Assessment()
    if isinstance(value, bool):
        return Assessment(EvidenceState.KNOWN, value)
    if isinstance(value, Mapping):
        data = dict(value)
        data["state"] = EvidenceState(data.get("state", "KNOWN" if "value" in data else "UNKNOWN"))
        if "overlaps_criteria" in data:
            data["overlaps_criteria"] = tuple(data["overlaps_criteria"])
        return Assessment(**data)
    raise TypeError("assessment must be a boolean, mapping, Assessment, or None")


def _features(value: BiologicalFeatures | Mapping[str, object]) -> BiologicalFeatures:
    if isinstance(value, BiologicalFeatures):
        return value
    if not isinstance(value, Mapping):
        raise TypeError("biological_features must be BiologicalFeatures or a mapping")
    data = dict(value)
    for name in _ASSESSMENTS & data.keys():
        data[name] = _assessment(data[name])
    if data.get("requested_strength") is not None:
        data["requested_strength"] = Strength(data["requested_strength"])
    return BiologicalFeatures(**data)


def _curated(assessment: Assessment) -> bool:
    return (
        assessment.state == EvidenceState.KNOWN
        and assessment.source is not None
        and assessment.quality in {"HIGH", "MODERATE"}
        and assessment.independent_of_pp2 is not False
        and assessment.disease_match is not False
        and assessment.ascertainment_bias is not True
    )


def evaluate_pp2(biological_features: BiologicalFeatures | Mapping[str, object]) -> Decision:
    """Apply the pinned VCEP rule where scoped; otherwise apply qualitative ACMG PP2.

    Primary ACMG assessments must be curated. Counts, fractions and constraint are
    context only: no universal numeric PP2 thresholds were published in ACMG 2015.
    """
    f = _features(biological_features)
    rule = find_rule(f.gene_id, f.disease_id)
    guideline = rule.guideline if rule else ACMG
    path = "VCEP_DEFINED" if rule else "ESTABLISHED_MISSENSE_MECHANISM"
    used: list[str] = []
    missing: list[str] = []
    exclusions: list[str] = []

    def result(status: Status, reason: str) -> Decision:
        return Decision("PP2", status, Strength.SUPPORTING if status == Status.ACTIVATED else Strength.NONE,
                        path, tuple(used), tuple(missing), tuple(exclusions), guideline, reason)

    if rule and f.pp2_specification_id not in (None, guideline.rule_id):
        missing.append("applicable_pp2_specification")
        return result(Status.UNDETERMINED, "Supplied specification conflicts with the pinned gene/disease rule.")
    if not rule and f.pp2_specification_id:
        missing.append("verified_pp2_specification")
        return result(Status.UNDETERMINED, "The requested PP2 specification is not implemented or verified.")
    if rule and not rule.permitted:
        exclusions.append("vcep_pp2_not_applicable")
        used.extend(("gene_id", "disease_id"))
        return result(Status.NOT_APPLICABLE, "The scoped ClinGen VCEP specification makes PP2 inapplicable.")

    if f.variant_consequence is None:
        missing.append("variant_consequence")
    elif f.variant_consequence.upper() not in {"MISSENSE", "MISSENSE_VARIANT"}:
        exclusions.append("nonmissense_consequence")
        used.append("variant_consequence")
        return result(Status.NOT_MET, "The candidate is not a missense variant.")
    else:
        used.append("variant_consequence")

    for name in ("consequence_reliable", "mapping_reliable", "relevant_transcript", "gene_disease_valid", "mechanism_compatible"):
        item: Assessment = getattr(f, name)
        if item.state == EvidenceState.KNOWN and item.value is False:
            exclusions.append(name)
            used.append(name)
            return result(Status.NOT_MET, f"{name} is contradicted by supplied evidence.")
        if item.state == EvidenceState.CONFLICTING:
            missing.append(name)
        elif item.state == EvidenceState.KNOWN:
            used.append(name)
        elif name in {"relevant_transcript", "gene_disease_valid"}:
            missing.append(name)

    if f.alternative_splice_consequence.state == EvidenceState.KNOWN and f.alternative_splice_consequence.value:
        missing.append("missense_primary_effect_vs_splicing")
        used.append("alternative_splice_consequence")
    elif f.alternative_splice_consequence.state == EvidenceState.CONFLICTING:
        missing.append("alternative_splice_consequence")

    if rule:
        for name in ("missense_mechanism_common", "benign_missense_uncommon"):
            item: Assessment = getattr(f, name)
            if item.state == EvidenceState.CONFLICTING:
                missing.append(name)
            elif _curated(item) and item.value is False:
                exclusions.append(name)
                used.append(name)
                return result(Status.NOT_MET, f"New curated evidence contradicts the VCEP PP2 premise: {name}.")

        if f.inheritance_mode is None:
            missing.append("inheritance_mode")
        elif f.inheritance_mode.upper() != rule.inheritance_mode:
            exclusions.append("vcep_inheritance_mismatch")
            return result(Status.NOT_APPLICABLE, "This VCEP PP2 rule applies only to autosomal dominant disease.")
        else:
            used.append("inheritance_mode")
        if f.transcript_id_version is None:
            missing.append("transcript_id_version")
        elif f.transcript_id_version != rule.transcript_id_version:
            exclusions.append("vcep_transcript_mismatch")
            return result(Status.NOT_MET, "The supplied transcript does not match the pinned VCEP transcript.")
        else:
            used.append("transcript_id_version")
        if missing:
            return result(Status.UNDETERMINED, "Essential candidate, gene/disease, or VCEP scope evidence is missing or conflicting.")
        used.extend(("gene_id", "disease_id", "vcep_pp2_specification"))
        return result(Status.ACTIVATED, "The candidate meets the ACTA1 autosomal dominant PP2 specification at Supporting strength.")

    if not f.gene_id:
        missing.append("gene_id")
    if not f.disease_id:
        missing.append("disease_id")
    for name in ("missense_mechanism_common", "benign_missense_uncommon"):
        item: Assessment = getattr(f, name)
        if item.state == EvidenceState.CONFLICTING:
            missing.append(name)
        elif not _curated(item):
            missing.append(name)
        elif item.value is False:
            used.append(name)
            exclusions.append(name)
            return result(Status.NOT_MET, f"Curated evidence establishes that {name} is false for this disease.")
        else:
            used.append(name)
    if missing:
        return result(Status.UNDETERMINED, "Essential PP2 evidence is missing, unreliable, or conflicting.")
    used.extend(("gene_id", "disease_id"))
    return result(Status.ACTIVATED, "Curated evidence establishes a common disease-matched missense mechanism and uncommon benign missense variation.")
