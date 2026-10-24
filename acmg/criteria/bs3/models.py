"""Small, typed input and auditable output models for BS3."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any, Mapping


class EvidenceState(StrEnum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"
    NOT_APPLICABLE = "NOT_APPLICABLE"


@dataclass(frozen=True)
class Evidence:
    state: EvidenceState = EvidenceState.UNKNOWN
    value: Any = None
    source: str | None = None

    @classmethod
    def parse(cls, raw: Any) -> "Evidence":
        if isinstance(raw, cls):
            return raw
        if raw is None:
            return cls()
        if isinstance(raw, Mapping) and "state" in raw:
            state = EvidenceState(str(raw["state"]).upper())
            value = raw.get("value")
            if state == EvidenceState.KNOWN and value is None:
                raise ValueError("KNOWN evidence requires a value")
            return cls(state, value if state == EvidenceState.KNOWN else None, raw.get("source"))
        return cls(EvidenceState.KNOWN, raw)


@dataclass(frozen=True)
class Context:
    gene_id: Evidence = field(default_factory=Evidence)
    disease_id: Evidence = field(default_factory=Evidence)
    disease_mechanism: Evidence = field(default_factory=Evidence)
    variant_id: Evidence = field(default_factory=Evidence)
    variant_consequence: Evidence = field(default_factory=Evidence)
    gene_disease_validity: Evidence = field(default_factory=Evidence)
    splice_effect_possible: Evidence = field(default_factory=Evidence)
    variant_mechanism_compatibility: Evidence = field(default_factory=Evidence)


@dataclass(frozen=True)
class Assay:
    assay_id: Evidence = field(default_factory=Evidence)
    assay_type: Evidence = field(default_factory=Evidence)
    study_id: Evidence = field(default_factory=Evidence)
    detectable_mechanisms: Evidence = field(default_factory=Evidence)
    mechanism_compatible: Evidence = field(default_factory=Evidence)
    variant_class_applicable: Evidence = field(default_factory=Evidence)
    detects_relevant_defect: Evidence = field(default_factory=Evidence)
    assay_sensitivity_adequate: Evidence = field(default_factory=Evidence)
    validation_status: Evidence = field(default_factory=Evidence)
    result_class: Evidence = field(default_factory=Evidence)
    result_value: Evidence = field(default_factory=Evidence)
    validated_normal_range: Evidence = field(default_factory=Evidence)
    result_reliable: Evidence = field(default_factory=Evidence)
    replicate_consistent: Evidence = field(default_factory=Evidence)
    batch_effect: Evidence = field(default_factory=Evidence)
    artifact: Evidence = field(default_factory=Evidence)
    saturation: Evidence = field(default_factory=Evidence)
    floor_effect: Evidence = field(default_factory=Evidence)
    wildtype_control: Evidence = field(default_factory=Evidence)
    pathogenic_control_count: Evidence = field(default_factory=Evidence)
    benign_control_count: Evidence = field(default_factory=Evidence)
    controls_independently_classified: Evidence = field(default_factory=Evidence)
    control_separation: Evidence = field(default_factory=Evidence)
    control_classification_quality: Evidence = field(default_factory=Evidence)
    calibration_validated: Evidence = field(default_factory=Evidence)
    odds_path: Evidence = field(default_factory=Evidence)
    accepted_normal_interpretation: Evidence = field(default_factory=Evidence)
    historical_acceptance_documented: Evidence = field(default_factory=Evidence)
    tissue_relevant: Evidence = field(default_factory=Evidence)
    transcript_relevant: Evidence = field(default_factory=Evidence)
    rna_detection_adequate: Evidence = field(default_factory=Evidence)
    nmd_addressed: Evidence = field(default_factory=Evidence)
    variant_allele_expressed: Evidence = field(default_factory=Evidence)
    overlaps_ps3: Evidence = field(default_factory=Evidence)
    evidence_conflict: Evidence = field(default_factory=Evidence)
    evidence_source: Evidence = field(default_factory=Evidence)
    mechanism_coverage_complete: Evidence = field(default_factory=Evidence)
    expression_adequate: Evidence = field(default_factory=Evidence)
    result_uncertainty_acceptable: Evidence = field(default_factory=Evidence)
    ps3_claimed: Evidence = field(default_factory=Evidence)


@dataclass(frozen=True)
class BS3Input:
    context: Context
    assays: tuple[Assay, ...]


@dataclass(frozen=True)
class BS3Result:
    criterion: str
    status: str
    strength: str
    decision_path: str
    assays_used: tuple[str, ...]
    features_used: tuple[str, ...]
    missing_critical_features: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    conflicting_evidence: tuple[str, ...]
    guideline: Mapping[str, str]
    reason: str

    def to_dict(self) -> dict[str, Any]:
        return {name: list(getattr(self, name)) if name in {
            "assays_used", "features_used", "missing_critical_features",
            "exclusions_evaluated", "conflicting_evidence",
        } else dict(getattr(self, name)) if name == "guideline" else getattr(self, name)
                for name in self.__dataclass_fields__}


def parse_input(raw: BS3Input | Mapping[str, Any]) -> BS3Input:
    if isinstance(raw, BS3Input):
        return raw
    if not isinstance(raw, Mapping):
        raise TypeError("BS3 input must be a BS3Input or mapping")
    context_raw = raw.get("context", raw)
    if not isinstance(context_raw, Mapping):
        raise TypeError("context must be a mapping")
    aliases = {"disease_molecular_mechanism": "disease_mechanism"}
    context = Context(**{f: Evidence.parse(context_raw.get(f, context_raw.get(next((a for a, b in aliases.items() if b == f), ""))))
                         for f in Context.__dataclass_fields__})
    assays_raw = raw.get("assays", ())
    if not isinstance(assays_raw, (list, tuple)):
        raise TypeError("assays must be a sequence")
    assay_aliases = {
        "assay_validation_status": "validation_status",
        "variant_result_classification": "result_class",
        "variant_assay_result": "result_value",
        "assay_odds_pathogenicity": "odds_path",
        "assay_mechanism_compatibility": "mechanism_compatible",
        "assay_variant_class_applicability": "variant_class_applicable",
        "variant_result_replicate_consistency": "replicate_consistent",
        "experimental_batch_effect": "batch_effect",
        "experimental_artifact_evidence": "artifact",
        "assay_saturation": "saturation",
        "assay_floor_effect": "floor_effect",
        "benign_control_count": "benign_control_count",
        "control_classification_independence": "controls_independently_classified",
        "control_result_separation": "control_separation",
        "assay_control_reliability": "control_classification_quality",
        "mechanism_specific_assay_sensitivity": "assay_sensitivity_adequate",
        "conflicting_functional_results": "evidence_conflict",
        "rna_tissue_relevance": "tissue_relevant",
        "rna_transcript_relevance": "transcript_relevant",
        "assay_measures_disease_mechanism": "mechanism_compatible",
        "assay_measures_variant_effect": "detects_relevant_defect",
        "assay_sensitivity_sufficient": "assay_sensitivity_adequate",
        "variant_expression_adequate": "expression_adequate",
        "variant_result_uncertainty_acceptable": "result_uncertainty_acceptable",
        "control_classification_quality_adequate": "control_classification_quality",
    }
    assays = []
    for entry in assays_raw:
        if isinstance(entry, Assay):
            assays.append(entry)
        elif isinstance(entry, Mapping):
            normalized = {assay_aliases.get(k, k): v for k, v in entry.items()}
            assays.append(Assay(**{f: Evidence.parse(normalized.get(f)) for f in Assay.__dataclass_fields__}))
        else:
            raise TypeError("each assay must be an Assay or mapping")
    return BS3Input(context, tuple(assays))
