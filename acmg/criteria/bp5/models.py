"""Typed, already-curated case evidence. Missing values are UNKNOWN, never false."""

from dataclasses import asdict, dataclass, fields
from enum import StrEnum
from typing import Any, Mapping


class EvidenceState(StrEnum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"
    NOT_APPLICABLE = "NOT_APPLICABLE"


@dataclass(frozen=True)
class Evidence:
    value: Any = None
    state: EvidenceState = EvidenceState.UNKNOWN
    source: str | None = None
    version: str | None = None
    quality: str | None = None
    missing_reason: str | None = None

    @classmethod
    def parse(cls, raw: Any) -> "Evidence":
        if isinstance(raw, cls):
            result = raw
        elif isinstance(raw, Mapping) and "state" in raw:
            result = cls(raw.get("value"), EvidenceState(raw["state"]), raw.get("source"),
                         raw.get("version"), raw.get("quality"), raw.get("missing_reason"))
        elif raw is None:
            result = cls()
        else:
            result = cls(raw, EvidenceState.KNOWN)
        if (result.state == EvidenceState.KNOWN) != (result.value is not None):
            raise ValueError("KNOWN evidence needs a value; other states cannot have one")
        return result


@dataclass(frozen=True)
class BP5Features:
    # Candidate locus and its possible independent contribution.
    candidate_variant_id: Any = None
    candidate_gene_id: Any = None
    candidate_disease_id: Any = None
    candidate_variant_consequence: Any = None
    candidate_genotype: Any = None
    candidate_zygosity: Any = None
    candidate_gene_disease_validity: Any = None
    candidate_disease_inheritance: Any = None
    candidate_disease_mechanism: Any = None
    candidate_phenotype_consistency: Any = None
    candidate_independent_disease_evidence: Any = None
    candidate_specific_phenotype_evidence: Any = None
    candidate_additional_phenotype_contribution: Any = None
    candidate_disease_mechanism_independence: Any = None
    candidate_segregation_evidence: Any = None
    candidate_alternative_disease_overlap: Any = None
    candidate_alternative_causal_relationship: Any = None
    # Alternative molecular diagnosis. A classification alone is insufficient.
    alternative_diagnosis_id: Any = None
    alternative_gene_id: Any = None
    alternative_variant_id: Any = None
    alternative_variant_classification: Any = None
    alternative_variant_classification_source: Any = None
    alternative_variant_classification_evidence: Any = None
    alternative_variant_classification_independence: Any = None
    alternative_genotype: Any = None
    alternative_zygosity: Any = None
    alternative_inheritance_compatibility: Any = None
    alternative_gene_disease_validity: Any = None
    alternative_disease_mechanism: Any = None
    alternative_diagnosis_confidence: Any = None
    alternative_diagnosis_genotype_confirmation: Any = None
    alternative_diagnosis_clinical_confirmation: Any = None
    alternative_diagnosis_biochemical_confirmation: Any = None
    alternative_diagnosis_functional_confirmation: Any = None
    alternative_diagnosis_segregation_support: Any = None
    alternative_diagnosis_diagnostic_criteria_met: Any = None
    alternative_diagnosis_confirmation_quality: Any = None
    # Clinical interpretation is supplied by the caller, not calculated from HPO scores.
    patient_hpo_terms: Any = None
    patient_clinical_findings: Any = None
    patient_phenotype_present: Any = None
    patient_age_of_onset: Any = None
    patient_disease_progression: Any = None
    patient_phenotype_specificity: Any = None
    patient_phenotype_completeness: Any = None
    alternative_diagnosis_phenotype_match: Any = None
    alternative_diagnosis_explains_core_phenotypes: Any = None
    alternative_diagnosis_explains_all_relevant_phenotypes: Any = None
    candidate_disease_phenotype_match: Any = None
    unexplained_phenotype_features: Any = None
    phenotype_evidence_quality: Any = None
    # Known complications and family context.
    dual_diagnosis_evidence: Any = None
    blended_phenotype_evidence: Any = None
    multiple_molecular_diagnoses: Any = None
    alternative_diagnosis_exclusivity: Any = None
    family_history: Any = None
    family_phenotype_consistency: Any = None
    alternative_diagnosis_segregation_evidence: Any = None
    family_inheritance_pattern: Any = None
    family_evidence_quality: Any = None
    phenocopy_evidence: Any = None
    environmental_disease_explanation: Any = None
    non_genetic_diagnosis: Any = None
    somatic_variant_context: Any = None
    mosaicism_context: Any = None
    oligogenic_disease_evidence: Any = None
    modifier_variant_evidence: Any = None
    clinical_diagnostic_uncertainty: Any = None
    alternative_diagnosis_conflicts: Any = None
    # Pinned guideline selection and use of correlated evidence.
    vcep_specification_id: Any = None
    vcep_specification_version: Any = None
    clinical_likelihood_ratio: Any = None
    clinically_calibrated_lr: Any = None
    lr_evidence_description: Any = None
    breast_ovarian_coobservation_only: Any = None
    bp5_evidence_ids: Any = None
    other_criteria_evidence_ids: Any = None

    def __post_init__(self) -> None:
        boolean_fields = (
            "candidate_independent_disease_evidence", "candidate_specific_phenotype_evidence",
            "candidate_additional_phenotype_contribution", "alternative_variant_classification_independence",
            "alternative_inheritance_compatibility", "alternative_gene_disease_validity",
            "patient_phenotype_present", "alternative_diagnosis_explains_core_phenotypes",
            "alternative_diagnosis_explains_all_relevant_phenotypes", "dual_diagnosis_evidence",
            "blended_phenotype_evidence", "oligogenic_disease_evidence", "modifier_variant_evidence",
            "clinical_diagnostic_uncertainty", "alternative_diagnosis_conflicts",
            "clinically_calibrated_lr", "breast_ovarian_coobservation_only",
        )
        for name in boolean_fields:
            evidence = self.evidence(name)
            if evidence.state == EvidenceState.KNOWN and type(evidence.value) is not bool:
                raise ValueError(f"{name} must be Boolean when KNOWN")

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any]) -> "BP5Features":
        unknown = set(values) - {f.name for f in fields(cls)}
        if unknown:
            raise ValueError(f"Unsupported BP5 features: {', '.join(sorted(unknown))}")
        return cls(**values)

    def evidence(self, name: str) -> Evidence:
        return Evidence.parse(getattr(self, name))


@dataclass(frozen=True)
class BP5Result:
    criterion: str
    status: str
    strength: str
    decision_path: str | None
    alternative_diagnosis_id: str | None
    features_used: tuple[str, ...]
    missing_critical_features: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    guideline: dict[str, str]
    reason: str

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        for key in ("features_used", "missing_critical_features", "exclusions_evaluated"):
            result[key] = list(result[key])
        return result
