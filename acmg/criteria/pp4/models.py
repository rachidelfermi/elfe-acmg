"""Typed, already-extracted evidence for PP4. Unknown is never interpreted as false."""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import Generic, TypeVar

T = TypeVar("T")


class State(StrEnum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"
    NOT_APPLICABLE = "NOT_APPLICABLE"


@dataclass(frozen=True)
class Evidence(Generic[T]):
    state: State = State.UNKNOWN
    value: T | None = None
    source: str = ""
    method: str = ""
    version: str = ""
    quality: str = ""
    missing_reason: str = ""

    def __post_init__(self) -> None:
        if self.state == State.KNOWN and self.value is None:
            raise ValueError("KNOWN evidence needs a value")
        if self.state != State.KNOWN and self.value is not None:
            raise ValueError("Only KNOWN evidence may have a value")


def known(value: T, *, source: str = "") -> Evidence[T]:
    return Evidence(State.KNOWN, value, source=source)


@dataclass(frozen=True)
class PP4Features:
    variant_id: str = ""
    gene_id: str = ""
    disease_id: str = ""
    gene_disease_validity: Evidence[str] = field(default_factory=Evidence)
    inheritance_mode: Evidence[str] = field(default_factory=Evidence)
    genotype_compatible: Evidence[bool] = field(default_factory=Evidence)
    phenotype_compatible: Evidence[bool] = field(default_factory=Evidence)
    phenotype_highly_specific: Evidence[bool] = field(default_factory=Evidence)
    family_history_highly_specific: Evidence[bool] = field(default_factory=Evidence)
    single_genetic_etiology: Evidence[bool] = field(default_factory=Evidence)
    alternative_diagnosis_confirmed: Evidence[bool] = field(default_factory=Evidence)
    phenotype_conflict: Evidence[bool] = field(default_factory=Evidence)
    age_sex_applicable: Evidence[bool] = field(default_factory=Evidence)
    evidence_reliable: Evidence[bool] = field(default_factory=Evidence)
    family_evidence_used_for_pp1: Evidence[bool] = field(default_factory=Evidence)
    phenotype_evidence_used_for_ps2_pm6: Evidence[bool] = field(default_factory=Evidence)
    # Context and audit fields. These do not independently activate PP4.
    patient_hpo_terms: tuple[str, ...] = ()
    patient_absent_hpo_terms: tuple[str, ...] = ()
    clinical_findings: tuple[str, ...] = ()
    patient_age_years: float | None = None
    patient_sex: str = ""
    hpo_similarity_score: float | None = None
    phenotype_match_rank: int | None = None
    phenotype_similarity_method: str = ""
    evidence_overlap: tuple[str, ...] = ()
    # ADA VCEP v2.2: individual, documented components of its PP4 points table.
    scid_diagnostic_criteria_met: Evidence[bool] = field(default_factory=Evidence)
    scid_sequencing_no_alternative: Evidence[bool] = field(default_factory=Evidence)
    scid_family_history: Evidence[bool] = field(default_factory=Evidence)
    ada_enzyme_activity_percent: Evidence[float] = field(default_factory=Evidence)
    ada_dado_above_reference: Evidence[bool] = field(default_factory=Evidence)
    ada_dado_pre_treatment_or_nontransfused: Evidence[bool] = field(default_factory=Evidence)
    ada_phenotype_corrected_by_treatment: Evidence[bool] = field(default_factory=Evidence)
    ada_cnv_testing_done: Evidence[bool] = field(default_factory=Evidence)
    scid_t_minus_b_minus_nk_minus: Evidence[bool] = field(default_factory=Evidence)
    # Glanzmann thrombasthenia VCEP v2.1 clinical and laboratory requirements.
    gt_clinical_diagnosis: Evidence[bool] = field(default_factory=Evidence)
    gt_mucocutaneous_bleeding: Evidence[bool] = field(default_factory=Evidence)
    gt_normal_ristocetin: Evidence[bool] = field(default_factory=Evidence)
    gt_minimal_aggregation_at_least_two_agonists: Evidence[bool] = field(default_factory=Evidence)
    gt_aIIbbeta3_expression_below_25_percent: Evidence[bool] = field(default_factory=Evidence)
    gt_impaired_fibrinogen_binding: Evidence[bool] = field(default_factory=Evidence)
    gt_both_genes_fully_sequenced: Evidence[bool] = field(default_factory=Evidence)


@dataclass(frozen=True)
class Guideline:
    name: str
    version: str
    rule_id: str
    url: str


@dataclass(frozen=True)
class PP4Result:
    criterion: str
    status: str
    strength: str
    decision_path: str
    features_used: tuple[str, ...]
    missing_critical_features: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    guideline: Guideline | None
    reason: str
