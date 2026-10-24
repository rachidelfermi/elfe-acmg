"""Typed, already-curated evidence for BP2. Unknown is never treated as false."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum


class State(StrEnum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class Phase(StrEnum):
    CIS = "CIS"
    TRANS = "TRANS"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"


class Classification(StrEnum):
    PATHOGENIC = "PATHOGENIC"
    LIKELY_PATHOGENIC = "LIKELY_PATHOGENIC"
    VUS = "VUS"
    LIKELY_BENIGN = "LIKELY_BENIGN"
    BENIGN = "BENIGN"
    UNKNOWN = "UNKNOWN"


class Inheritance(StrEnum):
    AD = "AUTOSOMAL_DOMINANT"
    AR = "AUTOSOMAL_RECESSIVE"
    XL = "X_LINKED"
    OTHER = "OTHER"
    UNKNOWN = "UNKNOWN"


class DecisionStatus(StrEnum):
    ACTIVATED = "ACTIVATED"
    NOT_MET = "NOT_MET"
    UNDETERMINED = "UNDETERMINED"
    NOT_APPLICABLE = "NOT_APPLICABLE"


@dataclass(frozen=True)
class Finding:
    """A curated finding with its provenance. A known false exclusion is explicit."""

    state: State = State.UNKNOWN
    value: bool | None = None
    source: str | None = None
    method: str | None = None
    version: str | None = None
    missing_reason: str | None = None

    def __post_init__(self) -> None:
        if self.state == State.KNOWN and self.value is None:
            raise ValueError("KNOWN findings require a boolean value")
        if self.state != State.KNOWN and self.value is not None:
            raise ValueError("Only KNOWN findings may carry a value")

    @classmethod
    def yes(cls, source: str | None = None) -> Finding:
        return cls(State.KNOWN, True, source)

    @classmethod
    def no(cls, source: str | None = None) -> Finding:
        return cls(State.KNOWN, False, source)

    def is_yes(self) -> bool:
        return self.state == State.KNOWN and self.value is True

    def is_no(self) -> bool:
        return self.state == State.KNOWN and self.value is False


@dataclass(frozen=True)
class Comparison:
    variant_id: str
    gene_id: str | None = None
    classification: Classification = Classification.UNKNOWN
    classification_source: str | None = None
    review_status: str | None = None
    classification_evidence: str | None = None
    classification_independent: Finding = field(default_factory=Finding)
    observed_independently: Finding = field(default_factory=Finding)
    disease_match: Finding = field(default_factory=Finding)
    mechanism_compatible: Finding = field(default_factory=Finding)
    vcep_classified: Finding = field(default_factory=Finding)
    genotype: str | None = None
    call_reliable: Finding = field(default_factory=Finding)
    allele_balance_reliable: Finding = field(default_factory=Finding)


@dataclass(frozen=True)
class Observation:
    comparison: Comparison
    phase: Phase = Phase.UNKNOWN
    phase_method: str | None = None
    phase_reliable: Finding = field(default_factory=Finding)
    phase_conflict: Finding = field(default_factory=Finding)
    candidate_genotype: str | None = None
    case_id: str | None = None
    family_id: str | None = None
    study_id: str | None = None
    affected: Finding = field(default_factory=Finding)
    diagnosis_matches: Finding = field(default_factory=Finding)
    comparison_explains_disease: Finding = field(default_factory=Finding)
    clinical_evidence_reliable: Finding = field(default_factory=Finding)
    phenotype_not_more_severe: Finding = field(default_factory=Finding)
    relevant_biallelic_effect_excluded: Finding = field(default_factory=Finding)
    evidence_source: str | None = None


@dataclass(frozen=True)
class BiologicalFeatures:
    candidate_variant_id: str
    gene_id: str
    disease_id: str
    inheritance: Inheritance
    gene_disease_valid: Finding
    observations: tuple[Observation, ...]
    disease_mechanism_compatible: Finding = field(default_factory=Finding)
    full_penetrance: Finding = field(default_factory=Finding)
    candidate_call_reliable: Finding = field(default_factory=Finding)
    candidate_allele_balance_reliable: Finding = field(default_factory=Finding)
    synergy_or_complex_allele: Finding = field(default_factory=Finding)
    haplotype_specific_pathogenicity: Finding = field(default_factory=Finding)
    candidate_independent_disease_effect: Finding = field(default_factory=Finding)
    dual_diagnosis: Finding = field(default_factory=Finding)
    evidence_conflict: Finding = field(default_factory=Finding)
    other_criteria_from_same_observation: tuple[str, ...] = ()
    guideline_version: str | None = None


@dataclass(frozen=True)
class GuidelineRef:
    name: str
    version: str
    rule_id: str
    url: str


@dataclass(frozen=True)
class BP2Result:
    criterion: str
    status: DecisionStatus
    strength: str
    decision_path: str | None
    comparison_variant_id: str | None
    phase_status: str | None
    features_used: tuple[str, ...]
    missing_critical_features: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    guideline: GuidelineRef
    reason: str
    observation_ids: tuple[str, ...] = ()
    criteria_overlap: tuple[str, ...] = ()
