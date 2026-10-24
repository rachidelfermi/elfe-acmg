"""Typed input and output for the BS1 criterion only."""

from dataclasses import dataclass, field
from enum import StrEnum


class EvidenceState(StrEnum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class Status(StrEnum):
    ACTIVATED = "ACTIVATED"
    NOT_MET = "NOT_MET"
    UNDETERMINED = "UNDETERMINED"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class Strength(StrEnum):
    STRONG = "STRONG"
    MODERATE = "MODERATE"
    SUPPORTING = "SUPPORTING"
    NONE = "NONE"


@dataclass(frozen=True)
class VariantIdentity:
    assembly: str
    chromosome: str
    position: int
    reference: str
    alternate: str
    normalized: bool = True
    variant_type: str = "SNV"  # SNV or INDEL; other classes need a specific rule

    def __post_init__(self):
        if self.position < 1 or not self.reference or not self.alternate:
            raise ValueError("A complete positive genomic coordinate and alleles are required")


@dataclass(frozen=True)
class FrequencyRecord:
    variant: VariantIdentity
    database: str
    version: str
    metric: str  # POPMAX_FAF, ANCESTRY_FAF, GLOBAL_AF, ANCESTRY_AF
    frequency: float | None
    ancestry: str | None = None
    allele_count: int | None = None
    allele_number: int | None = None
    filter_status: str = "UNKNOWN"
    quality: EvidenceState = EvidenceState.UNKNOWN
    sampling: EvidenceState = EvidenceState.UNKNOWN
    founder_population: bool = False
    source: str | None = None
    platform: str | None = None  # exome or genome
    depth: float | None = None
    technical_limitation: bool = False

    def __post_init__(self):
        if self.frequency is not None and not 0 <= self.frequency <= 1:
            raise ValueError("Frequency must be in [0, 1]")
        if self.allele_count is not None and self.allele_count < 0:
            raise ValueError("Allele count must be nonnegative")
        if self.allele_number is not None and self.allele_number <= 0:
            raise ValueError("Allele number must be positive")
        if self.allele_count is not None and self.allele_number is not None and self.allele_count > self.allele_number:
            raise ValueError("Allele count exceeds allele number")


@dataclass(frozen=True)
class GuidelineRule:
    name: str
    version: str
    rule_id: str
    source_url: str
    gene: str
    disease: str
    inheritance: str
    threshold: float
    metric: str
    operator: str  # > or >=, exactly as published
    strength: Strength = Strength.STRONG
    database: str | None = None
    minimum_allele_number: int | None = None
    minimum_exome_depth: float | None = None
    eligible_ancestries: frozenset[str] | None = None
    excluded_variant_types: frozenset[str] = frozenset({"SV", "CNV", "MITOCHONDRIAL"})
    requires_non_founder: bool = True
    verified: bool = False  # custom rules require an explicit attestation

    def __post_init__(self):
        if not 0 < self.threshold <= 1 or self.operator not in {">", ">="}:
            raise ValueError("Invalid threshold or comparison operator")
        if not self.source_url or not self.gene or not self.disease:
            raise ValueError("A rule requires a source and exact gene/disease scope")


@dataclass(frozen=True)
class BS1Features:
    variant: VariantIdentity
    gene: str
    disease: str
    inheritance: str
    frequencies: tuple[FrequencyRecord, ...] = ()
    guideline: GuidelineRule | None = None
    variant_mapping: EvidenceState = EvidenceState.KNOWN
    gene_disease_validity: EvidenceState = EvidenceState.KNOWN
    exception_status: EvidenceState = EvidenceState.UNKNOWN
    exception_source: str | None = None
    exception_disease_scope: str | None = None
    conflicting_evidence: bool = False
    ba1_established: bool = False
    pm2_established: bool = False
    bs2_established: bool = False
    # These are context flags, not automatic vetoes. A verified exception is.
    founder_variant: bool = False
    reduced_penetrance: bool = False
    late_onset: bool = False
    risk_allele: bool = False
    hypomorphic_allele: bool = False
    prevalence: float | None = None
    penetrance: float | None = None


@dataclass(frozen=True)
class BS1Result:
    criterion: str
    status: Status
    strength: Strength
    decision_path: str | None
    frequency_used: float | None
    threshold_used: float | None
    features_used: tuple[str, ...] = ()
    missing_critical_features: tuple[str, ...] = ()
    exclusions_evaluated: tuple[str, ...] = ()
    guideline: dict[str, str] | None = None
    reason: str = ""
