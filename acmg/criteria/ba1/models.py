"""Typed, already-extracted evidence. Unknown is never interpreted as zero."""

from dataclasses import dataclass, field
from enum import StrEnum


class EvidenceState(StrEnum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"
    NOT_APPLICABLE = "NOT_APPLICABLE"


@dataclass(frozen=True)
class VariantIdentity:
    variant_id: str | None = None
    assembly: str | None = None
    chromosome: str | None = None
    position: int | None = None
    reference: str | None = None
    alternate: str | None = None
    normalized: bool | None = None
    variant_type: str = "SNV"  # SNV, INDEL, SV, CNV, MITOCHONDRIAL
    mapping_confidence: EvidenceState = EvidenceState.UNKNOWN
    clingen_allele_id: str | None = None
    sex_chromosome_ploidy: str | None = None
    pseudoautosomal_region: bool | None = None
    mitochondrial_heteroplasmy: float | None = None


@dataclass(frozen=True)
class PopulationRecord:
    database: str
    version: str
    assembly: str
    chromosome: str
    position: int
    reference: str
    alternate: str
    normalized: bool
    ancestry: str
    continental: bool | None = None
    allele_frequency: float | None = None
    allele_count: int | None = None
    allele_number: int | None = None
    filtering_allele_frequency: float | None = None
    filter_pass: bool | None = None
    quality: EvidenceState = EvidenceState.UNKNOWN
    mapping_ambiguous: bool | None = None
    callability_adequate: bool | None = None
    founder_or_bottlenecked: bool | None = None
    sampling_bias: bool | None = None
    technical_limitation: bool | None = None
    population_dataset: str | None = None  # e.g. gnomAD SV
    detection_sensitivity_adequate: bool | None = None


@dataclass(frozen=True)
class ExceptionCheck:
    state: EvidenceState = EvidenceState.UNKNOWN
    is_exception: bool | None = None
    source: str | None = None
    version: str | None = None
    disease_id: str | None = None
    reason: str | None = None


@dataclass(frozen=True)
class BA1Evidence:
    variant: VariantIdentity
    populations: tuple[PopulationRecord, ...] = ()
    gene_id: str | None = None
    disease_id: str | None = None
    exception_check: ExceptionCheck = field(default_factory=ExceptionCheck)
    guideline_modification_check: EvidenceState = EvidenceState.UNKNOWN
    # A KNOWN modification is selected only from bundled verified specifications.
    known_modification: bool | None = None
    pathogenic_evidence_conflict: bool = False
    frequency_conflict: bool = False
    reduced_penetrance: bool = False
    late_onset: bool = False
    risk_allele: bool = False
    hypomorphic_allele: bool = False
    founder_pathogenic_context: bool = False
    evidence_source: str | None = None
    evidence_version: str | None = None
    missing_reason: str | None = None


@dataclass(frozen=True)
class BA1Result:
    criterion: str
    status: str
    strength: str | None
    decision_path: str | None
    frequency_used: float | None
    threshold_used: float | None
    features_used: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    guideline: dict[str, str] | None
    reason: str
