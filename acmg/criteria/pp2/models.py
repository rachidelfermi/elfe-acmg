"""Typed inputs and deterministic PP2 output."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import StrEnum
from typing import Any


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
class Assessment:
    """Curated assertion, with enough provenance to reject circular or conflicting support."""

    state: EvidenceState = EvidenceState.UNKNOWN
    value: bool | None = None
    source: str | None = None
    version: str | None = None
    quality: str | None = None  # HIGH, MODERATE, LOW, or UNKNOWN
    independent_of_pp2: bool | None = None
    disease_match: bool | None = None
    ascertainment_bias: bool | None = None
    overlaps_criteria: tuple[str, ...] = ()
    missing_reason: str | None = None

    def __post_init__(self) -> None:
        if self.state == EvidenceState.KNOWN and self.value is None:
            raise ValueError("KNOWN assessments require a boolean value")
        if self.state != EvidenceState.KNOWN and self.value is not None:
            raise ValueError("Only KNOWN assessments may carry a value")


@dataclass(frozen=True)
class BiologicalFeatures:
    variant_id: str | None = None
    gene_id: str | None = None
    disease_id: str | None = None
    variant_consequence: str | None = None
    coding_hgvs: str | None = None
    protein_hgvs: str | None = None
    transcript_id_version: str | None = None
    protein_accession_version: str | None = None
    relevant_transcript: Assessment = field(default_factory=Assessment)
    consequence_reliable: Assessment = field(default_factory=Assessment)
    alternative_splice_consequence: Assessment = field(default_factory=Assessment)
    mapping_reliable: Assessment = field(default_factory=Assessment)
    gene_disease_valid: Assessment = field(default_factory=Assessment)
    inheritance_mode: str | None = None
    disease_molecular_mechanism: str | None = None
    missense_mechanism_common: Assessment = field(default_factory=Assessment)
    benign_missense_uncommon: Assessment = field(default_factory=Assessment)
    mechanism_compatible: Assessment = field(default_factory=Assessment)
    pathogenic_missense_fraction: float | None = None
    pathogenic_missense_variant_count: int | None = None
    independent_pathogenic_missense_positions: int | None = None
    benign_missense_variant_count: int | None = None
    benign_missense_fraction: float | None = None
    population_missense_variant_count: int | None = None
    missense_constraint_z_score: float | None = None
    missense_observed_expected_ratio: float | None = None
    regional_missense_constraint: float | None = None
    constraint_source: str | None = None
    constraint_version: str | None = None
    disease_specific_variant_spectrum: str | None = None
    pp2_specification_id: str | None = None
    requested_strength: Strength | None = None


@dataclass(frozen=True)
class Guideline:
    name: str
    version: str
    rule_id: str
    url: str


@dataclass(frozen=True)
class Decision:
    criterion: str
    status: Status
    strength: Strength
    decision_path: str
    features_used: tuple[str, ...]
    missing_critical_features: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    guideline: Guideline
    reason: str

    def to_dict(self) -> dict[str, Any]:
        def plain(value: Any) -> Any:
            if isinstance(value, StrEnum):
                return value.value
            if isinstance(value, dict):
                return {key: plain(item) for key, item in value.items()}
            if isinstance(value, (list, tuple)):
                return [plain(item) for item in value]
            return value

        return plain(asdict(self))
