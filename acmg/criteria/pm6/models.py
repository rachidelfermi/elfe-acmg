"""Typed, already-curated evidence for a PM6 decision."""

from __future__ import annotations

from dataclasses import dataclass, field
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
    NONE = "NONE"
    SUPPORTING = "SUPPORTING"
    MODERATE = "MODERATE"
    STRONG = "STRONG"
    VERY_STRONG = "VERY_STRONG"


class Phenotype(StrEnum):
    HIGHLY_SPECIFIC = "HIGHLY_SPECIFIC"
    CONSISTENT = "CONSISTENT"
    HETEROGENEOUS = "HETEROGENEOUS"
    INCONSISTENT = "INCONSISTENT"


@dataclass(frozen=True)
class Observation:
    case_id: str | None = None
    family_id: str | None = None
    study_ids: tuple[str, ...] = ()
    proband_variant_present: bool | None = None
    maternal_variant_absent: bool | None = None
    paternal_variant_absent: bool | None = None
    maternal_test_adequate: bool | None = None
    paternal_test_adequate: bool | None = None
    maternity_confirmed: bool | None = None
    paternity_confirmed: bool | None = None
    sample_identity_valid: bool | None = None
    affected: bool | None = None
    phenotype: Phenotype | None = None
    independent: bool | None = None
    related_to_other_case: bool = False
    duplicate_of: str | None = None
    family_history_compatible: bool | None = None
    parental_mosaicism: bool | None = None
    proband_mosaicism: bool | None = None
    apparent_germline_mosaicism: bool = False
    technical_artifact: bool | None = None
    mapping_ambiguity: bool | None = None
    contamination: bool | None = None
    conflicting_segregation: bool | None = None
    alternative_diagnosis: bool | None = None
    already_counted_ps2: bool = False
    evidence_conflict: bool = False
    # Curator-supplied technical conclusions; no invented numeric cutoffs.
    proband_call_adequate: bool | None = None
    evidence_source: str | None = None
    evidence_method: str | None = None
    evidence_version: str | None = None
    missing_reason: str | None = None
    vcep_facts: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class PM6Input:
    variant_id: str | None = None
    gene_id: str | None = None
    disease_id: str | None = None
    gene_disease_valid: bool | None = None
    inheritance_mode: str | None = None
    mechanism_compatible: bool | None = None
    second_pathogenic_allele: bool | None = None
    vcep_profile: str | None = None
    observations: tuple[Observation, ...] = ()


@dataclass(frozen=True)
class ObservationResult:
    case_id: str | None
    status: Status
    category: str | None
    phenotype: str | None
    points: float
    missing_critical_features: tuple[str, ...] = ()
    exclusions: tuple[str, ...] = ()
    rule_id: str = "SVI-v1.1-table-1"

    def to_dict(self) -> dict[str, Any]:
        return {"case_id": self.case_id, "status": self.status.value,
                "category": self.category, "phenotype": self.phenotype,
                "points": self.points,
                "missing_critical_features": list(self.missing_critical_features),
                "exclusions": list(self.exclusions), "rule_id": self.rule_id}


@dataclass(frozen=True)
class PM6Result:
    status: Status
    strength: Strength
    total_points: float
    eligible_observations: int
    observations: tuple[ObservationResult, ...]
    features_used: tuple[str, ...]
    missing_critical_features: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    guideline: dict[str, str]
    reason: str
    delegation: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {"criterion": "PM6", "status": self.status.value,
                "strength": self.strength.value, "total_points": self.total_points,
                "eligible_observations": self.eligible_observations,
                "observations": [item.to_dict() for item in self.observations],
                "features_used": list(self.features_used),
                "missing_critical_features": list(self.missing_critical_features),
                "exclusions_evaluated": list(self.exclusions_evaluated),
                "guideline": self.guideline, "reason": self.reason,
                "delegation": self.delegation}
