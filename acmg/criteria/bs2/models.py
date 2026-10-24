"""Typed, explicit input and output structures for BS2."""

from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Any


class EvidenceState(str, Enum):
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
            return raw
        if isinstance(raw, dict):
            state = EvidenceState(raw.get("state", "KNOWN" if "value" in raw else "UNKNOWN"))
            value = raw.get("value")
            if state == EvidenceState.KNOWN and value is None:
                raise ValueError("KNOWN evidence requires a value")
            return cls(value if state == EvidenceState.KNOWN else None, state,
                       raw.get("source"), raw.get("version"), raw.get("quality"),
                       raw.get("missing_reason"))
        if raw is None:
            return cls()
        return cls(raw, EvidenceState.KNOWN)


@dataclass(frozen=True)
class DiseaseContext:
    gene_id: str | None = None
    disease_id: str | None = None
    gene_disease_validity: Evidence = field(default_factory=Evidence)
    inheritance_mode: Evidence = field(default_factory=Evidence)
    molecular_mechanism_compatible: Evidence = field(default_factory=Evidence)
    full_penetrance: Evidence = field(default_factory=Evidence)
    early_onset: Evidence = field(default_factory=Evidence)
    age_dependent_penetrance: Evidence = field(default_factory=Evidence)
    variable_expressivity: Evidence = field(default_factory=Evidence)
    relevant_subclinical_features: Evidence = field(default_factory=Evidence)


@dataclass(frozen=True)
class HealthyObservation:
    individual_id: str
    genotype: Evidence = field(default_factory=Evidence)
    genotype_reliable: Evidence = field(default_factory=Evidence)
    confirmed_unaffected: Evidence = field(default_factory=Evidence)
    adult: Evidence = field(default_factory=Evidence)
    beyond_expected_onset: Evidence = field(default_factory=Evidence)
    clinical_assessment_adequate: Evidence = field(default_factory=Evidence)
    independent: Evidence = field(default_factory=Evidence)
    maternal_allele: Evidence = field(default_factory=Evidence)
    artifact: Evidence = field(default_factory=Evidence)
    deletion_mimicking_homozygosity: Evidence = field(default_factory=Evidence)
    mosaicism: Evidence = field(default_factory=Evidence)
    protective_modifier: Evidence = field(default_factory=Evidence)
    genetic_rescue: Evidence = field(default_factory=Evidence)
    subclinical_manifestation: Evidence = field(default_factory=Evidence)
    phenotype_misclassification: Evidence = field(default_factory=Evidence)
    source_type: str | None = None
    overlaps_population_evidence: bool = False
    related_group: str | None = None


@dataclass(frozen=True)
class BS2Input:
    variant_id: str | None = None
    disease: DiseaseContext = field(default_factory=DiseaseContext)
    observations: tuple[HealthyObservation, ...] = ()
    variant_reliable: Evidence = field(default_factory=Evidence)
    vcep_applicable: Evidence = field(default_factory=Evidence)
    vcep_specification_version: str | None = None


@dataclass(frozen=True)
class BS2Result:
    criterion: str
    status: str
    strength: str
    decision_path: str | None
    eligible_healthy_observations: int
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
