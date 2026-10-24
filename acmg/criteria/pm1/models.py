"""Typed, already-curated inputs for the PM1 decision engine."""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class State(str, Enum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"
    NOT_APPLICABLE = "NOT_APPLICABLE"


@dataclass(frozen=True)
class Evidence:
    state: State = State.UNKNOWN
    value: Any = None
    source: str | None = None

    @classmethod
    def known(cls, value: Any, source: str | None = None) -> "Evidence":
        return cls(State.KNOWN, value, source)


def evidence(value: Any) -> Evidence:
    if isinstance(value, Evidence):
        return value
    if isinstance(value, dict) and "state" in value:
        return Evidence(State(value["state"]), value.get("value"), value.get("source"))
    return Evidence() if value is None else Evidence.known(value)


@dataclass(frozen=True)
class Region:
    region_id: str
    kind: str  # HOTSPOT, CRITICAL_DOMAIN, GENERIC_DOMAIN, ENRICHMENT
    start: int
    end: int
    transcript_id_version: str
    protein_accession_version: str | None = None
    definition_source: str | None = None
    definition_version: str | None = None
    established: Evidence = field(default_factory=Evidence)
    disease_relevant: Evidence = field(default_factory=Evidence)
    mechanism_compatible: Evidence = field(default_factory=Evidence)
    pathogenic_clustering: Evidence = field(default_factory=Evidence)
    critical_function: Evidence = field(default_factory=Evidence)
    benign_variation_assessed: Evidence = field(default_factory=Evidence)
    lacks_benign_variation: Evidence = field(default_factory=Evidence)
    classification_independent: Evidence = field(default_factory=Evidence)
    classification_conflicts: Evidence = field(default_factory=Evidence)
    ascertainment_bias_resolved: Evidence = field(default_factory=Evidence)
    evidence_quality: Evidence = field(default_factory=Evidence)
    population_coverage_adequate: Evidence = field(default_factory=Evidence)
    excluded: Evidence = field(default_factory=Evidence)
    # Context is audit data; no numeric threshold is inferred from it.
    pathogenic_count: int | None = None
    benign_count: int | None = None
    independent_pathogenic_positions: int | None = None
    enrichment_statistic: float | None = None
    functional_validation_source: str | None = None
    evidence_overlap: tuple[str, ...] = ()

    @classmethod
    def from_dict(cls, data: dict) -> "Region":
        fields = dict(data)
        for name in _REGION_EVIDENCE:
            if name in fields:
                fields[name] = evidence(fields[name])
        return cls(**fields)


_REGION_EVIDENCE = (
    "established", "disease_relevant", "mechanism_compatible",
    "pathogenic_clustering", "critical_function", "benign_variation_assessed",
    "lacks_benign_variation", "classification_independent", "classification_conflicts",
    "ascertainment_bias_resolved", "evidence_quality",
    "population_coverage_adequate", "excluded",
)


@dataclass(frozen=True)
class BiologicalFeatures:
    gene_id: str | None = None
    disease_id: str | None = None
    variant_id: str | None = None
    variant_consequence: str | None = None
    transcript_id_version: str | None = None
    protein_accession_version: str | None = None
    protein_position: int | None = None
    amino_acid_change: str | None = None
    exon_number: int | None = None
    disease_molecular_mechanism: str | None = None
    variant_mechanism_compatibility: Evidence = field(default_factory=Evidence)
    gene_disease_validity: Evidence = field(default_factory=Evidence)
    mapping_confidence: Evidence = field(default_factory=Evidence)
    somatic_same_amino_acid_observations: Evidence = field(default_factory=Evidence)
    somatic_source_version: str | None = None
    pm5_strength: str | None = None
    regions: tuple[Region, ...] = ()

    @classmethod
    def from_dict(cls, data: dict) -> "BiologicalFeatures":
        fields = dict(data)
        fields["regions"] = tuple(
            r if isinstance(r, Region) else Region.from_dict(r)
            for r in fields.get("regions", ())
        )
        for name in ("variant_mechanism_compatibility", "gene_disease_validity",
                     "mapping_confidence", "somatic_same_amino_acid_observations"):
            if name in fields:
                fields[name] = evidence(fields[name])
        return cls(**fields)


@dataclass(frozen=True)
class Decision:
    criterion: str
    status: str
    strength: str
    decision_path: str | None
    region_id: str | None
    features_used: tuple[str, ...]
    missing_critical_features: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    guideline: dict[str, str] | None
    reason: str

    def to_dict(self) -> dict:
        return {
            **self.__dict__,
            "features_used": list(self.features_used),
            "missing_critical_features": list(self.missing_critical_features),
            "exclusions_evaluated": list(self.exclusions_evaluated),
        }
