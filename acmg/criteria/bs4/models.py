"""Typed, already-curated inputs for the BS4 criterion."""

from dataclasses import asdict, dataclass, field
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
    method: str | None = None
    version: str | None = None
    missing_reason: str | None = None


def evidence(value: Any) -> Evidence:
    if isinstance(value, Evidence):
        return value
    if isinstance(value, dict) and "state" in value:
        return Evidence(**{**value, "state": State(value["state"])})
    if value is None:
        return Evidence()
    return Evidence(State.KNOWN, value)


@dataclass(frozen=True)
class Observation:
    family_id: str
    individual_id: str
    affected: Evidence = field(default_factory=Evidence)
    variant_present: Evidence = field(default_factory=Evidence)
    genotype_reliable: Evidence = field(default_factory=Evidence)
    phenotype_reliable: Evidence = field(default_factory=Evidence)
    phenotype_consistent: Evidence = field(default_factory=Evidence)
    biological_relationship_confirmed: Evidence = field(default_factory=Evidence)
    phenocopy_excluded: Evidence = field(default_factory=Evidence)
    heterogeneity_excluded: Evidence = field(default_factory=Evidence)
    alternative_diagnosis: Evidence = field(default_factory=Evidence)
    phenocopy_evidence: Evidence = field(default_factory=Evidence)
    genetic_heterogeneity: Evidence = field(default_factory=Evidence)
    locus_heterogeneity: Evidence = field(default_factory=Evidence)
    phenotype_misclassification: Evidence = field(default_factory=Evidence)
    mosaicism_explanation: Evidence = field(default_factory=Evidence)
    genotyping_artifact: Evidence = field(default_factory=Evidence)
    copy_number_conflict: Evidence = field(default_factory=Evidence)
    phasing_conflict: Evidence = field(default_factory=Evidence)
    recombination_explanation: Evidence = field(default_factory=Evidence)
    sample_identity_conflict: Evidence = field(default_factory=Evidence)
    pedigree_conflict: Evidence = field(default_factory=Evidence)
    phenotype_score: Evidence = field(default_factory=Evidence)
    evidence_quality: Evidence = field(default_factory=Evidence)
    source: str | None = None
    study_id: str | None = None

    @classmethod
    def from_mapping(cls, data: dict) -> "Observation":
        fields = cls.__dataclass_fields__
        values = {k: evidence(v) if k in EVIDENCE_FIELDS else v
                  for k, v in data.items() if k in fields}
        return cls(**values)


EVIDENCE_FIELDS = {name for name, definition in Observation.__dataclass_fields__.items()
                   if definition.type == Evidence or str(definition.type) == "Evidence"}


@dataclass(frozen=True)
class BS4Input:
    variant_id: str | None = None
    gene_id: str | None = None
    disease_id: str | None = None
    vcep_specification_version: str | None = None
    gene_disease_validity: Evidence = field(default_factory=Evidence)
    inheritance_compatible: Evidence = field(default_factory=Evidence)
    variant_reliable: Evidence = field(default_factory=Evidence)
    variant_mapping_confidence: Evidence = field(default_factory=Evidence)
    variant_call_quality: Evidence = field(default_factory=Evidence)
    observations: tuple[Observation, ...] = ()
    quantitative_likelihood_ratio: Evidence = field(default_factory=Evidence)
    quantitative_validated: Evidence = field(default_factory=Evidence)
    quantitative_model_appropriate: Evidence = field(default_factory=Evidence)
    quantitative_independent: Evidence = field(default_factory=Evidence)
    family_overlap: Evidence = field(default_factory=Evidence)
    pp1_conflict: Evidence = field(default_factory=Evidence)
    other_criteria_overlap: tuple[str, ...] = ()

    @classmethod
    def from_mapping(cls, data: dict) -> "BS4Input":
        keys = cls.__dataclass_fields__
        values = {k: evidence(v) if k in CONTEXT_EVIDENCE_FIELDS else v
                  for k, v in data.items() if k in keys and k != "observations"}
        values["observations"] = tuple(
            item if isinstance(item, Observation) else Observation.from_mapping(item)
            for item in data.get("observations", ())
        )
        return cls(**values)


CONTEXT_EVIDENCE_FIELDS = {"gene_disease_validity", "inheritance_compatible",
                           "variant_reliable", "variant_mapping_confidence",
                           "variant_call_quality", "quantitative_likelihood_ratio",
                           "quantitative_validated", "quantitative_model_appropriate",
                           "quantitative_independent", "family_overlap", "pp1_conflict"}


@dataclass(frozen=True)
class Decision:
    criterion: str
    status: str
    strength: str
    decision_path: str
    informative_nonsegregation_events: int
    families_used: tuple[str, ...]
    features_used: tuple[str, ...]
    missing_critical_features: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    guideline: dict[str, str]
    reason: str
    conflicts: tuple[str, ...] = ()
    overlaps: tuple[str, ...] = ()

    def to_dict(self) -> dict:
        result = asdict(self)
        for key in ("families_used", "features_used", "missing_critical_features",
                    "exclusions_evaluated", "conflicts", "overlaps"):
            result[key] = list(result[key])
        return result
