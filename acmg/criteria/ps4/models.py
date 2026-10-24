"""Typed PS4 inputs. None and omitted fields mean unknown, never zero or false."""

from dataclasses import asdict, dataclass, field
from typing import Any, Literal

State = Literal["KNOWN", "UNKNOWN", "CONFLICTING", "NOT_APPLICABLE"]


@dataclass(frozen=True)
class Evidence:
    state: State = "UNKNOWN"
    value: Any = None
    source: str | None = None

    @classmethod
    def parse(cls, item: Any) -> "Evidence":
        if isinstance(item, cls):
            return item
        if isinstance(item, dict) and "state" in item:
            return cls(**item)
        return cls("UNKNOWN" if item is None else "KNOWN", item)

    def known(self, value: Any = True) -> bool:
        return self.state == "KNOWN" and self.value == value


@dataclass(frozen=True)
class Study:
    study_id: str
    cohort_id: str | None = None
    variant_id: str | None = None
    disease_id: str | None = None
    case_carriers: int | None = None
    case_noncarriers: int | None = None
    control_carriers: int | None = None
    control_noncarriers: int | None = None
    odds_ratio: float | None = None
    relative_risk: float | None = None
    confidence_interval: tuple[float, float] | None = None
    confidence_level: float | None = None
    p_value: float | None = None
    statistic_source: str | None = None
    statistic_method: str | None = None
    cases_affected: Evidence = field(default_factory=Evidence)
    controls_valid: Evidence = field(default_factory=Evidence)
    variant_calls_valid: Evidence = field(default_factory=Evidence)
    ancestry_comparable: Evidence = field(default_factory=Evidence)
    stratification_adjusted: Evidence = field(default_factory=Evidence)
    technical_comparable: Evidence = field(default_factory=Evidence)
    confounding_resolved: Evidence = field(default_factory=Evidence)
    independent_cases: Evidence = field(default_factory=Evidence)
    overlapping_cohort: Evidence = field(default_factory=Evidence)
    founder_effect: Evidence = field(default_factory=Evidence)
    batch_effect: Evidence = field(default_factory=Evidence)
    case_ids: tuple[str, ...] = ()
    control_source: str | None = None

    @classmethod
    def parse(cls, obj: "Study | dict") -> "Study":
        if isinstance(obj, cls):
            return obj
        d = dict(obj)
        for name in cls.__dataclass_fields__:
            if name in d and name in EVIDENCE_STUDY_FIELDS:
                d[name] = Evidence.parse(d[name])
        if "confidence_interval" in d and d["confidence_interval"] is not None:
            d["confidence_interval"] = tuple(d["confidence_interval"])
        if "case_ids" in d:
            d["case_ids"] = tuple(d["case_ids"])
        return cls(**d)


EVIDENCE_STUDY_FIELDS = {
    "cases_affected", "controls_valid", "variant_calls_valid", "ancestry_comparable",
    "stratification_adjusted", "technical_comparable", "confounding_resolved",
    "independent_cases", "overlapping_cohort", "founder_effect", "batch_effect",
}


@dataclass(frozen=True)
class Proband:
    case_id: str
    study_id: str
    family_id: str | None = None
    clinically_assessed: Evidence = field(default_factory=Evidence)
    phenotype: str | None = None
    has_other_myoc_vus_lp_p: Evidence = field(default_factory=Evidence)
    overlaps_other_criterion: Evidence = field(default_factory=Evidence)

    @classmethod
    def parse(cls, obj: "Proband | dict") -> "Proband":
        if isinstance(obj, cls):
            return obj
        d = dict(obj)
        for key in ("clinically_assessed", "has_other_myoc_vus_lp_p", "overlaps_other_criterion"):
            if key in d:
                d[key] = Evidence.parse(d[key])
        return cls(**d)


@dataclass(frozen=True)
class PS4Input:
    variant_id: str
    gene_id: str
    disease_id: str
    gene_disease_validity: Evidence = field(default_factory=Evidence)
    studies: tuple[Study, ...] = ()
    probands: tuple[Proband, ...] = ()
    pm2_met: Evidence = field(default_factory=Evidence)
    studies_independent: Evidence = field(default_factory=Evidence)
    specification: str | None = None

    @classmethod
    def parse(cls, obj: "PS4Input | dict") -> "PS4Input":
        if isinstance(obj, cls):
            return obj
        d = dict(obj)
        d["gene_disease_validity"] = Evidence.parse(d.get("gene_disease_validity"))
        d["pm2_met"] = Evidence.parse(d.get("pm2_met"))
        d["studies_independent"] = Evidence.parse(d.get("studies_independent"))
        d["studies"] = tuple(Study.parse(x) for x in d.get("studies", ()))
        d["probands"] = tuple(Proband.parse(x) for x in d.get("probands", ()))
        return cls(**d)


@dataclass(frozen=True)
class PS4Result:
    status: Literal["ACTIVATED", "NOT_MET", "UNDETERMINED", "NOT_APPLICABLE"]
    strength: Literal["VERY_STRONG", "STRONG", "MODERATE", "SUPPORTING", "NONE"] = "NONE"
    decision_path: str | None = None
    studies_used: tuple[str, ...] = ()
    eligible_case_count: int | None = None
    odds_ratio: float | None = None
    relative_risk: float | None = None
    confidence_interval: tuple[float, float] | None = None
    p_value: float | None = None
    features_used: tuple[str, ...] = ()
    missing_critical_features: tuple[str, ...] = ()
    exclusions_evaluated: tuple[str, ...] = ()
    conflicting_evidence: tuple[str, ...] = ()
    guideline: dict[str, str] = field(default_factory=dict)
    reason: str = ""
    criterion: str = "PS4"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
