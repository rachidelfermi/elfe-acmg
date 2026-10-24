"""Typed input and output for a supplied, disease-specific BP1 assessment."""

from dataclasses import asdict, dataclass, field
from enum import StrEnum
from typing import Any, Mapping


class EvidenceState(StrEnum):
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

    @classmethod
    def parse(cls, raw: Any) -> "Evidence":
        if isinstance(raw, cls):
            return raw
        if isinstance(raw, Mapping) and "state" in raw:
            state = EvidenceState(raw["state"])
            value = raw.get("value")
            if state != EvidenceState.KNOWN and value is not None:
                raise ValueError("Only KNOWN evidence may have a value")
            return cls(value, state, raw.get("source"), raw.get("version"))
        if raw is None:
            return cls()
        return cls(raw, EvidenceState.KNOWN)


@dataclass(frozen=True)
class BP1Features:
    """Already extracted features. Missing fields are explicitly unknown.

    Identifiers and versioned evidence provenance are kept separate from biological
    conclusions. Counts and constraint measures are context only; they cannot
    create a qualitative 'primarily truncating' finding by themselves.
    """

    variant_id: Any = None
    gene_id: Any = None
    disease_id: Any = None
    variant_consequence: Any = None
    variant_type: Any = None
    transcript_id_version: Any = None
    biologically_relevant_transcript: Any = None
    missense_consequence_confidence: Any = None
    alternative_splice_consequence: Any = None
    gene_disease_validity: Any = None
    disease_inheritance_mode: Any = None
    disease_molecular_mechanism: Any = None
    truncating_predominant: Any = None
    missense_common_mechanism: Any = None
    loss_of_function_disease_mechanism: Any = None
    missense_gain_of_function_mechanism: Any = None
    missense_dominant_negative_mechanism: Any = None
    missense_specific_disease_subtype: Any = None
    exceptions_assessed: Any = None
    candidate_region_missense_pathogenicity: Any = None
    candidate_residue_missense_pathogenicity: Any = None
    candidate_is_known_pathogenic_missense: Any = None
    candidate_splice_effect: Any = None
    candidate_splice_prediction: Any = None
    alternative_transcript_consequence: Any = None
    mechanism_evidence_quality: Any = None
    mechanism_evidence_source: Any = None
    mechanism_evidence_version: Any = None
    spectrum_truncating_predominant: Any = None
    spectrum_evidence_quality: Any = None
    spectrum_classification_independent: Any = None
    spectrum_unique_variants: Any = None
    spectrum_disease_match: Any = None
    spectrum_ascertainment_bias_addressed: Any = None
    pathogenic_truncating_variant_count: Any = None
    pathogenic_missense_variant_count: Any = None
    pathogenic_truncating_fraction: Any = None
    pathogenic_missense_fraction: Any = None
    missense_constraint: Any = None
    population_missense_tolerance: Any = None
    benign_missense_variant_count: Any = None
    prior_criteria: Any = None
    vcep_specification_id: Any = None
    vcep_specification_version: Any = None

    def __post_init__(self) -> None:
        boolean_fields = (
            "biologically_relevant_transcript", "gene_disease_validity", "truncating_predominant",
            "missense_common_mechanism", "loss_of_function_disease_mechanism",
            "missense_gain_of_function_mechanism", "missense_dominant_negative_mechanism",
            "missense_specific_disease_subtype", "exceptions_assessed",
            "candidate_region_missense_pathogenicity", "candidate_residue_missense_pathogenicity",
            "candidate_is_known_pathogenic_missense", "spectrum_truncating_predominant",
            "spectrum_classification_independent", "spectrum_unique_variants", "spectrum_disease_match",
            "spectrum_ascertainment_bias_addressed",
        )
        for name in boolean_fields:
            evidence = self.evidence(name)
            if evidence.state == EvidenceState.KNOWN and type(evidence.value) is not bool:
                raise ValueError(f"{name} must be a Boolean when KNOWN")

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any]) -> "BP1Features":
        allowed = cls.__dataclass_fields__
        unknown = set(values) - set(allowed)
        if unknown:
            raise ValueError(f"Unsupported BP1 features: {', '.join(sorted(unknown))}")
        return cls(**values)

    def evidence(self, name: str) -> Evidence:
        return Evidence.parse(getattr(self, name))


@dataclass(frozen=True)
class BP1Result:
    criterion: str
    status: str
    strength: str
    decision_path: str | None
    features_used: tuple[str, ...] = ()
    missing_critical_features: tuple[str, ...] = ()
    exclusions_evaluated: tuple[str, ...] = ()
    guideline: dict[str, str] = field(default_factory=dict)
    reason: str = ""

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        for key in ("features_used", "missing_critical_features", "exclusions_evaluated"):
            result[key] = list(result[key])
        return result
