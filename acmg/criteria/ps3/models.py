"""Typed, missing-data-aware models for the PS3-only evaluator."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Generic, Mapping, TypeVar


T = TypeVar("T")


class TextEnum(str, Enum):
    pass


class EvidenceState(TextEnum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class DecisionStatus(TextEnum):
    ACTIVATED = "ACTIVATED"
    NOT_MET = "NOT_MET"
    UNDETERMINED = "UNDETERMINED"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class Strength(TextEnum):
    VERY_STRONG = "VERY_STRONG"
    STRONG = "STRONG"
    MODERATE = "MODERATE"
    SUPPORTING = "SUPPORTING"
    NONE = "NONE"


class FunctionalResult(TextEnum):
    ABNORMAL = "ABNORMAL"
    NORMAL = "NORMAL"
    INTERMEDIATE = "INTERMEDIATE"
    INDETERMINATE = "INDETERMINATE"


class DiseaseMechanism(TextEnum):
    LOSS_OF_FUNCTION = "LOSS_OF_FUNCTION"
    GAIN_OF_FUNCTION = "GAIN_OF_FUNCTION"
    DOMINANT_NEGATIVE = "DOMINANT_NEGATIVE"
    ALTERED_FUNCTION = "ALTERED_FUNCTION"
    OTHER = "OTHER"


class ValidationStatus(TextEnum):
    FORMALLY_CALIBRATED = "FORMALLY_CALIBRATED"
    VALIDATED = "VALIDATED"
    LIMITED = "LIMITED"
    HISTORICALLY_ACCEPTED = "HISTORICALLY_ACCEPTED"
    UNVALIDATED = "UNVALIDATED"


class MechanismRelevance(TextEnum):
    DIRECT = "DIRECT"
    RELEVANT = "RELEVANT"
    INDIRECT = "INDIRECT"
    NOT_RELEVANT = "NOT_RELEVANT"


class SchemaValidationError(ValueError):
    """The supplied evidence cannot be represented without guessing."""


def token(value: Any) -> str:
    return str(value).strip().upper().replace("-", "_").replace(" ", "_")


@dataclass(frozen=True)
class EvidenceValue(Generic[T]):
    state: EvidenceState = EvidenceState.UNKNOWN
    value: T | None = None
    provenance: tuple[Mapping[str, Any], ...] = ()

    @property
    def known(self) -> bool:
        return self.state == EvidenceState.KNOWN

    @classmethod
    def parse(cls, raw: Any, name: str) -> "EvidenceValue[Any]":
        if raw is None:
            return cls()
        if isinstance(raw, EvidenceValue):
            return raw
        value = raw
        state: Any = EvidenceState.KNOWN
        provenance: Any = ()
        if isinstance(raw, Mapping) and any(
            key in raw for key in ("value", "state", "status", "provenance")
        ):
            value = raw.get("value")
            state = raw.get("state", raw.get("status", "KNOWN"))
            provenance = raw.get("provenance", ())
        elif isinstance(raw, str) and token(raw) in EvidenceState.__members__:
            state = token(raw)
            value = None
        try:
            state = state if isinstance(state, EvidenceState) else EvidenceState(token(state))
        except ValueError as exc:
            raise SchemaValidationError(f"{name}.state is invalid: {state!r}") from exc
        if state == EvidenceState.KNOWN and value is None:
            raise SchemaValidationError(f"{name} is KNOWN but has no value")
        if provenance is None:
            provenance = ()
        if isinstance(provenance, Mapping):
            provenance = (provenance,)
        if not isinstance(provenance, (list, tuple)):
            raise SchemaValidationError(f"{name}.provenance must be a list")
        return cls(state=state, value=value if state == EvidenceState.KNOWN else None,
                   provenance=tuple(dict(item) for item in provenance))


CONTEXT_FIELDS = (
    "variant_id", "gene_id", "disease_id", "variant_consequence",
    "disease_mechanism", "gene_disease_validity", "assay_mechanism_relevance",
)

# These fields are retained because they affect applicability, validity,
# reliability, interpretation, strength, exclusions, or auditability. Counts
# that can be derived across studies are intentionally not separate inputs.
ASSAY_FIELDS = (
    "assay_id", "assay_type", "assay_name", "assay_protocol", "assay_system",
    "assay_biological_model", "assay_endpoint", "assay_endpoint_direction",
    "assay_measures_variant_effect", "assay_measures_disease_mechanism",
    "assay_validation_status", "wildtype_control_result",
    "known_pathogenic_control_results", "known_benign_control_results",
    "positive_control_count", "negative_control_count",
    "control_variant_classification_quality", "control_classification_independence",
    "control_result_separation", "assay_dynamic_range", "assay_calibration_method",
    "assay_sensitivity", "assay_specificity", "assay_odds_pathogenicity",
    "assay_validation_sample_size", "assay_validation_method",
    "assay_reproducibility", "inter_assay_variability", "intra_assay_variability",
    "assay_quality_control", "assay_limitations", "variant_assay_result",
    "variant_result_normalized", "variant_effect_direction", "variant_effect_size",
    "variant_result_uncertainty", "variant_result_confidence_interval",
    "variant_result_classification", "variant_result_relative_to_controls",
    "variant_result_within_assay_range", "variant_result_replicate_count",
    "variant_result_replicate_consistency", "technical_replicates",
    "biological_replicates", "experimental_batch_effect", "measurement_error",
    "experimental_blinding", "assay_independent_replication",
    "experimental_artifact_evidence", "variant_expression_level",
    "protein_expression_level", "assay_saturation_or_floor_effect",
    "loss_of_function_effect", "gain_of_function_effect", "dominant_negative_effect",
    "altered_protein_localization", "altered_protein_stability",
    "altered_binding_activity", "altered_enzyme_activity", "altered_channel_function",
    "altered_signaling_activity", "functional_rescue_evidence",
    "mechanism_compatibility", "study_id", "study_independence",
    "evidence_overlap", "evidence_source", "evidence_version", "evidence_quality",
    "evidence_conflicts",
)


@dataclass(frozen=True)
class DiseaseContext:
    values: Mapping[str, EvidenceValue[Any]]

    def get(self, name: str) -> EvidenceValue[Any]:
        return self.values.get(name, EvidenceValue())


@dataclass(frozen=True)
class FunctionalAssay:
    values: Mapping[str, EvidenceValue[Any]]
    index: int = 0

    def get(self, name: str) -> EvidenceValue[Any]:
        return self.values.get(name, EvidenceValue())

    @property
    def identifier(self) -> str:
        for name in ("assay_id", "assay_name", "evidence_source"):
            item = self.get(name)
            if item.known:
                return str(item.value)
        return f"assay-{self.index + 1}"


@dataclass(frozen=True)
class GuidelineReference:
    name: str
    version: str
    rule_id: str
    source_url: str

    def to_dict(self) -> dict[str, str]:
        return {
            "name": self.name,
            "version": self.version,
            "rule_id": self.rule_id,
            "source_url": self.source_url,
        }


@dataclass(frozen=True)
class VCEPRule:
    rule_id: str
    when: Mapping[str, Any]
    status: DecisionStatus
    strength: Strength
    reason: str
    decision_path: str


@dataclass(frozen=True)
class VCEPSpecification:
    name: str
    version: str
    specification_id: str
    source_url: str
    gene_id: str
    disease_id: str | None
    rules: tuple[VCEPRule, ...]
    exclusive: bool = True


@dataclass(frozen=True)
class PS3Evidence:
    context: DiseaseContext
    functional_assays: tuple[FunctionalAssay, ...]
    vcep_specifications: tuple[VCEPSpecification, ...] = ()
    global_values: Mapping[str, EvidenceValue[Any]] = field(default_factory=dict)


@dataclass(frozen=True)
class AssayEvaluation:
    assay_id: str
    status: DecisionStatus
    strength: Strength = Strength.NONE
    decision_path: str = ""
    reason: str = ""
    features_used: tuple[str, ...] = ()
    missing: tuple[str, ...] = ()
    exclusions: tuple[str, ...] = ()
    guideline: GuidelineReference | None = None
    relevance_rank: int = 0
    validation_rank: int = 0
    effect: FunctionalResult | None = None
    bs3_relevant: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "assay_id": self.assay_id,
            "status": self.status.value,
            "strength": self.strength.value,
            "decision_path": self.decision_path,
            "reason": self.reason,
            "features_used": list(self.features_used),
            "missing_critical_features": list(self.missing),
            "exclusions_evaluated": list(self.exclusions),
        }


@dataclass(frozen=True)
class PS3Decision:
    status: DecisionStatus
    strength: Strength
    decision_path: str
    assays_used: tuple[str, ...]
    features_used: tuple[str, ...]
    missing_critical_features: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    conflicting_evidence: tuple[str, ...]
    guideline: GuidelineReference
    reason: str
    assay_evaluations: tuple[AssayEvaluation, ...] = ()
    bs3_relevance: str = "NOT_INDICATED"

    @property
    def features_missing(self) -> tuple[str, ...]:
        """Compatibility alias used by other evidence engines in this workspace."""
        return self.missing_critical_features

    def to_dict(self) -> dict[str, Any]:
        return {
            "criterion": "PS3",
            "status": self.status.value,
            "strength": self.strength.value,
            "decision_path": self.decision_path,
            "assays_used": list(self.assays_used),
            "features_used": list(self.features_used),
            "missing_critical_features": list(self.missing_critical_features),
            "exclusions_evaluated": list(self.exclusions_evaluated),
            "conflicting_evidence": list(self.conflicting_evidence),
            "guideline": self.guideline.to_dict(),
            "reason": self.reason,
            "bs3_relevance": self.bs3_relevance,
            "assay_evaluations": [item.to_dict() for item in self.assay_evaluations],
        }


def unique(items: list[str] | tuple[str, ...]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(item for item in items if item))
