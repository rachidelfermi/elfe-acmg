"""Typed input and output models for PVS1 evaluation.

The models intentionally keep *knowledge state* separate from a value.  A
missing Boolean is not a negative Boolean, and a conflicting annotation is not
silently reduced to either.  This is the central safety property of the public
input schema.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
import math
import re
from typing import Any, Callable, Generic, Iterable, Mapping, Optional, Sequence, TypeVar


T = TypeVar("T")


class SchemaValidationError(ValueError):
    """Raised when an input feature set cannot be interpreted safely."""


class TextEnum(str, Enum):
    """String enum with values that serialize directly to JSON."""


class EvidenceStatus(TextEnum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"


class EvidenceKind(TextEnum):
    OBSERVED = "OBSERVED"
    PREDICTED = "PREDICTED"
    DERIVED = "DERIVED"
    CURATED = "CURATED"


class DecisionStatus(TextEnum):
    ACTIVATED = "ACTIVATED"
    NOT_MET = "NOT_MET"
    UNDETERMINED = "UNDETERMINED"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class PVS1Strength(TextEnum):
    VERY_STRONG = "VERY_STRONG"
    STRONG = "STRONG"
    MODERATE = "MODERATE"
    SUPPORTING = "SUPPORTING"
    NONE = "NONE"


class VariantConsequence(TextEnum):
    STOP_GAIN = "STOP_GAIN"
    FRAMESHIFT = "FRAMESHIFT"
    CANONICAL_SPLICE = "CANONICAL_SPLICE"
    SPLICE_SITE = "SPLICE_SITE"
    START_LOST = "START_LOST"
    EXON_DELETION = "EXON_DELETION"
    WHOLE_GENE_DELETION = "WHOLE_GENE_DELETION"
    DUPLICATION = "DUPLICATION"
    OTHER = "OTHER"


class StructuralEvent(TextEnum):
    DELETION = "DELETION"
    DUPLICATION = "DUPLICATION"


class ReadingFrame(TextEnum):
    DISRUPTED = "DISRUPTED"
    PRESERVED = "PRESERVED"


class TandemDuplicationStatus(TextEnum):
    PROVEN_TANDEM = "PROVEN_TANDEM"
    PRESUMED_TANDEM = "PRESUMED_TANDEM"
    PROVEN_NOT_TANDEM = "PROVEN_NOT_TANDEM"


class SpliceOutcome(TextEnum):
    OUT_OF_FRAME = "OUT_OF_FRAME"
    IN_FRAME = "IN_FRAME"
    NO_ABERRANT_SPLICING = "NO_ABERRANT_SPLICING"


class GeneDiseaseValidity(TextEnum):
    DEFINITIVE = "DEFINITIVE"
    STRONG = "STRONG"
    MODERATE = "MODERATE"
    LIMITED = "LIMITED"
    DISPUTED = "DISPUTED"
    REFUTED = "REFUTED"


class LoFMechanismLevel(TextEnum):
    """A curated conclusion equivalent to the SVI 2018 Table 1 evidence tier.

    ``ESTABLISHED`` permits the Figure 1 branch strength.  The two ``CAP_*``
    values retain the published Table 1 one- and two-level caps.  They are not
    feature-count scores.
    """

    ESTABLISHED = "ESTABLISHED"
    CAP_STRONG = "CAP_STRONG"
    CAP_MODERATE = "CAP_MODERATE"
    NOT_ESTABLISHED = "NOT_ESTABLISHED"


class TranscriptSelectionPolicy(TextEnum):
    DISEASE_CURATED = "DISEASE_CURATED"
    MANE_SELECT = "MANE_SELECT"
    EXPLICIT_DOCUMENTED = "EXPLICIT_DOCUMENTED"


def _token(value: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "_", value.upper()).strip("_")


def _enum_value(
    raw: Any,
    enum_type: type[TextEnum],
    field_name: str,
    aliases: Optional[Mapping[str, TextEnum]] = None,
) -> TextEnum:
    if isinstance(raw, enum_type):
        return raw
    if not isinstance(raw, str):
        raise SchemaValidationError(
            f"{field_name} must be a {enum_type.__name__} string, not {raw!r}"
        )
    token = _token(raw)
    if aliases and token in aliases:
        return aliases[token]
    for member in enum_type:
        if token in {_token(member.name), _token(member.value)}:
            return member
    choices = ", ".join(member.value for member in enum_type)
    raise SchemaValidationError(
        f"{field_name} has unsupported value {raw!r}; expected one of {choices}"
    )


def _parse_bool(raw: Any, field_name: str) -> bool:
    if isinstance(raw, bool):
        return raw
    if isinstance(raw, str):
        normalized = raw.strip().lower()
        if normalized in {"true", "yes"}:
            return True
        if normalized in {"false", "no"}:
            return False
    raise SchemaValidationError(f"{field_name} must be true or false, not {raw!r}")


def _parse_nmd(raw: Any, field_name: str) -> bool:
    """Accept the published NMD wording as well as a strict Boolean."""

    if isinstance(raw, str):
        normalized = _token(raw)
        if normalized in {"PREDICTED", "NMD_PREDICTED", "YES"}:
            return True
        if normalized in {
            "NOT_PREDICTED",
            "NMD_NOT_PREDICTED",
            "NOT_PREDICTED_TO_UNDERGO_NMD",
            "NO",
        }:
            return False
    return _parse_bool(raw, field_name)


def _parse_string(raw: Any, field_name: str) -> str:
    if not isinstance(raw, str) or not raw.strip():
        raise SchemaValidationError(f"{field_name} must be a non-empty string")
    return raw.strip()


def _parse_int(raw: Any, field_name: str) -> int:
    if isinstance(raw, bool) or not isinstance(raw, int):
        raise SchemaValidationError(f"{field_name} must be an integer, not {raw!r}")
    return raw


def _parse_splice_site_position(raw: Any, field_name: str) -> int:
    """Parse a signed splice offset without accepting the nonexistent zero site."""

    value = _parse_int(raw, field_name)
    if value == 0:
        raise SchemaValidationError(
            f"{field_name} must be a non-zero signed splice-site offset"
        )
    return value


def _parse_nonnegative_int(raw: Any, field_name: str) -> int:
    value = _parse_int(raw, field_name)
    if value < 0:
        raise SchemaValidationError(f"{field_name} must be greater than or equal to zero")
    return value


def _parse_percent(raw: Any, field_name: str) -> float:
    if isinstance(raw, bool) or not isinstance(raw, (int, float)):
        raise SchemaValidationError(f"{field_name} must be a number from 0 to 100")
    result = float(raw)
    if not math.isfinite(result) or result < 0 or result > 100:
        raise SchemaValidationError(f"{field_name} must be between 0 and 100")
    return result


def _parse_string_tuple(raw: Any, field_name: str) -> tuple[str, ...]:
    if not isinstance(raw, (list, tuple)) or not raw:
        raise SchemaValidationError(f"{field_name} must be a non-empty list")
    result = tuple(_parse_string(value, field_name) for value in raw)
    return result


def _parse_splice_outcomes(raw: Any, field_name: str) -> tuple[SpliceOutcome, ...]:
    values: Sequence[Any]
    if isinstance(raw, (list, tuple)):
        values = raw
    else:
        values = (raw,)
    if not values:
        raise SchemaValidationError(f"{field_name} must contain at least one outcome")
    aliases = {
        "FRAMESHIFT": SpliceOutcome.OUT_OF_FRAME,
        "FRAMESHIFT_NMD": SpliceOutcome.OUT_OF_FRAME,
        "FRAME_DISRUPTING": SpliceOutcome.OUT_OF_FRAME,
        "OUT_OF_FRAME": SpliceOutcome.OUT_OF_FRAME,
        "INFRAME": SpliceOutcome.IN_FRAME,
        "PRESERVES_READING_FRAME": SpliceOutcome.IN_FRAME,
        "NO_EFFECT": SpliceOutcome.NO_ABERRANT_SPLICING,
        "NORMAL_SPLICING": SpliceOutcome.NO_ABERRANT_SPLICING,
    }
    parsed = tuple(
        _enum_value(value, SpliceOutcome, field_name, aliases) for value in values
    )
    # Preserve caller order for provenance but remove exact duplicates.
    return tuple(dict.fromkeys(parsed))


def _json_value(value: Any) -> Any:
    """Return a deterministic JSON-compatible representation for audit output."""

    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    if isinstance(value, set):
        return [_json_value(item) for item in sorted(value, key=repr)]
    return value


@dataclass(frozen=True)
class Provenance:
    """Source information for observed, predicted, derived, or curated evidence."""

    source: str
    kind: EvidenceKind = EvidenceKind.CURATED
    identifier: Optional[str] = None
    method: Optional[str] = None
    note: Optional[str] = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "Provenance":
        if not isinstance(raw, Mapping):
            raise SchemaValidationError("provenance entries must be objects")
        source = _parse_string(raw.get("source"), "provenance.source")
        kind = _enum_value(raw.get("kind", "CURATED"), EvidenceKind, "provenance.kind")
        optional: dict[str, Optional[str]] = {}
        for name in ("identifier", "method", "note"):
            value = raw.get(name)
            if value is not None and not isinstance(value, str):
                raise SchemaValidationError(f"provenance.{name} must be a string")
            optional[name] = value
        return cls(source=source, kind=kind, **optional)

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {"source": self.source, "kind": self.kind.value}
        for name in ("identifier", "method", "note"):
            value = getattr(self, name)
            if value is not None:
                result[name] = value
        return result


@dataclass(frozen=True)
class EvidenceValue(Generic[T]):
    """A value with explicit epistemic status.

    A known negative Boolean has ``status=KNOWN, value=False``.  A missing or
    unavailable value has ``status=UNKNOWN, value=None``.  Conflicting sources
    use ``status=CONFLICTING`` and are deliberately not collapsed.
    """

    status: EvidenceStatus = EvidenceStatus.UNKNOWN
    value: Optional[T] = None
    provenance: tuple[Provenance, ...] = ()

    @classmethod
    def unknown(cls) -> "EvidenceValue[Any]":
        return cls()

    @property
    def is_known(self) -> bool:
        return self.status == EvidenceStatus.KNOWN

    @property
    def is_unknown(self) -> bool:
        return self.status == EvidenceStatus.UNKNOWN

    @property
    def is_conflicting(self) -> bool:
        return self.status == EvidenceStatus.CONFLICTING

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {"state": self.status.value}
        if self.is_known:
            result["value"] = _json_value(self.value)
        if self.provenance:
            result["provenance"] = [item.to_dict() for item in self.provenance]
        return result


def _parse_evidence(
    raw: Any,
    parser: Callable[[Any, str], T],
    field_name: str,
) -> EvidenceValue[T]:
    """Parse a scalar shorthand or the explicit evidence-envelope form."""

    if raw is _MISSING or raw is None:
        return EvidenceValue()

    state_raw: Any = "KNOWN"
    provenance_raw: Any = ()
    value_raw = raw
    wrapped = isinstance(raw, Mapping) and any(
        key in raw for key in ("value", "state", "status", "provenance")
    )
    if wrapped:
        if "state" in raw and "status" in raw:
            stated = _enum_value(raw["state"], EvidenceStatus, f"{field_name}.state")
            status = _enum_value(raw["status"], EvidenceStatus, f"{field_name}.status")
            if stated != status:
                raise SchemaValidationError(
                    f"{field_name}.state and {field_name}.status disagree"
                )
        state_raw = raw.get("state", raw.get("status", "KNOWN"))
        provenance_raw = raw.get("provenance", ())
        value_raw = raw.get("value", _MISSING)

    if isinstance(state_raw, EvidenceStatus):
        state = state_raw
    else:
        state = _enum_value(state_raw, EvidenceStatus, f"{field_name}.state")

    # Shorthand strings make unknown/conflicting ergonomic without compromising
    # fields whose meaningful values are themselves strings.
    if not wrapped and isinstance(raw, str):
        shorthand = raw.strip().lower()
        if shorthand in {"unknown", "unavailable", "not_assessed", "not assessed"}:
            state = EvidenceStatus.UNKNOWN
        elif shorthand in {"conflicting", "conflict"}:
            state = EvidenceStatus.CONFLICTING

    if provenance_raw is None:
        provenance_raw = ()
    if isinstance(provenance_raw, Mapping):
        provenance_raw = (provenance_raw,)
    if not isinstance(provenance_raw, (list, tuple)):
        raise SchemaValidationError(f"{field_name}.provenance must be a list of objects")
    provenance = tuple(Provenance.from_mapping(item) for item in provenance_raw)

    if state != EvidenceStatus.KNOWN:
        return EvidenceValue(status=state, provenance=provenance)
    if value_raw is _MISSING or value_raw is None:
        raise SchemaValidationError(f"{field_name} is known but has no value")
    return EvidenceValue(status=state, value=parser(value_raw, field_name), provenance=provenance)


_MISSING = object()


def _variant_consequence(raw: Any, field_name: str) -> VariantConsequence:
    aliases = {
        "NONSENSE": VariantConsequence.STOP_GAIN,
        "STOP_GAINED": VariantConsequence.STOP_GAIN,
        "STOP_GAIN": VariantConsequence.STOP_GAIN,
        "FRAMESHIFT_VARIANT": VariantConsequence.FRAMESHIFT,
        "CANONICAL_SPLICE": VariantConsequence.CANONICAL_SPLICE,
        "CANONICAL_SPLICE_SITE": VariantConsequence.CANONICAL_SPLICE,
        "CANONICAL_SPLICE_VARIANT": VariantConsequence.CANONICAL_SPLICE,
        "SPLICE_DONOR_VARIANT": VariantConsequence.CANONICAL_SPLICE,
        "SPLICE_ACCEPTOR_VARIANT": VariantConsequence.CANONICAL_SPLICE,
        "INITIATION_CODON": VariantConsequence.START_LOST,
        "INITIATOR_CODON_VARIANT": VariantConsequence.START_LOST,
        "START_LOSS": VariantConsequence.START_LOST,
        "SINGLE_EXON_DELETION": VariantConsequence.EXON_DELETION,
        "MULTI_EXON_DELETION": VariantConsequence.EXON_DELETION,
        "WHOLE_GENE_DEL": VariantConsequence.WHOLE_GENE_DELETION,
        "EXON_DUPLICATION": VariantConsequence.DUPLICATION,
        # Standard VEP/SO consequences that are outside generic PVS1 scope.
        # They remain valid *typed* inputs and produce NOT_APPLICABLE rather
        # than being mistaken for a known biological exclusion of PVS1.
        "MISSENSE_VARIANT": VariantConsequence.OTHER,
        "SYNONYMOUS_VARIANT": VariantConsequence.OTHER,
        "INFRAME_DELETION": VariantConsequence.OTHER,
        "INFRAME_DELETION_VARIANT": VariantConsequence.OTHER,
        "INFRAME_INSERTION": VariantConsequence.OTHER,
        "INFRAME_INSERTION_VARIANT": VariantConsequence.OTHER,
        "DISRUPTIVE_INFRAME_DELETION": VariantConsequence.OTHER,
        "DISRUPTIVE_INFRAME_INSERTION": VariantConsequence.OTHER,
        "CONSERVATIVE_INFRAME_DELETION": VariantConsequence.OTHER,
        "CONSERVATIVE_INFRAME_INSERTION": VariantConsequence.OTHER,
        "PROTEIN_ALTERING_VARIANT": VariantConsequence.OTHER,
        "CODING_SEQUENCE_VARIANT": VariantConsequence.OTHER,
        "INTRON_VARIANT": VariantConsequence.OTHER,
        "INTERGENIC_VARIANT": VariantConsequence.OTHER,
        "UPSTREAM_GENE_VARIANT": VariantConsequence.OTHER,
        "DOWNSTREAM_GENE_VARIANT": VariantConsequence.OTHER,
        "REGULATORY_REGION_VARIANT": VariantConsequence.OTHER,
        "5_PRIME_UTR_VARIANT": VariantConsequence.OTHER,
        "3_PRIME_UTR_VARIANT": VariantConsequence.OTHER,
        "NON_CODING_TRANSCRIPT_VARIANT": VariantConsequence.OTHER,
        "NON_CODING_TRANSCRIPT_EXON_VARIANT": VariantConsequence.OTHER,
        "MATURE_MIRNA_VARIANT": VariantConsequence.OTHER,
        "SPLICE_REGION_VARIANT": VariantConsequence.SPLICE_SITE,
        "STOP_LOST": VariantConsequence.OTHER,
        "STOP_LOSS": VariantConsequence.OTHER,
    }
    return _enum_value(raw, VariantConsequence, field_name, aliases)  # type: ignore[return-value]


def _lof_mechanism(raw: Any, field_name: str) -> LoFMechanismLevel:
    aliases = {
        "FLOWCHART_ELIGIBLE": LoFMechanismLevel.ESTABLISHED,
        "KNOWN": LoFMechanismLevel.ESTABLISHED,
        "TRUE": LoFMechanismLevel.ESTABLISHED,
        "SUPPORTED_STRONG": LoFMechanismLevel.CAP_STRONG,
        "ONE_LEVEL_DOWNGRADE": LoFMechanismLevel.CAP_STRONG,
        "SUPPORTED_MODERATE": LoFMechanismLevel.CAP_MODERATE,
        "TWO_LEVEL_DOWNGRADE": LoFMechanismLevel.CAP_MODERATE,
        "FALSE": LoFMechanismLevel.NOT_ESTABLISHED,
        "ABSENT": LoFMechanismLevel.NOT_ESTABLISHED,
    }
    return _enum_value(raw, LoFMechanismLevel, field_name, aliases)  # type: ignore[return-value]


def _structural_event(raw: Any, field_name: str) -> StructuralEvent:
    aliases = {
        "DEL": StructuralEvent.DELETION,
        "EXON_DELETION": StructuralEvent.DELETION,
        "DUP": StructuralEvent.DUPLICATION,
        "DUPLICATION_VARIANT": StructuralEvent.DUPLICATION,
    }
    return _enum_value(raw, StructuralEvent, field_name, aliases)  # type: ignore[return-value]


def _reading_frame(raw: Any, field_name: str) -> ReadingFrame:
    aliases = {
        "OUT_OF_FRAME": ReadingFrame.DISRUPTED,
        "FRAMESHIFT": ReadingFrame.DISRUPTED,
        "IN_FRAME": ReadingFrame.PRESERVED,
        "INFRAME": ReadingFrame.PRESERVED,
    }
    return _enum_value(raw, ReadingFrame, field_name, aliases)  # type: ignore[return-value]


def _tandem_status(raw: Any, field_name: str) -> TandemDuplicationStatus:
    aliases = {
        "TANDEM": TandemDuplicationStatus.PROVEN_TANDEM,
        "NOT_TANDEM": TandemDuplicationStatus.PROVEN_NOT_TANDEM,
        "PRESUMED": TandemDuplicationStatus.PRESUMED_TANDEM,
    }
    return _enum_value(raw, TandemDuplicationStatus, field_name, aliases)  # type: ignore[return-value]


def _validity(raw: Any, field_name: str) -> GeneDiseaseValidity:
    return _enum_value(raw, GeneDiseaseValidity, field_name)  # type: ignore[return-value]


def _policy(raw: Any, field_name: str) -> TranscriptSelectionPolicy:
    aliases = {
        "CURATED": TranscriptSelectionPolicy.DISEASE_CURATED,
        "DISEASE_RELEVANT": TranscriptSelectionPolicy.DISEASE_CURATED,
        "MANE": TranscriptSelectionPolicy.MANE_SELECT,
        "DOCUMENTED": TranscriptSelectionPolicy.EXPLICIT_DOCUMENTED,
    }
    return _enum_value(raw, TranscriptSelectionPolicy, field_name, aliases)  # type: ignore[return-value]


def _rule_list(raw: Any, field_name: str) -> tuple[Mapping[str, Any], ...]:
    values = raw if isinstance(raw, (list, tuple)) else (raw,)
    if not values:
        raise SchemaValidationError(f"{field_name} cannot be empty when supplied")
    result: list[Mapping[str, Any]] = []
    for value in values:
        if not isinstance(value, Mapping):
            raise SchemaValidationError(f"{field_name} must contain specification objects")
        result.append(dict(value))
    return tuple(result)


@dataclass(frozen=True)
class BiologicalFeatures:
    """Validated biological evidence accepted by :func:`evaluate_pvs1`.

    All members are optional evidence values.  The evaluator only asks for a
    member if the relevant decision-tree branch needs it.  Several supplied
    candidate fields are retained as provenance-bearing inputs even where the
    generic tree intentionally does not derive a decision from them (for
    example, ``distance_to_final_exon_junction`` never auto-creates an NMD
    prediction).
    """

    variant_consequence: EvidenceValue[VariantConsequence] = field(default_factory=EvidenceValue.unknown)
    variant_type: EvidenceValue[str] = field(default_factory=EvidenceValue.unknown)
    coding_position: EvidenceValue[int] = field(default_factory=EvidenceValue.unknown)
    premature_stop_position: EvidenceValue[int] = field(default_factory=EvidenceValue.unknown)
    splice_site_position: EvidenceValue[int] = field(default_factory=EvidenceValue.unknown)
    start_codon_affected: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)

    gene_id: EvidenceValue[str] = field(default_factory=EvidenceValue.unknown)
    disease_id: EvidenceValue[str] = field(default_factory=EvidenceValue.unknown)
    mode_of_inheritance: EvidenceValue[str] = field(default_factory=EvidenceValue.unknown)
    lof_disease_mechanism: EvidenceValue[LoFMechanismLevel] = field(default_factory=EvidenceValue.unknown)
    gene_disease_validity: EvidenceValue[GeneDiseaseValidity] = field(default_factory=EvidenceValue.unknown)
    pathogenic_lof_variant_count: EvidenceValue[int] = field(default_factory=EvidenceValue.unknown)
    pathogenic_lof_variant_fraction_percent: EvidenceValue[float] = field(default_factory=EvidenceValue.unknown)
    pathogenic_lof_variants_across_multiple_exons: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    null_model_recaptures_disease: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    single_exon_gene: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)

    transcript_id: EvidenceValue[str] = field(default_factory=EvidenceValue.unknown)
    transcript_version: EvidenceValue[str] = field(default_factory=EvidenceValue.unknown)
    transcript_selection_policy: EvidenceValue[TranscriptSelectionPolicy] = field(default_factory=EvidenceValue.unknown)
    biologically_relevant_transcript: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    exon_location: EvidenceValue[str] = field(default_factory=EvidenceValue.unknown)
    exon_biological_relevance: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    alternative_transcript_rescue: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    population_evidence_for_exon_lof_tolerance: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)

    nmd_prediction: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    distance_to_final_exon_junction: EvidenceValue[int] = field(default_factory=EvidenceValue.unknown)
    reading_frame_disruption: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    protein_truncation_percent: EvidenceValue[float] = field(default_factory=EvidenceValue.unknown)
    critical_protein_region: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    protein_functional_region_evidence: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)

    splice_outcomes: EvidenceValue[tuple[SpliceOutcome, ...]] = field(default_factory=EvidenceValue.unknown)
    predicted_exon_skipping: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    cryptic_splice_site_usage: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    rna_splicing_evidence: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    rna_assay_interpretable: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    rna_aberrant_splicing_sufficient_for_lof: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    residual_normal_transcript_sufficient_for_function: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    rna_in_frame_loss_of_clinically_relevant_residues: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    splice_rescue_evidence: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)

    cnv_event_type: EvidenceValue[StructuralEvent] = field(default_factory=EvidenceValue.unknown)
    affected_exons: EvidenceValue[tuple[str, ...]] = field(default_factory=EvidenceValue.unknown)
    whole_gene_deletion: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    whole_gene_duplication: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    cnv_reading_frame: EvidenceValue[ReadingFrame] = field(default_factory=EvidenceValue.unknown)
    tandem_duplication_status: EvidenceValue[TandemDuplicationStatus] = field(default_factory=EvidenceValue.unknown)
    duplication_contained_within_gene: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)

    alternative_start_codon: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    upstream_pathogenic_variant_evidence: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    gene_specific_pvs1_rules: EvidenceValue[tuple[Mapping[str, Any], ...]] = field(default_factory=EvidenceValue.unknown)

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "BiologicalFeatures":
        if not isinstance(raw, Mapping):
            raise SchemaValidationError("biological_features must be a mapping")

        aliases = {"splice_effect": "splice_outcomes"}
        normalized: dict[str, Any] = {}
        for name, value in raw.items():
            target = aliases.get(name, name)
            if target in normalized:
                raise SchemaValidationError(
                    f"{name} and {target} both supply the same biological feature"
                )
            normalized[target] = value

        parsers: Mapping[str, Callable[[Any, str], Any]] = {
            "variant_consequence": _variant_consequence,
            "variant_type": _parse_string,
            "coding_position": _parse_int,
            "premature_stop_position": _parse_int,
            "splice_site_position": _parse_splice_site_position,
            "start_codon_affected": _parse_bool,
            "gene_id": _parse_string,
            "disease_id": _parse_string,
            "mode_of_inheritance": _parse_string,
            "lof_disease_mechanism": _lof_mechanism,
            "gene_disease_validity": _validity,
            "pathogenic_lof_variant_count": _parse_nonnegative_int,
            "pathogenic_lof_variant_fraction_percent": _parse_percent,
            "pathogenic_lof_variants_across_multiple_exons": _parse_bool,
            "null_model_recaptures_disease": _parse_bool,
            "single_exon_gene": _parse_bool,
            "transcript_id": _parse_string,
            "transcript_version": _parse_string,
            "transcript_selection_policy": _policy,
            "biologically_relevant_transcript": _parse_bool,
            "exon_location": _parse_string,
            "exon_biological_relevance": _parse_bool,
            "alternative_transcript_rescue": _parse_bool,
            "population_evidence_for_exon_lof_tolerance": _parse_bool,
            "nmd_prediction": _parse_nmd,
            "distance_to_final_exon_junction": _parse_int,
            "reading_frame_disruption": _parse_bool,
            "protein_truncation_percent": _parse_percent,
            "critical_protein_region": _parse_bool,
            "protein_functional_region_evidence": _parse_bool,
            "splice_outcomes": _parse_splice_outcomes,
            "predicted_exon_skipping": _parse_bool,
            "cryptic_splice_site_usage": _parse_bool,
            "rna_splicing_evidence": _parse_bool,
            "rna_assay_interpretable": _parse_bool,
            "rna_aberrant_splicing_sufficient_for_lof": _parse_bool,
            "residual_normal_transcript_sufficient_for_function": _parse_bool,
            "rna_in_frame_loss_of_clinically_relevant_residues": _parse_bool,
            "splice_rescue_evidence": _parse_bool,
            "cnv_event_type": _structural_event,
            "affected_exons": _parse_string_tuple,
            "whole_gene_deletion": _parse_bool,
            "whole_gene_duplication": _parse_bool,
            "cnv_reading_frame": _reading_frame,
            "tandem_duplication_status": _tandem_status,
            "duplication_contained_within_gene": _parse_bool,
            "alternative_start_codon": _parse_bool,
            "upstream_pathogenic_variant_evidence": _parse_bool,
            "gene_specific_pvs1_rules": _rule_list,
        }
        unknown = sorted(set(normalized) - set(parsers))
        if unknown:
            raise SchemaValidationError(
                "unsupported biological feature(s): " + ", ".join(unknown)
            )
        values = {
            name: _parse_evidence(normalized.get(name, _MISSING), parser, name)
            for name, parser in parsers.items()
        }
        return cls(**values)

    def evidence(self, field_name: str) -> EvidenceValue[Any]:
        try:
            return getattr(self, field_name)
        except AttributeError as error:
            raise SchemaValidationError(f"unknown feature {field_name!r}") from error

    def normalized(self) -> "BiologicalFeatures":
        """Revalidate a direct dataclass input through the public schema.

        Public callers may construct BiologicalFeatures manually. Python
        dataclasses do not enforce annotations, so the entry point must not
        trust those values more than mapping input.
        """

        raw: dict[str, Any] = {}
        for name in self.__dataclass_fields__:
            value = getattr(self, name)
            if not isinstance(value, EvidenceValue):
                raise SchemaValidationError(
                    f"{name} must be an EvidenceValue when BiologicalFeatures is constructed directly"
                )
            if not isinstance(value.status, EvidenceStatus):
                raise SchemaValidationError(
                    f"{name}.status must be an EvidenceStatus when BiologicalFeatures is constructed directly"
                )
            # Do not delegate this validation to ``EvidenceValue.to_dict``:
            # manually constructed EvidenceValue instances can contain an
            # invalid provenance member, and that method would otherwise leak
            # an AttributeError instead of the schema-level error promised by
            # this public entry point.
            provenance: list[dict[str, Any]] = []
            if not isinstance(value.provenance, tuple):
                raise SchemaValidationError(
                    f"{name}.provenance must be a tuple when BiologicalFeatures is constructed directly"
                )
            for index, item in enumerate(value.provenance):
                if not isinstance(item, Provenance):
                    raise SchemaValidationError(
                        f"{name}.provenance[{index}] must be a Provenance when BiologicalFeatures is constructed directly"
                    )
                try:
                    reparsed = Provenance.from_mapping(
                        {
                            "source": item.source,
                            "kind": item.kind,
                            "identifier": item.identifier,
                            "method": item.method,
                            "note": item.note,
                        }
                    )
                except (SchemaValidationError, AttributeError, TypeError) as error:
                    raise SchemaValidationError(
                        f"{name}.provenance[{index}] is invalid when BiologicalFeatures is constructed directly"
                    ) from error
                provenance.append(reparsed.to_dict())
            envelope: dict[str, Any] = {"state": value.status.value}
            if value.status == EvidenceStatus.KNOWN:
                envelope["value"] = value.value
            if provenance:
                envelope["provenance"] = provenance
            raw[name] = envelope
        return BiologicalFeatures.from_mapping(raw)


@dataclass(frozen=True)
class GuidelineReference:
    """A pinned document/reference, not an invented branch identifier."""

    authority: str
    document: str
    version: Optional[str]
    source_url: str
    source_locator: str
    rule_id: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "authority": self.authority,
            "name": self.document,
            "version": self.version,
            "source_url": self.source_url,
            "source_locator": self.source_locator,
            "rule_id": self.rule_id,
        }


@dataclass(frozen=True)
class PVS1Decision:
    """Auditable result of evaluating exactly one PVS1 criterion."""

    status: DecisionStatus
    strength: PVS1Strength
    guideline: GuidelineReference
    decision_path: str
    reason: str
    features_used: tuple[str, ...] = ()
    features_missing: tuple[str, ...] = ()
    exclusions_evaluated: tuple[Mapping[str, str], ...] = ()
    evidence_provenance: tuple[Provenance, ...] = ()
    evidence_used: Mapping[str, Mapping[str, Any]] = field(default_factory=dict)
    guideline_conditions: Mapping[str, Any] = field(default_factory=dict)
    evidence_mode: Optional[str] = None

    def __post_init__(self) -> None:
        if self.status == DecisionStatus.ACTIVATED and self.strength == PVS1Strength.NONE:
            raise ValueError("an activated PVS1 decision requires a non-NONE strength")
        if self.status != DecisionStatus.ACTIVATED and self.strength != PVS1Strength.NONE:
            raise ValueError("a non-activated PVS1 decision must have NONE strength")

    @property
    def criterion(self) -> str:
        return "PVS1"

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "criterion": self.criterion,
            "status": self.status.value,
            "strength": self.strength.value,
            "guideline": self.guideline.to_dict(),
            "decision_path": self.decision_path,
            "features_used": list(self.features_used),
            "features_missing": list(self.features_missing),
            "exclusions_evaluated": [dict(item) for item in self.exclusions_evaluated],
            "evidence_provenance": [item.to_dict() for item in self.evidence_provenance],
            "evidence_used": {
                key: _json_value(value) for key, value in self.evidence_used.items()
            },
            "guideline_conditions": _json_value(self.guideline_conditions),
            "reason": self.reason,
        }
        if self.evidence_mode is not None:
            result["evidence_mode"] = self.evidence_mode
        return result


def stable_unique(values: Iterable[str]) -> tuple[str, ...]:
    """Stable deduplication used to make decisions reproducible."""

    return tuple(dict.fromkeys(value for value in values if value))
