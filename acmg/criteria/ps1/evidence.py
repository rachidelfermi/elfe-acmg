"""Evidence-state and provenance primitives used by the PS1 engine.

The evaluator deliberately separates a value from what is known about that
value.  In particular, an omitted value is ``UNKNOWN`` rather than ``FALSE``.
This makes it possible for a rule to distinguish a documented negative finding
from an assertion that has not been assessed.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Callable, Generic, Mapping, Optional, TypeVar


T = TypeVar("T")


class SchemaValidationError(ValueError):
    """Raised when supplied biological evidence cannot safely be interpreted."""


class EvidenceState(str, Enum):
    """Epistemic state for an evidence assertion.

    ``TRUE`` and ``FALSE`` are literal for Boolean features.  For a typed
    categorical feature, ``TRUE`` means that the supplied value is an asserted
    observation; its biological meaning is carried by the value itself.  This
    keeps every feature in the four-state model requested by the API while
    preserving a categorical value such as ``PATHOGENIC``.
    """

    TRUE = "TRUE"
    FALSE = "FALSE"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"


class EvidenceKind(str, Enum):
    OBSERVED = "OBSERVED"
    PREDICTED = "PREDICTED"
    DERIVED = "DERIVED"
    CURATED = "CURATED"


def _token(value: Any) -> str:
    return "_".join(
        piece
        for piece in "".join(
            character if character.isalnum() else " " for character in str(value)
        ).upper().split()
        if piece
    )


def _required_string(value: Any, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise SchemaValidationError(f"{field_name} must be a non-empty string")
    return value.strip()


@dataclass(frozen=True)
class Provenance:
    """Traceability metadata for an input assertion.

    The engine does not impose a universal minimum review status or database
    threshold.  Provenance is retained so a deployment can apply its own data
    governance policy and so that activated decisions remain auditable.
    """

    source: str
    kind: EvidenceKind = EvidenceKind.CURATED
    identifier: Optional[str] = None
    method: Optional[str] = None
    version: Optional[str] = None
    quality: Optional[str] = None
    note: Optional[str] = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "Provenance":
        if not isinstance(raw, Mapping):
            raise SchemaValidationError("provenance entries must be objects")
        source = _required_string(raw.get("source", raw.get("feature_source")), "provenance.source")
        raw_kind = raw.get("kind", "CURATED")
        try:
            kind = EvidenceKind(_token(raw_kind))
        except ValueError as error:
            choices = ", ".join(item.value for item in EvidenceKind)
            raise SchemaValidationError(
                f"provenance.kind must be one of {choices}"
            ) from error
        optional: dict[str, Optional[str]] = {}
        aliases = {
            "method": "feature_method",
            "version": "feature_version",
            "quality": "feature_quality",
        }
        for name in ("identifier", "method", "version", "quality", "note"):
            value = raw.get(name, raw.get(aliases.get(name, "")))
            if value is not None and not isinstance(value, str):
                raise SchemaValidationError(f"provenance.{name} must be a string")
            optional[name] = value.strip() if isinstance(value, str) else value
        return cls(source=source, kind=kind, **optional)

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {"source": self.source, "kind": self.kind.value}
        for name in ("identifier", "method", "version", "quality", "note"):
            value = getattr(self, name)
            if value is not None:
                result[name] = value
        return result


@dataclass(frozen=True)
class Evidence(Generic[T]):
    """A typed observation with an explicit certainty state."""

    state: EvidenceState = EvidenceState.UNKNOWN
    value: Optional[T] = None
    provenance: tuple[Provenance, ...] = ()
    missing_reason: Optional[str] = None

    @classmethod
    def unknown(cls, reason: Optional[str] = None) -> "Evidence[Any]":
        return cls(state=EvidenceState.UNKNOWN, missing_reason=reason)

    @property
    def is_known(self) -> bool:
        return self.state in {EvidenceState.TRUE, EvidenceState.FALSE}

    @property
    def is_unknown(self) -> bool:
        return self.state == EvidenceState.UNKNOWN

    @property
    def is_conflicting(self) -> bool:
        return self.state == EvidenceState.CONFLICTING

    @property
    def is_true(self) -> bool:
        return self.state == EvidenceState.TRUE and self.value is True

    @property
    def is_false(self) -> bool:
        return self.state == EvidenceState.FALSE and self.value is False

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {"state": self.state.value}
        if self.is_known and self.value is not None:
            value: Any = self.value
            if isinstance(value, Enum):
                value = value.value
            result["value"] = value
        if self.provenance:
            result["provenance"] = [entry.to_dict() for entry in self.provenance]
        if self.missing_reason:
            result["missing_reason"] = self.missing_reason
        return result


_MISSING = object()


def parse_bool(raw: Any, field_name: str) -> bool:
    if isinstance(raw, bool):
        return raw
    if isinstance(raw, str):
        token = raw.strip().upper()
        if token in {"TRUE", "YES"}:
            return True
        if token in {"FALSE", "NO"}:
            return False
    raise SchemaValidationError(f"{field_name} must be true or false")


def parse_string(raw: Any, field_name: str) -> str:
    return _required_string(raw, field_name)


def parse_number(raw: Any, field_name: str) -> float:
    if isinstance(raw, bool) or not isinstance(raw, (int, float)):
        raise SchemaValidationError(f"{field_name} must be numeric")
    return float(raw)


def parse_optional_percentage(raw: Any, field_name: str) -> float:
    result = parse_number(raw, field_name)
    if not 0 <= result <= 100:
        raise SchemaValidationError(f"{field_name} must be between 0 and 100")
    return result


def parse_enum(
    enum_type: type[Enum],
    raw: Any,
    field_name: str,
    aliases: Optional[Mapping[str, Enum]] = None,
) -> Enum:
    if isinstance(raw, enum_type):
        return raw
    token = _token(raw)
    if aliases and token in aliases:
        return aliases[token]
    for member in enum_type:
        if token in {_token(member.name), _token(member.value)}:
            return member
    choices = ", ".join(str(member.value) for member in enum_type)
    raise SchemaValidationError(f"{field_name} has unsupported value {raw!r}; expected {choices}")


def parse_evidence(
    raw: Any,
    parser: Callable[[Any, str], T],
    field_name: str,
    *,
    boolean: bool = False,
) -> Evidence[T]:
    """Parse shorthand or an explicit evidence envelope.

    Accepted examples::

        True
        {"state": "FALSE"}
        {"state": "TRUE", "value": "PATHOGENIC", "provenance": {...}}
        {"state": "UNKNOWN", "missing_reason": "not assessed"}

    ``KNOWN`` is accepted as a compatibility spelling for a categorical or
    Boolean assertion with an explicit value.
    """

    if raw is _MISSING or raw is None:
        return Evidence.unknown()

    wrapped = isinstance(raw, Mapping) and any(
        key in raw
        for key in (
            "state",
            "status",
            "value",
            "provenance",
            "feature_source",
            "missing_reason",
            "feature_missing_reason",
        )
    )
    value_raw = raw
    state_raw: Any = None
    provenance_raw: Any = ()
    missing_reason: Optional[str] = None
    if wrapped:
        value_raw = raw.get("value", _MISSING)
        state_raw = raw.get("state", raw.get("status"))
        provenance_raw = raw.get("provenance", ())
        missing_reason = raw.get("missing_reason", raw.get("feature_missing_reason"))
        if missing_reason is not None and not isinstance(missing_reason, str):
            raise SchemaValidationError(f"{field_name}.missing_reason must be a string")
        # Permit the compact feature_source/method/version/quality form.
        if not provenance_raw and raw.get("feature_source") is not None:
            provenance_raw = {
                key: raw[key]
                for key in (
                    "feature_source",
                    "kind",
                    "identifier",
                    "feature_method",
                    "feature_version",
                    "feature_quality",
                    "note",
                )
                if key in raw
            }

    if not wrapped and isinstance(raw, str):
        token = raw.strip().upper()
        if token in {"UNKNOWN", "UNAVAILABLE", "NOT_ASSESSED"}:
            return Evidence(state=EvidenceState.UNKNOWN)
        if token in {"CONFLICTING", "CONFLICT"}:
            return Evidence(state=EvidenceState.CONFLICTING)

    if state_raw is None:
        if boolean and isinstance(raw, bool):
            state = EvidenceState.TRUE if raw else EvidenceState.FALSE
        else:
            state = EvidenceState.TRUE
    else:
        state_token = _token(state_raw)
        if state_token == "KNOWN":
            state = EvidenceState.TRUE
        else:
            try:
                state = EvidenceState(state_token)
            except ValueError as error:
                choices = ", ".join(item.value for item in EvidenceState)
                raise SchemaValidationError(
                    f"{field_name}.state must be one of {choices}"
                ) from error

    if provenance_raw is None:
        provenance_raw = ()
    if isinstance(provenance_raw, Mapping):
        provenance_raw = (provenance_raw,)
    if not isinstance(provenance_raw, (tuple, list)):
        raise SchemaValidationError(f"{field_name}.provenance must be an object or list")
    provenance = tuple(Provenance.from_mapping(item) for item in provenance_raw)

    if state in {EvidenceState.UNKNOWN, EvidenceState.CONFLICTING}:
        return Evidence(state=state, provenance=provenance, missing_reason=missing_reason)

    # ``FALSE`` carries a meaningful value only for Boolean assertions.  A
    # categorical fact such as ``PATHOGENIC`` or ``PVS1`` is represented by
    # its enum value, while absence/conflict uses UNKNOWN/CONFLICTING.  Accepting
    # a FALSE envelope with a positive categorical value would make a
    # contradictory assertion look usable to a decision rule.
    if not boolean and state == EvidenceState.FALSE:
        raise SchemaValidationError(
            f"{field_name} cannot use FALSE state for categorical evidence"
        )

    if boolean and value_raw is _MISSING:
        value_raw = state == EvidenceState.TRUE
    if value_raw is _MISSING or value_raw is None:
        raise SchemaValidationError(f"{field_name} is asserted but has no value")
    value = parser(value_raw, field_name)
    if boolean:
        if state == EvidenceState.TRUE and value is not True:
            raise SchemaValidationError(f"{field_name} has TRUE state but false value")
        if state == EvidenceState.FALSE and value is not False:
            raise SchemaValidationError(f"{field_name} has FALSE state but true value")
    return Evidence(state=state, value=value, provenance=provenance, missing_reason=missing_reason)


def stable_unique(values: list[str] | tuple[str, ...]) -> tuple[str, ...]:
    """Stable de-duplication for reproducible outputs."""

    return tuple(dict.fromkeys(value for value in values if value))

