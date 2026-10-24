"""Small audit records for BP6 policy decisions."""

from dataclasses import asdict, dataclass, field
from typing import Literal


@dataclass(frozen=True)
class ExternalAssertion:
    classification: str | None = None
    source: str | None = None
    date: str | None = None
    review_status: str | None = None
    evidence_available: bool | None = None
    evidence_independently_assessable: bool | None = None
    conflicts: bool | None = None
    provenance: str | None = None


@dataclass(frozen=True)
class BP6Result:
    criterion: Literal["BP6"]
    status: Literal["NOT_APPLICABLE", "UNDETERMINED"]
    strength: Literal["NONE"]
    guideline: dict[str, str | None]
    reason: str
    variant: dict[str, str | None] = field(default_factory=dict)
    external_assertions: list[ExternalAssertion] = field(default_factory=list)
    features_used: list[str] = field(default_factory=list)
    missing_critical_features: list[str] = field(default_factory=list)
    recommendations: list[str] = field(default_factory=list)
    applicable_vcep_specification: dict[str, str] | None = None
    policy_conflicts: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Return a JSON-compatible audit record."""
        return asdict(self)
