"""Structured, JSON-friendly PS1 decisions and internal rule outcomes."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional

from .evidence import Provenance
from .models import DecisionStatus, PS1Strength


@dataclass(frozen=True)
class GuidelineReference:
    """Pinned source metadata for an evaluated rule."""

    authority: str
    name: str
    version: str
    source_url: str
    source_locator: str
    rule_id: str = "PS1"

    def to_dict(self) -> dict[str, str]:
        return {
            "authority": self.authority,
            "name": self.name,
            "version": self.version,
            "source_url": self.source_url,
            "source_locator": self.source_locator,
            "rule_id": self.rule_id,
        }


@dataclass(frozen=True)
class RuleEvaluation:
    """One reference/pathway result before cross-reference resolution."""

    status: DecisionStatus
    strength: PS1Strength
    decision_path: str
    guideline: GuidelineReference
    reason: str
    reference_index: Optional[int] = None
    features_used: tuple[str, ...] = ()
    missing_critical_features: tuple[str, ...] = ()
    exclusions_evaluated: tuple[Mapping[str, str], ...] = ()

    def __post_init__(self) -> None:
        if self.status == DecisionStatus.ACTIVATED and self.strength == PS1Strength.NONE:
            raise ValueError("an activated PS1 rule must have a non-NONE strength")
        if self.status != DecisionStatus.ACTIVATED and self.strength != PS1Strength.NONE:
            raise ValueError("a non-activated PS1 rule must have NONE strength")


@dataclass(frozen=True)
class PS1Decision:
    """Final auditable decision for exactly one PS1 assessment."""

    status: DecisionStatus
    strength: PS1Strength
    decision_path: str
    guideline: GuidelineReference
    reason: str
    reference_variants_used: tuple[Mapping[str, str], ...] = ()
    features_used: tuple[str, ...] = ()
    missing_critical_features: tuple[str, ...] = ()
    exclusions_evaluated: tuple[Mapping[str, str], ...] = ()
    evidence_provenance: tuple[Provenance, ...] = ()
    reference_evaluations: tuple[Mapping[str, Any], ...] = ()

    def __post_init__(self) -> None:
        if self.status == DecisionStatus.ACTIVATED and self.strength == PS1Strength.NONE:
            raise ValueError("an activated PS1 decision requires a non-NONE strength")
        if self.status != DecisionStatus.ACTIVATED and self.strength != PS1Strength.NONE:
            raise ValueError("a non-activated PS1 decision must have NONE strength")

    @property
    def criterion(self) -> str:
        return "PS1"

    def to_dict(self) -> dict[str, Any]:
        # `features_missing` is included as a compatibility alias for callers
        # using the PVS1-style audit vocabulary; the canonical PS1 field is
        # `missing_critical_features`.
        return {
            "criterion": self.criterion,
            "status": self.status.value,
            "strength": self.strength.value,
            "decision_path": self.decision_path,
            "reference_variants_used": [dict(item) for item in self.reference_variants_used],
            "features_used": list(self.features_used),
            "missing_critical_features": list(self.missing_critical_features),
            "features_missing": list(self.missing_critical_features),
            "exclusions_evaluated": [dict(item) for item in self.exclusions_evaluated],
            "guideline": self.guideline.to_dict(),
            "evidence_provenance": [entry.to_dict() for entry in self.evidence_provenance],
            "reference_evaluations": [dict(item) for item in self.reference_evaluations],
            "reason": self.reason,
        }

