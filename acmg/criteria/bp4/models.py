"""Typed, already-extracted inputs and deterministic BP4 outputs."""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum


class State(StrEnum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class Status(StrEnum):
    ACTIVATED = "ACTIVATED"
    NOT_MET = "NOT_MET"
    UNDETERMINED = "UNDETERMINED"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class Strength(StrEnum):
    STRONG = "STRONG"
    MODERATE = "MODERATE"
    SUPPORTING = "SUPPORTING"
    NONE = "NONE"


@dataclass(frozen=True)
class Prediction:
    name: str
    score: float | None = None
    version: str | None = None
    variant_match: State = State.UNKNOWN
    transcript_match: State = State.UNKNOWN
    genome_build_match: State = State.UNKNOWN
    model_applicability: State = State.UNKNOWN
    quality: State = State.UNKNOWN
    out_of_distribution: bool = False
    source: str | None = None


@dataclass(frozen=True)
class BiologicalFeatures:
    variant_id: str | None = None
    gene_id: str | None = None
    disease_id: str | None = None
    consequence: str | None = None  # missense, noncanonical_splice, intronic, synonymous, etc.
    transcript_id_version: str | None = None
    biologically_relevant_transcript: State = State.UNKNOWN
    protein_effect_relevant: bool | None = None
    splice_region_eligible: bool | None = None  # True only after site/location review
    predictions: tuple[Prediction, ...] = ()
    pp3_evidence_conflict: bool = False
    bp7_same_prediction_prohibited: bool = False
    splice_evidence_double_counting: bool = False
    conflicting_evidence: bool = False
    vcep_version: str | None = None


@dataclass(frozen=True)
class GuidelineRef:
    name: str
    version: str
    rule_id: str
    url: str


@dataclass(frozen=True)
class BP4Result:
    criterion: str
    status: Status
    strength: Strength
    decision_path: str | None
    predictor_used: str | None
    features_used: tuple[str, ...] = ()
    missing_critical_features: tuple[str, ...] = ()
    exclusions_evaluated: tuple[str, ...] = ()
    guideline: GuidelineRef | None = None
    reason: str = ""
