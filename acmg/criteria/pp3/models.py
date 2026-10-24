"""Typed, already-extracted evidence and deterministic PP3 result."""

from dataclasses import dataclass, field
from enum import Enum


class EvidenceState(str, Enum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class Decision(str, Enum):
    ACTIVATED = "ACTIVATED"
    NOT_MET = "NOT_MET"
    UNDETERMINED = "UNDETERMINED"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class Strength(str, Enum):
    VERY_STRONG = "VERY_STRONG"
    STRONG = "STRONG"
    MODERATE = "MODERATE"
    SUPPORTING = "SUPPORTING"
    NONE = "NONE"


@dataclass(frozen=True)
class Prediction:
    name: str
    score: float | None = None
    version: str | None = None
    calibration: str | None = None
    state: EvidenceState = EvidenceState.KNOWN
    source: str | None = None
    variant_id: str | None = None
    transcript_id_version: str | None = None
    genome_build: str | None = None
    quality_ok: bool | None = None
    applicable: bool | None = None
    out_of_distribution: bool | None = None
    interpretation: str | None = None


@dataclass(frozen=True)
class Features:
    variant_consequence: str | None = None
    gene_id: str | None = None
    disease_id: str | None = None
    variant_id: str | None = None
    variant_type: str | None = None
    coding_hgvs: str | None = None
    protein_hgvs: str | None = None
    transcript_id_version: str | None = None
    genome_build: str | None = None
    biologically_relevant_transcript: bool | None = None
    gene_disease_validity: bool | None = None
    mechanism_compatible: bool | None = None
    splice_context_established: bool | None = None
    canonical_splice_site: bool | None = None
    intron_number: int | None = None
    pvs1_splice_overlap: bool | None = None
    rna_splice_evidence_used: bool | None = None
    bp4_evidence_conflict: bool | None = None
    conflicting_computational_evidence: bool | None = None
    selected_missense_predictor: str | None = None
    predictions: tuple[Prediction, ...] = ()
    vcep_id: str | None = None
    vcep_version: str | None = None
    missing_states: dict[str, EvidenceState] = field(default_factory=dict)
    conservation: dict | None = None
    structural_predictions: dict | None = None


@dataclass(frozen=True)
class Result:
    status: Decision
    strength: Strength
    decision_path: str | None
    predictor_used: str | None
    features_used: tuple[str, ...]
    missing_critical_features: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    guideline: dict[str, str] | None
    reason: str
    bp4_candidate: bool = False

    @property
    def decision(self) -> Decision:
        return self.status

    def to_dict(self) -> dict:
        return {
            "criterion": "PP3", "status": self.status.value,
            "strength": self.strength.value, "decision_path": self.decision_path,
            "predictor_used": self.predictor_used,
            "features_used": list(self.features_used),
            "missing_critical_features": list(self.missing_critical_features),
            "exclusions_evaluated": list(self.exclusions_evaluated),
            "guideline": self.guideline, "reason": self.reason,
            "bp4_candidate": self.bp4_candidate,
        }
