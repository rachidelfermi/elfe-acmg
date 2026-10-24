"""Typed, already-curated biological inputs and deterministic result."""

from dataclasses import asdict, dataclass, field
from typing import Literal

EvidenceState = Literal["KNOWN", "UNKNOWN", "CONFLICTING", "NOT_APPLICABLE"]
Classification = Literal["PATHOGENIC", "LIKELY_PATHOGENIC", "VUS", "LIKELY_BENIGN", "BENIGN"]


@dataclass(frozen=True)
class Candidate:
    variant_id: str | None = None
    consequence: str | None = None
    transcript_id_version: str | None = None
    protein_accession_version: str | None = None
    amino_acid_position: int | None = None
    reference_amino_acid: str | None = None
    alternate_amino_acid: str | None = None
    transcript_relevant: bool | None = None
    splice_mechanism: Literal["PROTEIN", "SPLICING", "UNKNOWN", "CONFLICTING"] | None = None
    grantham_distance: int | None = None  # Curated distance from reference to alternate amino acid.


@dataclass(frozen=True)
class Reference(Candidate):
    classification: Classification | None = None
    classification_source: str | None = None
    classification_evidence_reviewed: bool | None = None
    classification_independent_of_pm5: bool | None = None
    classification_conflicting: bool | None = None
    disease_match: bool | None = None
    mechanism_compatible: bool | None = None
    classification_date: str | None = None


@dataclass(frozen=True)
class BiologicalFeatures:
    candidate: Candidate
    references: tuple[Reference, ...] = ()
    gene_id: str | None = None
    disease_id: str | None = None
    gene_disease_validity: Literal["ESTABLISHED", "DISPUTED", "REFUTED", "UNKNOWN"] | None = None
    candidate_mechanism_compatible: bool | None = None
    known_benign_substitution_at_residue: bool | None = None
    biological_evidence_conflicting: bool | None = None
    specification: str | None = None  # ACMG_AMP_2015, HNF1A_1.1, PLATELET_1
    evidence_source: str | None = None
    evidence_version: str | None = None
    evidence_states: dict[str, EvidenceState] = field(default_factory=dict)


@dataclass(frozen=True)
class PM5Result:
    criterion: str
    status: Literal["ACTIVATED", "NOT_MET", "UNDETERMINED", "NOT_APPLICABLE"]
    strength: Literal["VERY_STRONG", "STRONG", "MODERATE", "SUPPORTING", "NONE"]
    decision_path: str
    reference_variants_used: tuple[str, ...] = ()
    features_used: tuple[str, ...] = ()
    missing_critical_features: tuple[str, ...] = ()
    exclusions_evaluated: tuple[str, ...] = ()
    guideline: dict[str, str] = field(default_factory=dict)
    reason: str = ""

    def to_dict(self) -> dict:
        return asdict(self)
