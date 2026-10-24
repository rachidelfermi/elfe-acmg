"""Supplied evidence only; none of these fields triggers annotation or prediction."""
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Literal


class EvidenceState(StrEnum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"
    NOT_APPLICABLE = "NOT_APPLICABLE"


@dataclass(frozen=True)
class SplicePrediction:
    # A complete SpliceAI score is the maximum of all four delta scores, supplied upstream.
    predictor: str = "SpliceAI"
    score: float | None = None
    no_splice_effect: bool | None = None  # For original ACMG, based on a validated assessment.
    transcript_id_version: str | None = None
    genome_assembly: str | None = None
    variant_match: bool | None = None
    applicable: bool | None = None
    calibrated: bool | None = None
    quality: EvidenceState = EvidenceState.UNKNOWN
    conflict: bool = False
    source: str | None = None
    version: str | None = None


@dataclass(frozen=True)
class RNAEvidence:
    normal_splicing: bool | None = None
    abnormal_splicing: bool | None = None
    validated: bool | None = None
    controls: bool | None = None
    relevant_tissue: bool | None = None
    adequate_sensitivity: bool | None = None
    transcript_id_version: str | None = None
    variant_match: bool | None = None
    conflict: bool = False
    already_used_for_bs3: bool = False
    source: str | None = None


@dataclass(frozen=True)
class BP7Features:
    variant_id: str | None = None
    gene_id: str | None = None
    disease_id: str | None = None
    genome_assembly: str | None = None
    transcript_id_version: str | None = None
    variant_consequence: Literal["synonymous", "intronic", "missense", "other"] | None = None
    consequence_confident: bool | None = None
    biologically_relevant_transcript: bool | None = None
    variant_mapping_valid: bool | None = None
    # First exon base or one of the last three exon bases excludes SVI computational BP7.
    exon_core_splice_motif: bool | None = None
    # Signed intronic coordinate: donor +N or acceptor -N (not genomic strand).
    intronic_position: int | None = None
    canonical_splice_site_overlap: bool | None = None
    regulatory_overlap: bool | None = None
    established_regulatory_splice_effect: bool | None = None
    alternative_transcript_effect: bool | None = None
    known_pathogenic_synonymous_mechanism: bool | None = None
    established_alternative_pathogenic_effect: bool | None = None
    # Original ACMG: upstream curator's reliable categorical assessment; no raw-score cutoff inferred.
    conservation: Literal["not_highly_conserved", "highly_conserved", "unknown", "conflicting"] | None = None
    conservation_reliable: bool | None = None
    conservation_source: str | None = None
    prediction: SplicePrediction | None = None
    rna: RNAEvidence | None = None
    # True means independent application of the same RNA evidence under BS3.
    bp4_splicing_applied: bool = False
    guideline: Literal["ACMG_2015", "SVI_2023", "TP53_2.1.0"] = "SVI_2023"


@dataclass(frozen=True)
class BP7Result:
    criterion: str
    status: Literal["ACTIVATED", "NOT_MET", "UNDETERMINED", "NOT_APPLICABLE"]
    strength: Literal["STRONG", "MODERATE", "SUPPORTING", "NONE"]
    decision_path: str | None
    features_used: tuple[str, ...] = ()
    missing_critical_features: tuple[str, ...] = ()
    exclusions_evaluated: tuple[str, ...] = ()
    guideline: dict[str, str] = field(default_factory=dict)
    reason: str = ""
    overlap_note: str | None = None
