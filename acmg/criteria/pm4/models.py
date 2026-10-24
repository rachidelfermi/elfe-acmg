"""Typed input and output for PM4. None means UNKNOWN for optional evidence."""

from dataclasses import asdict, dataclass, field
from enum import StrEnum
from typing import Any, Generic, TypeVar


class EvidenceState(StrEnum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"
    NOT_APPLICABLE = "NOT_APPLICABLE"


T = TypeVar("T")


@dataclass(frozen=True)
class Evidence(Generic[T]):
    state: EvidenceState
    value: T | None = None
    source: str | None = None

    def __post_init__(self) -> None:
        if (self.state == EvidenceState.KNOWN) != (self.value is not None):
            raise ValueError("KNOWN evidence needs a value; other states cannot have one")


@dataclass(frozen=True)
class BiologicalFeatures:
    # Identity and consequence. Consequence values are INFRAME_DELETION,
    # INFRAME_INSERTION, STOP_LOSS, OTHER, or FRAMESHIFT.
    variant_id: str | None = None
    gene_id: str | None = None
    disease_id: str | None = None
    variant_type: str | None = None
    variant_consequence: Any = None
    genomic_variant_normalization: str | None = None
    coding_hgvs: str | None = None
    protein_hgvs: str | None = None
    transcript_id_version: str | None = None
    protein_accession_version: str | None = None
    transcript_biological_relevance: Any = None
    disease_relevant_transcript: str | None = None
    relevant_transcript_consequences: tuple[Any, ...] = ()
    # Protein consequence. Difference is alternate minus reference length.
    reference_protein_length: int | None = None
    alternate_protein_length: int | None = None
    protein_length_difference: int | None = None
    inserted_amino_acid_count: int | None = None
    deleted_amino_acid_count: int | None = None
    protein_change_start: int | None = None
    protein_change_end: int | None = None
    reading_frame_preserved: Any = None
    protein_extension_length: int | None = None
    native_stop_codon_disrupted: Any = None
    downstream_inframe_stop_codon: Any = None
    predicted_protein_extension: Any = None
    stop_loss_rescue_evidence: Any = None
    # Repeat and benign variation.
    repeat_region_status: Any = None
    repeat_region_type: str | None = None
    repeat_region_boundaries: tuple[int, int] | None = None
    known_benign_repeat_length_variation: Any = None
    regional_indel_tolerance: Any = None
    known_benign_indel_same_region: Any = None
    # Molecular reliability and disease context.
    variant_call_quality: Any = None
    indel_alignment_reliability: Any = None
    protein_consequence_prediction_quality: Any = None
    gene_disease_validity: Any = None
    protein_length_change_mechanism_compatibility: Any = None
    alternative_transcript_rescue: Any = None
    rna_splicing_evidence: Any = None
    protein_product_confirmation: Any = None
    pvs1_applied_to_same_consequence: Any = None
    # Context is retained for audit. It does not independently raise strength.
    critical_protein_domain_overlap: Any = None
    functional_region_overlap: Any = None
    affected_region_functional_evidence: Any = None
    evidence_source: str | None = None
    evidence_version: str | None = None
    evidence_conflicts: tuple[str, ...] = ()
    evidence_missing_reason: str | None = None
    vcep_specification_version: str | None = None


@dataclass(frozen=True)
class PM4Result:
    criterion: str
    status: str
    strength: str
    decision_path: str | None
    features_used: tuple[str, ...]
    missing_critical_features: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    guideline: dict[str, str]
    reason: str

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        for key in ("features_used", "missing_critical_features", "exclusions_evaluated"):
            data[key] = list(data[key])
        return data
