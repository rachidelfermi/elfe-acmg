"""Typed, already-extracted evidence for BP3. None denotes unknown evidence."""

from dataclasses import asdict, dataclass
from enum import StrEnum
from typing import Generic, TypeVar


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
    version: str | None = None
    quality: str | None = None

    def __post_init__(self) -> None:
        if (self.state == EvidenceState.KNOWN) != (self.value is not None):
            raise ValueError("Only KNOWN evidence has a value")


@dataclass(frozen=True)
class BiologicalFeatures:
    variant_id: str | None = None
    gene_id: str | None = None
    disease_id: str | None = None
    variant_type: str | None = None
    variant_consequence: str | Evidence[str] | None = None
    coding_hgvs: str | None = None
    protein_hgvs: str | None = None
    transcript_id_version: str | None = None
    protein_accession_version: str | None = None
    biologically_relevant_transcript: bool | Evidence[bool] | None = None
    reading_frame_preserved: bool | Evidence[bool] | None = None
    protein_length_change: int | Evidence[int] | None = None
    inserted_amino_acids: str | None = None
    deleted_amino_acids: str | None = None
    affected_protein_positions: tuple[int, int] | Evidence[tuple[int, int]] | None = None
    protein_consequence_confidence: bool | Evidence[bool] | None = None
    repeat_region_status: bool | Evidence[bool] | None = None
    repeat_region_type: str | None = None
    repeat_region_boundaries: tuple[int, int] | None = None
    repeat_unit_sequence: str | None = None
    repeat_unit_length: int | None = None
    repeat_copy_number: int | None = None
    repeat_length_variability: bool | Evidence[bool] | None = None
    tandem_repeat_status: bool | Evidence[bool] | None = None
    low_complexity_region_status: bool | Evidence[bool] | None = None
    repeat_region_mapping_confidence: bool | Evidence[bool] | None = None
    repeat_annotation_source: str | None = None
    # NO_KNOWN_FUNCTION requires an affirmative, sourced assessment; UNKNOWN does not qualify.
    repeat_region_functional_importance: str | Evidence[str] | None = None
    repeat_region_disease_relevance: bool | Evidence[bool] | None = None
    critical_protein_domain_overlap: bool | Evidence[bool] | None = None
    active_site_overlap: bool | Evidence[bool] | None = None
    binding_site_overlap: bool | Evidence[bool] | None = None
    protein_interaction_region_overlap: bool | Evidence[bool] | None = None
    repeat_region_functional_validation: str | Evidence[str] | None = None
    repeat_region_tolerance_evidence: bool | Evidence[bool] | None = None
    known_benign_repeat_length_variants: bool | Evidence[bool] | None = None
    known_pathogenic_repeat_length_variants: bool | Evidence[bool] | None = None
    repeat_length_population_variability: bool | Evidence[bool] | None = None
    repeat_region_conservation: str | None = None
    repeat_region_constraint: str | None = None
    repeat_region_clinical_evidence_quality: str | None = None
    gene_disease_validity: bool | Evidence[bool] | None = None
    disease_inheritance_mode: str | None = None
    disease_molecular_mechanism: str | None = None
    repeat_length_disease_mechanism: bool | Evidence[bool] | None = None
    protein_length_change_mechanism_compatibility: bool | Evidence[bool] | None = None
    gene_specific_repeat_pathogenicity: bool | Evidence[bool] | None = None
    disease_specific_repeat_exceptions: bool | Evidence[bool] | None = None
    indel_alignment_reliability: bool | Evidence[bool] | None = None
    variant_call_quality: bool | Evidence[bool] | None = None
    repeat_mapping_ambiguity: bool | Evidence[bool] | None = None
    protein_consequence_prediction_quality: bool | Evidence[bool] | None = None
    alternative_splice_consequence: bool | Evidence[bool] | None = None
    rna_splicing_evidence: bool | Evidence[bool] | None = None
    alternative_transcript_consequence: bool | Evidence[bool] | None = None
    complex_allele_evidence: bool | Evidence[bool] | None = None
    pm4_applied_to_same_consequence: bool | Evidence[bool] | None = None
    evidence_source: str | None = None
    evidence_version: str | None = None
    evidence_quality: str | None = None
    evidence_conflict_status: bool | Evidence[bool] | None = None
    evidence_missing_reason: str | None = None
    vcep_specification_version: str | None = None


@dataclass(frozen=True)
class BP3Result:
    criterion: str
    status: str
    strength: str
    decision_path: str | None
    features_used: tuple[str, ...]
    missing_critical_features: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    guideline: dict[str, str]
    reason: str

    def to_dict(self) -> dict:
        value = asdict(self)
        for field in ("features_used", "missing_critical_features", "exclusions_evaluated"):
            value[field] = list(value[field])
        return value
