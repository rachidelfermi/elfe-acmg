"""Typed input and output. None always means unknown, never negative."""

from dataclasses import asdict, dataclass, field
from typing import Any, Literal

State = Literal["KNOWN", "UNKNOWN", "CONFLICTING", "NOT_APPLICABLE"]
Status = Literal["ACTIVATED", "NOT_MET", "UNDETERMINED", "NOT_APPLICABLE"]
Strength = Literal["STRONG", "MODERATE", "SUPPORTING", "NONE"]


@dataclass
class Individual:
    individual_id: str
    affected: bool | None = None
    carries_variant: bool | None = None
    is_proband: bool = False
    informative: bool | None = None  # Curated independent segregation, excluding proband.
    relationship_confirmed: bool | None = None
    genotype_reliable: bool | None = None
    phenotype_consistent: bool | None = None
    biallelic_in_trans: bool | None = None
    unaffected_recessive_at_risk: bool | None = None
    age_at_examination: float | None = None
    expected_onset_age: float | None = None
    sex: str | None = None
    father_id: str | None = None
    mother_id: str | None = None
    genotype_quality: float | None = None
    read_depth: int | None = None
    allele_balance: float | None = None
    copy_number_context: str | None = None
    mosaicism: bool | None = None
    phenotype_specificity: str | None = None
    diagnosis: str | None = None
    phenotype_uncertainty: str | None = None
    evidence_source: str | None = None


@dataclass
class LodEvidence:
    score: float
    method: str
    disease_model_validated: bool | None = None
    penetrance_assumptions_validated: bool | None = None
    pedigree_validated: bool | None = None
    phenotype_model_validated: bool | None = None
    likelihood_ratio: float | None = None
    source: str | None = None


@dataclass
class Family:
    family_id: str
    individuals: list[Individual] = field(default_factory=list)
    affected_segregations: int | None = None
    unaffected_recessive_segregations: int | None = None
    count_validated: bool | None = None
    proband_excluded: bool | None = None
    relationships_confirmed: bool | None = None
    genotype_reliable: bool | None = None
    phenotype_relevant: bool | None = None
    biallelic_phase_confirmed: bool | None = None
    classic_ar_model_validated: bool | None = None
    affected_noncarrier_count: int | None = None
    phenocopy_explained: bool | None = None
    unaffected_carrier_count: int | None = None
    unaffected_carriers_young: bool | None = None
    penetrance_model_validated: bool | None = None
    reduced_penetrance: bool | None = None
    pedigree_conflict: bool | None = None
    sample_identity_conflict: bool | None = None
    genotyping_artifact: bool | None = None
    linked_variant_evidence: bool | None = None
    family_overlap: bool | None = None
    independent: bool | None = None
    study_id: str | None = None
    published_family_duplication: bool | None = None
    lod: LodEvidence | None = None


@dataclass
class PP1Input:
    variant_id: str | None = None
    gene_id: str | None = None
    disease_id: str | None = None
    gene_disease_validity: str | None = None
    inheritance_mode: str | None = None
    variant_reliable: bool | None = None
    disease_mechanism_compatible: bool | None = None
    families: list[Family] = field(default_factory=list)
    guideline_id: str | None = None
    alternative_causal_variant: bool | None = None
    alternative_molecular_diagnosis: bool | None = None
    haplotype_confounding: bool | None = None
    genetic_heterogeneity: bool | None = None
    evidence_overlap_prohibited: bool | None = None
    evidence_overlap_codes: list[str] = field(default_factory=list)
    evidence_states: dict[str, State] = field(default_factory=dict)
    evidence_source: str | None = None
    evidence_version: str | None = None


@dataclass
class PP1Result:
    criterion: str = "PP1"
    status: Status = "UNDETERMINED"
    strength: Strength = "NONE"
    decision_path: str | None = None
    informative_meioses: int | None = None
    lod_score: float | None = None
    families_used: list[str] = field(default_factory=list)
    features_used: list[str] = field(default_factory=list)
    missing_critical_features: list[str] = field(default_factory=list)
    exclusions_evaluated: list[str] = field(default_factory=list)
    guideline: dict[str, str] | None = None
    reason: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
