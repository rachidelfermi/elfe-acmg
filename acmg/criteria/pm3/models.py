"""Typed input records. None means unknown; explicit False means established absence."""

from dataclasses import dataclass, field
from typing import Literal

Classification = Literal["PATHOGENIC", "LIKELY_PATHOGENIC", "VUS", "LIKELY_BENIGN", "BENIGN"]
Phase = Literal["TRANS", "CIS", "UNKNOWN", "CONFLICTING"]
Genotype = Literal["REF", "HET", "HOM_ALT", "HEMIZYGOUS", "OTHER"]


@dataclass
class SecondVariant:
    variant_id: str | None = None
    gene_id: str | None = None
    genotype: Genotype | None = None
    classification: Classification | None = None
    classification_source: str | None = None
    classification_review_status: str | None = None
    classification_evidence: str | None = None
    classification_independent: bool | None = None
    classification_credible: bool | None = None
    vcep_classification_confirmed: bool | None = None
    classification_conflict: bool = False
    disease_match: bool | None = None
    mechanism_compatible: bool | None = None
    rare_for_disease: bool | None = None
    call_valid: bool | None = None
    allele_balance_valid: bool | None = None


@dataclass
class Case:
    case_id: str | None = None
    family_id: str | None = None
    study_ids: tuple[str, ...] = ()
    affected: bool | None = None
    clinical_disease_match: bool | None = None
    phenotype_consistent: bool | None = None
    diagnosis_confirmed: bool | None = None
    alternative_molecular_diagnosis: bool | None = None
    candidate_genotype: Genotype | None = None
    candidate_call_valid: bool | None = None
    candidate_allele_balance_valid: bool | None = None
    second_variant: SecondVariant | None = None
    phase: Phase | None = None
    phase_method: str | None = None
    phase_confident: bool | None = None
    # Parent calls are used only when both variants have informative results.
    maternal_candidate: Genotype | None = None
    paternal_candidate: Genotype | None = None
    maternal_second: Genotype | None = None
    paternal_second: Genotype | None = None
    homozygous_observation: bool | None = None
    homozygous_genotype_valid: bool | None = None
    copy_number_normal: bool | None = None
    deletion_mimicking_homozygosity: bool | None = None
    uniparental_disomy: bool | None = None
    independence_confirmed: bool | None = None
    case_overlap: bool = False
    pp1_overlap: bool = False
    pseudodeficiency_cis_context: bool = False
    functional_support: bool | None = None
    clinical_laboratory_support: bool | None = None
    evidence_conflict: bool = False
    evidence_source: str | None = None
    evidence_method: str | None = None
    evidence_version: str | None = None


@dataclass
class PM3Input:
    candidate_variant_id: str | None = None
    gene_id: str | None = None
    disease_id: str | None = None
    inheritance_mode: str | None = None
    gene_disease_valid: bool | None = None
    disease_mechanism_compatible: bool | None = None
    candidate_rare_for_disease: bool | None = None
    candidate_pm2_supporting: bool | None = None
    cases: list[Case] = field(default_factory=list)
    # Explicit selection of a versioned local specification is required.
    specification: str = "SVI_PM3_V1"
