"""Pinned, directly verified PP2 specifications. No inferred numeric cutoffs."""

from __future__ import annotations

from dataclasses import dataclass

from .models import Guideline


ACMG = Guideline(
    "ACMG/AMP standards and guidelines", "2015", "PP2",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/",
)


@dataclass(frozen=True)
class VcepRule:
    guideline: Guideline
    gene_id: str
    disease_id: str
    inheritance_mode: str
    transcript_id_version: str
    permitted: bool


ACTA1_AD = VcepRule(
    Guideline("ClinGen Congenital Myopathies VCEP", "2.0", "ACTA1_AD_PP2",
              "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN147"),
    "ACTA1", "MONDO:0100084", "AD", "NM_001100.4", True,
)
DICER1_AD = VcepRule(
    Guideline("ClinGen DICER1 VCEP", "1.4", "DICER1_PP2_NOT_APPLICABLE",
              "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN024"),
    "DICER1", "MONDO:0100216", "AD", "NM_177438.2", False,
)

PINNED_RULES = (ACTA1_AD, DICER1_AD)


def find_rule(gene_id: str | None, disease_id: str | None) -> VcepRule | None:
    return next((rule for rule in PINNED_RULES if gene_id == rule.gene_id and disease_id == rule.disease_id), None)
