"""Pinned, published BP7 rule definitions. See docs/bp7-guideline-matrix.md."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Guideline:
    name: str
    version: str
    rule_id: str
    source: str
    requires_conservation: bool
    permits_intronic: bool
    spliceai_cutoff: float | None
    permits_rna_strong: bool


GUIDELINES = {
    "ACMG_2015": Guideline("ACMG/AMP", "2015", "BP7", "https://pubmed.ncbi.nlm.nih.gov/25741868/", True, False, None, False),
    "SVI_2023": Guideline("ClinGen SVI Splicing", "2023", "BP7 / BP7_Strong (RNA)", "https://pubmed.ncbi.nlm.nih.gov/37352859/", False, True, 0.1, True),
    "TP53_2.1.0": Guideline("ClinGen TP53 VCEP", "2.1.0", "BP7 / BP7_Strong (RNA)", "https://cspec.genome.network/cspec/ui/svi/doc/GN009?version=2.1.0", False, True, 0.1, True),
}
