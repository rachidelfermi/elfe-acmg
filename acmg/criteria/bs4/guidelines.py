"""Pinned BS4 rules transcribed from the source-to-rule matrix."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Guideline:
    name: str
    version: str
    rule_id: str
    url: str
    gene_id: str | None = None
    disease_id: str | None = None
    phenotype_score_minimum: int | None = None
    one_event_strength: str = "STRONG"
    multiple_event_strength: str = "STRONG"
    quantitative: bool = False

    def citation(self) -> dict[str, str]:
        return {"name": self.name, "version": self.version,
                "rule_id": self.rule_id, "url": self.url}


ACMG = Guideline("ACMG/AMP", "2015", "BS4",
                 "https://pubmed.ncbi.nlm.nih.gov/25741868/")

PIK3R1 = Guideline("ClinGen Antibody Deficiencies PIK3R1", "1.0", "GN160:BS4",
                   "https://erepo.clinicalgenome.org/cspec/ui/svi/doc/GN160",
                   "PIK3R1", "MONDO:1060136", 6, "SUPPORTING", "STRONG")

CTLA4 = Guideline("ClinGen Antibody Deficiencies CTLA4", "1.0", "1567102241:BS4",
                 "https://cspec.genome.network/cspec/ui/svi/doc/1567102241?version=1.0.0",
                 "CTLA4", "MONDO:0014493", 6, "SUPPORTING", "STRONG")

BRCA2 = Guideline("ClinGen ENIGMA BRCA2", "1.2", "GN097:BS4",
                 "https://cspec.genome.network/cspec/ui/svi/doc/GN097?version=1.2",
                 "BRCA2", "MONDO:0700269", quantitative=True)

SPECIFICATIONS = (PIK3R1, CTLA4, BRCA2)


def select_guideline(gene_id: str | None, disease_id: str | None) -> Guideline:
    for guideline in SPECIFICATIONS:
        if (gene_id or "").upper() == guideline.gene_id and disease_id == guideline.disease_id:
            return guideline
    return ACMG


def brca2_lr_strength(lr: float) -> str | None:
    # GN097 v1.2 also defines Very Strong (<=0.00285), outside this API's strengths.
    if lr <= 0.05:
        return "STRONG"
    if lr <= 0.23:
        return "MODERATE"
    if lr <= 0.48:
        return "SUPPORTING"
    return None
