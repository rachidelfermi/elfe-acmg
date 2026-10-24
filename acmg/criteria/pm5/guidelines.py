"""Pinned primary-source PM5 specifications; no user-supplied rule execution."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Guideline:
    key: str
    name: str
    version: str
    rule_id: str
    url: str
    gene: str | None = None
    disease: str | None = None
    transcript: str | None = None

    def citation(self) -> dict[str, str]:
        return {"name": self.name, "version": self.version, "rule_id": self.rule_id, "url": self.url}


ACMG = Guideline("ACMG_AMP_2015", "ACMG/AMP", "2015", "PM5", "https://pubmed.ncbi.nlm.nih.gov/25741868/")
HNF1A = Guideline("HNF1A_1.1", "ClinGen Monogenic Diabetes VCEP", "1.1", "PM5", "https://www.clinicalgenome.org/site/assets/files/7133/clingen_diabetes_acmg_specifications_v1_1.pdf", "HNF1A", "MONDO:0015967", "NM_000545.8")
PLATELET_ITGA2B = Guideline("PLATELET_1", "ClinGen Platelet Disorders VCEP", "1", "PM5", "https://clinicalgenome.org/docs/clingen-platelet-disorders-expert-panel-specifications-to-the-acmg-amp-variant-interpretation-guidelines-version-1/", "ITGA2B", "MONDO:0010119", "NM_000419.4")
PLATELET_ITGB3 = Guideline("PLATELET_1", "ClinGen Platelet Disorders VCEP", "1", "PM5", "https://clinicalgenome.org/docs/clingen-platelet-disorders-expert-panel-specifications-to-the-acmg-amp-variant-interpretation-guidelines-version-1/", "ITGB3", "MONDO:0010119", "NM_000212.2")


def select_guideline(gene: str | None, disease: str | None, requested: str | None) -> Guideline | None:
    scoped = {"HNF1A": HNF1A, "ITGA2B": PLATELET_ITGA2B, "ITGB3": PLATELET_ITGB3}.get(gene or "")
    if scoped:
        return scoped if disease == scoped.disease and requested in (None, scoped.key) else None
    return ACMG if requested in (None, ACMG.key) else None
