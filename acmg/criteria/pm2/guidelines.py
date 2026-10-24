"""Pinned, source-backed PM2 rules. No user-supplied thresholds are accepted."""

from __future__ import annotations

from dataclasses import dataclass

from .models import PM2Input


@dataclass(frozen=True)
class Guideline:
    name: str
    version: str
    rule_id: str
    url: str
    gene: str | None = None
    disease: str | None = None
    inheritance: str | None = None
    frequency_limit: float | None = None
    minimum_ancestry_allele_number: int | None = None
    requires_male_absence: bool = False

    def citation(self) -> dict[str, str]:
        return {"name": self.name, "version": self.version, "rule_id": self.rule_id, "url": self.url}


SVI = Guideline(
    name="ClinGen SVI PM2 recommendation",
    version="1.0 (2020-09-04)",
    rule_id="PM2_Supporting",
    url="https://www.clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf",
)

ACADVL = Guideline(
    name="ClinGen ACADVL VCEP",
    version="2.2",
    rule_id="PM2_Supporting",
    url="https://erepo.clinicalgenome.org/cspec/ui/svi/doc/639509375",
    gene="ACADVL",
    disease="MONDO:0008723",
    inheritance="AUTOSOMAL_RECESSIVE",
    frequency_limit=0.001,
    minimum_ancestry_allele_number=2001,  # Published rule says >2000 alleles.
)

F9 = Guideline(
    name="ClinGen Coagulation Factor Deficiency VCEP",
    version="2.1",
    rule_id="PM2_Supporting",
    url="https://erepo.clinicalgenome.org/cspec/ui/svi/doc/1527859052",
    gene="F9",
    disease="MONDO:0010604",
    inheritance="X_LINKED",
    requires_male_absence=True,
)

PINNED = (ACADVL, F9)


def select_guideline(features: PM2Input) -> tuple[Guideline, str | None]:
    """Select only an exact supported VCEP scope; never silently use another scope."""
    gene = (features.gene_id or "").upper().removeprefix("HGNC:")
    if features.gene_id == "HGNC:92":
        gene = "ACADVL"
    elif features.gene_id == "HGNC:3551":
        gene = "F9"
    for rule in PINNED:
        if gene == rule.gene:
            if features.disease_id != rule.disease:
                return rule, "disease_id does not match the pinned VCEP scope"
            if features.inheritance_mode is not None and features.inheritance_mode.upper() != rule.inheritance:
                return rule, "inheritance_mode does not match the pinned VCEP scope"
            if features.vcep_specification_version is not None and features.vcep_specification_version != rule.version:
                return rule, "vcep_specification_version does not match the pinned rule"
            return rule, None
    if features.vcep_specification_version is not None:
        return SVI, "requested VCEP specification is not pinned in this engine"
    return SVI, None
