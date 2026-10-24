"""Pinned source rules. No runtime network access or inferred VCEP tables."""

from .models import Phenotype, Strength

SVI_URL = "https://clinicalgenome.org/site/assets/files/3461/svi_proposal_for_de_novo_criteria_v1_1.pdf"
HNF1A_URL = "https://calculator.clinicalgenome.org/cspec/ui/svi/doc/467903868?version=1.1.0"
SLC9A6_URL = "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN033"

POINTS = {
    "ASSUMED": {Phenotype.HIGHLY_SPECIFIC: 1.0, Phenotype.CONSISTENT: 0.5,
                Phenotype.HETEROGENEOUS: 0.25, Phenotype.INCONSISTENT: 0.0},
    "CONFIRMED": {Phenotype.HIGHLY_SPECIFIC: 2.0, Phenotype.CONSISTENT: 1.0,
                  Phenotype.HETEROGENEOUS: 0.5, Phenotype.INCONSISTENT: 0.0},
}
THRESHOLDS = ((4.0, Strength.VERY_STRONG), (2.0, Strength.STRONG),
              (1.0, Strength.MODERATE), (0.5, Strength.SUPPORTING))


def strength_for(points: float) -> Strength:
    return next((strength for threshold, strength in THRESHOLDS if points >= threshold), Strength.NONE)


def select_applicable_guideline(gene_id: str | None, disease_id: str | None,
                                profile: str | None) -> dict[str, str]:
    """The verified HNF1A scope is selected automatically when it matches."""
    hnf1a_scope = gene_id in {"HNF1A", "HGNC:11621"} and disease_id == "MONDO:0015967"
    slc9a6_scope = gene_id in {"SLC9A6", "HGNC:11079"} and disease_id == "MONDO:0010278"
    if hnf1a_scope and profile in {None, "HNF1A-v3.1"}:
        return {"name": "ClinGen Monogenic Diabetes VCEP", "version": "3.1",
                "rule_id": "HNF1A-v3.1-PS2", "source_url": HNF1A_URL}
    if slc9a6_scope and profile in {None, "SLC9A6-v6.0"}:
        return {"name": "ClinGen Rett and Angelman-like VCEP", "version": "6.0",
                "rule_id": "SLC9A6-v6.0-PM6", "source_url": SLC9A6_URL}
    if profile is None or profile == "SVI-v1.1":
        if hnf1a_scope or slc9a6_scope:
            raise ValueError("Applicable VCEP supersedes the general PM6 rule")
        return {"name": "ClinGen SVI de novo", "version": "1.1",
                "rule_id": "SVI-v1.1-tables-1-2", "source_url": SVI_URL}
    raise ValueError("Unsupported or out-of-scope VCEP profile; supply a verified applicable profile")
