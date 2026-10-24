"""Versioned rules transcribed from the linked primary sources in README."""

from dataclasses import dataclass


@dataclass(frozen=True)
class BA1Guideline:
    name: str
    version: str
    rule_id: str
    gene_id: str | None
    disease_id: str | None
    threshold: float
    operator: str
    metric: str
    database: str | None = None
    minimum_allele_number: int | None = None
    continental_only: bool = False
    allowed_variant_types: tuple[str, ...] = ("SNV", "INDEL")


GENERAL = BA1Guideline(
    "ClinGen SVI BA1 refinement", "2018", "SVI-BA1", None, None,
    0.05, ">", "allele_frequency", minimum_allele_number=2000,
    continental_only=True,
)

LDLR_FH = BA1Guideline(
    "ClinGen FH VCEP", "1.2", "LDLR-FH-BA1", "LDLR", "MONDO:0007750",
    0.005, ">=", "filtering_allele_frequency", database="gnomAD",
    allowed_variant_types=("SNV", "INDEL", "SV", "CNV"),
)

# ClinGen SVI BA1 exception list, July 2018. Allele Registry identifiers avoid
# ambiguous transcript and genome-build representations.
SVI_EXCEPTION_ALLELES = frozenset({
    "CA114709", "CA172210", "CA113797", "CA113795", "CA280114",
    "CA280116", "CA210261", "CA312214", "CA090886",
})


def select_guideline(gene_id: str | None, disease_id: str | None) -> BA1Guideline | None:
    if gene_id == "LDLR":
        return LDLR_FH if disease_id == LDLR_FH.disease_id else None
    return GENERAL
