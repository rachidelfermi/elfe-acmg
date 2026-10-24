"""Pinned, source-specific PP1 rules; no user-supplied thresholds."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Guideline:
    id: str
    name: str
    version: str
    rule_id: str
    source: str
    gene: str | None
    disease_scope: str | None
    modes: tuple[str, ...]
    thresholds: tuple[tuple[int, str], ...]  # descending affected segregations
    min_families_strong: int = 1
    lod_thresholds: tuple[tuple[float, str], ...] = ()
    ar_table: bool = False


GUIDELINES = {
    "SVI_HL_2018": Guideline(
        "SVI_HL_2018", "ClinGen SVI / Hearing Loss VCEP", "1.0 (2018)", "PP1 Table 4a/4b",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC6188673/", None, "HEARING_LOSS",
        ("AD", "AR"), ((5, "STRONG"), (4, "MODERATE"), (2, "SUPPORTING")),
        lod_thresholds=((1.5, "STRONG"), (1.2, "MODERATE"), (0.6, "SUPPORTING")), ar_table=True,
    ),
    "RUNX1_3.0.0": Guideline(
        "RUNX1_3.0.0", "ClinGen RUNX1 Myeloid Malignancy VCEP", "3.0.0", "PP1",
        "https://cspec.genome.network/cspec/ui/svi/doc/GN008?version=3.0.0",
        "RUNX1", "FPD_AML", ("AD",), ((7, "STRONG"), (5, "MODERATE"), (3, "SUPPORTING")),
    ),
    "PTEN_2.0.0": Guideline(
        "PTEN_2.0.0", "ClinGen PTEN VCEP", "2.0.0", "PP1",
        "https://cspec.genome.network/cspec/ui/svi/doc/GN003?version=2.0.0",
        "PTEN", "PHTS", ("AD",), ((7, "STRONG"), (5, "MODERATE"), (3, "SUPPORTING")),
        min_families_strong=2,
    ),
}


def select_guideline(gene: str | None, disease: str | None, requested: str | None) -> Guideline | None:
    gene = gene.upper() if gene else None
    disease = disease.upper() if disease else None
    if requested:
        rule = GUIDELINES.get(requested)
        if rule is None or (rule.gene and rule.gene != gene) or rule.disease_scope != disease:
            return None
        if rule.id == "SVI_HL_2018" and gene not in {"CDH23", "COCH", "GJB2", "KCNQ4", "MYO6", "MYO7A", "SLC26A4", "TECTA", "USH2A"}:
            return None
        return rule
    for rule in GUIDELINES.values():
        if rule.gene and rule.gene == gene and rule.disease_scope == disease:
            return rule
    if disease == "HEARING_LOSS" and gene in {"CDH23", "COCH", "GJB2", "KCNQ4", "MYO6", "MYO7A", "SLC26A4", "TECTA", "USH2A"}:
        return GUIDELINES["SVI_HL_2018"]
    return None
