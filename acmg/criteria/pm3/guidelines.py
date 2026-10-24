"""Source mapped, versioned PM3 specifications.

SVI v1.0, approved 2019-05-02, Tables 1 and 2:
https://www.clinicalgenome.org/site/assets/files/3717/svi_proposal_for_pm3_criterion_-_version_1.pdf
GAA/Pompe VCEP v2, approved 2021-06-02, PM3 section pp. 10-11:
https://www.clinicalgenome.org/site/assets/files/6993/clingen_lsd_acmg_specifications_v2-1.pdf
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Guideline:
    name: str
    version: str
    rule_id: str
    source: str
    gene: str | None = None
    disease: str | None = None
    allow_vus: bool = True
    repeated_unphased_pair_cap: float | None = None


SVI = Guideline(
    "ClinGen SVI PM3", "1.0 (2019-05-02)", "Tables 1-2",
    "https://www.clinicalgenome.org/site/assets/files/3717/svi_proposal_for_pm3_criterion_-_version_1.pdf",
)
GAA = Guideline(
    "ClinGen LSD VCEP GAA PM3", "2 (2021-06-02)", "PM3, pp. 10-11",
    "https://www.clinicalgenome.org/site/assets/files/6993/clingen_lsd_acmg_specifications_v2-1.pdf",
    gene="GAA", disease="Pompe disease", allow_vus=False,
    repeated_unphased_pair_cap=1.0,
)

POINTS = {
    ("TRANS", "PATHOGENIC"): 1.0,
    ("TRANS", "LIKELY_PATHOGENIC"): 1.0,
    ("UNKNOWN", "PATHOGENIC"): 0.5,
    ("UNKNOWN", "LIKELY_PATHOGENIC"): 0.25,
    ("TRANS", "VUS"): 0.25,
}
THRESHOLDS = ((4.0, "VERY_STRONG"), (2.0, "STRONG"), (1.0, "MODERATE"), (0.5, "SUPPORTING"))


def select_guideline(specification: str, gene_id: str | None, disease_id: str | None) -> Guideline:
    if specification == "SVI_PM3_V1":
        return SVI
    if specification == "GAA_POMPE_V2" and gene_id == "GAA" and disease_id == "Pompe disease":
        return GAA
    raise ValueError("Unknown or out-of-scope PM3 specification")
