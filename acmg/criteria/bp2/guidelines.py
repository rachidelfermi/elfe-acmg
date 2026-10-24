"""Source-pinned BP2 rules. Add VCEP rules only after primary-source review."""

from dataclasses import dataclass

from .models import GuidelineRef


@dataclass(frozen=True)
class Specification:
    gene: str | None
    disease: str | None
    ref: GuidelineRef
    allowed_phases: tuple[str, ...]
    comparison_lp_allowed: bool = False
    requires_vcep_classification: bool = False
    requires_affected: bool = False
    cis_min_cases: int = 1
    cis_min_unique_comparisons: int = 1
    unknown_phase_min_cases: int | None = None
    trans_min_cases: int = 1
    trans_requires_full_penetrance: bool = True
    trans_requires_phenotype_check: bool = False


ACMG = Specification(
    None, None,
    GuidelineRef("ACMG/AMP and ClinGen Variant Curation SOP", "2015; SOP 3.2 (2022)", "BP2",
                 "https://clinicalgenome.org/site/assets/files/7438/variant_curation_sop_v3_2_oct_2022.pdf"),
    ("CIS", "TRANS"),
)

DICER1 = Specification(
    "DICER1", "MONDO:0100216",
    GuidelineRef("ClinGen DICER1 VCEP", "1.4", "BP2",
                 "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN024"),
    ("CIS", "TRANS", "UNKNOWN"),
    comparison_lp_allowed=True,
    requires_vcep_classification=True,
    cis_min_cases=3,
    cis_min_unique_comparisons=2,
    unknown_phase_min_cases=3,
    trans_requires_full_penetrance=False,
)

BMPR2 = Specification(
    "BMPR2", "MONDO:0015924",
    GuidelineRef("ClinGen Pulmonary Hypertension VCEP", "2.0", "BP2",
                 "https://ldh.clinicalgenome.org/cspec/ui/svi/doc/1568416522"),
    ("CIS",),
    comparison_lp_allowed=True,
    requires_affected=True,
)


def applicable_specification(gene: str, disease: str) -> Specification:
    for spec in (DICER1, BMPR2):
        if gene == spec.gene and disease == spec.disease:
            return spec
    return ACMG
