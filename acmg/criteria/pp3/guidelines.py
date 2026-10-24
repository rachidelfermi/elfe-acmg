"""Pinned source-to-rule mapping. See docs/pp3-guideline-matrix.md."""

from dataclasses import dataclass


PEJAVER = {"name": "Pejaver et al., ClinGen SVI", "version": "2022", "rule_id": "Table 2 REVEL"}
WALKER = {"name": "Walker et al., ClinGen SVI", "version": "2023", "rule_id": "Table 1 SpliceAI"}
RUNX1 = {"name": "ClinGen Myeloid Malignancy VCEP", "version": "3.1", "rule_id": "GN008 PP3"}


@dataclass(frozen=True)
class VCEP:
    id: str
    version: str
    gene: str
    disease: str
    transcript: str
    missense_min: float
    splice_min: float


VCEPS = {("GN008", "3.1"): VCEP("GN008", "3.1", "RUNX1", "MONDO:0011071", "NM_001754.4", 0.88, 0.38)}


def revel_band(score: float) -> tuple[str, str]:
    if score >= 0.932:
        return "PP3", "STRONG"
    if score >= 0.773:
        return "PP3", "MODERATE"
    if score >= 0.644:
        return "PP3", "SUPPORTING"
    if score <= 0.290:
        return "BP4_CANDIDATE", "NONE"
    return "NEUTRAL", "NONE"


def spliceai_band(score: float) -> str:
    if score >= 0.2:
        return "PP3"
    if score <= 0.1:
        return "BP4_CANDIDATE"
    return "NEUTRAL"
