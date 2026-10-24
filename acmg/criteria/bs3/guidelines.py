"""Source-pinned BS3 strength and one published VCEP restriction."""

from __future__ import annotations

from math import isfinite

SVI = {"name": "ClinGen SVI functional evidence recommendations", "version": "2019", "rule_id": "Brnich-Table-3", "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6938631/"}
SVI_CONTROLS = {**SVI, "rule_id": "Brnich-Step-3-control-count"}
WALKER = {"name": "ClinGen SVI splicing recommendations", "version": "2023", "rule_id": "BP7-RNA-not-BS3", "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/pmid/37352859/"}
MYOC = {"name": "ClinGen Glaucoma VCEP MYOC specification", "version": "2022 publication", "rule_id": "MYOC-BS3-solubility-secretion", "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9771967/"}


def benign_strength_from_oddspath(value: float) -> str:
    """Brnich Table 3: strict comparisons, strong maximum."""
    if not isfinite(value) or value < 0:
        raise ValueError("OddsPath must be finite and nonnegative")
    if value < 0.053:
        return "STRONG"
    if value < 0.23:
        return "MODERATE"
    if value < 0.48:
        return "SUPPORTING"
    return "NONE"


def control_count_strength(benign: int, pathogenic: int) -> str:
    """Brnich Step 3: mixed independent controls; 11 total permit moderate."""
    if any(isinstance(n, bool) or not isinstance(n, int) or n < 0 for n in (benign, pathogenic)):
        raise ValueError("control counts must be nonnegative integers")
    if not benign or not pathogenic:
        return "NONE"
    return "MODERATE" if benign + pathogenic >= 11 else "SUPPORTING"
