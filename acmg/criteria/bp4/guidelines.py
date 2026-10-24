"""Versioned rules transcribed from cited primary sources; no runtime downloads."""
from __future__ import annotations

from dataclasses import dataclass
from .models import GuidelineRef, Strength


PEJAVER = GuidelineRef("ClinGen SVI missense calibration (Pejaver et al.)", "2022", "Table 2: REVEL", "https://pubmed.ncbi.nlm.nih.gov/36413997/")
WALKER = GuidelineRef("ClinGen SVI splicing (Walker et al.)", "2023", "SpliceAI noncanonical BP4", "https://pubmed.ncbi.nlm.nih.gov/37352859/")
RPE65 = GuidelineRef("ClinGen LCA/eoRD VCEP RPE65", "1.0", "BP4", "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/1566365955?version=1.0.0")


@dataclass(frozen=True)
class Rule:
    name: str
    consequence: str
    bands: tuple[tuple[float, Strength], ...]  # ascending inclusive upper limits
    guideline: GuidelineRef
    allowed_versions: tuple[str, ...] = ()  # only when a specification pins model versions

    def strength(self, score: float) -> Strength:
        for upper, strength in self.bands:
            if score <= upper:
                return strength
        return Strength.NONE


REVEL = Rule("REVEL", "missense", ((0.003, Strength.STRONG), (0.016, Strength.STRONG), (0.183, Strength.MODERATE), (0.290, Strength.SUPPORTING)), PEJAVER)
# Pejaver identifies <=0.003 as very strong. This API has no VERY_STRONG level;
# the conservative STRONG cap is explicit, never an invented upgrade.
SPLICEAI = Rule("SpliceAI", "noncanonical_splice", ((0.1, Strength.MODERATE),), WALKER)
BAYESDEL = Rule("BayesDel_noAF", "missense", ((-0.36, Strength.MODERATE), (-0.18, Strength.SUPPORTING)), GuidelineRef(PEJAVER.name, PEJAVER.version, "Table 2: BayesDel without allele frequency", PEJAVER.url))

RPE65_GENE = {"RPE65", "HGNC:10294"}
RPE65_DISEASE = {"MONDO:0100368", "RPE65-related recessive retinopathy"}
RPE65_TRANSCRIPT = "NM_000329.3"
PALB2 = GuidelineRef("ClinGen HBOPC VCEP PALB2", "1.3", "BP4", "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN077")
PALB2_GENE = {"PALB2", "HGNC:26144"}
PALB2_DISEASE = {"MONDO:0700272", "MONDO:0012565"}
PALB2_TRANSCRIPT = "NM_024675.3"
