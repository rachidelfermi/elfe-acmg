"""Pinned primary-source BP3 rules. No user-supplied strength overrides."""

GENERAL = {"name": "ACMG/AMP", "version": "Richards et al. 2015", "rule_id": "BP3", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/"}
FOXG1 = {"name": "ClinGen Rett/Angelman-like VCEP", "version": "2", "rule_id": "BP3", "url": "https://clinicalgenome.org/site/assets/files/7339/clingen_rettas_acmg_specifications_v2.pdf"}
RYR1 = {"name": "ClinGen Malignant Hyperthermia VCEP", "version": "2", "rule_id": "BP3", "url": "https://www.clinicalgenome.org/site/assets/files/7458/clingen_mhs_acmg_specifications_ryr1_v2.pdf"}
BMPR2 = {"name": "ClinGen Pulmonary Hypertension VCEP", "version": "1.0.0", "rule_id": "BP3", "url": "https://cspec.genome.network/cspec/ui/svi/doc/GN125?version=1.0.0"}
MYOC = {"name": "ClinGen Glaucoma VCEP", "version": "1", "rule_id": "BP3", "url": "https://clinicalgenome.org/site/assets/files/7186/clingen_glaucoma_acmg_specifications_v1.pdf"}

FOXG1_REGIONS = ((47, 57), (70, 73), (58, 61), (65, 69), (74, 80))
FOXG1_PROTEIN = "NP_005240.3"
FOXG1_TRANSCRIPT = "NM_005249.4"
BMPR2_TRANSCRIPT = "NM_001204.7"


def select_guideline(gene: str | None, disease: str | None) -> tuple[str, dict[str, str]]:
    gene = gene.upper() if gene else None
    disease = disease.upper() if disease else None
    if gene == "FOXG1":
        return ("FOXG1", FOXG1) if disease in {"FOXG1 DISORDER", "FOXG1 SYNDROME", "MONDO:0100040"} else ("SCOPE_UNKNOWN", FOXG1)
    if gene == "RYR1":
        return ("RYR1", RYR1) if disease in {"MALIGNANT HYPERTHERMIA", "MONDO:0018493"} else ("SCOPE_UNKNOWN", RYR1)
    if gene == "BMPR2":
        return ("BMPR2", BMPR2) if disease in {"PULMONARY ARTERIAL HYPERTENSION", "MONDO:0015924"} else ("SCOPE_UNKNOWN", BMPR2)
    if gene == "MYOC":
        return ("MYOC", MYOC) if disease in {"PRIMARY OPEN-ANGLE GLAUCOMA", "PRIMARY OPEN ANGLE GLAUCOMA", "JUVENILE ONSET OPEN-ANGLE GLAUCOMA"} else ("SCOPE_UNKNOWN", MYOC)
    return "GENERAL", GENERAL
