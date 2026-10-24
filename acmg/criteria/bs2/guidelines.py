"""Pinned published BS2 rules. No user-defined thresholds are executed."""

ACMG_2015 = {
    "name": "ACMG/AMP",
    "version": "2015",
    "rule_id": "BS2",
    "url": "https://pmc.ncbi.nlm.nih.gov/articles/pmid/25741868/",
}

RETT_ANGELMAN = {
    "name": "ClinGen Rett/Angelman-like VCEP",
    "version": "2022 publication",
    "rule_id": "BS2",
    "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9135956/",
}

RUNX1 = {
    "name": "ClinGen RUNX1 MM-VCEP",
    "version": "GN008 v3.0.0",
    "rule_id": "BS2 not applicable",
    "url": "https://cspec.genome.network/cspec/ui/svi/doc/GN008?version=3.0.0",
}

# The panel permits related observations for these counts. UBE3A requires
# maternal inheritance when related, and maternal allele is required here
# for every UBE3A observation to avoid inferring imprinting phase.
RETT_GENOTYPES = {
    "FOXG1": ("HETEROZYGOUS",),
    "TCF4": ("HETEROZYGOUS",),
    "SLC9A6": ("HEMIZYGOUS",),
    "CDKL5": ("HETEROZYGOUS", "HEMIZYGOUS"),
    "MECP2": ("HETEROZYGOUS", "HEMIZYGOUS"),
    "UBE3A": ("HETEROZYGOUS",),
}
