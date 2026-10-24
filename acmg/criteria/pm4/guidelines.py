"""Pinned PM4 sources and the one explicitly implemented VCEP profile."""

GENERAL = {
    "name": "ACMG/AMP",
    "version": "Richards et al. 2015",
    "rule_id": "PM4",
    "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/",
}

RETT = {
    "name": "ClinGen Rett/Angelman-like Expert Panel",
    "version": "published 2022 profile",
    "rule_id": "PM4",
    "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9135956/",
}

RETT_GENES = frozenset({"MECP2", "UBE3A", "CDKL5", "FOXG1", "TCF4", "SLC9A6"})

# Inclusive protein coordinates from Table 4. A variant's affected interval
# must be supplied in the cited protein reference sequence for automatic overlap.
RETT_REPEAT_REGIONS = {
    "MECP2": ((381, 405),),
    "FOXG1": ((37, 57), (58, 86), (105, 112)),
}
RETT_REPEAT_PROTEINS = {"MECP2": "NP_004983.1", "FOXG1": "NP_005240.3"}
