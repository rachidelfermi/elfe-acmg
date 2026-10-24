"""Pinned, narrowly scoped ClinGen rules and the ACMG/AMP base rule."""

ACMG = {
    "name": "ACMG/AMP", "version": "2015", "rule_id": "PM1",
    "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/",
}

VCEPS = {
    ("TP53", "MONDO:0018875"): {
        "name": "ClinGen TP53 VCEP", "version": "2.4", "rule_id": "PM1",
        "url": "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN009",
        "transcript": "NM_000546.5", "codons": (175, 245, 248, 249, 273, 282),
    },
    ("DICER1", "MONDO:0100216"): {
        "name": "ClinGen DICER1 VCEP", "version": "1.4", "rule_id": "PM1",
        "url": "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN024",
        "transcript": "NM_177438.2",
        "metal_codons": (1344, 1705, 1709, 1713, 1809, 1810, 1813),
        "domain": (1682, 1846),
    },
    ("F9", "MONDO:0010604"): {
        "name": "ClinGen Coagulation Factor Deficiency VCEP", "version": "2.1",
        "rule_id": "PM1",
        "url": "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN080",
        "transcript": "NM_000133.4",
        "critical_codons": (191, 192, 226, 227, 267, 315, 411),
        "exons": (3, 4, 5),
    },
}
