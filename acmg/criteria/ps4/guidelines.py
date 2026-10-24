"""Versioned, source-linked PS4 rules implemented by this package."""

ACMG = {
    "name": "ACMG/AMP standards for sequence variant interpretation",
    "version": "2015",
    "rule_id": "PS4 Note 1",
    "url": "https://pubmed.ncbi.nlm.nih.gov/25741868/",
}

MYOC = {
    "name": "ClinGen Glaucoma VCEP MYOC",
    "version": "2.1",
    "rule_id": "PS4 proband counting",
    "url": "https://actionability.clinicalgenome.org/cspec/ui/svi/doc/635003679?version=2.0.0",
}

MYOC_GENES = frozenset({"MYOC", "HGNC:7610"})
MYOC_DISEASES = frozenset({"MONDO:0005338", "open-angle glaucoma"})


def myoc_scope(gene: str, disease: str) -> bool:
    return gene in MYOC_GENES and disease in MYOC_DISEASES
