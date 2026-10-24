# Adaptive BP5 evaluator

`evaluate_bp5` evaluates **BP5 only** from already-curated case evidence. It does not classify a variant or infer a molecular diagnosis from a disease label, ClinVar assertion, or HPO similarity score. The output is a frozen `BP5Result`; call `.to_dict()` for JSON-ready data.

```python
from adaptive_bp5 import evaluate_bp5

result = evaluate_bp5({
    "candidate_variant_id": "case:variant-1",
    "candidate_gene_id": "GENE1",
    "candidate_disease_id": "MONDO:0000001",
    "candidate_independent_disease_evidence": False,
    "alternative_diagnosis_id": "confirmed-disorder",
    "alternative_variant_id": "case:variant-2",
    "alternative_variant_classification": "PATHOGENIC",
    "alternative_variant_classification_source": "reviewed laboratory report",
    "alternative_variant_classification_evidence": "reviewed primary variant evidence",
    "alternative_variant_classification_independence": True,
    "alternative_gene_disease_validity": True,
    "alternative_inheritance_compatibility": True,
    "patient_clinical_findings": ["relevant clinical finding"],
    "alternative_diagnosis_explains_all_relevant_phenotypes": True,
})
print(result.to_dict())
```

Inputs can be plain values or evidence records such as `{"state": "KNOWN", "value": true, "source": "case report", "version": "2026-01", "quality": "HIGH"}`. `UNKNOWN`, `CONFLICTING`, and `NOT_APPLICABLE` records have no value. Missing fields are unknown. The caller is responsible for supplying reviewed Boolean clinical conclusions and traceable evidence, not computed similarity scores.

## Source-to-rule matrix

| Source and scope | Implemented BP5 rule | Strength |
| --- | --- | --- |
| [Richards et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/4544753/), general ACMG/AMP | Established alternate **molecular** cause that explains the relevant phenotype. An independent candidate contribution, carrier interpretation, blended phenotype, or oligogenic explanation defeats the simple inference. | Supporting |
| [ClinGen SVI guidance index](https://clinicalgenome.org/working-groups/sequence-variant-interpretation/) and [Biesecker et al. guidance](https://clinicalgenome.org/docs/clingen-guidance-for-use-of-the-pp1-bs4-co-segregation-and-pp4-phenotype-specificity-criteria-for-sequence-variant-pathogenicity-classification/) | Phenotype and segregation observations are tracked for reuse across PP4, PP1, BS4, PS2, and PM6. A supplied shared evidence ID blocks BP5 pending review. No universal BP5 strength upgrade is specified here. | None added |
| [ClinGen CSpec GN011 v2.1](https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN011), ITGA2B/ITGB3, Glanzmann thrombasthenia | BP5 explicitly disabled because a patient can carry an unrelated recessive pathogenic allele. | Not applicable |
| [ClinGen CSpec GN092 v1.2](https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN092), BRCA1 cancer predisposition | Clinically calibrated combined LR, with a reviewable description. Co-observation in breast–ovarian cancer risk genes is excluded. LR ≤0.48, ≤0.23, ≤0.05 supports increasing strength. | Supporting, Moderate, Strong |
| [ClinGen CSpec GN127 v1.3](https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN127), RRAS2 RASopathy | Requires supplemental BP5/BP2 case-point scoring. The engine recognizes the scope and returns `UNDETERMINED` rather than substituting general BP5. | None assigned |

The BRCA1 specification also defines Very Strong at LR ≤0.00285. The requested output contract has no `VERY_STRONG` value, so the engine conservatively reports `STRONG` for that range. It does not compute an LR from raw clinical data. Other VCEP rules are not embedded. A caller who requests a VCEP ID or version outside the pinned scope gets `UNDETERMINED`; a gene–disease pair without a pinned rule uses the general ACMG pathway and requires external review for any additional applicable VCEP rule.

## Decision semantics

The general path needs a candidate variant/gene/disease, an independently supported pathogenic or likely pathogenic alternative variant, a valid alternative gene–disease relationship, a genotype assessed compatible with that disease, a relevant patient phenotype, an assessment that the alternative explains **all relevant** findings, and an explicit assessment that the candidate has no independent disease contribution. Biochemical, functional, segregation, age, and family details can inform those supplied conclusions but are not universally required. A clinical or biochemical diagnosis without an alternative molecular variant does not activate general BP5.

`NOT_MET` means supplied reliable evidence defeats the pathway. `UNDETERMINED` means a critical conclusion or provenance is missing or conflicting. `NOT_APPLICABLE` is used only for explicit VCEP exclusions. `ACTIVATED` is criterion-level evidence, not a benign variant classification. All tests use synthetic case fixtures, not published curated patient examples.

Run `python -m pytest tests/test_bp5_engine.py -q` with Python 3.12 or newer.
