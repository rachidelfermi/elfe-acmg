# Adaptive BS4

Offline BS4 decisions from pre-curated biological and family evidence. See the [source-to-rule matrix](../docs/bs4-guideline-matrix.md).

```python
from adaptive_bs4 import evaluate_bs4

result = evaluate_bs4({
    "variant_id": "synthetic-v1", "gene_id": "TEST1", "disease_id": "TEST:D1",
    "gene_disease_validity": True, "inheritance_compatible": True,
    "variant_reliable": True,
    "observations": [{
        "family_id": "family-1", "individual_id": "relative-1",
        "affected": True, "variant_present": False,
        "genotype_reliable": True, "phenotype_reliable": True,
        "phenotype_consistent": True,
        "biological_relationship_confirmed": True,
        "phenocopy_excluded": True, "heterogeneity_excluded": True,
    }],
})
print(result.to_dict())  # ACTIVATED, STRONG, synthetic example
```

Each evidence field also accepts `{"state": "KNOWN", "value": true, "source": "..."}`. Other states are `UNKNOWN`, `CONFLICTING`, and `NOT_APPLICABLE`. The input is **already extracted and clinically assessed**; the evaluator does not infer phenotype, genotypes, relationships, phenocopy exclusion, or a likelihood model. Unrelated clinical and variant annotations are optional. Unrecognized VCEP rules are not invented. If a newly published specification applies, review and add a pinned rule before clinical use.
