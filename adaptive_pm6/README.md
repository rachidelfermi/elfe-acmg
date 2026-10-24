# Adaptive PM6

Standalone Python 3.12+ evaluator for already-curated ACMG/ClinGen de novo evidence. It does not call variants, infer family relationships, extract phenotypes, download guidelines, or classify a variant overall. The [source-to-rule matrix](SOURCES.md) records each implemented rule.

```python
from adaptive_pm6 import evaluate_pm6

decision = evaluate_pm6({
    "variant_id": "example-variant",
    "gene_id": "example-gene",
    "disease_id": "example-disease",
    "gene_disease_valid": True,
    "proband_variant_present": True,
    "maternal_genotype": "ABSENT",
    "paternal_genotype": "ABSENT",
    "maternal_test_adequate": True,
    "paternal_test_adequate": True,
    "affected": True,
    "phenotype": "HIGHLY_SPECIFIC",
})
print(decision.to_dict())  # PM6, MODERATE, 1.0 de novo point
```

`phenotype` is one of `HIGHLY_SPECIFIC`, `CONSISTENT`, `HETEROGENEOUS`, or `INCONSISTENT`. Use `observations: [...]` for multiple cases, with `case_id` and distinct `family_id`, or curator-confirmed `independent: true`, to establish independence. Shared top-level fields are inherited by each case. A case with both maternity and paternity confirmed is referred to PS2. Unknown parentage remains eligible for assumed de novo scoring once both parental variant tests are reliably negative. Explicitly disproved relationships exclude the case.

For each parent, supply an `ABSENT`/`REFERENCE` genotype or `*_variant_absence: true`, plus `*_test_adequate: true` or an equivalent curated quality conclusion. Raw depth and allele balance have no invented numerical threshold. Explicit zero coverage, poor quality, mosaicism, contamination, mapping ambiguity, sample identity problems, and contradictions prevent scoring. An optional field omitted from a complete pathway does not block activation. Evidence envelopes may use `{"state": "KNOWN", "value": ...}` or `UNKNOWN`, `CONFLICTING`, `NOT_APPLICABLE`.

The built-in VCEP scopes are HNF1A/monogenic diabetes (v3.1, PM6 subsumed by PS2) and SLC9A6/Christianson syndrome (v6.0, gene-specific phenotype points). They are selected by gene and disease identifiers. An explicit unsupported or mismatched `vcep_profile` returns `UNDETERMINED`. The evaluator does not contain the full CSpec registry, so users must check whether another approved VCEP applies before relying on generic SVI scoring. This is the principal remaining scientific limitation.

Run the synthetic rule tests with `python -m pytest tests/test_adaptive_pm6.py -q`.
