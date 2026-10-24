# BP6 evaluator

`adaptive_bp6` is an offline, standalone policy evaluator for BP6. It accepts
already extracted external classification evidence and returns an auditable
decision. It does not retrieve classifications or score any other ACMG criterion.

```python
from adaptive_bp6 import evaluate_bp6

result = evaluate_bp6({
    "variant_id": "example",
    "external_classification": "BENIGN",
    "classification_source": "ClinVar",
    "classification_review_status": "criteria provided, single submitter",
})
assert result.status == "NOT_APPLICABLE"
assert result.strength == "NONE"
print(result.to_dict())
```

The default `bp6_usage_policy` is `CURRENT_CLINGEN`. Under that policy, BP6
is `NOT_APPLICABLE` with `NONE` strength even when a ClinVar, expert panel,
or laboratory benign assertion is present. Review stars are retained as
provenance, never converted into evidence strength. Missing classifications,
review status, underlying evidence, and other ACMG criteria are allowed.
External assertions may be supplied as flat `external_classification` and
`classification_*` fields or as an `external_assertions` list of mappings or
`ExternalAssertion` records. Variant identity fields are optional and retained
for audit. `features_used` names policy inputs, not biological evidence.

`HISTORICAL_ACMG_2015` selects the original BP6 rule for audit. It returns
`UNDETERMINED` and assigns no strength because source reputation, recency,
evidence availability, and independence cannot be established from a label or
review status. If a verified, matching VCEP explicitly excludes BP6, the
historical result is `NOT_APPLICABLE`. Thus this evaluator does not generate
historical BP6 supporting evidence.

## Verified policy sources

- [Richards et al. (2015), Table 4](https://pubmed.ncbi.nlm.nih.gov/25741868/): BP6 originally meant a reputable source recently reported a variant benign while the evidence was unavailable to the laboratory for independent evaluation. It was classified as benign supporting evidence. “2015” is a publication year, not a numbered guideline version.
- [Biesecker and Harrison, ClinGen SVI (2018)](https://pubmed.ncbi.nlm.nih.gov/29543229/): recommends discontinuing PP5 and BP6. The [ClinGen guidance index](https://clinicalgenome.org/tools/clingen-variant-classification-guidance/) continues to list this recommendation. The recommendation has no numbered version in this evaluator; `guideline.version` is `null`.
- [ClinGen Coagulation Factor Deficiency VCEP, F9, CSpec GN080 version 2.1](https://cspec.genome.network/cspec/ui/svi/doc/GN080?version=2.1): explicitly lists BP6 as Not Applicable for F9 (HGNC:3551) and hemophilia B (MONDO:0010604). The evaluator applies this VCEP record only with matching ID, version, gene and disease. Unknown or mismatched specifications are reported as unverified and cannot authorize BP6.

An external assertion may prompt examination of independently assessable
primary evidence under appropriate benign criteria. Its classification is
never counted as independent biological evidence or added to the strength of
criteria using the same underlying observations. Because the evaluator does
not inspect that primary evidence, it cannot assign those criteria or perform
a complete variant classification.
