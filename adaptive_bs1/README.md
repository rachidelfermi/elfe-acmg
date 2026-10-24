# Adaptive BS1 evaluator

`adaptive_bs1` evaluates only ACMG/AMP BS1 from supplied variant and population evidence. It does not annotate variants, download population data, classify a variant overall, or infer a disease threshold from missing assumptions.

## Source-to-rule matrix

| Source | Implemented rule | Scope and conditions |
| --- | --- | --- |
| [Richards et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/) | BS1 is **Strong** when allele frequency is **greater than** expected for the disorder. | General definition only; no universal numerical BS1 cutoff is encoded. BA1 and BS2 are separate criteria. |
| [ClinGen SVI/VCEP guidance on gnomAD v4](https://www.clinicalgenome.org/site/assets/files/9445/clingen_guidance_to_vceps_regarding_the_use_of_gnomad_v4_march_2024.pdf) | gnomAD v4 GroupMax FAF is the 95% lower confidence bound for the group with the highest FAF; bottlenecked groups are excluded. | The evaluator consumes a supplied FAF and does not recalculate it from AC/AN. It requires an explicitly quality-validated record. |
| [ClinGen TPM1 VCEP v1.0, GN100](https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN100) | BS1 Strong at **popmax FAF ≥0.0001**. | TPM1 / hypertrophic cardiomyopathy (MONDO:0005045) / AD; gnomAD FAF. Scope must match exactly. |
| [ClinGen RMRP VCEP v1.3](https://erepo.clinicalgenome.org/cspec/ui/svi/doc/1529230095) | BS1 Strong at **popmax FAF >0.00089**. | RMRP / cartilage-hair hypoplasia (MONDO:0009595) / AR; gnomAD FAF. Exome-derived AF requires ≥20× site depth. |

The TPM1 and RMRP thresholds are published rule values, not values calculated by this package. No generic prevalence/penetrance formula is applied: inheritance, allelic heterogeneity, and disease model must be sourced together. Other VCEP rules may be supplied as `GuidelineRule` with a source, exact scope, operator, metric, and `verified=True` attestation by the caller. That attestation is an input assertion, not independent verification by the package.

## Use

```python
from adaptive_bs1 import BS1Features, EvidenceState, FrequencyRecord, VariantIdentity, evaluate_bs1

variant = VariantIdentity("GRCh38", "15", 63000001, "A", "G")
record = FrequencyRecord(
    variant=variant, database="gnomAD", version="4.1", metric="POPMAX_FAF",
    frequency=0.0002, filter_status="PASS", quality=EvidenceState.KNOWN,
    platform="genome",
)
result = evaluate_bs1(BS1Features(
    variant=variant, gene="TPM1", disease="MONDO:0005045",
    inheritance="AD", frequencies=(record,),
))
assert result.status == "ACTIVATED"
```

`quality=KNOWN` means the caller has validated the record's call quality, callability, mapping and technical limitations. A reported frequency with unknown or failing quality returns `UNDETERMINED`. Exact variant identity includes assembly and representation; normalize indels before use. The evaluator does not treat missing frequency as zero. A reliable frequency below the applicable threshold returns `NOT_MET`; missing or conflicting essential data returns `UNDETERMINED`. BA1 suppresses separate BS1 counting. PM2 at the same time as BS1 raises a conflict. BS2 genotype evidence remains separate.

Founder, risk, and hypomorphic contexts require an explicit exception review (`exception_status=NOT_APPLICABLE`) before activation. A sourced known BS1 exception returns `NOT_APPLICABLE`. Reduced penetrance or late onset alone do not negate a published threshold that already incorporates disease assumptions. Structural, copy-number and mitochondrial variants require class-specific rules and are outside the bundled rules.

The tests use **synthetic** frequency and variant fixtures to exercise the published rule operators and gates; they are not externally curated clinical variant cases. Run `python -m pytest -q tests/test_bs1.py` with Python 3.12+.
