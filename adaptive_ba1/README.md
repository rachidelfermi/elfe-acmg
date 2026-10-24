# Adaptive BA1 evaluator

`adaptive_ba1` evaluates the ACMG/ClinGen **BA1** benign stand-alone criterion from
already-extracted evidence. It does not classify a variant or fetch population data.
All fixtures in `tests/test_ba1.py` are synthetic.

## Source-to-rule matrix

| Rule | Encoded decision | Primary source |
| --- | --- | --- |
| ACMG/AMP 2015 BA1 | Benign stand-alone; original population frequency >5% | [Richards et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/) |
| ClinGen SVI 2018 | AF strictly >0.05 in a general continental population with at least 2,000 observed alleles, unless a gene or variant modification applies | [ClinGen SVI recommendation](https://www.clinicalgenome.org/docs/updated-recommendation-for-the-benign-stand-alone-acmg-amp-criterion/), [ClinGen curation SOP](https://www.clinicalgenome.org/site/assets/files/7438/variant_curation_sop_v3_2_oct_2022.pdf) |
| ClinGen SVI exceptions | Nine Allele Registry IDs excluded from the general rule | [ClinGen July 2018 exception list](https://www.clinicalgenome.org/docs/ba1-exception-list/) |
| ClinGen FH VCEP v1.2 | LDLR / MONDO:0007750: gnomAD PopMax FAF ≥0.005 | [ClinGen specification v1.2](https://www.clinicalgenome.org/site/assets/files/7249/clingen_fh_acmg_specifications_v1_2.pdf) |
| FH structural variants | Same FH thresholds for whole-exon deletions/duplications using gnomAD SV | [FH VCEP publication](https://www.gimjournal.org/article/S1098-3600%2821%2904140-X/fulltext) |

The guideline definitions are versioned in `guidelines.py`. The general rule never
borrows the LDLR threshold. Other VCEP rules are **not bundled**. A known modification
without a bundled specification returns `UNDETERMINED`.

## Example

```python
from adaptive_ba1 import (
    BA1Evidence, EvidenceState, ExceptionCheck, PopulationRecord,
    VariantIdentity, evaluate_ba1,
)

variant = VariantIdentity(
    variant_id="example", assembly="GRCh38", chromosome="1", position=100,
    reference="A", alternate="G", normalized=True,
    mapping_confidence=EvidenceState.KNOWN,
)
population = PopulationRecord(
    database="gnomAD", version="4.1", assembly="GRCh38", chromosome="1",
    position=100, reference="A", alternate="G", normalized=True,
    ancestry="AFR", continental=True, allele_count=600, allele_number=10000,
    filter_pass=True, quality=EvidenceState.KNOWN,
    mapping_ambiguous=False, callability_adequate=True,
)
evidence = BA1Evidence(
    variant=variant, populations=(population,),
    exception_check=ExceptionCheck(
        state=EvidenceState.KNOWN, is_exception=False,
        source="ClinGen SVI exception list", version="July 2018",
    ),
    guideline_modification_check=EvidenceState.KNOWN,
    known_modification=False,
)
result = evaluate_ba1(evidence)
assert result.status == "ACTIVATED"
assert result.strength == "STAND_ALONE"
```

The caller must verify the variant identity, population group classification,
technical quality, applicable guideline modifications, and exception status. A
`KNOWN` exception check needs a source and version. The bundled nine Allele
Registry IDs are automatically excluded from the **general** rule; claiming a
variant is absent from an exception list still requires a documented check.
Unverified mapping, quality, applicability, or exception status yields
`UNDETERMINED`. A reliable frequency below the selected threshold yields
`NOT_MET`. An excluded variant class or verified exception yields
`NOT_APPLICABLE`.

The engine accepts no user-defined thresholds because their clinical validity
cannot be inferred from input alone. Disease prevalence, penetrance, onset,
founder effects, and risk-allele context do not produce a new cutoff. Supplied
special contexts instead flag review. Population homozygote and hemizygote
counts are outside BA1's frequency test; they may be relevant to other criteria.
Mitochondrial, sex chromosome ploidy-sensitive, and most structural variant
contexts need a separately verified specification. The LDLR SV pathway requires
a validated gnomAD SV record and adequate detection sensitivity.

Run: `python -m pytest -q tests/test_ba1.py`.
