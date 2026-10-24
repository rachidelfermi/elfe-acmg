# PS2/PM6 de novo evidence decision table

## Authority and scope

`ps2pm6` evaluates only the shared, mutually exclusive ACMG/AMP de novo
criteria. It accepts already-curated biological evidence; it does not call a
variant, infer a trio relationship, extract a phenotype, query a database, or
combine ACMG/AMP criteria into a classification.

The generic evaluator is pinned to:

- Richards et al. 2015, [ACMG/AMP PS2 and PM6 definitions](https://pubmed.ncbi.nlm.nih.gov/25741868/).
- ClinGen SVI, [Recommendation for de novo Criteria (PS2/PM6), v1.1](https://clinicalgenome.org/site/assets/files/3461/svi_proposal_for_de_novo_criteria_v1_1.pdf), approved 2018-03-18 and updated 2021-05-05.
- ClinGen, [Variant Curation SOP v3.2](https://clinicalgenome.org/site/assets/files/7438/variant_curation_sop_v3_2_oct_2022.pdf), for mutually exclusive use and de-duplication of observations.

The SVI v1.1 clarification matters: **confirmed** and **assumed** refer to
whether parental relationships are confirmed, not to whether a negative
parental genotype is merely guessed. Both parents must be tested for the
variant before an observation can receive a point. Unknown parental genotypes
never become `ABSENT`.

## Generic SVI table implemented

| Phenotype category | Confirmed parental relationships (PS2 points) | Unconfirmed parental relationships (PM6 points) |
| --- | ---: | ---: |
| Highly specific for the gene | 2 | 1 |
| Consistent, but not highly specific | 1 | 0.5 |
| Consistent, not highly specific, and highly genetically heterogeneous | 0.5 | 0.25 |
| Inconsistent with the disorder | 0 | 0 |

The combined contribution from the third row is capped at **1 point**, as
specified in the SVI table. Points are not Tavtigian/Bayesian classification
points.

| Pooled points | Generic strength |
| --- | --- |
| `>= 4` | Very Strong |
| `>= 2` | Strong |
| `>= 1` | Moderate |
| `>= 0.5` | Supporting |
| `< 0.5` | no criterion |

All eligible, independent observations are pooled once. The output always
assigns one code: it prefers `PS2` if any pooled case has confirmed parentage;
otherwise it assigns `PM6`. Thus a mixed confirmed/assumed series is not
reported as both PS2 and PM6, and no patient can be counted twice.

When more than one distinct representative is pooled, each must explicitly
document `case_independence: true`. The engine does not pool known conflicting
gene, disease, or inheritance contexts. Entirely omitted generic context stays
nonblocking, but partially supplied context must agree across pooled cases.

For an autosomal-recessive disorder without an additional pathogenic/likely
pathogenic allele, the final generic strength is reduced by one level. Apparent
germline mosaicism requires confirmed parentage. Credible parental mosaicism,
a parental variant call, an incompatible phenotype, a poor call, or a
non-biological relationship prevents generic activation.

## Input contract and adaptive requirements

Import the evaluator with either supported package spelling:

```python
from ps2pm6 import evaluate_ps2_pm6

result = evaluate_ps2_pm6({
    "proband_variant_present": True,
    "maternal_genotype": "ABSENT",
    "paternal_genotype": "ABSENT",
    "variant_call_quality": "ADEQUATE",
    "maternity_confirmed": True,
    "paternity_confirmed": True,
    "gene_disease_validity": "DEFINITIVE",
    "inheritance_pattern": "AUTOSOMAL_DOMINANT",
    "phenotype_category": "HIGHLY_SPECIFIC",
})

assert result.assigned_criterion == "PS2"
assert result.strength.value == "STRONG"
```

`ps2_pm6`, `adaptive_ps2pm6`, and `adaptive_ps2_pm6` re-export the same
public API for integration compatibility.

Every scalar can be compact (`true`, `"ABSENT"`) or an evidence envelope:

```json
{
  "state": "TRUE",
  "provenance": [{
    "source": "validated trio analysis",
    "kind": "OBSERVED",
    "identifier": "TRIO-42"
  }]
}
```

`UNKNOWN` and `CONFLICTING` remain distinct. Omission is unknown, not a known
negative value. A Boolean `TRUE`/`FALSE` envelope serializes internally as a
known Boolean so it remains compatible with the shared PVS1 evidence model.

For one case, fields can be flat. For multiple cases, use `observations` (or
the `independent_cases` alias); shared top-level fields are inherited only by
cases that omit them:

```json
{
  "gene_id": "HGNC:example",
  "disease_id": "MONDO:example",
  "observations": [
    {"case_id": "A", "case_independence": true, "...": "..."},
    {"case_id": "B", "case_independence": true, "...": "..."}
  ]
}
```

The schema supports the requested genotype, parentage, technical, mosaicism,
disease/phenotype, independent-case, VCEP, and provenance fields. The direct
`phenotype_category` enum is the smallest unambiguous way to choose an SVI
row. Alternatively, `phenotype_consistency`, `phenotype_specificity`, and
`phenotype_heterogeneity` derive the applicable SVI row; when both forms are
supplied, they are cross-checked rather than scored twice.

`independent_case_parentage_status` is a per-observation convenience assertion
(`CONFIRMED` or `ASSUMED`) and supports the same `KNOWN`/`UNKNOWN`/
`CONFLICTING` evidence envelopes as other fields. It never substitutes for
biological parent testing: it only prevents shared parentage context from
silently upgrading the local observation.

Only the active pathway's decision-critical facts are required. For example,
raw parental depth and allele balance are preserved for audit, but the engine
does not invent a universal numerical cut-off from them. An explicit
`LOW`/`INADEQUATE` testing conclusion or zero coverage blocks scoring; a
laboratory should supply `parental_testing_adequacy` or an absence-confidence
assessment when its own validated threshold is decisive.

The generic engine does not infer ClinGen's X-linked unaffected carrier-mother
exception from a child's `maternal_genotype: PRESENT` call or from
`family_history_consistent`. That exception needs a separately modeled de novo
observation in the mother, her own parentage/genotypes, and structured pedigree
facts about male relatives. A maternal-present call therefore remains
non-de-novo in this pathway until a future specialized schema/pathway supports
that evidence explicitly.

## Output contract

`PS2PM6Decision.to_dict()` returns JSON-safe data including:

- `criteria`: `['PS2', 'PM6']`
- `assigned_criterion`: exactly `PS2`, exactly `PM6`, or `null`
- status, strength, pooled total, confirmed/assumed counts
- per-case parentage category, phenotype category, raw and capped point
  contributions, plus `duplicate_of` for an audit-only duplicate report
- minimum missing critical features, exclusions, feature provenance, source,
  version, and deterministic decision path

The only statuses are `ACTIVATED`, `NOT_MET`, `UNDETERMINED`, and
`NOT_APPLICABLE`. Non-activated results always have strength `NONE`.

## VCEP declarations

The package deliberately does not ship a mutable or incomplete VCEP catalogue.
An applicable VCEP document is passed in `applicable_vcep_de_novo_rules` as a
pinned, scoped declaration. It must include a name, version, HTTP(S) source,
gene, disease, and one or more declarative rules. A rule can use a documented
`point_overrides`, `strength_thresholds`, `heterogeneous_points_cap`, final
`assigned_criterion` assertion, or an explicit exclusion.

```json
{
  "name": "Example VCEP PS2/PM6 specification",
  "version": "pinned-version",
  "source_url": "https://example.org/cspec",
  "specification_id": "document-id",
  "gene_id": "HGNC:example",
  "disease_id": "MONDO:example",
  "rules": [{
    "source_locator": "Published PS2/PM6 table",
    "when": {"phenotype_category": "HIGHLY_SPECIFIC"},
    "point_overrides": {
      "HIGHLY_SPECIFIC": {"CONFIRMED": 2, "ASSUMED": 1}
    }
  }]
}
```

An unresolved potentially applicable VCEP scope/condition produces
`UNDETERMINED` rather than bypassing a higher-precedence specification.
Matching VCEP provenance is returned in the output. `assigned_criterion` does
not relabel evidence: the engine derives `PS2` versus `PM6` from qualifying
parentage and rejects a mismatched VCEP assertion. VCEP manifests can restrict
or adjust already eligible de novo observations; they cannot bypass the
generic parental-absence gate to infer the X-linked carrier-mother exception.
The caller remains responsible for transcribing and governing the current
approved CSpec.
