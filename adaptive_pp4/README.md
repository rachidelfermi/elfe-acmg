# Adaptive PP4 evaluator

This offline package evaluates **PP4 only** from already extracted evidence. It does not infer findings from HPO terms, diagnose a patient, or classify a variant. `KNOWN` means that the caller has reviewed and supplied a finding; absent fields remain `UNKNOWN`.

```python
from adaptive_pp4 import PP4Features, known, evaluate_pp4

features = PP4Features(
    variant_id="example-1", gene_id="GENE", disease_id="MONDO:example",
    gene_disease_validity=known("DEFINITIVE"),
    single_genetic_etiology=known(True),
    phenotype_highly_specific=known(True, source="reviewed clinical record"),
)
result = evaluate_pp4(features)
print(result.status, result.strength, result.guideline.rule_id)
```

`PP4Result` is a frozen dataclass. `dataclasses.asdict(result)` gives a JSON-compatible mapping. All decisions include the rule source, used features, missing critical features, checked exclusions, and explanation. In the general pathway, `phenotype_highly_specific` and `single_genetic_etiology` must be expert-reviewed assertions. A similarity score, a disease rank, a generic HPO match, or a biomarker without a scoped rule cannot create those assertions.

## Source-to-rule matrix

| Source/version | Scope and pathway | Minimum sufficient evidence | Permitted strength | Exclusions and unsupported situations |
|---|---|---|---|---|
| [Richards et al., ACMG/AMP 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/) | General highly specific phenotype or family history | Established gene–disease relationship, reviewed highly specific phenotype **or** family history, and single genetic etiology | Supporting | No activation from similarity or consistency alone; heterogeneous genetic etiology, confirmed alternative cause, or conflicting evidence prevents activation. The evaluator does not infer specificity from raw HPO terms. |
| [ClinGen SVI PP1/PP4 guidance, 2023](https://clinicalgenome.org/docs/clingen-guidance-for-use-of-the-pp1-bs4-co-segregation-and-pp4-phenotype-specificity-criteria-for-sequence-variant-pathogenicity-classification/) | General PP4/PP1 overlap | PP4 family evidence must be distinct from evidence assigned to PP1 | No generic strength escalation implemented | The paper describes a coupled quantitative framework; its disease-specific parameters are not supplied by this engine, so it does not calculate those points. |
| [ClinGen SCID VCEP, ADA v2.2, GN114](https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN114) | ADA, MONDO:0007064 | At least one point from documented SCID diagnosis, sequencing without alternative cause, qualifying ADA biochemical finding, phenotype correction, or lymphocyte profile | Supporting 1–<2; Moderate 2–<6; Strong ≥6 points | Family-history points require proband sequencing without an alternative; treatment correction is capped below Strong when CNV testing is absent and its points are necessary for Strong. ADA activity <1% is the unambiguous subset of the published <1–2% range; values in the boundary range require clinical review. |
| [ClinGen Platelet Disorders VCEP v2.1, GN011](https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN011) | ITGA2B/ITGB3, MONDO:0100326 | Clinical GT diagnosis, mucocutaneous bleeding, normal ristocetin response, minimal/no aggregation with at least two other agonists | Moderate; Strong with qualifying αIIbβ3 expression or fibrinogen-binding confirmation and full sequencing of both genes | Missing required clinical/laboratory evidence blocks activation; known absent requirement means not met. Supporting-only applications are not automated because the rule's generic wording needs clinical interpretation. |
| [ClinGen CDH1 VCEP v3.1](https://www.clinicalgenome.org/site/assets/files/7580/clingen_cdh1_acmg_specifications_v3_1.pdf) | CDH1, MONDO:0007648 | PP4 prohibited; qualifying case evidence is considered under PS4 | None | Returns `NOT_APPLICABLE`. |

## Scientific limits

The implementation has explicit coverage for the general ACMG definition and the three scoped VCEP specifications above. Other VCEP specifications may override the general rule; the caller must check for a current applicable specification before treating a general result as final. This package deliberately does not fetch or update CSpec. The general PP4 pathway relies on expert-reviewed specificity and etiology assertions supplied by the caller. It does not derive those assertions from disease profiles, penetrance, age, sex, diagnostic test performance, or HPO similarity. `age_sex_applicable` can block a known inapplicable presentation; an unknown value does not block an otherwise complete pathway. Component scores are synthetic test inputs, not externally validated patient cases.

The `phenotype_evidence_used_for_ps2_pm6` field is retained for audit. It does not create a universal PP4 exclusion; an applicable VCEP restriction would need its own scoped implementation. Family evidence explicitly marked as reused for PP1 cannot by itself satisfy this engine's general family-history pathway.
