# Adaptive BP2 evaluator

This package evaluates **BP2 only** from curated observations. It does not phase variants, classify comparison variants, or compute a final ACMG classification. All examples and tests are synthetic rule fixtures, not validated clinical cases.

## Source-to-rule matrix

| Source | BP2 rule implemented | Strength |
| --- | --- | --- |
| [Richards et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/) | Confirmed cis with a pathogenic variant in any inheritance pattern, or confirmed trans with a pathogenic variant in a fully penetrant dominant disorder. | Supporting |
| [ClinGen Variant Curation SOP 3.2, p. 49](https://clinicalgenome.org/site/assets/files/7438/variant_curation_sop_v3_2_oct_2022.pdf) | Comparison pathogenic variant must be observed independently of candidate; do not apply if there is evidence of synergy. | Supporting |
| [DICER1 VCEP 1.4, BP2](https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN024) | One confirmed trans observation with a VCEP classified P/LP DICER1 variant, or at least three cis and/or unphased observations involving at least two different VCEP classified P/LP comparisons. | Supporting |
| [BMPR2 VCEP 2.0, BP2](https://ldh.clinicalgenome.org/cspec/ui/svi/doc/1568416522) | Cis with P/LP BMPR2 variant in an individual with pulmonary arterial hypertension; trans excluded because penetrance is incomplete. | Supporting |

The DICER1 and BMPR2 specifications apply only to their exact gene and MONDO disease IDs. Those genes with a different disease ID return `UNDETERMINED` rather than falling back to the base rule. Other gene/disease contexts use the base ACMG rule. **Before clinical use, check whether a current VCEP specification exists for the exact gene and disease.** This package does not claim exhaustive VCEP coverage or infer scope from a gene symbol alone. It does not accept unverified user-defined strength upgrades.

## Evidence contract

- `Finding` carries `KNOWN`, `UNKNOWN`, `CONFLICTING`, or `NOT_APPLICABLE`, plus optional provenance. A known negative is distinct from unknown.
- The comparison classification requires a source, supporting evidence summary, an independent classification assessment, independent observation of the comparison pathogenic variant, and a disease/mechanism match. A ClinVar assertion alone does not establish these facts.
- Confirmed phase requires a method and a curated reliable-phase finding. Molecular phasing needs no parental genotypes. Unknown phase is accepted only by the DICER1 multi-observation exception and must include both variant genotypes.
- A known complex-allele interaction, haplotype-specific pathogenicity, independent candidate disease effect, or dual diagnosis blocks activation. Cis activation also requires an explicit negative synergy assessment. Unknown ancillary assays and absent call-quality metrics do not block a complete pathway; a known unreliable call does.
- Observation counts use case IDs to avoid counting duplicate publications. The DICER1 count requires case IDs. Family or study IDs alone do not prove independent cases; curators must assign stable case IDs and resolve relatedness.
- Generic BP2 never upgrades above Supporting based on counts. Overlap with PM3, PP1, PS2/PM6, or BP5 is reported but those criteria are not evaluated.

## Example

```python
from dataclasses import asdict
from adaptive_bp2 import BiologicalFeatures, Classification, Comparison, Finding, Inheritance, Observation, Phase, evaluate_bp2

yes, no = Finding.yes(), Finding.no()
comparison = Comparison(
    "ref-variant", gene_id="GENE1", classification=Classification.PATHOGENIC,
    classification_source="curated laboratory", classification_evidence="independent case evidence",
    classification_independent=yes, observed_independently=yes,
    disease_match=yes, mechanism_compatible=yes,
)
observation = Observation(comparison, Phase.CIS, "long-read sequencing", yes, case_id="case-1")
features = BiologicalFeatures(
    "candidate", "GENE1", "MONDO:example", Inheritance.AR, yes, (observation,),
    synergy_or_complex_allele=no,
)
result = evaluate_bp2(features)
print(asdict(result))
```

`result.status` is `ACTIVATED`, `NOT_MET`, `UNDETERMINED`, or `NOT_APPLICABLE` (when the gene-disease relationship is known invalid). A complete result includes the source version, decisive features, missing critical features, and exclusions.

Run tests with `python -m pytest tests/test_bp2.py -q` from the repository root.
