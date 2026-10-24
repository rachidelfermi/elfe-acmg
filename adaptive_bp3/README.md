# Adaptive BP3 evaluator

`adaptive_bp3` evaluates **only BP3** from supplied, already-extracted protein and
biological evidence. It never annotates a variant or infers a repeat from HGVS or
low-complexity labels. The returned strength is `SUPPORTING` only when BP3 is
activated; otherwise it is `NONE`.

```python
from adaptive_bp3 import BiologicalFeatures, evaluate_bp3

features = BiologicalFeatures(
    gene_id="EXAMPLE_GENE", disease_id="EXAMPLE_DISEASE",
    variant_consequence="INFRAME_DELETION",
    biologically_relevant_transcript=True,
    protein_consequence_confidence=True,
    repeat_region_status=True,
    repeat_region_functional_importance="NO_KNOWN_FUNCTION",
    repeat_length_disease_mechanism=False,
)
print(evaluate_bp3(features).to_dict())
```

The example is synthetic. A caller must substantiate each affirmative assertion.
`NO_KNOWN_FUNCTION` means an affirmative assessment under the relevant guideline;
an absent annotation is `None`, not that value. `Evidence(state=...)` permits
`KNOWN`, `UNKNOWN`, `CONFLICTING`, and `NOT_APPLICABLE` input states. Raw `None`
means unknown. `features_used`, `missing_critical_features`, and
`exclusions_evaluated` explain the decision. No optional evidence field is a
universal requirement.

## Source-to-rule matrix

| Primary source | Implemented BP3 rule |
| --- | --- |
| [Richards et al. 2015, ACMG/AMP](https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/) | In-frame insertion/deletion in a genuine repeat without a known function; Supporting. PM4 and BP3 are not assigned together from the same unresolved repeat evidence. |
| [ClinGen variant classification guidance](https://www.clinicalgenome.org/tools/clingen-variant-classification-guidance/) | No separate SVI BP3 strength or tolerance threshold is listed; none is invented here. |
| [ClinGen Rett/Angelman-like VCEP v2](https://clinicalgenome.org/site/assets/files/7339/clingen_rettas_acmg_specifications_v2.pdf) | For FOXG1 disorder, BP3 applies to wholly contained in-frame indels in poly-His 47–57, poly-Gln 70–73, or poly-Pro 58–61, 65–69, 74–80 on `NP_005240.3`; `NM_005249.4` is the transcript. Supporting. Other repeat-rich FOXG1 regions alone do not qualify. |
| [ClinGen RYR1/MHS VCEP v2](https://www.clinicalgenome.org/site/assets/files/7458/clingen_mhs_acmg_specifications_ryr1_v2.pdf) | BP3 is not applicable for RYR1-related malignant hyperthermia. |
| [ClinGen BMPR2/PAH VCEP v1.0.0](https://cspec.genome.network/cspec/ui/svi/doc/GN125?version=1.0.0) | General BP3 definition retained at Supporting strength for BMPR2-related pulmonary arterial hypertension, on the specified `NM_001204.7` transcript. |
| [ClinGen MYOC/POAG VCEP v1](https://clinicalgenome.org/site/assets/files/7186/clingen_glaucoma_acmg_specifications_v1.pdf) | BP3 is not applicable for MYOC-related open-angle glaucoma. |

The general path needs an explicit disease-specific assessment that a pathogenic
repeat-length mechanism does not apply. Known critical function or pathogenic
repeat-length evidence yields `NOT_MET`; unreliable alignment or conflicting
consequences yields `UNDETERMINED`. A known nonrepeat indel yields `NOT_MET`.
For the curated FOXG1 repeat path, the pinned VCEP definition itself supplies
the repeat and population-tolerance assessment; unrelated population or RNA
fields are not required. A benign repeat variant alone does not create a
separate general BP3 pathway. No universal repeat-copy or length threshold is
used.

## Scope and validation

The pinned VCEP profiles are a documented subset of all ClinGen specifications.
Unrecognized gene–disease pairs use the general ACMG rule. For a recognized
VCEP gene with unresolved disease scope, the result is `UNDETERMINED` rather
than applying the gene-specific rule elsewhere. Coordinates use inclusive
protein positions and must use the pinned accession for FOXG1. The engine
does not validate external evidence sources or adjudicate clinical relevance;
the caller supplies those assessments. Tests use synthetic fixtures, not
published variant-level curations.

Run the BP3 tests with `python -m pytest -q tests/test_bp3_*.py`.
