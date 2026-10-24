# Unified ACMG/ClinGen validation audit — 2026-10-23

## Scope and conclusion

The repository contains all 28 original criterion packages, 26 active registry
entries (PS2 and PM6 share one de novo evaluator), and disabled PP5/BP6. The
unified classifier accepts pre-extracted criterion-specific evidence. This is a
**code and synthetic-rule validation**, not clinical validation. A full
criterion-by-criterion comparison with expert-curated variant dossiers remains
open; the repository does not contain those source evidence dossiers.

## Criterion coverage matrix

`Tested` means synthetic evaluator tests exist; it does not mean clinical
concordance. Every listed active evaluator distinguishes an activated result
from missing/undetermined input in its own API. `Scoped` means the documented
generic or pinned paths only; it does not imply every current VCEP rule.

| Criterion | Actual status | Guideline / pathway scope | Tests | Main boundary |
|---|---|---|---|---|
| PVS1 | Active, scoped | ACMG 2015; SVI 2018 LoF tree; SVI 2023 RNA; caller-pinned VCEP | `test_pvs1_engine` | No universal VCEP catalogue |
| PS1 | Active, scoped | ACMG same amino acid; curated splice comparison; pinned VCEP | `test_adaptive_ps1_engine` | Reference classification independence required |
| PS2 | Active via PS2_PM6 | SVI de novo v1.1, confirmed parentage | `test_ps2pm6_engine` | One pooled de novo code |
| PS3 | Active, scoped | ACMG functional assays; SVI assay calibration | `test_ps3`, `test_assay_validation` | Assay validity must be supplied |
| PS4 | Active, scoped | ACMG case-control/case evidence; selected VCEP | `test_ps4` | No universal case threshold |
| PM1 | Active, scoped | ACMG hotspot/domain; selected VCEP | `test_pm1_*` | Curated domain evidence required |
| PM2 | Active, scoped | ACMG rarity; SVI PM2 Supporting v1.0; selected VCEP | `test_pm2` | Population quality and scope matter |
| PM3 | Active, scoped | ACMG trans phase; SVI points | `test_pm3_missing_evidence`, phase tests | No phase inferred from co-occurrence |
| PM4 | Active, scoped | ACMG length change/nonrepeat; selected VCEP | `test_pm4_missing_evidence`, indel tests | Repeat/LoF overlap matters |
| PM5 | Active, scoped | ACMG different missense at same residue; selected VCEP | `test_pm5_engine` | Independent reference classification |
| PM6 | Active via PS2_PM6 | SVI de novo v1.1, assumed parentage | `test_ps2pm6_engine` | Never scored alongside PS2 from same pool |
| PP1 | Active, scoped | ACMG segregation; selected VCEP | `test_pp1` | Family observations must be valid |
| PP2 | Active, scoped | ACMG missense mechanism; selected VCEP | `test_pp2` | Disease-specific mechanism required |
| PP3 | Active, scoped | SVI calibrated REVEL/SpliceAI; RUNX1 | `test_pp3_*` | Prediction overlap with PVS1 excluded when known |
| PP4 | Active, scoped | ACMG phenotype specificity | `test_pp4_engine` | No universal de novo phenotype reuse ban |
| PP5 | Disabled | SVI 2018 discontinue | `test_pp5` | External label is not independent evidence |
| BA1 | Active, scoped | ACMG/SVI frequency and exceptions; selected VCEP | `test_ba1` | Exceptions and validated population record required |
| BS1 | Active, scoped | Disease-specific frequency; selected VCEP | `test_bs1` | Not a universal numeric threshold |
| BS2 | Active, scoped | ACMG healthy adult genotype; selected VCEP | `test_bs2` | Penetrance/onset and clinical assessment required |
| BS3 | Active, scoped | ACMG/SVI functional assays; selected VCEP | `test_bs3_*` | Independent validated normal assay required |
| BS4 | Active, scoped | ACMG nonsegregation; selected VCEP | `test_bs4` | Phenocopies/penetrance assessed |
| BP1 | Active, scoped | ACMG truncating mechanism; selected VCEP | `test_bp1_*` | Curated disease mechanism required |
| BP2 | Active, scoped | ACMG cis/trans with pathogenic variant | `test_bp2` | Phase and comparison classification required |
| BP3 | Active, scoped | ACMG repeat-region in-frame indel; selected VCEP | `test_bp3_*` | Functional-region exception matters |
| BP4 | Active, scoped | SVI calibrated computational benign evidence; selected VCEP | `test_bp4_*` | PP3 contradiction needs review |
| BP5 | Active, scoped | ACMG alternate molecular diagnosis | `test_bp5_engine` | Alternative diagnosis must be established |
| BP6 | Disabled | SVI 2018 discontinue | `test_bp6` | External label is not independent evidence |
| BP7 | Active, scoped | ACMG synonymous/no splice effect; selected VCEP | `test_bp7_engine` | Splice assessment required |

## Guideline-to-code verification matrix

| Published rule | Code and checked behavior | Result / limit |
|---|---|---|
| [ACMG/AMP 2015 Tables 3–5](https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/) | `acmg/classification_rules.py`; all Table 5 branches have explicit boundary tests | Synthetic combinations pass; conflicting pathogenic/benign evidence returns VUS |
| [SVI PM2 v1.0 (2020)](https://www.clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf) | `adaptive_pm2`; PVS1 + one Supporting rule in combiner | Implemented; generic Table 5 alone would leave this VUS |
| [SVI de novo v1.1](https://www.clinicalgenome.org/site/assets/files/3461/svi_proposal_for_de_novo_criteria_v1_1.pdf) | `ps2pm6`; one pooled PS2/PM6 result | Synthetic pathway tests pass |
| [SVI PVS1](https://pmc.ncbi.nlm.nih.gov/articles/PMC6185798/) and [splicing](https://pmc.ncbi.nlm.nih.gov/articles/PMC9980257/) | `pvs1`; `adaptive_pp3` | Implemented branches documented in `coverage-matrix.md` and `pp3-guideline-matrix.md`; no complete VCEP catalogue |
| [SVI PP5/BP6 discontinuation](https://pubmed.ncbi.nlm.nih.gov/29543229/) | Registry and classifier | Both disabled |
| [ClinGen CSpec](https://cspec.clinicalgenome.org/cspec/ui/svi/summary) | Pinned rules in criterion packages | Only selected gene/disease/version scopes; no universal VCEP combination engine |

Criterion-specific source/threshold/exception mappings live in the existing
`docs/*guideline-matrix.md` and decision-table files. This audit did not
independently recalculate every threshold in those files against each current
CSpec revision; those paths must not be described as exhaustively verified.

## Defects corrected and dependency audit

- The unified package was not included in setuptools discovery. It is now
  included. The classifier now rejects explicitly conflicting variant IDs
  across supplied criterion evidence rather than combining different variants.
- The 2015 Very Strong combination is anchored to PVS1; a non-PVS1 Very
  Strong result has no automatic 2015 Table 5 rule. ClinGen's published
  PVS1 + Supporting Likely Pathogenic combination is implemented separately.
- PS2/PM6: one pooled evaluation. PS3/BS3, PP3/BP4, PP1/BS4, PM4/BP3,
  BA1/BS1/PM2: `evidence_dependencies.py` records or suppresses known
  contradictions/reuse. PVS1/PP3 overlap is suppressed when explicitly
  identified. PM3/BP2 and PP4/de novo can involve distinct observations;
  no universal suppression is scientifically justified. Opposing activated
  evidence yields VUS pending review.
- Cross-evaluator observation/source IDs are not standardized. The resolver
  cannot prove independence when IDs are omitted. Every activated result is
  preserved in `criteria_activated`; suppressed rows are marked `counted=False`.

## Published-reference and missing-data results

| Published case | Reference evidence available | Engine comparison |
|---|---|---|
| [SLC6A1 NM_003042.4:c.1070C>A](https://pmc.ncbi.nlm.nih.gov/articles/PMC7431429/), de novo example | Confirmed parentage, nonspecific seizures, one de novo point; no complete variant classification dossier | Incomplete; no concordance claim |
| [ABCA4 c.5603A>T, p.Asn1868Ile](https://pmc.ncbi.nlm.nih.gov/articles/PMC8700904/) | Hypomorphic BA1 exception described; no complete validated population-record/criterion dossier | Incomplete; no forced benign call |
| [TYR c.1205G>A, p.Arg402Gln](https://pmc.ncbi.nlm.nih.gov/articles/PMC8700904/) | High frequency plus context-dependent hypomorphic effect; no complete criterion dossier | Incomplete; no forced benign call |

The tests cover empty input, population-only BA1, computational-only PP3,
missing essential PVS1 evidence, deprecated codes, de novo pooling, opposing
evidence, and unsupported VCEP declarations. They do not provide measured
clinical sensitivity, specificity, or concordance across P/LP/VUS/LB/B
expert-curated variants.

## Regression and remaining limitations

`python -m pytest -q`: **1083 passed, 66 subtests passed** on the audit
snapshot. The unified classification and dependency tests passed as well.

Remaining scientific work: obtain source evidence dossiers and run blinded
criterion-level comparisons; verify each selected VCEP specification and
version against CSpec; standardize evidence observation IDs and provenance
across evaluators; audit all 28 positive/negative/conflicting/gene-exception
pathways against published sources. Until then, this engine is decision
support with synthetic-rule coverage, not a clinically validated classifier.
