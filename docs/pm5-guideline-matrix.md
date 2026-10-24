# PM5 guideline to code coverage

Inputs are curated facts. This package does not infer transcript mappings, clinical classifications, splice predictions, or Grantham distances.

The [ClinGen SVI guidance index](https://clinicalgenome.org/working-groups/sequence-variant-interpretation/) was reviewed for general PM5 recommendations. It lists no standalone PM5 strength algorithm, so the general pathway remains the 2015 ACMG/AMP rule; only the two cited VCEP specifications modify it here.

| Path | Primary source / pinned version | Minimum sufficient evidence | Exclusions and missing evidence | Permitted strength | Code |
| --- | --- | --- | --- | --- | --- |
| Standard same residue, different missense | [Richards et al. ACMG/AMP 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/), Table 4 and PM5 discussion | Missense candidate and independently established Pathogenic missense reference; same relevant transcript, residue and reference amino acid; different alternate amino acid; matching disease | Same substitution routes to PS1 review; splice-mediated reference, incompatible disease/mechanism, circular classification, or invalid mapping excludes; absent reference review is unresolved | Moderate | `evaluator._reference_check`, `evaluate_pm5` |
| Multiple references under original ACMG | [Richards et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/) | At least one qualifying distinct Pathogenic reference | Duplicate reports of the same amino acid substitution count once; no original rule for automatic increase | Moderate | `evaluate_pm5` |
| HNF1A monogenic diabetes | [ClinGen Monogenic Diabetes VCEP v1.1](https://www.clinicalgenome.org/site/assets/files/7133/clingen_diabetes_acmg_specifications_v1_1.pdf), PM5 rows | HNF1A, MONDO:0015967, NM_000545.8; one qualifying Pathogenic reference with candidate Grantham distance at least as great as reference | One Pathogenic reference with greater Grantham distance, or a Likely Pathogenic reference, supports the VCEP Supporting route; unknown required distances leave one Pathogenic route unresolved | Moderate or Supporting | `evaluate_pm5` |
| HNF1A two changes | Same VCEP v1.1, PM5_Strong row | Two distinct qualifying Pathogenic amino acid changes at the residue | Two reports or nucleotide variants producing the same amino acid change count once | Strong | `evaluate_pm5` |
| ITGA2B/ITGB3 Glanzmann thrombasthenia | [ClinGen Platelet Disorders VCEP v1](https://clinicalgenome.org/docs/clingen-platelet-disorders-expert-panel-specifications-to-the-acmg-amp-variant-interpretation-guidelines-version-1/), PM5 and PM5_Supporting | Matching gene, MONDO:0010119, specified transcript and one qualifying reference | Pathogenic reference uses PM5 as written; Likely Pathogenic reference uses Supporting | Moderate or Supporting | `evaluate_pm5` |

## Scientific boundaries

- The source classification and its underlying evidence must be reviewed outside this package. `classification_evidence_reviewed` and `classification_independent_of_pm5` are assertions by the curator, not inferred from ClinVar stars. A bare ClinVar label is insufficient.
- The source ACMG paper warns that a missense label may mask a splice mechanism and that different substitutions can yield different phenotypes. Known splice-mediated reference pathogenicity blocks PM5; unknown splice effects are not interpreted as normal splicing. A known incompatible candidate mechanism also blocks PM5.
- For general ACMG, a Likely Pathogenic reference alone does not satisfy the wording “pathogenic”; it remains undetermined for PM5. The two pinned VCEP rules explicitly permit weaker routes.
- Protein accession mismatch and transcript mismatch are not silently reconciled. A separate, verified mapping must be supplied as a common transcript and protein context.
- A known benign substitution at the residue or conflicting biological evidence returns `UNDETERMINED` for manual review; it is not converted to a numerical penalty.
- No other VCEP is encoded. The default ACMG rule applies to genes outside the three pinned VCEP scopes; users must select and verify any other applicable disease-specific VCEP before clinical use. The included VCEP versions are historical pinned specifications and should be checked against current versions for clinical work.
- `features_used` and `exclusions_evaluated` describe decision predicates, not independent ACMG criteria. PM5 is returned alone and no overall variant classification is made.
