# BP1 source-to-rule matrix

| Rule | Primary source | Implemented decision |
| --- | --- | --- |
| General BP1 | [Richards et al., 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/) | Missense in a gene for which primarily truncating variants cause the relevant disease; Supporting only. |
| Gene-centric assessment | [ClinGen variant curation SOP v2.0, BP1 section](https://www.clinicalgenome.org/site/assets/files/5933/variant_curation_sop_2_0_jan_2021.pdf) | Require a disease-specific curated assessment. Neither LOF intolerance nor raw counts alone establish BP1. |
| VCEP use of BP1 | [ClinGen standardized VCEP text v3.0](https://clinicalgenome.org/site/assets/files/6773/vcep_general_standardized_summary_text_v3_dec_2025.pdf) | A scoped VCEP can define BP1 applicability; do not extend a gene rule to unrelated diseases. |
| DICER1 tumor predisposition, GN024 v1.4 | [CSpec GN024](https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN024) | BP1 Not Applicable for DICER1 / MONDO:0100216. |
| RYR1 malignant hyperthermia, GN012 v2.0 | [CSpec GN012](https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN012) | BP1 Not Applicable for RYR1 / MONDO:0018493. |
| MAP2K1 RASopathy, GN045 v2.3 | [CSpec GN045](https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN045) | Disease-specific reversal: specified truncating/LOF classes in MAP2K1 / MONDO:0021060, with missense GOF mechanism and no established LOF correlation; Supporting. |

The general mechanism route accepts a *supplied* documented qualitative conclusion
that truncating variants predominate. The spectrum route likewise accepts a
curated qualitative conclusion and checks disease match, independent pathogenic
classifications, unique variants, and addressed ascertainment bias. The evaluator does not infer
predominance from counts, fractions, constraint, population tolerance, or rarity:
the cited sources provide no universal quantitative cutoff for BP1. A pathogenic
missense variant, affected region or residue, relevant missense GOF/DN mechanism,
or established damaging splice effect prevents general activation. High-risk
splice prediction leaves the decision undetermined pending resolution. A prior
PP2 assignment also leaves BP1 undetermined until the contradictory mechanism
premises are resolved. PP3 and PM1 are recorded but do not independently change
BP1 because no general combination prohibition was verified.

All tests using TEST1 and SYNTHETIC disease IDs are synthetic rule fixtures. The
three pinned CSpec tests check published applicability statements; they are not
externally validated variant classifications. No general BP1 strength upgrade is
implemented because none was verified for these scoped rules. Other VCEP
specifications are not preloaded, so users must verify applicability independently
before using this evaluator for a gene–disease pair outside the pinned set.
