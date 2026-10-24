# BP4 source-to-rule matrix

| Source | Encoded rule | Scope and limits |
| --- | --- | --- |
| [Richards et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/) | BP4 computational evidence is counted once, with applicability checked. | The 2015 text gives no universal numerical cutoffs. |
| [Pejaver et al. 2022, Table 2](https://pubmed.ncbi.nlm.nih.gov/36413997/) | REVEL: ≤0.290 supporting, ≤0.183 moderate, ≤0.016 strong. BayesDel without allele frequency: ≤−0.18 supporting, ≤−0.36 moderate. | Missense only; individual calibrated tools may suffice. REVEL ≤0.003 is labelled very strong in the paper but capped at strong by this API. Other tools remain unsupported until explicitly transcribed and checked. |
| [Walker et al. 2023](https://pubmed.ncbi.nlm.nih.gov/37352859/) | SpliceAI ≤0.1 gives moderate evidence against spliceogenicity for eligible noncanonical splice variants. | Caller must verify site eligibility, relevant transcript and absence of unresolved protein effects. The score says nothing about protein function. |
| [RPE65 VCEP 1.0 BP4](https://cspec.clinicalgenome.org/cspec/ui/svi/doc/1566365955?version=1.0.0) | Missense: REVEL ≤0.183 and SpliceAI ≤0.1 gives moderate; REVEL >0.183 to ≤0.290 and SpliceAI ≤0.1 gives supporting. Splice-only rule supports BP4 at ≤0.1. | Restricted to RPE65, MONDO:0100368, NM_000329.3 and version 1.0. Splice-only site/location must be reviewed by caller. |
| [PALB2 VCEP 1.3 BP4](https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN077) | Missense BP4 is prohibited. Splice-only SpliceAI ≤0.1 gives supporting BP4. | Restricted to PALB2, the two listed MONDO diseases, NM_024675.3 and version 1.3. |

No variant examples in the tests are published clinical cases; all fixtures are synthetic boundary and failure examples. This evaluator does not determine gene–disease validity, run predictors, discover VCEP rules, or assert that a prediction is experimentally validated. Inputs asserting identity, transcript mapping, model applicability, and site eligibility must be independently checked.
