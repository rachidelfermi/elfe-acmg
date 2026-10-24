# PP3 guideline-to-code map

This offline evaluator receives already-extracted predictions. It does not classify variants, run predictors, or activate BP4.

| Rule | Source | Eligible input and minimum evidence | PP3 result |
| --- | --- | --- | --- |
| General missense REVEL | [Pejaver et al. 2022, Table 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/) | Missense; relevant versioned transcript; matched, valid REVEL prediction and predictor version; no known technical exclusion | REVEL ≥0.644 Supporting, ≥0.773 Moderate, ≥0.932 Strong. Boundary inclusive. |
| General noncanonical splicing | [Walker et al. 2023, Table 1 and recommendations](https://pmc.ncbi.nlm.nih.gov/articles/PMC10357475/) | Eligible noncanonical variant; relevant versioned transcript and splice context; matched, valid SpliceAI prediction; no known splice/PVS1 RNA overlap | SpliceAI maximum delta ≥0.20: PP3 Supporting. Walker's Table 1 calibrates **spliceogenicity**, not pathogenic OddsPath; this implementation conservatively retains PP3 Supporting. |
| RUNX1 VCEP | [ClinGen CSpec GN008 v3.1, PP3](https://cspec.genome.network/cspec/ui/svi/svi/GN008) | RUNX1, MONDO:0011071, NM_001754.4; explicitly selected GN008 v3.1; missense for REVEL; missense, synonymous, or intronic intron 4–8 for SpliceAI; no canonical ±1/2 site or RNA proven splice effect | REVEL ≥0.88 or SpliceAI ≥0.38: Supporting only. VCEP overrides general thresholds and strengths. |

The original [ACMG/AMP PP3 criterion](https://pubmed.ncbi.nlm.nih.gov/25741868/) supplies the criterion definition. The calibrated Pejaver recommendation supports one prospectively selected tool instead of voting across predictors. Set `selected_missense_predictor="REVEL"` before evaluating scores; REVEL is the only implemented missense policy. Prediction records require `calibration="PEJAVER_2022"` (or `WALKER_2023` for SpliceAI). The engine does not sum missense and splice predictions. An explicit applicable PP3/BP4 conflict returns `UNDETERMINED`. Low scores are reported as BP4 candidates only and are never activated as BP4.

Scores from other listed tools, conservation, and structural predictions are accepted as audit context but cannot activate PP3 under this implemented rule set. Additional VCEP versions require source verification and a pinned rule implementation. The engine does not turn a high SpliceAI score into PVS1 or claim a specific RNA consequence.

Prediction records may provide variant, transcript, and genome build identifiers. When an input identifier is supplied, a prediction identifier must match before activation. Known poor quality, inapplicability, out-of-distribution status, absent score, and version mismatch prevent activation. Omitted optional predictors are irrelevant to the selected pathway.
