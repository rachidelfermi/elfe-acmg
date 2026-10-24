# BP7 source-to-rule matrix

| Rule | Verified source | Implemented decision | Strength |
| --- | --- | --- | --- |
| ACMG/AMP 2015 BP7 | [Richards et al.](https://pubmed.ncbi.nlm.nih.gov/25741868/) | Synonymous on a relevant transcript; validated prediction of no splice impact; reliable determination that nucleotide is not highly conserved. No numeric conservation or predictor cutoff is supplied by the source, so this path accepts only upstream categorical evidence. | Supporting |
| ClinGen SVI 2023 BP7 | [Walker et al.](https://pubmed.ncbi.nlm.nih.gov/37352859/) | Synonymous outside the first and last three exon bases, or intronic at +7/−21 or farther; a calibrated SpliceAI maximum delta score ≤0.1 meets the BP4 splicing prerequisite. Conservation is not required. | Supporting |
| ClinGen SVI 2023 BP7_Strong (RNA) | [Walker et al.](https://pubmed.ncbi.nlm.nih.gov/37352859/) | Validated, sensitive, relevant RNA assay compared with controls shows no abnormal splicing. The same RNA assay must not also be counted as BS3. | Strong |
| TP53 VCEP v2.1.0 | [ClinGen CSpec GN009](https://cspec.genome.network/cspec/ui/svi/doc/GN009?version=2.1.0) | Limited to TP53 and Li-Fraumeni syndrome; same permitted positions; SpliceAI ≤0.1 and BP4 prerequisite; no conservation requirement. Its BP7_Strong (RNA) is also supported, with the CSpec's BP4 prerequisite retained. | Supporting / Strong |

The published SVI analysis reports a SpliceAI **≤0.1** non-spliceogenic threshold; this is evidence for the BP4 splicing prerequisite, not a Moderate BP7 strength. The engine never derives a predictor score, transcript consequence, conservation status, or RNA result.

The caller supplies transcript-specific, variant-matched evidence. `calibrated=True` means the supplied SpliceAI maximum score was produced and interpreted under the selected published rule; four sub-scores are not recomputed here. A non-SpliceAI tool can supply a validated categorical no-effect assessment for the original ACMG path, which has no numeric threshold. The default is SVI 2023; choosing `ACMG_2015` preserves the original conservation condition. Only TP53 has a built-in VCEP rule. Other VCEP specifications require review and explicit implementation before use. All tests are synthetic rule fixtures; no clinically curated variants are claimed.
