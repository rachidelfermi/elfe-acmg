# BS3 guideline to code matrix

| Published recommendation | Code | Minimum evidence and limit |
| --- | --- | --- |
| ACMG/AMP 2015 BS3: well-established in vitro or in vivo studies show no damaging effect | `evaluator._evaluate_assay` | Normal disease-relevant function in a suitable assay; a normal observation alone is insufficient. |
| Brnich et al. 2019, four-step functional assay framework | `evaluator._evaluate_assay` | Disease mechanism, assay applicability, validation, variant result; known artifact, saturation, floor effects and incompatible readouts exclude activation. |
| Brnich et al. 2019, Table 3 | `guidelines.benign_strength_from_oddspath` | Validated formal OddsPath: `<0.48` supporting, `<0.23` moderate, `<0.053` strong. Boundaries are strict. No automatic very strong. |
| Brnich et al. 2019, Step 3 control-count path | `guidelines.control_count_strength` | Mixed independently classified benign and pathogenic variant controls, WT control, separation and replicates; 11 total permit moderate, fewer permit supporting. Never strong by control count alone. |
| Brnich et al. 2019, historically accepted assay path | `evaluator._evaluate_assay` | Documented historical acceptance plus an accepted normal interpretation permits supporting only. |
| Walker et al. 2023, direct RNA-splicing assays | `evaluator._evaluate_assay` | Report outside BS3 under BP7(RNA); tissue, transcript, detection, NMD and variant-allele expression are surfaced as adequacy considerations. |
| Glaucoma VCEP MYOC specification | `evaluator._evaluate_assay` | Scope MYOC and primary open-angle glaucoma; secretion or solubility, calibrated OddsPath, and a moderate maximum. |
| Brnich et al. 2019, multiple studies | `evaluator.evaluate_bs3` | Evaluate studies individually; never add strengths; unresolved credible normal/damaging conflict gives `UNDETERMINED`. |

Sources: [Richards et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/), [Brnich et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6938631/), [Walker et al. 2023](https://pmc.ncbi.nlm.nih.gov/articles/pmid/37352859/), [MYOC Glaucoma VCEP publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC9771967/).

The MYOC rule is the only bundled VCEP rule. Other VCEP specifications are not inferred from a gene name. The engine evaluates supplied evidence offline; it does not authenticate a caller's source claims or calculate OddsPath from raw control results.
