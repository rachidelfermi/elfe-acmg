# Evidence feature contract

The unified API accepts already extracted evidence under `criteria`, with one
payload per criterion. Each evaluator owns its biological fields. Shared names
below describe comparable facts, not a universal input schema or permission to
count one observation twice. The [criterion coverage and limits](unified-scientific-validation.md)
and criterion-specific matrices define the actual rule conditions.

| Evidence concept | Existing fields and consumers | Decision role |
| --- | --- | --- |
| Variant identity and consequence | `variant_id`, `variant_consequence` or criterion-specific `consequence`; PVS1, PS1, PM4, PM5, PP2, BP1, BP3, BP7 | Establishes the eligible variant effect and guards against cross-criterion identity mismatch. Consequence labels have pathway-specific meanings. |
| Gene, disease, mechanism and transcript | `gene_id`, `disease_id`, `disease_mechanism`, `biologically_relevant_transcript` and criterion-specific validity fields; PVS1, PS1, PS3, PP2, BP1, BP4 and others | Selects applicable disease and VCEP scope; supports or excludes the relevant mechanism. |
| Population frequency | `population_records` (PM2), `populations` (BA1), `frequencies` (BS1); `allele_frequency` and `filtering_allele_frequency` have different uses | PM2 rarity and BA1/BS1 benign frequency require their own source, quality, variant identity and threshold checks. These collections are not interchangeable aliases. |
| Computational prediction | `predictions`, REVEL or SpliceAI scores; PP3, BP4, BP7, and splice-related PVS1 | A calibrated result can be sufficient where its rule permits; reuse and contradictory directions need review. Predictor names alone do not prove observation identity. |
| Functional assay | `functional_assays` (PS3), `assays` (BS3), assay identifiers and calibration | One qualifying assay can meet a supported pathway; opposing use of the same assay is resolved before combination. |
| De novo and family observations | `observations`, case and family identifiers; PS2/PM6, PP1, BS4 | Parentage, phenotype and case independence affect the applicable criterion and strength; family observations retain criterion-specific interpretations. |
| Phase | `phase`, `phase_method`, reliability and family facts; PM3, BP2 | Reliable supplied molecular phase can be used without inventing parental calls. Contradictory phase remains unresolved. |
| Protein length and repeat region | `protein_length_change`, `repeat_region_status` and functional context; PM4, BP3 | A repeat-region interpretation can change the direction of evidence; shared events are reconciled. |
| Phenotype | `phenotype_category` or criterion-specific phenotype evidence; PS2/PM6, PP4 | Changes de novo points or phenotype specificity; the same phenotype observation cannot be presumed independent. |

Supported input aliases remain local to their evaluator. For example, the
shared de novo parser accepts `cases` and `independent_cases` for
`observations`; BS3 accepts `disease_molecular_mechanism` for
`disease_mechanism`. Do not apply those aliases across unrelated criteria.
Evidence wrappers retain `KNOWN`, `UNKNOWN`, `CONFLICTING`, and
`NOT_APPLICABLE` where supported, plus source and version when supplied. The
unified boundary converts typed nested wrappers and ignores fields unrelated
to the selected typed evaluator. An unknown essential fact cannot activate a
criterion; missing optional facts do not block a complete pathway.

## Minimum evidence and dependencies

The minimum is pathway-specific: PS1 may use one qualifying reference variant;
the shared PS2/PM6 engine may use one eligible independent de novo case; PS3
may use one sufficiently validated disease-relevant assay; PP3 may use one
calibrated predictor; and BS1 may use a verified disease-specific frequency
threshold without a new prevalence calculation. Each pathway still enforces
its own essential eligibility, exclusion, quality and strength conditions.
The [source matrices](unified-classification-source-matrix.md) and the
criterion-specific decision tables give those conditions; this summary adds
no thresholds.

The [dependency resolver](../acmg/evidence_dependencies.py) records shared
PS3/BS3 assays, PP3/BP4 predictors, PP1/BS4 families, PVS1/PP3 splice reuse,
PM4/BP3 length conflicts, and BA1/BS1/PM2 population conflicts. The shared
de novo evaluator assigns either PS2 or PM6. Criterion evaluators also check
their own phase, case and phenotype overlap. Where independence cannot be
established from the supplied evidence, the engine does not invent it.

## Supported boundary

PP5 and BP6 remain inactive. The classifier supports its current ACMG/AMP
2015 combinations and ClinGen SVI PM2 update; complete VCEP-wide and Bayesian
combining remain unsupported. The standalone PM6 interface retains separately
verified VCEP behavior and is not used by the classifier. No variant
annotation, database extraction or clinical validation is provided here.
