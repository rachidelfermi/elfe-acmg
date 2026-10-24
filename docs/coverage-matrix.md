# PVS1 guideline-to-code coverage matrix

## Scope and authority

This is a PVS1-only decision-support component. It does not annotate variants,
query external databases, or combine ACMG/AMP criteria into a variant
classification. The generic branches are deliberately limited to the published
ACMG/AMP and ClinGen SVI decision trees:

- Richards et al. 2015, [ACMG/AMP PVS1 eligibility](https://pubmed.ncbi.nlm.nih.gov/25741868/).
- Abou Tayoun et al. 2018, [ClinGen SVI PVS1 Figure 1 and Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC6185798/).
- Walker et al. 2023, [ClinGen SVI splicing/PVS1(RNA) guidance](https://pmc.ncbi.nlm.nih.gov/articles/PMC10357475/).
- ClinGen CSpec documents, [VCEP specification catalogue](https://cspec.genome.network/cspec/ui/svi/summary), when a caller supplies a pinned declaration.

`ACTIVATED`, `NOT_MET`, `UNDETERMINED`, and `NOT_APPLICABLE` are deliberately
distinct. An omitted feature is unknown, not a negative finding.

## Feature-impact audit

| Feature group | Present | Decision impact | Applicable pathway | Action / current handling |
| --- | --- | --- | --- | --- |
| Curated LoF mechanism; gene-disease validity; SVI Table 1 facts | Yes | Critical | All generic branches | Keep. A direct curated conclusion is sufficient; Table 1 is a fallback and can cap strength. A known non-LoF mechanism is `NOT_MET`; incomplete fallback facts are not negative evidence. |
| Gene, disease, transcript/version, inheritance | Yes | Critical for scoped rules; high otherwise | VCEP; transcript-sensitive generic paths | Keep. Versioned VCEP scope matches these facts exactly; generic transcript relevance is separately assessed. |
| Relevant transcript/exon, selection policy, alternative-transcript rescue | Yes | Critical | N/F, splice, exon-level CNV, duplication, initiation | Keep. Any direct relevance finding or documented policy can establish relevance. Known irrelevance or credible rescue excludes the applicable generic branch; contradiction is `UNDETERMINED`. |
| Population evidence that exon LoF is tolerated | Yes | High | NMD-escape/in-frame role-unknown percentage ladder | Keep, but only in this branch. It does **not** veto predicted-NMD or critical-region branches. |
| NMD prediction | Yes | Critical | N/F; out-of-frame splice/deletion/duplication | Keep. It must be supplied; position or junction-distance inputs never manufacture an NMD conclusion. |
| Protein amount/region affected | Yes | High | NMD escape; in-frame splice; in-frame exon deletion | Keep. `critical_protein_region` supports the applicable Strong branch; the generic role-unknown ladder uses only published `>10%` / `<10%` outcomes and leaves exactly 10% unresolved. In-frame genomic exon deletion has no generic percentage ladder. |
| Canonical splice location, predicted skipping/cryptic use, frame, resolved outcomes | Yes | Critical | Canonical splice | Keep. A canonical position alone never assumes exon skipping or an out-of-frame product. Generic prediction needs an explicit splice signal plus an explicit frame result. |
| Direct RNA splice evidence, assay interpretability, LoF sufficiency, residual normal transcript, functional rescue | Yes | Critical for non-canonical/RNA-confirmed; high for canonical | PVS1(RNA) splice | Keep. Positive RNA mode requires observed/curated assay provenance, an interpretable assay, an LoF-sufficient profile, and an observed/curated allele-relevant splice outcome. It can route an unpositioned splice-site call or an otherwise non-null sequence annotation when the direct RNA evidence establishes the splice event. Mixed normal/aberrant products require a gene-disease-specific residual-function conclusion. |
| RNA-confirmed in-frame loss of clinically relevant residues | Yes | High | PVS1(RNA) in-frame branch | Keep as a distinct RNA-specific route; it is not inferred from the generic critical-region Boolean. |
| Initiation-codon effect, alternative start, upstream pathogenic variants | Yes | Critical | Initiation codon | Keep. Functional alternative initiation excludes generic PVS1; upstream pathogenic evidence changes Moderate versus Supporting strength. |
| Whole-gene versus exon deletion, affected exons, frame | Yes | Critical | Deletions | Keep. Whole-gene deletion is Very Strong after LoF eligibility; exon deletion requires actual exon extent and a supplied frame result. A conflicting full-gene assertion alongside a known small-variant, start-codon, splice-position, or affected-exon anchor is a material preflight conflict, not an ignored extra field. |
| Duplication containment, whole-gene status, affected exons, tandem status, frame, NMD | Yes | Critical | Duplications | Keep. The generic route requires explicitly non-whole-gene, intragenic structure and a proven/presumed tandem assessment. Unknown structure is `UNDETERMINED`, not presumed null. |
| Versioned VCEP rules, source, scope, decision conditions | Yes | Critical when supplied | Any supported VCEP pathway | Keep. Caller-supplied pinned manifests take precedence over generic SVI logic; conflicting or potentially applicable unresolved rules return `UNDETERMINED`. A path-specific rule outranks `ANY` only inside the same specification; separate supplied specifications have no implicit scope/version ranking. If generic routing is incomplete, an `ANY` rule may apply, while a path-specific rule must match a compatible `variant_consequence`, CNV/whole-gene, start-codon, canonical/non-canonical `splice_site_position`, or structural `affected_exons` anchor so an unrelated rule cannot hijack the case. |
| Provenance and conflict state | Yes | High | All paths | Keep. Evidence has `KNOWN`, `UNKNOWN`, or `CONFLICTING` state plus observed/predicted/derived/curated provenance; the result retains only decision-critical `evidence_used` snapshots and provenance rather than adding a score. Generic output also records its framework/branch/source/evaluated conditions, while VCEP output records its exact `scope` and `when` match. |
| Coordinates, exact PTC/junction distance, exon location, raw domain evidence | Retained, not generically derived | Moderate context | Upstream annotation / VCEP | Retain for audit and VCEP conditions. The engine intentionally does not create NMD, transcript relevance, or a critical-region conclusion from these raw values. |
| Tissue expression/pext, constitutive exon use, natural exon skipping, translation reinitiation, assay tissue quality | No separate numerical model | Potentially high, context-specific | Transcript/splice/NMD escape | Do not invent a universal threshold. Supply a curated relevance/rescue/RNA-interpretable conclusion or a scoped VCEP condition where the published gene-disease specification defines how it changes PVS1. |
| Breakpoint sequence, duplication orientation/insertion topology, complex rearrangement | No full structural model | Potentially high | Complex CNV/duplication | Do not infer a null effect. The generic duplication route accepts only the documented categorical structure facts; unresolved or more complex cases need a VCEP rule or expert assessment. |

## Guideline branch-to-code coverage

| Guideline branch | Implemented behavior | Main code path / evidence | Coverage boundary |
| --- | --- | --- | --- |
| ACMG/AMP PVS1 eligibility | Evaluates null-variant evidence only in a gene-disease assertion with an eligible LoF mechanism. | `_mechanism_assessment`; `lof_disease_mechanism`, SVI Table 1 fields | Does not classify the variant or combine other ACMG criteria. |
| SVI Table 1 mechanism tiers | Direct established/capped conclusion, or the applicable published Table 1 combinations; applies a strength cap where needed. | `_mechanism_assessment` | No additive point scoring and no inference from partial Table 1 facts. |
| SVI Figure 1 N/F, predicted NMD | Relevant N/F event with supplied predicted NMD activates Very Strong. | `_evaluate_null_variant`, `_nmd_or_escape` | NMD is supplied evidence, not calculated from coordinates. |
| SVI Figure 1 N/F, NMD escape | Critical region gives Strong; otherwise the non-tolerant, role-unknown `>10%` / `<10%` ladder gives Strong/Moderate. | `_non_nmd_ladder` | Exactly 10% is deliberately unresolved; no invented threshold fills the gap. |
| SVI Figure 1 canonical +/-1/2 splice | Evaluates the actual/predicted resolved outcome; out-of-frame uses the NMD branch and in-frame uses the in-frame branch. | `_identify_pathway`, `_resolve_splice_outcomes`, `_evaluate_splice` | Canonical position alone is insufficient to assume a product or frame. |
| Walker et al. PVS1(RNA), non-canonical/direct-RNA-confirmed splice | A usable observed RNA profile can enter the resolved-outcome tree even when a splice-site position is unavailable or the upstream sequence annotation is otherwise non-null; direct normal/no-LoF findings return `NOT_MET`. | `_assess_rna_profile`, `_evaluate_splice` | Prediction-only or unprovenanced RNA never activates RNA mode. An explicit positional conflict remains unresolved rather than being silently overridden. RNA-mode output cites 2023 direct-splice guidance and records any 2018 Figure 1 strength row as a separate source. No generic RNA-abundance cutoff is claimed. |
| RNA mixed products and in-frame products | Assesses residual normal transcript sufficiency; chooses the lowest strength across credible damaging outcomes. A specific RNA-confirmed clinically relevant in-frame loss has its separate route. | `_evaluate_splice`, `_evaluate_in_frame_splice` | The caller must supply the gene-disease functional conclusion; the engine does not estimate functional sufficiency from an RNA percentage. |
| SVI Figure 1 initiation codon | No functional alternative start plus upstream pathogenic evidence gives Moderate; absent such evidence gives Supporting. | `_evaluate_initiation_codon` | Alternative start and upstream evidence are required only for this path. |
| SVI Figure 1 exon deletion | Out-of-frame deletion follows NMD/NMD-escape logic; in-frame deletion needs a critical region for generic Strong. | `_evaluate_exon_deletion` | No generic in-frame-deletion percentage ladder. |
| SVI Figure 1 full-gene deletion | Eligible whole-gene deletion gives Very Strong without a transcript relevance requirement. | `_evaluate_whole_gene_deletion` | Complex copy-number calls must first establish whole-gene status. |
| SVI Figure 1 intragenic duplication | Proven tandem/disrupted-frame/NMD yields Very Strong; presumed tandem yields Strong. | `_evaluate_duplication` | Non-tandem, in-frame, non-NMD, whole-gene, or non-intragenic events do not receive generic PVS1. |
| ClinGen VCEP specifications | Evaluates a versioned, gene/disease-scoped declarative rule before the generic route and reports that source. | `VCEPSpecification`, `select_specific_rule` | No mutable or universal VCEP database is bundled. If generic routing is unknown, `ANY` remains eligible; a path-specific exception needs an explicit pathway anchor in its conditions. |

## Remaining supported-boundary limitations

- The package has no annotation, transcript-selection, NMD-prediction, RNA-quantification, population-frequency, or database-query pipeline. It evaluates supplied, traceable conclusions.
- There is no built-in catalogue of PTEN or other VCEP rules. A published VCEP specification must be represented by a caller-pinned, versioned manifest with its original source URL and exact scope; this package does not attest that a caller's manifest transcribes a source correctly.
- It intentionally has no universal pext/exon-expression, residual-normal-RNA, assay-quality, downstream-truncation, or translation-reinitiation numerical cutoff. Such thresholds are gene/disease-specific when valid and belong in a governed VCEP declaration or a curated input assertion.
- Multiple transcripts/products are represented by supplied relevance, rescue, splice-outcome, and residual-function conclusions rather than a complete isoform-expression model. Contradictory material input is returned as `UNDETERMINED`.
- Detailed CNV breakpoints, duplication orientation/insertion site, complex rearrangements, and allele phase are not structurally reconstructed. The generic duplication/deletion paths require the limited facts published for their branches and otherwise stop safely.
- `variant_type`, coding/termination positions, junction distance, exon location, and raw functional-region inputs are retained for provenance and VCEP conditions but do not independently alter a generic decision. This avoids silently inventing NMD or functional relevance.
- Tests are synthetic decision-engine tests, not clinical validation cases. A production deployment still needs controlled source snapshots, review by domain experts, release governance, and validation on its intended use population.
