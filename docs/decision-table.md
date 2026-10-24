# PVS1 decision table and evidence contract

## Authority and scope

This package evaluates only PVS1. It does not assign a pathogenicity
classification or combine criteria. The generic implementation is pinned to:

- Richards et al., ACMG/AMP 2015: PVS1 is a null variant in a gene where LoF
  is a known disease mechanism.
- Abou Tayoun et al., ClinGen SVI 2018, PMID 30192042,
  [Figure 1 and Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC6185798/).
- For a non-canonical or otherwise RNA-confirmed splice event, Walker et al.,
  ClinGen SVI Splicing Subgroup 2023, PMID 37352859,
  [PVS1(RNA) guidance](https://pmc.ncbi.nlm.nih.gov/articles/PMC10357475/).

The SVI 2018 paper labels branches in Figure 1 but does not publish formal
per-branch IDs or a semantic-versioned rules API. Accordingly, output uses the
real criterion ID `PVS1`, a pinned document version/citation, a source
locator, and a separate deterministic `decision_path` engine identifier. It
does **not** present local path identifiers as ClinGen rule IDs.

An applicable versioned ClinGen VCEP rule supplied by the caller has precedence
over the generic branch. This package deliberately embeds no incomplete,
mutable VCEP catalogue. A path-specific rule outranks an `ANY` fallback only
within the same supplied specification. Separate matching specifications (or a
potentially applicable unresolved one) have no implicit authority, version, or
scope ranking and return `UNDETERMINED` for review.

## Evidence envelopes

Every feature may be supplied in a compact form:

```json
{"nmd_prediction": true}
```

or with an explicit knowledge state and provenance:

```json
{
  "nmd_prediction": {
    "state": "KNOWN",
    "value": true,
    "provenance": [{
      "source": "validated transcript/NMD annotation",
      "kind": "PREDICTED",
      "identifier": "annotation-release-2026-01"
    }]
  }
}
```

`KNOWN`, `UNKNOWN`, and `CONFLICTING` are distinct. An omitted field is
`UNKNOWN`; it is never coerced to `false`. `kind` is one of `OBSERVED`,
`PREDICTED`, `DERIVED`, or `CURATED`.

The decision output preserves the decision-critical `evidence_used` snapshots,
`features_missing`, evaluated exclusions, and their provenance. Every generic
result records `guideline_conditions` with its SVI/RNA framework, branch,
source locator, and normalized decision-critical values; a matching VCEP rule
instead records its exact scoped `guideline_conditions` (`scope` and `when`).
Unrelated supplied inputs are not copied into the decision audit.

The engine accepts coordinate-like input such as `premature_stop_position` or
`distance_to_final_exon_junction`, but intentionally does not manufacture an
NMD conclusion from it. Supply `nmd_prediction` with the method/provenance
used to derive it. Likewise, it does not infer a splice outcome from a
canonical position or a reading-frame consequence from a duplication.

## Decision-critical schema

Only the relevant subset is needed for an individual branch. All documented
fields are validated; unsupported keys are rejected to prevent misspellings
from silently removing clinical evidence.

| Role | Core fields | Notes |
| --- | --- | --- |
| Eligibility | `gene_id`, `disease_id`, `lof_disease_mechanism`, `gene_disease_validity`, `mode_of_inheritance` | The direct curated LoF conclusion is the smallest sufficient eligibility input. |
| Table 1 fallback | `pathogenic_lof_variant_count`, `pathogenic_lof_variant_fraction_percent`, `pathogenic_lof_variants_across_multiple_exons`, `null_model_recaptures_disease`, `single_exon_gene` | Used only when no curated LoF conclusion is supplied. Count means prior pathogenic LoF variants assessed without PVS1. |
| Variant class | `variant_consequence`, `cnv_event_type`, `whole_gene_deletion`, `whole_gene_duplication`, `start_codon_affected` | Identifies a pathway; it does not itself activate PVS1. |
| Transcript/exon | `biologically_relevant_transcript`, `exon_biological_relevance`, `transcript_selection_policy`, `alternative_transcript_rescue`, `population_evidence_for_exon_lof_tolerance` | Either a relevance finding or documented selection policy can establish relevance. Known alternative-transcript rescue excludes the applicable generic path. Population exon LoF tolerance is evaluated only in the NMD-escape/in-frame role-unknown percentage branch; it does not veto an NMD or critical-region branch. |
| NMD/protein | `nmd_prediction`, `protein_truncation_percent`, `critical_protein_region` | Protein percentage is a percentage removed/altered on the relevant protein. Exactly 10% is deliberately unresolved by the generic rule. |
| Splicing | `splice_site_position`, `splice_outcomes` (or `splice_effect` alias), `rna_splicing_evidence`, `rna_assay_interpretable`, `rna_aberrant_splicing_sufficient_for_lof`, `residual_normal_transcript_sufficient_for_function`, `rna_in_frame_loss_of_clinically_relevant_residues`, `splice_rescue_evidence` | `splice_site_position` is a non-zero signed relative splice offset; `splice_outcomes` is one or more resolved credible outcomes: `OUT_OF_FRAME`, `IN_FRAME`, or `NO_ABERRANT_SPLICING`. Positive PVS1(RNA) use requires observed/curated RNA provenance, an interpretable assay, a LoF-sufficient profile, and an observed/curated allele-relevant outcome; a prediction alone is not RNA evidence. |
| CNV | `affected_exons`, `cnv_reading_frame` (or `reading_frame_disruption`), `duplication_contained_within_gene`, `tandem_duplication_status` | A duplication must be demonstrated intragenic and given an explicit proven/presumed tandem assessment. |
| Start loss | `alternative_start_codon`, `upstream_pathogenic_variant_evidence` | The latter is specifically evidence of pathogenic variants 5′ of the closest downstream potential in-frame start. |
| VCEP override | `gene_specific_pvs1_rules` | Validated versioned snapshot scoped to gene, disease, and optionally transcript/version/MOI. |

The remaining candidate input fields (`variant_type`, `coding_position`,
`premature_stop_position`, `distance_to_final_exon_junction`, `exon_location`,
and `protein_functional_region_evidence`) are retained with provenance for
upstream derivation/audit but do not add arbitrary scoring requirements.
`predicted_exon_skipping` and `cryptic_splice_site_usage` can support the
canonical-splice route only when an explicit reading-frame consequence is also
supplied; neither creates a splice product or NMD conclusion by itself.

## Generic decision table

The table below is a direct implementation map, not an additive feature score.
`PVS1` in a source column means Very Strong in the output.

| Pathway | Minimum sufficient positive conditions | Known exclusion / negative conclusion | Generic strength | Source locator |
| --- | --- | --- | --- | --- |
| Nonsense or frameshift, NMD | Eligible LoF mechanism; NMD predicted; affected exon/transcript biologically relevant | Alternative biologically relevant transcript rescues the exon; irrelevant exon/transcript | Very Strong | SVI 2018 Fig. 1, N/F → predicted NMD |
| Nonsense or frameshift, no NMD | Eligible mechanism; relevant exon; critical region | Alternative-transcript rescue or irrelevant exon/transcript | Strong | SVI 2018 Fig. 1, N/F → not predicted NMD → critical region |
| Nonsense or frameshift, no NMD, role unknown | Eligible mechanism; relevant non-tolerant exon; protein altered `>10%` | Alternative-transcript rescue, irrelevant exon/transcript, or exon LoF tolerance in this percentage branch | Strong | SVI 2018 Fig. 1, N/F → not predicted NMD → `>10%` |
| Same, role unknown | As above but protein altered `<10%` | Same branch-specific exclusions | Moderate | SVI 2018 Fig. 1, N/F → not predicted NMD → `<10%` |
| Canonical ±1/2 splice, out of frame | Resolved outcome is out of frame; then use the matching NMD/non-NMD N/F row | Resolved no-aberrant-splicing, functional rescue, or irrelevant exon/transcript | Very Strong / Strong / Moderate | SVI 2018 Fig. 1, GT–AG ±1/2 |
| Canonical ±1/2 splice, in frame | Relevant exon plus critical region, or role-unknown percentage ladder | Functional rescue or irrelevant exon/transcript; tolerance only in the percentage route | Strong / Moderate | SVI 2018 Fig. 1, GT–AG ±1/2 → preserves reading frame |
| Single/multi-exon deletion, out of frame | Resolved disrupted frame; then use NMD/non-NMD N/F row | Alternative-transcript rescue or irrelevant exon/transcript | Very Strong / Strong / Moderate | SVI 2018 Fig. 1, deletion → disrupts frame |
| Single/multi-exon deletion, in frame | Relevant exon and critical protein region | No generic percentage ladder for this branch | Strong | SVI 2018 Fig. 1, deletion → preserves reading frame |
| Whole-gene deletion | Eligible LoF mechanism | No transcript details needed for this PVS1-only branch | Very Strong | SVI 2018 Fig. 1, full gene deletion |
| Intragenic duplication, proven tandem | Completely within gene; explicitly not whole-gene; frame disrupted; NMD predicted; relevant transcript/exon | Not tandem, in frame, no NMD, whole-gene duplication, or known exclusion | Very Strong | SVI 2018 Fig. 1, duplication → proven tandem |
| Intragenic duplication, presumed tandem | Same facts, with explicit `PRESUMED_TANDEM` assessment | Same | Strong | SVI 2018 Fig. 1, duplication → presumed tandem |
| Start loss | Relevant transcript; no functional alternative start; pathogenic variant(s) upstream of closest downstream in-frame start | Functional alternative start | Moderate | SVI 2018 Fig. 1, initiation codon |
| Start loss | Relevant transcript; no functional alternative start; no such upstream pathogenic variants | Functional alternative start | Supporting | SVI 2018 Fig. 1, initiation codon |
| Non-canonical or direct-RNA-confirmed splice | Observed/curated, interpretable RNA evidence of an aberrant LoF-sufficient profile and observed/curated allele-relevant outcome, followed by the resolved outcome branch | Interpretable direct RNA shows no relevant aberrant splicing; a non-LoF-sufficient profile; no usable RNA/VCEP rule | Per resolved outcome, marked `evidence_mode=RNA` | Walker et al. 2023 PVS1(RNA) |

For multiple credible splice outcomes, every outcome is evaluated and the
lowest activated strength is used. A credible normal/non-LoF outcome with an
unresolved rescue proportion returns `UNDETERMINED`, not a fabricated PVS1
result. The engine never combines PVS1 with PP3 or PS3; those criteria are
outside its scope and should not be double-counted by a calling classifier.

For direct RNA evaluation, `rna_splicing_evidence` is not a generic proxy for
assay quality. The RNA branch requires observed or curated provenance,
`rna_assay_interpretable=true`, and
`rna_aberrant_splicing_sufficient_for_lof=true`. Where credible normal and
aberrant products coexist, the caller must also provide the gene-disease
assessment `residual_normal_transcript_sufficient_for_function`; this engine
does not invent an RNA abundance threshold. A known functional rescue excludes
PVS1. A positive canonical-position prediction can still be evaluated under
the generic SVI route when no usable direct RNA result is supplied. In
particular, an otherwise positive RNA assertion missing assay interpretability
or gene-disease LoF sufficiency remains insufficient for the non-canonical
RNA-only path, but is not an extra requirement for the canonical generic path.
An RNA/profile conflict or an interpretable provenanced negative RNA result is
material and remains terminal rather than being bypassed by a prediction.
`splice_site_position` is needed to classify a generic predicted splice-site
route, but a qualified direct RNA assertion can establish the PVS1(RNA) route
when that position is absent or an upstream sequence annotator supplied an
otherwise non-null consequence (for example, synonymous). An explicitly
conflicting canonical/non-canonical positional assertion remains
`UNDETERMINED`; the engine does not select one coordinate interpretation over
the other.

RNA-mode output cites the 2023 PVS1(RNA) guidance as its primary guideline. If
the observed product then uses the 2018 Figure 1 NMD or in-frame strength
ladder, `guideline_conditions.strength_resolution_source` records that second
document and exact Figure 1 row separately; a Figure 1 locator is never
misrepresented as a locator inside the RNA guidance.

## Quantitative limits and contextual evidence

The generic engine uses only the SVI Figure 1 `>10%` / `<10%` protein-altered
branches where that flowchart explicitly calls for them. It has no universal
pext, tissue-expression, RNA-abundance, residual-normal-transcript, assay
quality, or translation-reinitiation threshold. Those observations may support
the supplied relevance, rescue, RNA-interpretability, or functional-sufficiency
conclusions, or an applicable pinned VCEP rule; they are not independently
scored.

## LoF mechanism cap

The direct `lof_disease_mechanism` value may be:

- `ESTABLISHED`: permits the Figure 1 strength.
- `CAP_STRONG`: SVI Table 1 one-level cap.
- `CAP_MODERATE`: SVI Table 1 two-level cap.
- `NOT_ESTABLISHED`: PVS1 is `NOT_MET`.

If it is omitted, the engine recognizes the SVI 2018 Table 1 combinations:

1. Strong/Definitive validity, at least three prior pathogenic LoF variants
   without PVS1, `>10%` phenotype-associated variants LoF, and distribution
   across more than one exon (unless single-exon gene) → Figure 1 strength.
2. At least Moderate validity, at least two qualifying prior LoF variants plus
   null-model phenocopy → cap at Strong.
3. At least Moderate validity, either of those two supporting findings → cap
   at Moderate.

Partial Table 1 inputs do not prove that LoF is absent. The smallest requested
missing evidence is normally a curated `lof_disease_mechanism` conclusion.

## VCEP specification declaration

The following is a compact, synthetic example. A clinical deployment should
use the actual CSpec document identifier, approved version, source URL, and
its exact documented conditions.

```json
{
  "name": "Example VCEP PVS1 specification",
  "version": "pinned-version",
  "source_url": "https://example.org/pinned-specification",
  "specification_id": "published-document-id-if-available",
  "gene_id": "HGNC:example",
  "disease_id": "MONDO:example",
  "transcript_id": "NM_example.1",
  "mode_of_inheritance": "AD",
  "rules": [{
    "source_locator": "Published PVS1 table, initiation-codon row",
    "pathway": "INITIATION_CODON",
    "when": {
      "alternative_start_codon": false,
      "upstream_pathogenic_variant_evidence": false
    },
    "decision": {
      "status": "ACTIVATED",
      "strength": "VERY_STRONG",
      "decision_path": "EXAMPLE_VCEP_START_LOSS"
    }
  }]
}
```

Matching is exact and deterministic. Within one specification, a matching
path-specific rule takes precedence over that specification's `ANY` fallback.
Across separately supplied specifications, two matches, a conflicting source,
or a potentially applicable declaration whose scope/conditions are unknown
yields `UNDETERMINED`; the generic SVI rule is not used to bypass it.

VCEP specifications are caller-supplied, versioned snapshots, not an embedded
or live VCEP catalogue. Their `when` clauses use a deliberately small,
non-executable condition language: exact values/alternatives, numeric
comparisons and inclusive ranges, and explicit collection membership tests for
affected exons or splice outcomes. This permits published, scoped
gene/disease-specific refinements without claiming universal VCEP coverage.
Unordered set operands are canonicalized before matching and audit output, so
the same pinned declaration produces reproducible payloads across processes.
When generic pathway identification is incomplete, an `ANY` rule remains
eligible; a path-specific VCEP exception must carry and match an explicit
pathway-compatible anchor in its conditions: `variant_consequence`, a
CNV/whole-gene fact, a start-codon fact, a canonical/non-canonical
`splice_site_position`, or `affected_exons` for an exon-level structural path.
This allows a published positional/exon rule to resolve a missing generic
class while preventing an unrelated rule from taking over the case.

## Status mapping

- `ACTIVATED`: all requirements of one branch are established; strength is the
  strongest result justified after applicable caps/exclusions.
- `NOT_MET`: the supported path was assessed and a known biological condition
  precludes PVS1 (for example, LoF mechanism not established, rescue, or
  non-disruptive duplication).
- `UNDETERMINED`: a decision-critical condition is unknown or conflicting.
- `NOT_APPLICABLE`: the consequence is outside the generic PVS1 scope (for
  example, missense) or a non-canonical splice variant lacks RNA/VCEP support.

All non-activated decisions have strength `NONE`; `NONE` never means an
unknown result by itself.
