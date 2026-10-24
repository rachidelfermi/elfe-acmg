# Adaptive PS1 decision table

This document is the versioned source-to-rule map implemented by
`adaptive_ps1`. It is deliberately local and static: evaluation does not query
ClinVar, the ClinGen registry, transcript databases, or splice predictors.
All biological assertions are inputs supplied by an upstream annotation and
curation workflow.

## Implemented generic rules

| Engine path | Source | Minimum sufficient evidence | Material exclusions / unresolved evidence | Permitted result |
|---|---|---|---|---|
| `SAME_AMINO_ACID_CHANGE` | Richards et al., ACMG/AMP 2015, Table 3, PS1 ([PMID 25741868](https://pubmed.ncbi.nlm.nih.gov/25741868/)) | A candidate and independently classified **Pathogenic** reference have a curated established classification record, a validated identical amino-acid substitution, distinct alleles, the same verified protein reference or a curated equivalence assertion, compatible disease context, and documented no-splice-impact conclusions (or a curated compatibility conclusion without contrary evidence). | Same allele; different protein effects; transcript/protein-reference mismatch without established equivalence; VUS/benign/LP generic reference; circular or conflicting reference classification; missing classification independence; predicted/observed splice mechanism or missing required splice assessment. | `ACTIVATED / STRONG` only. |
| `SPLICING_SAME_PREDICTED_EVENT` | Walker et al., ClinGen SVI Splicing Subgroup 2023, Table 2 ([PMID 37352859](https://pubmed.ncbi.nlm.nih.gov/37352859/)) | Candidate and P/LP reference are distinct alleles on the same curation transcript (or have a curated transcript-equivalence assertion), affect the same donor/acceptor motif, have precisely the same predicted splice event, the candidate prediction is curated as similar to or stronger than the reference prediction, and a Table 2 baseline/position row is fully established. The reference pathogenic mechanism must be established as splicing. Encoded exonic protein effects must be reconciled when relevant. | Different event/motif; unknown or conflicting transcript, event, position, baseline, encoded protein-effect, reference mechanism, or prediction-strength relation; a reliable RNA result demonstrating normal splicing, a different product, or a different residual-normal-splicing profile; reference circularity/conflict. | Exactly the strength in Table 2 below. |
| `VCEP_*` | Caller-supplied pinned ClinGen VCEP configuration | Matching gene/disease (and optional inheritance/transcript) scope, explicit reference-classification condition, exact declarative rule, independent P/LP reference preconditions, and the non-overridable same-AA or same-splice-event safety checks for its declared pathway. | Unknown potentially applicable VCEP scope/rule blocks fallback; multiple matching VCEP rules conflict; an unresolved or contradictory core comparison cannot be bypassed by configuration. | Only the strength explicitly declared in the pinned configuration. |

The standard path intentionally does not treat a missing splice field as
normal splicing. RNA experiments are optional only after a non-splice-impact
conclusion has otherwise been established; their absence never provides that
conclusion.

## Walker et al. 2023 Table 2 mapping

All rows additionally require precisely matched predicted splice events, a
candidate prediction of similar or greater strength, and a common curation
transcript or a curated `splice_transcript_equivalence` assertion. The engine
accepts the curated relation `splice_prediction_strength_comparable`; it does
not invent a numeric SpliceAI delta. Available reliable RNA evidence is checked
for contradictions but is not required to activate the predicted-event path.

| Candidate location | Candidate baseline supplied to engine | Reference relation | P reference | LP reference |
|---|---|---|---|---|
| Outside canonical ±1,2 but within the same motif | `PP3` | Same nucleotide | Strong | Moderate |
| Outside canonical ±1,2 but within the same motif | `PP3` | Different nucleotide in same donor/acceptor motif | Moderate | Supporting |
| Canonical ±1,2 | `PVS1` | Within the same donor/acceptor ±1,2 dinucleotide | Supporting | Not applicable |
| Canonical ±1,2 | `PVS1` | Reference in same motif outside canonical ±1,2 | Supporting | Supporting |
| Canonical ±1,2 | `PVS1_STRONG`, `PVS1_MODERATE`, or `PVS1_SUPPORTING` | Within the same donor/acceptor ±1,2 dinucleotide | Strong | Not applicable |
| Canonical ±1,2 | `PVS1_STRONG`, `PVS1_MODERATE`, or `PVS1_SUPPORTING` | Reference in same motif outside canonical ±1,2 | Moderate | Supporting |

`applicable_splicing_pvs1_strength` is a pre-curated input. The PS1 package
does not invoke or duplicate PVS1 logic. RNA splicing results are treated as
potentially material corroboration or contradiction; they are not re-labelled
as PS3 by this engine.

## VCEP configuration contract

No comprehensive VCEP catalogue is bundled. To apply a VCEP modification,
provide a reviewed immutable snapshot such as:

```python
{
    "name": "Example VCEP PS1 specification",
    "version": "1.2.0",
    "source_url": "https://cspec.genome.network/cspec/ui/svi/doc/GN000?version=1.2.0",
    "specification_id": "GN000-PS1-1.2.0",
    "scope": {"gene_id": "GENE", "disease_id": "MONDO:0000000"},
    "rules": [{
        "source_locator": "PS1 section",
        "pathway": "SAME_AMINO_ACID_CHANGE",
        "when": {
            "protein_substitution_equivalence": True,
            "candidate_splice_impact": "NO_IMPACT",
        },
        "reference_when": {
            "classification": "LIKELY_PATHOGENIC",
            "classification_evidence": "ESTABLISHED",
        },
        "unless": {"protein_effect_mechanism_compatibility": False},
        "decision": {
            "status": "ACTIVATED",
            "strength": "MODERATE",
            "decision_path": "EXAMPLE_VCEP_PS1_MODERATE",
        },
    }],
}
```

This is configuration syntax, not a scientific claim about an actual VCEP.
The caller is responsible for validating the declared source and translating
the source wording faithfully. The engine refuses a VCEP configuration without
gene/disease scope, version, source URL, an explicit reference-classification
condition for activation, and explicit eligibility conditions. `pathway` must
be `SAME_AMINO_ACID_CHANGE` or `SPLICING_SAME_PREDICTED_EVENT`; use
`decision_path` for a VCEP-specific label. A VCEP rule can set its documented
strength or LP/P requirement, but cannot bypass a different protein effect,
same allele, unresolved transcript mapping, non-comparable splice event, or a
known RNA contradiction. If more than one supplied specification version
matches the candidate scope, `vcep_specification_version` must select one;
otherwise the result is undetermined. Likewise, a potentially applicable but
unresolved scoped rule blocks a favorable VCEP rule because this compact
configuration contract has no implicit rule-priority system. A matched
`unless` or `reference_unless` clause is treated as an explicit VCEP PS1
exclusion and cannot fall through to a generic path.

Rule conditions support scalar equality, a membership form such as
`{"in": ["PATHOGENIC", "LIKELY_PATHOGENIC"]}`, and numeric comparison
operators (for example, `{"<=": 0.2}` or `{">=": 0.8}`). Numeric conditions
are executed only when they appear in the caller-pinned VCEP snapshot; the
engine neither derives nor applies a generic predictor cutoff. For a VCEP
splice rule, only a matched candidate prediction-magnitude bound can replace
the generic candidate-versus-reference strength assertion; a reference-only
score bound remains an additional condition, not evidence that the candidate
prediction is sufficiently strong.

## Supported input assertions

Candidate/reference identifiers, genome/coding/protein annotations,
transcript and protein-accession versions, amino-acid components, reference
classification provenance, independence, disease association, molecular
mechanism, splice motif/event/prediction/RNA evidence, transcript/event/motif
pairwise equivalence, splice-product functional consequence, and VCEP metadata
are typed fields in `BiologicalFeatures`. Pairwise assertions can be supplied
inside an individual `references[]` entry so multiple comparison variants do
not share an accidental global assertion. Each accepts either a scalar
shorthand or an evidence envelope:

```python
{
    "state": "TRUE",
    "value": "PATHOGENIC",
    "provenance": [{
        "source": "approved VCEP curation record",
        "kind": "CURATED",
        "identifier": "record-123",
    }],
}
```

Boolean evidence uses all four states: `TRUE`, `FALSE`, `UNKNOWN`, and
`CONFLICTING`. Categorical values use typed enums plus `UNKNOWN` or
`CONFLICTING` evidence states. Missing data remains unknown and is never
converted to a favorable finding. An explicitly supplied reference-level
`UNKNOWN` pairwise assertion overrides a legacy top-level default for that
comparison rather than inheriting a favorable assertion intended for another
reference.

## Deliberate limits

- The package does not validate HGVS syntax, normalize genome alleles, align
  transcripts, query databases, or calculate splice prediction scores.
- It does not decide PVS1, PP3, PS3, or an overall ACMG classification.
- It does not turn reference review stars, a raw SpliceAI magnitude, or a
  count of matching references into an invented PS1 threshold.
- A splice PS1 comparison remains undetermined if the reference pathogenic
  mechanism has not been established as splicing. A known different mechanism
  excludes that reference from both generic and VCEP splice paths.
- It does not bundle unverified or incomplete VCEP specifications. A supplied
  VCEP snapshot is applied only within its exact declared scope and version.
- No VCEP specification is bundled in this release: every VCEP remains
  unsupported until a caller supplies its reviewed, immutable snapshot. The
  standalone `vcep_gene_scope`, `vcep_ps1_eligibility_rules`, and related
  inventory fields are retained as typed intake metadata, but cannot activate
  PS1 without the complete `vcep_specifications` configuration and will block
  unsafe generic fallback when they imply an incomplete VCEP declaration.
- The ANO5 two-LP-reference exception and other multi-reference VCEP
  exceptions are not inferred by aggregation. They require an explicit pinned
  VCEP configuration rule.

