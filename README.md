# Adaptive ACMG/ClinGen evidence engines

The installable `adaptive-acmg-evidence` distribution has one primary `acmg`
package. Scientific evaluators live under `acmg.criteria`; the earlier
top-level package names remain thin compatibility imports. Inputs are already
extracted biological evidence, supplied separately for each criterion.

```python
from acmg import classify_variant
from acmg.criteria.ps1 import evaluate_ps1
from acmg.criteria.ps2_pm6 import evaluate_de_novo
from acmg.criteria.pm3 import evaluate_pm3

result = classify_variant({"variant_id": "example", "criteria": {}})
```

The classifier uses the shared `ps2_pm6` evaluator for mutually exclusive PS2
and PM6 assignment. `acmg.criteria.pm6.evaluate_pm6` retains the older
standalone PM6 interface and its separately verified VCEP behavior for existing
callers. It is not used by the unified classifier. PP5 and BP6 remain inactive.
See [package architecture](docs/package-architecture.md) for the canonical
paths, compatibility imports, and known limits.
The [evidence feature contract](docs/feature-contract.md) maps shared facts,
local aliases, minimum-evidence examples, and dependency handling.

The [unified scientific validation audit](docs/unified-scientific-validation.md)
records criterion coverage, source checks, combination tests, dependency
limits, and the remaining clinical-validation work.

## Unified classification

`acmg.classifier.classify_variant` runs the existing criterion evaluators over
already extracted evidence and applies ACMG/AMP 2015 Table 5 with the ClinGen
SVI PM2 Supporting update. Give each criterion only its native evaluator input;
missing criterion inputs remain undetermined.

```python
from acmg import classify_variant

result = classify_variant({
    "variant_id": "example",
    "criteria": {
        "PVS1": {
            "gene_id": "GENE-TEST", "disease_id": "DISEASE-TEST",
            "transcript_id": "TX-TEST.1", "variant_consequence": "frameshift_variant",
            "lof_disease_mechanism": "established", "gene_disease_validity": "definitive",
            "biologically_relevant_transcript": True, "nmd_prediction": True,
        }
    },
})
```

Use the single `PS2_PM6` input for de novo observations. `PP5` and `BP6` are
audited but inactive. An explicitly requested VCEP-wide or Bayesian
classification remains unsupported until complete verified combining and
prohibition rules are registered. See the
[source-to-code matrix](docs/unified-classification-source-matrix.md).

## BP6 external benign classifications

The standalone [BP6 evaluator](adaptive_bp6/README.md) preserves external
classifications as context. Under current ClinGen guidance it always returns
`NOT_APPLICABLE` with `NONE` strength. Its historical mode documents the 2015
rule without assigning BP6 evidence.

## BP3 in-frame repeat indels

The standalone [BP3 evaluator](adaptive_bp3/README.md) applies the general
ACMG/AMP rule and pinned ClinGen VCEP profiles for FOXG1, RYR1, BMPR2, and
MYOC. It accepts extracted protein and repeat evidence, returns a structured
decision, and runs offline.

## BP1 missense in a primarily truncating disease context

`adaptive_bp1` evaluates supplied biological features offline. Its general ACMG/AMP
pathway requires a relevant missense consequence, valid gene–disease relationship,
a cited curated assessment that truncating pathogenic variants predominate,
an assessment that missense is not a common disease mechanism, and a completed
exception assessment. A separate curated spectrum pathway requires disease-matched,
independent classifications and addressed ascertainment bias. The package also
pins three ClinGen VCEP rules: DICER1 and RYR1 BP1 exclusions, and the MAP2K1
RASopathy truncating-variant rule. See the [BP1 source-to-rule matrix](docs/bp1-guideline-matrix.md).

```python
from adaptive_bp1 import evaluate_bp1

# Synthetic fixture only; TEST1 is not a clinical gene assertion.
result = evaluate_bp1({
    "gene_id": "TEST1", "disease_id": "SYNTHETIC:1",
    "variant_consequence": "MISSENSE", "biologically_relevant_transcript": True,
    "gene_disease_validity": True, "truncating_predominant": True,
    "missense_common_mechanism": False, "exceptions_assessed": True,
    "mechanism_evidence_quality": "CURATED",
    "mechanism_evidence_source": "synthetic curated assessment",
})
print(result.to_dict())  # ACTIVATED, SUPPORTING
```

Inputs may use explicit evidence wrappers such as
`{"state": "KNOWN", "value": true, "source": "...", "version": "..."}`.
Other states are `UNKNOWN`, `CONFLICTING`, and `NOT_APPLICABLE`.

## BS4 family nonsegregation

`adaptive_bs4` evaluates supplied family observations or a validated BRCA2 likelihood ratio offline. It includes pinned general ACMG, CTLA4, PIK3R1, and BRCA2 BS4 rules. See the [BS4 API](adaptive_bs4/README.md) and [source-to-rule matrix](docs/bs4-guideline-matrix.md).

## BS2 healthy-individual genotype evidence

`adaptive_bs2` evaluates supplied genotype and clinical evidence offline. It implements the 2015 general BS2 rule, the pinned Rett/Angelman-like VCEP count rule, and the RUNX1 FPD/AML exclusion. See the [BS2 source-to-rule matrix](docs/bs2-guideline-matrix.md) for exact scope and limitations.

```python
from adaptive_bs2 import evaluate_bs2

# Synthetic example; TEST1 is not a curated clinical gene.
result = evaluate_bs2({
    "variant_id": "synthetic-variant", "variant_reliable": True,
    "disease": {
        "gene_id": "TEST1", "disease_id": "SYNTHETIC",
        "gene_disease_validity": True, "inheritance_mode": "AR",
        "molecular_mechanism_compatible": True,
        "full_penetrance": True, "early_onset": True,
    },
    "observations": [{
        "individual_id": "person-1", "genotype": "HOMOZYGOUS",
        "genotype_reliable": True, "confirmed_unaffected": True,
        "adult": True, "beyond_expected_onset": True,
        "clinical_assessment_adequate": True,
    }],
})
print(result.to_dict())  # ACTIVATED, STRONG
```

Any evidence value may instead be supplied as `{"state": "KNOWN", "value": true, "source": "..."}`. Other states are `UNKNOWN`, `CONFLICTING`, and `NOT_APPLICABLE`. Missing essential evidence yields `UNDETERMINED`; known failed requirements yield `NOT_MET`. Optional biological exclusions only affect a result when supplied as present or conflicting. Population database membership alone never establishes `confirmed_unaffected`.

## BS1 population-frequency criterion

`adaptive_bs1` evaluates supplied, validated population frequencies against
versioned ClinGen TPM1 and RMRP BS1 specifications. It returns an auditable
BS1 status and strength. See the [BS1 source-to-rule matrix and API](adaptive_bs1/README.md).

## PP5 reputable source criterion

`adaptive_pp5` accepts already extracted external classifications and returns an
auditable PP5 decision. The default follows the ClinGen SVI 2018 recommendation
to discontinue PP5. See [PP5 policy and sources](docs/pp5-policy.md).

```python
from adaptive_pp5 import evaluate_pp5

result = evaluate_pp5({
    "external_classification": "PATHOGENIC",
    "classification_source": "ClinVar",
    "classification_evidence_available": True,
})
print(result.to_dict())  # NOT_APPLICABLE, NONE
```

The `HISTORICAL_ACMG_2015` selector documents the original definition but does
not award PP5 without verified source-specific conditions.

## PP3 computational evidence evaluator

`adaptive_pp3` evaluates supplied REVEL and SpliceAI scores offline using pinned ClinGen SVI rules and RUNX1 GN008 v3.1. The [PP3 rule matrix](docs/pp3-guideline-matrix.md) lists exact sources, thresholds, scope, and limitations.

```python
from adaptive_pp3 import evaluate_pp3

result = evaluate_pp3({
    "variant_id": "example", "variant_consequence": "MISSENSE",
    "transcript_id_version": "NM_EXAMPLE.1",
    "biologically_relevant_transcript": True,
    "selected_missense_predictor": "REVEL",
    "predictions": [{"name": "REVEL", "score": 0.80, "version": "1",
                     "calibration": "PEJAVER_2022",
                     "variant_id": "example", "transcript_id_version": "NM_EXAMPLE.1"}],
})
print(result.to_dict())  # PP3 ACTIVATED, MODERATE
```

The evaluator emits one PP3 decision with used features, missing critical features, exclusions, source rule, and BP4 candidacy. It does not apply BP4 or combine ACMG criteria.

## PP2 missense mechanism evaluator

`adaptive_pp2` evaluates supplied evidence offline. It implements qualitative ACMG/AMP 2015 PP2 at Supporting strength and pinned ClinGen ACTA1 v2.0 and DICER1 v1.4 specifications. See the [PP2 source-to-rule matrix](docs/pp2-guideline-matrix.md) for scope, sources, and limitations.

```python
from adaptive_pp2 import evaluate_pp2

# Synthetic example; TEST1 is not a real curated gene.
result = evaluate_pp2({
    "gene_id": "TEST1", "disease_id": "MONDO:synthetic",
    "variant_consequence": "MISSENSE", "relevant_transcript": True,
    "gene_disease_valid": True,
    "missense_mechanism_common": {
        "value": True, "source": "synthetic expert review", "quality": "HIGH",
        "independent_of_pp2": True, "disease_match": True,
    },
    "benign_missense_uncommon": {
        "value": True, "source": "synthetic expert review", "quality": "HIGH",
        "independent_of_pp2": True, "disease_match": True,
    },
})
print(result.to_dict())  # ACTIVATED, SUPPORTING
```

Assessments accept `KNOWN`, `UNKNOWN`, `CONFLICTING`, and `NOT_APPLICABLE` states. A missing primary assessment returns `UNDETERMINED`; reliable contradictory evidence returns `NOT_MET`. Counts and missense constraint are audit context, never stand-alone activation or strength upgrades.

## PP1 family segregation evaluator

`adaptive_pp1` evaluates curated PP1 evidence offline for the pinned ClinGen Hearing Loss, RUNX1 and PTEN scopes. See the [PP1 rule matrix](docs/pp1-guideline-matrix.md) for source versions, evidence requirements and limitations. The package does not infer genotypes or pedigree relationships.

```python
from adaptive_pp1 import evaluate_pp1

result = evaluate_pp1({
    "variant_id": "example-variant", "gene_id": "RUNX1", "disease_id": "FPD_AML",
    "gene_disease_validity": "DEFINITIVE", "inheritance_mode": "AD",
    "variant_reliable": True, "disease_mechanism_compatible": True,
    "families": [{"family_id": "F1", "affected_segregations": 3,
                  "count_validated": True, "proband_excluded": True,
                  "relationships_confirmed": True, "genotype_reliable": True,
                  "phenotype_relevant": True}],
})
print(result.to_dict())  # ACTIVATED, SUPPORTING; illustrative synthetic evidence
```

## PM5 same-residue missense evaluator

`adaptive_pm5` evaluates curated candidate and reference missense substitutions offline. It implements ACMG/AMP 2015 PM5 at Moderate strength and pinned HNF1A v1.1 and ITGA2B/ITGB3 v1 VCEP modifications. The [PM5 guideline matrix](docs/pm5-guideline-matrix.md) gives source links, scope, minimum evidence, and limitations. All package tests use synthetic examples.

```python
from adaptive_pm5 import evaluate_pm5

result = evaluate_pm5({
    "gene_id": "EXAMPLE", "disease_id": "MONDO:synthetic",
    "gene_disease_validity": "ESTABLISHED",
    "candidate": {
        "variant_id": "candidate", "consequence": "MISSENSE",
        "transcript_id_version": "NM_TEST.1", "amino_acid_position": 132,
        "reference_amino_acid": "Arg", "alternate_amino_acid": "His",
        "transcript_relevant": True,
    },
    "references": [{
        "variant_id": "reference", "consequence": "MISSENSE",
        "transcript_id_version": "NM_TEST.1", "amino_acid_position": 132,
        "reference_amino_acid": "Arg", "alternate_amino_acid": "Cys",
        "classification": "PATHOGENIC", "classification_source": "curated dossier",
        "classification_evidence_reviewed": True,
        "classification_independent_of_pm5": True, "disease_match": True,
    }],
})
print(result.to_dict())
```

`KNOWN`, `UNKNOWN`, `CONFLICTING`, and `NOT_APPLICABLE` can be supplied in `evidence_states` for `candidate` and `gene_disease`; omitted values remain `None` in the typed input model. A result has one of `ACTIVATED`, `NOT_MET`, `UNDETERMINED`, or `NOT_APPLICABLE`. Optional provenance never counts as additional PM5 evidence.

## PM4 protein length change evaluator

`adaptive_pm4` evaluates supplied in-frame deletion, in-frame insertion, and stop-loss
evidence offline. It returns a PM4 status and strength, not an overall variant
classification. The [PM4 guideline matrix](docs/pm4-guideline-matrix.md) lists the
implemented sources and scientific limits.

```python
from adaptive_pm4 import evaluate_pm4

result = evaluate_pm4({
    "gene_id": "EXAMPLE", "disease_id": "EXAMPLE_DISEASE",
    "transcript_id_version": "NM_EXAMPLE.1",
    "transcript_biological_relevance": True,
    "variant_consequence": "INFRAME_DELETION",
    "reading_frame_preserved": True,
    "protein_length_difference": -2,
    "repeat_region_status": False,
})
print(result.to_dict())
```

The optional `vcep_specification_version="rett_angelman_2022"` selects the
published Rett/Angelman-like panel PM4 profile for its six scoped genes. Other
VCEP rules are not built in. Known, unknown, conflicting, and not-applicable
evidence can be supplied with `Evidence(EvidenceState.KNOWN, value)` or a
`{"state": "KNOWN", "value": ...}` mapping; omitted fields are unknown. The
caller is responsible for correct disease context and molecular annotations.

## PM2 population rarity evaluator

`adaptive_pm2` evaluates supplied population evidence offline. It applies general ClinGen SVI PM2 at Supporting strength for confirmed absence, and pinned ACADVL v2.2 and F9 v2.1 VCEP rules. It does not query population databases or calculate a general rarity cutoff. See the [PM2 guideline matrix](docs/pm2-guideline-matrix.md) for scope and scientific limitations.

```python
from adaptive_pm2 import evaluate_pm2

result = evaluate_pm2({
    "variant_id": "example", "genome_assembly": "GRCh38",
    "normalized_variant_representation": "NC_000001.11:g.100A>G",
    "gene_id": "EXAMPLE", "disease_id": "MONDO:example",
    "population_records": [{
        "name": "gnomAD", "version": "4.1", "genome_assembly": "GRCh38",
        "normalized_variant_representation": "NC_000001.11:g.100A>G",
        "allele_count": 0, "allele_number": 100000,
        "identity_confirmed": True, "filter_pass": True,
        "site_callable": True, "technically_reliable": True,
    }],
})
print(result.to_dict())
```

## PS2/PM6 de novo evaluator

`ps2pm6` evaluates already-extracted de novo observations under ClinGen SVI
PS2/PM6 v1.1. It pools only eligible independent cases, applies the published
2/1, 1/0.5, and 0.5/0.25 point table (including the heterogeneous-phenotype
cap), and emits exactly one criterion: PS2 when any pooled observation has
confirmed parentage, otherwise PM6. Missing parental calls are never treated
as negative calls. See the [PS2/PM6 decision table](docs/ps2-pm6-decision-table.md)
for the versioned sources, schema, VCEP declaration format, and examples.

```python
from ps2pm6 import evaluate_ps2_pm6

decision = evaluate_ps2_pm6({
    "proband_variant_present": True,
    "maternal_genotype": "ABSENT",
    "paternal_genotype": "ABSENT",
    "variant_call_quality": "ADEQUATE",
    "maternity_confirmed": True,
    "paternity_confirmed": True,
    "gene_disease_validity": "DEFINITIVE",
    "inheritance_pattern": "AUTOSOMAL_DOMINANT",
    "phenotype_category": "HIGHLY_SPECIFIC",
})
```

## PM1 regional evidence evaluator

`adaptive_pm1` evaluates supplied, already-curated protein hotspot/domain evidence and three pinned ClinGen VCEP specifications. It operates offline and does not annotate variants or classify overall pathogenicity.

```python
from adaptive_pm1 import evaluate_pm1

decision = evaluate_pm1({
    "gene_id": "DICER1",
    "disease_id": "MONDO:0100216",
    "variant_consequence": "missense",
    "transcript_id_version": "NM_177438.2",
    "protein_position": 1705,
})
print(decision.to_dict())  # PM1, ACTIVATED, MODERATE
```

For general ACMG/AMP pathways, pass one or more `Region` records with curated positive/negative/unknown evidence. Omission means unknown; numerical regional counts alone never establish a hotspot. The [PM1 source-to-rule matrix](docs/pm1-guideline-matrix.md) documents minimum evidence, scope, and scientific limits.

## PS4 association evaluator

`adaptive_ps4` evaluates already-curated PS4 evidence without network access. It
supports ACMG/AMP 2015 case-control OR/RR evidence and the scoped ClinGen
Glaucoma VCEP MYOC v2.1 proband-count rule. It does not infer missing evidence,
pool studies, or classify a variant. See the [PS4 guideline matrix](docs/ps4-guideline-matrix.md).

```python
from adaptive_ps4 import evaluate_ps4

result = evaluate_ps4({
    "variant_id": "example", "gene_id": "GENE", "disease_id": "DISEASE",
    "gene_disease_validity": True,
    "studies": [{
        "study_id": "study-1", "variant_id": "example", "disease_id": "DISEASE",
        "case_carriers": 30, "case_noncarriers": 70,
        "control_carriers": 1, "control_noncarriers": 999,
        "cases_affected": True, "controls_valid": True,
        "variant_calls_valid": True, "technical_comparable": True,
        "independent_cases": True, "overlapping_cohort": False,
        "ancestry_comparable": True,
    }],
})
print(result.to_dict())
```

For published estimates without counts, provide `odds_ratio` or
`relative_risk`, `confidence_interval`, `confidence_level: 0.95`,
`statistic_source`, and `statistic_method` alongside study-validity evidence.
Omitted evidence is unknown; envelopes such as `{"state": "CONFLICTING"}`
are also accepted. MYOC proband counting requires `pm2_met`,
`studies_independent`, and individually identified probands with study and
family IDs. The returned `status` is `ACTIVATED`, `NOT_MET`,
`UNDETERMINED`, or `NOT_APPLICABLE`.

## PS3 functional-evidence evaluator

`adaptive_ps3` is a standalone, deterministic evaluator for the ACMG/AMP PS3
criterion. It implements the four-step ClinGen SVI functional-evidence
framework from Brnich et al. (2019): disease mechanism, assay-class
applicability, validation of the specific assay instance, and interpretation
of the tested variant. It does not annotate variants, download databases,
classify a variant, or activate BS3.

```python
from adaptive_ps3 import evaluate_ps3

decision = evaluate_ps3({
    "variant_id": "example",
    "gene_id": "GENE",
    "disease_id": "DISEASE",
    "disease_mechanism": "loss_of_function",
    "functional_assays": [{
        "assay_name": "example validated assay",
        "assay_measures_variant_effect": True,
        "assay_measures_disease_mechanism": True,
        "variant_result_classification": "abnormal",
        "positive_control_count": 1,
        "negative_control_count": 1,
        "technical_replicates": 3,
        "known_pathogenic_control_results": ["P1", "P2", "P3", "P4", "P5"],
        "known_benign_control_results": ["B1", "B2", "B3", "B4", "B5", "B6"],
        "assay_validation_sample_size": 11,
        "control_variant_classification_quality": True,
        "control_classification_independence": True,
    }],
})

print(decision.to_dict())
```

Inputs may use evidence envelopes such as
`{"state": "UNKNOWN"}` or `{"state": "KNOWN", "value": true,
"provenance": [...]}`. Omission means unknown, never false.

Formal calibration uses the exact pathogenic OddsPath cutoffs published in
Brnich et al. Table 3: `>2.1` Supporting, `>4.3` Moderate, `>18.7` Strong,
and `>350` Very Strong. Without formal statistics, mixed independently
classified pathogenic and benign validation variants permit Supporting with
10 or fewer controls and Moderate with at least 11. Assay type, effect size,
sensitivity, specificity, or statistical significance never creates a
strength by itself.

VCEP rules are accepted as pinned, versioned, gene/disease-scoped
configuration. This package deliberately does not freeze the changing CSpec
Registry into source code. See [the PS3 decision table](docs/ps3-decision-table.md)
for the exact implemented branches and sources.

Run all tests with:

```text
python -m pytest -q
```

## PVS1 evaluator

This package evaluates the ACMG/AMP PVS1 criterion only. It is deliberately
not a complete variant classifier and never combines PVS1 with other ACMG/AMP
criteria.

The generic decision tree implements the ClinGen SVI PVS1 recommendations
published by Abou Tayoun et al. (2018; PMID 30192042), with the original
ACMG/AMP PVS1 eligibility requirement. Gene/disease-specific ClinGen VCEP
rules are supplied by the caller as versioned, scoped declarations; the
package does not embed a mutable or incomplete VCEP database.

See [the decision table](docs/decision-table.md) for sources, branch scope,
and the minimum decision-critical evidence for every generic pathway.

## Quick start

```python
from pvs1 import evaluate_pvs1

decision = evaluate_pvs1({
    "variant_consequence": "frameshift_variant",
    "lof_disease_mechanism": "established",
    "biologically_relevant_transcript": True,
    "nmd_prediction": True,
})

assert decision.status.value == "ACTIVATED"
assert decision.strength.value == "VERY_STRONG"
print(decision.to_dict())
```

All values can instead be supplied as an evidence envelope with a state and
provenance:

```python
{
    "value": True,
    "state": "known",
    "provenance": [{
        "source": "curated transcript assessment",
        "kind": "CURATED",
        "identifier": "internal-record-42"
    }]
}
```

Omitted values are `UNKNOWN`, never `False`. The evaluator returns
`UNDETERMINED` rather than treating them as negative evidence when they are
needed for the active pathway.

## Validation and test suite

The implementation has no runtime dependencies. With Python 3.12+ installed:

```text
python -m unittest discover -s tests -t . -v
```

The engine is decision support software, not a substitute for expert clinical
review. A production deployment should pin and govern its VCEP specifications,
input provenance, release process, and clinical validation.

## PM3 evaluator

`adaptive_pm3` evaluates supplied affected-person observations under ClinGen
SVI PM3 v1.0 or the explicitly selected GAA/Pompe VCEP v2. It has no network
calls or feature extraction. See [the PM3 source-to-rule map](adaptive_pm3/SOURCES.md).

```python
from adaptive_pm3 import evaluate_pm3

result = evaluate_pm3({
    "candidate_variant_id": "candidate-1",
    "gene_id": "PAH", "disease_id": "PKU",
    "inheritance_mode": "AUTOSOMAL_RECESSIVE",
    "candidate_rare_for_disease": True,
    "cases": [{
        "case_id": "proband-1", "affected": True,
        "clinical_disease_match": True,
        "candidate_genotype": "HET", "candidate_call_valid": True,
        "second_variant": {
            "variant_id": "partner-1", "gene_id": "PAH", "genotype": "HET",
            "classification": "PATHOGENIC", "classification_credible": True,
            "classification_independent": True, "disease_match": True,
            "rare_for_disease": True, "call_valid": True,
        },
        "phase": "TRANS", "phase_method": "molecular",
    }],
})
assert result["status"] == "ACTIVATED"
assert result["strength"] == "MODERATE"
```

For GAA/Pompe, set `specification="GAA_POMPE_V2"`, use `disease_id="Pompe
disease"`, and attest `diagnosis_confirmed` and
`second_variant.vcep_classification_confirmed`. Missing required attestations
produce `UNDETERMINED`. `None` means unknown; a supported exclusion uses
`False` or its explicit exclusion field. Classification credibility and
independence must be curated upstream. The evaluator cannot verify ClinVar
review status, population thresholds, or genotype assay quality itself.
For the GAA pseudodeficiency context, attest `candidate_pm2_supporting`,
`functional_support`, and `clinical_laboratory_support`.

The tests in `tests/test_*pm3*` and the six PM3 pathway test modules use
synthetic cases; they are not published case-level validation examples.

## PS1 evaluator

`adaptive_ps1` is a standalone, deterministic evaluator for ACMG/AMP PS1. It
accepts already-curated annotations and does not perform variant annotation,
database lookup, transcript alignment, splice prediction, PVS1 evaluation, or
an overall ACMG classification. The generic implementation covers ACMG/AMP
2015 same-amino-acid PS1 and the ClinGen SVI 2023 same-predicted-splice-event
matrix. It applies VCEP differences only from caller-supplied, pinned
gene/disease/version-scoped configuration. It requires a common curated
transcript (or an explicit mapping assertion) for splice-event comparisons and
does not let a VCEP configuration bypass core allele, protein-effect,
splice-event, or known-RNA-contradiction checks.

```python
from adaptive_ps1 import evaluate_ps1

result = evaluate_ps1({
    "gene_id": "EXAMPLE", "disease_id": "MONDO:example",
    "candidate_protein_hgvs": "p.Val12Leu",
    "candidate_genomic_allele": "NC_000001.11:g.100G>C",
    "candidate_transcript_id_version": "NM_000001.2",
    "candidate_protein_accession_version": "NP_000001.1",
    "candidate_splice_impact": "NO_IMPACT",
    "references": [{
        "variant_id": "reference-1",
        "genomic_allele": "NC_000001.11:g.100G>T",
        "protein_hgvs": "p.Val12Leu",
        "transcript_id_version": "NM_000001.2",
        "protein_accession_version": "NP_000001.1",
        "splice_impact": "NO_IMPACT",
        "classification": "PATHOGENIC",
        "classification_evidence": "ESTABLISHED",
        "classification_independence": True,
        "disease_association": True,
        "molecular_mechanism": "PROTEIN_ALTERATION",
    }],
})

assert result.status.value == "ACTIVATED"
assert result.strength.value == "STRONG"
```

Every decision contains the source locator, features used, critical missing
fields, exclusions, reference-level audit results, and input provenance. See
the [PS1 source-to-rule matrix](docs/ps1-decision-table.md) for exact rule
coverage, the six Walker et al. splice rows, VCEP configuration syntax, and
scientific limitations.
# BA1 evaluator

The standalone BA1 implementation, usage example, scientific source-to-rule
matrix, and limitations are in [adaptive_ba1/README.md](adaptive_ba1/README.md).

# Adaptive BS3 decision engine

`adaptive_bs3.evaluate_bs3` evaluates already extracted functional evidence for BS3 only. It returns a typed `BS3Result`; call `to_dict()` for a JSON-ready mapping. The evaluator runs offline and never estimates assay metrics.

```python
from adaptive_bs3 import evaluate_bs3

result = evaluate_bs3({
    "context": {"gene_id": "GENE", "disease_id": "DISEASE", "disease_mechanism": "LOSS_OF_FUNCTION"},
    "assays": [{
        "assay_id": "study-1", "assay_type": "ENZYME",
        "mechanism_compatible": True, "variant_class_applicable": True,
        "detects_relevant_defect": True, "validation_status": "FORMALLY_CALIBRATED",
        "result_class": "NORMAL", "accepted_normal_interpretation": True,
        "calibration_validated": True, "odds_path": 0.04,
    }],
})
print(result.to_dict())
```

For quantitative results, supply `result_value` and `validated_normal_range`. A categorical result instead needs `accepted_normal_interpretation=True`. For the control-count route, supply wild-type and independently classified benign/pathogenic controls, separation and replicate consistency. Optional reliability fields can exclude activation when they reveal a problem. Any evidence field accepts `{"state": "UNKNOWN"}`, `CONFLICTING`, or `NOT_APPLICABLE` instead of a value.

Direct RNA-splicing results are routed to BP7(RNA) under the 2023 SVI guidance and never activate BS3 by themselves. See the [guideline matrix](docs/bs3-guideline-matrix.md) for source-linked pathways and limits.

# Adaptive BP4 evaluator

`adaptive_bp4` is an offline evaluator for supplied missense and noncanonical splice predictions. It is independent of the other criterion packages in this repository. The [source-to-rule matrix](docs/BP4_SOURCE_RULE_MATRIX.md) records the exact rules and scopes.

```python
from adaptive_bp4 import BiologicalFeatures, Prediction, State, evaluate_bp4

verified = dict(variant_match=State.KNOWN, transcript_match=State.KNOWN,
                genome_build_match=State.KNOWN, model_applicability=State.KNOWN,
                quality=State.KNOWN)
features = BiologicalFeatures(
    variant_id="example", gene_id="GENE", consequence="missense",
    transcript_id_version="NM_123.1", biologically_relevant_transcript=State.KNOWN,
    predictions=(Prediction("REVEL", 0.18, **verified),),
)
result = evaluate_bp4(features)
print(result.status, result.strength, result.reason)
```

`KNOWN` mapping checks are caller assertions; omitted checks remain `UNKNOWN` and prevent activation. For splice-only assessment, set `splice_region_eligible=True` only after reviewing the site, and `protein_effect_relevant=False` only after excluding a relevant protein effect. RPE65 VCEP 1.0 is selected only for scoped gene, disease, and transcript inputs. Unsupported predictors are ignored; several benign scores are never treated as a vote.

Run focused tests with `python -m pytest tests/test_bp4_*.py -q`.

## BP2

The standalone `adaptive_bp2` package evaluates curated cis/trans observations under the base ACMG/ClinGen BP2 rule and source-pinned DICER1 and BMPR2 VCEP exceptions. See the [BP2 source matrix, evidence contract, and example](adaptive_bp2/README.md). Run focused tests with `python -m pytest tests/test_bp2.py -q`.
# BP5

The standalone, offline BP5 evaluator is documented in [adaptive_bp5/README.md](adaptive_bp5/README.md). Import it with `from adaptive_bp5 import evaluate_bp5`.

# Adaptive BP7

`adaptive_bp7` evaluates only the BP7 criterion from supplied transcript-specific evidence. It defaults to the ClinGen SVI 2023 splicing recommendation and also includes pinned ACMG/AMP 2015 and TP53 VCEP v2.1.0 paths. It performs no annotation, scoring, or database access.

```python
from adaptive_bp7 import evaluate_bp7

result = evaluate_bp7({
    "gene_id": "GENE", "genome_assembly": "GRCh38",
    "transcript_id_version": "NM_EXAMPLE.1",
    "variant_consequence": "synonymous", "consequence_confident": True,
    "biologically_relevant_transcript": True, "variant_mapping_valid": True,
    "exon_core_splice_motif": False,
    "prediction": {
        "predictor": "SpliceAI", "score": 0.05,
        "transcript_id_version": "NM_EXAMPLE.1", "genome_assembly": "GRCh38",
        "variant_match": True, "applicable": True, "calibrated": True,
        "quality": "KNOWN",
    },
})
assert result.status == "ACTIVATED"
```

`ACTIVATED` means the selected BP7 rule is met, not that the variant is benign or likely benign overall. Missing essential evidence returns `UNDETERMINED`. Select `guideline="ACMG_2015"` to require a reliable upstream conservation assessment. Select `guideline="TP53_2.1.0"` only for TP53/Li-Fraumeni syndrome. See [the BP7 guideline matrix](docs/bp7-guideline-matrix.md) for source links and scope.

