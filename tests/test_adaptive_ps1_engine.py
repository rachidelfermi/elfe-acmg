"""Public-API acceptance tests for the standalone adaptive PS1 evaluator.

All examples are synthetic.  They intentionally exercise only the evidence
that a particular ACMG/ClinGen path needs, so optional annotation fields do
not accidentally become requirements.
"""

from __future__ import annotations

from copy import deepcopy

import pytest

from adaptive_ps1 import SchemaValidationError, evaluate_ps1


def pathogenic_reference(**overrides):
    """A credible, independently classified reference for the standard path."""

    reference = {
        "variant_id": "REF-P-1",
        "gene_id": "GENE-TEST",
        "disease_id": "DISEASE-TEST",
        "genomic_allele": "chr1:g.100A>C",
        "coding_hgvs": "c.100A>C",
        "protein_hgvs": "p.Arg34Gln",
        "transcript_id_version": "NM_000001.1",
        "protein_accession_version": "NP_000001.1",
        "classification": "PATHOGENIC",
        "classification_evidence": "ESTABLISHED",
        "classification_independence": True,
        "disease_association": True,
        "molecular_mechanism": "PROTEIN_ALTERATION",
        "splice_prediction": "NO_IMPACT",
    }
    reference.update(overrides)
    return reference


def standard_features(**overrides):
    """Smallest ordinary protein-equivalence fixture that can activate PS1."""

    evidence = {
        "gene_id": "GENE-TEST",
        "disease_id": "DISEASE-TEST",
        "candidate_genomic_allele": "chr1:g.100A>G",
        "candidate_coding_hgvs": "c.100A>G",
        "candidate_protein_hgvs": "p.Arg34Gln",
        "candidate_transcript_id_version": "NM_000001.1",
        "candidate_protein_accession_version": "NP_000001.1",
        "candidate_variant_consequence": "MISSENSE",
        "candidate_splice_prediction": "NO_IMPACT",
        "reference_variants": [pathogenic_reference()],
    }
    evidence.update(overrides)
    return evidence


def table2_features(**overrides):
    """Walker et al. Table 2 non-canonical, same-nucleotide PS1 fixture."""

    evidence = standard_features(
        candidate_genomic_allele="chr1:g.100+5G>A",
        candidate_coding_hgvs="c.100+5G>A",
        candidate_variant_consequence="SPLICE_REGION",
        candidate_splice_prediction="SPLICE_ALTERING",
        candidate_splice_site_type="DONOR",
        candidate_splice_motif_position="WITHIN_MOTIF",
        candidate_splice_motif_id="intron-3-donor",
        candidate_splice_nucleotide_position="c.100+5",
        candidate_spliceai_event_type="EXON_SKIP",
        candidate_spliceai_event_position="exon-3",
        splice_event_equivalence=True,
        splice_prediction_strength_comparable=True,
        applicable_splicing_pvs1_strength="PP3",
        reference_variants=[
            pathogenic_reference(
                variant_id="REF-SPLICE-1",
                genomic_allele="chr1:g.100+5G>C",
                coding_hgvs="c.100+5G>C",
                molecular_mechanism="SPLICING",
                splice_prediction="SPLICE_ALTERING",
                splice_site_type="DONOR",
                splice_motif_position="WITHIN_MOTIF",
                splice_motif_id="intron-3-donor",
                splice_nucleotide_position="c.100+5",
                spliceai_event_type="EXON_SKIP",
                spliceai_event_position="exon-3",
            )
        ],
    )
    evidence.update(overrides)
    return evidence


def canonical_table2_features(*, baseline):
    """Walker Table 2 canonical-dinucleotide fixture for a P reference."""

    evidence = table2_features(
        candidate_genomic_allele="chr1:g.100+1G>A",
        candidate_coding_hgvs="c.100+1G>A",
        candidate_splice_motif_position="CANONICAL_DINUCLEOTIDE",
        candidate_splice_nucleotide_position="c.100+1",
        applicable_splicing_pvs1_strength=baseline,
    )
    evidence["reference_variants"][0].update(
        {
            "genomic_allele": "chr1:g.100+1G>C",
            "coding_hgvs": "c.100+1G>C",
            "splice_motif_position": "CANONICAL_DINUCLEOTIDE",
            "splice_nucleotide_position": "c.100+1",
        }
    )
    return evidence


def synthetic_vcep_specification(*, disease="DISEASE-TEST", transcript="NM_000001.1"):
    return {
        "name": "Synthetic PS1 VCEP specification",
        "version": "1.0.0",
        "source_url": "https://example.invalid/ps1-vcep-1",
        "specification_id": "SYNTHETIC-PS1-1",
        "scope": {
            "gene_id": "GENE-TEST",
            "disease_id": disease,
            "candidate_transcript_id_version": transcript,
        },
        "rules": [
            {
                "source_locator": "Synthetic VCEP PS1 table, row 1",
                "pathway": "SAME_AMINO_ACID_CHANGE",
                "when": {"candidate_variant_consequence": "MISSENSE"},
                "reference_when": {"classification": "PATHOGENIC"},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "MODERATE",
                    "decision_path": "SYNTHETIC_VCEP_PS1_MODERATE",
                },
            }
        ],
    }


def assert_outcome(evidence, status, strength):
    decision = evaluate_ps1(evidence)
    assert (decision.status.value, decision.strength.value) == (status, strength)
    assert decision.to_dict()["criterion"] == "PS1"
    return decision


def test_standard_path_activates_with_minimum_sufficient_evidence_and_audit_output():
    decision = assert_outcome(standard_features(), "ACTIVATED", "STRONG")
    payload = decision.to_dict()

    assert decision.decision_path == "SAME_AMINO_ACID_CHANGE"
    assert payload["missing_critical_features"] == []
    assert payload["reference_variants_used"] == [
        {
            "identifier": "REF-P-1",
            "classification": "PATHOGENIC",
            "status": "ACTIVATED",
            "strength": "STRONG",
            "decision_path": "SAME_AMINO_ACID_CHANGE",
        }
    ]
    assert payload["guideline"]["rule_id"] == "PS1"
    assert payload["guideline"]["source_url"]
    assert payload["features_used"]


def test_missing_unrelated_optional_annotations_do_not_block_standard_activation():
    evidence = standard_features()
    evidence.pop("candidate_coding_hgvs")
    evidence["reference_variants"][0].pop("coding_hgvs")

    assert_outcome(evidence, "ACTIVATED", "STRONG")


def test_different_amino_acid_substitution_at_same_residue_is_not_met():
    evidence = standard_features()
    evidence["reference_variants"][0]["protein_hgvs"] = "p.Arg34Trp"

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    assert any(
        entry["outcome"] == "DIFFERENT_PROTEIN_EFFECT"
        for entry in decision.to_dict()["exclusions_evaluated"]
    )


def test_same_underlying_allele_is_not_an_independent_ps1_comparison():
    evidence = standard_features()
    reference = evidence["reference_variants"][0]
    reference["genomic_allele"] = evidence["candidate_genomic_allele"]
    reference["coding_hgvs"] = evidence["candidate_coding_hgvs"]

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    assert any(
        entry["outcome"] == "SAME_ALLELE"
        for entry in decision.to_dict()["exclusions_evaluated"]
    )


def test_same_hgvs_text_on_different_transcript_versions_is_undetermined():
    evidence = standard_features()
    evidence["reference_variants"][0]["transcript_id_version"] = "NM_000001.2"

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "protein_substitution_equivalence" in decision.missing_critical_features


def test_generic_standard_path_rejects_likely_pathogenic_reference():
    evidence = standard_features()
    evidence["reference_variants"][0]["classification"] = "LIKELY_PATHOGENIC"

    assert_outcome(evidence, "NOT_MET", "NONE")


def test_reference_classification_circularity_is_an_exclusion():
    evidence = standard_features()
    evidence["reference_variants"][0]["classification_independence"] = False

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    assert any(
        entry["outcome"] == "CIRCULAR_PS1_EVIDENCE"
        for entry in decision.to_dict()["exclusions_evaluated"]
    )


def test_demonstrated_candidate_splice_effect_blocks_standard_protein_path():
    evidence = standard_features(candidate_splice_prediction="SPLICE_ALTERING")

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    assert any(
        entry["outcome"] == "SPLICE_MEDIATED_OR_POSSIBLE"
        for entry in decision.to_dict()["exclusions_evaluated"]
    )


def test_missing_splice_conclusions_are_not_assumed_to_be_negative():
    evidence = standard_features()
    evidence.pop("candidate_splice_prediction")
    evidence["reference_variants"][0].pop("splice_prediction")

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "candidate_splice_impact" in decision.missing_critical_features
    assert "candidate_splice_prediction" in decision.missing_critical_features


def test_walker_table2_noncanonical_same_nucleotide_pathogenic_reference_is_strong():
    decision = assert_outcome(table2_features(), "ACTIVATED", "STRONG")

    assert decision.decision_path == "SPLICING_SAME_PREDICTED_EVENT"
    assert decision.guideline.source_locator.startswith("Table 2")


def test_walker_table2_noncanonical_different_nucleotides_is_moderate():
    evidence = table2_features(candidate_splice_nucleotide_position="c.100+6")

    assert_outcome(evidence, "ACTIVATED", "MODERATE")


@pytest.mark.parametrize(
    ("same_nucleotide", "expected_strength"),
    ((True, "MODERATE"), (False, "SUPPORTING")),
)
def test_walker_table2_noncanonical_likely_pathogenic_reference_rows(
    same_nucleotide, expected_strength
):
    evidence = table2_features()
    evidence["reference_variants"][0]["classification"] = "LIKELY_PATHOGENIC"
    if not same_nucleotide:
        evidence["candidate_splice_nucleotide_position"] = "c.100+6"

    assert_outcome(evidence, "ACTIVATED", expected_strength)


def test_different_predicted_splice_events_are_not_met():
    evidence = table2_features()
    evidence.pop("splice_event_equivalence")
    evidence["reference_variants"][0]["spliceai_event_position"] = "exon-4"

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    splice_path = decision.to_dict()["reference_evaluations"][0]["pathways"][1]
    assert splice_path["decision_path"] == "SPLICING_SAME_PREDICTED_EVENT"
    assert splice_path["status"] == "NOT_MET"


def test_table2_does_not_invent_prediction_strength_when_comparison_is_missing():
    evidence = table2_features()
    evidence.pop("splice_prediction_strength_comparable")

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "splice_prediction_strength_comparable" in decision.missing_critical_features


def test_scoped_vcep_rule_can_apply_a_documented_strength_adjustment():
    evidence = standard_features(
        applicable_vcep_ps1_specification=[synthetic_vcep_specification()]
    )

    decision = assert_outcome(evidence, "ACTIVATED", "MODERATE")
    payload = decision.to_dict()
    assert decision.decision_path == "SYNTHETIC_VCEP_PS1_MODERATE"
    assert payload["guideline"]["name"] == "Synthetic PS1 VCEP specification"
    assert payload["guideline"]["version"] == "1.0.0"
    assert payload["guideline"]["rule_id"] == "SYNTHETIC-PS1-1"


def test_wrong_disease_scope_does_not_apply_vcep_rule():
    evidence = standard_features(
        applicable_vcep_ps1_specification=[
            synthetic_vcep_specification(disease="OTHER-DISEASE")
        ]
    )

    decision = assert_outcome(evidence, "ACTIVATED", "STRONG")
    assert decision.decision_path == "SAME_AMINO_ACID_CHANGE"


def test_unresolved_applicable_vcep_scope_blocks_generic_fallback():
    specification = synthetic_vcep_specification()
    specification["scope"]["inheritance_context"] = "AUTOSOMAL_DOMINANT"
    decision = assert_outcome(
        standard_features(applicable_vcep_ps1_specification=[specification]),
        "UNDETERMINED",
        "NONE",
    )

    assert "inheritance_context" in decision.missing_critical_features


def test_multiple_references_are_audited_but_never_summed_into_a_higher_strength():
    evidence = standard_features()
    second = deepcopy(evidence["reference_variants"][0])
    second.update(
        {
            "variant_id": "REF-P-2",
            "genomic_allele": "chr1:g.100A>T",
            "coding_hgvs": "c.100A>T",
        }
    )
    evidence["reference_variants"].append(second)

    decision = assert_outcome(evidence, "ACTIVATED", "STRONG")
    assert {row["identifier"] for row in decision.to_dict()["reference_variants_used"]} == {
        "REF-P-1",
        "REF-P-2",
    }
    assert "not summed" in decision.reason
    assert "references[1].classification" in decision.features_used


def test_missing_reference_classification_is_undetermined_not_negative():
    evidence = standard_features()
    evidence["reference_variants"][0].pop("classification")

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "references[0].classification" in decision.missing_critical_features


def test_missing_disease_context_is_undetermined():
    evidence = standard_features()
    evidence.pop("disease_id")

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "disease_id" in decision.missing_critical_features


def test_conflicting_reference_classification_is_undetermined():
    evidence = standard_features()
    evidence["reference_variants"][0]["classification_conflicts"] = True

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "references[0].classification_conflicts" in decision.missing_critical_features


def test_unsupported_synonymous_candidate_is_not_applicable():
    evidence = standard_features(candidate_variant_consequence="SYNONYMOUS")
    evidence.pop("candidate_splice_prediction")
    evidence["reference_variants"][0].pop("splice_prediction")

    assert_outcome(evidence, "NOT_APPLICABLE", "NONE")


def test_conflicting_derived_protein_equivalence_is_not_resolved_by_choosing_a_positive_path():
    evidence = standard_features(
        protein_change_match={"state": "CONFLICTING"}
    )

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "protein_substitution_equivalence" in decision.missing_critical_features


def test_provenance_is_emitted_and_input_is_not_mutated_or_order_dependent():
    evidence = standard_features()
    evidence["reference_variants"][0]["classification"] = {
        "state": "TRUE",
        "value": "PATHOGENIC",
        "provenance": {
            "source": "synthetic clinical assertion",
            "kind": "CURATED",
            "identifier": "REF-ASSERTION-1",
        },
    }
    original = deepcopy(evidence)

    first = evaluate_ps1(evidence).to_dict()
    second = evaluate_ps1(evidence).to_dict()

    assert evidence == original
    assert first == second
    assert any(
        item.get("identifier") == "REF-ASSERTION-1"
        for item in first["evidence_provenance"]
    )


def test_schema_rejects_unknown_biological_features_instead_of_silently_ignoring_them():
    evidence = standard_features(not_a_ps1_feature=True)

    with pytest.raises(SchemaValidationError, match="unsupported biological feature"):
        evaluate_ps1(evidence)


def test_reference_for_a_different_disease_context_is_not_met():
    evidence = standard_features()
    evidence["reference_variants"][0].update(
        {
            "disease_id": "OTHER-DISEASE",
            "disease_association": False,
        }
    )

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    assert any(
        entry["outcome"] == "MISMATCH"
        for entry in decision.to_dict()["exclusions_evaluated"]
    )


def test_reference_splice_mediated_mechanism_excludes_standard_protein_path():
    evidence = standard_features()
    evidence["reference_variants"][0].update(
        {
            "molecular_mechanism": "SPLICING",
            "splice_prediction": "SPLICE_ALTERING",
        }
    )

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    assert "non-protein-alteration mechanism" in decision.reason


@pytest.mark.parametrize(
    ("baseline", "expected_strength"),
    (
        ("PVS1", "SUPPORTING"),
        ("PVS1_STRONG", "STRONG"),
    ),
)
def test_walker_table2_canonical_rows_respect_the_applicable_pvs1_baseline(
    baseline, expected_strength
):
    decision = assert_outcome(
        canonical_table2_features(baseline=baseline),
        "ACTIVATED",
        expected_strength,
    )

    assert decision.decision_path == "SPLICING_SAME_PREDICTED_EVENT"
    assert decision.guideline.source_locator.startswith("Table 2")


@pytest.mark.parametrize(
    ("baseline", "classification", "expected_strength"),
    (
        ("PVS1", "PATHOGENIC", "SUPPORTING"),
        ("PVS1", "LIKELY_PATHOGENIC", "SUPPORTING"),
        ("PVS1_STRONG", "PATHOGENIC", "MODERATE"),
        ("PVS1_STRONG", "LIKELY_PATHOGENIC", "SUPPORTING"),
    ),
)
def test_walker_table2_canonical_to_noncanonical_rows(
    baseline, classification, expected_strength
):
    evidence = canonical_table2_features(baseline=baseline)
    evidence["reference_variants"][0].update(
        {
            "classification": classification,
            "genomic_allele": "chr1:g.100+5G>C",
            "coding_hgvs": "c.100+5G>C",
            "splice_motif_position": "WITHIN_MOTIF",
            "splice_nucleotide_position": "c.100+5",
        }
    )

    assert_outcome(evidence, "ACTIVATED", expected_strength)


def test_walker_table2_canonical_pvs1_likely_pathogenic_same_nucleotide_is_na():
    evidence = canonical_table2_features(baseline="PVS1")
    evidence["reference_variants"][0]["classification"] = "LIKELY_PATHOGENIC"

    decision = evaluate_ps1(evidence).to_dict()
    assert decision["reference_evaluations"][0]["pathways"][1]["status"] == "NOT_APPLICABLE"


def test_equivalent_reliable_rna_results_preserve_table2_activation():
    evidence = table2_features(
        candidate_rna_splicing_evidence=True,
        candidate_observed_splice_product="exon-3-skipped",
        splice_rna_assay_reliability="ESTABLISHED",
    )
    evidence["reference_variants"][0].update(
        {
            "rna_splicing_evidence": True,
            "observed_splice_product": "exon-3-skipped",
        }
    )

    assert_outcome(evidence, "ACTIVATED", "STRONG")


def test_reliable_rna_product_contradiction_prevents_splice_ps1_activation():
    evidence = table2_features(
        candidate_rna_splicing_evidence=True,
        candidate_observed_splice_product="exon-3-skipped",
        splice_rna_assay_reliability="ESTABLISHED",
    )
    evidence["reference_variants"][0].update(
        {
            "rna_splicing_evidence": True,
            "observed_splice_product": "exon-4-skipped",
        }
    )

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    splice_path = decision.to_dict()["reference_evaluations"][0]["pathways"][1]
    assert splice_path["status"] == "NOT_MET"


def test_different_residual_normal_splicing_profiles_prevent_splice_ps1_activation():
    evidence = table2_features(
        candidate_rna_splicing_evidence=True,
        candidate_observed_splice_product="exon-3-skipped",
        candidate_residual_normal_splicing=True,
        splice_rna_assay_reliability="ESTABLISHED",
    )
    evidence["reference_variants"][0].update(
        {
            "rna_splicing_evidence": True,
            "observed_splice_product": "exon-3-skipped",
            "residual_normal_splicing": False,
        }
    )

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    splice_path = decision.to_dict()["reference_evaluations"][0]["pathways"][1]
    assert splice_path["status"] == "NOT_MET"


def test_scoped_vcep_not_met_rule_takes_precedence_over_generic_activation():
    specification = synthetic_vcep_specification()
    specification["rules"][0]["decision"] = {
        "status": "NOT_MET",
        "strength": "NONE",
        "decision_path": "SYNTHETIC_VCEP_PS1_EXCLUDED",
    }
    evidence = standard_features(applicable_vcep_ps1_specification=[specification])

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    assert decision.decision_path == "SYNTHETIC_VCEP_PS1_EXCLUDED"
    assert decision.guideline.name == "Synthetic PS1 VCEP specification"


def test_irrelevant_optional_annotations_do_not_change_standard_decision():
    baseline = standard_features()
    enriched = deepcopy(baseline)
    enriched.update(
        {
            "candidate_variant_id": "CANDIDATE-1",
            "candidate_amino_acid_position": 34,
            "candidate_splice_prediction_magnitude": 0.01,
        }
    )
    enriched["reference_variants"][0].update(
        {
            "amino_acid_position": 34,
            "splice_prediction_magnitude": 0.01,
        }
    )

    first = evaluate_ps1(baseline)
    second = evaluate_ps1(enriched)
    assert (first.status.value, first.strength.value, first.decision_path) == (
        second.status.value,
        second.strength.value,
        second.decision_path,
    )


def test_protein_accession_mismatch_without_derived_mapping_is_undetermined():
    evidence = standard_features()
    evidence["reference_variants"][0]["protein_accession_version"] = "NP_000001.2"

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "protein_substitution_equivalence" in decision.missing_critical_features


def test_conflicting_genomic_and_coding_representations_do_not_prove_distinct_alleles():
    evidence = standard_features()
    reference = evidence["reference_variants"][0]
    reference["genomic_allele"] = "1:g.100A>G"  # alternate chromosome spelling
    reference["coding_hgvs"] = evidence["candidate_coding_hgvs"]

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "underlying_allele_equivalence" in decision.missing_critical_features


def test_categorical_false_evidence_envelope_is_rejected():
    evidence = standard_features()
    evidence["reference_variants"][0]["classification"] = {
        "state": "FALSE",
        "value": "PATHOGENIC",
    }

    with pytest.raises(SchemaValidationError, match="categorical evidence"):
        evaluate_ps1(evidence)


def test_splice_table2_requires_a_common_or_curated_equivalent_transcript_mapping():
    evidence = table2_features()
    evidence["reference_variants"][0]["transcript_id_version"] = "NM_999999.1"

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "splice_transcript_equivalence" in decision.missing_critical_features


def test_curated_splice_transcript_mapping_can_resolve_different_labels():
    evidence = table2_features()
    evidence["reference_variants"][0].update(
        {
            "transcript_id_version": "NM_999999.1",
            "splice_transcript_equivalence": True,
        }
    )

    assert_outcome(evidence, "ACTIVATED", "STRONG")


def test_reference_specific_splice_strength_comparison_overrides_legacy_default():
    evidence = table2_features(splice_prediction_strength_comparable=False)
    evidence["reference_variants"][0]["splice_prediction_strength_comparable"] = True

    decision = assert_outcome(evidence, "ACTIVATED", "STRONG")
    assert "references[0].splice_prediction_strength_comparable" in decision.features_used


def test_one_sided_reliable_normal_rna_result_contradicts_splice_prediction():
    evidence = table2_features(
        candidate_rna_splicing_evidence=True,
        candidate_observed_splice_product="NORMAL",
        splice_rna_assay_reliability="ESTABLISHED",
    )

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    splice_path = decision.to_dict()["reference_evaluations"][0]["pathways"][1]
    assert any(
        entry["outcome"] == "RNA_CONTRADICTS_SPLICE_PREDICTION"
        for entry in splice_path["exclusions_evaluated"]
    )


def test_exonic_protein_effect_disagreement_blocks_generic_splice_activation():
    evidence = table2_features(candidate_variant_consequence="MISSENSE")
    reference = evidence["reference_variants"][0]
    reference["variant_consequence"] = "MISSENSE"
    reference["protein_hgvs"] = "p.Arg34Trp"

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "protein_effect_mechanism_compatibility" in decision.missing_critical_features


def test_non_loss_of_function_product_conflicts_with_pvs1_splice_baseline():
    evidence = canonical_table2_features(baseline="PVS1")
    evidence["splice_product_functional_consequence"] = "NOT_LOSS_OF_FUNCTION"

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "splice_product_functional_consequence" in decision.missing_critical_features


def test_vcep_activation_cannot_bypass_same_amino_acid_safety_checks():
    evidence = standard_features(
        applicable_vcep_ps1_specification=[synthetic_vcep_specification()]
    )
    evidence["reference_variants"][0]["protein_hgvs"] = "p.Arg34Trp"

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    assert decision.decision_path == "SYNTHETIC_VCEP_PS1_MODERATE"


def test_incomplete_legacy_vcep_metadata_blocks_generic_fallback():
    decision = assert_outcome(
        standard_features(vcep_gene_scope="GENE-TEST"), "UNDETERMINED", "NONE"
    )
    assert "vcep_specifications" in decision.missing_critical_features


def test_vcep_activation_requires_an_explicit_reference_classification_condition():
    specification = synthetic_vcep_specification()
    specification["rules"][0]["reference_when"] = {}

    with pytest.raises(SchemaValidationError, match="reference-classification"):
        evaluate_ps1(
            standard_features(applicable_vcep_ps1_specification=[specification])
        )


def test_incompatible_same_nucleotide_motif_categories_are_undetermined():
    evidence = table2_features()
    evidence["reference_variants"][0]["splice_motif_position"] = "CANONICAL_DINUCLEOTIDE"

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "candidate_splice_motif_position" in decision.missing_critical_features


def test_compact_feature_method_provenance_is_preserved():
    evidence = standard_features()
    evidence["candidate_splice_prediction"] = {
        "state": "TRUE",
        "value": "NO_IMPACT",
        "feature_source": "synthetic predictor",
        "feature_method": "test method",
    }

    payload = evaluate_ps1(evidence).to_dict()
    assert any(item.get("method") == "test method" for item in payload["evidence_provenance"])


def test_redundant_vcep_version_metadata_must_match_the_pinned_snapshot():
    evidence = standard_features(
        applicable_vcep_ps1_specification=[synthetic_vcep_specification()],
        vcep_specification_version="9.9.9",
    )

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "VCEP scope/version metadata conflicts" in decision.reason


def test_pinned_vcep_numeric_splice_threshold_is_applied_without_a_generic_cutoff():
    specification = synthetic_vcep_specification()
    rule = specification["rules"][0]
    rule["pathway"] = "SPLICING_SAME_PREDICTED_EVENT"
    rule["when"] = {"candidate_splice_prediction_magnitude": {">=": 0.2}}
    rule["decision"] = {
        "status": "ACTIVATED",
        "strength": "MODERATE",
        "decision_path": "SYNTHETIC_VCEP_SPLICE_THRESHOLD",
    }

    above_threshold = table2_features(
        candidate_splice_prediction_magnitude=0.3,
        applicable_vcep_ps1_specification=[specification],
    )
    at_or_below_generic_path = table2_features(
        candidate_splice_prediction_magnitude=0.1,
        applicable_vcep_ps1_specification=[specification],
    )
    # The first result proves that the documented VCEP threshold, rather than
    # the generic Walker pairwise-strength assertion, is the sufficient
    # prediction-strength evidence for this scoped rule.
    above_threshold.pop("splice_prediction_strength_comparable")

    assert_outcome(above_threshold, "ACTIVATED", "MODERATE")
    assert_outcome(at_or_below_generic_path, "ACTIVATED", "STRONG")


@pytest.mark.parametrize(
    ("baseline", "classification", "status", "strength"),
    (
        ("PVS1", "PATHOGENIC", "ACTIVATED", "SUPPORTING"),
        ("PVS1", "LIKELY_PATHOGENIC", "NOT_APPLICABLE", "NONE"),
        ("PVS1_STRONG", "PATHOGENIC", "ACTIVATED", "STRONG"),
        ("PVS1_STRONG", "LIKELY_PATHOGENIC", "NOT_APPLICABLE", "NONE"),
    ),
)
def test_walker_canonical_rows_use_the_same_dinucleotide_not_the_same_coordinate(
    baseline, classification, status, strength
):
    evidence = canonical_table2_features(baseline=baseline)
    evidence["reference_variants"][0].update(
        {
            "classification": classification,
            "genomic_allele": "chr1:g.100+2G>C",
            "coding_hgvs": "c.100+2G>C",
            "splice_nucleotide_position": "c.100+2",
        }
    )

    splice_path = evaluate_ps1(evidence).to_dict()["reference_evaluations"][0]["pathways"][1]
    assert (splice_path["status"], splice_path["strength"]) == (status, strength)


def test_walker_canonical_row_does_not_require_an_exact_nucleotide_position():
    evidence = canonical_table2_features(baseline="PVS1_STRONG")
    evidence.pop("candidate_splice_nucleotide_position")
    evidence["reference_variants"][0].pop("splice_nucleotide_position")

    assert_outcome(evidence, "ACTIVATED", "STRONG")


def test_walker_same_coordinate_row_does_not_require_arbitrary_motif_identifiers():
    evidence = table2_features()
    evidence.pop("candidate_splice_motif_id")
    evidence["reference_variants"][0].pop("splice_motif_id")

    assert_outcome(evidence, "ACTIVATED", "STRONG")


def test_identical_genomic_allele_cannot_be_overridden_by_a_false_pairwise_assertion():
    evidence = standard_features(underlying_allele_equivalence=False)
    evidence["reference_variants"][0].update(
        {
            "genomic_allele": evidence["candidate_genomic_allele"],
            "coding_hgvs": "c.100A>T",
        }
    )

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "underlying_allele_equivalence" in decision.missing_critical_features


@pytest.mark.parametrize("fixture", (standard_features, table2_features))
def test_false_rna_flag_with_an_observed_product_is_never_silently_ignored(fixture):
    evidence = fixture(
        candidate_rna_splicing_evidence=False,
        candidate_observed_splice_product="EXON_SKIP",
    )

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "candidate_rna_splicing_evidence" in decision.missing_critical_features


def test_residual_normal_splicing_profiles_are_not_silently_ignored():
    evidence = table2_features(candidate_residual_normal_splicing=True)
    evidence["reference_variants"][0]["residual_normal_splicing"] = False

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "splice_rna_assay_reliability" in decision.missing_critical_features


def test_residual_normal_splicing_without_a_complete_result_blocks_standard_ps1():
    evidence = standard_features(candidate_residual_normal_splicing=True)

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "candidate_observed_splice_product" in decision.missing_critical_features


def test_uncertain_splice_prediction_is_undetermined_not_a_negative_finding():
    evidence = table2_features(candidate_splice_prediction="UNCERTAIN")

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "candidate_splice_prediction" in decision.missing_critical_features


def test_vcep_pairwise_when_uses_the_current_reference_override():
    specification = synthetic_vcep_specification()
    rule = specification["rules"][0]
    rule["when"] = {"splice_event_equivalence": True}
    rule["decision"] = {
        "status": "NOT_APPLICABLE",
        "strength": "NONE",
        "decision_path": "SYNTHETIC_VCEP_PAIRWISE_NOT_APPLICABLE",
    }
    evidence = standard_features(
        splice_event_equivalence=True,
        applicable_vcep_ps1_specification=[specification],
    )
    evidence["reference_variants"][0]["splice_event_equivalence"] = False

    decision = assert_outcome(evidence, "ACTIVATED", "STRONG")
    assert decision.decision_path == "SAME_AMINO_ACID_CHANGE"


def test_vcep_pairwise_unless_uses_the_current_reference_override():
    specification = synthetic_vcep_specification()
    rule = specification["rules"][0]
    rule["unless"] = {"splice_event_equivalence": True}
    evidence = standard_features(
        splice_event_equivalence=True,
        applicable_vcep_ps1_specification=[specification],
    )
    evidence["reference_variants"][0]["splice_event_equivalence"] = False

    decision = assert_outcome(evidence, "ACTIVATED", "MODERATE")
    assert decision.decision_path == "SYNTHETIC_VCEP_PS1_MODERATE"


def test_legacy_vcep_scope_metadata_without_a_matching_pinned_snapshot_blocks_fallback():
    specification = synthetic_vcep_specification()
    specification["scope"]["gene_id"] = "GENE-OTHER"
    evidence = standard_features(
        applicable_vcep_ps1_specification=[specification],
        vcep_gene_scope="GENE-TEST",
        vcep_disease_scope="DISEASE-TEST",
        vcep_specification_version="1.0.0",
    )

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "vcep_specifications" in decision.missing_critical_features


def test_vcep_condition_objects_cannot_be_empty():
    specification = synthetic_vcep_specification()
    specification["rules"][0]["when"] = {"candidate_variant_consequence": {}}

    with pytest.raises(SchemaValidationError, match="empty condition"):
        evaluate_ps1(standard_features(applicable_vcep_ps1_specification=[specification]))


def test_reference_only_vcep_score_threshold_cannot_replace_candidate_comparability():
    specification = synthetic_vcep_specification()
    rule = specification["rules"][0]
    rule["pathway"] = "SPLICING_SAME_PREDICTED_EVENT"
    rule["when"] = {"candidate_variant_consequence": "SPLICE_REGION"}
    rule["reference_when"] = {
        "classification": "PATHOGENIC",
        "splice_prediction_magnitude": {">=": 0.2},
    }
    rule["decision"] = {
        "status": "ACTIVATED",
        "strength": "MODERATE",
        "decision_path": "SYNTHETIC_VCEP_REFERENCE_THRESHOLD",
    }
    evidence = table2_features(applicable_vcep_ps1_specification=[specification])
    evidence.pop("splice_prediction_strength_comparable")
    evidence["reference_variants"][0]["splice_prediction_magnitude"] = 0.3

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "splice_prediction_strength_comparable" in decision.missing_critical_features


def test_vcep_metadata_provenance_is_preserved_in_the_final_audit():
    evidence = standard_features(
        applicable_vcep_ps1_specification=[synthetic_vcep_specification()],
        vcep_gene_scope={
            "state": "TRUE",
            "value": "GENE-TEST",
            "provenance": {
                "source": "synthetic VCEP metadata",
                "identifier": "VCEP-META-1",
            },
        },
    )

    payload = evaluate_ps1(evidence).to_dict()
    assert any(item.get("identifier") == "VCEP-META-1" for item in payload["evidence_provenance"])


def test_multiple_matching_vcep_versions_require_an_explicit_version_pin():
    version_one = synthetic_vcep_specification()
    version_two = deepcopy(version_one)
    version_two["version"] = "2.0.0"
    version_two["specification_id"] = "SYNTHETIC-PS1-2"
    version_two["rules"][0]["when"] = {"candidate_variant_consequence": "SPLICE_REGION"}

    decision = assert_outcome(
        standard_features(applicable_vcep_ps1_specification=[version_one, version_two]),
        "UNDETERMINED",
        "NONE",
    )
    assert "Multiple VCEP specification versions" in decision.reason


def test_unresolved_scoped_vcep_rule_blocks_a_competing_favorable_rule():
    specification = synthetic_vcep_specification()
    specification["rules"].append(
        {
            "source_locator": "Synthetic VCEP PS1 table, unresolved exclusion row",
            "pathway": "SAME_AMINO_ACID_CHANGE",
            "when": {"candidate_amino_acid_position": 34},
            "reference_when": {"classification": "PATHOGENIC"},
            "decision": {
                "status": "NOT_MET",
                "strength": "NONE",
                "decision_path": "SYNTHETIC_VCEP_UNRESOLVED_EXCLUSION",
            },
        }
    )

    decision = assert_outcome(
        standard_features(applicable_vcep_ps1_specification=[specification]),
        "UNDETERMINED",
        "NONE",
    )
    assert "candidate_amino_acid_position" in decision.missing_critical_features


@pytest.mark.parametrize(
    ("exclusion_field", "exclusion_value", "candidate_fields"),
    (
        ("unless", {"candidate_amino_acid_position": 34}, {"candidate_amino_acid_position": 34}),
        ("reference_unless", {"classification": "PATHOGENIC"}, {}),
    ),
)
def test_vcep_unless_conditions_are_explicit_exclusions_not_generic_fallbacks(
    exclusion_field, exclusion_value, candidate_fields
):
    specification = synthetic_vcep_specification()
    specification["rules"][0][exclusion_field] = exclusion_value
    evidence = standard_features(
        applicable_vcep_ps1_specification=[specification],
        **candidate_fields,
    )

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    assert decision.decision_path == "SYNTHETIC_VCEP_PS1_MODERATE"


def test_explicit_unknown_pairwise_reference_assertion_does_not_inherit_global_true():
    evidence = table2_features()
    reference = evidence["reference_variants"][0]
    reference.pop("spliceai_event_type")
    reference.pop("spliceai_event_position")
    reference["splice_event_equivalence"] = {"state": "UNKNOWN"}

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "references[0].splice_event_equivalence" in decision.missing_critical_features


@pytest.mark.parametrize("use_vcep", (False, True))
def test_splice_ps1_rejects_reference_with_non_splicing_pathogenic_mechanism(use_vcep):
    evidence = table2_features()
    evidence["reference_variants"][0]["molecular_mechanism"] = "OTHER"
    if use_vcep:
        specification = synthetic_vcep_specification()
        rule = specification["rules"][0]
        rule["pathway"] = "SPLICING_SAME_PREDICTED_EVENT"
        rule["when"] = {"candidate_splice_prediction_magnitude": {">=": 0.2}}
        rule["decision"] = {
            "status": "ACTIVATED",
            "strength": "MODERATE",
            "decision_path": "SYNTHETIC_VCEP_SPLICE_THRESHOLD",
        }
        evidence["candidate_splice_prediction_magnitude"] = 0.3
        evidence["applicable_vcep_ps1_specification"] = [specification]

    decision = assert_outcome(evidence, "NOT_MET", "NONE")
    assert any(
        item["outcome"] == "INCOMPATIBLE_SPLICE_MECHANISM"
        for pathway in decision.to_dict()["reference_evaluations"][0]["pathways"]
        for item in pathway["exclusions_evaluated"]
    )


def test_splice_ps1_requires_established_reference_splicing_mechanism():
    evidence = table2_features()
    evidence["reference_variants"][0].pop("molecular_mechanism")

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert "references[0].molecular_mechanism" in decision.missing_critical_features


def test_duplicate_reference_identity_with_conflicting_classifications_is_undetermined():
    evidence = standard_features()
    duplicate = deepcopy(evidence["reference_variants"][0])
    duplicate["classification"] = "VUS"
    evidence["reference_variants"].append(duplicate)

    decision = assert_outcome(evidence, "UNDETERMINED", "NONE")
    assert decision.decision_path == "PS1_DUPLICATE_REFERENCE_CONFLICT"
    assert set(decision.missing_critical_features) == {
        "references[0].classification",
        "references[1].classification",
    }
