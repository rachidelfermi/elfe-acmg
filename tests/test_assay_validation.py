from adaptive_ps3 import evaluate_ps3

from .ps3_fixtures import assay, evidence, formal_assay, supporting_assay


def test_insufficient_control_classification_quality_rejects_calibration():
    result = evaluate_ps3(evidence(assay(control_variant_classification_quality=False)))
    assert result.status.value == "NOT_MET"


def test_pathogenic_and_benign_controls_reach_moderate_at_eleven_total():
    result = evaluate_ps3(evidence(assay()))
    assert result.strength.value == "MODERATE"


def test_incorrect_disease_mechanism_is_not_applicable():
    raw = assay(assay_measures_disease_mechanism=False)
    assert evaluate_ps3(evidence(raw)).status.value == "NOT_APPLICABLE"


def test_lof_assay_for_lof_disease_is_compatible():
    raw = assay(loss_of_function_effect=True)
    assert evaluate_ps3(evidence(raw)).status.value == "ACTIVATED"


def test_lof_assay_for_gof_disease_is_not_ps3_evidence():
    raw = assay(loss_of_function_effect=True)
    result = evaluate_ps3(evidence(raw, disease_mechanism="gain_of_function"))
    assert result.status.value == "NOT_MET"


def test_gof_assay_for_gof_disease_is_compatible():
    raw = assay(gain_of_function_effect=True)
    result = evaluate_ps3(evidence(raw, disease_mechanism="gain_of_function"))
    assert result.status.value == "ACTIVATED"


def test_functional_rescue_is_audited_but_does_not_inflate_strength():
    raw = supporting_assay(functional_rescue_evidence=True)
    result = evaluate_ps3(evidence(raw))
    assert result.strength.value == "SUPPORTING"
    assert "functional_rescue_evidence" in result.features_used


def test_replicate_inconsistency_blocks_ps3():
    raw = assay(variant_result_replicate_consistency=False)
    assert evaluate_ps3(evidence(raw)).status.value == "NOT_MET"


def test_experimental_batch_effect_blocks_ps3():
    raw = assay(experimental_batch_effect=True)
    assert evaluate_ps3(evidence(raw)).status.value == "NOT_MET"


def test_assay_saturation_blocks_ps3():
    raw = assay(assay_saturation_or_floor_effect=True)
    assert evaluate_ps3(evidence(raw)).status.value == "NOT_MET"


def test_failed_quality_control_or_control_separation_blocks_ps3():
    assert evaluate_ps3(evidence(assay(assay_quality_control=False))).status.value == "NOT_MET"
    assert evaluate_ps3(evidence(assay(control_result_separation=False))).status.value == "NOT_MET"


def test_formal_normal_result_is_only_flagged_for_separate_bs3_evaluation():
    raw = formal_assay(0.02, variant_result_classification="normal")
    result = evaluate_ps3(evidence(raw))
    assert result.status.value == "NOT_MET"
    assert result.strength.value == "NONE"
    assert result.bs3_relevance == "POTENTIALLY_APPLICABLE"


def test_multiple_consistent_independent_studies_do_not_stack():
    one = supporting_assay(assay_id="a", study_id="s1", study_independence=True)
    two = supporting_assay(assay_id="b", study_id="s2", study_independence=True)
    result = evaluate_ps3(evidence(one, two))
    assert result.strength.value == "SUPPORTING"
    assert len(result.assays_used) == 1


def test_multiple_equally_credible_conflicting_studies_are_undetermined():
    abnormal = formal_assay(20, assay_id="abnormal")
    normal = formal_assay(0.02, assay_id="normal", variant_result_classification="normal")
    result = evaluate_ps3(evidence(abnormal, normal))
    assert result.status.value == "UNDETERMINED"
    assert set(result.conflicting_evidence) == {"abnormal", "normal"}


def test_duplicate_overlapping_studies_never_stack():
    one = supporting_assay(assay_id="a", study_id="same", evidence_overlap=True)
    two = supporting_assay(assay_id="b", study_id="same", evidence_overlap=True)
    result = evaluate_ps3(evidence(one, two))
    assert result.strength.value == "SUPPORTING"


def test_better_validated_mechanism_relevant_assay_can_resolve_conflict():
    abnormal = formal_assay(20, assay_id="best", mechanism_compatibility=True)
    normal = supporting_assay(assay_id="weaker", variant_result_classification="normal")
    result = evaluate_ps3(evidence(abnormal, normal))
    assert result.status.value == "ACTIVATED"
    assert result.assays_used == ("best",)
