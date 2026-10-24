from adaptive_ps3 import evaluate_ps3

from .ps3_fixtures import assay, evidence, formal_assay, supporting_assay


def test_validated_damaging_assay():
    result = evaluate_ps3(evidence())
    assert (result.status.value, result.strength.value) == ("ACTIVATED", "MODERATE")


def test_validated_normal_function_assay_reports_bs3_relevance_only():
    result = evaluate_ps3(evidence(assay(variant_result_classification="normal")))
    assert result.status.value == "NOT_MET"
    assert result.strength.value == "NONE"
    assert result.bs3_relevance == "POTENTIALLY_APPLICABLE"


def test_strong_effect_with_poor_validation_is_not_strength_evidence():
    weak = assay(
        assay_validation_status="unvalidated",
        variant_effect_size=1000,
        known_pathogenic_control_results=None,
        known_benign_control_results=None,
        assay_validation_sample_size=None,
    )
    result = evaluate_ps3(evidence(weak))
    assert result.status.value != "ACTIVATED"


def test_limited_validation_can_support_published_reduced_strength():
    result = evaluate_ps3(evidence(supporting_assay()))
    assert (result.status.value, result.strength.value) == ("ACTIVATED", "SUPPORTING")


def test_quantitative_assay_uses_formal_oddspath():
    result = evaluate_ps3(evidence(formal_assay(20)))
    assert result.strength.value == "STRONG"
    assert result.decision_path == "VALIDATED_QUANTITATIVE_ASSAY"


def test_qualitative_assay_uses_same_validation_principles():
    result = evaluate_ps3(evidence(assay(assay_type="qualitative")))
    assert result.status.value == "ACTIVATED"
    assert result.decision_path == "VALIDATED_QUALITATIVE_ASSAY"


def test_high_throughput_is_not_automatically_stronger():
    high = supporting_assay(assay_type="high_throughput", assay_name="MAVE")
    low = supporting_assay(assay_type="low_throughput", assay_name="single-variant")
    assert evaluate_ps3(evidence(high)).strength.value == "SUPPORTING"
    assert evaluate_ps3(evidence(low)).strength.value == "SUPPORTING"


def test_no_arbitrary_strength_from_effect_size_sensitivity_or_specificity():
    raw = supporting_assay(variant_effect_size=999, assay_sensitivity=1.0, assay_specificity=1.0)
    assert evaluate_ps3(evidence(raw)).strength.value == "SUPPORTING"


def test_formal_oddspath_thresholds_are_strict_and_very_strong_is_explicit():
    assert evaluate_ps3(evidence(formal_assay(350))).strength.value == "STRONG"
    assert evaluate_ps3(evidence(formal_assay(351))).strength.value == "VERY_STRONG"
