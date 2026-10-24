from adaptive_bs3 import evaluate_bs3
from tests.test_bs3_functional_assays import case


def test_missing_context_and_assay():
    assert evaluate_bs3({"assays": []}).status == "UNDETERMINED"
    assert evaluate_bs3({"context": {"gene_id": "GENE", "disease_id": "D", "disease_mechanism": "LOSS_OF_FUNCTION"}}).status == "UNDETERMINED"


def test_explicit_evidence_states():
    data = case()
    data["assays"][0]["result_class"] = {"state": "CONFLICTING"}
    assert evaluate_bs3(data).status == "UNDETERMINED"
    assert "a1.result_class" in evaluate_bs3(data).missing_critical_features


def test_optional_evidence_is_not_required():
    assert evaluate_bs3(case()).status == "ACTIVATED"


def test_known_reliability_failures_override_calibrated_normal_result():
    for field in ("assay_sensitivity_adequate", "mechanism_coverage_complete",
                  "expression_adequate", "result_uncertainty_acceptable",
                  "control_classification_quality"):
        assert evaluate_bs3(case(**{field: False})).status == "NOT_MET"


def test_variant_mechanism_and_ps3_overlap():
    data = case()
    data["context"]["variant_mechanism_compatibility"] = False
    assert evaluate_bs3(data).status == "NOT_MET"
    assert evaluate_bs3(case(ps3_claimed=True)).status == "UNDETERMINED"


def test_nonsignificant_comparison_is_not_normal_function():
    result = evaluate_bs3(case(result_class="NON_SIGNIFICANT"))
    assert (result.status, result.strength) == ("UNDETERMINED", "NONE")
