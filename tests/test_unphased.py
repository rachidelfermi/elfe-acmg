from adaptive_pm3 import evaluate_pm3
from tests.pm3_fixtures import case, data, second


def test_unphased_scores():
    for classification, points in (("PATHOGENIC", 0.5), ("LIKELY_PATHOGENIC", 0.25)):
        result = evaluate_pm3(data(case(phase="UNKNOWN", phase_method=None, second_variant=second(classification=classification))))
        assert result["total_points"] == points
        assert result["strength"] == ("SUPPORTING" if points == 0.5 else "NONE")


def test_vus_limited_points_and_benign_zero():
    vus = evaluate_pm3(data(case(second_variant=second(classification="VUS"))))
    assert vus["total_points"] == 0.25 and vus["status"] == "NOT_MET"
    assert evaluate_pm3(data(case(second_variant=second(classification="BENIGN"))))["total_points"] == 0
    assert evaluate_pm3(data(case(phase="UNKNOWN", second_variant=second(classification="VUS"))))["total_points"] == 0


def test_missing_classification_and_circularity():
    assert evaluate_pm3(data(case(second_variant=second(classification=None))))["status"] == "UNDETERMINED"
    assert evaluate_pm3(data(case(second_variant=second(classification_independent=False))))["status"] == "NOT_MET"
    assert evaluate_pm3(data(case(second_variant=second(classification_conflict=True))))["status"] == "UNDETERMINED"
