from adaptive_pm3 import evaluate_pm3
from tests.pm3_fixtures import case, data, second


def test_confirmed_pathogenic_and_likely_pathogenic_both_moderate():
    for classification in ("PATHOGENIC", "LIKELY_PATHOGENIC"):
        result = evaluate_pm3(data(case(second_variant=second(classification=classification))))
        assert (result["status"], result["strength"], result["total_points"]) == ("ACTIVATED", "MODERATE", 1.0)


def test_parent_and_molecular_phase():
    parent = case(phase=None, phase_method=None, maternal_candidate="REF", maternal_second="HET")
    molecular = case(phase="TRANS", phase_method="long_read")
    assert evaluate_pm3(data(parent))["total_points"] == 1.0
    assert evaluate_pm3(data(molecular))["total_points"] == 1.0


def test_cis_and_conflict():
    assert evaluate_pm3(data(case(phase="CIS")))["status"] == "NOT_MET"
    assert evaluate_pm3(data(case(phase="CONFLICTING")))["status"] == "UNDETERMINED"
    assert evaluate_pm3(data(case(phase="TRANS", maternal_candidate="HET", maternal_second="REF")))["total_points"] == 1.0
    assert evaluate_pm3(data(case(phase="CIS", maternal_candidate="HET", maternal_second="REF")))["status"] == "UNDETERMINED"


def test_no_unproved_confirmed_phase():
    result = evaluate_pm3(data(case(phase="TRANS", phase_method=None)))
    assert result["total_points"] == 0
    assert result["status"] == "UNDETERMINED"
