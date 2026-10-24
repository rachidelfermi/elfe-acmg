from adaptive_pm3 import evaluate_pm3
from tests.pm3_fixtures import case, data, second


def test_independent_cases_reach_strong_and_boundaries():
    cases = [case(case_id=f"c{i}", family_id=f"f{i}", independence_confirmed=True,
                  second_variant=second(variant_id=f"partner{i}")) for i in range(4)]
    for n, strength in ((1, "MODERATE"), (2, "STRONG"), (4, "VERY_STRONG")):
        result = evaluate_pm3(data(*cases[:n]))
        assert result["strength"] == strength


def test_duplicate_and_related_cases_do_not_accumulate():
    first = case(case_id="c1", family_id="f1", independence_confirmed=True)
    duplicate = case(case_id="c1", family_id="f2", independence_confirmed=True)
    related = case(case_id="c2", family_id="f1", independence_confirmed=True)
    assert evaluate_pm3(data(first, duplicate, related))["total_points"] == 1.0


def test_missing_independence_keeps_existing_activation():
    first = case(case_id="c1", family_id="f1", independence_confirmed=True)
    unknown = case(case_id="c2", family_id=None)
    result = evaluate_pm3(data(first, unknown))
    assert result["status"] == "ACTIVATED" and result["total_points"] == 1.0


def test_reported_case_overlap_is_not_counted():
    first = case(case_id="c1", family_id="f1", independence_confirmed=True)
    overlap = case(case_id="c2", family_id="f2", independence_confirmed=True, case_overlap=True)
    result = evaluate_pm3(data(first, overlap))
    assert result["total_points"] == 1.0
    assert "overlapping_case_report" in result["exclusions_evaluated"]
