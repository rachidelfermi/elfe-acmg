import pytest

from adaptive_bs3.guidelines import benign_strength_from_oddspath, control_count_strength
from tests.test_bs3_functional_assays import case
from adaptive_bs3 import evaluate_bs3


@pytest.mark.parametrize("value,strength", [
    (0.0529, "STRONG"), (0.053, "MODERATE"),
    (0.2299, "MODERATE"), (0.23, "SUPPORTING"),
    (0.4799, "SUPPORTING"), (0.48, "NONE"),
])
def test_published_strict_oddspath_boundaries(value, strength):
    assert benign_strength_from_oddspath(value) == strength


def test_control_count_route_and_bad_controls():
    changes = dict(validation_status="VALIDATED", calibration_validated=False,
                   odds_path=None, wildtype_control=True,
                   pathogenic_control_count=5, benign_control_count=6,
                   controls_independently_classified=True,
                   control_separation=True, replicate_consistent=True)
    assert evaluate_bs3(case(**changes)).strength == "MODERATE"
    assert evaluate_bs3(case(**{**changes, "pathogenic_control_count": 0})).status == "NOT_MET"
    assert evaluate_bs3(case(**{**changes, "benign_control_count": 0})).status == "NOT_MET"
    assert evaluate_bs3(case(**{**changes, "controls_independently_classified": False})).status == "NOT_MET"
    assert control_count_strength(1, 1) == "SUPPORTING"


def test_unverified_oddspath_does_not_bypass_controls():
    assert evaluate_bs3(case(calibration_validated=False)).status == "UNDETERMINED"
