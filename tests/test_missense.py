import pytest

from adaptive_pp3 import evaluate_pp3
from .pp3_fixtures import missense


@pytest.mark.parametrize("score,strength", [(0.644, "SUPPORTING"), (0.773, "MODERATE"),
                                          (0.932, "STRONG"), (0.95, "STRONG")])
def test_revel_activation_boundaries(score, strength):
    result = evaluate_pp3(missense(score)).to_dict()
    assert (result["status"], result["strength"]) == ("ACTIVATED", strength)
    assert result["predictor_used"] == "REVEL"


@pytest.mark.parametrize("score", [0.643, 0.3, 0.290])
def test_revel_below_threshold(score):
    result = evaluate_pp3(missense(score)).to_dict()
    assert result["status"] == "NOT_MET"
    assert result["bp4_candidate"] == (score <= 0.290)


def test_missing_optional_predictors():
    assert evaluate_pp3(missense(0.9)).status.value == "ACTIVATED"


def test_wrong_variant_class():
    assert evaluate_pp3(missense(0.9, variant_consequence="FRAMESHIFT")).status.value == "UNDETERMINED"
