from adaptive_pp3 import evaluate_pp3
from .pp3_fixtures import missense


def test_no_scores():
    data = missense()
    data["predictions"] = []
    result = evaluate_pp3(data)
    assert result.status.value == "UNDETERMINED"
    assert "calibrated_prediction" in result.missing_critical_features


def test_conservation_or_structure_alone():
    for field in ("conservation", "structural_predictions"):
        data = missense()
        data["predictions"] = []
        data[field] = {"score": 999}
        assert evaluate_pp3(data).status.value == "UNDETERMINED"


def test_relevant_transcript_missing():
    assert evaluate_pp3(missense(0.9, biologically_relevant_transcript=None)).status.value == "UNDETERMINED"


def test_mechanism_incompatible():
    assert evaluate_pp3(missense(0.9, mechanism_compatible=False)).status.value == "NOT_APPLICABLE"


def test_bp4_conflict():
    assert evaluate_pp3(missense(0.9, bp4_evidence_conflict=True)).status.value == "UNDETERMINED"
