import pytest

from adaptive_pp3 import Decision, evaluate_pp3


def test_invalid_score_rejected():
    result = evaluate_pp3({"variant_consequence": "MISSENSE", "transcript_id_version": "NM_TEST.1",
                           "biologically_relevant_transcript": True, "selected_missense_predictor": "REVEL",
                           "predictions": [{"name": "REVEL", "score": 1.1, "version": "1", "calibration": "PEJAVER_2022"}]})
    assert result.decision == Decision.UNDETERMINED


def test_unsupported_predictor_cannot_activate():
    result = evaluate_pp3({"variant_consequence": "MISSENSE", "selected_missense_predictor": "AlphaMissense",
                           "predictions": [{"name": "AlphaMissense", "score": 0.99, "version": "1"}]})
    assert result.decision == Decision.UNDETERMINED
