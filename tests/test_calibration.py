import pytest

from adaptive_pp3 import evaluate_pp3
from .pp3_fixtures import missense, splicing


@pytest.mark.parametrize("field,value", [("score", None), ("score", 1.1),
                                         ("version", None), ("applicable", False),
                                         ("quality_ok", False), ("out_of_distribution", True),
                                         ("variant_id", "wrong"), ("transcript_id_version", "NM_BAD.1"),
                                         ("genome_build", "GRCh37")])
def test_invalid_prediction_cannot_activate(field, value):
    data = missense(0.9)
    data["predictions"][0][field] = value
    assert evaluate_pp3(data).status.value == "UNDETERMINED"


def test_calibration_version_mismatch():
    data = missense(0.9)
    data["predictions"][0]["calibration"] = "OTHER"
    assert evaluate_pp3(data).status.value == "UNDETERMINED"


def test_uncalibrated_context_does_not_vote():
    data = missense(0.8)
    data["predictions"] = [{"name": "SIFT", "score": 0.01, "interpretation": "DAMAGING"},
                           {"name": "POLYPHEN2", "score": 0.99, "interpretation": "DAMAGING"},
                           {"name": "CADD", "score": 35, "interpretation": "DAMAGING"}]
    assert evaluate_pp3(data).status.value == "UNDETERMINED"


def test_concordant_context_does_not_stack():
    data = missense(0.8)
    data["predictions"].append({"name": "SIFT", "score": 0.01, "interpretation": "DAMAGING"})
    assert evaluate_pp3(data).strength.value == "MODERATE"


def test_two_calibrated_positive_pathways_still_one_pp3():
    data = missense(0.8, canonical_splice_site=False,
                    splice_context_established=True, pvs1_splice_overlap=False,
                    rna_splice_evidence_used=False)
    data["predictions"].append({"name": "SPLICEAI", "score": 0.7, "version": "1.3",
                                "calibration": "WALKER_2023", "variant_id": "v1",
                                "transcript_id_version": "NM_TEST.1", "genome_build": "GRCh38"})
    result = evaluate_pp3(data)
    assert (result.status.value, result.strength.value, result.predictor_used) == (
        "ACTIVATED", "MODERATE", "REVEL")


def test_explicit_conflict_blocks():
    assert evaluate_pp3(missense(0.9, conflicting_computational_evidence=True)).status.value == "UNDETERMINED"


def test_splice_version_mismatch():
    data = splicing(0.9)
    data["predictions"][0]["calibration"] = "OTHER"
    assert evaluate_pp3(data).status.value == "UNDETERMINED"
