from adaptive_pp3 import evaluate_pp3
from .pp3_fixtures import missense, splicing


SCOPE = {"gene_id": "RUNX1", "disease_id": "MONDO:0011071",
         "vcep_id": "GN008", "vcep_version": "3.1"}


def _runx1(data):
    data.update(SCOPE)
    data["transcript_id_version"] = "NM_001754.4"
    for p in data["predictions"]:
        p["transcript_id_version"] = "NM_001754.4"
    return data


def test_vcep_revel_specific_threshold_and_strength():
    assert evaluate_pp3(_runx1(missense(0.879))).status.value == "NOT_MET"
    result = evaluate_pp3(_runx1(missense(0.932)))
    assert (result.status.value, result.strength.value) == ("ACTIVATED", "SUPPORTING")


def test_vcep_splice_specific_threshold():
    data = _runx1(splicing(0.38, intron_number=4))
    assert evaluate_pp3(data).status.value == "ACTIVATED"
    data["predictions"][0]["score"] = 0.379
    assert evaluate_pp3(data).status.value == "NOT_MET"


def test_vcep_scope_and_exclusions():
    assert evaluate_pp3(_runx1(splicing(0.9, intron_number=3))).status.value == "UNDETERMINED"
    data = _runx1(missense(0.9))
    data["gene_id"] = "OTHER"
    assert evaluate_pp3(data).status.value == "NOT_APPLICABLE"


def test_unpinned_runx1_does_not_use_generic_threshold():
    assert evaluate_pp3(missense(0.7, gene_id="RUNX1")).status.value == "UNDETERMINED"
