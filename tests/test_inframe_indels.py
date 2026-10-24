import pytest

from adaptive_pm4 import evaluate_pm4


BASE = dict(gene_id="GENE", disease_id="DISEASE", transcript_id_version="NM_1.1",
            transcript_biological_relevance=True, reading_frame_preserved=True,
            repeat_region_status=False)


@pytest.mark.parametrize("consequence,difference", [
    ("INFRAME_DELETION", -2), ("INFRAME_INSERTION", 2)])
def test_nonrepeat_activation(consequence, difference):
    result = evaluate_pm4({**BASE, "variant_consequence": consequence,
                           "protein_length_difference": difference})
    assert (result.status, result.strength, result.decision_path) == (
        "ACTIVATED", "MODERATE", consequence)


def test_no_length_change_and_frameshift():
    zero = evaluate_pm4({**BASE, "variant_consequence": "INFRAME_DELETION",
                         "protein_length_difference": 0})
    frame = evaluate_pm4({**BASE, "variant_consequence": "FRAMESHIFT"})
    assert zero.status == "NOT_MET"
    assert frame.status == "NOT_APPLICABLE"


def test_derived_length_and_irrelevant_features():
    base = {**BASE, "variant_consequence": "INFRAME_INSERTION",
            "reference_protein_length": 100, "alternate_protein_length": 101}
    first = evaluate_pm4(base)
    second = evaluate_pm4({**base, "critical_protein_domain_overlap": True,
                           "evidence_source": "synthetic"})
    assert first.status == second.status == "ACTIVATED"
    assert first.strength == second.strength == "MODERATE"


def test_conflicting_length():
    r = evaluate_pm4({**BASE, "variant_consequence": "INFRAME_INSERTION",
                      "reference_protein_length": 100, "alternate_protein_length": 101,
                      "protein_length_difference": 2})
    assert r.status == "UNDETERMINED"
