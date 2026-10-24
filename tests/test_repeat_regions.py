from adaptive_pm4 import evaluate_pm4


BASE = dict(gene_id="GENE", disease_id="DISEASE", transcript_id_version="NM_1.1",
            transcript_biological_relevance=True, variant_consequence="INFRAME_DELETION",
            reading_frame_preserved=True, protein_length_difference=-1)


def test_repeat_exclusion():
    r = evaluate_pm4({**BASE, "repeat_region_status": True})
    assert (r.status, r.strength) == ("NOT_MET", "NONE")


def test_benign_repeat_polymorphism():
    r = evaluate_pm4({**BASE, "repeat_region_status": False,
                      "known_benign_repeat_length_variation": True})
    assert r.status == "NOT_MET"


def test_unknown_repeat_status():
    r = evaluate_pm4(BASE)
    assert r.status == "UNDETERMINED"
    assert "repeat_region_status" in r.missing_critical_features


def test_regional_benign_indel_requires_review():
    r = evaluate_pm4({**BASE, "repeat_region_status": False,
                      "known_benign_indel_same_region": True})
    assert r.status == "UNDETERMINED"
