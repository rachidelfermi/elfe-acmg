from adaptive_pm4 import evaluate_pm4


BASE = dict(gene_id="GENE", disease_id="DISEASE", transcript_id_version="NM_1.1",
            transcript_biological_relevance=True, variant_consequence="STOP_LOSS",
            native_stop_codon_disrupted=True)


def test_extension_and_downstream_stop():
    r = evaluate_pm4({**BASE, "protein_extension_length": 1,
                      "downstream_inframe_stop_codon": True})
    assert (r.status, r.strength) == ("ACTIVATED", "MODERATE")


def test_alternative_termination_with_no_extension():
    r = evaluate_pm4({**BASE, "protein_extension_length": 0,
                      "downstream_inframe_stop_codon": True})
    assert r.status == "NOT_MET"


def test_extension_without_confirmed_stop_loss():
    data = {**BASE, "protein_extension_length": 5}
    del data["native_stop_codon_disrupted"]
    r = evaluate_pm4(data)
    assert r.status == "UNDETERMINED"
    assert "native_stop_codon_disrupted" in r.missing_critical_features


def test_stop_retained_or_rescued():
    assert evaluate_pm4({**BASE, "native_stop_codon_disrupted": False,
                         "protein_extension_length": 5}).status == "NOT_MET"
    assert evaluate_pm4({**BASE, "protein_extension_length": 5,
                         "stop_loss_rescue_evidence": True}).status == "NOT_MET"
