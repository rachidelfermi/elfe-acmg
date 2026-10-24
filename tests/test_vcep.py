from adaptive_pm4 import evaluate_pm4


BASE = dict(gene_id="MECP2", disease_id="Rett syndrome",
            transcript_id_version="NM_004992.4", transcript_biological_relevance=True,
            vcep_specification_version="rett_angelman_2022")


def test_panel_indel_size_strength():
    base = {**BASE, "variant_consequence": "INFRAME_DELETION",
            "reading_frame_preserved": True, "repeat_region_status": False}
    small = evaluate_pm4({**base, "protein_length_difference": -2})
    large = evaluate_pm4({**base, "protein_length_difference": -3})
    assert small.strength == "SUPPORTING"
    assert large.strength == "MODERATE"


def test_panel_small_functional_region_exception():
    r = evaluate_pm4({**BASE, "variant_consequence": "INFRAME_DELETION",
                      "reading_frame_preserved": True, "repeat_region_status": False,
                      "protein_length_difference": -1,
                      "critical_protein_domain_overlap": True,
                      "affected_region_functional_evidence": True})
    assert r.strength == "MODERATE"


def test_panel_repeat_and_stop_loss():
    repeat = evaluate_pm4({**BASE, "variant_consequence": "INFRAME_DELETION",
                           "reading_frame_preserved": True,
                           "protein_length_difference": -3,
                           "protein_accession_version": "NP_004983.1",
                           "protein_change_start": 390, "protein_change_end": 392})
    stop = evaluate_pm4({**BASE, "variant_consequence": "STOP_LOSS",
                         "native_stop_codon_disrupted": True,
                         "protein_extension_length": 4})
    assert repeat.status == "NOT_MET"
    assert (stop.status, stop.strength) == ("ACTIVATED", "STRONG")


def test_profile_cannot_transfer_to_another_gene():
    r = evaluate_pm4({**BASE, "gene_id": "OTHER", "variant_consequence": "STOP_LOSS"})
    assert r.status == "UNDETERMINED"
