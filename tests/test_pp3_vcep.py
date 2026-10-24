from adaptive_pp3 import Decision, Strength, evaluate_pp3


def test_runx1_missense_override():
    base = {"gene_id": "RUNX1", "disease_id": "MONDO:0011071", "variant_consequence": "MISSENSE", "transcript_id_version": "NM_001754.4", "biologically_relevant_transcript": True, "vcep_id": "GN008", "vcep_version": "3.1", "selected_missense_predictor": "REVEL"}
    base["predictions"] = [{"name": "REVEL", "score": 0.879, "version": "1.3", "calibration": "PEJAVER_2022", "transcript_id_version": "NM_001754.4"}]
    assert evaluate_pp3(base).decision == Decision.NOT_MET
    base["predictions"][0]["score"] = 0.95
    assert evaluate_pp3(base).strength == Strength.SUPPORTING


def test_runx1_splice_scope():
    base = {"gene_id": "RUNX1", "disease_id": "MONDO:0011071", "variant_consequence": "INTRONIC", "intron_number": 4,
            "vcep_id": "GN008", "vcep_version": "3.1", "transcript_id_version": "NM_001754.4",
            "biologically_relevant_transcript": True, "splice_context_established": True,
            "canonical_splice_site": False, "pvs1_splice_overlap": False, "rna_splice_evidence_used": False,
            "predictions": [{"name": "SpliceAI", "score": 0.38, "version": "1.3", "calibration": "WALKER_2023", "transcript_id_version": "NM_001754.4"}]}
    assert evaluate_pp3(base).decision == Decision.ACTIVATED
    base["intron_number"] = 3
    assert evaluate_pp3(base).decision == Decision.UNDETERMINED


def test_unpinned_vcep_fails_closed():
    assert evaluate_pp3({"gene_id": "RUNX1", "variant_consequence": "MISSENSE"}).decision == Decision.UNDETERMINED
