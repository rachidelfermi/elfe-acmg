from adaptive_pm4 import Evidence, EvidenceState, evaluate_pm4


BASE = dict(gene_id="GENE", disease_id="DISEASE", transcript_id_version="NM_1.1",
            transcript_biological_relevance=True, variant_consequence="INFRAME_DELETION",
            reading_frame_preserved=True, protein_length_difference=-1,
            repeat_region_status=False)


def test_missing_protein_consequence_and_frame():
    missing = dict(BASE)
    del missing["protein_length_difference"]
    assert evaluate_pm4(missing).status == "UNDETERMINED"
    assert evaluate_pm4({**BASE, "reading_frame_preserved": None}).status == "UNDETERMINED"


def test_transcript_mismatch_and_rescue():
    assert evaluate_pm4({**BASE, "disease_relevant_transcript": "NM_2.1"}).status == "NOT_MET"
    assert evaluate_pm4({**BASE, "alternative_transcript_rescue": True}).status == "NOT_MET"


def test_rna_contradiction_and_multiple_transcripts():
    assert evaluate_pm4({**BASE, "rna_splicing_evidence": "CONTRADICTS_LENGTH_CHANGE"}).status == "NOT_MET"
    assert evaluate_pm4({**BASE, "relevant_transcript_consequences": (
        "INFRAME_DELETION", "FRAMESHIFT")}).status == "UNDETERMINED"


def test_alignment_and_conflicting_annotations():
    assert evaluate_pm4({**BASE, "indel_alignment_reliability": False}).status == "UNDETERMINED"
    assert evaluate_pm4({**BASE, "variant_consequence": Evidence(EvidenceState.CONFLICTING)}).status == "UNDETERMINED"


def test_optional_absence_and_no_arbitrary_upgrade():
    result = evaluate_pm4(BASE)
    assert (result.status, result.strength) == ("ACTIVATED", "MODERATE")
    assert result.missing_critical_features == ()
    assert evaluate_pm4({**BASE, "protein_length_difference": -50,
                         "critical_protein_domain_overlap": True}).strength == "MODERATE"


def test_same_consequence_pvs1_overlap():
    assert evaluate_pm4({**BASE, "pvs1_applied_to_same_consequence": True}).status == "NOT_MET"
