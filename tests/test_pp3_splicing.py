from adaptive_pp3 import Decision, Strength, evaluate_pp3


def case(score, **extra):
    return {
        "gene_id": "TEST1", "variant_consequence": "INTRONIC",
        "transcript_id_version": "NM_000001.1", "biologically_relevant_transcript": True,
        "splice_context_established": True, "canonical_splice_site": False,
        "pvs1_splice_overlap": False, "rna_splice_evidence_used": False,
        "predictions": [{"name": "SpliceAI", "score": score, "version": "1.3", "calibration": "WALKER_2023", "transcript_id_version": "NM_000001.1"}],
        **extra,
    }


def test_splice_thresholds():
    result = evaluate_pp3(case(0.2))
    assert (result.decision, result.strength) == (Decision.ACTIVATED, Strength.SUPPORTING)
    assert evaluate_pp3(case(0.199)).decision == Decision.NOT_MET
    assert evaluate_pp3(case(0.1)).bp4_candidate


def test_splice_exclusions():
    assert evaluate_pp3(case(0.8, canonical_splice_site=True)).decision == Decision.UNDETERMINED
    assert evaluate_pp3(case(0.8, pvs1_splice_overlap=True)).decision == Decision.UNDETERMINED
    assert evaluate_pp3(case(0.8, rna_splice_evidence_used=True)).decision == Decision.UNDETERMINED
    assert evaluate_pp3(case(0.8, transcript_id_version=None)).decision == Decision.UNDETERMINED


def test_splice_score_does_not_clear_missense_effect():
    assert not evaluate_pp3(case(0.05, variant_consequence="MISSENSE")).bp4_candidate
