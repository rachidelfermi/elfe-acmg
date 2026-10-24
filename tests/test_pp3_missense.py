from adaptive_pp3 import Decision, Strength, evaluate_pp3


def case(score, **extra):
    return {
        "gene_id": "TEST1", "variant_consequence": "MISSENSE",
        "transcript_id_version": "NM_TEST.1", "biologically_relevant_transcript": True,
        "selected_missense_predictor": "REVEL",
        "predictions": [{"name": "REVEL", "score": score, "version": "1.3", "calibration": "PEJAVER_2022", "transcript_id_version": "NM_TEST.1"}],
        **extra,
    }


def test_revel_boundaries():
    for score, strength in [(0.644, Strength.SUPPORTING), (0.773, Strength.MODERATE), (0.932, Strength.STRONG)]:
        result = evaluate_pp3(case(score))
        assert (result.decision, result.strength) == (Decision.ACTIVATED, strength)


def test_revel_gap_and_bp4():
    assert evaluate_pp3(case(0.643)).decision == Decision.NOT_MET
    result = evaluate_pp3(case(0.290))
    assert result.decision == Decision.NOT_MET and result.bp4_candidate


def test_extra_predictor_does_not_raise_strength():
    input_ = case(0.7)
    input_["predictions"].append({"name": "CADD", "score": 99, "version": "1"})
    assert evaluate_pp3(input_).strength == Strength.SUPPORTING


def test_missing_calibration_and_selection():
    assert evaluate_pp3(case(0.9, selected_missense_predictor=None)).decision == Decision.UNDETERMINED
    input_ = case(0.9)
    input_["predictions"][0].pop("calibration")
    assert evaluate_pp3(input_).decision == Decision.UNDETERMINED


def test_conflict():
    assert evaluate_pp3(case(0.9, conflicting_computational_evidence=True)).decision == Decision.UNDETERMINED
