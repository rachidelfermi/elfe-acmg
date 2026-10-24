from dataclasses import replace

import pytest

from adaptive_bp4 import BiologicalFeatures, Prediction as RawPrediction, State, Status, Strength, evaluate_bp4


def Prediction(name, score=None, **kwargs):
    checks = dict.fromkeys(("variant_match", "transcript_match", "genome_build_match", "model_applicability", "quality"), State.KNOWN)
    return RawPrediction(name, score, **(checks | kwargs))


def case(*predictions):
    return BiologicalFeatures("v1", "GENE", consequence="missense", transcript_id_version="NM_1.1", biologically_relevant_transcript=State.KNOWN, predictions=predictions)


@pytest.mark.parametrize("score,status,strength", [
    (0.1, Status.ACTIVATED, Strength.MODERATE),
    (0.183, Status.ACTIVATED, Strength.MODERATE),
    (0.290, Status.ACTIVATED, Strength.SUPPORTING),
    (0.291, Status.NOT_MET, Strength.NONE),
    (0.016, Status.ACTIVATED, Strength.STRONG),
])
def test_revel_boundaries(score, status, strength):
    result = evaluate_bp4(case(Prediction("REVEL", score)))
    assert (result.status, result.strength) == (status, strength)


@pytest.mark.parametrize("prediction", [
    Prediction("REVEL"), Prediction("REVEL", 1.01),
    Prediction("REVEL", 0.1, variant_match=State.CONFLICTING),
    Prediction("REVEL", 0.1, transcript_match=State.UNKNOWN),
    Prediction("REVEL", 0.1, genome_build_match=State.CONFLICTING),
    Prediction("REVEL", 0.1, model_applicability=State.NOT_APPLICABLE),
    Prediction("REVEL", 0.1, out_of_distribution=True),
])
def test_unreliable_scores_do_not_activate(prediction):
    assert evaluate_bp4(case(prediction)).status != Status.ACTIVATED


def test_bayesdel_single_predictor_and_conflict():
    assert evaluate_bp4(case(Prediction("BayesDel_noAF", -0.4))).strength == Strength.MODERATE
    assert evaluate_bp4(case(Prediction("REVEL", 0.1), Prediction("BayesDel_noAF", 0.3))).status == Status.UNDETERMINED
    assert evaluate_bp4(case(Prediction("REVEL", 0.1), Prediction("BayesDel_noAF", -0.4))).strength == Strength.MODERATE
    assert evaluate_bp4(case(Prediction("REVEL", 0.4), Prediction("BayesDel_noAF", -0.4))).strength == Strength.MODERATE
    assert evaluate_bp4(case(Prediction("REVEL", 0.1, quality=State.UNKNOWN), Prediction("BayesDel_noAF", -0.4))).strength == Strength.MODERATE


def test_unsupported_predictors_do_not_vote():
    result = evaluate_bp4(case(Prediction("SIFT", 0.9), Prediction("PolyPhen2", 0.01), Prediction("CADD", 1.0)))
    assert result.status == Status.UNDETERMINED


def test_irrelevant_predictions_do_not_change_revel():
    base = evaluate_bp4(case(Prediction("REVEL", 0.25)))
    extra = evaluate_bp4(case(Prediction("REVEL", 0.25), Prediction("AlphaMissense", 0.01)))
    assert (base.status, base.strength, base.predictor_used) == (extra.status, extra.strength, extra.predictor_used)


def test_overlap_and_missing_identity():
    assert evaluate_bp4(replace(case(Prediction("REVEL", 0.1)), pp3_evidence_conflict=True)).status == Status.UNDETERMINED
    assert evaluate_bp4(replace(case(Prediction("REVEL", 0.1)), variant_id=None)).status == Status.UNDETERMINED
    assert evaluate_bp4(replace(case(Prediction("REVEL", 0.1)), biologically_relevant_transcript=State.UNKNOWN)).status == Status.UNDETERMINED


def test_conservation_and_structure_alone():
    assert evaluate_bp4(case(Prediction("phyloP", -1.0), Prediction("AlphaFold", 0.99))).status == Status.UNDETERMINED
