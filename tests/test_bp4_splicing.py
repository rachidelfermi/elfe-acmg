from dataclasses import replace

import pytest

from adaptive_bp4 import BiologicalFeatures, Prediction as RawPrediction, State, Status, Strength, evaluate_bp4


def Prediction(name, score=None, **kwargs):
    checks = dict.fromkeys(("variant_match", "transcript_match", "genome_build_match", "model_applicability", "quality"), State.KNOWN)
    return RawPrediction(name, score, **(checks | kwargs))


def case(score):
    return BiologicalFeatures("v2", "GENE", consequence="noncanonical_splice", transcript_id_version="NM_2.1", biologically_relevant_transcript=State.KNOWN, protein_effect_relevant=False, splice_region_eligible=True, predictions=(Prediction("SpliceAI", score),))


@pytest.mark.parametrize("score,status", [(0.0, Status.ACTIVATED), (0.1, Status.ACTIVATED), (0.1001, Status.NOT_MET)])
def test_spliceai_boundaries(score, status):
    result = evaluate_bp4(case(score))
    assert result.status == status
    if status == Status.ACTIVATED:
        assert result.strength == Strength.MODERATE


def test_splicing_eligibility_and_overlap():
    assert evaluate_bp4(replace(case(0.01), protein_effect_relevant=True)).status == Status.UNDETERMINED
    assert evaluate_bp4(replace(case(0.01), splice_region_eligible=None)).status == Status.UNDETERMINED
    assert evaluate_bp4(replace(case(0.01), bp7_same_prediction_prohibited=True)).status == Status.UNDETERMINED
    assert evaluate_bp4(replace(case(0.01), splice_evidence_double_counting=True)).status == Status.UNDETERMINED
    assert evaluate_bp4(replace(case(0.01), predictions=(Prediction("SpliceAI", 0.01, transcript_match=State.CONFLICTING),))).status == Status.UNDETERMINED


def test_low_spliceai_does_not_resolve_missense_protein_effect():
    assert evaluate_bp4(replace(case(0.01), consequence="missense")).status == Status.UNDETERMINED
