from dataclasses import replace

import pytest

from adaptive_bp4 import BiologicalFeatures, Prediction as RawPrediction, State, Status, Strength, evaluate_bp4


def Prediction(name, score=None, **kwargs):
    checks = dict.fromkeys(("variant_match", "transcript_match", "genome_build_match", "model_applicability", "quality"), State.KNOWN)
    return RawPrediction(name, score, **(checks | kwargs))


def case(revel=0.1, splice=0.05):
    return BiologicalFeatures("v3", "RPE65", "MONDO:0100368", "missense", "NM_000329.3", biologically_relevant_transcript=State.KNOWN, predictions=(Prediction("REVEL", revel), Prediction("SpliceAI", splice)), vcep_version="1.0")


@pytest.mark.parametrize("revel,splice,status,strength", [
    (0.183, 0.1, Status.ACTIVATED, Strength.MODERATE),
    (0.184, 0.1, Status.ACTIVATED, Strength.SUPPORTING),
    (0.290, 0.1, Status.ACTIVATED, Strength.SUPPORTING),
    (0.291, 0.1, Status.NOT_MET, Strength.NONE),
    (0.1, 0.101, Status.NOT_MET, Strength.NONE),
])
def test_rpe65_joint_rule(revel, splice, status, strength):
    result = evaluate_bp4(case(revel, splice))
    assert (result.status, result.strength) == (status, strength)


def test_rpe65_scope_and_required_evidence():
    assert evaluate_bp4(replace(case(), disease_id="other")).status == Status.UNDETERMINED
    assert evaluate_bp4(replace(case(), vcep_version="2.0")).status == Status.UNDETERMINED
    assert evaluate_bp4(replace(case(), predictions=(Prediction("REVEL", 0.1),))).status == Status.UNDETERMINED
    assert evaluate_bp4(replace(case(), transcript_id_version="NM_000329.2")).status == Status.UNDETERMINED


def test_rpe65_splice_only():
    f = replace(case(), consequence="noncanonical_splice", protein_effect_relevant=False, splice_region_eligible=True, predictions=(Prediction("SpliceAI", 0.05),))
    assert evaluate_bp4(f).strength == Strength.SUPPORTING


def test_palb2_missense_prohibition_and_splice_rule():
    f = BiologicalFeatures("v4", "PALB2", "MONDO:0700272", "missense", "NM_024675.3", biologically_relevant_transcript=State.KNOWN, predictions=(Prediction("REVEL", 0.001),), vcep_version="1.3")
    assert evaluate_bp4(f).status == Status.NOT_APPLICABLE
    splice = replace(f, consequence="noncanonical_splice", protein_effect_relevant=False, splice_region_eligible=True, predictions=(Prediction("SpliceAI", 0.1),))
    assert evaluate_bp4(splice).strength == Strength.SUPPORTING
