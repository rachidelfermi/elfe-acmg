from dataclasses import replace

from adaptive_bp4 import BiologicalFeatures, Prediction, State, Status, evaluate_bp4
from adaptive_bp4.evaluator import _valid
from adaptive_bp4.guidelines import REVEL, Rule


def test_explicit_version_pinning_rejects_mismatch():
    pinned = replace(REVEL, allowed_versions=("1.0",))
    prediction = Prediction("REVEL", 0.1, version="2.0", variant_match=State.KNOWN, transcript_match=State.KNOWN, genome_build_match=State.KNOWN, model_applicability=State.KNOWN, quality=State.KNOWN)
    assert _valid(prediction, pinned) == (False, "predictor version does not match calibration", True)


def test_no_invented_activation_and_unrelated_predictors():
    f = BiologicalFeatures("v", "G", consequence="missense", transcript_id_version="NM_1.1")
    assert evaluate_bp4(f).status == Status.UNDETERMINED
    assert evaluate_bp4(replace(f, predictions=(Prediction("unknown", 0.0),))).status == Status.UNDETERMINED
