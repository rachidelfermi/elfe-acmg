import pytest

from adaptive_bs3 import evaluate_bs3


def case(**changes):
    assay = dict(assay_id="a1", assay_type="ENZYME", mechanism_compatible=True,
                 variant_class_applicable=True, detects_relevant_defect=True,
                 validation_status="FORMALLY_CALIBRATED", result_class="NORMAL",
                 accepted_normal_interpretation=True, calibration_validated=True,
                 odds_path=0.04)
    assay.update(changes)
    return {"context": {"gene_id": "GENE", "disease_id": "DISEASE", "disease_mechanism": "LOSS_OF_FUNCTION"}, "assays": [assay]}


def test_validated_normal_assay():
    result = evaluate_bs3(case())
    assert (result.status, result.strength) == ("ACTIVATED", "STRONG")
    assert result.to_dict()["guideline"]["rule_id"] == "Brnich-Table-3"


@pytest.mark.parametrize("change,status", [
    ({"result_class": "DAMAGING"}, "NOT_MET"),
    ({"validation_status": "UNVALIDATED"}, "NOT_MET"),
    ({"mechanism_compatible": False}, "NOT_MET"),
    ({"detects_relevant_defect": False}, "NOT_MET"),
    ({"saturation": True}, "NOT_MET"),
    ({"floor_effect": True}, "NOT_MET"),
    ({"batch_effect": True}, "NOT_MET"),
    ({"replicate_consistent": False}, "NOT_MET"),
    ({"result_reliable": False}, "NOT_MET"),
    ({"artifact": True}, "NOT_MET"),
    ({"assay_sensitivity_adequate": False}, "NOT_MET"),
    ({"result_class": "INTERMEDIATE"}, "NOT_MET"),
    ({"accepted_normal_interpretation": False}, "UNDETERMINED"),
])
def test_known_exclusions(change, status):
    assert evaluate_bs3(case(**change)).status == status


@pytest.mark.parametrize("assay_type", ["ENZYME", "CHANNEL", "RECEPTOR_SIGNALING"])
def test_mechanism_specific_assays_are_alternatives(assay_type):
    assert evaluate_bs3(case(assay_type=assay_type)).status == "ACTIVATED"


def test_quantitative_assay_requires_validated_normal_range():
    assert evaluate_bs3(case(result_value=0.9)).status == "UNDETERMINED"
    assert evaluate_bs3(case(result_value=0.9, validated_normal_range=[0.8, 1.2])).status == "ACTIVATED"
    assert evaluate_bs3(case(result_value=0.7, validated_normal_range=[0.8, 1.2])).status == "NOT_MET"


def test_qualitative_assay_requires_accepted_interpretation():
    assert evaluate_bs3(case(accepted_normal_interpretation=None)).status == "UNDETERMINED"


def test_gain_of_function_not_excluded_by_unrelated_enzyme_readout():
    data = case(detectable_mechanisms=["LOSS_OF_FUNCTION"])
    data["context"]["disease_mechanism"] = "GAIN_OF_FUNCTION"
    assert evaluate_bs3(data).status == "NOT_MET"


def test_credible_conflict_prevents_activation():
    data = case()
    data["assays"].append({**data["assays"][0], "assay_id": "a2", "result_class": "DAMAGING"})
    result = evaluate_bs3(data)
    assert (result.status, result.decision_path) == ("UNDETERMINED", "CONFLICTING_STUDIES")


def test_multiple_normal_assays_do_not_stack():
    data = case()
    data["assays"].append({**data["assays"][0], "assay_id": "a2", "odds_path": 0.3})
    assert evaluate_bs3(data).strength == "STRONG"


def test_possible_splice_defect_is_not_cleared_by_protein_assay():
    data = case()
    data["context"]["splice_effect_possible"] = True
    assert evaluate_bs3(data).status == "UNDETERMINED"
