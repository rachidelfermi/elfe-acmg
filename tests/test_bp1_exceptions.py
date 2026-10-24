import pytest

from adaptive_bp1 import evaluate_bp1
from test_bp1_eligibility import base


@pytest.mark.parametrize("field", [
    "candidate_is_known_pathogenic_missense", "candidate_region_missense_pathogenicity",
    "candidate_residue_missense_pathogenicity", "missense_specific_disease_subtype",
    "missense_gain_of_function_mechanism", "missense_dominant_negative_mechanism",
])
def test_known_missense_exceptions(field):
    assert evaluate_bp1(base(**{field: True})).status == "NOT_MET"


def test_damaging_splicing_and_prediction():
    assert evaluate_bp1(base(candidate_splice_effect="DAMAGING")).status == "NOT_MET"
    assert evaluate_bp1(base(candidate_splice_prediction="HIGH_RISK")).status == "UNDETERMINED"


def test_conflicting_exception():
    assert evaluate_bp1(base(candidate_residue_missense_pathogenicity={"state": "CONFLICTING"})).status == "UNDETERMINED"


def test_prior_pp2_premise_conflict():
    assert evaluate_bp1(base(prior_criteria=["PP2"])).status == "UNDETERMINED"
