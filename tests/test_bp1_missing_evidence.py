import pytest

from adaptive_bp1 import Evidence, EvidenceState, evaluate_bp1
from test_bp1_eligibility import base


@pytest.mark.parametrize("field", ["gene_id", "disease_id", "variant_consequence",
                                     "exceptions_assessed", "missense_common_mechanism",
                                     "truncating_predominant"])
def test_missing_essential(field):
    r = evaluate_bp1(base(**{field: None}))
    assert r.status == "UNDETERMINED"


def test_explicit_evidence_states():
    assert Evidence.parse({"state": "NOT_APPLICABLE"}).state == EvidenceState.NOT_APPLICABLE
    assert evaluate_bp1(base(truncating_predominant={"state": "CONFLICTING"})).status == "UNDETERMINED"


def test_unknown_mechanism_never_inferred_from_tolerance():
    data = base(truncating_predominant=None, missense_constraint=0.1,
                population_missense_tolerance=True, benign_missense_variant_count=100)
    assert evaluate_bp1(data).status == "UNDETERMINED"
