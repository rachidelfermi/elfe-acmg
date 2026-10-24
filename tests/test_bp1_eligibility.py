"""All scenarios use synthetic biological features unless a pinned VCEP is named."""

import pytest

from adaptive_bp1 import evaluate_bp1


def base(**changes):
    values = dict(gene_id="TEST1", disease_id="SYNTHETIC:1", variant_consequence="MISSENSE",
                  biologically_relevant_transcript=True, gene_disease_validity=True,
                  exceptions_assessed=True, missense_common_mechanism=False,
                  truncating_predominant=True, mechanism_evidence_quality="CURATED",
                  mechanism_evidence_source="synthetic curated assessment")
    values.update(changes)
    return values


def test_minimum_sufficient_and_optional_missing():
    r = evaluate_bp1(base())
    assert (r.status, r.strength, r.decision_path) == ("ACTIVATED", "SUPPORTING", "TRUNCATING_DISEASE_MECHANISM")
    assert r.missing_critical_features == ()
    assert r.to_dict()["criterion"] == "BP1"


@pytest.mark.parametrize("change,status", [
    ({"variant_consequence": "SYNONYMOUS"}, "NOT_MET"),
    ({"biologically_relevant_transcript": False}, "NOT_MET"),
    ({"biologically_relevant_transcript": None}, "UNDETERMINED"),
    ({"gene_disease_validity": False}, "NOT_MET"),
    ({"gene_disease_validity": None}, "UNDETERMINED"),
    ({"missense_consequence_confidence": "LOW"}, "UNDETERMINED"),
])
def test_eligibility(change, status):
    assert evaluate_bp1(base(**change)).status == status


def test_irrelevant_context_does_not_change_decision():
    original = evaluate_bp1(base())
    extra = evaluate_bp1(base(pathogenic_truncating_variant_count=4,
                              benign_missense_variant_count=200, missense_constraint=0.2,
                              prior_criteria=["PP3"]))
    assert (extra.status, extra.strength, extra.decision_path) == (original.status, original.strength, original.decision_path)
