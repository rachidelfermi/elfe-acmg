"""BP6 policy and provenance regression tests."""

import pytest

from adaptive_bp6 import ExternalAssertion, evaluate_bp6


@pytest.mark.parametrize(
    "assertion",
    [
        {"classification": "BENIGN", "source": "ClinVar", "review_status": "criteria provided, single submitter"},
        {"classification": "LIKELY_BENIGN", "source": "ClinVar", "review_status": "reviewed by expert panel"},
        {"classification": "BENIGN", "source": "expert panel", "review_status": "reviewed by expert panel"},
        {"classification": "BENIGN", "source": "clinical laboratory"},
        {"classification": "CONFLICTING", "source": "ClinVar", "conflicts": True},
        {"classification": "BENIGN", "source": "ClinVar", "evidence_available": False},
        {"classification": "BENIGN", "source": "ClinVar", "evidence_available": True,
         "evidence_independently_assessable": True},
        {"classification": "BENIGN", "source": "ClinVar", "review_status": None,
         "evidence_available": None, "provenance": None},
    ],
)
def test_all_external_assertions_are_context_only(assertion):
    result = evaluate_bp6({"external_assertions": [assertion]})
    assert result.status == "NOT_APPLICABLE"
    assert result.strength == "NONE"
    assert result.external_assertions[0].classification == assertion["classification"]
    assert result.external_assertions[0].source == assertion["source"]
    assert "external_classification" not in result.features_used
    assert result.missing_critical_features == []


def test_current_policy_needs_no_external_evidence():
    result = evaluate_bp6()
    assert result.to_dict()["criterion"] == "BP6"
    assert result.status == "NOT_APPLICABLE"
    assert result.strength == "NONE"
    assert result.guideline["version"] is None
    assert result.external_assertions == []
    assert result.missing_critical_features == []


def test_flat_assertion_preserves_provenance_and_variant_identity():
    result = evaluate_bp6({
        "variant_id": "v1", "normalized_variant_representation": "NC_000023.11:g.1A>G",
        "external_classification": "BENIGN", "classification_source": "ClinVar",
        "classification_review_status": "two stars", "classification_provenance": "submission 123",
    })
    assert result.variant["variant_id"] == "v1"
    assert result.external_assertions[0].provenance == "submission 123"
    assert result.external_assertions[0].review_status == "two stars"
    assert result.strength == "NONE"
    assert result.features_used == ["bp6_usage_policy"]


def test_conflicting_classifications_remain_separate_and_unscored():
    assertions = [ExternalAssertion("BENIGN", "lab A"), ExternalAssertion("PATHOGENIC", "lab B", conflicts=True)]
    result = evaluate_bp6({"external_assertions": assertions})
    assert len(result.external_assertions) == 2
    assert result.strength == "NONE"
    assert any("conflict" in text.lower() for text in result.recommendations)


def test_independently_assessable_evidence_is_routed_elsewhere():
    result = evaluate_bp6({"external_assertions": [{
        "classification": "BENIGN", "source": "ClinVar",
        "evidence_independently_assessable": True,
    }]})
    assert result.strength == "NONE"
    assert any("appropriate ACMG benign criteria" in text for text in result.recommendations)


def test_historical_policy_is_explicit_and_audit_only():
    result = evaluate_bp6({"bp6_usage_policy": "HISTORICAL_ACMG_2015",
                           "external_classification": "BENIGN",
                           "classification_evidence_available": False})
    assert result.status == "UNDETERMINED"
    assert result.strength == "NONE"
    assert result.guideline["name"] == "ACMG/AMP 2015"
    assert result.missing_critical_features == ["verified_historical_source_conditions"]


def test_applicable_vcep_excludes_bp6_even_under_historical_selection():
    result = evaluate_bp6({
        "bp6_usage_policy": "HISTORICAL_ACMG_2015", "gene_id": "HGNC:3551",
        "disease_id": "MONDO:0010604",
        "applicable_vcep_specification": {"id": "GN080", "version": "2.1"},
    })
    assert result.status == "NOT_APPLICABLE"
    assert result.strength == "NONE"
    assert result.applicable_vcep_specification["bp6"] == "EXCLUDED"
    assert result.missing_critical_features == []


@pytest.mark.parametrize("spec,gene,disease", [
    ({"id": "GN999", "version": "1"}, "HGNC:3551", "MONDO:0010604"),
    ({"id": "GN080"}, "HGNC:3551", "MONDO:0010604"),
    ({"id": "GN080", "version": "2.1"}, "HGNC:3551", "MONDO:9999999"),
])
def test_unverified_or_mismatched_vcep_never_changes_policy(spec, gene, disease):
    result = evaluate_bp6({"gene_id": gene, "disease_id": disease,
                           "applicable_vcep_specification": spec})
    assert result.status == "NOT_APPLICABLE"
    assert result.applicable_vcep_specification is None
    assert result.policy_conflicts


def test_unverified_claimed_exception_is_reported_not_applied():
    result = evaluate_bp6({"applicable_vcep_specification": {
        "id": "unknown", "version": "1", "bp6": "ALLOWED"}})
    assert result.status == "NOT_APPLICABLE"
    assert any("conflicts" in item for item in result.policy_conflicts)


def test_invalid_policy_and_invented_version_are_rejected():
    with pytest.raises(ValueError):
        evaluate_bp6({"bp6_usage_policy": "ALLOW"})
    with pytest.raises(ValueError):
        evaluate_bp6({"guideline_version": "2025"})
