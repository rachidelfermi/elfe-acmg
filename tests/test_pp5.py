import pytest

from adaptive_pp5 import evaluate_pp5


@pytest.mark.parametrize("classification,source", [
    ("PATHOGENIC", "ClinVar"),
    ("LIKELY_PATHOGENIC", "ClinVar"),
    ("CONFLICTING", "ClinVar"),
    ("PATHOGENIC", "expert panel"),
    ("PATHOGENIC", "laboratory"),
])
def test_external_classifications_are_context_only(classification, source):
    result = evaluate_pp5({
        "external_classification": classification,
        "classification_source": source,
        "classification_review_status": "reviewed by expert panel",
    })
    assert (result.status, result.strength) == ("NOT_APPLICABLE", "NONE")
    assert result.external_assertions[0].classification == classification
    assert result.external_assertions[0].source == source
    assert "classification_review_status" not in result.features_used


def test_missing_evidence_and_metadata_are_not_required():
    result = evaluate_pp5({})
    assert result.status == "NOT_APPLICABLE"
    assert result.external_assertions == []
    assert result.missing_critical_features == []
    assert result.guideline["version"] is None


@pytest.mark.parametrize("available", [None, False, True])
def test_evidence_availability_never_activates_pp5(available):
    result = evaluate_pp5({
        "external_classification": "PATHOGENIC",
        "classification_evidence_available": available,
        "classification_evidence_independently_assessable": available,
    })
    assert result.status == "NOT_APPLICABLE"
    assert result.strength == "NONE"
    assert bool(result.recommendations) == (available is True)


def test_multiple_conflicting_assertions_are_preserved():
    result = evaluate_pp5({"external_assertions": [
        {"classification": "PATHOGENIC", "source": "lab A", "conflicts": True},
        {"classification": "BENIGN", "source": "lab B", "conflicts": True},
    ]})
    assert len(result.external_assertions) == 2
    assert result.status == "NOT_APPLICABLE"
    assert any("conflicting" in note for note in result.recommendations)


def test_explicit_historical_policy_is_not_silently_activated():
    result = evaluate_pp5({
        "pp5_usage_policy": "HISTORICAL_ACMG_2015",
        "external_classification": "PATHOGENIC",
        "classification_evidence_available": False,
    })
    assert result.guideline["name"] == "ACMG/AMP 2015"
    assert result.guideline["version"] is None
    assert (result.status, result.strength) == ("UNDETERMINED", "NONE")
    assert result.missing_critical_features == ["verified_source_specific_conditions"]


def test_versioned_applicable_vcep_excludes_pp5():
    result = evaluate_pp5({
        "gene_id": "HGNC:7989", "disease_id": "MONDO:0021060",
        "applicable_vcep_specification": {"id": "GN039", "version": "2.2"},
        "pp5_usage_policy": "HISTORICAL_ACMG_2015",
    })
    assert result.status == "NOT_APPLICABLE"
    assert result.applicable_vcep_specification["pp5"] == "EXCLUDED"


def test_unknown_or_mismatched_vcep_cannot_change_policy():
    for spec, gene in [({"id": "GN999", "version": "1"}, "NRAS"),
                       ({"id": "GN039", "version": "2.2"}, "OTHER")]:
        result = evaluate_pp5({
            "gene_id": gene, "disease_id": "MONDO:0021060",
            "applicable_vcep_specification": spec,
        })
        assert result.status == "NOT_APPLICABLE"
        assert result.applicable_vcep_specification is None
        assert result.policy_conflicts


def test_unverified_permissive_vcep_is_reported_and_not_applied():
    result = evaluate_pp5({"applicable_vcep_specification": {
        "id": "GN999", "version": "1", "pp5": "ALLOWED",
    }})
    assert result.status == "NOT_APPLICABLE"
    assert len(result.policy_conflicts) == 2


def test_result_serializes_for_audit():
    result = evaluate_pp5({
        "variant_id": "test", "variant_normalization": "GRCh38 normalized",
        "external_classification": "PATHOGENIC",
        "classification_provenance": "SCV:test",
        "classification_evidence_available": True,
    }).to_dict()
    assert result["external_assertions"][0]["provenance"] == "SCV:test"
    assert result["strength"] == "NONE"
    assert result["variant"]["variant_normalization"] == "GRCh38 normalized"
    assert "variant_id" not in result["features_used"]


def test_invalid_policy_is_rejected():
    with pytest.raises(ValueError):
        evaluate_pp5({"pp5_usage_policy": "USE_PP5"})
    with pytest.raises(ValueError):
        evaluate_pp5({"guideline_version": "2018"})
