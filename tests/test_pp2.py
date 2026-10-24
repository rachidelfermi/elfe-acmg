"""Synthetic decision fixtures; these are not clinically validated variant cases."""

import pytest

from adaptive_pp2 import Status, Strength, evaluate_pp2


def curated(value=True, **changes):
    return {"state": "KNOWN", "value": value, "source": "synthetic expert review",
            "version": "fixture-1", "quality": "HIGH", "independent_of_pp2": True,
            "disease_match": True, **changes}


@pytest.fixture
def general():
    return {
        "variant_id": "synthetic-1", "gene_id": "TEST1", "disease_id": "MONDO:synthetic",
        "variant_consequence": "MISSENSE", "relevant_transcript": True,
        "gene_disease_valid": True, "missense_mechanism_common": curated(),
        "benign_missense_uncommon": curated(),
    }


@pytest.fixture
def acta1():
    return {
        "gene_id": "ACTA1", "disease_id": "MONDO:0100084",
        "variant_consequence": "MISSENSE", "relevant_transcript": True,
        "gene_disease_valid": True, "inheritance_mode": "AD",
        "transcript_id_version": "NM_001100.4",
    }


def test_minimum_sufficient_general(general):
    result = evaluate_pp2(general)
    assert (result.status, result.strength) == (Status.ACTIVATED, Strength.SUPPORTING)
    assert result.missing_critical_features == ()
    assert result.guideline.version == "2015"
    assert result.to_dict()["status"] == "ACTIVATED"


@pytest.mark.parametrize("change,expected", [
    ({"variant_consequence": "SYNONYMOUS"}, Status.NOT_MET),
    ({"gene_disease_valid": False}, Status.NOT_MET),
    ({"relevant_transcript": False}, Status.NOT_MET),
    ({"missense_mechanism_common": curated(False)}, Status.NOT_MET),
    ({"benign_missense_uncommon": curated(False)}, Status.NOT_MET),
    ({"mechanism_compatible": False}, Status.NOT_MET),
    ({"mapping_reliable": False}, Status.NOT_MET),
    ({"alternative_splice_consequence": True}, Status.UNDETERMINED),
    ({"missense_mechanism_common": {"state": "CONFLICTING"}}, Status.UNDETERMINED),
    ({"missense_mechanism_common": curated(True, independent_of_pp2=False)}, Status.UNDETERMINED),
    ({"missense_mechanism_common": curated(True, disease_match=False)}, Status.UNDETERMINED),
    ({"missense_mechanism_common": curated(True, ascertainment_bias=True)}, Status.UNDETERMINED),
    ({"missense_mechanism_common": {"state": "UNKNOWN"}}, Status.UNDETERMINED),
    ({"benign_missense_uncommon": {"state": "UNKNOWN"}}, Status.UNDETERMINED),
    ({"benign_missense_uncommon": curated(True, quality="LOW")}, Status.UNDETERMINED),
])
def test_general_failures(general, change, expected):
    general.update(change)
    assert evaluate_pp2(general).status == expected


def test_optional_fields_do_not_change_decision(general):
    baseline = evaluate_pp2(general)
    general.update({"pathogenic_missense_variant_count": 100, "benign_missense_variant_count": 0,
                    "missense_constraint_z_score": 8.0, "population_missense_variant_count": 1000,
                    "regional_missense_constraint": 0.1, "requested_strength": "STRONG"})
    result = evaluate_pp2(general)
    assert (result.status, result.strength) == (baseline.status, baseline.strength)


def test_constraint_and_counts_do_not_invent_assessments(general):
    general.pop("missense_mechanism_common")
    general.pop("benign_missense_uncommon")
    general.update({"missense_constraint_z_score": 10.0,
                    "pathogenic_missense_fraction": 0.99, "pathogenic_missense_variant_count": 500,
                    "benign_missense_variant_count": 0})
    assert evaluate_pp2(general).status == Status.UNDETERMINED


def test_duplicate_submissions_do_not_create_mechanism(general):
    general.pop("missense_mechanism_common")
    general["pathogenic_missense_variant_count"] = 100
    general["independent_pathogenic_missense_positions"] = 1
    assert evaluate_pp2(general).status == Status.UNDETERMINED


def test_population_variants_are_not_benign_classifications(general):
    general.pop("benign_missense_uncommon")
    general["population_missense_variant_count"] = 300
    assert evaluate_pp2(general).status == Status.UNDETERMINED


def test_acta1_pinned_rule(acta1):
    result = evaluate_pp2(acta1)
    assert (result.status, result.strength) == (Status.ACTIVATED, Strength.SUPPORTING)
    assert result.guideline.version == "2.0"
    assert result.guideline.rule_id == "ACTA1_AD_PP2"


@pytest.mark.parametrize("change,expected", [
    ({"inheritance_mode": "AR"}, Status.NOT_APPLICABLE),
    ({"inheritance_mode": None}, Status.UNDETERMINED),
    ({"transcript_id_version": "NM_001100.3"}, Status.NOT_MET),
    ({"transcript_id_version": None}, Status.UNDETERMINED),
    ({"variant_consequence": "SPLICE_DONOR"}, Status.NOT_MET),
    ({"pp2_specification_id": "other-version"}, Status.UNDETERMINED),
    ({"benign_missense_uncommon": curated(False)}, Status.NOT_MET),
    ({"missense_mechanism_common": {"state": "CONFLICTING"}}, Status.UNDETERMINED),
])
def test_acta1_scope(acta1, change, expected):
    acta1.update(change)
    assert evaluate_pp2(acta1).status == expected


def test_dicer1_vcep_blocks_high_constraint():
    result = evaluate_pp2({"gene_id": "DICER1", "disease_id": "MONDO:0100216",
                           "variant_consequence": "MISSENSE", "missense_constraint_z_score": 4.23})
    assert (result.status, result.strength) == (Status.NOT_APPLICABLE, Strength.NONE)
    assert result.guideline.version == "1.4"


def test_unknown_vcep_cannot_supply_a_threshold(general):
    general["pp2_specification_id"] = "invented-threshold"
    assert evaluate_pp2(general).status == Status.UNDETERMINED


def test_no_arbitrary_strength_upgrade(acta1):
    acta1["requested_strength"] = "STRONG"
    assert evaluate_pp2(acta1).strength == Strength.SUPPORTING


def test_missing_essential_eligibility(general):
    general.pop("relevant_transcript")
    result = evaluate_pp2(general)
    assert result.status == Status.UNDETERMINED
    assert "relevant_transcript" in result.missing_critical_features
