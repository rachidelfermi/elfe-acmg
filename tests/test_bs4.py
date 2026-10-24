"""Synthetic fixtures; no externally curated clinical cases are asserted."""

from copy import deepcopy

import pytest

from adaptive_bs4 import evaluate_bs4


def case():
    return {
        "variant_id": "synthetic-v1", "gene_id": "TEST1", "disease_id": "TEST:D1",
        "gene_disease_validity": True, "inheritance_compatible": True,
        "variant_reliable": True,
        "observations": [{
            "family_id": "f1", "individual_id": "p1", "affected": True,
            "variant_present": False, "genotype_reliable": True,
            "phenotype_reliable": True, "phenotype_consistent": True,
            "biological_relationship_confirmed": True,
            "phenocopy_excluded": True, "heterogeneity_excluded": True,
        }],
    }


def status(data):
    result = evaluate_bs4(data)
    return result.status, result.strength


def test_minimum_evidence_and_irrelevant_features():
    data = case()
    assert status(data) == ("ACTIVATED", "STRONG")
    data["population_frequency"] = 0.3
    data["observations"][0]["age"] = 12
    assert status(data) == ("ACTIVATED", "STRONG")


def test_multiple_noncarriers_and_families():
    data = case()
    second = deepcopy(data["observations"][0])
    second.update(family_id="f2", individual_id="p2")
    data["observations"].append(second)
    result = evaluate_bs4(data)
    assert (result.status, result.informative_nonsegregation_events) == ("ACTIVATED", 2)
    assert result.families_used == ("f1", "f2")


@pytest.mark.parametrize("change,expected", [
    ({"affected": False}, "NOT_MET"),
    ({"variant_present": True}, "NOT_MET"),
    ({"affected": None}, "UNDETERMINED"),
    ({"variant_present": None}, "UNDETERMINED"),
    ({"genotype_reliable": None}, "UNDETERMINED"),
    ({"phenotype_reliable": None}, "UNDETERMINED"),
    ({"biological_relationship_confirmed": None}, "UNDETERMINED"),
    ({"phenocopy_excluded": None}, "UNDETERMINED"),
    ({"heterogeneity_excluded": None}, "UNDETERMINED"),
    ({"phenocopy_excluded": False}, "NOT_MET"),
    ({"heterogeneity_excluded": False}, "NOT_MET"),
    ({"alternative_diagnosis": True}, "NOT_MET"),
    ({"phenocopy_evidence": True}, "NOT_MET"),
    ({"genetic_heterogeneity": True}, "NOT_MET"),
    ({"locus_heterogeneity": True}, "NOT_MET"),
    ({"phenotype_misclassification": True}, "NOT_MET"),
    ({"mosaicism_explanation": True}, "NOT_MET"),
    ({"genotyping_artifact": True}, "NOT_MET"),
    ({"copy_number_conflict": True}, "NOT_MET"),
    ({"phasing_conflict": True}, "NOT_MET"),
    ({"recombination_explanation": True}, "NOT_MET"),
    ({"sample_identity_conflict": True}, "NOT_MET"),
    ({"pedigree_conflict": True}, "NOT_MET"),
])
def test_observation_exclusions_and_missing(change, expected):
    data = case()
    data["observations"][0].update(change)
    assert status(data)[0] == expected


def test_unaffected_carrier_is_not_bs4():
    data = case()
    data["observations"][0].update(affected=False, variant_present=True)
    data["other_criteria_overlap"] = ["BS2"]
    result = evaluate_bs4(data)
    assert result.status == "NOT_MET" and result.overlaps == ("BS2",)


def test_duplicate_report_does_not_upgrade():
    data = case()
    data.update(gene_id="CTLA4", disease_id="MONDO:0014493")
    data["observations"][0]["phenotype_score"] = 6
    data["observations"].append(deepcopy(data["observations"][0]))
    result = evaluate_bs4(data)
    assert (result.informative_nonsegregation_events, result.strength) == (1, "SUPPORTING")


def test_conflicting_duplicate_and_family_overlap():
    data = case()
    duplicate = deepcopy(data["observations"][0])
    duplicate["variant_present"] = True
    data["observations"].append(duplicate)
    assert status(data)[0] == "UNDETERMINED"
    data = case()
    data["family_overlap"] = True
    assert status(data)[0] == "UNDETERMINED"


@pytest.mark.parametrize("gene,disease", [
    ("CTLA4", "MONDO:0014493"),
    ("PIK3R1", "MONDO:1060136"),
])
def test_vcep_count_and_phenotype_rule(gene, disease):
    data = case()
    data.update(gene_id=gene, disease_id=disease)
    assert status(data)[0] == "UNDETERMINED"
    data["observations"][0]["phenotype_score"] = 5
    assert status(data)[0] == "NOT_MET"
    data["observations"][0]["phenotype_score"] = 6
    assert status(data) == ("ACTIVATED", "SUPPORTING")
    second = deepcopy(data["observations"][0])
    second.update(individual_id="p2", family_id="f2")
    data["observations"].append(second)
    assert status(data) == ("ACTIVATED", "STRONG")
    data["disease_id"] = "OTHER:D1"
    assert status(data) == ("ACTIVATED", "STRONG")
    assert evaluate_bs4(data).guideline["rule_id"] == "BS4"


@pytest.mark.parametrize("lr,strength", [
    (0.49, "NONE"), (0.48, "SUPPORTING"), (0.23, "MODERATE"),
    (0.05, "STRONG"), (0.002, "STRONG"),
])
def test_brca2_likelihood_thresholds(lr, strength):
    data = case()
    data.update(gene_id="BRCA2", disease_id="MONDO:0700269",
                observations=[], quantitative_likelihood_ratio=lr,
                quantitative_validated=True, quantitative_model_appropriate=True,
                quantitative_independent=True)
    result = evaluate_bs4(data)
    assert result.strength == strength
    assert result.status == ("NOT_MET" if strength == "NONE" else "ACTIVATED")


def test_brca2_invalid_model_and_missing_lr():
    data = case()
    data.update(gene_id="BRCA2", disease_id="MONDO:0700269", observations=[],
                quantitative_likelihood_ratio=0.01, quantitative_validated=True,
                quantitative_model_appropriate=False, quantitative_independent=True)
    assert status(data)[0] == "NOT_MET"
    data["quantitative_model_appropriate"] = True
    data["quantitative_likelihood_ratio"] = None
    assert status(data)[0] == "UNDETERMINED"


def test_unpinned_vcep_version_is_not_applied():
    data = case()
    data.update(gene_id="CTLA4", disease_id="MONDO:0014493",
                vcep_specification_version="9.9")
    assert status(data)[0] == "UNDETERMINED"


def test_conflicting_pp1_and_explicit_states():
    data = case()
    data["pp1_conflict"] = True
    assert status(data)[0] == "UNDETERMINED"
    data["pp1_conflict"] = False
    data["observations"][0]["variant_present"] = {"state": "CONFLICTING"}
    assert status(data)[0] == "UNDETERMINED"


def test_no_invented_observations_and_context():
    data = case()
    data["observations"] = []
    assert status(data)[0] == "UNDETERMINED"
    data = case()
    data["gene_disease_validity"] = None
    assert status(data)[0] == "UNDETERMINED"
    data["gene_disease_validity"] = False
    assert status(data)[0] == "NOT_APPLICABLE"
