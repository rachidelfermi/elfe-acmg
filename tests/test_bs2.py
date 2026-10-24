"""Synthetic fixtures: these are rule checks, not clinically curated cases."""

from copy import deepcopy

import pytest

from adaptive_bs2 import evaluate_bs2


def case(mode="AD", genotype="HETEROZYGOUS", gene="TEST1", disease="SYNTHETIC"):
    return {
        "variant_id": "synthetic-variant", "variant_reliable": True,
        "disease": {"gene_id": gene, "disease_id": disease,
                    "gene_disease_validity": True, "inheritance_mode": mode,
                    "molecular_mechanism_compatible": True, "full_penetrance": True,
                    "early_onset": True},
        "observations": [{"individual_id": "person-1", "genotype": genotype,
                          "genotype_reliable": True, "confirmed_unaffected": True,
                          "adult": True, "beyond_expected_onset": True,
                          "clinical_assessment_adequate": True}],
    }


@pytest.mark.parametrize("mode,genotype,status", [
    ("AD", "HETEROZYGOUS", "ACTIVATED"),
    ("AR", "HOMOZYGOUS", "ACTIVATED"),
    ("AR", "HETEROZYGOUS", "NOT_MET"),
    ("XL", "HEMIZYGOUS", "ACTIVATED"),
    ("XL", "HETEROZYGOUS", "NOT_MET"),
])
def test_inheritance(mode, genotype, status):
    result = evaluate_bs2(case(mode, genotype))
    assert result.status == status
    assert result.strength == ("STRONG" if status == "ACTIVATED" else "NONE")


@pytest.mark.parametrize("field,value,status", [
    ("adult", False, "UNDETERMINED"),
    ("beyond_expected_onset", False, "UNDETERMINED"),
    ("confirmed_unaffected", False, "NOT_MET"),
    ("clinical_assessment_adequate", None, "UNDETERMINED"),
    ("genotype_reliable", False, "NOT_MET"),
    ("artifact", True, "NOT_MET"),
    ("deletion_mimicking_homozygosity", True, "NOT_MET"),
    ("mosaicism", True, "NOT_MET"),
    ("protective_modifier", True, "NOT_MET"),
    ("genetic_rescue", True, "NOT_MET"),
    ("subclinical_manifestation", True, "NOT_MET"),
    ("phenotype_misclassification", True, "NOT_MET"),
])
def test_observation_failures(field, value, status):
    data = case()
    data["observations"][0][field] = value
    assert evaluate_bs2(data).status == status


@pytest.mark.parametrize("field,value,status", [
    ("full_penetrance", False, "NOT_MET"),
    ("early_onset", False, "NOT_MET"),
    ("full_penetrance", None, "UNDETERMINED"),
    ("early_onset", None, "UNDETERMINED"),
    ("age_dependent_penetrance", True, "NOT_MET"),
    ("variable_expressivity", True, "UNDETERMINED"),
])
def test_disease_context(field, value, status):
    data = case()
    data["disease"][field] = value
    assert evaluate_bs2(data).status == status


def test_variable_expressivity_assessed():
    data = case()
    data["disease"].update(variable_expressivity=True, relevant_subclinical_features=False)
    assert evaluate_bs2(data).status == "ACTIVATED"


def test_population_genotype_without_clinical_status():
    data = case("AR", "HOMOZYGOUS")
    data["observations"][0].update(source_type="POPULATION_DATABASE", confirmed_unaffected=None)
    result = evaluate_bs2(data)
    assert result.status == "UNDETERMINED"
    assert "person-1.confirmed_unaffected" in result.missing_critical_features


def test_missing_essential_and_optional():
    data = case()
    assert evaluate_bs2(data).status == "ACTIVATED"
    data["observations"][0]["genotype"] = {"state": "CONFLICTING"}
    assert evaluate_bs2(data).status == "UNDETERMINED"


def test_extra_uncertain_observation_does_not_block_complete_path():
    data = case()
    data["observations"].append({"individual_id": "person-2"})
    assert evaluate_bs2(data).status == "ACTIVATED"


def test_duplicate_and_related_general_observations_no_strength_upgrade():
    data = case()
    duplicate = deepcopy(data["observations"][0])
    data["observations"].append(duplicate)
    result = evaluate_bs2(data)
    assert result.status == "ACTIVATED"
    assert result.eligible_healthy_observations == 1
    assert result.strength == "STRONG"
    data["observations"][1]["individual_id"] = "person-2"
    data["observations"][1]["related_group"] = "family-1"
    assert evaluate_bs2(data).strength == "STRONG"


def test_runx1_vcep_exclusion():
    data = case(gene="RUNX1", disease="FPD_AML")
    assert evaluate_bs2(data).status == "NOT_APPLICABLE"


def rett_case(gene="FOXG1", genotype="HETEROZYGOUS"):
    data = case("AD", genotype, gene, "RETT_ANGELMAN_LIKE")
    data["vcep_applicable"] = True
    data["vcep_specification_version"] = "2022 publication"
    return data


def test_rett_count_and_genotype():
    data = rett_case()
    assert evaluate_bs2(data).status == "NOT_MET"
    second = deepcopy(data["observations"][0])
    second["individual_id"] = "person-2"
    data["observations"].append(second)
    assert evaluate_bs2(data).status == "ACTIVATED"
    data["observations"][1]["genotype"] = "HOMOZYGOUS"
    assert evaluate_bs2(data).status == "NOT_MET"


def test_rett_age_and_exclusion():
    data = rett_case()
    second = deepcopy(data["observations"][0])
    second["individual_id"] = "person-2"
    data["observations"].append(second)
    data["observations"][1]["adult"] = False
    assert evaluate_bs2(data).status == "UNDETERMINED"
    data["observations"][1]["adult"] = True
    data["observations"][1]["mosaicism"] = True
    assert evaluate_bs2(data).status == "NOT_MET"


def test_ube3a_maternal_and_four_observations():
    data = rett_case("UBE3A")
    data["observations"] = [deepcopy(data["observations"][0]) for _ in range(4)]
    for i, o in enumerate(data["observations"]):
        o["individual_id"] = f"person-{i}"
        o["maternal_allele"] = True
    assert evaluate_bs2(data).status == "ACTIVATED"
    data["observations"][3]["maternal_allele"] = None
    assert evaluate_bs2(data).status == "UNDETERMINED"
    for o in data["observations"]:
        o["maternal_allele"] = None
        o["independent"] = True
    assert evaluate_bs2(data).status == "ACTIVATED"


def test_vcep_requires_pinned_version_and_scope():
    data = rett_case()
    data["vcep_specification_version"] = None
    assert evaluate_bs2(data).status == "UNDETERMINED"
    data["disease"]["disease_id"] = "OTHER_DISEASE"
    assert evaluate_bs2(data).guideline["name"] == "ACMG/AMP"


def test_explicit_states_and_irrelevant_feature():
    data = case()
    data["observations"][0]["artifact"] = {"state": "NOT_APPLICABLE"}
    assert evaluate_bs2(data).status == "ACTIVATED"
    data["observations"][0]["artifact"] = {"state": "CONFLICTING"}
    assert evaluate_bs2(data).status == "UNDETERMINED"
    assert evaluate_bs2(case()).to_dict()["criterion"] == "BS2"
