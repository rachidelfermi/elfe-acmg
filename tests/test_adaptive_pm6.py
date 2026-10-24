"""Synthetic fixtures derived from the cited rule tables, not validated cases."""

from copy import deepcopy

import pytest

from adaptive_pm6 import evaluate_pm6


def base(**changes):
    data = {"variant_id": "synthetic-1", "gene_id": "TEST", "disease_id": "TEST:1",
            "gene_disease_valid": True, "inheritance_mode": "AUTOSOMAL_DOMINANT",
            "proband_variant_present": True, "maternal_genotype": "ABSENT",
            "paternal_genotype": "ABSENT", "maternal_test_adequate": True,
            "paternal_test_adequate": True, "affected": True,
            "phenotype": "HIGHLY_SPECIFIC", "case_id": "case-1", "family_id": "family-1"}
    data.update(changes)
    return data


@pytest.mark.parametrize("phenotype,points,strength", [
    ("HIGHLY_SPECIFIC", 1.0, "MODERATE"),
    ("CONSISTENT", 0.5, "SUPPORTING"),
    ("HETEROGENEOUS", 0.25, "NONE"),
])
def test_single_assumed_table(phenotype, points, strength):
    result = evaluate_pm6(base(phenotype=phenotype)).to_dict()
    assert result["total_points"] == points
    assert result["strength"] == strength
    assert result["status"] == ("NOT_MET" if strength == "NONE" else "ACTIVATED")


def test_confirmed_and_mixed_delegate_to_ps2():
    confirmed = base(maternity_confirmed=True, paternity_confirmed=True)
    result = evaluate_pm6(confirmed).to_dict()
    assert (result["status"], result["delegation"], result["total_points"]) == ("NOT_APPLICABLE", "PS2", 2)
    second = base(case_id="case-2", family_id="family-2", phenotype="CONSISTENT")
    mixed = evaluate_pm6({**base(), "observations": [confirmed, second]}).to_dict()
    assert (mixed["delegation"], mixed["total_points"], mixed["eligible_observations"]) == ("PS2", 2.5, 2)


@pytest.mark.parametrize("missing", ["maternal_genotype", "paternal_genotype", "proband_variant_present", "gene_disease_valid"])
def test_missing_critical_evidence(missing):
    data = base()
    del data[missing]
    assert evaluate_pm6(data).status.value == "UNDETERMINED"


def test_parentage_unconfirmed_is_assumed_and_disproved_is_excluded():
    assert evaluate_pm6(base(maternity_confirmation_status="UNCONFIRMED", paternity_confirmation_status="UNCONFIRMED")).status.value == "ACTIVATED"
    assert evaluate_pm6(base(maternity_confirmation_status="DISPROVED")).status.value == "NOT_MET"


@pytest.mark.parametrize("change", [
    {"maternal_genotype": "PRESENT"}, {"parental_mosaicism": True},
    {"maternal_test_adequate": False}, {"proband_call_adequate": False},
    {"sample_identity_valid": False}, {"technical_artifact": True},
    {"conflicting_segregation": True}, {"phenotype": "INCONSISTENT"},
    {"affected": False}, {"alternative_diagnosis": True},
    {"already_counted_ps2": True},
])
def test_explicit_exclusions(change):
    assert evaluate_pm6(base(**change)).status.value == "NOT_MET"


def test_raw_depth_zero_blocks_but_no_arbitrary_positive_cutoff():
    assert evaluate_pm6(base(maternal_read_depth=0)).status.value == "NOT_MET"
    assert evaluate_pm6(base(maternal_read_depth=3)).status.value == "ACTIVATED"


def test_proband_mosaicism_can_be_scored_when_presence_is_reliable():
    assert evaluate_pm6(base(proband_mosaicism=True)).status.value == "ACTIVATED"


def test_multiple_independent_cases_and_boundaries():
    cases = [base(case_id=f"c{i}", family_id=f"f{i}", phenotype="CONSISTENT") for i in range(8)]
    expected = [(1, "SUPPORTING"), (2, "MODERATE"), (4, "STRONG"), (8, "VERY_STRONG")]
    for count, strength in expected:
        assert evaluate_pm6({**base(), "observations": cases[:count]}).strength.value == strength


def test_heterogeneous_cap_and_duplicate_family():
    cases = [base(case_id=f"c{i}", family_id=f"f{i}", phenotype="HETEROGENEOUS") for i in range(8)]
    assert evaluate_pm6({**base(), "observations": cases}).total_points == 1
    cases[1]["family_id"] = cases[0]["family_id"]
    result = evaluate_pm6({**base(), "observations": cases[:2]})
    assert result.eligible_observations == 1


def test_unknown_case_independence_cannot_multiply_evidence():
    cases = [base(case_id=f"c{i}", family_id=None) for i in range(2)]
    result = evaluate_pm6({**base(), "observations": cases})
    assert result.status.value == "UNDETERMINED"


def test_recessive_strength_reduction_and_germline_mosaicism():
    data = base(inheritance_mode="AUTOSOMAL_RECESSIVE", second_pathogenic_allele=False)
    assert evaluate_pm6(data).strength.value == "SUPPORTING"
    del data["second_pathogenic_allele"]
    assert evaluate_pm6(data).status.value == "UNDETERMINED"
    assert evaluate_pm6(base(apparent_germline_mosaicism=True)).status.value == "NOT_MET"


def test_vcep_scope_and_ps2_reassignment():
    data = base(gene_id="HNF1A", disease_id="MONDO:0015967",
                vcep_facts={"mody_probability_at_least_50_percent": True, "hnf4a_negative": True})
    result = evaluate_pm6(data).to_dict()
    assert (result["status"], result["delegation"], result["total_points"]) == ("NOT_APPLICABLE", "PS2", 1)
    data["vcep_facts"]["autoimmune_diabetes"] = True
    assert evaluate_pm6(data).total_points == 0
    data["gene_id"] = "OTHER"
    data["vcep_profile"] = "HNF1A-v3.1"
    assert evaluate_pm6(data).status.value == "UNDETERMINED"


def test_slc9a6_verified_point_modification_and_strength_threshold():
    data = base(gene_id="SLC9A6", disease_id="MONDO:0010278", phenotype="CONSISTENT",
                vcep_facts={"neurodevelopmental_phenotype_consistent": True})
    single = evaluate_pm6(data)
    assert (single.total_points, single.strength.value, single.guideline["version"]) == (1, "MODERATE", "6.0")
    second = base(case_id="case-2", family_id="family-2", phenotype="CONSISTENT",
                  vcep_facts={"neurodevelopmental_phenotype_consistent": True})
    assert evaluate_pm6({**data, "observations": [data, second]}).strength.value == "STRONG"
    del data["vcep_facts"]
    assert evaluate_pm6(data).status.value == "UNDETERMINED"


def test_evidence_states_and_no_optional_field_requirement():
    data = base(maternal_genotype={"state": "UNKNOWN"})
    assert evaluate_pm6(data).status.value == "UNDETERMINED"
    data = base()
    assert evaluate_pm6(data).status.value == "ACTIVATED"
    assert evaluate_pm6(deepcopy(data)).to_dict()["missing_critical_features"] == []
    data["maternal_genotype"] = {"state": "CONFLICTING"}
    assert evaluate_pm6(data).status.value == "UNDETERMINED"
    data = base(phenotype={"state": "CONFLICTING"})
    assert evaluate_pm6(data).status.value == "UNDETERMINED"
