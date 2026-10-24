"""Synthetic PP1 fixtures; no patient data or externally validated classifications."""

from copy import deepcopy

import pytest

from adaptive_pp1 import evaluate_pp1


def case(gene="RUNX1", disease="FPD_AML", mode="AD", count=3):
    return {
        "variant_id": "synthetic-1", "gene_id": gene, "disease_id": disease,
        "gene_disease_validity": "DEFINITIVE", "inheritance_mode": mode,
        "variant_reliable": True, "disease_mechanism_compatible": True,
        "families": [{
            "family_id": "F1", "affected_segregations": count,
            "count_validated": True, "proband_excluded": True,
            "relationships_confirmed": True, "genotype_reliable": True,
            "phenotype_relevant": True,
        }],
    }


def lod(score=0.6):
    return {"score": score, "method": "curated segregation likelihood",
            "disease_model_validated": True, "penetrance_assumptions_validated": True,
            "pedigree_validated": True, "phenotype_model_validated": True}


@pytest.mark.parametrize("count,status,strength", [
    (0, "NOT_MET", "NONE"), (2, "NOT_MET", "NONE"),
    (3, "ACTIVATED", "SUPPORTING"), (5, "ACTIVATED", "MODERATE"),
    (7, "ACTIVATED", "STRONG"),
])
def test_runx1_meioses(count, status, strength):
    result = evaluate_pp1(case(count=count))
    assert (result.status, result.strength) == (status, strength)
    assert result.lod_score is None


def test_pten_two_families_for_strong():
    one = case("PTEN", "PHTS", count=7)
    assert evaluate_pp1(one).strength == "MODERATE"
    two = deepcopy(one)
    two["families"][0].update(independent=True, family_overlap=False, published_family_duplication=False)
    extra = deepcopy(two["families"][0])
    extra.update(family_id="F2", affected_segregations=1)
    two["families"].append(extra)
    assert evaluate_pp1(two).strength == "STRONG"


@pytest.mark.parametrize("issue", ["duplicate", "related", "overlap", "publication"])
def test_multiple_families_must_be_independent(issue):
    data = case(count=2)
    first = data["families"][0]
    first.update(independent=True, family_overlap=False, published_family_duplication=False)
    second = deepcopy(first)
    second["family_id"] = "F2"
    data["families"].append(second)
    assert evaluate_pp1(data).strength == "SUPPORTING"
    if issue == "duplicate":
        second["family_id"] = "F1"
    elif issue == "related":
        second["independent"] = False
    elif issue == "overlap":
        second["family_overlap"] = True
    else:
        second["published_family_duplication"] = True
    assert evaluate_pp1(data).status == "UNDETERMINED"


@pytest.mark.parametrize("field,value,expected", [
    ("affected_noncarrier_count", 1, "UNDETERMINED"),
    ("linked_variant_evidence", True, "UNDETERMINED"),
    ("genotyping_artifact", True, "UNDETERMINED"),
    ("pedigree_conflict", True, "UNDETERMINED"),
    ("sample_identity_conflict", True, "UNDETERMINED"),
])
def test_exclusions(field, value, expected):
    data = case()
    data["families"][0][field] = value
    assert evaluate_pp1(data).status == expected


def test_phenocopy_explains_affected_noncarrier():
    data = case()
    data["families"][0].update(affected_noncarrier_count=1, phenocopy_explained=True)
    assert evaluate_pp1(data).status == "ACTIVATED"
    data["families"][0]["phenocopy_explained"] = False
    assert evaluate_pp1(data).status == "NOT_MET"


def test_unaffected_carrier_penetrance():
    data = case()
    family = data["families"][0]
    family["unaffected_carrier_count"] = 1
    assert evaluate_pp1(data).status == "UNDETERMINED"
    family["unaffected_carriers_young"] = True
    assert evaluate_pp1(data).status == "UNDETERMINED"
    family["penetrance_model_validated"] = True
    assert evaluate_pp1(data).status == "ACTIVATED"
    family["unaffected_carriers_young"] = False
    family.update(reduced_penetrance=True, penetrance_model_validated=True)
    assert evaluate_pp1(data).status == "ACTIVATED"
    family["reduced_penetrance"] = False
    assert evaluate_pp1(data).status == "NOT_MET"


def test_missing_family_genotypes_and_phenotype():
    data = case()
    data["families"][0].pop("affected_segregations")
    assert evaluate_pp1(data).status == "UNDETERMINED"
    data["families"][0]["individuals"] = [{
        "individual_id": "R1", "informative": True, "relationship_confirmed": True,
        "genotype_reliable": True, "carries_variant": True,
    }]
    assert "affected status" in evaluate_pp1(data).missing_critical_features[0]


def test_individuals_count_once_and_exclude_proband():
    data = case()
    family = data["families"][0]
    family.pop("affected_segregations")
    family["individuals"] = [
        {"individual_id": "P", "is_proband": True, "affected": True, "carries_variant": True},
        *[{"individual_id": f"R{i}", "informative": True, "affected": True,
           "carries_variant": True, "relationship_confirmed": True,
           "genotype_reliable": True, "phenotype_consistent": True} for i in range(3)],
    ]
    assert evaluate_pp1(data).informative_meioses == 3
    family["individuals"].append(deepcopy(family["individuals"][1]))
    assert evaluate_pp1(data).status == "UNDETERMINED"


def test_hl_autosomal_dominant_and_recessive():
    ad = case("KCNQ4", "HEARING_LOSS", count=2)
    assert evaluate_pp1(ad).strength == "SUPPORTING"
    ar = case("GJB2", "HEARING_LOSS", mode="AR", count=1)
    ar["families"][0].update(biallelic_phase_confirmed=True, classic_ar_model_validated=True)
    result = evaluate_pp1(ar)
    assert (result.status, result.strength) == ("ACTIVATED", "SUPPORTING")
    assert result.lod_score is not None
    ar["families"][0]["biallelic_phase_confirmed"] = False
    assert evaluate_pp1(ar).status == "UNDETERMINED"


def test_valid_lod_without_raw_genotypes_and_invalid_assumptions():
    data = case("KCNQ4", "HEARING_LOSS")
    data["families"][0].pop("affected_segregations")
    data["families"][0]["lod"] = lod(1.5)
    assert evaluate_pp1(data).strength == "STRONG"
    data["families"][0]["lod"]["penetrance_assumptions_validated"] = False
    assert evaluate_pp1(data).status == "UNDETERMINED"


def test_lod_threshold_is_not_rounded_up():
    data = case("KCNQ4", "HEARING_LOSS")
    data["families"][0].pop("affected_segregations")
    data["families"][0]["lod"] = lod(1.499)
    assert evaluate_pp1(data).strength == "MODERATE"


def test_no_invented_lod_or_ar_risk_model():
    data = case("GJB2", "HEARING_LOSS", mode="AR", count=2)
    data["families"][0]["biallelic_phase_confirmed"] = True
    assert evaluate_pp1(data).status == "UNDETERMINED"
    assert evaluate_pp1(data).lod_score is None


def test_missing_optional_fields_and_irrelevant_values():
    data = case()
    assert evaluate_pp1(data).status == "ACTIVATED"
    data["evidence_source"] = "synthetic fixture"
    data["families"][0]["study_id"] = "test"
    assert evaluate_pp1(data).strength == "SUPPORTING"


def test_no_arbitrary_activation_or_scope_leakage():
    assert evaluate_pp1(case("EXAMPLE", "UNSPECIFIED")).status == "UNDETERMINED"
    assert evaluate_pp1(case("RUNX1", "FPD_AML", mode="AR")).status == "NOT_APPLICABLE"
    assert evaluate_pp1(case("KCNQ4", "HEARING_LOSS", mode="XL")).status == "NOT_APPLICABLE"
    data = case("EXAMPLE", "HEARING_LOSS")
    data["guideline_id"] = "SVI_HL_2018"
    assert evaluate_pp1(data).status == "UNDETERMINED"


def test_alternative_cause_and_overlap():
    data = case()
    data["alternative_causal_variant"] = True
    assert evaluate_pp1(data).status == "NOT_MET"
    data.pop("alternative_causal_variant")
    data["evidence_overlap_prohibited"] = True
    assert evaluate_pp1(data).status == "UNDETERMINED"


def test_missing_essential_and_conflicting_state():
    data = case()
    data.pop("variant_reliable")
    assert evaluate_pp1(data).status == "UNDETERMINED"
    data["variant_reliable"] = True
    data["evidence_states"] = {"genotypes": "CONFLICTING"}
    assert evaluate_pp1(data).status == "UNDETERMINED"
    data["evidence_states"] = {"variant": "UNKNOWN"}
    assert evaluate_pp1(data).status == "UNDETERMINED"
    data["evidence_states"] = {"variant": "NOT_APPLICABLE"}
    assert evaluate_pp1(data).status == "NOT_APPLICABLE"


def test_invalid_optional_lod_does_not_block_valid_count():
    data = case()
    data["families"][0]["lod"] = lod()
    data["families"][0]["lod"]["pedigree_validated"] = False
    assert evaluate_pp1(data).status == "ACTIVATED"


def test_two_valid_independent_lods_combine():
    data = case("KCNQ4", "HEARING_LOSS")
    first = data["families"][0]
    first.pop("affected_segregations")
    first.update(lod=lod(0.61), independent=True, family_overlap=False,
                 published_family_duplication=False)
    second = deepcopy(first)
    second["family_id"] = "F2"
    data["families"].append(second)
    assert evaluate_pp1(data).strength == "MODERATE"


def test_recessive_unaffected_at_risk_counts_but_carrier_parent_does_not():
    data = case("GJB2", "HEARING_LOSS", mode="AR", count=0)
    family = data["families"][0]
    family.update(biallelic_phase_confirmed=True, classic_ar_model_validated=True,
                  unaffected_recessive_segregations=5)
    assert evaluate_pp1(data).strength == "SUPPORTING"
    family["unaffected_recessive_segregations"] = 0
    assert evaluate_pp1(data).status == "NOT_MET"


def test_confirmed_individual_overlap_across_families():
    data = case(count=0)
    first = data["families"][0]
    first.update(independent=True, family_overlap=False, published_family_duplication=False,
                 individuals=[{"individual_id": "R1"}])
    second = deepcopy(first)
    second["family_id"] = "F2"
    data["families"].append(second)
    assert evaluate_pp1(data).status == "UNDETERMINED"
