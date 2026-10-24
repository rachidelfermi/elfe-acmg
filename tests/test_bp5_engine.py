"""Synthetic BP5 fixtures; these are not externally curated patient cases."""

import pytest

from adaptive_bp5 import BP5Features, Evidence, EvidenceState, evaluate_bp5


def case(**changes):
    values = dict(
        candidate_variant_id="candidate-1", candidate_gene_id="GENE1",
        candidate_disease_id="MONDO:0000001", candidate_independent_disease_evidence=False,
        alternative_diagnosis_id="diagnosis-2", alternative_variant_id="alternative-2",
        alternative_variant_classification="PATHOGENIC",
        alternative_variant_classification_source="curated laboratory record",
        alternative_variant_classification_evidence="independent variant-level assessment",
        alternative_variant_classification_independence=True,
        alternative_gene_disease_validity=True,
        alternative_inheritance_compatibility=True,
        patient_clinical_findings=["specific finding"],
        alternative_diagnosis_explains_all_relevant_phenotypes=True,
    )
    values.update(changes)
    return values


def status(**changes):
    return evaluate_bp5(case(**changes))


def test_minimum_sufficient_and_optional_fields():
    r = status()
    assert (r.status, r.strength, r.decision_path) == (
        "ACTIVATED", "SUPPORTING", "ESTABLISHED_ALTERNATIVE_DIAGNOSIS")
    assert r.alternative_diagnosis_id == "diagnosis-2"
    assert "alternative_diagnosis_biochemical_confirmation" not in r.missing_critical_features
    assert status(alternative_diagnosis_functional_confirmation=True,
                  patient_age_of_onset="childhood").to_dict()["status"] == "ACTIVATED"


@pytest.mark.parametrize("changes,expected", [
    ({"alternative_diagnosis_explains_core_phenotypes": True}, "ACTIVATED"),
    ({"alternative_diagnosis_explains_all_relevant_phenotypes": False}, "NOT_MET"),
    ({"alternative_diagnosis_explains_core_phenotypes": False}, "NOT_MET"),
    ({"alternative_inheritance_compatibility": False}, "NOT_MET"),
    ({"alternative_variant_classification": "VUS"}, "NOT_MET"),
    ({"alternative_variant_classification_evidence": None}, "UNDETERMINED"),
    ({"alternative_variant_classification_independence": None}, "UNDETERMINED"),
    ({"alternative_variant_classification_source": None}, "UNDETERMINED"),
    ({"alternative_gene_disease_validity": False}, "NOT_MET"),
    ({"candidate_independent_disease_evidence": True}, "NOT_MET"),
    ({"candidate_independent_disease_evidence": None}, "UNDETERMINED"),
    ({"dual_diagnosis_evidence": True}, "NOT_MET"),
    ({"blended_phenotype_evidence": True}, "NOT_MET"),
    ({"candidate_specific_phenotype_evidence": True}, "NOT_MET"),
    ({"candidate_additional_phenotype_contribution": True}, "NOT_MET"),
    ({"oligogenic_disease_evidence": True}, "NOT_MET"),
    ({"modifier_variant_evidence": True}, "NOT_MET"),
    ({"alternative_diagnosis_conflicts": True}, "UNDETERMINED"),
    ({"clinical_diagnostic_uncertainty": True}, "UNDETERMINED"),
    ({"unexplained_phenotype_features": ["second syndrome"]}, "UNDETERMINED"),
    ({"alternative_diagnosis_confirmation_quality": "LOW"}, "UNDETERMINED"),
    ({"phenotype_evidence_quality": "LOW"}, "UNDETERMINED"),
    ({"candidate_disease_id": None}, "UNDETERMINED"),
    ({"alternative_diagnosis_id": None}, "UNDETERMINED"),
    ({"alternative_variant_id": None}, "UNDETERMINED"),
    ({"patient_clinical_findings": None}, "UNDETERMINED"),
    ({"alternative_diagnosis_explains_all_relevant_phenotypes": None}, "UNDETERMINED"),
    ({"alternative_diagnosis_conflicts": {"state": "CONFLICTING"}}, "UNDETERMINED"),
    ({"candidate_disease_inheritance": "RECESSIVE", "candidate_zygosity": "HET"}, "NOT_MET"),
])
def test_decisive_features(changes, expected):
    assert status(**changes).status == expected


def test_phenotype_only_is_not_a_molecular_diagnosis():
    r = evaluate_bp5({"candidate_variant_id": "v", "candidate_gene_id": "G",
                      "candidate_disease_id": "D", "candidate_independent_disease_evidence": False,
                      "alternative_diagnosis_phenotype_match": "HIGH",
                      "patient_hpo_terms": ["HP:0000001"]})
    assert r.status == "UNDETERMINED"
    assert r.alternative_diagnosis_id is None


def test_clinical_or_biochemical_confirmation_is_not_molecular_confirmation():
    assert status(alternative_variant_id=None,
                  alternative_diagnosis_clinical_confirmation=True,
                  alternative_diagnosis_biochemical_confirmation=True).status == "UNDETERMINED"


def test_environmental_explanation_alone_is_not_bp5():
    assert evaluate_bp5({"candidate_variant_id": "v", "candidate_gene_id": "G",
                         "candidate_disease_id": "D", "candidate_independent_disease_evidence": False,
                         "environmental_disease_explanation": True}).status == "UNDETERMINED"


def test_overlapping_criterion_evidence_is_not_reused():
    r = status(bp5_evidence_ids=["case-observation-1"],
               other_criteria_evidence_ids={"PP4": ["case-observation-1"]})
    assert r.status == "UNDETERMINED"
    assert "CORRELATED_EVIDENCE" in r.exclusions_evaluated


def test_vcep_disabled_and_scope():
    disabled = case(candidate_gene_id="ITGA2B", candidate_disease_id="MONDO:0100326")
    r = evaluate_bp5(disabled)
    assert r.status == "NOT_APPLICABLE"
    assert r.guideline["rule_id"] == "GN011:BP5"
    disabled["candidate_disease_id"] = "OTHER"
    assert evaluate_bp5(disabled).status == "ACTIVATED"
    assert status(vcep_specification_id="GN011").status == "UNDETERMINED"


def test_unimplemented_rasopathy_scoring_does_not_fall_back_to_general_bp5():
    r = status(candidate_gene_id="RRAS2", candidate_disease_id="MONDO:0021060")
    assert r.status == "UNDETERMINED"
    assert r.guideline["rule_id"] == "GN127:BP5"


def brca(**changes):
    values = {"candidate_variant_id": "BRCA1-v", "candidate_gene_id": "BRCA1",
              "candidate_disease_id": "MONDO:0700268", "clinically_calibrated_lr": True,
              "lr_evidence_description": "reviewed multifactorial clinical LR",
              "clinical_likelihood_ratio": 0.4}
    values.update(changes)
    return evaluate_bp5(values)


@pytest.mark.parametrize("lr,strength", [(0.48, "SUPPORTING"), (0.23, "MODERATE"),
                                           (0.05, "STRONG"), (0.001, "STRONG")])
def test_brca1_verified_strengths(lr, strength):
    r = brca(clinical_likelihood_ratio=lr)
    assert (r.status, r.strength, r.decision_path) == ("ACTIVATED", strength, "VCEP_CLINICAL_LR")


def test_brca1_exclusions_and_missing():
    assert brca(breast_ovarian_coobservation_only=True).status == "NOT_APPLICABLE"
    assert brca(clinical_likelihood_ratio=0.5).status == "NOT_MET"
    assert brca(clinical_likelihood_ratio=None).status == "UNDETERMINED"
    assert brca(clinically_calibrated_lr=False).status == "NOT_MET"
    assert brca(vcep_specification_version="0.1").status == "UNDETERMINED"


def test_explicit_states_and_type_validation():
    r = status(alternative_gene_disease_validity=Evidence(state=EvidenceState.CONFLICTING))
    assert r.status == "UNDETERMINED"
    with pytest.raises(ValueError):
        BP5Features.from_mapping({"invented": True})
    with pytest.raises(ValueError):
        Evidence.parse({"state": "UNKNOWN", "value": False})
    with pytest.raises(ValueError):
        brca(clinical_likelihood_ratio=-1)
