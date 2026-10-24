"""Synthetic fixtures exercising source-defined rules; not validated patient cases."""
from dataclasses import replace

import pytest

from adaptive_pp4 import Evidence, PP4Features, State, evaluate_pp4, known


def base(**changes):
    f = PP4Features(gene_id="GENE1", disease_id="MONDO:synthetic",
                    gene_disease_validity=known("DEFINITIVE"),
                    single_genetic_etiology=known(True),
                    phenotype_highly_specific=known(True))
    return replace(f, **changes)


def ada(**changes):
    return replace(base(gene_id="ADA", disease_id="MONDO:0007064",
                        single_genetic_etiology=Evidence(), phenotype_highly_specific=Evidence()), **changes)


def gt(**changes):
    return replace(base(gene_id="ITGA2B", disease_id="MONDO:0100326",
                        single_genetic_etiology=Evidence(), phenotype_highly_specific=Evidence(),
                        gt_clinical_diagnosis=known(True), gt_mucocutaneous_bleeding=known(True),
                        gt_normal_ristocetin=known(True),
                        gt_minimal_aggregation_at_least_two_agonists=known(True)), **changes)


def test_general_minimum_sufficient_evidence():
    r = evaluate_pp4(base())
    assert (r.status, r.strength, r.decision_path) == ("ACTIVATED", "SUPPORTING", "HIGHLY_SPECIFIC_PHENOTYPE")
    assert r.guideline.rule_id == "PP4"


@pytest.mark.parametrize("change,status", [
    ({"phenotype_highly_specific": known(False), "family_history_highly_specific": known(False)}, "NOT_MET"),
    ({"single_genetic_etiology": known(False)}, "NOT_MET"),
    ({"single_genetic_etiology": Evidence()}, "UNDETERMINED"),
    ({"gene_disease_validity": known("REFUTED")}, "NOT_MET"),
    ({"alternative_diagnosis_confirmed": known(True)}, "NOT_MET"),
    ({"phenotype_conflict": known(True)}, "UNDETERMINED"),
    ({"phenotype_conflict": Evidence(State.CONFLICTING)}, "UNDETERMINED"),
    ({"age_sex_applicable": known(False)}, "NOT_MET"),
    ({"evidence_reliable": known(False)}, "UNDETERMINED"),
])
def test_general_blockers(change, status):
    assert evaluate_pp4(base(**change)).status == status


def test_hpo_and_rank_never_activate_without_specificity():
    f = base(phenotype_highly_specific=Evidence(), patient_hpo_terms=("HP:0000001",),
             hpo_similarity_score=1.0, phenotype_match_rank=1,
             phenotype_similarity_method="Exomiser")
    assert evaluate_pp4(f).status == "UNDETERMINED"
    assert evaluate_pp4(replace(f, phenotype_highly_specific=known(False),
                                family_history_highly_specific=known(False))).status == "NOT_MET"


def test_missing_phenotype_and_clinical_diagnosis():
    assert evaluate_pp4(base(phenotype_highly_specific=Evidence())).status == "UNDETERMINED"
    assert evaluate_pp4(base(patient_hpo_terms=())).status == "ACTIVATED"


def test_family_history_without_reusing_pp1():
    f = base(phenotype_highly_specific=Evidence(), family_history_highly_specific=known(True))
    assert evaluate_pp4(f).decision_path == "HIGHLY_SPECIFIC_FAMILY_HISTORY"
    assert evaluate_pp4(replace(f, family_evidence_used_for_pp1=known(True))).status == "UNDETERMINED"


def test_ps2_pm6_overlap_is_not_universal_exclusion():
    assert evaluate_pp4(base(phenotype_evidence_used_for_ps2_pm6=known(True))).status == "ACTIVATED"


def test_irrelevant_features_do_not_change_result():
    original = evaluate_pp4(base())
    added = evaluate_pp4(base(patient_hpo_terms=("HP:1234567",), patient_age_years=3,
                              hpo_similarity_score=.99, phenotype_match_rank=1))
    assert (original.status, original.strength, original.decision_path) == (added.status, added.strength, added.decision_path)


def test_ada_single_verified_component_is_enough():
    r = evaluate_pp4(ada(scid_sequencing_no_alternative=known(True)))
    assert (r.status, r.strength, r.guideline.version) == ("ACTIVATED", "SUPPORTING", "2.2")


@pytest.mark.parametrize("features,strength", [
    ({"ada_enzyme_activity_percent": known(0.5)}, "MODERATE"),
    ({"ada_dado_above_reference": known(True), "ada_dado_pre_treatment_or_nontransfused": known(True)}, "MODERATE"),
    ({"ada_enzyme_activity_percent": known(0.5), "scid_sequencing_no_alternative": known(True)}, "STRONG"),
    ({"ada_phenotype_corrected_by_treatment": known(True), "ada_cnv_testing_done": known(True)}, "STRONG"),
    ({"ada_phenotype_corrected_by_treatment": known(True)}, "MODERATE"),
])
def test_ada_strengths(features, strength):
    assert evaluate_pp4(ada(**features)).strength == strength


def test_ada_biomarker_below_and_missing():
    assert evaluate_pp4(ada(ada_enzyme_activity_percent=known(2.0))).status == "UNDETERMINED"
    assert evaluate_pp4(ada()).status == "UNDETERMINED"
    assert evaluate_pp4(ada(ada_dado_above_reference=known(True))).status == "UNDETERMINED"


def test_ada_family_requires_sequencing_and_no_pp1_overlap():
    f = ada(scid_family_history=known(True), scid_sequencing_no_alternative=known(True))
    assert evaluate_pp4(f).strength == "SUPPORTING"
    assert "scid_family_history" in evaluate_pp4(f).features_used
    assert "scid_family_history" not in evaluate_pp4(replace(f, family_evidence_used_for_pp1=known(True))).features_used


def test_ada_scope_does_not_transfer_to_other_disease():
    r = evaluate_pp4(replace(ada(ada_enzyme_activity_percent=known(0.2)), disease_id="MONDO:other"))
    assert r.status == "UNDETERMINED" and r.guideline.rule_id == "PP4"


def test_gt_vcep_moderate_and_strong():
    assert evaluate_pp4(gt()).strength == "MODERATE"
    assert evaluate_pp4(gt(gt_aIIbbeta3_expression_below_25_percent=known(True),
                           gt_both_genes_fully_sequenced=known(True))).strength == "STRONG"
    assert evaluate_pp4(gt(gt_impaired_fibrinogen_binding=known(True))).strength == "MODERATE"


def test_gt_missing_diagnosis_and_failed_lab():
    assert evaluate_pp4(gt(gt_clinical_diagnosis=Evidence())).status == "UNDETERMINED"
    assert evaluate_pp4(gt(gt_normal_ristocetin=known(False))).status == "NOT_MET"
    assert evaluate_pp4(gt(gt_mucocutaneous_bleeding=known(False))).status == "NOT_MET"


def test_gt_scope_is_gene_and_disease_specific():
    assert evaluate_pp4(replace(gt(), disease_id="MONDO:other")).status == "UNDETERMINED"
    assert evaluate_pp4(replace(gt(), gene_id="OTHER")).status == "UNDETERMINED"


def test_cdh1_vcep_prohibits_pp4():
    f = base(gene_id="CDH1", disease_id="MONDO:0007648")
    r = evaluate_pp4(f)
    assert (r.status, r.guideline.version) == ("NOT_APPLICABLE", "3.1")


def test_invalid_input_rejected():
    with pytest.raises(ValueError):
        Evidence(State.KNOWN)
    with pytest.raises(TypeError):
        evaluate_pp4({})
