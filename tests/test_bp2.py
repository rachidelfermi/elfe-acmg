"""Synthetic BP2 fixtures; these are rule tests, not validated patient cases."""

from dataclasses import replace

import pytest

from adaptive_bp2 import (BiologicalFeatures, Classification, Comparison, DecisionStatus,
                          Finding, Inheritance, Observation, Phase, evaluate_bp2)

Y, N = Finding.yes(), Finding.no()


def fixture(phase=Phase.CIS, *, gene="GENE1", disease="MONDO:TEST", inheritance=Inheritance.AR,
            classification=Classification.PATHOGENIC, case="case-1", comparison_id="ref-1"):
    comparison = Comparison(comparison_id, gene, classification, "curated laboratory",
                            classification_evidence="independent case and functional evidence",
                            classification_independent=Y, observed_independently=Y,
                            disease_match=Y, mechanism_compatible=Y, vcep_classified=Y,
                            genotype="HET")
    obs = Observation(comparison, phase, "long-read sequencing" if phase in (Phase.CIS, Phase.TRANS) else None,
                      Y if phase in (Phase.CIS, Phase.TRANS) else Finding(), candidate_genotype="HET",
                      case_id=case, affected=Y, diagnosis_matches=Y, comparison_explains_disease=Y)
    return BiologicalFeatures("candidate", gene, disease, inheritance, Y, (obs,),
                              synergy_or_complex_allele=N, haplotype_specific_pathogenicity=N,
                              candidate_independent_disease_effect=N, dual_diagnosis=N,
                              full_penetrance=Y)


def change_obs(f, **changes):
    return replace(f, observations=(replace(f.observations[0], **changes),))


def change_comp(f, **changes):
    obs = f.observations[0]
    return change_obs(f, comparison=replace(obs.comparison, **changes))


def test_confirmed_cis_minimum():
    result = evaluate_bp2(fixture())
    assert (result.status, result.strength, result.decision_path) == (
        DecisionStatus.ACTIVATED, "SUPPORTING", "CONFIRMED_CIS_WITH_PATHOGENIC")


def test_confirmed_trans_dominant():
    f = fixture(Phase.TRANS, inheritance=Inheritance.AD)
    assert evaluate_bp2(f).status == DecisionStatus.ACTIVATED


def test_trans_recessive_is_not_bp2():
    f = fixture(Phase.TRANS)
    assert evaluate_bp2(f).status == DecisionStatus.NOT_MET


@pytest.mark.parametrize("phase", [Phase.UNKNOWN, Phase.CONFLICTING])
def test_unresolved_phase(phase):
    assert evaluate_bp2(fixture(phase)).status == DecisionStatus.UNDETERMINED


@pytest.mark.parametrize("classification,status", [
    (Classification.PATHOGENIC, DecisionStatus.ACTIVATED),
    (Classification.LIKELY_PATHOGENIC, DecisionStatus.NOT_MET),
    (Classification.VUS, DecisionStatus.NOT_MET),
    (Classification.BENIGN, DecisionStatus.NOT_MET),
    (Classification.UNKNOWN, DecisionStatus.UNDETERMINED),
])
def test_comparison_classification(classification, status):
    assert evaluate_bp2(fixture(classification=classification)).status == status


def test_circular_classification():
    assert evaluate_bp2(change_comp(fixture(), classification_independent=N)).status == DecisionStatus.NOT_MET


def test_comparison_not_independently_observed():
    assert evaluate_bp2(change_comp(fixture(), observed_independently=N)).status == DecisionStatus.NOT_MET


@pytest.mark.parametrize("changes", [
    {"gene_id": "OTHER"}, {"disease_match": N}, {"mechanism_compatible": N},
])
def test_comparison_context_mismatch(changes):
    assert evaluate_bp2(change_comp(fixture(), **changes)).status == DecisionStatus.NOT_MET


@pytest.mark.parametrize("field", ["synergy_or_complex_allele", "haplotype_specific_pathogenicity",
                                   "candidate_independent_disease_effect", "dual_diagnosis"])
def test_known_exclusions(field):
    assert evaluate_bp2(replace(fixture(), **{field: Y})).status == DecisionStatus.NOT_MET


@pytest.mark.parametrize("penetrance,status", [(N, DecisionStatus.NOT_MET),
                                               (Finding(), DecisionStatus.UNDETERMINED)])
def test_trans_penetrance(penetrance, status):
    f = replace(fixture(Phase.TRANS, inheritance=Inheritance.AD), full_penetrance=penetrance)
    assert evaluate_bp2(f).status == status


def test_molecular_phase_does_not_need_parental_genotypes():
    assert evaluate_bp2(fixture()).status == DecisionStatus.ACTIVATED


def test_missing_phase_method():
    assert evaluate_bp2(change_obs(fixture(), phase_method=None)).status == DecisionStatus.UNDETERMINED


def test_unreliable_call_blocks():
    assert evaluate_bp2(replace(fixture(), candidate_call_reliable=N)).status == DecisionStatus.NOT_MET


def test_optional_features_do_not_change_decision():
    f = fixture()
    enriched = replace(f, candidate_call_reliable=Y, candidate_allele_balance_reliable=Y)
    assert evaluate_bp2(f).status == evaluate_bp2(enriched).status == DecisionStatus.ACTIVATED


def test_missing_essential_context():
    assert evaluate_bp2(replace(fixture(), gene_disease_valid=Finding())).status == DecisionStatus.UNDETERMINED


def test_phenotype_not_required_for_base_cis():
    f = change_obs(fixture(), affected=Finding(), diagnosis_matches=Finding())
    assert evaluate_bp2(f).status == DecisionStatus.ACTIVATED


def test_dicer1_single_trans():
    f = fixture(Phase.TRANS, gene="DICER1", disease="MONDO:0100216",
                inheritance=Inheritance.AD, classification=Classification.LIKELY_PATHOGENIC)
    f = replace(f, full_penetrance=Finding())
    result = evaluate_bp2(f)
    assert result.status == DecisionStatus.ACTIVATED
    assert result.guideline.version == "1.4"


def test_dicer1_cis_requires_three_cases_two_comparisons():
    f = fixture(gene="DICER1", disease="MONDO:0100216", inheritance=Inheritance.AD)
    a = f.observations[0]
    two = replace(a, case_id="case-2")
    three = replace(a, case_id="case-3", comparison=replace(a.comparison, variant_id="ref-2"))
    assert evaluate_bp2(replace(f, observations=(a, two))).status == DecisionStatus.UNDETERMINED
    assert evaluate_bp2(replace(f, observations=(a, two, three))).status == DecisionStatus.ACTIVATED


def test_dicer1_unknown_phase_exception():
    f = fixture(Phase.UNKNOWN, gene="DICER1", disease="MONDO:0100216", inheritance=Inheritance.AD)
    a = f.observations[0]
    b = replace(a, case_id="case-2")
    c = replace(a, case_id="case-3", comparison=replace(a.comparison, variant_id="ref-2"))
    assert evaluate_bp2(replace(f, observations=(a, b, c))).status == DecisionStatus.ACTIVATED


def test_duplicate_case_reports_not_counted():
    f = fixture(gene="DICER1", disease="MONDO:0100216", inheritance=Inheritance.AD)
    a = f.observations[0]
    duplicate = replace(a, study_id="second publication")
    third = replace(a, case_id="case-2", comparison=replace(a.comparison, variant_id="ref-2"))
    assert evaluate_bp2(replace(f, observations=(a, duplicate, third))).status == DecisionStatus.UNDETERMINED


def test_dicer1_requires_vcep_comparison_classification():
    f = fixture(Phase.TRANS, gene="DICER1", disease="MONDO:0100216", inheritance=Inheritance.AD)
    assert evaluate_bp2(change_comp(f, vcep_classified=N)).status == DecisionStatus.NOT_MET


def test_bmpr2_cis_only_and_affected():
    f = fixture(gene="BMPR2", disease="MONDO:0015924", inheritance=Inheritance.AD,
                classification=Classification.LIKELY_PATHOGENIC)
    assert evaluate_bp2(f).status == DecisionStatus.ACTIVATED
    assert evaluate_bp2(change_obs(f, affected=N)).status == DecisionStatus.NOT_MET
    assert evaluate_bp2(change_obs(f, diagnosis_matches=N)).status == DecisionStatus.NOT_MET
    trans = fixture(Phase.TRANS, gene="BMPR2", disease="MONDO:0015924", inheritance=Inheritance.AD)
    assert evaluate_bp2(trans).status == DecisionStatus.NOT_MET


def test_version_mismatch():
    f = replace(fixture(gene="BMPR2", disease="MONDO:0015924", inheritance=Inheritance.AD),
                guideline_version="1.0")
    assert evaluate_bp2(f).status == DecisionStatus.UNDETERMINED


def test_vcep_gene_with_wrong_disease_cannot_fall_back_to_base():
    f = fixture(gene="DICER1", disease="MONDO:OTHER", inheritance=Inheritance.AD)
    assert evaluate_bp2(f).status == DecisionStatus.UNDETERMINED


def test_pm3_overlap_recorded_without_bp2_in_recessive_trans():
    f = replace(fixture(Phase.TRANS), other_criteria_from_same_observation=("PM3",))
    r = evaluate_bp2(f)
    assert r.status == DecisionStatus.NOT_MET
    assert r.criteria_overlap == ("PM3",)


def test_no_arbitrary_strength_upgrade():
    f = fixture()
    more = replace(f, observations=(f.observations[0], replace(f.observations[0], case_id="case-2")))
    assert evaluate_bp2(more).strength == "SUPPORTING"


def test_conflicting_phase_flag():
    assert evaluate_bp2(change_obs(fixture(), phase_conflict=Y)).status == DecisionStatus.UNDETERMINED


def test_invalid_gene_disease_relationship_is_not_applicable():
    assert evaluate_bp2(replace(fixture(), gene_disease_valid=N)).status == DecisionStatus.NOT_APPLICABLE


def test_contradictory_duplicate_phase_reports_are_unresolved():
    f = fixture()
    a = f.observations[0]
    contradictory = replace(a, phase=Phase.TRANS)
    assert evaluate_bp2(replace(f, observations=(a, contradictory))).status == DecisionStatus.UNDETERMINED


def test_two_comparisons_in_one_case_still_count_one_case():
    f = fixture(gene="DICER1", disease="MONDO:0100216", inheritance=Inheritance.AD)
    a = f.observations[0]
    second_comparison = replace(a, comparison=replace(a.comparison, variant_id="ref-2"))
    b = replace(a, case_id="case-2")
    c = replace(a, case_id="case-3")
    assert evaluate_bp2(replace(f, observations=(a, second_comparison, b, c))).status == DecisionStatus.ACTIVATED
