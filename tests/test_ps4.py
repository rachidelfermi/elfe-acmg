"""PS4 decision boundaries and evidence-integrity checks."""

from adaptive_ps4 import evaluate_ps4
from adaptive_ps4.statistics import fisher_or_ci


BASE = {"variant_id": "v1", "gene_id": "GENE", "disease_id": "D",
        "gene_disease_validity": True}
VALID = {"study_id": "s1", "variant_id": "v1", "disease_id": "D",
         "cases_affected": True, "controls_valid": True, "variant_calls_valid": True,
         "technical_comparable": True, "independent_cases": True,
         "overlapping_cohort": False, "ancestry_comparable": True}


def decision(study, **top):
    return evaluate_ps4({**BASE, "studies": [{**VALID, **study}], **top})


def test_exact_interval_and_significant_enrichment():
    r = decision({"case_carriers": 30, "case_noncarriers": 70,
                  "control_carriers": 1, "control_noncarriers": 999})
    assert r.status == "ACTIVATED" and r.strength == "STRONG"
    assert r.odds_ratio > 5 and r.confidence_interval[0] > 1
    assert r.p_value < .05


def test_no_enrichment_and_protective():
    for a, b, c, d in [(10, 90, 100, 900), (1, 99, 10, 90), (0, 100, 5, 100)]:
        assert decision({"case_carriers": a, "case_noncarriers": b,
                         "control_carriers": c, "control_noncarriers": d}).status == "NOT_MET"


def test_small_sample_and_zero_control_cell():
    r = decision({"case_carriers": 1, "case_noncarriers": 9,
                  "control_carriers": 0, "control_noncarriers": 10})
    assert r.status == "NOT_MET" and r.confidence_interval[0] <= 1
    assert fisher_or_ci(5, 95, 0, 100)[0] == float("inf")


def test_missing_counts_never_become_zero():
    r = decision({"case_carriers": 5, "case_noncarriers": 95})
    assert r.status == "UNDETERMINED"
    assert "complete_case_control_2x2_counts" in r.missing_critical_features


def test_published_estimate_with_provenance():
    r = decision({"odds_ratio": 8, "confidence_interval": [2, 30], "confidence_level": .95,
                  "statistic_source": "PMID:example", "statistic_method": "logistic_regression"})
    assert r.status == "ACTIVATED" and r.decision_path == "PUBLISHED_ASSOCIATION"
    assert decision({"odds_ratio": 8}).status == "UNDETERMINED"
    assert decision({"relative_risk": 8, "confidence_interval": [2, 30],
                     "confidence_level": .95, "statistic_source": "PMID:example",
                     "statistic_method": "cohort_model"}).status == "ACTIVATED"
    assert decision({"odds_ratio": 8, "confidence_interval": [2, 30],
                     "statistic_source": "PMID:example", "statistic_method": "logistic_regression"}).status == "UNDETERMINED"


def test_inconsistent_statistics():
    r = decision({"case_carriers": 30, "case_noncarriers": 70,
                  "control_carriers": 1, "control_noncarriers": 999, "odds_ratio": 1.2})
    assert r.status == "UNDETERMINED"
    assert "reported_odds_ratio_vs_raw_counts" in r.conflicting_evidence
    zero = decision({"case_carriers": 5, "case_noncarriers": 95,
                     "control_carriers": 0, "control_noncarriers": 100, "odds_ratio": 3})
    assert zero.status == "UNDETERMINED"


def test_ancestry_founder_batch_and_duplicates_block():
    counts = {"case_carriers": 30, "case_noncarriers": 70,
              "control_carriers": 1, "control_noncarriers": 999}
    for extra in ({"ancestry_comparable": False}, {"founder_effect": True},
                  {"batch_effect": True}, {"overlapping_cohort": True},
                  {"independent_cases": False}):
        assert decision({**counts, **extra}).status == "UNDETERMINED"
    assert decision({**counts, "ancestry_comparable": False,
                     "stratification_adjusted": True}).status == "ACTIVATED"


def test_myoc_scoped_proband_thresholds():
    top = {"variant_id": "v1", "gene_id": "MYOC", "disease_id": "MONDO:0005338",
           "gene_disease_validity": True, "pm2_met": True, "studies_independent": True}
    def cases(n):
        return [{"case_id": f"c{i}", "study_id": f"s{i % 2}", "family_id": f"f{i}",
                 "clinically_assessed": True, "phenotype": "JOAG",
                 "has_other_myoc_vus_lp_p": False} for i in range(n)]
    for n, strength in [(2, "SUPPORTING"), (6, "MODERATE"), (15, "STRONG")]:
        r = evaluate_ps4({**top, "probands": cases(n)})
        assert r.status == "ACTIVATED" and r.strength == strength
    assert evaluate_ps4({**top, "probands": cases(2), "pm2_met": False}).status == "NOT_MET"
    assert evaluate_ps4({**top, "probands": cases(2), "studies_independent": None}).status == "UNDETERMINED"
    assert evaluate_ps4({**top, "probands": cases(2), "gene_id": "OTHER",
                         "specification": "MYOC_2.1"}).status == "NOT_APPLICABLE"
    duplicate = cases(2)
    duplicate[1]["family_id"] = duplicate[0]["family_id"]
    assert evaluate_ps4({**top, "probands": duplicate}).status == "UNDETERMINED"


def test_opposing_studies_require_review():
    a = {**VALID, "case_carriers": 30, "case_noncarriers": 70,
         "control_carriers": 1, "control_noncarriers": 999}
    b = {**VALID, "study_id": "s2", "case_carriers": 1, "case_noncarriers": 99,
         "control_carriers": 10, "control_noncarriers": 90}
    r = evaluate_ps4({**BASE, "studies": [a, b]})
    assert r.status == "UNDETERMINED" and r.conflicting_evidence
