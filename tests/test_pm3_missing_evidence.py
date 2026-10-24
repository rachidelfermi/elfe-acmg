from adaptive_pm3 import evaluate_pm3
from tests.pm3_fixtures import case, data, second


def test_global_scope_and_rarity():
    assert evaluate_pm3(data(case(), inheritance_mode="AUTOSOMAL_DOMINANT"))["status"] == "NOT_APPLICABLE"
    assert evaluate_pm3(data(case(), candidate_rare_for_disease=None))["status"] == "UNDETERMINED"
    assert evaluate_pm3(data(case(second_variant=second(rare_for_disease=False))))["status"] == "NOT_MET"


def test_technical_and_disease_mismatch():
    assert evaluate_pm3(data(case(candidate_call_valid=False)))["status"] == "NOT_MET"
    assert evaluate_pm3(data(case(second_variant=second(disease_match=False))))["status"] == "NOT_MET"
    assert evaluate_pm3(data(case(affected=None)))["status"] == "UNDETERMINED"
    assert evaluate_pm3(data(case(affected=False)))["status"] == "NOT_MET"


def test_optional_fields_absent_and_dict_input():
    result = evaluate_pm3(data(case()))
    assert result["status"] == "ACTIVATED"
    raw = {"candidate_variant_id": "v1", "gene_id": "PAH", "disease_id": "PKU",
           "inheritance_mode": "AUTOSOMAL_RECESSIVE", "candidate_rare_for_disease": True,
           "cases": [{"case_id": "c1", "affected": True, "clinical_disease_match": True,
                      "candidate_genotype": "HET", "candidate_call_valid": True,
                      "second_variant": vars(second()), "phase": "TRANS", "phase_method": "molecular"}]}
    assert evaluate_pm3(raw)["strength"] == "MODERATE"
