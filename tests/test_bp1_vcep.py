from adaptive_bp1 import evaluate_bp1
from test_bp1_eligibility import base


def test_dicer1_vcep_prohibits_bp1():
    r = evaluate_bp1(base(gene_id="DICER1", disease_id="MONDO:0100216"))
    assert (r.status, r.strength, r.guideline["rule_id"]) == ("NOT_APPLICABLE", "NONE", "GN024:BP1")


def test_ryr1_scope_is_disease_specific():
    assert evaluate_bp1(base(gene_id="RYR1", disease_id="MONDO:0018493")).status == "NOT_APPLICABLE"
    assert evaluate_bp1(base(gene_id="RYR1", disease_id="SYNTHETIC:OTHER")).status == "ACTIVATED"


def test_map2k1_vcep_reverses_variant_class():
    data = base(gene_id="MAP2K1", disease_id="MONDO:0021060", variant_consequence="FRAMESHIFT",
                loss_of_function_disease_mechanism=False, missense_gain_of_function_mechanism=True)
    r = evaluate_bp1(data)
    assert (r.status, r.strength, r.decision_path) == ("ACTIVATED", "SUPPORTING", "VCEP_RASOPATHY_LOF")
    assert evaluate_bp1({**data, "variant_consequence": "MISSENSE"}).status == "NOT_MET"
    assert evaluate_bp1({**data, "loss_of_function_disease_mechanism": True}).status == "NOT_MET"


def test_unverified_vcep_version_cannot_adjust_strength():
    data = base(gene_id="MAP2K1", disease_id="MONDO:0021060", variant_consequence="FRAMESHIFT",
                loss_of_function_disease_mechanism=False, missense_gain_of_function_mechanism=True,
                vcep_specification_version="999")
    assert evaluate_bp1(data).status == "UNDETERMINED"
