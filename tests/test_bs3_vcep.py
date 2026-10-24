from adaptive_bs3 import evaluate_bs3
from tests.test_bs3_functional_assays import case


def test_myoc_vcep_caps_strength_and_limits_assay_type():
    data = case(assay_type="SECRETION")
    data["context"]["gene_id"] = "MYOC"
    data["context"]["disease_id"] = "PRIMARY_OPEN_ANGLE_GLAUCOMA"
    assert evaluate_bs3(data).strength == "MODERATE"
    data["assays"][0]["assay_type"] = "ENZYME"
    assert evaluate_bs3(data).status == "NOT_APPLICABLE"
