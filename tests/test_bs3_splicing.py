from adaptive_bs3 import evaluate_bs3
from tests.test_bs3_functional_assays import case


def test_direct_normal_rna_is_routed_to_bp7_not_bs3():
    data = case(assay_type="RNA", tissue_relevant=True, transcript_relevant=True,
                rna_detection_adequate=True, nmd_addressed=True,
                variant_allele_expressed=True)
    result = evaluate_bs3(data)
    assert result.status == "NOT_APPLICABLE"
    assert result.guideline["rule_id"] == "BP7-RNA-not-BS3"


def test_rna_with_unknown_tissue_is_unresolved():
    assert evaluate_bs3(case(assay_type="RNA")).status == "UNDETERMINED"
