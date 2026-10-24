import pytest

from adaptive_pp3 import evaluate_pp3
from .pp3_fixtures import splicing


@pytest.mark.parametrize("score,status", [(0.199, "NOT_MET"), (0.2, "ACTIVATED"),
                                          (0.8, "ACTIVATED"), (0.1, "NOT_MET")])
def test_spliceai_boundary(score, status):
    assert evaluate_pp3(splicing(score)).status.value == status


def test_canonical_site_excluded():
    assert evaluate_pp3(splicing(0.9, canonical_splice_site=True)).status.value == "UNDETERMINED"


def test_pvs1_and_rna_overlap():
    for field in ("pvs1_splice_overlap", "rna_splice_evidence_used"):
        assert evaluate_pp3(splicing(0.8, **{field: True})).status.value == "UNDETERMINED"


def test_low_splice_score_is_only_bp4_candidate():
    result = evaluate_pp3(splicing(0.1)).to_dict()
    assert result["status"] == "NOT_MET" and result["bp4_candidate"]
