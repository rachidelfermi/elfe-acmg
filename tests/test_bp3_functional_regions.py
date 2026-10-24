import pytest
from adaptive_bp3 import evaluate_bp3
from test_bp3_repeat_indels import base


@pytest.mark.parametrize("name", ["critical_protein_domain_overlap", "active_site_overlap",
                                  "binding_site_overlap", "protein_interaction_region_overlap"])
def test_critical_region(name):
    assert evaluate_bp3(base(**{name: True})).status == "NOT_MET"


@pytest.mark.parametrize("name", ["known_pathogenic_repeat_length_variants", "repeat_length_disease_mechanism",
                                  "gene_specific_repeat_pathogenicity", "disease_specific_repeat_exceptions"])
def test_pathogenic_repeat_evidence(name):
    assert evaluate_bp3(base(**{name: True})).status == "NOT_MET"


def test_established_function():
    assert evaluate_bp3(base(repeat_region_functional_importance="ESTABLISHED_FUNCTION")).status == "NOT_MET"


def test_damaging_validation():
    assert evaluate_bp3(base(repeat_region_functional_validation="DAMAGING")).status == "NOT_MET"


def test_low_complexity_only():
    r = evaluate_bp3(base(repeat_region_status=None, repeat_region_functional_importance=None,
                          low_complexity_region_status=True))
    assert r.status == "UNDETERMINED" and "repeat_region_status" in r.missing_critical_features
