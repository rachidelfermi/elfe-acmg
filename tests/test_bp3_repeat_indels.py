from dataclasses import replace
import pytest
from adaptive_bp3 import BiologicalFeatures, evaluate_bp3


def base(**changes):
    data = dict(gene_id="GENE1", disease_id="DISEASE1", variant_consequence="INFRAME_DELETION",
                reading_frame_preserved=True, biologically_relevant_transcript=True,
                protein_consequence_confidence=True, repeat_region_status=True,
                repeat_region_functional_importance="NO_KNOWN_FUNCTION",
                repeat_length_disease_mechanism=False)
    data.update(changes)
    return BiologicalFeatures(**data)


def test_minimum_deletion():
    r = evaluate_bp3(base())
    assert (r.status, r.strength, r.decision_path) == ("ACTIVATED", "SUPPORTING", "STANDARD_NONFUNCTIONAL_REPEAT")
    assert not r.missing_critical_features


def test_insertion():
    assert evaluate_bp3(base(variant_consequence="INFRAME_INSERTION")).status == "ACTIVATED"


def test_outside_repeat():
    assert evaluate_bp3(base(repeat_region_status=False)).status == "NOT_MET"


def test_frameshift():
    assert evaluate_bp3(base(variant_consequence="FRAMESHIFT", reading_frame_preserved=False)).status == "NOT_MET"


def test_benign_polymorphism_needs_function_assessment():
    assert evaluate_bp3(base(known_benign_repeat_length_variants=True,
                             repeat_region_functional_importance=None)).status == "UNDETERMINED"


def test_benign_polymorphism_with_complete_pathway():
    assert evaluate_bp3(base(known_benign_repeat_length_variants=True)).status == "ACTIVATED"


def test_no_invented_activation():
    assert evaluate_bp3(BiologicalFeatures()).status == "UNDETERMINED"


def test_irrelevant_fields_do_not_change_decision():
    a = evaluate_bp3(base())
    b = evaluate_bp3(replace(base(), coding_hgvs="c.10_12del", repeat_unit_length=3))
    assert (a.status, a.strength, a.decision_path) == (b.status, b.strength, b.decision_path)


def test_determinism():
    assert evaluate_bp3(base()).to_dict() == evaluate_bp3(base()).to_dict()


@pytest.mark.parametrize("change,status", [(0, "NOT_MET"), (1, "UNDETERMINED")])
def test_length_change(change, status):
    assert evaluate_bp3(base(protein_length_change=change)).status == status
