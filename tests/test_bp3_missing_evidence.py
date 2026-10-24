import pytest
from adaptive_bp3 import Evidence, EvidenceState, evaluate_bp3
from test_bp3_repeat_indels import base


@pytest.mark.parametrize("name", ["repeat_region_status", "repeat_region_functional_importance",
                                  "variant_consequence",
                                  "biologically_relevant_transcript", "protein_consequence_confidence",
                                  "repeat_length_disease_mechanism"])
def test_missing_essential(name):
    assert evaluate_bp3(base(**{name: None})).status == "UNDETERMINED"


def test_missing_optional():
    assert evaluate_bp3(base(reading_frame_preserved=None, indel_alignment_reliability=None,
                             rna_splicing_evidence=None)).status == "ACTIVATED"


def test_conflicting_state():
    assert evaluate_bp3(base(repeat_region_status=Evidence(EvidenceState.CONFLICTING))).status == "UNDETERMINED"


def test_unreliable_alignment():
    assert evaluate_bp3(base(indel_alignment_reliability=False)).status == "UNDETERMINED"


@pytest.mark.parametrize("name", ["complex_allele_evidence", "alternative_splice_consequence",
                                  "alternative_transcript_consequence", "pm4_applied_to_same_consequence"])
def test_conflicting_consequence(name):
    assert evaluate_bp3(base(**{name: True})).status == "UNDETERMINED"


def test_coordinate_conflict():
    assert evaluate_bp3(base(affected_protein_positions=(50, 51),
                             repeat_region_boundaries=(100, 120))).status == "UNDETERMINED"


def test_transcript_mismatch():
    assert evaluate_bp3(base(biologically_relevant_transcript=False)).status == "NOT_MET"


def test_no_unverified_strength_upgrade():
    assert evaluate_bp3(base(known_benign_repeat_length_variants=True)).strength == "SUPPORTING"


def test_explicitly_unreliable_evidence_quality():
    assert evaluate_bp3(base(evidence_quality="UNRELIABLE")).status == "UNDETERMINED"


def test_irrelevant_conflicting_annotation_does_not_block():
    assert evaluate_bp3(base(repeat_length_variability=Evidence(EvidenceState.CONFLICTING))).status == "ACTIVATED"


def test_frame_label_conflict():
    assert evaluate_bp3(base(reading_frame_preserved=False)).status == "UNDETERMINED"
