from adaptive_bp3 import BiologicalFeatures, evaluate_bp3
from test_bp3_repeat_indels import base


def foxg1(**changes):
    data = dict(gene_id="FOXG1", disease_id="MONDO:0100040", variant_consequence="INFRAME_DELETION",
                reading_frame_preserved=True, biologically_relevant_transcript=True,
                protein_consequence_confidence=True, transcript_id_version="NM_005249.4",
                protein_accession_version="NP_005240.3", affected_protein_positions=(50, 51))
    data.update(changes)
    return BiologicalFeatures(**data)


def test_foxg1_curated_repeat():
    r = evaluate_bp3(foxg1())
    assert (r.status, r.strength, r.guideline["version"]) == ("ACTIVATED", "SUPPORTING", "2")


def test_foxg1_outside_region():
    assert evaluate_bp3(foxg1(affected_protein_positions=(90, 91))).status == "NOT_MET"


def test_foxg1_protein_mismatch():
    assert evaluate_bp3(foxg1(protein_accession_version="NP_005240.2")).status == "UNDETERMINED"


def test_foxg1_transcript_mismatch():
    assert evaluate_bp3(foxg1(transcript_id_version="NM_005249.3")).status == "UNDETERMINED"


def test_foxg1_exclusion():
    assert evaluate_bp3(foxg1(known_pathogenic_repeat_length_variants=True)).status == "NOT_MET"


def test_foxg1_strength_stays_supporting():
    assert evaluate_bp3(foxg1(known_benign_repeat_length_variants=True)).strength == "SUPPORTING"


def test_version_mismatch():
    assert evaluate_bp3(foxg1(vcep_specification_version="1")).status == "UNDETERMINED"


def test_ryr1_not_applicable():
    assert evaluate_bp3(base(gene_id="RYR1", disease_id="MONDO:0018493")).status == "NOT_APPLICABLE"


def test_myoc_not_applicable():
    assert evaluate_bp3(base(gene_id="MYOC", disease_id="PRIMARY OPEN-ANGLE GLAUCOMA")).status == "NOT_APPLICABLE"


def test_bmpr2_supporting():
    r = evaluate_bp3(base(gene_id="BMPR2", disease_id="MONDO:0015924",
                          transcript_id_version="NM_001204.7"))
    assert (r.status, r.strength) == ("ACTIVATED", "SUPPORTING")


def test_scope_unknown():
    assert evaluate_bp3(base(gene_id="FOXG1", disease_id="unrelated disease")).status == "UNDETERMINED"
