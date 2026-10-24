from adaptive_pm1 import evaluate_pm1


def vcep(gene, disease, transcript, position, **extra):
    return evaluate_pm1(dict(gene_id=gene, disease_id=disease,
                             transcript_id_version=transcript, protein_position=position,
                             variant_consequence="missense", **extra))


def test_tp53_curated_codon():
    d = vcep("TP53", "MONDO:0018875", "NM_000546.5", 175)
    assert (d.status, d.strength) == ("ACTIVATED", "MODERATE")


def test_tp53_somatic_supporting_and_moderate():
    common = dict(somatic_source_version="cancerhotspots.org:v2", amino_acid_change="p.Ala100Val")
    assert vcep("TP53", "MONDO:0018875", "NM_000546.5", 100,
                somatic_same_amino_acid_observations=2, **common).strength == "SUPPORTING"
    assert vcep("TP53", "MONDO:0018875", "NM_000546.5", 100,
                somatic_same_amino_acid_observations=10, **common).strength == "MODERATE"


def test_tp53_wrong_somatic_source_is_unresolved():
    assert vcep("TP53", "MONDO:0018875", "NM_000546.5", 100,
                somatic_same_amino_acid_observations=10).status == "UNDETERMINED"


def test_dicer1_metal_codon_and_domain_strength():
    args = ("DICER1", "MONDO:0100216", "NM_177438.2")
    assert vcep(*args, 1705).strength == "MODERATE"
    assert vcep(*args, 1720).strength == "SUPPORTING"


def test_dicer1_pm5_exclusion():
    assert vcep("DICER1", "MONDO:0100216", "NM_177438.2", 1705,
                pm5_strength="MODERATE").status == "NOT_MET"


def test_f9_residue_and_exon_strength():
    args = ("F9", "MONDO:0010604", "NM_000133.4")
    assert vcep(*args, 267).strength == "STRONG"
    assert vcep(*args, 100, exon_number=4).strength == "MODERATE"


def test_f9_combined_pm5_requires_coordination():
    assert vcep("F9", "MONDO:0010604", "NM_000133.4", 267,
                pm5_strength="MODERATE").status == "UNDETERMINED"


def test_vcep_transcript_mismatch():
    assert vcep("TP53", "MONDO:0018875", "NM_000546.4", 175).status == "NOT_MET"


def test_vcep_disease_scope():
    assert vcep("TP53", "OTHER", "NM_000546.5", 175).status == "UNDETERMINED"
