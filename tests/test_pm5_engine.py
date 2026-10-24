"""All cases are synthetic, not externally validated variant classifications."""
from copy import deepcopy
import pytest
from adaptive_pm5 import evaluate_pm5


def case():
    return {
        "gene_id": "EXAMPLE", "disease_id": "MONDO:synthetic", "gene_disease_validity": "ESTABLISHED",
        "candidate": {"variant_id": "candidate", "consequence": "MISSENSE", "transcript_id_version": "NM_TEST.1", "protein_accession_version": "NP_TEST.1", "amino_acid_position": 132, "reference_amino_acid": "Arg", "alternate_amino_acid": "His", "transcript_relevant": True},
        "references": [{"variant_id": "reference-1", "consequence": "MISSENSE", "transcript_id_version": "NM_TEST.1", "protein_accession_version": "NP_TEST.1", "amino_acid_position": 132, "reference_amino_acid": "Arg", "alternate_amino_acid": "Cys", "classification": "PATHOGENIC", "classification_source": "curated synthetic dossier", "classification_evidence_reviewed": True, "classification_independent_of_pm5": True, "disease_match": True}],
    }


def result(data):
    return evaluate_pm5(data).to_dict()


def test_minimum_and_irrelevant_features():
    d = case()
    assert result(d)["strength"] == "MODERATE"
    d["candidate"]["splice_mechanism"] = "UNKNOWN"
    d["evidence_source"] = "extra provenance"
    assert result(d)["status"] == "ACTIVATED"


@pytest.mark.parametrize("change,status", [
    ({"alternate_amino_acid": "His"}, "NOT_MET"),
    ({"amino_acid_position": 133}, "NOT_MET"),
    ({"transcript_id_version": "NM_OTHER.1"}, "NOT_MET"),
    ({"protein_accession_version": "NP_OTHER.1"}, "NOT_MET"),
    ({"classification": "VUS"}, "NOT_MET"),
    ({"classification": "BENIGN"}, "NOT_MET"),
    ({"classification": "LIKELY_PATHOGENIC"}, "UNDETERMINED"),
    ({"classification": None}, "UNDETERMINED"),
    ({"classification_conflicting": True}, "UNDETERMINED"),
    ({"classification_independent_of_pm5": False}, "NOT_MET"),
    ({"disease_match": False}, "NOT_MET"),
    ({"splice_mechanism": "SPLICING"}, "NOT_MET"),
    ({"amino_acid_position": None}, "UNDETERMINED"),
])
def test_reference_screen(change, status):
    d = case(); d["references"][0].update(change)
    assert result(d)["status"] == status


@pytest.mark.parametrize("change,status", [
    ({"consequence": "SYNONYMOUS"}, "NOT_APPLICABLE"),
    ({"transcript_relevant": False}, "NOT_MET"),
    ({"amino_acid_position": None}, "UNDETERMINED"),
    ({"splice_mechanism": "SPLICING"}, "NOT_MET"),
])
def test_candidate_screen(change, status):
    d = case(); d["candidate"].update(change)
    assert result(d)["status"] == status


def test_missing_splice_prediction_and_exclusions():
    d = case()
    assert result(d)["status"] == "ACTIVATED"
    d["biological_evidence_conflicting"] = True
    assert result(d)["status"] == "UNDETERMINED"
    d.pop("biological_evidence_conflicting")
    d["known_benign_substitution_at_residue"] = True
    assert result(d)["status"] == "UNDETERMINED"


def test_multiple_refs_and_duplicate_reports():
    d = case()
    second = deepcopy(d["references"][0]); second.update(variant_id="reference-2", alternate_amino_acid="Ser")
    d["references"].append(second)
    assert result(d)["strength"] == "MODERATE"
    second["alternate_amino_acid"] = "Cys"
    assert result(d)["strength"] == "MODERATE"
    d["references"].append({"classification": None})
    assert result(d)["status"] == "ACTIVATED"


def test_platelet_vcep():
    d = case(); d.update(gene_id="ITGA2B", disease_id="MONDO:0010119")
    d["candidate"]["transcript_id_version"] = "NM_000419.4"
    d["references"][0].update(transcript_id_version="NM_000419.4", classification="LIKELY_PATHOGENIC")
    assert result(d)["strength"] == "SUPPORTING"
    d["references"][0]["classification"] = "PATHOGENIC"
    assert result(d)["strength"] == "MODERATE"
    d["disease_id"] = "MONDO:other"
    assert result(d)["status"] == "UNDETERMINED"


def test_hnf1a_vcep():
    d = case(); d.update(gene_id="HNF1A", disease_id="MONDO:0015967")
    d["candidate"].update(transcript_id_version="NM_000545.8", grantham_distance=100)
    d["references"][0].update(transcript_id_version="NM_000545.8", grantham_distance=80)
    assert result(d)["strength"] == "MODERATE"
    d["candidate"]["grantham_distance"] = 60
    assert result(d)["strength"] == "SUPPORTING"
    second = deepcopy(d["references"][0]); second.update(variant_id="reference-2", alternate_amino_acid="Ser")
    d["references"].append(second)
    assert result(d)["strength"] == "STRONG"
    second["alternate_amino_acid"] = "Cys"
    assert result(d)["strength"] == "SUPPORTING"


def test_missing_essentials_and_unknown_rule():
    d = case(); d["references"] = []
    assert result(d)["status"] == "UNDETERMINED"
    d = case(); d.pop("gene_disease_validity")
    assert result(d)["status"] == "UNDETERMINED"
    d = case(); d["specification"] = "UNVERIFIED"
    assert result(d)["status"] == "UNDETERMINED"


@pytest.mark.parametrize("state,status", [("UNKNOWN", "UNDETERMINED"), ("CONFLICTING", "UNDETERMINED"), ("NOT_APPLICABLE", "NOT_APPLICABLE")])
def test_explicit_core_states(state, status):
    d = case(); d["evidence_states"] = {"candidate": state}
    assert result(d)["status"] == status
