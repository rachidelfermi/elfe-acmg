"""Synthetic fixtures: rule behavior, not externally validated variant classifications."""
from dataclasses import replace
import pytest

from adaptive_bp7 import BP7Features, RNAEvidence, SplicePrediction, EvidenceState, evaluate_bp7


def base(**changes):
    f = BP7Features(
        variant_id="synthetic", gene_id="GENE", disease_id="disease", genome_assembly="GRCh38",
        transcript_id_version="NM_TEST.1", variant_consequence="synonymous", consequence_confident=True,
        biologically_relevant_transcript=True, variant_mapping_valid=True, exon_core_splice_motif=False,
        prediction=SplicePrediction(score=0.05, transcript_id_version="NM_TEST.1", genome_assembly="GRCh38",
                                    variant_match=True, applicable=True, calibrated=True, quality=EvidenceState.KNOWN),
    )
    return replace(f, **changes)


def original(**changes):
    f = base(guideline="ACMG_2015", conservation="not_highly_conserved", conservation_reliable=True,
             prediction=replace(base().prediction, no_splice_effect=True))
    return replace(f, **changes)


def tp53(**changes):
    f = base(guideline="TP53_2.1.0", gene_id="TP53", disease_id="MONDO:0018875")
    return replace(f, **changes)


def rna(**changes):
    x = RNAEvidence(normal_splicing=True, abnormal_splicing=False, validated=True, controls=True,
                    relevant_tissue=True, adequate_sensitivity=True, transcript_id_version="NM_TEST.1",
                    variant_match=True, source="synthetic assay")
    return replace(x, **changes)


@pytest.mark.parametrize("case,status", [
    (base(), "ACTIVATED"),
    (base(prediction=replace(base().prediction, score=0.2)), "NOT_MET"),
    (base(prediction=replace(base().prediction, score=0.1)), "ACTIVATED"),
    (base(prediction=None), "UNDETERMINED"),
    (base(prediction=replace(base().prediction, transcript_id_version="NM_OTHER.1")), "UNDETERMINED"),
    (base(prediction=replace(base().prediction, genome_assembly="GRCh37")), "UNDETERMINED"),
    (base(exon_core_splice_motif=True), "NOT_MET"),
    (base(regulatory_overlap=True), "UNDETERMINED"),
    (base(regulatory_overlap=True, established_regulatory_splice_effect=False), "ACTIVATED"),
    (base(rna=rna()), "ACTIVATED"),
    (base(rna=rna(normal_splicing=False, abnormal_splicing=True)), "NOT_MET"),
    (base(rna=rna(), prediction=replace(base().prediction, score=0.2)), "UNDETERMINED"),
    (base(known_pathogenic_synonymous_mechanism=True), "NOT_MET"),
    (base(alternative_transcript_effect=True), "NOT_MET"),
    (base(variant_consequence="missense"), "NOT_APPLICABLE"),
    (base(consequence_confident=None), "UNDETERMINED"),
    (base(prediction=replace(base().prediction, calibrated=False)), "UNDETERMINED"),
    (base(prediction=replace(base().prediction, conflict=True)), "UNDETERMINED"),
    (base(rna=rna(already_used_for_bs3=True)), "ACTIVATED"),
    (original(), "ACTIVATED"),
    (original(conservation="highly_conserved"), "NOT_MET"),
    (original(conservation=None), "UNDETERMINED"),
    (original(prediction=replace(original().prediction, no_splice_effect=False)), "NOT_MET"),
    (tp53(), "ACTIVATED"),
    (tp53(gene_id="OTHER"), "NOT_APPLICABLE"),
    (tp53(variant_consequence="intronic", intronic_position=7), "ACTIVATED"),
    (tp53(variant_consequence="intronic", intronic_position=6), "NOT_MET"),
    (tp53(variant_consequence="intronic", intronic_position=-21), "ACTIVATED"),
    (tp53(variant_consequence="intronic", intronic_position=-20), "NOT_MET"),
    (tp53(rna=rna()), "ACTIVATED"),
    (base(rna=rna(), prediction=None), "ACTIVATED"),
    (tp53(rna=rna(), prediction=None), "UNDETERMINED"),
    (base(rna=rna(already_used_for_bs3=True), prediction=None), "UNDETERMINED"),
])
def test_synthetic_rules(case, status):
    assert evaluate_bp7(case).status == status


def test_strength_and_overlap():
    result = evaluate_bp7(base())
    assert result.strength == "SUPPORTING"
    assert result.overlap_note and "BP4" in result.overlap_note
    result = evaluate_bp7(base(rna=rna()))
    assert result.strength == "STRONG"
    assert "BS3" in result.overlap_note


def test_missing_optional_features_do_not_change_decision():
    a = evaluate_bp7(base())
    b = evaluate_bp7(base(variant_id="another", conservation="highly_conserved", conservation_reliable=True,
                          bp4_splicing_applied=True))
    assert (a.status, a.strength, a.decision_path) == (b.status, b.strength, b.decision_path)


def test_dict_api_and_no_invented_activation():
    data = {"variant_consequence": "synonymous", "consequence_confident": True,
            "biologically_relevant_transcript": True, "variant_mapping_valid": True,
            "transcript_id_version": "NM_TEST.1", "exon_core_splice_motif": False}
    result = evaluate_bp7(data)
    assert result.status == "UNDETERMINED" and "prediction" in result.missing_critical_features


def test_invalid_score_rejected():
    with pytest.raises(ValueError):
        evaluate_bp7(base(prediction=replace(base().prediction, score=1.5)))
