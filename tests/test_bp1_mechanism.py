from adaptive_bp1 import evaluate_bp1
from test_bp1_eligibility import base


def test_missense_and_mixed_mechanism():
    assert evaluate_bp1(base(missense_common_mechanism=True)).status == "NOT_MET"
    assert evaluate_bp1(base(truncating_predominant=False, missense_common_mechanism=True)).status == "NOT_MET"


def test_same_gene_different_disease():
    first = evaluate_bp1(base())
    second = evaluate_bp1(base(disease_id="SYNTHETIC:2", truncating_predominant=False))
    assert first.status == "ACTIVATED" and second.status == "NOT_MET"


def test_curated_spectrum_without_numeric_threshold():
    data = base(truncating_predominant=None, spectrum_truncating_predominant=True,
                spectrum_evidence_quality="CURATED", spectrum_classification_independent=True,
                spectrum_unique_variants=True,
                spectrum_disease_match=True, spectrum_ascertainment_bias_addressed=True)
    result = evaluate_bp1(data)
    assert (result.status, result.decision_path) == ("ACTIVATED", "CURATED_PATHOGENIC_SPECTRUM")


def test_counts_alone_never_activate():
    data = base(truncating_predominant=None, pathogenic_truncating_variant_count=100,
                pathogenic_missense_variant_count=0)
    assert evaluate_bp1(data).status == "UNDETERMINED"


def test_circular_or_duplicated_spectrum_unresolved():
    data = base(truncating_predominant=None, spectrum_truncating_predominant=True,
                spectrum_evidence_quality="CURATED", spectrum_classification_independent=False,
                spectrum_unique_variants=True,
                spectrum_disease_match=True, spectrum_ascertainment_bias_addressed=True)
    assert evaluate_bp1(data).status == "UNDETERMINED"
    data["spectrum_classification_independent"] = None
    assert evaluate_bp1(data).status == "UNDETERMINED"


def test_duplicate_submissions_do_not_activate():
    data = base(truncating_predominant=None, spectrum_truncating_predominant=True,
                spectrum_evidence_quality="CURATED", spectrum_classification_independent=True,
                spectrum_unique_variants=False, spectrum_disease_match=True,
                spectrum_ascertainment_bias_addressed=True)
    assert evaluate_bp1(data).status == "UNDETERMINED"


def test_uncurated_assessment_does_not_activate():
    assert evaluate_bp1(base(mechanism_evidence_quality="LOW")).status == "UNDETERMINED"
    assert evaluate_bp1(base(mechanism_evidence_source=None)).status == "UNDETERMINED"
