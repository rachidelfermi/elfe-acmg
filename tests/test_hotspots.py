from dataclasses import replace

from adaptive_pm1 import BiologicalFeatures, Evidence, Region, evaluate_pm1


def hotspot(**changes):
    base = Region(
        "curated-hotspot", "HOTSPOT", 10, 20, "NM_TEST.1",
        definition_source="curated published regional review",
        established=Evidence.known(True), disease_relevant=Evidence.known(True),
        mechanism_compatible=Evidence.known(True),
        pathogenic_clustering=Evidence.known(True),
        benign_variation_assessed=Evidence.known(True),
        lacks_benign_variation=Evidence.known(True),
        classification_independent=Evidence.known(True),
        evidence_quality=Evidence.known(True),
    )
    return replace(base, **changes)


def case(*regions, **changes):
    base = BiologicalFeatures(gene_id="TEST", disease_id="TEST-DISEASE",
                              variant_consequence="missense",
                              transcript_id_version="NM_TEST.1", protein_position=15,
                              regions=regions)
    return replace(base, **changes)


def test_validated_hotspot_minimum_evidence():
    result = evaluate_pm1(case(hotspot()))
    assert (result.status, result.strength, result.decision_path) == (
        "ACTIVATED", "MODERATE", "VALIDATED_MUTATIONAL_HOTSPOT")


def test_outside_hotspot():
    assert evaluate_pm1(case(hotspot(), protein_position=21)).status == "NOT_MET"


def test_insufficient_pathogenic_evidence():
    assert evaluate_pm1(case(hotspot(pathogenic_clustering=Evidence.known(False)))).status == "NOT_MET"


def test_benign_variation_excludes():
    assert evaluate_pm1(case(hotspot(lacks_benign_variation=Evidence.known(False)))).status == "NOT_MET"


def test_circularity_blocks():
    assert evaluate_pm1(case(hotspot(classification_independent=Evidence.known(False)))).status == "UNDETERMINED"


def test_conflicting_reference_classifications_block():
    assert evaluate_pm1(case(hotspot(classification_conflicts=Evidence.known(True)))).status == "UNDETERMINED"


def test_missing_benign_assessment_blocks():
    assert evaluate_pm1(case(hotspot(benign_variation_assessed=Evidence()))).status == "UNDETERMINED"


def test_mechanism_mismatch():
    assert evaluate_pm1(case(hotspot(mechanism_compatible=Evidence.known(False)))).status == "NOT_MET"


def test_isoform_mismatch():
    assert evaluate_pm1(case(hotspot(), transcript_id_version="NM_TEST.2")).status == "NOT_MET"


def test_poor_population_coverage_blocks_claim_of_no_benign_variation():
    assert evaluate_pm1(case(hotspot(population_coverage_adequate=Evidence.known(False)))).status == "UNDETERMINED"


def test_irrelevant_audit_features_do_not_change_strength():
    first = evaluate_pm1(case(hotspot()))
    second = evaluate_pm1(case(hotspot(pathogenic_count=40, enrichment_statistic=20.0)))
    assert (first.status, first.strength) == (second.status, second.strength)
