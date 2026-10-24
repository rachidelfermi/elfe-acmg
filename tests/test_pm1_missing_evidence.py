from adaptive_pm1 import Evidence, State, evaluate_pm1
from .test_hotspots import case, hotspot


def test_no_regions_is_unknown():
    assert evaluate_pm1(case()).status == "UNDETERMINED"


def test_missing_position():
    assert evaluate_pm1(case(hotspot(), protein_position=None)).status == "UNDETERMINED"


def test_conflicting_regional_evidence():
    r = hotspot(established=Evidence(State.CONFLICTING))
    assert evaluate_pm1(case(r)).status == "UNDETERMINED"


def test_conflicting_duplicate_region_annotations():
    assert evaluate_pm1(case(hotspot(), hotspot(established=Evidence.known(False)))).status == "UNDETERMINED"


def test_missing_protein_accession_for_accession_scoped_region():
    assert evaluate_pm1(case(hotspot(protein_accession_version="NP_TEST.1"))).status == "UNDETERMINED"


def test_missing_optional_data_does_not_block():
    assert evaluate_pm1(case(hotspot(population_coverage_adequate=Evidence()))).status == "ACTIVATED"


def test_input_envelopes():
    d = case(hotspot())
    payload = {**d.__dict__, "regions": [hotspot()]}
    payload["gene_disease_validity"] = {"state": "KNOWN", "value": True}
    assert evaluate_pm1(payload).status == "ACTIVATED"


def test_no_unpublished_strength_adjustment():
    assert evaluate_pm1(case(hotspot(pathogenic_count=1000))).strength == "MODERATE"
