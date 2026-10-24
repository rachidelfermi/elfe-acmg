from dataclasses import replace

from adaptive_pm1 import Evidence, evaluate_pm1
from .test_hotspots import case, hotspot


def domain(**changes):
    r = replace(hotspot(), region_id="critical-domain", kind="CRITICAL_DOMAIN",
                pathogenic_clustering=Evidence(), critical_function=Evidence.known(True))
    return replace(r, **changes)


def test_established_critical_domain():
    assert evaluate_pm1(case(domain())).decision_path == "ESTABLISHED_CRITICAL_DOMAIN"


def test_generic_domain_annotation_is_insufficient():
    assert evaluate_pm1(case(replace(domain(), kind="GENERIC_DOMAIN"))).status == "NOT_MET"


def test_domain_with_benign_variation():
    assert evaluate_pm1(case(domain(lacks_benign_variation=Evidence.known(False)))).status == "NOT_MET"


def test_domain_without_criticality_is_unresolved():
    assert evaluate_pm1(case(domain(critical_function=Evidence()))).status == "UNDETERMINED"


def test_other_complete_region_can_activate():
    result = evaluate_pm1(case(domain(critical_function=Evidence()), hotspot()))
    assert result.status == "ACTIVATED"


def test_overlapping_generic_domain_does_not_block():
    result = evaluate_pm1(case(replace(domain(), region_id="generic-domain", kind="GENERIC_DOMAIN"), domain()))
    assert result.status == "ACTIVATED"


def test_enrichment_has_no_invented_universal_threshold():
    r = replace(hotspot(), kind="ENRICHMENT", enrichment_statistic=999.0)
    assert evaluate_pm1(case(r)).status == "UNDETERMINED"
