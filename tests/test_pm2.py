"""Synthetic rule-boundary fixtures; these are not curated clinical variants."""

from copy import deepcopy

import pytest

from adaptive_pm2 import evaluate_pm2


def base(**changes):
    case = {
        "variant_id": "synthetic-1", "genome_assembly": "GRCh38",
        "normalized_variant_representation": "NC_000001.11:g.100A>G",
        "gene_id": "EXAMPLE", "disease_id": "MONDO:synthetic",
        "inheritance_mode": "AUTOSOMAL_DOMINANT",
        "population_records": [{
            "name": "gnomAD", "version": "4.1", "genome_assembly": "GRCh38",
            "normalized_variant_representation": "NC_000001.11:g.100A>G",
            "allele_count": 0, "allele_number": 100000,
            "identity_confirmed": True, "filter_pass": True,
            "site_callable": True, "technically_reliable": True,
        }],
    }
    case.update(changes)
    return case


def record_case(**changes):
    case = base()
    case["population_records"][0].update(changes)
    return case


def acadvl(**record_changes):
    case = base(gene_id="ACADVL", disease_id="MONDO:0008723", inheritance_mode="AUTOSOMAL_RECESSIVE")
    case["population_records"][0].update({
        "allele_count": 4, "allele_number": 100000,
        "ancestries_complete": True,
        "ancestries": [
            {"name": "African", "allele_count": 1, "allele_number": 10000},
            {"name": "European", "allele_count": 3, "allele_number": 90000},
        ],
    })
    case["population_records"][0].update(record_changes)
    return case


def f9(**record_changes):
    case = base(gene_id="F9", disease_id="MONDO:0010604", inheritance_mode="X_LINKED")
    case["population_records"][0].update({"allele_count": 2, "hemizygote_count": 0, "male_allele_number": 25000})
    case["population_records"][0].update(record_changes)
    return case


def test_confirmed_absence_and_optional_data():
    result = evaluate_pm2(base())
    assert (result.status, result.strength, result.decision_path) == ("ACTIVATED", "SUPPORTING", "CONFIRMED_ABSENCE")
    assert result.to_dict()["guideline"]["version"]
    assert evaluate_pm2(base(disease_prevalence=0.00001)).status == result.status


@pytest.mark.parametrize("change", [
    {"site_callable": False}, {"filter_pass": False}, {"technically_reliable": False},
    {"identity_confirmed": None}, {"sample_size_adequate": False},
    {"genome_assembly": "GRCh37"},
    {"normalized_variant_representation": "NC_000001.11:g.101A>G"},
    {"mapping_ambiguity": True}, {"low_complexity_region": True},
    {"segmental_duplication_region": True}, {"detection_sensitivity_adequate": False},
    {"sampling_bias_material": True}, {"population_substructure_material": True},
])
def test_unreliable_absence_is_undetermined(change):
    assert evaluate_pm2(record_case(**change)).status == "UNDETERMINED"


@pytest.mark.parametrize("change", [
    {"allele_count": None}, {"allele_number": None},
    {"allele_count": None, "allele_frequency": None},
])
def test_missing_population_quantities_do_not_become_zero(change):
    assert evaluate_pm2(record_case(**change)).status == "UNDETERMINED"


def test_missing_population_record():
    assert evaluate_pm2(base(population_records=[])).status == "UNDETERMINED"


def test_general_nonzero_frequency_is_not_absence():
    assert evaluate_pm2(record_case(allele_count=1)).status == "NOT_MET"
    recessive = record_case(allele_count=1)
    recessive["inheritance_mode"] = "AUTOSOMAL_RECESSIVE"
    assert evaluate_pm2(recessive).status == "UNDETERMINED"
    recessive["population_records"][0]["allele_count"] = 0
    assert evaluate_pm2(recessive).status == "ACTIVATED"


def test_two_sources_with_discordant_calls_require_reconciliation():
    case = base()
    second = deepcopy(case["population_records"][0])
    second.update(name="other", allele_count=1)
    case["population_records"].append(second)
    assert evaluate_pm2(case).status == "UNDETERMINED"


@pytest.mark.parametrize("flag", ["ba1_established", "bs1_established"])
def test_frequency_benign_criteria_block_pm2(flag):
    assert evaluate_pm2(base(**{flag: True})).status == "NOT_MET"


def test_bs2_requires_reconciliation():
    assert evaluate_pm2(base(bs2_established=True)).status == "UNDETERMINED"


def test_acadvl_nonzero_rarity_and_ancestry_exceedance():
    assert evaluate_pm2(acadvl()).status == "ACTIVATED"
    high = acadvl(ancestries=[
        {"name": "African", "allele_frequency": 0.002, "allele_number": 10000},
        {"name": "European", "allele_frequency": 0.00001, "allele_number": 90000},
    ])
    assert evaluate_pm2(high).status == "NOT_MET"


def test_acadvl_ancestry_completeness_and_allele_number():
    assert evaluate_pm2(acadvl(ancestries_complete=False)).status == "UNDETERMINED"
    assert evaluate_pm2(acadvl(ancestries=[{"name": "African", "allele_frequency": 0, "allele_number": 1000}])).status == "UNDETERMINED"


def test_founder_group_is_not_mistaken_for_a_continental_group():
    case = acadvl(allele_count=21, allele_number=20000, ancestries=[
        {"name": "European", "allele_count": 1, "allele_number": 10000},
        {"name": "Founder cohort", "continental": False, "founder_or_bottlenecked": True,
         "allele_count": 20, "allele_number": 10000},
    ])
    assert evaluate_pm2(case).status == "ACTIVATED"


def test_acadvl_homozygotes_need_expert_review():
    assert evaluate_pm2(acadvl(homozygote_count=1)).status == "UNDETERMINED"


def test_acadvl_zero_count_and_large_allele_number():
    assert evaluate_pm2(acadvl(allele_count=0, allele_number=100000, ancestries=[
        {"name": "European", "allele_count": 0, "allele_number": 90000},
    ])).status == "ACTIVATED"
    assert evaluate_pm2(acadvl(allele_count=0, allele_number=100000, ancestries=[])).status == "UNDETERMINED"


def test_f9_male_specific_requirement():
    assert evaluate_pm2(f9()).status == "ACTIVATED"
    assert evaluate_pm2(f9(hemizygote_count=1)).status == "NOT_MET"
    assert evaluate_pm2(f9(hemizygote_count=None)).status == "UNDETERMINED"
    assert evaluate_pm2(f9(male_allele_number=None)).status == "UNDETERMINED"


def test_scope_and_version_are_strict():
    assert evaluate_pm2(base(gene_id="F9", disease_id="MONDO:other")).status == "NOT_APPLICABLE"
    assert evaluate_pm2(acadvl()).guideline["version"] == "2.2"
    case = acadvl()
    case["vcep_specification_version"] = "1.0"
    assert evaluate_pm2(case).status == "NOT_APPLICABLE"


@pytest.mark.parametrize("changes", [
    {"variant_type": "INDEL"}, {"variant_type": "MULTIALLELIC"},
    {"variant_type": "SV"}, {"variant_type": "CNV"},
    {"variant_type": "MITOCHONDRIAL"},
])
def test_special_variant_class_needs_matching_validation(changes):
    assert evaluate_pm2(base(**changes)).status == "UNDETERMINED"


def test_indel_can_be_assessed_after_normalization():
    assert evaluate_pm2(base(variant_type="INDEL", indel_normalization_adequate=True)).status == "ACTIVATED"


def test_multiallelic_and_structural_classes_need_validated_matching():
    assert evaluate_pm2(base(variant_type="MULTIALLELIC", multiallelic_allele_resolved=True)).status == "ACTIVATED"
    assert evaluate_pm2(base(variant_type="SV", structural_variant_frequency_resolved=True)).status == "ACTIVATED"


def test_explicit_frequency_conflict_is_unresolved():
    assert evaluate_pm2(base(frequency_conflict=True)).status == "UNDETERMINED"
    assert evaluate_pm2(record_case(frequency_evidence_state="CONFLICTING")).status == "UNDETERMINED"
    assert evaluate_pm2(record_case(frequency_evidence_state="UNKNOWN")).status == "UNDETERMINED"


def test_invalid_counts_fail_fast():
    with pytest.raises(ValueError):
        evaluate_pm2(record_case(allele_count=-1))
