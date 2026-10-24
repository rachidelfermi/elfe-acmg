"""Synthetic fixtures; source thresholds are independently documented in README."""

from dataclasses import replace

import pytest

from adaptive_bs1 import (
    BS1Features, EvidenceState, FrequencyRecord, GuidelineRule, Status,
    Strength, VariantIdentity, evaluate_bs1,
)


V = VariantIdentity("GRCh38", "15", 63000001, "A", "G")


def tpm1(frequency=0.0002, **record_changes):
    record = FrequencyRecord(V, "gnomAD", "4.1", "POPMAX_FAF", frequency,
                             filter_status="PASS", quality=EvidenceState.KNOWN,
                             platform="genome")
    return BS1Features(V, "TPM1", "MONDO:0005045", "AD",
                       (replace(record, **record_changes),))


def rmrp(frequency=0.001, **record_changes):
    record = FrequencyRecord(V, "gnomAD", "4.1", "POPMAX_FAF", frequency,
                             filter_status="PASS", quality=EvidenceState.KNOWN,
                             platform="exome", depth=25)
    return BS1Features(V, "RMRP", "MONDO:0009595", "AR",
                       (replace(record, **record_changes),))


@pytest.mark.parametrize("frequency,status", [
    (0.0002, Status.ACTIVATED), (0.00009, Status.NOT_MET),
    (0.0001, Status.ACTIVATED),
])
def test_tpm1_inclusive_operator(frequency, status):
    result = evaluate_bs1(tpm1(frequency))
    assert result.status == status
    assert result.strength == (Strength.STRONG if status == Status.ACTIVATED else Strength.NONE)


@pytest.mark.parametrize("frequency,status", [
    (0.001, Status.ACTIVATED), (0.00089, Status.NOT_MET),
    (0.0005, Status.NOT_MET),
])
def test_rmrp_strict_operator(frequency, status):
    assert evaluate_bs1(rmrp(frequency)).status == status


def test_no_threshold_and_no_invention():
    assert evaluate_bs1(replace(tpm1(), gene="OTHER")).status == Status.UNDETERMINED
    assert evaluate_bs1(replace(tpm1(), gene="OTHER", prevalence=0.01,
                                penetrance=1)).status == Status.UNDETERMINED


def test_missing_prevalence_penetrance_and_unrelated_features():
    assert evaluate_bs1(tpm1()).status == Status.ACTIVATED
    assert evaluate_bs1(replace(tpm1(), reduced_penetrance=True,
                                late_onset=True, bs2_established=True)).status == Status.ACTIVATED


@pytest.mark.parametrize("changes", [
    {"frequency": None}, {"quality": EvidenceState.UNKNOWN},
    {"filter_status": "FAIL"}, {"technical_limitation": True},
    {"metric": "GLOBAL_AF"}, {"database": "other"},
])
def test_incomplete_population_evidence(changes):
    assert evaluate_bs1(tpm1(**changes)).status == Status.UNDETERMINED


def test_variant_and_build_identity():
    other_build = replace(V, assembly="GRCh37")
    assert evaluate_bs1(tpm1(variant=other_build)).status == Status.UNDETERMINED
    indel = replace(V, reference="AT", alternate="A", variant_type="INDEL", normalized=False)
    case = replace(tpm1(), variant=indel,
                   frequencies=(replace(tpm1().frequencies[0], variant=indel),))
    assert evaluate_bs1(case).status == Status.UNDETERMINED


def test_rmrp_exome_depth():
    assert evaluate_bs1(rmrp(depth=19)).status == Status.UNDETERMINED
    assert evaluate_bs1(rmrp(depth=20)).status == Status.ACTIVATED
    assert evaluate_bs1(rmrp(platform="genome", depth=None)).status == Status.ACTIVATED
    assert evaluate_bs1(rmrp(platform=None)).status == Status.UNDETERMINED


@pytest.mark.parametrize("field", ["founder_variant", "risk_allele", "hypomorphic_allele"])
def test_special_context_requires_exception_review(field):
    case = replace(tpm1(), **{field: True})
    assert evaluate_bs1(case).status == Status.UNDETERMINED
    assert evaluate_bs1(replace(case, exception_status=EvidenceState.NOT_APPLICABLE)).status == Status.ACTIVATED


def test_sourced_exception():
    case = replace(tpm1(), exception_status=EvidenceState.KNOWN,
                   exception_source="https://example.org/curated-exception",
                   exception_disease_scope="MONDO:0005045")
    assert evaluate_bs1(case).status == Status.NOT_APPLICABLE
    assert evaluate_bs1(replace(case, exception_source=None)).status == Status.UNDETERMINED
    assert evaluate_bs1(replace(case, exception_disease_scope="OTHER")).status == Status.UNDETERMINED


def test_founder_population_record_not_used_for_popmax():
    assert evaluate_bs1(tpm1(founder_population=True)).status == Status.UNDETERMINED


def test_conflicting_sources_and_other_population_criteria():
    case = tpm1()
    low = replace(case.frequencies[0], frequency=0.00005, version="2.1.1")
    assert evaluate_bs1(replace(case, frequencies=case.frequencies + (low,))).status == Status.UNDETERMINED
    assert evaluate_bs1(replace(case, ba1_established=True)).status == Status.NOT_APPLICABLE
    assert evaluate_bs1(replace(case, pm2_established=True)).status == Status.UNDETERMINED
    assert evaluate_bs1(replace(case, bs2_established=True)).status == Status.ACTIVATED


def test_custom_attested_rule_strength_and_ancestry_requirements():
    rule = GuidelineRule("Sourced example", "1", "EXAMPLE:BS1", "https://example.org/rule",
                         "GENE", "DISEASE", "AR", 0.002, "ANCESTRY_FAF", ">",
                         strength=Strength.SUPPORTING, minimum_allele_number=10000,
                         eligible_ancestries=frozenset({"EUR"}), verified=True)
    record = FrequencyRecord(V, "cohort", "1", "ANCESTRY_FAF", 0.003,
                             ancestry="EUR", allele_number=20000,
                             filter_status="PASS", quality=EvidenceState.KNOWN,
                             sampling=EvidenceState.KNOWN)
    case = BS1Features(V, "GENE", "DISEASE", "AR", (record,), guideline=rule)
    assert evaluate_bs1(case).strength == Strength.SUPPORTING
    assert evaluate_bs1(replace(case, frequencies=(replace(record, allele_number=20),))).status == Status.UNDETERMINED
    assert evaluate_bs1(replace(case, frequencies=(replace(record, sampling=EvidenceState.UNKNOWN),))).status == Status.UNDETERMINED
    assert evaluate_bs1(replace(case, frequencies=(replace(record, ancestry="AFR"),))).status == Status.UNDETERMINED
    assert evaluate_bs1(replace(case, gene="OTHER")).status == Status.NOT_APPLICABLE
    assert evaluate_bs1(replace(case, guideline=replace(rule, verified=False))).status == Status.UNDETERMINED


@pytest.mark.parametrize("changes", [
    {"variant_mapping": EvidenceState.UNKNOWN},
    {"gene_disease_validity": EvidenceState.CONFLICTING},
    {"conflicting_evidence": True},
    {"exception_status": EvidenceState.CONFLICTING},
])
def test_unresolved_essential_evidence(changes):
    assert evaluate_bs1(replace(tpm1(), **changes)).status == Status.UNDETERMINED


def test_unsupported_variant_class_and_scope():
    sv = replace(V, variant_type="SV")
    case = replace(tpm1(), variant=sv, frequencies=(replace(tpm1().frequencies[0], variant=sv),))
    assert evaluate_bs1(case).status == Status.NOT_APPLICABLE
    assert evaluate_bs1(replace(tpm1(), disease="MONDO:other")).status == Status.UNDETERMINED
