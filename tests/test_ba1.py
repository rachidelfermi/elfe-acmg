"""Synthetic decision fixtures; no case here is claimed as clinically curated."""

from dataclasses import replace

import pytest

from adaptive_ba1 import (BA1Evidence, EvidenceState, ExceptionCheck,
                          PopulationRecord, VariantIdentity, evaluate_ba1)


def base(*, frequency=0.06, count=600, number=10000, faf=None,
         gene=None, disease=None, variant_type="SNV"):
    variant = VariantIdentity("synthetic", "GRCh38", "1", 100, "A", "G", True,
                              variant_type, EvidenceState.KNOWN)
    record = PopulationRecord("gnomAD", "4.1", "GRCh38", "1", 100, "A", "G", True,
                              "AFR", True, frequency, count, number, faf, True,
                              EvidenceState.KNOWN, False, True, False, False, False)
    return BA1Evidence(variant, (record,), gene, disease,
                       ExceptionCheck(EvidenceState.KNOWN, False, "ClinGen SVI list", "2018"),
                       EvidenceState.KNOWN, gene == "LDLR")


@pytest.mark.parametrize("frequency,count,expected", [
    (0.06, 600, "ACTIVATED"), (0.04, 400, "NOT_MET"),
    (0.05, 500, "NOT_MET"), (0.0501, 501, "ACTIVATED"),
])
def test_general_frequency_boundary(frequency, count, expected):
    result = evaluate_ba1(base(frequency=frequency, count=count))
    assert result.status == expected
    assert result.strength == ("STAND_ALONE" if expected == "ACTIVATED" else None)


def test_exception_list():
    e = base()
    e = replace(e, variant=replace(e.variant, clingen_allele_id="CA172210"))
    assert evaluate_ba1(e).status == "NOT_APPLICABLE"


def test_verified_disease_exception():
    e = base(gene="LDLR", disease="MONDO:0007750", faf=0.01)
    e = replace(e, exception_check=ExceptionCheck(EvidenceState.KNOWN, True,
                "ClinGen FH VCEP", "1.2", "MONDO:0007750"))
    assert evaluate_ba1(e).status == "NOT_APPLICABLE"


def test_ldlr_uses_faf_and_inclusive_boundary():
    e = base(frequency=0.03, count=300, faf=0.005,
             gene="LDLR", disease="MONDO:0007750")
    r = evaluate_ba1(e)
    assert (r.status, r.frequency_used, r.threshold_used) == ("ACTIVATED", 0.005, 0.005)
    assert r.guideline["version"] == "1.2"


def test_ldlr_bundled_rule_needs_no_separate_modification_check():
    e = base(gene="LDLR", disease="MONDO:0007750", faf=0.005)
    e = replace(e, guideline_modification_check=EvidenceState.UNKNOWN,
                known_modification=None)
    assert evaluate_ba1(e).status == "ACTIVATED"


def test_ldlr_does_not_use_raw_af():
    e = base(frequency=0.06, count=600, faf=0.004,
             gene="LDLR", disease="MONDO:0007750")
    assert evaluate_ba1(e).status == "NOT_MET"


@pytest.mark.parametrize("change", [
    lambda e: replace(e, populations=()),
    lambda e: replace(e, populations=(replace(e.populations[0], allele_frequency=None,
                                         allele_count=None, allele_number=None),)),
    lambda e: replace(e, populations=(replace(e.populations[0], allele_count=None),)),
    lambda e: replace(e, populations=(replace(e.populations[0], allele_number=None),)),
    lambda e: replace(e, populations=(replace(e.populations[0], filter_pass=False),)),
    lambda e: replace(e, populations=(replace(e.populations[0], quality=EvidenceState.UNKNOWN),)),
    lambda e: replace(e, populations=(replace(e.populations[0], assembly="GRCh37"),)),
    lambda e: replace(e, populations=(replace(e.populations[0], reference="AT", normalized=False),)),
    lambda e: replace(e, populations=(replace(e.populations[0], allele_number=1000),)),
    lambda e: replace(e, populations=(replace(e.populations[0], founder_or_bottlenecked=True),)),
    lambda e: replace(e, populations=(replace(e.populations[0], mapping_ambiguous=True),)),
    lambda e: replace(e, populations=(replace(e.populations[0], callability_adequate=False),)),
    lambda e: replace(e, populations=(replace(e.populations[0], technical_limitation=True),)),
    lambda e: replace(e, exception_check=ExceptionCheck()),
    lambda e: replace(e, guideline_modification_check=EvidenceState.UNKNOWN),
    lambda e: replace(e, frequency_conflict=True),
    lambda e: replace(e, pathogenic_evidence_conflict=True),
    lambda e: replace(e, reduced_penetrance=True),
    lambda e: replace(e, late_onset=True),
    lambda e: replace(e, risk_allele=True),
    lambda e: replace(e, hypomorphic_allele=True),
    lambda e: replace(e, founder_pathogenic_context=True),
])
def test_unresolved_evidence_never_activates(change):
    assert evaluate_ba1(change(base())).status == "UNDETERMINED"


def test_missing_optional_fields_do_not_block():
    e = replace(base(), evidence_source=None, evidence_version=None,
                missing_reason=None)
    assert evaluate_ba1(e).status == "ACTIVATED"


def test_irrelevant_provenance_does_not_change_decision():
    e = base()
    assert evaluate_ba1(e).status == evaluate_ba1(replace(e, evidence_source="lab")).status


@pytest.mark.parametrize("variant_type", ["SV", "CNV", "MITOCHONDRIAL"])
def test_general_class_exclusion(variant_type):
    e = base(variant_type=variant_type)
    assert evaluate_ba1(e).status == "NOT_APPLICABLE"


def test_ldlr_sv_requires_validated_sv_data():
    e = base(gene="LDLR", disease="MONDO:0007750", faf=0.006, variant_type="SV")
    assert evaluate_ba1(e).status == "UNDETERMINED"
    rec = replace(e.populations[0], population_dataset="gnomAD SV",
                  detection_sensitivity_adequate=True)
    assert evaluate_ba1(replace(e, populations=(rec,))).status == "ACTIVATED"


def test_sex_chromosome_requires_ploidy_pathway():
    e = base()
    variant = replace(e.variant, chromosome="X", sex_chromosome_ploidy="male")
    record = replace(e.populations[0], chromosome="X")
    assert evaluate_ba1(replace(e, variant=variant, populations=(record,))).status == "UNDETERMINED"
    variant = replace(variant, pseudoautosomal_region=True)
    assert evaluate_ba1(replace(e, variant=variant, populations=(record,))).status == "ACTIVATED"


def test_mitochondrial_chromosome_is_never_general_ba1():
    e = base()
    variant = replace(e.variant, chromosome="MT", mitochondrial_heteroplasmy=0.5)
    record = replace(e.populations[0], chromosome="MT")
    assert evaluate_ba1(replace(e, variant=variant, populations=(record,))).status == "NOT_APPLICABLE"


def test_vcep_database_requirement():
    e = base(gene="LDLR", disease="MONDO:0007750", faf=0.006)
    rec = replace(e.populations[0], database="OtherDB")
    assert evaluate_ba1(replace(e, populations=(rec,))).status == "UNDETERMINED"


def test_wrong_disease_never_borrows_ldlr_threshold():
    e = base(gene="LDLR", disease="MONDO:other", faf=0.006)
    assert evaluate_ba1(e).status == "UNDETERMINED"


def test_high_af_with_no_faf_does_not_activate_ldlr():
    e = base(gene="LDLR", disease="MONDO:0007750")
    assert evaluate_ba1(e).status == "UNDETERMINED"


def test_general_uses_qualifying_ancestry_not_global_only():
    e = base(frequency=0.01, count=100)
    high = replace(e.populations[0], ancestry="EAS", allele_frequency=0.06,
                   allele_count=600)
    assert evaluate_ba1(replace(e, populations=e.populations + (high,))).status == "ACTIVATED"


def test_no_invented_frequency_from_missing_counts():
    e = base()
    rec = replace(e.populations[0], allele_frequency=None,
                  allele_count=None, allele_number=10000)
    assert evaluate_ba1(replace(e, populations=(rec,))).status == "UNDETERMINED"
