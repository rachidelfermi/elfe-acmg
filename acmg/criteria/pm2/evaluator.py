"""Deterministic PM2 decision paths over supplied evidence."""

from __future__ import annotations

from .guidelines import ACADVL, F9, Guideline, select_guideline
from .models import EvidenceState, PM2Input, PM2Result, PopulationRecord, parse_input


def _result(
    status: str, rule: Guideline, reason: str, *, path: str | None = None,
    records: tuple[PopulationRecord, ...] = (), used: tuple[str, ...] = (),
    missing: tuple[str, ...] = (), exclusions: tuple[str, ...] = (),
) -> PM2Result:
    return PM2Result(
        criterion="PM2", status=status,
        strength="SUPPORTING" if status == "ACTIVATED" else "NONE",
        decision_path=path, population_sources_used=tuple(f"{r.name} {r.version}" for r in records),
        features_used=used, missing_critical_features=missing,
        exclusions_evaluated=exclusions, guideline=rule.citation(), reason=reason,
    )


def _valid_record(record: PopulationRecord, features: PM2Input) -> tuple[bool, tuple[str, ...]]:
    problems = []
    if record.genome_assembly != features.genome_assembly:
        problems.append("genome_assembly mismatch")
    if record.normalized_variant_representation != features.normalized_variant_representation:
        problems.append("normalized_variant_representation mismatch")
    for name in ("identity_confirmed", "filter_pass", "site_callable", "technically_reliable"):
        if getattr(record, name) is not True:
            problems.append(name)
    for name in (
        "detection_sensitivity_adequate", "genotype_quality_adequate", "sample_size_adequate",
    ):
        if getattr(record, name) is False:
            problems.append(name)
    for name in (
        "low_complexity_region", "segmental_duplication_region", "mapping_ambiguity",
        "sampling_bias_material", "case_enrichment_material", "population_substructure_material",
        "frequency_conflict",
    ):
        if getattr(record, name) is True:
            problems.append(name)
    return not problems, tuple(problems)


def _validate_numbers(record: PopulationRecord) -> None:
    for name in ("allele_count", "allele_number", "homozygote_count", "hemizygote_count", "male_allele_number"):
        value = getattr(record, name)
        if value is not None and (not isinstance(value, int) or isinstance(value, bool) or value < 0):
            raise ValueError(f"{name} must be a nonnegative integer")
    for name in ("allele_frequency", "filtering_allele_frequency"):
        value = getattr(record, name)
        if value is not None and not 0 <= value <= 1:
            raise ValueError(f"{name} must be between zero and one")
    if record.allele_count is not None and record.allele_number is not None and record.allele_count > record.allele_number:
        raise ValueError("allele_count exceeds allele_number")
    if record.allele_count == 0 and record.allele_frequency is not None and record.allele_frequency > 0:
        raise ValueError("zero allele count conflicts with nonzero allele frequency")
    for ancestry in record.ancestries:
        if ancestry.allele_number is not None and ancestry.allele_number < 0:
            raise ValueError("ancestry allele_number must be nonnegative")
        if ancestry.allele_count is not None and ancestry.allele_count < 0:
            raise ValueError("ancestry allele_count must be nonnegative")
        if ancestry.allele_count is not None and ancestry.allele_number is not None and ancestry.allele_count > ancestry.allele_number:
            raise ValueError("ancestry allele_count exceeds allele_number")
        if ancestry.allele_frequency is not None and not 0 <= ancestry.allele_frequency <= 1:
            raise ValueError("ancestry allele_frequency must be between zero and one")


def _observed_frequency(record: PopulationRecord) -> float | None:
    if record.allele_frequency is not None:
        return record.allele_frequency
    if record.allele_count is not None and record.allele_number:
        return record.allele_count / record.allele_number
    return None


def evaluate_pm2(biological_features: PM2Input | dict) -> PM2Result:
    """Evaluate only PM2; unknown evidence never becomes a zero observation."""
    features = parse_input(biological_features)
    rule, scope_problem = select_guideline(features)
    if scope_problem:
        return _result("NOT_APPLICABLE", rule, scope_problem)
    if not features.variant_id or not features.genome_assembly or not features.normalized_variant_representation:
        return _result("UNDETERMINED", rule, "Exact variant identity is incomplete.",
                       missing=("variant_id", "genome_assembly", "normalized_variant_representation"))
    if features.gene_disease_validity is False:
        return _result("NOT_APPLICABLE", rule, "Gene-disease relationship is explicitly invalid.")
    if features.ba1_established or features.bs1_established:
        return _result("NOT_MET", rule, "BA1 or BS1 is established for the same variant and disease context.",
                       exclusions=("BA1", "BS1"))
    if features.bs2_established:
        return _result("UNDETERMINED", rule, "BS2 evidence needs expert reconciliation with PM2.",
                       exclusions=("BS2",))
    if features.frequency_conflict:
        return _result("UNDETERMINED", rule, "Population frequency evidence is conflicting.",
                       missing=("frequency_conflict_resolution",))
    if features.variant_mapping_confidence in ("LOW", "CONFLICTING"):
        return _result("UNDETERMINED", rule, "Variant mapping is not reliable.", missing=("variant_mapping_confidence",))
    if features.variant_type.upper() in {"INDEL", "INSERTION", "DELETION"} and features.indel_normalization_adequate is not True:
        return _result("UNDETERMINED", rule, "Indel normalization needs confirmation.", missing=("indel_normalization_adequate",))
    if features.variant_type.upper() == "MULTIALLELIC" and features.multiallelic_allele_resolved is not True:
        return _result("UNDETERMINED", rule, "The alternate allele must be resolved at this multiallelic site.",
                       missing=("multiallelic_allele_resolved",))
    if features.variant_type.upper() in {"SV", "CNV"} and features.structural_variant_frequency_resolved is not True:
        return _result("UNDETERMINED", rule, "Structural variant frequency and detection require validation.",
                       missing=("structural_variant_frequency_resolved",))
    if features.variant_type.upper() == "MITOCHONDRIAL" and features.mitochondrial_heteroplasmy_resolved is not True:
        return _result("UNDETERMINED", rule, "Mitochondrial heteroplasmy interpretation is unresolved.",
                       missing=("mitochondrial_heteroplasmy_resolved",))
    if rule is F9 and features.ploidy_context_resolved is False:
        return _result("UNDETERMINED", rule, "X-linked ploidy context is unresolved.",
                       missing=("ploidy_context_resolved",))

    valid = []
    problems = []
    for record in features.population_records:
        _validate_numbers(record)
        if record.frequency_evidence_state == EvidenceState.CONFLICTING:
            return _result("UNDETERMINED", rule, "A population source explicitly reports conflicting frequency evidence.",
                           records=(record,), missing=("frequency_conflict_resolution",))
        if record.frequency_evidence_state in (EvidenceState.UNKNOWN, EvidenceState.NOT_APPLICABLE):
            problems.append("frequency_evidence_state")
            continue
        good, why = _valid_record(record, features)
        if good:
            valid.append(record)
        else:
            problems.extend(why)
    if not valid:
        return _result("UNDETERMINED", rule, "No matched, callable, reliable population record is available.",
                       missing=tuple(sorted(set(problems))) if problems else ("population_records",))
    records = tuple(valid)
    if rule is ACADVL:
        return _evaluate_acadvl(records, rule)
    if rule is F9:
        return _evaluate_f9(records, rule)
    return _evaluate_general(records, features, rule)


def _evaluate_general(records: tuple[PopulationRecord, ...], features: PM2Input, rule: Guideline) -> PM2Result:
    zero = []
    observed = []
    for record in records:
        if record.allele_count is not None and record.allele_number is not None and record.allele_number > 0:
            if record.allele_count == 0:
                zero.append(record)
            else:
                observed.append(record)
        elif _observed_frequency(record) is not None and _observed_frequency(record) > 0:
            observed.append(record)
    if zero and observed:
        return _result("UNDETERMINED", rule, "One source reports absence while another observes the variant; reconcile cohorts and calls.",
                       records=records, missing=("population_source_reconciliation",))
    if observed:
        if (features.inheritance_mode or "").upper() == "AUTOSOMAL_RECESSIVE":
            return _result("UNDETERMINED", rule, "A nonzero recessive frequency needs a verified disease-specific rarity rule.",
                           records=tuple(observed), missing=("verified_recessive_frequency_rule",))
        return _result("NOT_MET", rule, "The variant is observed; the general absence pathway is not met.",
                       records=tuple(observed), used=("allele_count_or_frequency",))
    if zero:
        return _result("ACTIVATED", rule, "Confirmed absence (AC=0, AN>0) in a matched, callable, reliable population source.",
                       path="CONFIRMED_ABSENCE", records=tuple(zero),
                       used=("normalized_variant_representation", "genome_assembly", "allele_count", "allele_number",
                             "identity_confirmed", "filter_pass", "site_callable", "technically_reliable"),
                       exclusions=("BA1", "BS1", "BS2"))
    return _result("UNDETERMINED", rule, "Neither confirmed absence nor a usable frequency observation is available.",
                   records=records, missing=("allele_count", "allele_number"))


def _evaluate_acadvl(records: tuple[PopulationRecord, ...], rule: Guideline) -> PM2Result:
    gnomad = tuple(r for r in records if r.name.lower() == "gnomad")
    if not gnomad:
        return _result("UNDETERMINED", rule, "The pinned ACADVL rule requires gnomAD population evidence.",
                       missing=("gnomAD",))
    if len(gnomad) > 1:
        return _result("UNDETERMINED", rule, "Select one gnomAD version to avoid counting overlapping samples.",
                       records=gnomad, missing=("selected_gnomAD_version",))
    record = gnomad[0]
    if record.ancestries_complete is not True:
        return _result("UNDETERMINED", rule, "The highest eligible continental ancestry frequency cannot be established.",
                       records=(record,), missing=("ancestries_complete",))
    eligible = [a for a in record.ancestries if a.continental and not a.founder_or_bottlenecked
                and a.allele_number is not None and a.allele_number > 2000]
    if not eligible:
        return _result("UNDETERMINED", rule, "No ancestry group has more than 2000 alleles.",
                       records=(record,), missing=("eligible_ancestry_allele_number",))
    if record.allele_count == 0 and record.allele_number is not None and record.allele_number > 0:
        if any(a.allele_count not in (None, 0) or a.allele_frequency not in (None, 0) for a in eligible):
            return _result("UNDETERMINED", rule, "Ancestry observations conflict with the zero global allele count.",
                           records=(record,), missing=("population_source_reconciliation",))
        return _result("ACTIVATED", rule, "Confirmed gnomAD absence with an eligible continental ancestry group.",
                       path="VCEP_ACADVL_ABSENCE", records=(record,),
                       used=("allele_count", "allele_number", "ancestry_allele_number", "site_callable"),
                       exclusions=("BA1", "BS1", "BS2"))
    frequencies = []
    for ancestry in eligible:
        if ancestry.technically_reliable is False:
            return _result("UNDETERMINED", rule, f"The {ancestry.name} ancestry call is unreliable.", records=(record,),
                           missing=("ancestry_call_reliability",))
        frequency = ancestry.allele_frequency
        if frequency is None and ancestry.allele_count is not None:
            frequency = ancestry.allele_count / ancestry.allele_number
        if frequency is None:
            return _result("UNDETERMINED", rule, f"Frequency is missing for {ancestry.name}.", records=(record,),
                           missing=("ancestry_allele_frequency",))
        frequencies.append(frequency)
    highest = max(frequencies)
    if highest >= 0.001:
        return _result("NOT_MET", rule, f"Highest eligible ancestry frequency {highest:g} is not below 0.001.",
                       records=(record,), used=("ancestry_allele_frequency", "ancestry_allele_number"))
    if record.homozygote_count:
        return _result("UNDETERMINED", rule, "ACADVL homozygotes require expert discussion under the VCEP specification.",
                       records=(record,), missing=("homozygote_expert_review",))
    return _result("ACTIVATED", rule, f"Highest eligible ancestry frequency {highest:g} is below 0.001.",
                   path="VCEP_ACADVL_RARITY", records=(record,),
                   used=("ancestry_allele_frequency", "ancestry_allele_number", "ancestries_complete"),
                   exclusions=("BA1", "BS1", "BS2"))


def _evaluate_f9(records: tuple[PopulationRecord, ...], rule: Guideline) -> PM2Result:
    if any(r.hemizygote_count is not None and r.hemizygote_count > 0 for r in records):
        return _result("NOT_MET", rule, "The variant is observed in males; F9 PM2 requires male absence.",
                       records=records, used=("hemizygote_count",))
    confirmed = tuple(r for r in records if r.hemizygote_count == 0 and r.male_allele_number is not None and r.male_allele_number > 0)
    if not confirmed:
        return _result("UNDETERMINED", rule, "Male absence is unconfirmed; female or global counts cannot substitute.",
                       records=records, missing=("hemizygote_count", "male_allele_number"))
    return _result("ACTIVATED", rule, "No male hemizygotes among callable male alleles in a reliable source.",
                   path="VCEP_F9_MALE_ABSENCE", records=confirmed,
                   used=("hemizygote_count", "male_allele_number", "site_callable"),
                   exclusions=("BA1", "BS1", "BS2"))
