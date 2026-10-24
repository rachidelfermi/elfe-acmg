"""Deterministic, conservative PS4 decisions from already-curated evidence."""

from .guidelines import ACMG, MYOC, myoc_scope
from .models import PS4Input, PS4Result, Study
from .statistics import fisher_or_ci


def _result(status: str, reason: str, **kwargs) -> PS4Result:
    return PS4Result(status=status, reason=reason, **kwargs)


def _study_result(study: Study, variant: str, disease: str) -> PS4Result:
    missing = []
    conflicts = []
    exclusions = []
    if study.variant_id is None:
        missing.append("study.variant_id")
    elif study.variant_id != variant:
        return _result("NOT_APPLICABLE", "Study concerns another variant.", guideline=ACMG)
    if study.disease_id is None:
        missing.append("study.disease_id")
    elif study.disease_id != disease:
        return _result("NOT_APPLICABLE", "Study concerns another disease.", guideline=ACMG)

    for name in ("cases_affected", "controls_valid", "variant_calls_valid", "technical_comparable", "independent_cases"):
        evidence = getattr(study, name)
        if evidence.state == "CONFLICTING":
            conflicts.append(name)
        elif evidence.state != "KNOWN":
            missing.append(name)
        elif evidence.value is False:
            exclusions.append(name)
    if study.overlapping_cohort.state == "CONFLICTING":
        conflicts.append("overlapping_cohort")
    elif study.overlapping_cohort.known(True):
        exclusions.append("overlapping_cohort")
    elif study.overlapping_cohort.state != "KNOWN":
        missing.append("overlapping_cohort")

    if study.ancestry_comparable.state == "CONFLICTING":
        conflicts.append("ancestry_comparable")
    elif study.ancestry_comparable.known(False) and not study.stratification_adjusted.known(True):
        exclusions.append("ancestry_mismatch_unadjusted")
    elif not (study.ancestry_comparable.known(True) or study.stratification_adjusted.known(True)):
        missing.append("ancestry_comparable_or_stratification_adjusted")
    for name in ("founder_effect", "batch_effect"):
        evidence = getattr(study, name)
        if evidence.state == "CONFLICTING":
            conflicts.append(name)
        elif evidence.known(True) and not study.confounding_resolved.known(True):
            exclusions.append(name + "_unresolved")
    if study.confounding_resolved.state == "CONFLICTING":
        conflicts.append("confounding_resolved")

    counts = (study.case_carriers, study.case_noncarriers,
              study.control_carriers, study.control_noncarriers)
    has_counts = all(x is not None for x in counts)
    any_counts = any(x is not None for x in counts)
    if any_counts and not has_counts:
        missing.append("complete_case_control_2x2_counts")
    raw_or = raw_ci = p_value = None
    adjusted = study.statistic_method == "adjusted" and study.odds_ratio is not None
    if study.odds_ratio is not None and study.relative_risk is not None:
        conflicts.append("mixed_OR_and_RR_estimates")
    if has_counts and study.relative_risk is not None:
        conflicts.append("RR_with_case_control_sampling")
    if has_counts:
        raw_or, raw_ci, p_value = fisher_or_ci(*counts)
        if study.odds_ratio is not None and not adjusted:
            if (raw_or == float("inf") and study.odds_ratio != float("inf")
                    or raw_or == 0 and study.odds_ratio != 0
                    or raw_or not in (0, float("inf")) and abs(study.odds_ratio / raw_or - 1) > 0.05):
                conflicts.append("reported_odds_ratio_vs_raw_counts")
        if study.confidence_interval is not None and not adjusted:
            reported_lo, reported_hi = study.confidence_interval
            if not reported_lo <= raw_or <= reported_hi:
                conflicts.append("reported_interval_excludes_raw_odds_ratio")
        if adjusted and (study.confidence_interval is None or study.confidence_level != 0.95
                         or not study.statistic_source
                         or not study.confounding_resolved.known(True)):
            missing.append("validated_adjusted_estimate")
    else:
        if study.odds_ratio is None and study.relative_risk is None:
            missing.append("odds_ratio_or_relative_risk")
        if study.confidence_interval is None:
            missing.append("confidence_interval")
        if study.confidence_level != 0.95:
            missing.append("95_percent_confidence_level")
        if not study.statistic_source:
            missing.append("statistic_source")
        if not study.statistic_method:
            missing.append("statistic_method")
        if study.statistic_method == "adjusted" and not study.confounding_resolved.known(True):
            missing.append("adjustment_validity")

    odds = study.odds_ratio if adjusted or not has_counts else raw_or
    rr = study.relative_risk if not has_counts else None
    ci = study.confidence_interval if adjusted or not has_counts else raw_ci
    effect = odds if odds is not None else rr
    if effect is not None and effect < 0:
        raise ValueError("association estimate must be nonnegative")
    if ci is not None:
        if len(ci) != 2 or ci[0] < 0 or ci[0] > ci[1]:
            raise ValueError("invalid confidence interval")
        if effect is not None and not (ci[0] <= effect <= ci[1]):
            conflicts.append("effect_outside_confidence_interval")
    if study.p_value is not None and not 0 <= study.p_value <= 1:
        raise ValueError("p_value must be between zero and one")
    if has_counts and study.statistic_method == "fisher_exact" and study.p_value is not None and abs(study.p_value - p_value) > 1e-4:
        conflicts.append("reported_p_value_vs_fisher_exact")

    common = dict(decision_path="CASE_CONTROL_ASSOCIATION" if has_counts else "PUBLISHED_ASSOCIATION",
                  studies_used=(study.study_id,), odds_ratio=odds, relative_risk=rr, confidence_interval=ci,
                  p_value=study.p_value if adjusted else p_value if has_counts else study.p_value,
                  guideline=ACMG, missing_critical_features=tuple(dict.fromkeys(missing)),
                  conflicting_evidence=tuple(conflicts), exclusions_evaluated=tuple(exclusions))
    if conflicts:
        return _result("UNDETERMINED", "Material evidence conflicts require review.", **common)
    if exclusions:
        return _result("UNDETERMINED", "Study validity or confounding prevents using the association.", **common)
    if missing:
        return _result("UNDETERMINED", "Essential study evidence is unavailable.", **common)
    if effect is None or ci is None:
        return _result("UNDETERMINED", "Effect estimate and uncertainty are required.", **common)
    if effect <= 1 or ci[1] <= 1:
        return _result("NOT_MET", "No positive case enrichment was established.", **common)
    if effect > 5 and ci[0] > 1:
        kind = "RR" if rr is not None else "OR"
        return _result("ACTIVATED", f"{kind} exceeds 5 and its 95% confidence interval excludes 1.",
                       strength="STRONG", features_used=("relative_risk" if rr is not None else "odds_ratio", "confidence_interval", "cohort_validity"), **common)
    return _result("NOT_MET", "Complete valid study does not satisfy ACMG/AMP PS4 Note 1.", **common)


def _myoc_result(data: PS4Input) -> PS4Result:
    common = dict(decision_path="MYOC_PROBAND_COUNT", guideline=MYOC)
    if not data.probands:
        return _result("UNDETERMINED", "No individually identified probands were provided.",
                       missing_critical_features=("probands",), **common)
    if data.pm2_met.state == "CONFLICTING":
        return _result("UNDETERMINED", "PM2 evidence conflicts.", conflicting_evidence=("pm2_met",), **common)
    if data.pm2_met.known(False):
        return _result("NOT_MET", "The MYOC PS4 proband pathway requires PM2.", **common)
    if not data.pm2_met.known(True):
        return _result("UNDETERMINED", "MYOC PM2 status is unknown.", missing_critical_features=("pm2_met",), **common)
    if data.studies_independent.state == "CONFLICTING":
        return _result("UNDETERMINED", "Study independence conflicts.", conflicting_evidence=("studies_independent",), **common)
    if not data.studies_independent.known(True):
        return _result("UNDETERMINED", "Independent study provenance is required.",
                       missing_critical_features=("studies_independent",), **common)
    seen_cases, families, studies = set(), set(), set()
    missing, conflicts, exclusions = [], [], []
    for p in data.probands:
        if p.case_id in seen_cases:
            conflicts.append("duplicate_case_id:" + p.case_id)
            continue
        seen_cases.add(p.case_id)
        if not p.study_id or not p.family_id:
            missing.append("study_id_or_family_id:" + p.case_id)
            continue
        if p.family_id in families:
            conflicts.append("related_probands:" + p.family_id)
            continue
        if p.clinically_assessed.state != "KNOWN" or p.phenotype is None:
            missing.append("clinical_assessment_or_phenotype:" + p.case_id)
            continue
        if not p.clinically_assessed.known(True) or p.phenotype not in ("JOAG", "POAG"):
            exclusions.append("ineligible_phenotype:" + p.case_id)
            continue
        if p.has_other_myoc_vus_lp_p.state != "KNOWN":
            missing.append("other_MYOC_variants:" + p.case_id)
            continue
        if p.has_other_myoc_vus_lp_p.known(True):
            exclusions.append("multiple_MYOC_variants:" + p.case_id)
            continue
        if p.overlaps_other_criterion.known(True) or p.overlaps_other_criterion.state == "CONFLICTING":
            conflicts.append("criterion_overlap:" + p.case_id)
            continue
        families.add(p.family_id)
        studies.add(p.study_id)
    n = len(families)
    common.update(eligible_case_count=n, studies_used=tuple(sorted(studies)),
                  missing_critical_features=tuple(missing), conflicting_evidence=tuple(conflicts),
                  exclusions_evaluated=tuple(exclusions))
    if conflicts:
        return _result("UNDETERMINED", "Case independence or criterion overlap is unresolved.", **common)
    if missing and (n < 15 or len(studies) < 2):
        return _result("UNDETERMINED", "Some probands lack essential eligibility evidence.", **common)
    if len(studies) < 2:
        return _result("UNDETERMINED", "MYOC PS4 requires probands from multiple independent studies.",
                       **common)
    strength = "STRONG" if n >= 15 else "MODERATE" if n >= 6 else "SUPPORTING" if n >= 2 else "NONE"
    if strength != "NONE":
        return _result("ACTIVATED", f"{n} eligible unrelated MYOC probands from multiple studies meet the VCEP threshold.",
                       strength=strength, features_used=("probands", "pm2_met", "clinical_phenotype", "study_independence"), **common)
    return _result("NOT_MET", "Fewer than two eligible unrelated probands.", **common)


def evaluate_ps4(biological_features: PS4Input | dict) -> PS4Result:
    """Evaluate PS4 only; no variant annotation, external data, or ACMG aggregation."""
    data = PS4Input.parse(biological_features)
    if data.gene_disease_validity.state == "CONFLICTING":
        return _result("UNDETERMINED", "Gene-disease validity conflicts.",
                       conflicting_evidence=("gene_disease_validity",), guideline=ACMG)
    if data.gene_disease_validity.known(False):
        return _result("NOT_APPLICABLE", "Gene-disease association is invalid.", guideline=ACMG)
    if not data.gene_disease_validity.known(True):
        return _result("UNDETERMINED", "Gene-disease validity is unknown.",
                       missing_critical_features=("gene_disease_validity",), guideline=ACMG)
    if data.specification is not None and data.specification != "MYOC_2.1":
        return _result("UNDETERMINED", "Unsupported specification; no rules were imported.",
                       missing_critical_features=("supported_specification",), guideline=ACMG)
    if data.specification == "MYOC_2.1" and not myoc_scope(data.gene_id, data.disease_id):
        return _result("NOT_APPLICABLE", "MYOC VCEP rule is outside its gene-disease scope.", guideline=MYOC)
    if myoc_scope(data.gene_id, data.disease_id):
        return _myoc_result(data)
    if not data.studies:
        return _result("UNDETERMINED", "No applicable case-control study or VCEP case-count rule.",
                       missing_critical_features=("case_control_study",), guideline=ACMG)
    results = [_study_result(s, data.variant_id, data.disease_id) for s in data.studies]
    active = [r for r in results if r.status == "ACTIVATED"]
    if active:
        opposing = [r for r in results if r.status == "NOT_MET" and r.odds_ratio is not None and r.odds_ratio <= 1]
        if opposing:
            return _result("UNDETERMINED", "Independent studies have opposing association directions.",
                           conflicting_evidence=tuple(r.studies_used[0] for r in opposing),
                           studies_used=tuple(r.studies_used[0] for r in active + opposing), guideline=ACMG)
        # Each complete study stands alone; studies are never pooled or counted twice.
        return active[0]
    uncertain = [r for r in results if r.status == "UNDETERMINED"]
    if uncertain:
        return uncertain[0]
    applicable = [r for r in results if r.status != "NOT_APPLICABLE"]
    return applicable[0] if applicable else _result("NOT_APPLICABLE", "No study matches the target variant and disease.", guideline=ACMG)
