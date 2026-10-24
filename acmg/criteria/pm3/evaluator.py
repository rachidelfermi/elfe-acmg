"""Conservative PM3 scorer; no annotation or external lookup is performed."""

from .guidelines import POINTS, THRESHOLDS, select_guideline
from .models import Case, PM3Input, SecondVariant


def _coerce(raw: PM3Input | dict) -> PM3Input:
    if isinstance(raw, PM3Input):
        return raw
    data = dict(raw)
    data["cases"] = [
        case if isinstance(case, Case) else Case(**{
            **case,
            "second_variant": (
                SecondVariant(**case["second_variant"])
                if isinstance(case.get("second_variant"), dict)
                else case.get("second_variant")
            ),
        }) for case in data.get("cases", [])
    ]
    return PM3Input(**data)


def _phase(case: Case) -> tuple[str, str | None]:
    if case.phase == "CONFLICTING" or case.phase_confident is False:
        return "CONFLICTING", "Conflicting or unreliable phase evidence"
    parent_calls = (
        (case.maternal_candidate, case.maternal_second),
        (case.paternal_candidate, case.paternal_second),
    )
    inferred = None
    for candidate, second in parent_calls:
        if candidate is not None and second is not None:
            if {candidate, second} == {"HET", "REF"}:
                inferred = "TRANS"
    if inferred and case.phase in ("TRANS", "CIS") and case.phase != inferred:
        return "CONFLICTING", "Parent calls conflict with reported phase"
    if case.phase in ("TRANS", "CIS"):
        if not case.phase_method and not inferred:
            return "UNKNOWN", "Phase confirmation method missing"
        return case.phase, None
    if inferred:
        return inferred, None
    return "UNKNOWN", None


def _case(case: Case, data: PM3Input, guideline) -> dict:
    missing: list[str] = []
    exclusions: list[str] = []
    used: list[str] = []

    def required(value, field):
        used.append(field)
        if value is None:
            missing.append(field)
        elif value is False:
            exclusions.append(field)

    required(case.affected, "affected")
    required(case.clinical_disease_match, "clinical_disease_match")
    required(case.candidate_call_valid, "candidate_call_valid")
    if case.candidate_allele_balance_valid is False:
        exclusions.append("candidate_allele_balance_valid")
    if case.alternative_molecular_diagnosis is True:
        exclusions.append("alternative_molecular_diagnosis")
    if case.evidence_conflict:
        missing.append("evidence_conflict_resolution")
    if case.case_overlap:
        exclusions.append("overlapping_case_report")
    if guideline.gene == "GAA":
        required(case.diagnosis_confirmed, "diagnosis_confirmed")
        if case.pseudodeficiency_cis_context:
            required(data.candidate_pm2_supporting, "candidate_pm2_supporting")
            required(case.functional_support, "functional_support")
            required(case.clinical_laboratory_support, "clinical_laboratory_support")

    category = "UNRESOLVED"
    phase = "NOT_APPLICABLE"
    qualification = "NOT_APPLICABLE"
    points = 0.0
    if case.candidate_genotype == "HOM_ALT" or case.homozygous_observation:
        category = "HOMOZYGOUS"
        used.append("candidate_genotype")
        if case.candidate_genotype != "HOM_ALT" or case.homozygous_observation is False:
            exclusions.append("inconsistent_homozygous_genotype")
        required(case.homozygous_genotype_valid, "homozygous_genotype_valid")
        if case.deletion_mimicking_homozygosity or case.uniparental_disomy:
            exclusions.append("apparent_homozygosity")
        if case.copy_number_normal is False:
            exclusions.append("copy_number_status")
        # A validated true homozygous genotype is sufficient; a separate CN assay is optional.
        points = 0.5
    elif case.candidate_genotype == "HET":
        category = "COMPOUND_HET"
        second = case.second_variant
        if second is None:
            missing.append("second_variant")
        else:
            used.extend(("second_variant", "second_variant_classification", "phase"))
            if not second.variant_id:
                missing.append("second_variant_id")
            if second.variant_id == data.candidate_variant_id:
                exclusions.append("same_variant")
            if second.gene_id is None:
                missing.append("second_variant_gene_id")
            elif second.gene_id != data.gene_id:
                exclusions.append("second_variant_gene_mismatch")
            if second.genotype is None:
                missing.append("second_variant_genotype")
            elif second.genotype != "HET":
                exclusions.append("second_variant_not_heterozygous")
            required(second.call_valid, "second_variant_call_valid")
            required(second.rare_for_disease, "second_variant_rare_for_disease")
            required(second.disease_match, "second_variant_disease_match")
            if second.mechanism_compatible is False:
                exclusions.append("second_variant_mechanism")
            if second.allele_balance_valid is False:
                exclusions.append("second_variant_allele_balance_valid")
            if second.classification_conflict:
                missing.append("second_variant_classification_conflict")
            if second.classification is None:
                missing.append("second_variant_classification")
            elif second.classification in ("BENIGN", "LIKELY_BENIGN"):
                exclusions.append("nonqualifying_second_variant_classification")
            if second.classification_independent is None:
                missing.append("second_variant_classification_independence")
            elif second.classification_independent is False:
                exclusions.append("circular_second_variant_classification")
            required(second.classification_credible, "second_variant_classification_credible")
            if guideline.gene == "GAA" and second.classification in ("PATHOGENIC", "LIKELY_PATHOGENIC"):
                required(second.vcep_classification_confirmed, "vcep_second_variant_classification")
            if guideline.gene == "GAA" and second.classification == "VUS":
                exclusions.append("GAA_VCEP_VUS_not_in_point_table")
            phase, phase_error = _phase(case)
            if phase_error:
                missing.append(phase_error)
            if phase == "CIS":
                exclusions.append("confirmed_cis")
            qualification = second.classification or "UNKNOWN"
            points = POINTS.get((phase, qualification), 0.0)
    else:
        missing.append("candidate_genotype")

    if exclusions or missing:
        points = 0.0
    return {
        "case_id": case.case_id,
        "family_id": case.family_id,
        "eligibility": "EXCLUDED" if exclusions else ("UNDETERMINED" if missing else "ELIGIBLE"),
        "evidence_state": "CONFLICTING" if (case.evidence_conflict or phase == "CONFLICTING" or (case.second_variant and case.second_variant.classification_conflict)) else ("UNKNOWN" if missing else "KNOWN"),
        "observation_category": category,
        "phase_status": phase,
        "second_variant_qualification": qualification,
        "points_awarded": points,
        "missing_critical_features": sorted(set(missing)),
        "exclusions": sorted(set(exclusions)),
        "features_used": sorted(set(used)),
        "pp1_overlap_review": case.pp1_overlap,
        "guideline_source": guideline.source,
        "second_variant_id": case.second_variant.variant_id if case.second_variant else None,
    }


def evaluate_pm3(biological_features: PM3Input | dict) -> dict:
    data = _coerce(biological_features)
    guideline = select_guideline(data.specification, data.gene_id, data.disease_id)
    for case in data.cases:
        if case.candidate_genotype not in (None, "REF", "HET", "HOM_ALT", "HEMIZYGOUS", "OTHER"):
            raise ValueError("Unknown candidate genotype")
        if case.phase not in (None, "TRANS", "CIS", "UNKNOWN", "CONFLICTING"):
            raise ValueError("Unknown phase state")
        if case.second_variant:
            if case.second_variant.classification not in (None, "PATHOGENIC", "LIKELY_PATHOGENIC", "VUS", "LIKELY_BENIGN", "BENIGN"):
                raise ValueError("Unknown second-variant classification")
            if case.second_variant.genotype not in (None, "REF", "HET", "HOM_ALT", "HEMIZYGOUS", "OTHER"):
                raise ValueError("Unknown second-variant genotype")
    global_missing = []
    global_exclusions = []
    for name in ("candidate_variant_id", "gene_id", "disease_id", "inheritance_mode", "candidate_rare_for_disease"):
        if getattr(data, name) is None:
            global_missing.append(name)
    if data.inheritance_mode is not None and data.inheritance_mode != "AUTOSOMAL_RECESSIVE":
        global_exclusions.append("inheritance_mode")
    for name in ("gene_disease_valid", "disease_mechanism_compatible", "candidate_rare_for_disease"):
        if getattr(data, name) is False:
            global_exclusions.append(name)
    if not data.cases:
        global_missing.append("cases")

    observations = [_case(case, data, guideline) for case in data.cases]
    seen_cases = set()
    seen_families = set()
    homo_points = vus_points = 0.0
    pair_points: dict[str, float] = {}
    for case, obs in zip(data.cases, observations):
        if obs["eligibility"] != "ELIGIBLE":
            continue
        if not case.case_id:
            obs["missing_critical_features"].append("case_id")
            obs["eligibility"] = "UNDETERMINED"
            obs["points_awarded"] = 0.0
            continue
        if case.case_id in seen_cases:
            obs["eligibility"] = "EXCLUDED"
            obs["exclusions"].append("duplicate_case")
            obs["points_awarded"] = 0.0
            continue
        seen_cases.add(case.case_id)
        if len(data.cases) > 1:
            if not case.family_id or case.independence_confirmed is None:
                obs["eligibility"] = "UNDETERMINED"
                obs["missing_critical_features"].append("family_id_or_independence")
                obs["points_awarded"] = 0.0
                continue
            if case.independence_confirmed is False or case.family_id in seen_families:
                obs["eligibility"] = "EXCLUDED"
                obs["exclusions"].append("related_or_overlapping_case")
                obs["points_awarded"] = 0.0
                continue
            seen_families.add(case.family_id)
        if obs["observation_category"] == "HOMOZYGOUS":
            awarded = min(obs["points_awarded"], max(0.0, 1.0 - homo_points))
            homo_points += awarded
            obs["points_awarded"] = awarded
        elif obs["second_variant_qualification"] == "VUS":
            awarded = min(obs["points_awarded"], max(0.0, 0.5 - vus_points))
            vus_points += awarded
            obs["points_awarded"] = awarded
        if guideline.repeated_unphased_pair_cap is not None and obs["phase_status"] == "UNKNOWN":
            pair = obs["second_variant_id"]
            awarded = min(obs["points_awarded"], max(0.0, guideline.repeated_unphased_pair_cap - pair_points.get(pair, 0.0)))
            pair_points[pair] = pair_points.get(pair, 0.0) + awarded
            obs["points_awarded"] = awarded
    if global_missing or global_exclusions:
        for obs in observations:
            obs["points_awarded"] = 0.0
    total = round(sum(obs["points_awarded"] for obs in observations), 2)
    strength = next((label for threshold, label in THRESHOLDS if total >= threshold), "NONE")
    status = ("NOT_APPLICABLE" if global_exclusions else
              "UNDETERMINED" if global_missing else
              "ACTIVATED" if strength != "NONE" else
              "UNDETERMINED" if any(obs["eligibility"] == "UNDETERMINED" for obs in observations) else "NOT_MET")
    return {
        "criterion": "PM3", "status": status, "strength": strength if status == "ACTIVATED" else "NONE",
        "total_points": total, "eligible_observations": sum(o["points_awarded"] > 0 for o in observations),
        "observations": observations,
        "features_used": sorted({f for o in observations for f in o["features_used"]}),
        "missing_critical_features": sorted(set(global_missing + [m for o in observations if o["eligibility"] == "UNDETERMINED" for m in o["missing_critical_features"]])),
        "exclusions_evaluated": sorted(set(global_exclusions + [e for o in observations for e in o["exclusions"]])),
        "guideline": {"name": guideline.name, "version": guideline.version, "rule_id": guideline.rule_id, "source": guideline.source},
        "reason": f"{status}: {total:g} verified PM3 points; {strength if status == 'ACTIVATED' else 'no applicable strength'}.",
    }
