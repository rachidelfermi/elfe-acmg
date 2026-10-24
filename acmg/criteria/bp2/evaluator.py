"""Deterministic BP2 evaluation from curated evidence."""

from .guidelines import ACMG, applicable_specification
from .models import BP2Result, BiologicalFeatures, Classification, DecisionStatus, Inheritance, Phase, State


def _result(f, spec, status, reason, path=None, observations=(), used=(), missing=(), exclusions=()):
    first = observations[0] if observations else None
    return BP2Result("BP2", status, "SUPPORTING" if status == DecisionStatus.ACTIVATED else "NONE",
                     path, first.comparison.variant_id if first else None,
                     first.phase.value if first else None, tuple(used), tuple(dict.fromkeys(missing)),
                     tuple(dict.fromkeys(exclusions)), spec.ref, reason,
                     tuple(o.case_id or f"observation:{f.observations.index(o)+1}" for o in observations),
                     tuple(sorted(set(f.other_criteria_from_same_observation) &
                                  {"PM3", "PP1", "PS2", "PM6", "BP5"})))


def _required(flag, label, missing, failures, value=True):
    if flag.state == State.CONFLICTING:
        missing.append(label + " (conflicting)")
    elif flag.state != State.KNOWN:
        missing.append(label)
    elif flag.value != value:
        failures.append(label)


def _dedupe(observations):
    seen = set()
    result = []
    for index, obs in enumerate(observations):
        key = (obs.case_id, obs.comparison.variant_id) if obs.case_id else f"unidentified:{index}"
        if key not in seen:
            seen.add(key)
            result.append(obs)
    return result


def _dedupe_cases(observations):
    seen = set()
    result = []
    for index, obs in enumerate(observations):
        key = obs.case_id or f"unidentified:{index}"
        if key not in seen:
            seen.add(key)
            result.append(obs)
    return result


def _evaluated_exclusions(features):
    return tuple(label for flag, label in (
        (features.synergy_or_complex_allele, "synergy or complex allele"),
        (features.haplotype_specific_pathogenicity, "haplotype-specific pathogenicity"),
        (features.candidate_independent_disease_effect, "independent candidate disease effect"),
        (features.dual_diagnosis, "dual diagnosis"),
    ) if flag.is_no())


def evaluate_bp2(features: BiologicalFeatures) -> BP2Result:
    """Apply base ACMG or built-in versioned VCEP BP2 rules, offline."""
    if not isinstance(features, BiologicalFeatures):
        raise TypeError("evaluate_bp2 expects BiologicalFeatures")
    spec = applicable_specification(features.gene_id, features.disease_id)
    if spec is ACMG and features.gene_id in {"DICER1", "BMPR2"}:
        return _result(features, spec, DecisionStatus.UNDETERMINED,
                       "Gene has a built-in VCEP rule, but the supplied disease ID is outside its verified scope.",
                       missing=("verified gene-disease specification scope",))
    if features.guideline_version and features.guideline_version != spec.ref.version:
        return _result(features, spec, DecisionStatus.UNDETERMINED,
                       "Supplied guideline version differs from the implemented version.",
                       missing=("matching guideline version",))
    if not features.observations:
        return _result(features, spec, DecisionStatus.UNDETERMINED,
                       "No comparison observation supplied.", missing=("comparison observation",))
    if features.gene_disease_valid.is_no():
        return _result(features, spec, DecisionStatus.NOT_APPLICABLE,
                       "The supplied gene-disease relationship is not valid for this disease.",
                       exclusions=("gene-disease relationship",))
    reported_phases = {}
    for obs in features.observations:
        if obs.case_id:
            key = (obs.case_id, obs.comparison.variant_id)
            if key in reported_phases and reported_phases[key] != obs.phase:
                return _result(features, spec, DecisionStatus.UNDETERMINED,
                               "Contradictory phase reports for the same case and comparison.",
                               missing=("reconciled phase reports",))
            reported_phases[key] = obs.phase
    missing, failures = [], []
    _required(features.gene_disease_valid, "gene-disease validity", missing, failures)
    if features.inheritance == Inheritance.UNKNOWN:
        missing.append("inheritance mode")
    if features.disease_mechanism_compatible.state == State.CONFLICTING:
        missing.append("disease mechanism compatibility (conflicting)")
    elif features.disease_mechanism_compatible.is_no():
        failures.append("disease mechanism incompatible")
    for flag, label in ((features.synergy_or_complex_allele, "synergy or complex allele"),
                        (features.haplotype_specific_pathogenicity, "haplotype-specific pathogenicity"),
                        (features.candidate_independent_disease_effect, "independent candidate disease effect"),
                        (features.dual_diagnosis, "dual diagnosis"),
                        (features.evidence_conflict, "evidence conflict")):
        if flag.state == State.CONFLICTING:
            missing.append(label + " (conflicting)")
        elif flag.is_yes():
            if label == "evidence conflict":
                missing.append(label)
            else:
                failures.append(label)
    if failures:
        return _result(features, spec, DecisionStatus.NOT_MET,
                       "BP2 excluded: " + ", ".join(failures) + ".", exclusions=failures)
    if missing:
        return _result(features, spec, DecisionStatus.UNDETERMINED,
                       "Essential context is unresolved.", missing=missing)

    cis, trans, unknown = [], [], []
    all_missing, all_failures = [], []
    for obs in _dedupe(features.observations):
        m, bad = [], []
        c = obs.comparison
        if c.variant_id == features.candidate_variant_id:
            bad.append("comparison is candidate")
        if c.gene_id is None:
            m.append("comparison gene")
        elif c.gene_id != features.gene_id:
            bad.append("comparison gene mismatch")
        if c.classification == Classification.UNKNOWN:
            m.append("comparison classification")
        elif c.classification != Classification.PATHOGENIC and not (
            spec.comparison_lp_allowed and c.classification == Classification.LIKELY_PATHOGENIC):
            bad.append("comparison classification not qualifying")
        if not c.classification_source:
            m.append("comparison classification source")
        if not c.classification_evidence:
            m.append("comparison classification evidence")
        for flag, label in ((c.classification_independent, "independent classification"),
                            (c.observed_independently, "independent observation of comparison"),
                            (c.disease_match, "comparison disease match"),
                            (c.mechanism_compatible, "comparison mechanism compatibility")):
            _required(flag, label, m, bad)
        if spec.requires_vcep_classification:
            _required(c.vcep_classified, "applicable VCEP classification", m, bad)
        for flag, label in ((features.candidate_call_reliable, "candidate call"),
                            (features.candidate_allele_balance_reliable, "candidate allele balance"),
                            (c.call_reliable, "comparison call"),
                            (c.allele_balance_reliable, "comparison allele balance")):
            if flag.state == State.CONFLICTING:
                m.append(label + " (conflicting)")
            elif flag.is_no():
                bad.append(label + " unreliable")
        if obs.phase_conflict.state == State.CONFLICTING or obs.phase_conflict.is_yes():
            m.append("phase conflict")
        for genotype, label in ((obs.candidate_genotype, "candidate genotype"),
                                (c.genotype, "comparison genotype")):
            if genotype and genotype.upper() not in ("HET", "HETEROZYGOUS"):
                bad.append(label + " incompatible")
        if obs.phase == Phase.CONFLICTING:
            m.append("phase (conflicting)")
        elif obs.phase == Phase.UNKNOWN:
            if spec.unknown_phase_min_cases is None:
                m.append("confirmed phase")
            elif not (obs.candidate_genotype and c.genotype):
                m.append("both genotypes for unphased co-occurrence")
        else:
            if not obs.phase_method:
                m.append("phase confirmation method")
            _required(obs.phase_reliable, "reliable confirmed phase", m, bad)
        if obs.phase not in (Phase.UNKNOWN, Phase.CONFLICTING) and obs.phase.value not in spec.allowed_phases:
            bad.append("phase excluded by specification")
        if obs.phase in (Phase.CIS, Phase.UNKNOWN):
            _required(features.synergy_or_complex_allele, "complex-allele interaction excluded", m, bad,
                      value=False)
        if obs.phase == Phase.TRANS:
            if features.inheritance != Inheritance.AD:
                bad.append("trans requires dominant disease")
            if spec.trans_requires_full_penetrance:
                _required(features.full_penetrance, "full penetrance", m, bad)
            _required(obs.comparison_explains_disease, "comparison explains disease", m, bad)
            _required(features.candidate_independent_disease_effect,
                      "independent candidate disease effect excluded", m, bad, value=False)
            _required(features.dual_diagnosis, "dual diagnosis excluded", m, bad, value=False)
            if spec.trans_requires_phenotype_check:
                _required(obs.phenotype_not_more_severe, "phenotype not more severe", m, bad)
        if spec.requires_affected:
            _required(obs.affected, "affected individual", m, bad)
            _required(obs.diagnosis_matches, "matching diagnosis", m, bad)
        if spec.cis_min_cases > 1 and obs.phase in (Phase.CIS, Phase.UNKNOWN) and not obs.case_id:
            m.append("case ID for independent count")
        if bad:
            all_failures.extend(bad)
        elif m:
            all_missing.extend(m)
        elif obs.phase == Phase.CIS:
            cis.append(obs)
        elif obs.phase == Phase.TRANS:
            trans.append(obs)
        elif obs.phase == Phase.UNKNOWN:
            unknown.append(obs)
    if len(trans) >= spec.trans_min_cases:
        return _result(features, spec, DecisionStatus.ACTIVATED,
                       "Verified trans observation satisfies the applicable BP2 rule.",
                       "CONFIRMED_TRANS_DOMINANT", trans[:spec.trans_min_cases],
                       ("comparison classification", "comparison independence", "confirmed trans phase",
                        "dominant inheritance", "penetrance or VCEP exception"),
                       exclusions=_evaluated_exclusions(features))
    all_allelic = cis + unknown
    allelic = _dedupe_cases(all_allelic)
    if len(allelic) >= spec.cis_min_cases and len({o.comparison.variant_id for o in all_allelic}) >= spec.cis_min_unique_comparisons:
        return _result(features, spec, DecisionStatus.ACTIVATED,
                       "Verified allelic observations satisfy the applicable BP2 rule.",
                       "VCEP_MULTIPLE_ALLELIC_OBSERVATIONS" if spec.cis_min_cases > 1 else "CONFIRMED_CIS_WITH_PATHOGENIC",
                       allelic[:spec.cis_min_cases],
                       ("comparison classification", "comparison independence", "allelic observations", "gene-disease context"),
                       exclusions=_evaluated_exclusions(features))
    if allelic and spec.cis_min_cases > 1:
        all_missing.append(f"{spec.cis_min_cases} independent cases with {spec.cis_min_unique_comparisons} unique comparisons")
    if all_missing:
        return _result(features, spec, DecisionStatus.UNDETERMINED,
                       "No complete pathway; essential evidence is missing or conflicting.",
                       missing=all_missing, exclusions=all_failures)
    return _result(features, spec, DecisionStatus.NOT_MET,
                   "Supplied observations do not satisfy an applicable BP2 pathway.",
                   exclusions=all_failures)
