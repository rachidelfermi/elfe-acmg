"""Evaluate already-curated de novo observations for PM6."""

from __future__ import annotations

from dataclasses import fields
from typing import Any, Mapping

from .guidelines import POINTS, select_applicable_guideline, strength_for
from .models import Observation, ObservationResult, PM6Input, PM6Result, Phenotype, Status, Strength


def _value(value: Any) -> Any:
    if isinstance(value, Mapping) and "state" in value:
        state = str(value["state"]).upper()
        if state == "KNOWN":
            return value.get("value")
        if state in {"UNKNOWN", "NOT_APPLICABLE"}:
            return None
        if state == "CONFLICTING":
            raise ValueError("Conflicting evidence")
        raise ValueError(f"Invalid evidence state: {state}")
    return value


def _boolean(value: Any, name: str) -> bool | None:
    value = _value(value)
    if value is None:
        return None
    if type(value) is bool:
        return value
    if name in {"maternity_confirmed", "paternity_confirmed"} and isinstance(value, str):
        if value.upper() == "CONFIRMED":
            return True
        if value.upper() in {"UNCONFIRMED", "UNKNOWN", "NOT_TESTED"}:
            return None
        if value.upper() in {"DISPROVED", "CONFLICTING"}:
            return False
    raise ValueError(f"{name} must be true, false, or an explicit evidence state")


def _absence(raw: Mapping[str, Any], parent: str) -> bool | None:
    absence = _value(raw.get(f"{parent}_variant_absence"))
    genotype = _value(raw.get(f"{parent}_genotype"))
    if absence == "CONFLICTING" or genotype == "CONFLICTING":
        raise ValueError(f"Conflicting {parent} variant evidence")
    if genotype is not None:
        if not isinstance(genotype, str) or genotype.upper() not in {"ABSENT", "PRESENT", "REFERENCE", "HET", "HOM_ALT", "MOSAIC"}:
            raise ValueError(f"Unsupported {parent}_genotype")
        genotype_absent = genotype.upper() in {"ABSENT", "REFERENCE"}
        if absence is not None and absence is not genotype_absent:
            raise ValueError(f"Conflicting {parent} genotype and absence")
        return genotype_absent
    return _boolean(absence, f"{parent}_variant_absence")


def _adequacy(raw: Mapping[str, Any], parent: str) -> bool | None:
    depth = _value(raw.get(f"{parent}_read_depth"))
    if depth is not None and (not isinstance(depth, (int, float)) or depth < 0):
        raise ValueError(f"Invalid {parent}_read_depth")
    if depth == 0:
        return False
    quality = _value(raw.get(f"{parent}_variant_call_quality"))
    if quality is not None:
        if str(quality).upper() in {"POOR", "LOW", "INADEQUATE"}:
            return False
    explicit = raw.get(f"{parent}_test_adequate", raw.get("parental_testing_reliability"))
    if explicit is not None:
        if isinstance(explicit, str):
            if explicit.upper() in {"ADEQUATE", "HIGH", "RELIABLE"}:
                return True
            if explicit.upper() in {"POOR", "LOW", "INADEQUATE"}:
                return False
        return _boolean(explicit, f"{parent}_test_adequate")
    if quality is not None:
        if str(quality).upper() in {"ADEQUATE", "HIGH", "RELIABLE"}:
            return True
    return None


def _phenotype(raw: Mapping[str, Any]) -> Phenotype | None:
    value = _value(raw.get("phenotype", raw.get("phenotype_category")))
    if value is None:
        if _value(raw.get("phenotype_consistency")) is False:
            return Phenotype.INCONSISTENT
        if _value(raw.get("phenotype_consistency")) is True:
            if _value(raw.get("phenotype_specificity")) in {"HIGH", "HIGHLY_SPECIFIC"}:
                return Phenotype.HIGHLY_SPECIFIC
            if _value(raw.get("phenotype_heterogeneity")) in {"HIGH", True}:
                return Phenotype.HETEROGENEOUS
            if _value(raw.get("phenotype_specificity")) in {"CONSISTENT", "NOT_HIGH"}:
                return Phenotype.CONSISTENT
        return None
    try:
        return Phenotype(str(value).upper())
    except ValueError as error:
        raise ValueError(f"Invalid phenotype category: {value}") from error


def _observation(raw: Mapping[str, Any]) -> Observation:
    if not isinstance(raw, Mapping):
        raise TypeError("Each observation must be a mapping")
    names = {item.name for item in fields(Observation)}
    data = {key: _value(value) for key, value in raw.items() if key in names}
    for parent in ("maternal", "paternal"):
        data[f"{parent}_variant_absent"] = _absence(raw, parent) if f"{parent}_variant_absence" in raw or f"{parent}_genotype" in raw else _boolean(raw.get(f"{parent}_variant_absent"), f"{parent}_variant_absent")
        data[f"{parent}_test_adequate"] = _adequacy(raw, parent)
    data["phenotype"] = _phenotype(raw)
    aliases = {"maternity_confirmed": "maternity_confirmation_status", "paternity_confirmed": "paternity_confirmation_status",
               "affected": "affected_status", "sample_identity_valid": "sample_identity_confirmation",
               "independent": "case_independence", "parental_mosaicism": "parental_low_level_mosaicism",
               "technical_artifact": "technical_artifact_evidence", "conflicting_segregation": "conflicting_segregation_evidence",
               "proband_call_adequate": "proband_variant_call_quality"}
    for target, source in aliases.items():
        if target not in raw and source in raw:
            value = _value(raw[source])
            if target == "proband_call_adequate" and isinstance(value, str):
                value = value.upper() in {"ADEQUATE", "HIGH", "RELIABLE"}
            data[target] = value
    if "proband_variant_present" in raw:
        data["proband_variant_present"] = _boolean(raw["proband_variant_present"], "proband_variant_present")
    if "study_ids" in data and data["study_ids"] is not None:
        data["study_ids"] = tuple(data["study_ids"])
    for key in ("maternity_confirmed", "paternity_confirmed", "sample_identity_valid", "affected", "independent", "parental_mosaicism", "proband_mosaicism", "technical_artifact", "mapping_ambiguity", "contamination", "conflicting_segregation", "alternative_diagnosis", "proband_call_adequate", "family_history_compatible"):
        if key in data:
            data[key] = _boolean(data[key], key)
    return Observation(**data)


def validate_features(raw: PM6Input | Mapping[str, Any]) -> PM6Input:
    if isinstance(raw, PM6Input):
        return raw
    if not isinstance(raw, Mapping):
        raise TypeError("biological_features must be a mapping or PM6Input")
    cases = raw.get("observations", raw.get("cases"))
    if cases is None:
        cases = [raw]
    if not isinstance(cases, (list, tuple)) or not cases:
        raise ValueError("observations must contain at least one case")
    shared = {key: value for key, value in raw.items() if key not in {"observations", "cases"}}
    observations = tuple(_observation({**shared, **case}) for case in cases)
    return PM6Input(variant_id=_value(raw.get("variant_id")), gene_id=_value(raw.get("gene_id")),
                    disease_id=_value(raw.get("disease_id")), gene_disease_valid=_boolean(raw.get("gene_disease_valid", raw.get("gene_disease_validity")), "gene_disease_valid"),
                    inheritance_mode=_value(raw.get("inheritance_mode", raw.get("disease_inheritance_mode"))),
                    mechanism_compatible=_boolean(raw.get("mechanism_compatible"), "mechanism_compatible"),
                    second_pathogenic_allele=_boolean(raw.get("second_pathogenic_allele"), "second_pathogenic_allele"),
                    vcep_profile=_value(raw.get("vcep_profile")), observations=observations)


def evaluate_de_novo_case(case: Observation, guideline: dict[str, str]) -> ObservationResult:
    missing: list[str] = []
    exclusions: list[str] = []
    if case.proband_variant_present is False:
        exclusions.append("proband_variant_absent")
    if case.maternal_variant_absent is False or case.paternal_variant_absent is False:
        exclusions.append("parent_carries_variant")
    if case.phenotype is Phenotype.INCONSISTENT or case.affected is False:
        exclusions.append("phenotype_or_affected_status_incompatible")
    if case.family_history_compatible is False:
        exclusions.append("family_history_incompatible")
    for name in ("sample_identity_valid", "proband_call_adequate", "maternal_test_adequate", "paternal_test_adequate"):
        if getattr(case, name) is False:
            exclusions.append(f"{name}_false")
    for name in ("parental_mosaicism", "technical_artifact", "mapping_ambiguity", "contamination", "conflicting_segregation", "alternative_diagnosis", "evidence_conflict"):
        if getattr(case, name):
            exclusions.append(name)
    if case.duplicate_of or case.related_to_other_case or case.independent is False:
        exclusions.append("non_independent_case")
    if case.already_counted_ps2:
        exclusions.append("already_counted_under_PS2")
    if case.apparent_germline_mosaicism and not (case.maternity_confirmed and case.paternity_confirmed):
        exclusions.append("germline_mosaicism_requires_confirmed_parentage")
    if exclusions:
        return ObservationResult(case.case_id, Status.NOT_MET, None, case.phenotype.value if case.phenotype else None, 0.0, exclusions=tuple(exclusions))
    for name in ("proband_variant_present", "maternal_variant_absent", "paternal_variant_absent", "affected", "phenotype"):
        if getattr(case, name) is None:
            missing.append(name)
    # A documented adequate negative test may be represented by an absence assertion;
    # explicit poor quality overrides it, while raw metrics are never mandatory.
    if case.maternal_variant_absent and case.maternal_test_adequate is None:
        missing.append("maternal_test_adequate")
    if case.paternal_variant_absent and case.paternal_test_adequate is None:
        missing.append("paternal_test_adequate")
    if missing:
        return ObservationResult(case.case_id, Status.UNDETERMINED, None, case.phenotype.value if case.phenotype else None, 0.0, tuple(missing))
    category = "CONFIRMED" if case.maternity_confirmed is True and case.paternity_confirmed is True else "ASSUMED"
    # An explicit relationship contradiction is different from lack of confirmation.
    if case.maternity_confirmed is False or case.paternity_confirmed is False:
        return ObservationResult(case.case_id, Status.NOT_MET, None, case.phenotype.value, 0.0, exclusions=("parentage_disproved",))
    if guideline["rule_id"] == "HNF1A-v3.1-PS2":
        facts = case.vcep_facts
        if facts.get("affected_parent_meets_PP4") or facts.get("affected_parent_similar_atypical_diabetes") or facts.get("affected_parent_nonautoimmune_diabetes_before_30"):
            return ObservationResult(case.case_id, Status.NOT_MET, None, case.phenotype.value, 0.0, exclusions=("HNF1A_affected_parent_exclusion",), rule_id="HNF1A-v3.1-PS2")
        if facts.get("autoimmune_diabetes") or facts.get("absolute_insulin_deficiency"):
            return ObservationResult(case.case_id, Status.NOT_MET, None, case.phenotype.value, 0.0, exclusions=("HNF1A_phenotype_exclusion",), rule_id="HNF1A-v3.1-PS2")
        if case.phenotype is Phenotype.HIGHLY_SPECIFIC and not (facts.get("mody_probability_at_least_50_percent") and facts.get("hnf4a_negative")):
            return ObservationResult(case.case_id, Status.UNDETERMINED, None, case.phenotype.value, 0.0, missing_critical_features=("HNF1A_PP4_criteria",), rule_id="HNF1A-v3.1-PS2")
    if guideline["rule_id"] == "SLC9A6-v6.0-PM6":
        consistent = case.vcep_facts.get("neurodevelopmental_phenotype_consistent")
        if consistent is None:
            return ObservationResult(case.case_id, Status.UNDETERMINED, None, case.phenotype.value, 0.0,
                                     missing_critical_features=("SLC9A6_neurodevelopmental_phenotype",), rule_id="SLC9A6-v6.0-PM6")
        if consistent is False:
            return ObservationResult(case.case_id, Status.NOT_MET, None, case.phenotype.value, 0.0,
                                     exclusions=("SLC9A6_neurodevelopmental_phenotype_incompatible",), rule_id="SLC9A6-v6.0-PM6")
        if case.proband_mosaicism is True:
            category = "CONFIRMED"
        return ObservationResult(case.case_id, Status.ACTIVATED, category, case.phenotype.value,
                                 2.0 if category == "CONFIRMED" else 1.0, rule_id="SLC9A6-v6.0-PM6")
    return ObservationResult(case.case_id, Status.ACTIVATED, category, case.phenotype.value,
                             POINTS[category][case.phenotype], rule_id=guideline["rule_id"] if guideline["rule_id"] != "SVI-v1.1-tables-1-2" else "SVI-v1.1-table-1")


def evaluate_pm6(biological_features: PM6Input | Mapping[str, Any]) -> PM6Result:
    try:
        evidence = validate_features(biological_features)
    except ValueError as error:
        if "conflict" not in str(error).lower():
            raise
        return PM6Result(Status.UNDETERMINED, Strength.NONE, 0.0, 0, (), (), ("resolve_conflicting_evidence",),
                         ("conflicting_evidence",),
                         {"name": "ClinGen SVI de novo", "version": "1.1", "rule_id": "SVI-v1.1-tables-1-2"},
                         str(error))
    try:
        guideline = select_applicable_guideline(evidence.gene_id, evidence.disease_id, evidence.vcep_profile)
    except ValueError as error:
        return PM6Result(Status.UNDETERMINED, Strength.NONE, 0.0, 0, (), (), ("verified_applicable_vcep_profile",), (),
                         {"name": "unresolved VCEP", "version": "unknown", "rule_id": "unresolved"}, str(error))
    seen_cases: set[str] = set()
    seen_families: set[str] = set()
    results: list[ObservationResult] = []
    for case in evidence.observations:
        if len(evidence.observations) > 1 and case.independent is not True and not case.family_id:
            results.append(ObservationResult(case.case_id, Status.UNDETERMINED, None,
                                             case.phenotype.value if case.phenotype else None, 0.0,
                                             missing_critical_features=("case_independence",)))
            continue
        if (case.case_id and case.case_id in seen_cases) or (case.family_id and case.family_id in seen_families):
            results.append(ObservationResult(case.case_id, Status.NOT_MET, None, case.phenotype.value if case.phenotype else None, 0.0, exclusions=("duplicate_or_related_case",)))
            continue
        result = evaluate_de_novo_case(case, guideline)
        results.append(result)
        if result.status is Status.ACTIVATED:
            if case.case_id:
                seen_cases.add(case.case_id)
            if case.family_id:
                seen_families.add(case.family_id)
    eligible = [r for r in results if r.status is Status.ACTIVATED and r.points > 0]
    missing = sorted({x for r in results for x in r.missing_critical_features})
    exclusions = sorted({x for r in results for x in r.exclusions})
    case_fields = ("proband_variant_present", "maternal_variant_absent", "paternal_variant_absent",
                   "maternal_test_adequate", "paternal_test_adequate", "sample_identity_valid",
                   "proband_call_adequate", "affected", "phenotype", "maternity_confirmed",
                   "paternity_confirmed", "parental_mosaicism", "proband_mosaicism",
                   "technical_artifact", "mapping_ambiguity", "contamination",
                   "conflicting_segregation", "alternative_diagnosis", "family_history_compatible",
                   "family_id", "independent", "vcep_facts")
    features_used = tuple(name for name in ("gene_id", "disease_id", "gene_disease_valid", "inheritance_mode", "mechanism_compatible", "second_pathogenic_allele")
                          if getattr(evidence, name) is not None) + tuple(
                              name for name in case_fields
                              if any(getattr(case, name) not in (None, {}) for case in evidence.observations))
    hetero = sum(r.points for r in eligible if r.phenotype == Phenotype.HETEROGENEOUS.value)
    total = round(sum(r.points for r in eligible) - max(0.0, hetero - 1.0), 2)
    confirmed = any(r.category == "CONFIRMED" for r in eligible)
    strength = strength_for(total)
    if evidence.inheritance_mode == "AUTOSOMAL_RECESSIVE" and evidence.second_pathogenic_allele is False and strength is not Strength.NONE:
        order = (Strength.NONE, Strength.SUPPORTING, Strength.MODERATE, Strength.STRONG, Strength.VERY_STRONG)
        strength = order[order.index(strength) - 1]
    if evidence.gene_disease_valid is False or evidence.mechanism_compatible is False:
        return PM6Result(Status.NOT_APPLICABLE, Strength.NONE, 0.0, 0, tuple(results), features_used, (), tuple(exclusions), guideline,
                         "Gene–disease relationship or molecular mechanism is incompatible.")
    if evidence.gene_disease_valid is None:
        needed = ("gene_disease_valid",)
        return PM6Result(Status.UNDETERMINED, Strength.NONE, total, len(eligible), tuple(results), features_used, tuple(sorted(set(missing + list(needed)))), tuple(exclusions), guideline,
                         "Essential gene–disease or mechanism evidence is missing.")
    if evidence.inheritance_mode == "AUTOSOMAL_RECESSIVE" and evidence.second_pathogenic_allele is None and total >= 0.5:
        return PM6Result(Status.UNDETERMINED, Strength.NONE, total, len(eligible), tuple(results), features_used,
                         tuple(sorted(set(missing + ["second_pathogenic_allele"]))), tuple(exclusions), guideline,
                         "Second pathogenic allele status is needed to apply the published AR strength adjustment.")
    if confirmed or guideline["rule_id"] == "HNF1A-v3.1-PS2":
        return PM6Result(Status.NOT_APPLICABLE, Strength.NONE, total, len(eligible), tuple(results), features_used, tuple(missing), tuple(exclusions), guideline,
                         "Pooled de novo evidence belongs under PS2; evaluate it once there.", "PS2")
    if strength is not Strength.NONE:
        return PM6Result(Status.ACTIVATED, strength, total, len(eligible), tuple(results), features_used, tuple(missing), tuple(exclusions), guideline,
                         f"{len(eligible)} independent assumed de novo observation(s) contribute {total:g} capped SVI points; PM6 {strength.value.lower()} threshold met.")
    if missing or any(r.status is Status.UNDETERMINED for r in results):
        return PM6Result(Status.UNDETERMINED, Strength.NONE, total, len(eligible), tuple(results), features_used, tuple(missing), tuple(exclusions), guideline,
                         "Available evidence does not reach a PM6 threshold; critical case evidence remains unresolved.")
    return PM6Result(Status.NOT_MET, Strength.NONE, total, len(eligible), tuple(results), features_used, (), tuple(exclusions), guideline,
                     "No eligible assumed de novo evidence reaches the 0.5-point supporting threshold.")
