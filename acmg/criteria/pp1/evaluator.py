"""Conservative PP1 evaluation with explicit minimum evidence per pathway."""

from dataclasses import fields
from math import isfinite, log10
from typing import Any

from .guidelines import Guideline, select_guideline
from .models import Family, Individual, LodEvidence, PP1Input, PP1Result

RANK = {"NONE": 0, "SUPPORTING": 1, "MODERATE": 2, "STRONG": 3}


def _convert(cls: type, value: Any):
    if isinstance(value, cls):
        return value
    if not isinstance(value, dict):
        raise TypeError(f"{cls.__name__} must be a mapping or {cls.__name__}")
    allowed = {f.name for f in fields(cls)}
    unknown = set(value) - allowed
    if unknown:
        raise ValueError(f"Unsupported {cls.__name__} fields: {', '.join(sorted(unknown))}")
    data = dict(value)
    if cls is PP1Input:
        data["families"] = [_convert(Family, f) for f in data.get("families", [])]
    elif cls is Family:
        data["individuals"] = [_convert(Individual, i) for i in data.get("individuals", [])]
        if data.get("lod") is not None:
            data["lod"] = _convert(LodEvidence, data["lod"])
    return cls(**data)


def _result(status: str, reason: str, rule: Guideline | None, **kwargs) -> PP1Result:
    guideline = None if rule is None else {"name": rule.name, "version": rule.version, "rule_id": rule.rule_id, "source": rule.source}
    return PP1Result(status=status, reason=reason, guideline=guideline, **kwargs)


def _strength(value: float, thresholds: tuple[tuple[float, str], ...]) -> str:
    return next((s for threshold, s in thresholds if value >= threshold), "NONE")


def _count(f: Family, mode: str) -> tuple[int | None, int | None, str | None]:
    """Use either validated aggregate counts or individually curated independent events."""
    if f.affected_segregations is not None:
        if f.affected_segregations < 0 or (f.unaffected_recessive_segregations or 0) < 0:
            return None, None, "negative segregation count"
        if not all((f.count_validated, f.proband_excluded, f.relationships_confirmed, f.genotype_reliable, f.phenotype_relevant)):
            return None, None, "validated counts, proband exclusion, relationships, genotypes and phenotype"
        if mode == "AR" and f.biallelic_phase_confirmed is not True:
            return None, None, "biallelic trans phase"
        return f.affected_segregations, f.unaffected_recessive_segregations or 0, None
    if not f.individuals:
        return None, None, "family segregation observations"
    seen: set[str] = set()
    affected = unaffected = 0
    for person in f.individuals:
        if person.individual_id in seen:
            return None, None, "unique individual identifiers"
        seen.add(person.individual_id)
        if person.is_proband or person.informative is False:
            continue
        if person.informative is not True:
            continue
        if person.relationship_confirmed is not True or person.genotype_reliable is not True:
            return None, None, "confirmed relationship and reliable genotype"
        if person.affected is None or person.carries_variant is None:
            return None, None, "affected status and genotype"
        if person.affected and person.carries_variant:
            if person.phenotype_consistent is not True:
                return None, None, "consistent affected phenotype"
            if mode == "AR" and person.biallelic_in_trans is not True:
                return None, None, "biallelic trans genotype"
            affected += 1
        elif mode == "AR" and not person.affected and person.unaffected_recessive_at_risk is True:
            if person.biallelic_in_trans is False:
                unaffected += 1
            else:
                return None, None, "unaffected relative genotype"
        elif person.affected and not person.carries_variant:
            return None, None, "affected noncarrier explanation"
    return affected, unaffected, None


def _lod_valid(lod: LodEvidence) -> bool:
    return (isfinite(lod.score) and bool(lod.method) and lod.disease_model_validated is True
            and lod.penetrance_assumptions_validated is True and lod.pedigree_validated is True
            and lod.phenotype_model_validated is True)


def evaluate_pp1(biological_features: PP1Input | dict[str, Any]) -> PP1Result:
    data = _convert(PP1Input, biological_features)
    rule = select_guideline(data.gene_id, data.disease_id, data.guideline_id)
    if rule is None:
        return _result("UNDETERMINED", "No implemented verified PP1 specification matches gene, disease and requested version.", None,
                       missing_critical_features=["applicable verified PP1 specification"])
    invalid_states = {k: v for k, v in data.evidence_states.items()
                      if v not in {"KNOWN", "UNKNOWN", "CONFLICTING", "NOT_APPLICABLE"}}
    if invalid_states:
        raise ValueError(f"Invalid evidence states: {invalid_states}")
    critical_states = {k: v for k, v in data.evidence_states.items()
                       if k in {"gene_disease", "variant", "inheritance", "segregation"}}
    if any(v == "NOT_APPLICABLE" for v in critical_states.values()):
        return _result("NOT_APPLICABLE", "Essential PP1 evidence is explicitly not applicable.", rule,
                       exclusions_evaluated=[k for k, v in critical_states.items() if v == "NOT_APPLICABLE"])
    if any(v == "UNKNOWN" for v in critical_states.values()):
        return _result("UNDETERMINED", "Essential PP1 evidence is explicitly unknown.", rule,
                       missing_critical_features=[k for k, v in critical_states.items() if v == "UNKNOWN"])
    if data.evidence_states and any(v == "CONFLICTING" for v in data.evidence_states.values()):
        return _result("UNDETERMINED", "Critical evidence has conflicting states.", rule,
                       missing_critical_features=[k for k, v in data.evidence_states.items() if v == "CONFLICTING"])
    if data.gene_disease_validity == "DISPUTED" or data.disease_mechanism_compatible is False:
        return _result("NOT_APPLICABLE", "Gene–disease relationship or mechanism does not support this PP1 scope.", rule)
    required = []
    if data.gene_disease_validity not in {"DEFINITIVE", "STRONG", "ESTABLISHED"}:
        required.append("established gene_disease_validity")
    if data.inheritance_mode is None:
        required.append("inheritance_mode")
    elif data.inheritance_mode not in rule.modes:
        return _result("NOT_APPLICABLE", "Inheritance mode is outside the specification scope.", rule)
    if data.variant_reliable is not True:
        required.append("reliable variant identity")
    if data.disease_mechanism_compatible is not True:
        required.append("disease mechanism compatibility")
    if required:
        return _result("UNDETERMINED", "Essential gene, variant or inheritance evidence is missing.", rule,
                       missing_critical_features=required)
    if data.alternative_causal_variant is True or data.alternative_molecular_diagnosis is True:
        return _result("NOT_MET", "An established alternative cause explains the family phenotype.", rule,
                       exclusions_evaluated=["alternative cause"])
    if data.haplotype_confounding is True or data.evidence_overlap_prohibited is True:
        return _result("UNDETERMINED", "Linked variant or prohibited use of the same family evidence is unresolved.", rule,
                       exclusions_evaluated=["haplotype confounding", "evidence overlap"])
    if not data.families:
        return _result("UNDETERMINED", "No genotyped family segregation evidence was supplied.", rule,
                       missing_critical_features=["family segregation evidence"])

    ids = [f.family_id for f in data.families]
    if len(set(ids)) != len(ids) or any(not x for x in ids):
        return _result("UNDETERMINED", "Family identifiers are duplicated or missing.", rule,
                       missing_critical_features=["unique family identifiers"])
    if len(ids) > 1 and any(f.independent is not True or f.family_overlap is not False or f.published_family_duplication is not False for f in data.families):
        return _result("UNDETERMINED", "Cross-family independence or publication duplication is unresolved.", rule,
                       missing_critical_features=["independent, nonoverlapping families"])
    individual_ids = [i.individual_id for f in data.families for i in f.individuals]
    if len(set(individual_ids)) != len(individual_ids):
        return _result("UNDETERMINED", "The same individual appears in multiple family evidence sets.", rule,
                       missing_critical_features=["nonoverlapping individual observations"])

    count_total = 0
    count_families: list[str] = []
    lod_values: list[tuple[str, float]] = []
    count_missing: list[str] = []
    lod_missing: list[str] = []
    exclusions: list[str] = []
    for family in data.families:
        for field in ("pedigree_conflict", "sample_identity_conflict", "genotyping_artifact"):
            if getattr(family, field) is True:
                return _result("UNDETERMINED", f"{family.family_id}: {field} compromises segregation evidence.", rule,
                               exclusions_evaluated=[field])
        if family.linked_variant_evidence is True:
            return _result("UNDETERMINED", f"{family.family_id}: a linked alternative variant may explain segregation.", rule,
                           exclusions_evaluated=["linked_variant_evidence"])
        observed_affected_noncarriers = family.affected_noncarrier_count or sum(
            1 for p in family.individuals if p.affected is True and p.carries_variant is False
            and p.genotype_reliable is True and p.phenotype_consistent is True
        )
        if observed_affected_noncarriers > 0 and family.phenocopy_explained is not True:
            if family.phenocopy_explained is False:
                return _result("NOT_MET", f"{family.family_id}: confirmed affected noncarrier contradicts PP1 and phenocopy has been excluded.", rule,
                               exclusions_evaluated=["affected noncarrier"])
            return _result("UNDETERMINED", f"{family.family_id}: affected noncarrier requires phenocopy assessment.", rule,
                           missing_critical_features=["phenocopy assessment"], exclusions_evaluated=["affected noncarrier"])
        observed_unaffected_carriers = family.unaffected_carrier_count or (
            sum(1 for p in family.individuals if p.affected is False and p.carries_variant is True
                and p.genotype_reliable is True) if data.inheritance_mode == "AD" else 0
        )
        if observed_unaffected_carriers > 0:
            if family.unaffected_carriers_young is True and family.penetrance_model_validated is True:
                exclusions.append(f"{family.family_id}: young unaffected carrier not counted as nonsegregation")
            elif family.reduced_penetrance is True and family.penetrance_model_validated is True:
                exclusions.append(f"{family.family_id}: unaffected carrier reviewed under validated reduced penetrance")
            elif family.reduced_penetrance is False and family.penetrance_model_validated is True:
                return _result("NOT_MET", f"{family.family_id}: unaffected carrier conflicts with validated fully penetrant model.", rule,
                               exclusions_evaluated=["unaffected carrier"])
            else:
                return _result("UNDETERMINED", f"{family.family_id}: unaffected carrier requires age and penetrance assessment.", rule,
                               missing_critical_features=["unaffected carrier penetrance assessment"])
        count, unaffected, error = _count(family, data.inheritance_mode)
        if error is None:
            if data.inheritance_mode == "AR":
                if family.classic_ar_model_validated is not True:
                    count_missing.append(f"{family.family_id}: validated classic AR risk model")
                else:
                    # Published SVI/HL Table 4b formula; proband and carrier parents excluded.
                    value = -log10(0.25 ** count * 0.75 ** unaffected)
                    family_strength = _strength(value, rule.lod_thresholds)
                    if family_strength != "NONE":
                        lod_values.append((family.family_id, value))
                    count_total += count
                    count_families.append(family.family_id)
            else:
                count_total += count
                count_families.append(family.family_id)
        else:
            count_missing.append(f"{family.family_id}: {error}")
        if family.lod is not None:
            if _lod_valid(family.lod) and rule.lod_thresholds:
                lod_values.append((family.family_id, family.lod.score))
            else:
                lod_missing.append(f"{family.family_id}: validated LOD model and assumptions")

    # A LOD for a family supersedes a count-derived LOD from that same family.
    if rule.lod_thresholds and lod_values:
        by_family: dict[str, float] = {}
        for family_id, score in lod_values:
            by_family[family_id] = score
        total_lod = sum(by_family.values())
        strength = _strength(total_lod, rule.lod_thresholds)
        if strength != "NONE" and len(by_family) == len(data.families) and not lod_missing:
            return _result("ACTIVATED", f"Validated segregation LOD {total_lod:.6g} meets {strength.lower()} threshold.", rule,
                           strength=strength, decision_path="FORMAL_LOD" if any(f.lod for f in data.families) else "INFORMATIVE_SEGREGATION",
                           lod_score=total_lod, informative_meioses=count_total if count_families else None,
                           families_used=list(by_family), features_used=["validated segregation LOD", "family independence"],
                           exclusions_evaluated=exclusions)
    if data.inheritance_mode != "AR" and len(count_families) == len(data.families) and not count_missing:
        strength = _strength(count_total, rule.thresholds)
        if strength == "STRONG" and len(count_families) < rule.min_families_strong:
            strength = "MODERATE"
        if strength != "NONE":
            return _result("ACTIVATED", f"{count_total} validated affected segregations meet {strength.lower()} threshold.", rule,
                           strength=strength, decision_path="INFORMATIVE_SEGREGATION", informative_meioses=count_total,
                           families_used=count_families, features_used=["affected segregations", "validated count", "family independence"],
                           exclusions_evaluated=exclusions)
    missing = count_missing + lod_missing
    if missing:
        return _result("UNDETERMINED", "Essential segregation evidence or LOD validation is missing.", rule,
                       missing_critical_features=missing, exclusions_evaluated=exclusions)
    return _result("NOT_MET", "Validated segregation evidence is below the applicable PP1 threshold.", rule,
                   informative_meioses=count_total if count_families else None, exclusions_evaluated=exclusions)
