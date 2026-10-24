"""Evaluate only BS2 using observed clinical and genotype evidence."""

from dataclasses import fields
from typing import Any

from .guidelines import ACMG_2015, RETT_ANGELMAN, RETT_GENOTYPES, RUNX1
from .models import BS2Input, BS2Result, DiseaseContext, Evidence, EvidenceState, HealthyObservation


def _fact(value: Evidence) -> Any:
    return value.value if value.state == EvidenceState.KNOWN else None


def _model(cls: type, raw: dict[str, Any]):
    names = {f.name for f in fields(cls)}
    unknown = set(raw) - names
    if unknown:
        raise ValueError(f"Unknown {cls.__name__} fields: {sorted(unknown)}")
    data = dict(raw)
    for f in fields(cls):
        if f.name in data and (f.type is Evidence or f.name in _EVIDENCE_FIELDS):
            data[f.name] = Evidence.parse(data[f.name])
    return cls(**data)


_EVIDENCE_FIELDS = {
    "gene_disease_validity", "inheritance_mode", "molecular_mechanism_compatible",
    "full_penetrance", "early_onset", "age_dependent_penetrance",
    "variable_expressivity", "relevant_subclinical_features", "genotype",
    "genotype_reliable", "confirmed_unaffected", "adult", "beyond_expected_onset",
    "clinical_assessment_adequate", "independent", "maternal_allele", "artifact",
    "deletion_mimicking_homozygosity", "mosaicism", "protective_modifier",
    "genetic_rescue", "subclinical_manifestation", "phenotype_misclassification",
    "variant_reliable", "vcep_applicable",
}


def _parse(raw: BS2Input | dict[str, Any]) -> BS2Input:
    if isinstance(raw, BS2Input):
        return raw
    if not isinstance(raw, dict):
        raise TypeError("BS2 input must be BS2Input or a mapping")
    data = dict(raw)
    disease = data.pop("disease", {})
    observations = data.pop("observations", ())
    data["disease"] = disease if isinstance(disease, DiseaseContext) else _model(DiseaseContext, disease)
    data["observations"] = tuple(o if isinstance(o, HealthyObservation) else _model(HealthyObservation, o)
                                 for o in observations)
    return _model(BS2Input, data)


def _result(status: str, reason: str, guideline: dict[str, str], *, strength: str = "NONE",
            path: str | None = None, count: int = 0, used: set[str] | None = None,
            missing: set[str] | None = None, exclusions: set[str] | None = None) -> BS2Result:
    return BS2Result("BS2", status, strength, path, count, tuple(sorted(used or ())),
                     tuple(sorted(missing or ())), tuple(sorted(exclusions or ())),
                     guideline, reason)


def _genotype_eligible(mode: str, genotype: str) -> bool:
    return (mode == "AD" and genotype == "HETEROZYGOUS") or (mode == "AR" and genotype == "HOMOZYGOUS") or (mode == "XL" and genotype == "HEMIZYGOUS")


def _assess_observation(o: HealthyObservation, *, mode: str | None,
                        allowed: tuple[str, ...] | None = None,
                        require_maternal: bool = False) -> tuple[str, set[str], set[str], set[str]]:
    """Return eligible, uncertain or failed; never infer a clinical phenotype."""
    missing: set[str] = set()
    used: set[str] = set()
    exclusions: set[str] = set()
    genotype = _fact(o.genotype)
    if genotype is None:
        missing.add(f"{o.individual_id}.genotype")
    elif not (genotype in allowed if allowed is not None else _genotype_eligible(mode or "", genotype)):
        return "failed", missing, used, exclusions
    else:
        used.add(f"{o.individual_id}.genotype")

    for name in ("genotype_reliable", "confirmed_unaffected", "adult",
                 "beyond_expected_onset", "clinical_assessment_adequate"):
        value = _fact(getattr(o, name))
        if value is False:
            return "failed" if name in ("genotype_reliable", "confirmed_unaffected", "clinical_assessment_adequate") else "uncertain", missing, used, exclusions
        if value is None:
            missing.add(f"{o.individual_id}.{name}")
        else:
            used.add(f"{o.individual_id}.{name}")
    if o.source_type == "POPULATION_DATABASE" and _fact(o.confirmed_unaffected) is not True:
        missing.add(f"{o.individual_id}.confirmed_unaffected")
    if require_maternal:
        if _fact(o.maternal_allele) is False:
            return "failed", missing, used, exclusions
        if _fact(o.maternal_allele) is None:
            missing.add(f"{o.individual_id}.maternal_allele")
        else:
            used.add(f"{o.individual_id}.maternal_allele")
    for name in ("artifact", "deletion_mimicking_homozygosity", "mosaicism",
                 "protective_modifier", "genetic_rescue", "subclinical_manifestation",
                 "phenotype_misclassification"):
        evidence = getattr(o, name)
        value = _fact(evidence)
        if value is True:
            exclusions.add(f"{o.individual_id}.{name}")
            return "failed", missing, used, exclusions
        if value is False:
            exclusions.add(f"{o.individual_id}.{name}")
        if evidence.state == EvidenceState.CONFLICTING:
            missing.add(f"{o.individual_id}.{name}")
    return ("uncertain" if missing else "eligible"), missing, used, exclusions


def evaluate_bs2(biological_features: BS2Input | dict[str, Any]) -> BS2Result:
    """Return a deterministic BS2 decision from supplied evidence.

    General ACMG applies only to fully penetrant early-onset AD/AR/XL disease.
    Pinned VCEP rules supersede it in their attested gene/disease scope.
    """
    data = _parse(biological_features)
    disease = data.disease
    gene = (disease.gene_id or "").upper()
    used: set[str] = set()
    missing: set[str] = set()
    exclusions: set[str] = set()
    if gene == "RUNX1" and disease.disease_id == "FPD_AML":
        return _result("NOT_APPLICABLE", "RUNX1 FPD/AML VCEP explicitly excludes BS2 because penetrance is incomplete and malignancy onset is later.", RUNX1)

    rett = (gene in RETT_GENOTYPES and disease.disease_id == "RETT_ANGELMAN_LIKE"
            and _fact(data.vcep_applicable) is True
            and data.vcep_specification_version == "2022 publication")
    rett_scope = gene in RETT_GENOTYPES and disease.disease_id == "RETT_ANGELMAN_LIKE"
    if gene in RETT_GENOTYPES and disease.disease_id == "RETT_ANGELMAN_LIKE" and _fact(data.vcep_applicable) is None:
        missing.add("vcep_applicable")
    if gene in RETT_GENOTYPES and disease.disease_id == "RETT_ANGELMAN_LIKE" and _fact(data.vcep_applicable) is True and not rett:
        missing.add("vcep_specification_version")
    if rett_scope and not rett:
        return _result("UNDETERMINED", "The in-scope VCEP rule requires confirmed applicability and its pinned version.", RETT_ANGELMAN,
                       missing=missing or {"vcep_applicable"})
    if _fact(data.variant_reliable) is False:
        return _result("NOT_MET", "Variant identity or call is unreliable.", RETT_ANGELMAN if rett else ACMG_2015)
    if _fact(data.variant_reliable) is None:
        missing.add("variant_reliable")
    else:
        used.add("variant_reliable")
    if not gene:
        missing.add("gene_id")
    if not disease.disease_id:
        missing.add("disease_id")
    if _fact(disease.gene_disease_validity) is False:
        return _result("NOT_APPLICABLE", "The gene-disease relationship is invalid for this context.", ACMG_2015)
    if _fact(disease.gene_disease_validity) is None:
        missing.add("gene_disease_validity")
    else:
        used.add("gene_disease_validity")
    if _fact(disease.molecular_mechanism_compatible) is False:
        return _result("NOT_MET", "The candidate variant is incompatible with the disease mechanism.", ACMG_2015)
    if _fact(disease.molecular_mechanism_compatible) is None:
        missing.add("molecular_mechanism_compatible")
    else:
        used.add("molecular_mechanism_compatible")
    mode = _fact(disease.inheritance_mode)
    if mode not in ("AD", "AR", "XL"):
        if mode is None:
            missing.add("inheritance_mode")
        else:
            return _result("NOT_APPLICABLE", "No verified BS2 pathway for this inheritance mode.", ACMG_2015)
    else:
        used.add("inheritance_mode")

    guideline = RETT_ANGELMAN if rett else ACMG_2015
    if not rett:
        for name in ("full_penetrance", "early_onset"):
            value = _fact(getattr(disease, name))
            if value is False:
                return _result("NOT_MET", f"The general BS2 rule requires {name.replace('_', ' ')}; no applicable alternative rule is pinned.", guideline)
            if value is None:
                missing.add(name)
            else:
                used.add(name)
        if _fact(disease.age_dependent_penetrance) is True:
            return _result("NOT_MET", "Age-dependent penetrance prevents use of the general fully penetrant early-onset rule.", guideline)
        if _fact(disease.variable_expressivity) is True and _fact(disease.relevant_subclinical_features) is not False:
            missing.add("relevant_subclinical_features")
    if not data.observations:
        missing.add("observations")
    eligible: list[HealthyObservation] = []
    uncertain = False
    observation_missing: set[str] = set()
    seen: set[str] = set()
    for o in data.observations:
        if not o.individual_id or o.individual_id in seen:
            exclusions.add("duplicate_or_missing_individual_id")
            continue
        seen.add(o.individual_id)
        state, omissing, oused, oexclusions = _assess_observation(
            o, mode=mode, allowed=RETT_GENOTYPES[gene] if rett else None,
            require_maternal=False)
        if rett and gene == "UBE3A" and state == "eligible":
            if _fact(o.maternal_allele) is True:
                oused.add(f"{o.individual_id}.maternal_allele")
            elif _fact(o.independent) is True and o.related_group is None:
                oused.add(f"{o.individual_id}.independent")
            elif _fact(o.maternal_allele) is False or _fact(o.independent) is False:
                state = "failed"
            else:
                state = "uncertain"
                omissing.add(f"{o.individual_id}.maternal_allele_or_unrelated")
        exclusions.update(oexclusions)
        if state == "eligible":
            eligible.append(o)
            used.update(oused)
        elif state == "uncertain":
            uncertain = True
            observation_missing.update(omissing)

    if rett:
        # Published VCEP table: two for five genes, four for UBE3A.
        threshold = 4 if gene == "UBE3A" else 2
        if len(eligible) >= threshold and not missing:
            return _result("ACTIVATED", f"{len(eligible)} documented unaffected adults satisfy the {gene} VCEP BS2 count of {threshold}.", guideline,
                           strength="STRONG", path="VCEP_RETT_ANGELMAN", count=len(eligible), used=used, exclusions=exclusions)
        missing.update(observation_missing)
        status = "UNDETERMINED" if missing or uncertain else "NOT_MET"
        return _result(status, f"{len(eligible)} eligible observations; {threshold} required by the {gene} VCEP rule.", guideline,
                       count=len(eligible), used=used, missing=missing, exclusions=exclusions)

    if eligible and not missing:
        return _result("ACTIVATED", "A reliably genotyped, clinically confirmed unaffected adult passed the expected early onset of a fully penetrant disorder.", guideline,
                       strength="STRONG", path="UNAFFECTED_EXPECTED_PATHOGENIC_GENOTYPE", count=len(eligible), used=used, exclusions=exclusions)
    if not eligible:
        missing.update(observation_missing)
    status = "UNDETERMINED" if missing or (uncertain and not eligible) else "NOT_MET"
    return _result(status, "No complete eligible BS2 observation is established." if status == "NOT_MET" else "Essential BS2 evidence is missing or conflicting.",
                   guideline, count=len(eligible), used=used, missing=missing, exclusions=exclusions)
