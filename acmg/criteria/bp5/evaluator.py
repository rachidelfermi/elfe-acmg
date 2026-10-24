"""Conservative BP5-only evaluator; no variant or phenotype inference."""

from math import isfinite
from typing import Any, Mapping

from .guidelines import ACMG_2015, Specification, find_specification
from .models import BP5Features, BP5Result, EvidenceState


_EXCLUSIONS = (
    "dual_diagnosis_evidence", "blended_phenotype_evidence",
    "candidate_specific_phenotype_evidence", "candidate_additional_phenotype_contribution",
    "oligogenic_disease_evidence", "modifier_variant_evidence",
    "clinical_diagnostic_uncertainty", "alternative_diagnosis_conflicts",
)


def evaluate_bp5(biological_features: BP5Features | Mapping[str, Any]) -> BP5Result:
    if isinstance(biological_features, Mapping):
        features = BP5Features.from_mapping(biological_features)
    elif isinstance(biological_features, BP5Features):
        features = biological_features
    else:
        raise TypeError("biological_features must be BP5Features or a mapping")

    used: list[str] = []
    missing: list[str] = []
    exclusions: list[str] = []
    spec: Specification | None = None

    def read(name: str, *, required: bool = False) -> Any:
        e = features.evidence(name)
        if e.state == EvidenceState.KNOWN:
            used.append(name)
            return e.value
        if required or e.state == EvidenceState.CONFLICTING:
            missing.append(name)
        return None

    def result(status: str, reason: str, path: str | None = None,
               strength: str = "NONE") -> BP5Result:
        diagnosis = features.evidence("alternative_diagnosis_id")
        return BP5Result("BP5", status, strength, path,
                         str(diagnosis.value) if diagnosis.state == EvidenceState.KNOWN else None,
                         tuple(dict.fromkeys(used)), tuple(dict.fromkeys(missing)),
                         tuple(dict.fromkeys(exclusions)),
                         spec.guideline if spec else ACMG_2015.copy(), reason)

    variant = read("candidate_variant_id", required=True)
    gene = read("candidate_gene_id", required=True)
    disease = read("candidate_disease_id", required=True)
    if not all((variant, gene, disease)):
        return result("UNDETERMINED", "Candidate variant, gene, and disease context are needed to assess BP5.")
    spec = find_specification(str(gene).upper(), str(disease).upper())
    requested = read("vcep_specification_id")
    version = read("vcep_specification_version")
    if requested is not None and (spec is None or requested != spec.spec_id):
        return result("UNDETERMINED", "Requested VCEP specification is not pinned for this gene and disease.")
    if version is not None and (spec is None or str(version) != spec.version):
        return result("UNDETERMINED", "Requested VCEP version is not the pinned version.")
    if spec and spec.mode == "DISABLED":
        exclusions.append("VCEP_BP5_DISABLED")
        return result("NOT_APPLICABLE", "This VCEP explicitly prohibits BP5 for the scoped disorder.")
    if spec and spec.mode == "REQUIRES_CASE_SCORING":
        missing.append("vcep_bp5_case_score")
        return result("UNDETERMINED", "This VCEP requires a disease-specific case-point table that is not implemented.")

    # Check observed complications before accepting a minimal positive pathway.
    for name in _EXCLUSIONS:
        e = features.evidence(name)
        if e.state == EvidenceState.CONFLICTING:
            missing.append(name)
            return result("UNDETERMINED", f"Conflicting evidence about {name.replace('_', ' ')}.")
        if e.state == EvidenceState.KNOWN:
            used.append(name)
            exclusions.append(name)
            if e.value is True:
                return result("NOT_MET" if name not in ("clinical_diagnostic_uncertainty", "alternative_diagnosis_conflicts")
                              else "UNDETERMINED",
                              f"{name.replace('_', ' ').capitalize()} prevents a BP5 inference.")

    independent = read("candidate_independent_disease_evidence", required=not (spec and spec.mode == "CLINICAL_LR"))
    if independent is True:
        exclusions.append("candidate_independent_disease_evidence")
        return result("NOT_MET", "The candidate has evidence for an independent disease contribution.")

    # Avoid counting the same observation twice. The caller supplies stable evidence IDs.
    bp5_ids = read("bp5_evidence_ids")
    other_ids = read("other_criteria_evidence_ids")
    if bp5_ids is not None and other_ids is not None:
        if not isinstance(bp5_ids, (list, tuple, set, frozenset)) or not isinstance(other_ids, Mapping):
            raise ValueError("bp5_evidence_ids must be a sequence and other_criteria_evidence_ids a mapping")
        overlaps = {"PP4", "PP1", "BS4", "PS2", "PM6"} & set(other_ids)
        if any(set(bp5_ids) & set(other_ids[code]) for code in overlaps):
            exclusions.append("CORRELATED_EVIDENCE")
            return result("UNDETERMINED", "The same observation is assigned to BP5 and another phenotype or family criterion.")

    if spec and spec.mode == "CLINICAL_LR":
        return _evaluate_clinical_lr(features, read, result, exclusions)

    if independent is None:
        return result("UNDETERMINED", "Whether the candidate contributes independently is unresolved.")

    # An isolated recessive carrier allele is not made benign by another diagnosis.
    inheritance = read("candidate_disease_inheritance")
    zygosity = read("candidate_zygosity")
    if str(inheritance).upper() in ("RECESSIVE", "AUTOSOMAL_RECESSIVE") and \
            str(zygosity).upper() in ("HETEROZYGOUS", "HET"):
        exclusions.append("RECESSIVE_CARRIER")
        return result("NOT_MET", "Another diagnosis does not establish benignity of a recessive carrier allele.")

    alt = read("alternative_diagnosis_id", required=True)
    alt_variant = read("alternative_variant_id", required=True)
    classification = read("alternative_variant_classification", required=True)
    classification_source = read("alternative_variant_classification_source")
    evidence = read("alternative_variant_classification_evidence", required=True)
    validity = read("alternative_gene_disease_validity", required=True)
    compatible = read("alternative_inheritance_compatibility", required=True)
    if not alt or not alt_variant:
        return result("UNDETERMINED", "A specific alternative molecular diagnosis and variant are required.")
    if classification is None or evidence is None or validity is None or compatible is None:
        return result("UNDETERMINED", "Alternative pathogenicity, disease validity, or genotype compatibility is unresolved.")
    if not classification_source and not features.evidence("alternative_variant_classification_evidence").source:
        missing.append("alternative_variant_classification_source")
        return result("UNDETERMINED", "A traceable source for the alternative classification is required.")
    if str(classification).upper().replace(" ", "_") not in ("PATHOGENIC", "LIKELY_PATHOGENIC"):
        return result("NOT_MET", "The alternative variant is not classified pathogenic or likely pathogenic.")
    if not evidence or validity is False or compatible is False:
        return result("NOT_MET", "The alternative molecular diagnosis is not established by the supplied evidence.")
    if validity is not True or compatible is not True:
        raise ValueError("Alternative validity and inheritance compatibility must be Boolean")

    independence = read("alternative_variant_classification_independence", required=True)
    if independence is False:
        exclusions.append("CIRCULAR_ALTERNATIVE_CLASSIFICATION")
        return result("UNDETERMINED", "Alternative pathogenicity depends on the BP5 conclusion being assessed.")
    if independence is None:
        return result("UNDETERMINED", "The alternative classification needs independently reviewable evidence; an assertion alone is insufficient.")

    confidence = read("alternative_diagnosis_confidence")
    if confidence is not None and str(confidence).upper() in ("LOW", "SUSPECTED", "UNCERTAIN"):
        return result("UNDETERMINED", "The alternative diagnosis is not established.")
    quality = read("alternative_diagnosis_confirmation_quality")
    if quality is not None and str(quality).upper() in ("LOW", "UNRELIABLE"):
        return result("UNDETERMINED", "Alternative diagnosis confirmation has insufficient quality.")

    phenotype = read("patient_phenotype_present")
    if phenotype is None:
        hpo = read("patient_hpo_terms")
        findings = read("patient_clinical_findings")
        phenotype = bool(hpo or findings) if hpo is not None or findings is not None else None
    if phenotype is False:
        return result("NOT_MET", "No relevant patient phenotype is supplied.")
    if phenotype is None:
        missing.append("patient_phenotype_present")
        return result("UNDETERMINED", "The patient's relevant phenotype is unknown.")

    core = read("alternative_diagnosis_explains_core_phenotypes")
    all_relevant = read("alternative_diagnosis_explains_all_relevant_phenotypes", required=True)
    if core is False or all_relevant is False:
        return result("NOT_MET", "The alternative diagnosis does not explain the relevant clinical presentation.")
    if all_relevant is None:
        return result("UNDETERMINED", "Whether the alternative diagnosis explains all relevant features is unresolved.")
    unexplained = read("unexplained_phenotype_features")
    if unexplained:
        exclusions.append("UNEXPLAINED_PHENOTYPE")
        return result("UNDETERMINED", "Unexplained findings conflict with a complete phenotype explanation.")
    if core is None and features.evidence("alternative_diagnosis_explains_core_phenotypes").state == EvidenceState.CONFLICTING:
        return result("UNDETERMINED", "Core phenotype explanation is conflicting.")
    if features.evidence("phenotype_evidence_quality").state == EvidenceState.KNOWN:
        if str(read("phenotype_evidence_quality")).upper() in ("LOW", "UNRELIABLE"):
            return result("UNDETERMINED", "Phenotype evidence has insufficient quality.")
    return result("ACTIVATED", "An established alternative molecular diagnosis explains the relevant phenotype; no independent candidate contribution is established in the supplied assessment.",
                  "ESTABLISHED_ALTERNATIVE_DIAGNOSIS", "SUPPORTING")


def _evaluate_clinical_lr(features: BP5Features, read: Any, result: Any,
                          exclusions: list[str]) -> BP5Result:
    if read("breast_ovarian_coobservation_only") is True:
        exclusions.append("BREAST_OVARIAN_COOBSERVATION")
        return result("NOT_APPLICABLE", "BRCA1 BP5 does not use co-observation of variants in breast–ovarian risk genes.")
    calibrated = read("clinically_calibrated_lr", required=True)
    description = read("lr_evidence_description", required=True)
    raw_lr = read("clinical_likelihood_ratio", required=True)
    if calibrated is None or not description or raw_lr is None:
        return result("UNDETERMINED", "The VCEP requires a reviewable, clinically calibrated likelihood ratio.")
    if calibrated is not True:
        return result("NOT_MET", "The supplied likelihood ratio is not clinically calibrated.")
    try:
        lr = float(raw_lr)
    except (TypeError, ValueError) as exc:
        raise ValueError("clinical_likelihood_ratio must be numeric") from exc
    if not isfinite(lr) or lr <= 0:
        raise ValueError("clinical_likelihood_ratio must be finite and positive")
    if lr > 0.48:
        return result("NOT_MET", "The combined clinical likelihood ratio does not reach BRCA1 BP5 Supporting.")
    # The published VCEP also permits Very Strong at <=0.00285. This API caps
    # output at STRONG, the highest strength requested by the public contract.
    strength = "STRONG" if lr <= 0.05 else "MODERATE" if lr <= 0.23 else "SUPPORTING"
    return result("ACTIVATED", f"ClinGen BRCA1 v1.2 calibrated clinical LR {lr:g} meets BP5 {strength}; evidence must not be reused in other criteria.",
                  "VCEP_CLINICAL_LR", strength)
