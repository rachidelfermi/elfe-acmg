"""Conservative PP4 evaluator: a complete pathway activates without unrelated fields."""
from __future__ import annotations

from .guidelines import ACMG, ADA, ADA_DISEASE, ADA_GENE, GT, GT_DISEASE, GT_GENES, CDH1, CDH1_DISEASE, CDH1_GENE
from .models import Evidence, PP4Features, PP4Result, State


def _yes(item: Evidence[bool]) -> bool:
    return item.state == State.KNOWN and item.value is True


def _no(item: Evidence[bool]) -> bool:
    return item.state == State.KNOWN and item.value is False


def _result(status: str, reason: str, *, strength: str = "NONE", path: str = "NONE",
            used: tuple[str, ...] = (), missing: tuple[str, ...] = (),
            exclusions: tuple[str, ...] = (), guideline=ACMG) -> PP4Result:
    return PP4Result("PP4", status, strength, path, used, missing, exclusions, guideline, reason)


def evaluate_pp4(features: PP4Features) -> PP4Result:
    if not isinstance(features, PP4Features):
        raise TypeError("evaluate_pp4 expects PP4Features")
    if not features.gene_id or not features.disease_id:
        return _result("UNDETERMINED", "Gene and disease identifiers are required.",
                       missing=tuple(n for n in ("gene_id", "disease_id") if not getattr(features, n)))
    if features.gene_id in ADA_GENE and features.disease_id in ADA_DISEASE:
        guideline = ADA
    elif features.gene_id in GT_GENES and features.disease_id in GT_DISEASE:
        guideline = GT
    elif features.gene_id in CDH1_GENE and features.disease_id in CDH1_DISEASE:
        guideline = CDH1
    else:
        guideline = ACMG
    if guideline == CDH1:
        return _result("NOT_APPLICABLE", "The CDH1 VCEP assigns qualifying phenotype evidence through PS4 and does not apply PP4.",
                       guideline=CDH1)
    if features.gene_disease_validity.state == State.KNOWN and features.gene_disease_validity.value in {"DISPUTED", "REFUTED", "NO_KNOWN_RELATIONSHIP"}:
        return _result("NOT_MET", "Gene–disease relationship is unsupported.", guideline=guideline,
                       used=("gene_disease_validity",))
    if _no(features.genotype_compatible) or _no(features.phenotype_compatible) or _yes(features.alternative_diagnosis_confirmed):
        return _result("NOT_MET", "A documented incompatibility or confirmed alternative diagnosis excludes this PP4 application.",
                       guideline=guideline, exclusions=("genotype_compatible", "phenotype_compatible", "alternative_diagnosis_confirmed"))
    if _yes(features.phenotype_conflict) or features.phenotype_conflict.state == State.CONFLICTING or _no(features.evidence_reliable):
        return _result("UNDETERMINED", "Clinical evidence conflicts or is unreliable.", guideline=guideline,
                       missing=("clinical_evidence_review",))
    if features.gene_disease_validity.state != State.KNOWN or features.gene_disease_validity.value not in {"DEFINITIVE", "STRONG", "MODERATE"}:
        return _result("UNDETERMINED", "Established gene–disease validity is not documented.", guideline=guideline,
                       missing=("gene_disease_validity",))
    if _no(features.age_sex_applicable):
        return _result("NOT_MET", "The documented manifestation is not applicable at the patient's age or sex.", guideline=guideline,
                       exclusions=("age_sex_applicable",))
    if guideline == ADA:
        return _evaluate_ada(features)
    if guideline == GT:
        return _evaluate_gt(features)
    return _evaluate_general(features)


def _evaluate_general(f: PP4Features) -> PP4Result:
    if _no(f.single_genetic_etiology):
        return _result("NOT_MET", "This presentation has multiple genetic etiologies and no scoped exception is encoded.",
                       used=("single_genetic_etiology",))
    personal = _yes(f.phenotype_highly_specific)
    family = _yes(f.family_history_highly_specific) and not _yes(f.family_evidence_used_for_pp1)
    if not personal and _yes(f.family_history_highly_specific) and _yes(f.family_evidence_used_for_pp1):
        return _result("UNDETERMINED", "Family evidence is already assigned to PP1; independent PP4 evidence is needed.",
                       missing=("independent_family_specificity_evidence",))
    if _no(f.phenotype_highly_specific) and _no(f.family_history_highly_specific):
        return _result("NOT_MET", "Neither patient phenotype nor family history is highly specific.",
                       used=("phenotype_highly_specific", "family_history_highly_specific"))
    missing = []
    if f.single_genetic_etiology.state != State.KNOWN:
        missing.append("single_genetic_etiology")
    if not personal and not family:
        missing.append("phenotype_or_family_history_highly_specific")
    if missing:
        return _result("UNDETERMINED", "Essential specificity evidence is missing.", missing=tuple(missing))
    field = "phenotype_highly_specific" if personal else "family_history_highly_specific"
    return _result("ACTIVATED", "A documented highly specific presentation and single genetic etiology meet ACMG/AMP PP4.",
                   strength="SUPPORTING", path="HIGHLY_SPECIFIC_PHENOTYPE" if personal else "HIGHLY_SPECIFIC_FAMILY_HISTORY",
                   used=("gene_disease_validity", "single_genetic_etiology", field),
                   exclusions=("alternative_diagnosis_confirmed", "phenotype_conflict"))


def _evaluate_ada(f: PP4Features) -> PP4Result:
    # ClinGen GN114 v2.2. Each component is counted at most once per patient.
    points = 0.0
    used = []
    def add(field: str, value: float) -> None:
        nonlocal points
        points += value
        used.append(field)
    if _yes(f.scid_diagnostic_criteria_met):
        add("scid_diagnostic_criteria_met", 0.5)
    if _yes(f.scid_sequencing_no_alternative):
        add("scid_sequencing_no_alternative", 1.0)
        if _yes(f.scid_family_history) and not _yes(f.family_evidence_used_for_pp1):
            add("scid_family_history", 0.5)
    enzyme = (f.ada_enzyme_activity_percent.state == State.KNOWN
              and 0 <= f.ada_enzyme_activity_percent.value < 1)
    dado = _yes(f.ada_dado_above_reference) and _yes(f.ada_dado_pre_treatment_or_nontransfused)
    if enzyme or dado:
        add("ada_enzyme_activity_percent" if enzyme else "ada_dado_above_reference", 5.0)
    treatment = _yes(f.ada_phenotype_corrected_by_treatment)
    if treatment:
        add("ada_phenotype_corrected_by_treatment", 6.0 if _yes(f.ada_cnv_testing_done) else 4.5)
    if _yes(f.scid_t_minus_b_minus_nk_minus):
        add("scid_t_minus_b_minus_nk_minus", 0.5)
    # The VCEP permits 4.5 treatment points without CNV testing, but requires
    # CNV testing if those points are necessary to reach Strong.
    independent = points - (6.0 if treatment and _yes(f.ada_cnv_testing_done) else 4.5 if treatment else 0.0)
    if points >= 6 and treatment and not _yes(f.ada_cnv_testing_done) and independent < 6:
        points = min(points, 5.5)
    if points >= 6:
        strength = "STRONG"
    elif points >= 2:
        strength = "MODERATE"
    elif points >= 1:
        strength = "SUPPORTING"
    else:
        strength = "NONE"
    if strength != "NONE":
        return _result("ACTIVATED", f"ADA VCEP v2.2 documented score: {points:g} points.", strength=strength,
                       path="ADA_VCEP", used=("gene_disease_validity", *used), guideline=ADA,
                       exclusions=("alternative_diagnosis_confirmed", "phenotype_conflict"))
    relevant = (f.scid_diagnostic_criteria_met, f.scid_sequencing_no_alternative,
                f.ada_enzyme_activity_percent, f.ada_dado_above_reference,
                f.ada_phenotype_corrected_by_treatment, f.scid_t_minus_b_minus_nk_minus)
    if any(x.state in {State.UNKNOWN, State.CONFLICTING} for x in relevant):
        return _result("UNDETERMINED", f"Documented ADA score is {points:g}; unassessed components could change it.",
                       path="ADA_VCEP", used=tuple(used), guideline=ADA,
                       missing=("ada_pp4_point_components",))
    return _result("NOT_MET", f"Documented ADA score {points:g} is below the VCEP's 1-point minimum.",
                   path="ADA_VCEP", used=tuple(used), guideline=ADA)


def _evaluate_gt(f: PP4Features) -> PP4Result:
    essential = ("gt_clinical_diagnosis", "gt_mucocutaneous_bleeding", "gt_normal_ristocetin",
                 "gt_minimal_aggregation_at_least_two_agonists")
    if any(_no(getattr(f, name)) for name in essential):
        return _result("NOT_MET", "A required clinical or platelet aggregation finding is absent.",
                       path="GT_VCEP", used=essential, guideline=GT)
    missing = tuple(name for name in essential if not _yes(getattr(f, name)))
    if missing:
        return _result("UNDETERMINED", "The GT clinical diagnosis and aggregation requirements are incomplete.",
                       path="GT_VCEP", missing=missing, guideline=GT)
    confirmation = (_yes(f.gt_aIIbbeta3_expression_below_25_percent)
                    or _yes(f.gt_impaired_fibrinogen_binding))
    if confirmation and _yes(f.gt_both_genes_fully_sequenced):
        strength = "STRONG"
    else:
        strength = "MODERATE"
    used = essential + (("gt_aIIbbeta3_expression_below_25_percent",) if _yes(f.gt_aIIbbeta3_expression_below_25_percent)
                        else ("gt_impaired_fibrinogen_binding",) if _yes(f.gt_impaired_fibrinogen_binding) else ())
    if strength == "STRONG":
        used += ("gt_both_genes_fully_sequenced",)
    return _result("ACTIVATED", "Clinician-diagnosed GT with mucocutaneous bleeding and qualifying aggregation findings meets the VCEP PP4 rule.",
                   strength=strength, path="GT_VCEP", used=("gene_disease_validity", *used), guideline=GT,
                   exclusions=("alternative_diagnosis_confirmed", "phenotype_conflict"))
