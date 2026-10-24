"""Conservative, pathway-based BS3 evaluation; no annotation or prediction."""

from __future__ import annotations

from typing import Any, Mapping

from .guidelines import MYOC, SVI, SVI_CONTROLS, WALKER, benign_strength_from_oddspath, control_count_strength
from .models import Assay, BS3Input, BS3Result, Evidence, EvidenceState, parse_input


def _value(e: Evidence) -> Any:
    return e.value if e.state == EvidenceState.KNOWN else None


def _yes(e: Evidence) -> bool:
    return _value(e) is True


def _label(assay: Assay, index: int) -> str:
    return str(_value(assay.assay_id) or f"assay_{index + 1}")


def _result(status: str, reason: str, *, strength: str = "NONE", path: str = "INSUFFICIENT_EVIDENCE",
            assays: tuple[str, ...] = (), used: tuple[str, ...] = (), missing: tuple[str, ...] = (),
            excluded: tuple[str, ...] = (), conflicts: tuple[str, ...] = (),
            guideline: Mapping[str, str] = SVI) -> BS3Result:
    return BS3Result("BS3", status, strength, path, assays, used, missing, excluded, conflicts, guideline, reason)


def _evaluate_assay(assay: Assay, mechanism: str, gene: str, disease: str, index: int) -> BS3Result:
    label = _label(assay, index)
    used: list[str] = []
    missing: list[str] = []
    excluded: list[str] = []
    for field in ("assay_type", "mechanism_compatible", "variant_class_applicable", "detects_relevant_defect", "result_class"):
        item = getattr(assay, field)
        if item.state != EvidenceState.KNOWN:
            missing.append(f"{label}.{field}")
        else:
            used.append(f"{label}.{field}")
    if _value(assay.mechanism_compatible) is False or _value(assay.variant_class_applicable) is False or _value(assay.detects_relevant_defect) is False:
        return _result("NOT_MET", "The assay cannot exclude the relevant disease mechanism or variant effect.", path="INAPPLICABLE_ASSAY", assays=(label,), used=tuple(used), excluded=("mechanism_or_variant_incompatibility",))
    if missing:
        return _result("UNDETERMINED", "Essential assay applicability or result evidence is missing or conflicting.", assays=(label,), used=tuple(used), missing=tuple(missing))
    assay_type = str(_value(assay.assay_type)).upper()
    if assay_type in ("RNA", "SPLICING", "RNA_SPLICING"):
        # Walker 2023 assigns direct normal splice observations to BP7(RNA).
        rna_missing = [f"{label}.{f}" for f in ("tissue_relevant", "transcript_relevant", "rna_detection_adequate", "nmd_addressed", "variant_allele_expressed")
                       if getattr(assay, f).state != EvidenceState.KNOWN]
        if rna_missing:
            return _result("UNDETERMINED", "RNA adequacy is unresolved; direct splicing findings belong in the BP7(RNA) framework.", path="RNA_ROUTED_TO_BP7", assays=(label,), missing=tuple(rna_missing), guideline=WALKER)
        return _result("NOT_APPLICABLE", "Direct RNA splicing evidence is evaluated under BP7(RNA), not BS3, per Walker et al. 2023.", path="RNA_ROUTED_TO_BP7", assays=(label,), used=tuple(used), guideline=WALKER)
    detectable = _value(assay.detectable_mechanisms)
    if detectable is not None and mechanism not in {str(x).upper() for x in (detectable if isinstance(detectable, (list, tuple, set)) else [detectable])}:
        return _result("NOT_MET", "The assay does not detect the disease mechanism.", path="INAPPLICABLE_ASSAY", assays=(label,), excluded=("mechanism_not_detected",))
    result_class = str(_value(assay.result_class)).upper()
    if result_class in ("NONSIGNIFICANT", "NON_SIGNIFICANT", "NO_SIGNIFICANT_DIFFERENCE"):
        return _result("UNDETERMINED", "A nonsignificant comparison does not establish normal function.", path="INDETERMINATE_RESULT", assays=(label,), missing=(f"{label}.validated_normal_result",))
    if result_class in ("DAMAGING", "ABNORMAL"):
        return _result("NOT_MET", "The observed result is damaging, so it does not support BS3.", path="DAMAGING_RESULT", assays=(label,), used=tuple(used))
    if result_class not in ("NORMAL", "INTERMEDIATE", "INDETERMINATE"):
        raise ValueError(f"{label}.result_class must be NORMAL, DAMAGING, INTERMEDIATE, or INDETERMINATE")
    if result_class != "NORMAL":
        return _result("NOT_MET", "The result is not in an accepted normal-function category.", path="NON_NORMAL_RESULT", assays=(label,), used=tuple(used))
    for field in ("replicate_consistent", "result_reliable", "assay_sensitivity_adequate",
                  "mechanism_coverage_complete", "expression_adequate",
                  "result_uncertainty_acceptable", "control_classification_quality"):
        item = getattr(assay, field)
        if _value(item) is False:
            excluded.append(field)
        elif item.state == EvidenceState.CONFLICTING:
            missing.append(f"{label}.{field}")
    for field in ("batch_effect", "artifact", "saturation", "floor_effect"):
        item = getattr(assay, field)
        if _yes(item):
            excluded.append(field)
        elif item.state == EvidenceState.CONFLICTING:
            missing.append(f"{label}.{field}")
    if excluded:
        return _result("NOT_MET", "A known experimental reliability problem invalidates the normal result.", path="RELIABILITY_EXCLUSION", assays=(label,), excluded=tuple(excluded))
    if missing:
        return _result("UNDETERMINED", "Experimental reliability has conflicting evidence.", assays=(label,), missing=tuple(missing))
    if _value(assay.assay_sensitivity_adequate) is False:
        return _result("NOT_MET", "Assay sensitivity is inadequate to exclude the relevant defect.", path="RELIABILITY_EXCLUSION", assays=(label,), excluded=("inadequate_sensitivity",))
    if _yes(assay.evidence_conflict):
        return _result("UNDETERMINED", "The functional evidence for this assay is conflicting.", path="CONFLICTING_STUDIES", assays=(label,), conflicts=(label,))
    if _yes(assay.overlaps_ps3) or _yes(assay.ps3_claimed):
        return _result("UNDETERMINED", "The same result is claimed for PS3; resolve the evidence overlap.", path="PS3_OVERLAP", assays=(label,), conflicts=(label,))
    validation = str(_value(assay.validation_status) or "UNKNOWN").upper()
    if validation == "UNVALIDATED":
        return _result("NOT_MET", "A normal result from an unvalidated assay does not establish BS3.", path="UNVALIDATED_ASSAY", assays=(label,), excluded=("unvalidated_assay",))
    if validation == "UNKNOWN":
        return _result("UNDETERMINED", "Assay validation is missing.", assays=(label,), missing=(f"{label}.validation_status",))
    if validation not in ("FORMALLY_CALIBRATED", "VALIDATED", "LIMITED", "HISTORICALLY_ACCEPTED"):
        raise ValueError(f"{label}.validation_status is invalid")
    # A numeric readout needs its own validated normal boundary. A categorical
    # readout needs explicit acceptance of the normal interpretation.
    numeric = assay.result_value.state == EvidenceState.KNOWN and isinstance(_value(assay.result_value), (int, float)) and not isinstance(_value(assay.result_value), bool)
    if numeric:
        bounds = _value(assay.validated_normal_range)
        if not isinstance(bounds, (tuple, list)) or len(bounds) != 2:
            return _result("UNDETERMINED", "A quantitative result needs a validated normal-function interval.", assays=(label,), missing=(f"{label}.validated_normal_range",))
        low, high = bounds
        if not (isinstance(low, (int, float)) and isinstance(high, (int, float)) and low <= high):
            raise ValueError("validated_normal_range must contain ordered numeric bounds")
        if not low <= _value(assay.result_value) <= high:
            return _result("NOT_MET", "Result lies outside the validated normal-function interval.", path="NON_NORMAL_RESULT", assays=(label,))
        used.extend((f"{label}.result_value", f"{label}.validated_normal_range"))
    elif not _yes(assay.accepted_normal_interpretation):
        return _result("UNDETERMINED", "A categorical normal call needs a validated accepted interpretation.", assays=(label,), missing=(f"{label}.accepted_normal_interpretation",))
    # An independently validated calibration is a complete path; per-study
    # control counts need not be repeated when the calibration is supplied.
    if _yes(assay.calibration_validated) and assay.odds_path.state == EvidenceState.KNOWN:
        strength = benign_strength_from_oddspath(float(_value(assay.odds_path)))
        guideline = SVI
        path = "VALIDATED_NORMAL_FUNCTION_ASSAY"
        used.extend((f"{label}.calibration_validated", f"{label}.odds_path"))
    elif validation == "HISTORICALLY_ACCEPTED" and _yes(assay.historical_acceptance_documented):
        strength = "SUPPORTING"
        guideline = SVI_CONTROLS
        path = "HISTORICALLY_ACCEPTED_ASSAY"
    else:
        required = ("wildtype_control", "pathogenic_control_count", "benign_control_count", "controls_independently_classified", "control_separation", "replicate_consistent")
        missing = [f"{label}.{f}" for f in required if getattr(assay, f).state != EvidenceState.KNOWN]
        if missing:
            return _result("UNDETERMINED", "The control-count pathway lacks essential validation evidence.", assays=(label,), missing=tuple(missing))
        if not all(_yes(getattr(assay, f)) for f in ("wildtype_control", "controls_independently_classified", "control_separation", "replicate_consistent")) or _value(assay.control_classification_quality) is False:
            return _result("NOT_MET", "Controls or replicates do not establish a reliable normal readout.", path="INVALID_CONTROLS", assays=(label,), excluded=("invalid_controls_or_replicates",))
        strength = control_count_strength(_value(assay.benign_control_count), _value(assay.pathogenic_control_count))
        guideline = SVI_CONTROLS
        path = "CONTROL_COUNT_NORMAL_ASSAY"
        used.extend(f"{label}.{f}" for f in required)
    if strength == "NONE":
        return _result("NOT_MET", "Validation does not support benign functional evidence at a published strength.", path=path, assays=(label,), used=tuple(used), guideline=guideline)
    if gene == "MYOC" and disease in ("PRIMARY_OPEN_ANGLE_GLAUCOMA", "GLAUCOMA", "MONDO:0005414"):
        # The Glaucoma VCEP accepts secretion/solubility and explicitly caps
        # BS3 at moderate because other mechanisms cannot be ruled out.
        if assay_type not in ("SECRETION", "SOLUBILITY"):
            return _result("NOT_APPLICABLE", "MYOC Glaucoma VCEP accepts secretion or solubility for this BS3 rule.", path="MYOC_VCEP", assays=(label,), guideline=MYOC)
        if not (_yes(assay.calibration_validated) and assay.odds_path.state == EvidenceState.KNOWN):
            return _result("UNDETERMINED", "The MYOC VCEP secretion/solubility rule requires a calibrated OddsPath.", path="MYOC_VCEP", assays=(label,), missing=(f"{label}.odds_path",), guideline=MYOC)
        if strength == "STRONG":
            strength = "MODERATE"
        guideline = MYOC
        path = "MYOC_VCEP"
    return _result("ACTIVATED", f"{label} shows validated preserved disease-relevant function at {strength.lower()} strength.", strength=strength, path=path, assays=(label,), used=tuple(used), excluded=("artifact", "saturation", "floor_effect", "batch_effect", "ps3_overlap"), guideline=guideline)


def evaluate_bs3(raw: BS3Input | Mapping[str, Any]) -> BS3Result:
    """Evaluate supplied evidence only; never fetch data or infer assay results."""
    data = parse_input(raw)
    context = data.context
    missing = [f"context.{f}" for f in ("gene_id", "disease_id", "disease_mechanism") if getattr(context, f).state != EvidenceState.KNOWN]
    if missing:
        return _result("UNDETERMINED", "Disease context is insufficient to judge functional relevance.", missing=tuple(missing))
    mechanism = str(_value(context.disease_mechanism)).upper()
    gene = str(_value(context.gene_id)).upper()
    if context.gene_disease_validity.state == EvidenceState.KNOWN and (
        _value(context.gene_disease_validity) is False or
        str(_value(context.gene_disease_validity)).upper() in ("REFUTED", "DISPUTED", "NO_KNOWN_RELATIONSHIP")
    ):
        return _result("NOT_APPLICABLE", "The supplied gene-disease relationship is not established.", path="DISEASE_CONTEXT_EXCLUSION")
    if _value(context.gene_disease_validity) is False:
        return _result("NOT_APPLICABLE", "The supplied gene-disease relationship is not established.", path="DISEASE_CONTEXT_EXCLUSION")
    if _yes(context.splice_effect_possible):
        rna_resolves_splicing = any(
            str(_value(a.assay_type) or "").upper() in ("RNA", "SPLICING", "RNA_SPLICING")
            and str(_value(a.result_class) or "").upper() == "NORMAL"
            and all(_yes(getattr(a, field)) for field in
                    ("tissue_relevant", "transcript_relevant", "rna_detection_adequate", "nmd_addressed", "variant_allele_expressed"))
            for a in data.assays
        )
        if not rna_resolves_splicing:
            return _result("UNDETERMINED", "A credible splice mechanism remains unassessed by the supplied protein assays.", path="UNRESOLVED_SPLICE_MECHANISM", missing=("relevant_splicing_evidence",))
    if _value(context.variant_mechanism_compatibility) is False:
        return _result("NOT_MET", "The variant is incompatible with the assessed disease mechanism.", path="DISEASE_CONTEXT_EXCLUSION")
    if not data.assays:
        return _result("UNDETERMINED", "No functional study was supplied.", missing=("assays",))
    disease = str(_value(context.disease_id)).upper()
    results = tuple(_evaluate_assay(a, mechanism, gene, disease, i) for i, a in enumerate(data.assays))
    active = [r for r in results if r.status == "ACTIVATED"]
    damaging = [r for i, r in enumerate(results) if r.decision_path == "DAMAGING_RESULT"
                and str(_value(data.assays[i].validation_status) or "").upper() in
                ("FORMALLY_CALIBRATED", "VALIDATED", "HISTORICALLY_ACCEPTED")
                and not any(_yes(getattr(data.assays[i], f)) for f in ("artifact", "batch_effect", "saturation", "floor_effect"))]
    explicit_conflicts = [r for r in results if r.conflicting_evidence]
    if active and (damaging or explicit_conflicts):
        labels = tuple(r.assays_used[0] for r in active + damaging + explicit_conflicts)
        return _result("UNDETERMINED", "Credible normal and damaging assays conflict; no benign strength is assigned.", path="CONFLICTING_STUDIES", assays=labels, conflicts=labels)
    if active:
        rank = {"SUPPORTING": 1, "MODERATE": 2, "STRONG": 3}
        best = max(active, key=lambda r: rank[r.strength])
        return BS3Result(best.criterion, best.status, best.strength, best.decision_path,
                         tuple(r.assays_used[0] for r in active), best.features_used,
                         best.missing_critical_features, best.exclusions_evaluated,
                         best.conflicting_evidence, best.guideline,
                         best.reason + " Independent strengths were not added.")
    unresolved = [r for r in results if r.status == "UNDETERMINED"]
    if unresolved:
        return _result("UNDETERMINED", "No complete BS3 pathway; essential evidence remains unresolved.",
                       assays=tuple(r.assays_used[0] for r in results if r.assays_used),
                       missing=tuple(x for r in unresolved for x in r.missing_critical_features),
                       conflicts=tuple(x for r in unresolved for x in r.conflicting_evidence))
    if all(r.status == "NOT_APPLICABLE" for r in results):
        return _result("NOT_APPLICABLE", "All supplied evidence belongs outside BS3.", path="NO_APPLICABLE_ASSAY", guideline=results[0].guideline)
    return _result("NOT_MET", "No supplied reliable assay establishes preserved disease-relevant function.", path="NO_QUALIFYING_ASSAY", assays=tuple(r.assays_used[0] for r in results if r.assays_used))
