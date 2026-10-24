"""Evaluation of one functional assay under the ClinGen SVI framework."""

from __future__ import annotations

from typing import Any

from .guidelines import (
    CLINGEN_SVI_2019,
    benign_strength_from_oddspath,
    strength_from_oddspath,
)
from .models import (
    AssayEvaluation,
    DecisionStatus,
    DiseaseContext,
    DiseaseMechanism,
    EvidenceState,
    FunctionalAssay,
    FunctionalResult,
    MechanismRelevance,
    Strength,
    ValidationStatus,
    token,
    unique,
)


def _bool(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    normalized = token(value)
    if normalized in {"TRUE", "YES", "PRESENT", "SUPPORTED", "ADEQUATE", "GOOD", "HIGH"}:
        return True
    if normalized in {"FALSE", "NO", "ABSENT", "UNSUPPORTED", "INADEQUATE", "POOR", "LOW"}:
        return False
    return None


def _count(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    try:
        result = int(value)
    except (TypeError, ValueError):
        return None
    return result if result >= 0 else None


def _items(value: Any) -> tuple[Any, ...]:
    if isinstance(value, (list, tuple)):
        return tuple(value)
    return (value,)


def _classification(raw: Any) -> FunctionalResult | None:
    aliases = {
        "DAMAGING": FunctionalResult.ABNORMAL,
        "FUNCTIONALLY_ABNORMAL": FunctionalResult.ABNORMAL,
        "PATHOGENIC_RANGE": FunctionalResult.ABNORMAL,
        "LOSS_OF_FUNCTION": FunctionalResult.ABNORMAL,
        "GAIN_OF_FUNCTION": FunctionalResult.ABNORMAL,
        "WILD_TYPE_LIKE": FunctionalResult.NORMAL,
        "FUNCTIONALLY_NORMAL": FunctionalResult.NORMAL,
        "BENIGN_RANGE": FunctionalResult.NORMAL,
        "PARTIAL": FunctionalResult.INTERMEDIATE,
        "HYPOMORPHIC": FunctionalResult.INTERMEDIATE,
        "UNCERTAIN": FunctionalResult.INDETERMINATE,
    }
    key = token(raw)
    if key in aliases:
        return aliases[key]
    try:
        return FunctionalResult(key)
    except ValueError:
        return None


def _mechanism(raw: Any) -> DiseaseMechanism | None:
    aliases = {
        "LOF": DiseaseMechanism.LOSS_OF_FUNCTION,
        "LOSS_OF_FUNCTION": DiseaseMechanism.LOSS_OF_FUNCTION,
        "HAPLOINSUFFICIENCY": DiseaseMechanism.LOSS_OF_FUNCTION,
        "GOF": DiseaseMechanism.GAIN_OF_FUNCTION,
        "GAIN_OF_FUNCTION": DiseaseMechanism.GAIN_OF_FUNCTION,
        "DN": DiseaseMechanism.DOMINANT_NEGATIVE,
        "DOMINANT_NEGATIVE": DiseaseMechanism.DOMINANT_NEGATIVE,
    }
    return aliases.get(token(raw), DiseaseMechanism.OTHER if raw is not None else None)


def _validation_status(raw: Any) -> ValidationStatus | None:
    aliases = {
        "FORMAL": ValidationStatus.FORMALLY_CALIBRATED,
        "ODDSPATH": ValidationStatus.FORMALLY_CALIBRATED,
        "WELL_ESTABLISHED": ValidationStatus.VALIDATED,
        "ACCEPTED": ValidationStatus.HISTORICALLY_ACCEPTED,
        "KIT": ValidationStatus.HISTORICALLY_ACCEPTED,
        "POOR": ValidationStatus.UNVALIDATED,
    }
    key = token(raw)
    if key in aliases:
        return aliases[key]
    try:
        return ValidationStatus(key)
    except ValueError:
        return None


def _result(
    assay: FunctionalAssay,
) -> tuple[FunctionalResult | None, tuple[str, ...], tuple[str, ...]]:
    used: list[str] = []
    conflicts: list[str] = []
    found: list[FunctionalResult] = []
    for name in ("variant_result_classification", "variant_result_relative_to_controls"):
        item = assay.get(name)
        if item.state == EvidenceState.CONFLICTING:
            conflicts.append(name)
        elif item.known:
            used.append(name)
            parsed = _classification(item.value)
            if parsed:
                found.append(parsed)
    if conflicts:
        return None, unique(used), unique(conflicts)
    if len(set(found)) > 1:
        return None, unique(used), ("variant functional result classifications disagree",)
    return (found[0] if found else None), unique(used), ()


def _applicability(
    assay: FunctionalAssay, context: DiseaseContext
) -> AssayEvaluation | tuple[int, tuple[str, ...], DiseaseMechanism]:
    used: list[str] = []
    missing: list[str] = []
    for name in ("gene_id", "disease_id", "disease_mechanism"):
        item = context.get(name)
        if item.state == EvidenceState.CONFLICTING:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.UNDETERMINED,
                decision_path="ASSAY_APPLICABILITY",
                reason=f"{name} is conflicting.", missing=(name,), guideline=CLINGEN_SVI_2019,
            )
        if not item.known:
            missing.append(name)
        else:
            used.append(name)
    if missing:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="ASSAY_APPLICABILITY",
            reason="The gene, disease, and established disease mechanism are required to assess assay relevance.",
            features_used=unique(used), missing=unique(missing), guideline=CLINGEN_SVI_2019,
        )
    mechanism = _mechanism(context.get("disease_mechanism").value)
    if mechanism is None:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="ASSAY_APPLICABILITY", reason="The disease mechanism is not interpretable.",
            features_used=unique(used), missing=("disease_mechanism",), guideline=CLINGEN_SVI_2019,
        )

    measures_variant = assay.get("assay_measures_variant_effect")
    if measures_variant.state == EvidenceState.CONFLICTING:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="ASSAY_APPLICABILITY", reason="Whether the assay measures the variant is conflicting.",
            features_used=unique(used), missing=("assay_measures_variant_effect",), guideline=CLINGEN_SVI_2019,
        )
    if not measures_variant.known:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="ASSAY_APPLICABILITY", reason="Variant-level interpretability is not documented.",
            features_used=unique(used), missing=("assay_measures_variant_effect",), guideline=CLINGEN_SVI_2019,
        )
    used.append("assay_measures_variant_effect")
    if _bool(measures_variant.value) is not True:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.NOT_APPLICABLE,
            decision_path="ASSAY_APPLICABILITY", reason="The assay does not measure a variant-level effect.",
            features_used=unique(used), exclusions=("not a variant-level assay",), guideline=CLINGEN_SVI_2019,
        )

    explicit = assay.get("mechanism_compatibility")
    measures = assay.get("assay_measures_disease_mechanism")
    contextual = context.get("assay_mechanism_relevance")
    for name, item in (
        ("mechanism_compatibility", explicit),
        ("assay_measures_disease_mechanism", measures),
        ("assay_mechanism_relevance", contextual),
    ):
        if item.state == EvidenceState.CONFLICTING:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.UNDETERMINED,
                decision_path="ASSAY_APPLICABILITY", reason=f"{name} is conflicting.",
                features_used=unique(used), missing=(name,), guideline=CLINGEN_SVI_2019,
            )

    # A mechanism-specific readout can expose a contradiction even if a broad
    # relevance flag was mistakenly asserted.
    mechanism_fields = {
        DiseaseMechanism.LOSS_OF_FUNCTION: "loss_of_function_effect",
        DiseaseMechanism.GAIN_OF_FUNCTION: "gain_of_function_effect",
        DiseaseMechanism.DOMINANT_NEGATIVE: "dominant_negative_effect",
    }
    opposite = [
        name for name in ("loss_of_function_effect", "gain_of_function_effect", "dominant_negative_effect")
        if name != mechanism_fields.get(mechanism) and _bool(assay.get(name).value) is True
    ]
    expected_name = mechanism_fields.get(mechanism)
    expected = assay.get(expected_name) if expected_name else None
    if opposite and (expected is None or not expected.known or _bool(expected.value) is not True):
        return AssayEvaluation(
            assay.identifier, DecisionStatus.NOT_MET,
            decision_path="MECHANISM_SPECIFIC_ASSAY",
            reason="The observed functional direction is incompatible with the established disease mechanism.",
            features_used=unique(used + opposite + ([expected_name] if expected_name else [])),
            exclusions=("disease-mechanism mismatch",), guideline=CLINGEN_SVI_2019,
        )

    relevance_rank = 0
    if explicit.known:
        used.append("mechanism_compatibility")
        if _bool(explicit.value) is False:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.NOT_APPLICABLE,
                decision_path="ASSAY_APPLICABILITY", reason="The assay is explicitly incompatible with the disease mechanism.",
                features_used=unique(used), exclusions=("disease-mechanism mismatch",), guideline=CLINGEN_SVI_2019,
            )
        if _bool(explicit.value) is True:
            relevance_rank = 3
    if measures.known:
        used.append("assay_measures_disease_mechanism")
        if _bool(measures.value) is False:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.NOT_APPLICABLE,
                decision_path="ASSAY_APPLICABILITY", reason="The assay does not measure the established disease mechanism.",
                features_used=unique(used), exclusions=("assay not relevant to disease mechanism",), guideline=CLINGEN_SVI_2019,
            )
        if _bool(measures.value) is True:
            relevance_rank = max(relevance_rank, 2)
    if contextual.known:
        used.append("assay_mechanism_relevance")
        key = token(contextual.value)
        ranks = {"DIRECT": 3, "RELEVANT": 2, "INDIRECT": 1, "NOT_RELEVANT": 0}
        relevance_rank = max(relevance_rank, ranks.get(key, 0))
        if key == MechanismRelevance.NOT_RELEVANT.value:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.NOT_APPLICABLE,
                decision_path="ASSAY_APPLICABILITY", reason="The assay is not relevant to the disease mechanism.",
                features_used=unique(used), exclusions=("assay not relevant to disease mechanism",), guideline=CLINGEN_SVI_2019,
            )
    if relevance_rank == 0:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="ASSAY_APPLICABILITY", reason="Assay relevance to the disease mechanism is not established.",
            features_used=unique(used), missing=("assay_measures_disease_mechanism",), guideline=CLINGEN_SVI_2019,
        )
    if expected and expected.known:
        used.append(expected_name)  # type: ignore[arg-type]
    return relevance_rank, unique(used), mechanism


def _control_presence(assay: FunctionalAssay) -> tuple[bool | None, bool | None, tuple[str, ...]]:
    used: list[str] = []
    positive: bool | None = None
    negative: bool | None = None
    for name, target in (("positive_control_count", "positive"), ("negative_control_count", "negative")):
        item = assay.get(name)
        if item.known:
            used.append(name)
            count = _count(item.value)
            present = count is not None and count > 0
            if target == "positive":
                positive = present
            else:
                negative = present
    for name, target in (
        ("known_pathogenic_control_results", "positive"),
        ("known_benign_control_results", "negative"),
        ("wildtype_control_result", "negative"),
    ):
        item = assay.get(name)
        if item.known:
            used.append(name)
            present = bool(_items(item.value))
            if target == "positive" and positive is None:
                positive = present
            if target == "negative" and negative is None:
                negative = present
    return positive, negative, unique(used)


def _replicates(assay: FunctionalAssay) -> tuple[bool | None, tuple[str, ...]]:
    used: list[str] = []
    answers: list[bool] = []
    consistency = assay.get("variant_result_replicate_consistency")
    if consistency.state == EvidenceState.CONFLICTING:
        return False, ("variant_result_replicate_consistency",)
    if consistency.known:
        used.append("variant_result_replicate_consistency")
        parsed = _bool(consistency.value)
        if parsed is not None:
            answers.append(parsed)
    for name in ("variant_result_replicate_count", "technical_replicates", "biological_replicates"):
        item = assay.get(name)
        if item.known:
            used.append(name)
            count = _count(item.value)
            if count is not None:
                answers.append(count >= 2)
    if False in answers:
        return False, unique(used)
    if True in answers:
        return True, unique(used)
    return None, unique(used)


def _well_understood(assay: FunctionalAssay) -> tuple[bool, tuple[str, ...]]:
    used: list[str] = []
    separation = assay.get("control_result_separation")
    dynamic = assay.get("assay_dynamic_range")
    calibration = assay.get("assay_calibration_method")
    for name, item in (
        ("control_result_separation", separation),
        ("assay_dynamic_range", dynamic),
        ("assay_calibration_method", calibration),
    ):
        if item.known:
            used.append(name)
    return (
        separation.known and _bool(separation.value) is True
        and dynamic.known and calibration.known,
        unique(used),
    )


def _validation_strength(
    assay: FunctionalAssay,
    effect: FunctionalResult,
) -> AssayEvaluation | tuple[Strength, int, str, tuple[str, ...]]:
    used: list[str] = []
    status_item = assay.get("assay_validation_status")
    if status_item.state == EvidenceState.CONFLICTING:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="ASSAY_VALIDATION", reason="Assay validation status is conflicting.",
            missing=("assay_validation_status",), guideline=CLINGEN_SVI_2019,
        )
    status = _validation_status(status_item.value) if status_item.known else None
    if status_item.known:
        used.append("assay_validation_status")

    positive, negative, control_used = _control_presence(assay)
    used.extend(control_used)
    replicate_consistency = assay.get("variant_result_replicate_consistency")
    if replicate_consistency.state == EvidenceState.CONFLICTING:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="EXPERIMENTAL_RELIABILITY",
            reason="Variant-result replicate consistency is conflicting.",
            features_used=unique(used),
            missing=("variant_result_replicate_consistency",), guideline=CLINGEN_SVI_2019,
        )
    replicate_ok, replicate_used = _replicates(assay)
    used.extend(replicate_used)
    understood, understood_used = _well_understood(assay)
    used.extend(understood_used)

    separation = assay.get("control_result_separation")
    if separation.state == EvidenceState.CONFLICTING:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="ASSAY_VALIDATION", reason="Control-result separation is conflicting.",
            features_used=unique(used), missing=("control_result_separation",), guideline=CLINGEN_SVI_2019,
        )
    if separation.known and _bool(separation.value) is False:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.NOT_MET,
            decision_path="ASSAY_VALIDATION",
            reason="The assay does not separate its normal and abnormal controls.",
            features_used=unique(used), exclusions=("controls do not separate",), guideline=CLINGEN_SVI_2019,
        )

    for name in ("assay_quality_control", "evidence_quality"):
        item = assay.get(name)
        if item.state == EvidenceState.CONFLICTING:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.UNDETERMINED,
                decision_path="EXPERIMENTAL_RELIABILITY", reason=f"{name} is conflicting.",
                features_used=unique(used), missing=(name,), guideline=CLINGEN_SVI_2019,
            )
        if item.known:
            used.append(name)
            if _bool(item.value) is False:
                return AssayEvaluation(
                    assay.identifier, DecisionStatus.NOT_MET,
                    decision_path="EXPERIMENTAL_RELIABILITY", reason=f"{name} is documented as inadequate.",
                    features_used=unique(used), exclusions=(f"inadequate {name}",), guideline=CLINGEN_SVI_2019,
                )

    for name, harmful_when_true in (
        ("experimental_artifact_evidence", True),
        ("experimental_batch_effect", True),
        ("assay_saturation_or_floor_effect", True),
        ("evidence_conflicts", True),
    ):
        item = assay.get(name)
        if item.state == EvidenceState.CONFLICTING:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.UNDETERMINED,
                decision_path="EXPERIMENTAL_RELIABILITY", reason=f"{name} is conflicting.",
                features_used=unique(used), missing=(name,), guideline=CLINGEN_SVI_2019,
            )
        if item.known:
            used.append(name)
            if _bool(item.value) is harmful_when_true:
                return AssayEvaluation(
                    assay.identifier, DecisionStatus.NOT_MET,
                    decision_path="EXPERIMENTAL_RELIABILITY",
                    reason=f"{name} prevents reliable use of this result for PS3.",
                    features_used=unique(used), exclusions=(name,), guideline=CLINGEN_SVI_2019,
                )
    reproducibility = assay.get("assay_reproducibility")
    if reproducibility.state == EvidenceState.CONFLICTING:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="EXPERIMENTAL_RELIABILITY", reason="Assay reproducibility is conflicting.",
            features_used=unique(used), missing=("assay_reproducibility",), guideline=CLINGEN_SVI_2019,
        )
    if reproducibility.known:
        used.append("assay_reproducibility")
        if _bool(reproducibility.value) is False:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.NOT_MET,
                decision_path="EXPERIMENTAL_RELIABILITY", reason="The assay is documented as not reproducible.",
                features_used=unique(used), exclusions=("assay not reproducible",), guideline=CLINGEN_SVI_2019,
            )
    if replicate_ok is False:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.NOT_MET,
            decision_path="EXPERIMENTAL_RELIABILITY", reason="Variant-result replicates are absent or inconsistent.",
            features_used=unique(used), exclusions=("replicate inconsistency",), guideline=CLINGEN_SVI_2019,
        )

    # Brnich permits Supporting for a broadly accepted, previously validated,
    # or kit assay even when the specific paper does not document controls and
    # replicates. This is the only generic bypass here.
    if status == ValidationStatus.HISTORICALLY_ACCEPTED:
        return Strength.SUPPORTING, 1, "LIMITED_VALIDATION_ASSAY", unique(used)

    if positive is False or negative is False:
        if not understood:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.NOT_MET,
                decision_path="ASSAY_VALIDATION",
                reason="The assay lacks both required abnormal/null and normal/wild-type controls.",
                features_used=unique(used), exclusions=("inadequate experimental controls",), guideline=CLINGEN_SVI_2019,
            )
    if (positive is None or negative is None or replicate_ok is None) and not understood:
        missing = []
        if positive is None:
            missing.append("positive_control_count")
        if negative is None:
            missing.append("negative_control_count")
        if replicate_ok is None:
            missing.append("technical_replicates or biological_replicates")
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="ASSAY_VALIDATION",
            reason="Experimental controls and replicates are not sufficiently documented for PS3.",
            features_used=unique(used), missing=unique(missing), guideline=CLINGEN_SVI_2019,
        )

    odds = assay.get("assay_odds_pathogenicity")
    if odds.state == EvidenceState.CONFLICTING:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="FORMAL_ODDSPATH", reason="OddsPath is conflicting.",
            features_used=unique(used), missing=("assay_odds_pathogenicity",), guideline=CLINGEN_SVI_2019,
        )
    if odds.known:
        used.append("assay_odds_pathogenicity")
        method = assay.get("assay_validation_method")
        if method.known:
            used.append("assay_validation_method")
        formal = status == ValidationStatus.FORMALLY_CALIBRATED or (
            method.known and "ODDS" in token(method.value)
        )
        if not formal:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.UNDETERMINED,
                decision_path="FORMAL_ODDSPATH",
                reason="An OddsPath value was supplied without a documented rigorous statistical calibration.",
                features_used=unique(used), missing=("assay_validation_method",), guideline=CLINGEN_SVI_2019,
            )
        try:
            value = float(odds.value)
        except (TypeError, ValueError):
            value = -1
        if value < 0:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.UNDETERMINED,
                decision_path="FORMAL_ODDSPATH", reason="OddsPath must be a non-negative number.",
                features_used=unique(used), missing=("assay_odds_pathogenicity",), guideline=CLINGEN_SVI_2019,
            )
        strength = (
            benign_strength_from_oddspath(value)
            if effect == FunctionalResult.NORMAL
            else strength_from_oddspath(value)
        )
        if strength == Strength.NONE:
            return AssayEvaluation(
                assay.identifier,
                DecisionStatus.UNDETERMINED if effect == FunctionalResult.NORMAL else DecisionStatus.NOT_MET,
                decision_path="FORMAL_ODDSPATH",
                reason=(
                    "The formal OddsPath is indeterminate for a result described as normal."
                    if effect == FunctionalResult.NORMAL
                    else "The formal OddsPath does not exceed the PS3 Supporting threshold (>2.1)."
                ),
                features_used=unique(used),
                exclusions=("OddsPath is not pathogenic evidence",), guideline=CLINGEN_SVI_2019,
            )
        ranks = {Strength.SUPPORTING: 1, Strength.MODERATE: 2, Strength.STRONG: 3, Strength.VERY_STRONG: 4}
        return strength, ranks[strength], "VALIDATED_QUANTITATIVE_ASSAY", unique(used)

    pathogenic = assay.get("known_pathogenic_control_results")
    benign = assay.get("known_benign_control_results")
    sample = assay.get("assay_validation_sample_size")
    path_n = len(_items(pathogenic.value)) if pathogenic.known else 0
    benign_n = len(_items(benign.value)) if benign.known else 0
    total = _count(sample.value) if sample.known else None
    if total is None and (path_n or benign_n):
        total = path_n + benign_n
    for name, item in (
        ("known_pathogenic_control_results", pathogenic),
        ("known_benign_control_results", benign),
        ("assay_validation_sample_size", sample),
    ):
        if item.known:
            used.append(name)

    if total and path_n and benign_n:
        quality = assay.get("control_variant_classification_quality")
        independence = assay.get("control_classification_independence")
        for name, item in (
            ("control_variant_classification_quality", quality),
            ("control_classification_independence", independence),
        ):
            if item.known:
                used.append(name)
        if quality.state == EvidenceState.CONFLICTING or independence.state == EvidenceState.CONFLICTING:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.UNDETERMINED,
                decision_path="CONTROL_COUNT_CALIBRATION", reason="Validation-control classification is conflicting.",
                features_used=unique(used), missing=("control validation",), guideline=CLINGEN_SVI_2019,
            )
        if quality.known and _bool(quality.value) is False:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.NOT_MET,
                decision_path="CONTROL_COUNT_CALIBRATION", reason="Validation-control classifications are not sufficiently high quality.",
                features_used=unique(used), exclusions=("poor validation-control classification",), guideline=CLINGEN_SVI_2019,
            )
        if independence.known and _bool(independence.value) is False:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.NOT_MET,
                decision_path="CONTROL_COUNT_CALIBRATION", reason="Validation-control classifications are circular or not independent.",
                features_used=unique(used), exclusions=("non-independent validation controls",), guideline=CLINGEN_SVI_2019,
            )
        if not quality.known or not independence.known:
            missing = []
            if not quality.known:
                missing.append("control_variant_classification_quality")
            if not independence.known:
                missing.append("control_classification_independence")
            return AssayEvaluation(
                assay.identifier, DecisionStatus.UNDETERMINED,
                decision_path="CONTROL_COUNT_CALIBRATION",
                reason="The validation controls cannot be used until classification quality and independence are documented.",
                features_used=unique(used), missing=unique(missing), guideline=CLINGEN_SVI_2019,
            )
        if total >= 11:
            return Strength.MODERATE, 2, "VALIDATED_FUNCTIONAL_ASSAY", unique(used)
        return Strength.SUPPORTING, 1, "LIMITED_VALIDATION_ASSAY", unique(used)

    if status == ValidationStatus.UNVALIDATED:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.NOT_MET,
            decision_path="ASSAY_VALIDATION", reason="The assay is explicitly unvalidated and has no accepted strength calibration.",
            features_used=unique(used), exclusions=("unvalidated assay",), guideline=CLINGEN_SVI_2019,
        )
    return AssayEvaluation(
        assay.identifier, DecisionStatus.UNDETERMINED,
        decision_path="ASSAY_VALIDATION",
        reason="No published generic strength pathway can be established: provide formal OddsPath, mixed independently classified validation controls, or documented historical acceptance.",
        features_used=unique(used),
        missing=("assay strength calibration",), guideline=CLINGEN_SVI_2019,
    )


def evaluate_assay(assay: FunctionalAssay, context: DiseaseContext) -> AssayEvaluation:
    applicability = _applicability(assay, context)
    if isinstance(applicability, AssayEvaluation):
        return applicability
    relevance_rank, applicable_used, _ = applicability

    result, result_used, result_conflicts = _result(assay)
    used = list(applicable_used + result_used)
    if result_conflicts:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="VARIANT_FUNCTIONAL_EFFECT",
            reason="The variant-specific functional result is conflicting.",
            features_used=unique(used), missing=result_conflicts,
            guideline=CLINGEN_SVI_2019, relevance_rank=relevance_rank,
        )
    if result is None:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="VARIANT_FUNCTIONAL_EFFECT",
            reason="A categorical result relative to validated normal/abnormal ranges is required; effect size or statistical significance alone is insufficient.",
            features_used=unique(used), missing=("variant_result_classification",),
            guideline=CLINGEN_SVI_2019, relevance_rank=relevance_rank,
        )
    within = assay.get("variant_result_within_assay_range")
    if within.state == EvidenceState.CONFLICTING:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="VARIANT_FUNCTIONAL_EFFECT", reason="Whether the result lies in the validated assay range is conflicting.",
            features_used=unique(used), missing=("variant_result_within_assay_range",),
            guideline=CLINGEN_SVI_2019, relevance_rank=relevance_rank, effect=result,
        )
    if within.known:
        used.append("variant_result_within_assay_range")
        if _bool(within.value) is False:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.UNDETERMINED,
                decision_path="VARIANT_FUNCTIONAL_EFFECT", reason="The result is outside the assay's validated interpretable range.",
                features_used=unique(used), exclusions=("result outside validated range",),
                guideline=CLINGEN_SVI_2019, relevance_rank=relevance_rank, effect=result,
            )

    validation = _validation_strength(assay, result)
    if isinstance(validation, AssayEvaluation):
        return AssayEvaluation(
            assay_id=validation.assay_id, status=validation.status, strength=validation.strength,
            decision_path=validation.decision_path, reason=validation.reason,
            features_used=unique(used + list(validation.features_used)), missing=validation.missing,
            exclusions=validation.exclusions, guideline=validation.guideline,
            relevance_rank=relevance_rank, validation_rank=validation.validation_rank,
            effect=result, bs3_relevant=False,
        )
    strength, validation_rank, path, validation_used = validation
    used.extend(validation_used)
    rescue = assay.get("functional_rescue_evidence")
    if rescue.known:
        used.append("functional_rescue_evidence")
    if result == FunctionalResult.NORMAL:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.NOT_MET,
            decision_path=path, reason="The validated assay result is functionally normal; PS3 is not met. The result may be evaluated separately for BS3.",
            features_used=unique(used), guideline=CLINGEN_SVI_2019,
            relevance_rank=relevance_rank, validation_rank=validation_rank,
            effect=result, bs3_relevant=True,
        )
    if result in {FunctionalResult.INTERMEDIATE, FunctionalResult.INDETERMINATE}:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path=path,
            reason="An intermediate or indeterminate functional result needs mechanism-specific expert interpretation; no generic PS3 strength is assigned.",
            features_used=unique(used), guideline=CLINGEN_SVI_2019,
            relevance_rank=relevance_rank, validation_rank=validation_rank, effect=result,
        )
    assay_type = assay.get("assay_type")
    if assay_type.known and token(assay_type.value) in {"QUALITATIVE", "CATEGORICAL"}:
        path = "VALIDATED_QUALITATIVE_ASSAY" if validation_rank >= 2 else path
        used.append("assay_type")
    return AssayEvaluation(
        assay.identifier, DecisionStatus.ACTIVATED, strength=strength,
        decision_path=path,
        reason="The functionally abnormal result is compatible with the disease mechanism and the assigned strength is supported by the assay's clinical validation.",
        features_used=unique(used), guideline=CLINGEN_SVI_2019,
        relevance_rank=relevance_rank, validation_rank=validation_rank, effect=result,
    )
