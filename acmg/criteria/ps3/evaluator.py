"""Public orchestration and multi-assay resolution for PS3."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from .functional_assays import evaluate_assay
from .guidelines import CLINGEN_SVI_2019, parse_vcep_specification, vcep_reference
from .models import (
    ASSAY_FIELDS,
    CONTEXT_FIELDS,
    AssayEvaluation,
    DecisionStatus,
    DiseaseContext,
    EvidenceState,
    EvidenceValue,
    FunctionalAssay,
    FunctionalResult,
    PS3Decision,
    PS3Evidence,
    SchemaValidationError,
    Strength,
    VCEPRule,
    VCEPSpecification,
    token,
    unique,
)


GLOBAL_FIELDS = (
    "independent_assay_count", "independent_study_count", "study_independence",
    "cross_assay_consistency", "conflicting_functional_results",
    "highest_validated_assay_strength", "evidence_overlap", "evidence_conflicts",
)


def _parse_values(raw: Mapping[str, Any], fields: Sequence[str]) -> dict[str, EvidenceValue[Any]]:
    return {name: EvidenceValue.parse(raw.get(name), name) for name in fields if name in raw}


def validate_features(raw: PS3Evidence | Mapping[str, Any]) -> PS3Evidence:
    if isinstance(raw, PS3Evidence):
        return raw
    if not isinstance(raw, Mapping):
        raise SchemaValidationError("biological_features must be a mapping or PS3Evidence")
    context = DiseaseContext(_parse_values(raw, CONTEXT_FIELDS))
    assays_raw = raw.get("functional_assays", raw.get("assays"))
    if assays_raw is None:
        flat = {name: raw[name] for name in ASSAY_FIELDS if name in raw}
        assays_raw = [flat] if flat else []
    if not isinstance(assays_raw, (list, tuple)):
        raise SchemaValidationError("functional_assays must be a list")
    assays: list[FunctionalAssay] = []
    for index, item in enumerate(assays_raw):
        if not isinstance(item, Mapping):
            raise SchemaValidationError(f"functional_assays[{index}] must be an object")
        assays.append(FunctionalAssay(_parse_values(item, ASSAY_FIELDS), index=index))

    specs_raw = raw.get("vcep_ps3_specification", raw.get("vcep_specifications", ()))
    if specs_raw is None:
        specs_raw = ()
    if isinstance(specs_raw, Mapping):
        specs_raw = (specs_raw,)
    if not isinstance(specs_raw, (list, tuple)):
        raise SchemaValidationError("vcep_ps3_specification must be an object or list")
    specs = tuple(parse_vcep_specification(item) for item in specs_raw)
    return PS3Evidence(
        context=context,
        functional_assays=tuple(assays),
        vcep_specifications=specs,
        global_values=_parse_values(raw, GLOBAL_FIELDS),
    )


def _condition_value(
    name: str, assay: FunctionalAssay, context: DiseaseContext
) -> EvidenceValue[Any]:
    assay_value = assay.get(name)
    return assay_value if assay_value.state != EvidenceState.UNKNOWN else context.get(name)


def _condition_matches(actual: Any, expected: Any) -> bool:
    if isinstance(expected, Mapping) and "operator" in expected:
        try:
            left = float(actual)
            right = float(expected["value"])
        except (TypeError, ValueError, KeyError):
            return False
        operator = str(expected["operator"])
        return {
            "<": left < right,
            "<=": left <= right,
            ">": left > right,
            ">=": left >= right,
            "==": left == right,
        }.get(operator, False)
    if isinstance(expected, (list, tuple, set)):
        return any(_condition_matches(actual, item) for item in expected)
    if isinstance(actual, str) or isinstance(expected, str):
        return token(actual) == token(expected)
    return actual == expected


def _evaluate_vcep(
    evidence: PS3Evidence, assay: FunctionalAssay
) -> AssayEvaluation | None:
    if not evidence.vcep_specifications:
        return None
    gene = evidence.context.get("gene_id")
    disease = evidence.context.get("disease_id")
    potential: list[VCEPSpecification] = []
    scoped: list[VCEPSpecification] = []
    for spec in evidence.vcep_specifications:
        if gene.known and token(gene.value) != token(spec.gene_id):
            continue
        if not gene.known:
            potential.append(spec)
            continue
        if spec.disease_id:
            if not disease.known:
                potential.append(spec)
                continue
            if token(disease.value) != token(spec.disease_id):
                continue
        scoped.append(spec)
    if potential and not scoped:
        missing = ["gene_id"] if not gene.known else ["disease_id"]
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="VCEP_SCOPE",
            reason="A supplied higher-precedence VCEP specification may apply, but its scope cannot be resolved.",
            missing=tuple(missing), guideline=CLINGEN_SVI_2019,
        )
    if not scoped:
        return None

    matches: list[tuple[VCEPSpecification, VCEPRule, tuple[str, ...]]] = []
    unresolved: list[str] = []
    for spec in scoped:
        for rule in spec.rules:
            used: list[str] = ["gene_id"]
            if spec.disease_id:
                used.append("disease_id")
            failed = False
            unknown = False
            for name, expected in rule.when.items():
                actual = _condition_value(name, assay, evidence.context)
                if actual.state in {EvidenceState.UNKNOWN, EvidenceState.CONFLICTING}:
                    unresolved.append(name)
                    unknown = True
                    break
                used.append(name)
                if not _condition_matches(actual.value, expected):
                    failed = True
                    break
            if not failed and not unknown:
                matches.append((spec, rule, unique(used)))
    if matches:
        signatures = {(rule.status, rule.strength) for _, rule, _ in matches}
        if len(signatures) > 1:
            return AssayEvaluation(
                assay.identifier, DecisionStatus.UNDETERMINED,
                decision_path="VCEP_CONFLICT",
                reason="Multiple equally applicable supplied VCEP rules conflict.",
                missing=(), exclusions=("conflicting VCEP rules",), guideline=CLINGEN_SVI_2019,
            )
        spec, rule, used = matches[0]
        return AssayEvaluation(
            assay.identifier, rule.status, strength=rule.strength,
            decision_path=rule.decision_path, reason=rule.reason,
            features_used=used, guideline=vcep_reference(spec, rule),
            relevance_rank=4,
            validation_rank={Strength.NONE: 0, Strength.SUPPORTING: 1, Strength.MODERATE: 2,
                             Strength.STRONG: 3, Strength.VERY_STRONG: 4}[rule.strength],
        )
    if unresolved:
        return AssayEvaluation(
            assay.identifier, DecisionStatus.UNDETERMINED,
            decision_path="VCEP_RULE_SELECTION",
            reason="A scoped VCEP rule may apply, but a rule condition is unresolved.",
            missing=unique(unresolved), guideline=CLINGEN_SVI_2019,
        )
    if any(spec.exclusive for spec in scoped):
        return AssayEvaluation(
            assay.identifier, DecisionStatus.NOT_MET,
            decision_path="VCEP_ASSAY_REJECTED",
            reason="The assay/result does not meet any accepted rule in the applicable exclusive VCEP specification.",
            exclusions=("not accepted by applicable VCEP specification",), guideline=CLINGEN_SVI_2019,
        )
    return None


def select_applicable_guideline(
    evidence: PS3Evidence, assay: FunctionalAssay
) -> AssayEvaluation:
    vcep = _evaluate_vcep(evidence, assay)
    return vcep if vcep is not None else evaluate_assay(assay, evidence.context)


def _dominates(left: AssayEvaluation, right: AssayEvaluation) -> bool:
    return (
        left.relevance_rank >= right.relevance_rank
        and left.validation_rank >= right.validation_rank
        and (left.relevance_rank > right.relevance_rank or left.validation_rank > right.validation_rank)
    )


def _decision(
    status: DecisionStatus,
    *,
    strength: Strength = Strength.NONE,
    path: str,
    reason: str,
    results: tuple[AssayEvaluation, ...],
    selected: tuple[AssayEvaluation, ...] = (),
    missing: tuple[str, ...] = (),
    conflicts: tuple[str, ...] = (),
    bs3: str = "NOT_INDICATED",
) -> PS3Decision:
    used_results = selected or results
    guideline = next((item.guideline for item in used_results if item.guideline), CLINGEN_SVI_2019)
    return PS3Decision(
        status=status,
        strength=strength if status == DecisionStatus.ACTIVATED else Strength.NONE,
        decision_path=path,
        assays_used=tuple(item.assay_id for item in selected),
        features_used=unique([name for item in used_results for name in item.features_used]),
        missing_critical_features=unique(list(missing)),
        exclusions_evaluated=unique([name for item in results for name in item.exclusions]),
        conflicting_evidence=unique(list(conflicts)),
        guideline=guideline,
        reason=reason,
        assay_evaluations=results,
        bs3_relevance=bs3,
    )


def resolve_ps3_decision(
    evidence: PS3Evidence, results: tuple[AssayEvaluation, ...]
) -> PS3Decision:
    if not results:
        return _decision(
            DecisionStatus.UNDETERMINED, path="NO_FUNCTIONAL_ASSAY",
            reason="No functional assay evidence was supplied.", results=(),
            missing=("functional_assays",),
        )

    for name in ("conflicting_functional_results", "evidence_conflicts"):
        value = evidence.global_values.get(name, EvidenceValue())
        declared_true = value.value is True or (
            isinstance(value.value, str) and token(value.value) in {"TRUE", "YES", "PRESENT"}
        )
        if value.state == EvidenceState.CONFLICTING or (value.known and declared_true):
            return _decision(
                DecisionStatus.UNDETERMINED, path="MULTIPLE_ASSAYS_CONFLICT",
                reason="The supplied evidence declares an unresolved functional-evidence conflict.",
                results=results, conflicts=(name,),
            )
    consistency = evidence.global_values.get("cross_assay_consistency", EvidenceValue())
    declared_false = consistency.value is False or (
        isinstance(consistency.value, str) and token(consistency.value) in {"FALSE", "NO", "ABSENT"}
    )
    if consistency.known and declared_false:
        return _decision(
            DecisionStatus.UNDETERMINED, path="MULTIPLE_ASSAYS_CONFLICT",
            reason="Cross-assay consistency is explicitly false and no justified override is supplied.",
            results=results, conflicts=("cross_assay_consistency",),
        )

    pathogenic = tuple(item for item in results if item.status == DecisionStatus.ACTIVATED)
    normal = tuple(item for item in results if item.effect == FunctionalResult.NORMAL and item.bs3_relevant)
    if pathogenic and normal:
        if all(any(_dominates(p, n) for p in pathogenic) for n in normal):
            best = max(pathogenic, key=lambda item: (item.validation_rank, item.relevance_rank))
            return _decision(
                DecisionStatus.ACTIVATED, strength=best.strength,
                path="CONFLICT_RESOLVED_BY_SUPERIOR_ASSAY",
                reason="Abnormal and normal results conflict, but the pathogenic result is from an assay that is both more mechanism-relevant and at least as well validated (or vice versa), as permitted by ClinGen SVI.",
                results=results, selected=(best,), conflicts=tuple(item.assay_id for item in normal),
            )
        if all(any(_dominates(n, p) for n in normal) for p in pathogenic):
            best = max(normal, key=lambda item: (item.validation_rank, item.relevance_rank))
            return _decision(
                DecisionStatus.NOT_MET, path="CONFLICT_RESOLVED_BY_SUPERIOR_ASSAY",
                reason="A superior validated, mechanism-relevant assay shows normal function; PS3 is not met. Consider the result separately for BS3.",
                results=results, selected=(best,), conflicts=tuple(item.assay_id for item in pathogenic),
                bs3="POTENTIALLY_APPLICABLE",
            )
        return _decision(
            DecisionStatus.UNDETERMINED, path="MULTIPLE_ASSAYS_CONFLICT",
            reason="Credible abnormal and normal functional results cannot be resolved by mechanism relevance and validation; ClinGen SVI says not to use the conflicting evidence.",
            results=results, conflicts=tuple(item.assay_id for item in pathogenic + normal),
        )

    if pathogenic:
        best = max(pathogenic, key=lambda item: (item.validation_rank, item.relevance_rank))
        return _decision(
            DecisionStatus.ACTIVATED, strength=best.strength,
            path=best.decision_path,
            reason=best.reason + " Concordant studies are not stacked; the best-validated applicable assay determines strength.",
            results=results, selected=(best,),
        )
    if normal:
        best = max(normal, key=lambda item: (item.validation_rank, item.relevance_rank))
        return _decision(
            DecisionStatus.NOT_MET, path=best.decision_path,
            reason=best.reason, results=results, selected=(best,),
            bs3="POTENTIALLY_APPLICABLE",
        )

    unresolved = tuple(item for item in results if item.status == DecisionStatus.UNDETERMINED)
    if unresolved:
        return _decision(
            DecisionStatus.UNDETERMINED, path="INSUFFICIENT_OR_UNRESOLVED_EVIDENCE",
            reason="No assay completes a verified PS3 pathway; decision-critical evidence is missing or unresolved.",
            results=results,
            missing=unique([name for item in unresolved for name in item.missing]),
        )
    if all(item.status == DecisionStatus.NOT_APPLICABLE for item in results):
        return _decision(
            DecisionStatus.NOT_APPLICABLE, path="NO_APPLICABLE_ASSAY",
            reason="No supplied assay measures variant function relevant to the established disease mechanism.",
            results=results,
        )
    return _decision(
        DecisionStatus.NOT_MET, path="FUNCTIONAL_EVIDENCE_NOT_SUPPORTIVE",
        reason="The available assay evidence does not support PS3.", results=results,
    )


class PS3Engine:
    def evaluate(self, biological_features: PS3Evidence | Mapping[str, Any]) -> PS3Decision:
        evidence = validate_features(biological_features)
        results = tuple(
            select_applicable_guideline(evidence, assay)
            for assay in evidence.functional_assays
        )
        return resolve_ps3_decision(evidence, results)


def evaluate_ps3(biological_features: PS3Evidence | Mapping[str, Any]) -> PS3Decision:
    """Evaluate PS3 only; no annotation, network access, or ACMG-code combining."""
    return PS3Engine().evaluate(biological_features)
