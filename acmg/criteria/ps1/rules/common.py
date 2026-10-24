"""Shared non-scoring checks for PS1 comparison pathways."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, Optional

from ..evidence import Evidence, EvidenceState, stable_unique
from ..guidelines import ACMG_PS1
from ..models import (
    BiologicalFeatures,
    ClassificationEvidence,
    DecisionStatus,
    PS1Strength,
    ReferenceClassification,
    ReferenceVariant,
)
from ..output import GuidelineReference, RuleEvaluation


def label(reference_index: int, field: str) -> str:
    return f"references[{reference_index}].{field}"


def pairwise_evidence(
    evidence: BiologicalFeatures,
    reference: ReferenceVariant,
    reference_index: int,
    field: str,
) -> tuple[Evidence[object], str]:
    """Return a reference-specific pairwise assertion or its legacy default.

    Pairwise comparisons cannot safely be shared across several reference
    variants.  Nested reference assertions therefore override the top-level
    compatibility spelling.  The latter remains supported for a single
    comparison and for backwards-compatible callers.
    """

    reference_value = getattr(reference, field, None)
    if isinstance(reference_value, Evidence) and (
        field in reference.pairwise_assertions_supplied
        or not reference_value.is_unknown
    ):
        return reference_value, label(reference_index, field)
    return getattr(evidence, field), field


def exclusion(feature: str, outcome: str, detail: str = "") -> dict[str, str]:
    result = {"feature": feature, "outcome": outcome}
    if detail:
        result["detail"] = detail
    return result


def undetermined(
    path: str,
    guideline: GuidelineReference,
    reason: str,
    missing: Iterable[str],
    *,
    reference_index: Optional[int] = None,
    used: Iterable[str] = (),
    exclusions: Iterable[dict[str, str]] = (),
) -> RuleEvaluation:
    return RuleEvaluation(
        status=DecisionStatus.UNDETERMINED,
        strength=PS1Strength.NONE,
        decision_path=path,
        guideline=guideline,
        reason=reason,
        reference_index=reference_index,
        features_used=stable_unique(tuple(used)),
        missing_critical_features=stable_unique(tuple(missing)),
        exclusions_evaluated=tuple(exclusions),
    )


def not_met(
    path: str,
    guideline: GuidelineReference,
    reason: str,
    *,
    reference_index: Optional[int] = None,
    used: Iterable[str] = (),
    exclusions: Iterable[dict[str, str]] = (),
) -> RuleEvaluation:
    return RuleEvaluation(
        status=DecisionStatus.NOT_MET,
        strength=PS1Strength.NONE,
        decision_path=path,
        guideline=guideline,
        reason=reason,
        reference_index=reference_index,
        features_used=stable_unique(tuple(used)),
        exclusions_evaluated=tuple(exclusions),
    )


def not_applicable(
    path: str,
    guideline: GuidelineReference,
    reason: str,
    *,
    reference_index: Optional[int] = None,
    used: Iterable[str] = (),
) -> RuleEvaluation:
    return RuleEvaluation(
        status=DecisionStatus.NOT_APPLICABLE,
        strength=PS1Strength.NONE,
        decision_path=path,
        guideline=guideline,
        reason=reason,
        reference_index=reference_index,
        features_used=stable_unique(tuple(used)),
    )


def activated(
    strength: PS1Strength,
    path: str,
    guideline: GuidelineReference,
    reason: str,
    *,
    reference_index: int,
    used: Iterable[str] = (),
    exclusions: Iterable[dict[str, str]] = (),
) -> RuleEvaluation:
    return RuleEvaluation(
        status=DecisionStatus.ACTIVATED,
        strength=strength,
        decision_path=path,
        guideline=guideline,
        reason=reason,
        reference_index=reference_index,
        features_used=stable_unique(tuple(used)),
        exclusions_evaluated=tuple(exclusions),
    )


def normalized(value: str) -> str:
    """Whitespace/case-insensitive comparison of an already-normalized annotation."""

    return re.sub(r"\s+", "", value).upper()


def evidence_relation(left: Evidence[object], right: Evidence[object]) -> Optional[bool]:
    """Return equality only when both observations are known and non-conflicting."""

    if left.is_conflicting or right.is_conflicting:
        return None
    if not left.is_known or not right.is_known:
        return None
    left_value = left.value
    right_value = right.value
    if isinstance(left_value, str) and isinstance(right_value, str):
        return normalized(left_value) == normalized(right_value)
    return left_value == right_value


@dataclass(frozen=True)
class ReferencePreconditions:
    """Common eligibility checks shared by generic and VCEP pathways."""

    terminal: Optional[RuleEvaluation] = None
    used: tuple[str, ...] = ()
    exclusions: tuple[dict[str, str], ...] = ()


def _conflict_or_unknown(
    value: Evidence[object], field: str, *, reference_index: int, path: str, guideline: GuidelineReference,
    context: str,
) -> Optional[RuleEvaluation]:
    if value.is_conflicting:
        return undetermined(
            path, guideline, f"{context} is conflicting.", (field,),
            reference_index=reference_index, used=(field,),
        )
    if value.is_unknown:
        return undetermined(
            path, guideline, f"{context} is required but has not been established.", (field,),
            reference_index=reference_index,
        )
    return None


def _same_underlying_allele(
    evidence: BiologicalFeatures, reference: ReferenceVariant, reference_index: int
) -> tuple[Optional[bool], tuple[str, ...], Optional[str]]:
    """Derive whether the two annotations describe the same allele.

    Upstream normalized genomic alleles are preferred.  An explicit derived
    `underlying_allele_equivalence` assertion can resolve representation
    differences, but is checked against an exact equality when available.
    """

    explicit, explicit_label = pairwise_evidence(
        evidence, reference, reference_index, "underlying_allele_equivalence"
    )
    used: list[str] = []
    if explicit.is_conflicting:
        return None, (explicit_label,), "allele-equivalence evidence is conflicting"

    candidate = evidence.candidate
    ref = reference.annotations
    comparisons: list[bool] = []
    if candidate.genomic_allele.is_known and ref.genomic_allele.is_known:
        comparisons.append(
            normalized(candidate.genomic_allele.value or "")
            == normalized(ref.genomic_allele.value or "")
        )
        used += ("candidate_genomic_allele", label(reference_index, "genomic_allele"))
    if candidate.coding_hgvs.is_known and ref.coding_hgvs.is_known:
        comparisons.append(
            normalized(candidate.coding_hgvs.value or "")
            == normalized(ref.coding_hgvs.value or "")
        )
        used += ("candidate_coding_hgvs", label(reference_index, "coding_hgvs"))

    if comparisons and all(comparisons):
        exact: Optional[bool] = True
        conflict_reason: Optional[str] = None
    elif comparisons and not any(comparisons):
        exact = False
        conflict_reason = None
    elif comparisons:
        # Different normalised labels can be alternate representations, but a
        # mixture of equal and unequal supplied identities cannot establish a
        # distinct allele without an upstream curated reconciliation.
        exact = None
        conflict_reason = "genomic and coding allele identifiers disagree about whether the variants are the same allele"
    else:
        exact = None
        conflict_reason = None

    if explicit.is_known:
        asserted = explicit.value is True
        used.append(explicit_label)
        # Different textual representations can denote the same allele, so a
        # curated affirmative equivalence may resolve unequal identifiers.  By
        # contrast, any *identical* normalized genomic or coding identifier is
        # direct evidence of the same underlying allele and cannot be negated
        # by a Boolean assertion.  This also covers a mixed genome/coding
        # input, which is internally contradictory rather than evidence of a
        # distinct comparison allele.
        if any(comparisons) and asserted is False:
            return None, tuple(used), "allele-equivalence assertion contradicts an identical normalized allele identifier"
        return asserted, tuple(used), None
    return exact, tuple(used), conflict_reason


def common_reference_preconditions(
    evidence: BiologicalFeatures,
    reference: ReferenceVariant,
    reference_index: int,
    *,
    path: str,
    guideline: GuidelineReference = ACMG_PS1,
) -> ReferencePreconditions:
    """Check shared non-negotiable PS1 comparison constraints.

    The generic pathways require a disease-scoped, independently established
    P/LP comparison.  Classification *strength* (P versus LP) is deliberately
    left to each pathway or VCEP rule.
    """

    used: list[str] = []
    exclusions: list[dict[str, str]] = []

    if evidence.feature_conflict_status.is_true:
        return ReferencePreconditions(
            terminal=undetermined(
                path, guideline,
                "The feature set declares unresolved material biological conflicts.",
                ("feature_conflict_status",), reference_index=reference_index,
                used=("feature_conflict_status",),
            )
        )
    if evidence.feature_conflict_status.is_conflicting:
        return ReferencePreconditions(
            terminal=undetermined(
                path, guideline, "The global feature-conflict assertion is conflicting.",
                ("feature_conflict_status",), reference_index=reference_index,
            )
        )

    for field, value, description in (
        ("gene_id", evidence.gene_id, "Candidate gene context"),
        ("disease_id", evidence.disease_id, "Candidate disease context"),
    ):
        terminal = _conflict_or_unknown(
            value, field, reference_index=reference_index, path=path, guideline=guideline, context=description
        )
        if terminal is not None:
            return ReferencePreconditions(terminal=terminal)
        used.append(field)

    # A reference's explicit disease identifier can establish/match disease
    # context; otherwise the curated association Boolean is required.
    disease_match: Optional[bool] = None
    ref_disease = reference.disease_id
    if ref_disease.is_conflicting:
        return ReferencePreconditions(
            terminal=undetermined(path, guideline, "Reference disease identity is conflicting.",
                (label(reference_index, "disease_id"),), reference_index=reference_index)
        )
    if ref_disease.is_known:
        disease_match = normalized(ref_disease.value or "") == normalized(evidence.disease_id.value or "")
        used += ("disease_id", label(reference_index, "disease_id"))
    association = reference.disease_association
    if association.is_conflicting:
        return ReferencePreconditions(
            terminal=undetermined(path, guideline, "Reference disease-association evidence is conflicting.",
                (label(reference_index, "disease_association"),), reference_index=reference_index)
        )
    if association.is_known:
        used.append(label(reference_index, "disease_association"))
        if disease_match is False and association.value is True:
            return ReferencePreconditions(
                terminal=undetermined(path, guideline, "Reference disease identifiers and disease-association assertion conflict.",
                    (label(reference_index, "disease_association"),), reference_index=reference_index,
                    used=tuple(used))
            )
        disease_match = association.value is True
    if disease_match is False:
        exclusions.append(exclusion(label(reference_index, "disease_association"), "MISMATCH"))
        return ReferencePreconditions(
            terminal=not_met(path, guideline, "The reference is not established for the candidate disease context.",
                reference_index=reference_index, used=tuple(used), exclusions=tuple(exclusions))
        )
    if disease_match is None:
        return ReferencePreconditions(
            terminal=undetermined(path, guideline, "Reference disease compatibility has not been established.",
                (label(reference_index, "disease_association"),), reference_index=reference_index, used=tuple(used))
        )

    ref_gene = reference.gene_id
    if ref_gene.is_conflicting:
        return ReferencePreconditions(
            terminal=undetermined(path, guideline, "Reference gene identity is conflicting.",
                (label(reference_index, "gene_id"),), reference_index=reference_index)
        )
    if ref_gene.is_known:
        used.append(label(reference_index, "gene_id"))
        if normalized(ref_gene.value or "") != normalized(evidence.gene_id.value or ""):
            exclusions.append(exclusion(label(reference_index, "gene_id"), "MISMATCH"))
            return ReferencePreconditions(
                terminal=not_met(path, guideline, "Candidate and reference gene identifiers differ.",
                    reference_index=reference_index, used=tuple(used), exclusions=tuple(exclusions))
            )

    conflict = reference.classification_conflicts
    if conflict.is_conflicting or conflict.is_true:
        return ReferencePreconditions(
            terminal=undetermined(
                path, guideline, "Reference classification evidence is conflicting.",
                (label(reference_index, "classification_conflicts"),), reference_index=reference_index,
                used=(label(reference_index, "classification_conflicts"),),
            )
        )
    if conflict.is_false:
        used.append(label(reference_index, "classification_conflicts"))

    classification = reference.classification
    terminal = _conflict_or_unknown(
        classification, label(reference_index, "classification"), reference_index=reference_index,
        path=path, guideline=guideline, context="Reference classification"
    )
    if terminal is not None:
        return ReferencePreconditions(terminal=terminal)
    used.append(label(reference_index, "classification"))
    if classification.value not in {ReferenceClassification.PATHOGENIC, ReferenceClassification.LIKELY_PATHOGENIC}:
        exclusions.append(exclusion(label(reference_index, "classification"), "NOT_PATHOGENIC_OR_LIKELY_PATHOGENIC"))
        return ReferencePreconditions(
            terminal=not_met(path, guideline, "The reference is not classified Pathogenic or Likely Pathogenic.",
                reference_index=reference_index, used=tuple(used), exclusions=tuple(exclusions))
        )

    independence = reference.classification_independence
    terminal = _conflict_or_unknown(
        independence, label(reference_index, "classification_independence"), reference_index=reference_index,
        path=path, guideline=guideline, context="Reference-classification independence"
    )
    if terminal is not None:
        return ReferencePreconditions(terminal=terminal)
    used.append(label(reference_index, "classification_independence"))
    if independence.value is False:
        exclusions.append(exclusion(label(reference_index, "classification_independence"), "CIRCULAR_PS1_EVIDENCE"))
        return ReferencePreconditions(
            terminal=not_met(path, guideline, "The reference classification is not independent of PS1 comparison evidence.",
                reference_index=reference_index, used=tuple(used), exclusions=tuple(exclusions))
        )

    same_allele, allele_used, conflict_reason = _same_underlying_allele(
        evidence, reference, reference_index
    )
    used.extend(allele_used)
    if conflict_reason:
        return ReferencePreconditions(
            terminal=undetermined(path, guideline, conflict_reason, ("underlying_allele_equivalence",),
                reference_index=reference_index, used=tuple(used))
        )
    if same_allele is None:
        return ReferencePreconditions(
            terminal=undetermined(path, guideline, "Distinct candidate and reference alleles have not been established.",
                ("candidate_genomic_allele", label(reference_index, "genomic_allele")),
                reference_index=reference_index, used=tuple(used))
        )
    if same_allele:
        _, allele_label = pairwise_evidence(
            evidence, reference, reference_index, "underlying_allele_equivalence"
        )
        exclusions.append(exclusion(allele_label, "SAME_ALLELE"))
        return ReferencePreconditions(
            terminal=not_met(path, guideline, "The comparison is the same underlying allele, not an independent PS1 reference.",
                reference_index=reference_index, used=tuple(used), exclusions=tuple(exclusions))
        )

    return ReferencePreconditions(used=stable_unique(tuple(used)), exclusions=tuple(exclusions))


def require_reference_evidence_quality(
    evidence: BiologicalFeatures,
    reference: ReferenceVariant,
    reference_index: int,
    *,
    path: str,
    guideline: GuidelineReference,
    used: Iterable[str] = (),
    exclusions: Iterable[dict[str, str]] = (),
) -> Optional[RuleEvaluation]:
    """Apply the generic requirement for a documented established reference."""

    value = reference.classification_evidence
    field = label(reference_index, "classification_evidence")
    accumulated = tuple(used)
    if value.is_conflicting:
        return undetermined(path, guideline, "Reference classification-evidence adequacy is conflicting.", (field,),
            reference_index=reference_index, used=accumulated + (field,), exclusions=exclusions)
    if value.is_unknown:
        return undetermined(path, guideline, "Reference classification evidence has not been shown to be established.", (field,),
            reference_index=reference_index, used=accumulated, exclusions=exclusions)
    if value.value != ClassificationEvidence.ESTABLISHED:
        return not_met(path, guideline, "Reference classification evidence is explicitly insufficient for generic PS1.",
            reference_index=reference_index, used=accumulated + (field,),
            exclusions=tuple(exclusions) + (exclusion(field, "INSUFFICIENT"),))
    return None

