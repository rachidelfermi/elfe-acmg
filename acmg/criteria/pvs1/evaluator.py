"""Guideline-driven, evidence-adaptive PVS1 evaluation.

The implementation mirrors decision branches rather than accumulating scores.
Every helper either proves one branch, records a known exclusion, or asks for
the smallest still-decision-critical evidence item.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Iterable, Mapping, Optional, Union

from .models import (
    BiologicalFeatures,
    DecisionStatus,
    EvidenceKind,
    EvidenceStatus,
    EvidenceValue,
    GeneDiseaseValidity,
    GuidelineReference,
    LoFMechanismLevel,
    PVS1Decision,
    PVS1Strength,
    Provenance,
    ReadingFrame,
    SpliceOutcome,
    StructuralEvent,
    TandemDuplicationStatus,
    VariantConsequence,
    stable_unique,
)
from .specifications import (
    guideline_for_specific_rule,
    select_specific_rule,
    validate_specific_rule_declarations,
)


SVI_2018_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6185798/"
RNA_SPLICE_2023_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC10357475/"

# "PVS1" is the authentic ACMG/AMP criterion identifier.  Figure 1 itself
# does not assign stable IDs to individual branches; `decision_path` is the
# deterministic engine path and `source_locator` remains the citation.
GENERAL_SVI_GUIDELINE = GuidelineReference(
    authority="ClinGen Sequence Variant Interpretation Working Group",
    document="Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion",
    version="2018 (PMID: 30192042; DOI: 10.1002/humu.23626)",
    source_url=SVI_2018_URL,
    source_locator="Figure 1; Table 1",
    rule_id="PVS1",
)

RNA_SPLICE_GUIDELINE = GuidelineReference(
    authority="ClinGen SVI Splicing Subgroup",
    document="Using the ACMG/AMP framework to capture evidence related to predicted and observed impact on splicing",
    version="2023 (PMID: 37352859; DOI: 10.1016/j.ajhg.2023.06.002)",
    source_url=RNA_SPLICE_2023_URL,
    source_locator="PVS1(RNA) direct-splice guidance",
    rule_id="PVS1",
)


class Pathway:
    NULL_VARIANT = "NONSENSE_OR_FRAMESHIFT"
    CANONICAL_SPLICE = "CANONICAL_SPLICE"
    NONCANONICAL_SPLICE = "NONCANONICAL_SPLICE"
    EXON_DELETION = "EXON_DELETION"
    WHOLE_GENE_DELETION = "WHOLE_GENE_DELETION"
    DUPLICATION = "DUPLICATION"
    INITIATION_CODON = "INITIATION_CODON"


_STRENGTH_RANK = {
    PVS1Strength.SUPPORTING: 0,
    PVS1Strength.MODERATE: 1,
    PVS1Strength.STRONG: 2,
    PVS1Strength.VERY_STRONG: 3,
    PVS1Strength.NONE: -1,
}


@dataclass(frozen=True)
class _Candidate:
    status: DecisionStatus
    base_strength: PVS1Strength
    path: str
    locator: str
    reason: str
    used: tuple[str, ...] = ()
    missing: tuple[str, ...] = ()
    exclusions: tuple[Mapping[str, str], ...] = ()
    guideline: GuidelineReference = GENERAL_SVI_GUIDELINE
    guideline_conditions: Mapping[str, Any] = field(default_factory=dict)
    evidence_mode: Optional[str] = None


@dataclass(frozen=True)
class _MechanismAssessment:
    cap: Optional[PVS1Strength] = None
    used: tuple[str, ...] = ()
    candidate: Optional[_Candidate] = None


@dataclass(frozen=True)
class _RelevanceAssessment:
    kind: str  # OK, EXCLUDED, UNDETERMINED
    used: tuple[str, ...] = ()
    missing: tuple[str, ...] = ()
    exclusions: tuple[Mapping[str, str], ...] = ()
    reason: str = ""


@dataclass(frozen=True)
class _PathwayAssessment:
    pathway: Optional[str] = None
    candidate: Optional[_Candidate] = None


@dataclass(frozen=True)
class _SpliceOutcomeAssessment:
    outcomes: Optional[tuple[SpliceOutcome, ...]] = None
    used: tuple[str, ...] = ()
    candidate: Optional[_Candidate] = None


@dataclass(frozen=True)
class _RNAAssessment:
    """Readiness of a direct RNA profile for a PVS1(RNA) evaluation."""

    state: str  # ABSENT, UNUSABLE, READY
    used: tuple[str, ...] = ()
    candidate: Optional[_Candidate] = None
    reason: str = ""
    # A canonical +/-1/2 consequence can be evaluated from its independent
    # generic prediction evidence.  An incomplete positive RNA assertion must
    # still stop the non-canonical RNA-only route, but it must not turn RNA
    # metadata that the generic route does not require into a universal gate.
    canonical_generic_fallback: bool = False


def _weaker(left: PVS1Strength, right: PVS1Strength) -> PVS1Strength:
    return left if _STRENGTH_RANK[left] <= _STRENGTH_RANK[right] else right


def _unique_provenance(values: Iterable[Provenance]) -> tuple[Provenance, ...]:
    return tuple(dict.fromkeys(values))


def _exclusion(feature: str, outcome: str, detail: str = "") -> Mapping[str, str]:
    result: dict[str, str] = {"feature": feature, "outcome": outcome}
    if detail:
        result["detail"] = detail
    return result


def _undetermined(
    path: str,
    locator: str,
    reason: str,
    missing: Iterable[str],
    *,
    used: Iterable[str] = (),
    exclusions: Iterable[Mapping[str, str]] = (),
    guideline: GuidelineReference = GENERAL_SVI_GUIDELINE,
    evidence_mode: Optional[str] = None,
) -> _Candidate:
    return _Candidate(
        status=DecisionStatus.UNDETERMINED,
        base_strength=PVS1Strength.NONE,
        path=path,
        locator=locator,
        reason=reason,
        used=stable_unique(used),
        missing=stable_unique(missing),
        exclusions=tuple(exclusions),
        guideline=guideline,
        evidence_mode=evidence_mode,
    )


def _not_met(
    path: str,
    locator: str,
    reason: str,
    *,
    used: Iterable[str] = (),
    exclusions: Iterable[Mapping[str, str]] = (),
    guideline: GuidelineReference = GENERAL_SVI_GUIDELINE,
    evidence_mode: Optional[str] = None,
) -> _Candidate:
    return _Candidate(
        status=DecisionStatus.NOT_MET,
        base_strength=PVS1Strength.NONE,
        path=path,
        locator=locator,
        reason=reason,
        used=stable_unique(used),
        exclusions=tuple(exclusions),
        guideline=guideline,
        evidence_mode=evidence_mode,
    )


def _not_applicable(path: str, reason: str, *, used: Iterable[str] = ()) -> _Candidate:
    return _Candidate(
        status=DecisionStatus.NOT_APPLICABLE,
        base_strength=PVS1Strength.NONE,
        path=path,
        locator="Figure 1 — supported PVS1 variant classes",
        reason=reason,
        used=stable_unique(used),
    )


def _activated(
    strength: PVS1Strength,
    path: str,
    locator: str,
    reason: str,
    *,
    used: Iterable[str] = (),
    exclusions: Iterable[Mapping[str, str]] = (),
    guideline: GuidelineReference = GENERAL_SVI_GUIDELINE,
    evidence_mode: Optional[str] = None,
) -> _Candidate:
    return _Candidate(
        status=DecisionStatus.ACTIVATED,
        base_strength=strength,
        path=path,
        locator=locator,
        reason=reason,
        used=stable_unique(used),
        exclusions=tuple(exclusions),
        guideline=guideline,
        evidence_mode=evidence_mode,
    )


def _is_known_true(value: EvidenceValue[Any]) -> bool:
    return value.status == EvidenceStatus.KNOWN and value.value is True


def _is_known_false(value: EvidenceValue[Any]) -> bool:
    return value.status == EvidenceStatus.KNOWN and value.value is False


def _has_direct_rna_splice_assertion(evidence: BiologicalFeatures) -> bool:
    """Whether the caller supplied an observed/curated RNA assay assertion.

    This is deliberately only a routing predicate.  `_assess_rna_profile`
    still requires assay interpretability, gene-disease LoF sufficiency, and
    an observed allele-relevant product before it can activate PVS1(RNA).
    """

    rna = evidence.rna_splicing_evidence
    assay = evidence.rna_assay_interpretable
    return (
        rna.status == EvidenceStatus.KNOWN
        # A known non-interpretable assay cannot establish a splice route or
        # displace a still-possible canonical route. Unknown/conflicting assay
        # quality remains a direct-RNA route so the evaluator can request that
        # exact missing/reconciliation fact rather than pretending it is
        # negative evidence.
        and not _is_known_false(assay)
        and any(
            item.kind in {EvidenceKind.OBSERVED, EvidenceKind.CURATED}
            for item in rna.provenance
        )
    )


def _mechanism_assessment(evidence: BiologicalFeatures) -> _MechanismAssessment:
    """Apply the general SVI 2018 Table 1 mechanism eligibility/cap logic."""

    direct = evidence.lof_disease_mechanism
    locator = "Table 1 — Criteria for LoF disease mechanism"
    if direct.status == EvidenceStatus.CONFLICTING:
        return _MechanismAssessment(
            candidate=_undetermined(
                "SVI2018.ELIGIBILITY.LOF_CONFLICT",
                locator,
                "The supplied gene–disease LoF-mechanism conclusion is conflicting.",
                ("lof_disease_mechanism",),
                used=("lof_disease_mechanism",),
            )
        )
    if direct.status == EvidenceStatus.KNOWN:
        # A direct curated conclusion is normally the minimum sufficient
        # evidence. Do not silently ignore an explicitly supplied
        # incompatible gene-disease-validity assertion.
        validity = evidence.gene_disease_validity
        if validity.status == EvidenceStatus.CONFLICTING:
            return _MechanismAssessment(
                candidate=_undetermined(
                    "SVI2018.ELIGIBILITY.DIRECT_VALIDITY_CONFLICT",
                    locator,
                    "The curated LoF-mechanism conclusion cannot be reconciled with conflicting gene-disease validity evidence.",
                    ("gene_disease_validity",),
                    used=("lof_disease_mechanism", "gene_disease_validity"),
                )
            )
        if (
            direct.value != LoFMechanismLevel.NOT_ESTABLISHED
            and validity.status == EvidenceStatus.KNOWN
            and validity.value in {GeneDiseaseValidity.DISPUTED, GeneDiseaseValidity.REFUTED}
        ):
            return _MechanismAssessment(
                candidate=_undetermined(
                    "SVI2018.ELIGIBILITY.DIRECT_VALIDITY_CONFLICT",
                    locator,
                    "The supplied LoF-mechanism conclusion conflicts with a disputed or refuted gene-disease validity assertion.",
                    ("gene_disease_validity",),
                    used=("lof_disease_mechanism", "gene_disease_validity"),
                )
            )
        caps = {
            LoFMechanismLevel.ESTABLISHED: PVS1Strength.VERY_STRONG,
            LoFMechanismLevel.CAP_STRONG: PVS1Strength.STRONG,
            LoFMechanismLevel.CAP_MODERATE: PVS1Strength.MODERATE,
        }
        if direct.value == LoFMechanismLevel.NOT_ESTABLISHED:
            return _MechanismAssessment(
                candidate=_not_met(
                    "SVI2018.ELIGIBILITY.LOF_NOT_ESTABLISHED",
                    locator,
                    "PVS1 cannot be applied because LoF is known not to be an established mechanism for this gene–disease assertion.",
                    used=("lof_disease_mechanism",),
                )
            )
        return _MechanismAssessment(cap=caps[direct.value], used=("lof_disease_mechanism",))

    # A fully supplied Table 1 evidence set is a valid alternative to a
    # curated conclusion.  Partial raw evidence never becomes negative proof.
    raw_names = (
        "gene_disease_validity",
        "pathogenic_lof_variant_count",
        "pathogenic_lof_variant_fraction_percent",
        "pathogenic_lof_variants_across_multiple_exons",
        "null_model_recaptures_disease",
        "single_exon_gene",
    )
    raw = {name: evidence.evidence(name) for name in raw_names}
    validity = raw["gene_disease_validity"]
    count = raw["pathogenic_lof_variant_count"]
    fraction = raw["pathogenic_lof_variant_fraction_percent"]
    distributed = raw["pathogenic_lof_variants_across_multiple_exons"]
    mouse = raw["null_model_recaptures_disease"]
    single_exon = raw["single_exon_gene"]

    def distributed_or_single() -> Optional[bool]:
        if _is_known_true(distributed) or _is_known_true(single_exon):
            return True
        if _is_known_false(distributed) and _is_known_false(single_exon):
            return False
        return None

    distribution = distributed_or_single()
    distribution_used = (
        ("pathogenic_lof_variants_across_multiple_exons",)
        if _is_known_true(distributed)
        else ("single_exon_gene",)
        if _is_known_true(single_exon)
        else ()
    )
    used = tuple(
        name for name, item in raw.items() if item.status == EvidenceStatus.KNOWN
    )
    if validity.status == EvidenceStatus.KNOWN and validity.value in {
        GeneDiseaseValidity.STRONG,
        GeneDiseaseValidity.DEFINITIVE,
    }:
        if (
            count.status == EvidenceStatus.KNOWN
            and count.value >= 3
            and fraction.status == EvidenceStatus.KNOWN
            and fraction.value > 10
            and distribution is True
        ):
            return _MechanismAssessment(
                cap=PVS1Strength.VERY_STRONG,
                used=(
                    "gene_disease_validity",
                    "pathogenic_lof_variant_count",
                    "pathogenic_lof_variant_fraction_percent",
                )
                + distribution_used,
            )

    if validity.status == EvidenceStatus.KNOWN and validity.value in {
        GeneDiseaseValidity.MODERATE,
        GeneDiseaseValidity.STRONG,
        GeneDiseaseValidity.DEFINITIVE,
    }:
        prior_lof = (
            count.status == EvidenceStatus.KNOWN and count.value >= 2 and distribution is True
        )
        mouse_positive = _is_known_true(mouse)
        if prior_lof and mouse_positive:
            return _MechanismAssessment(
                cap=PVS1Strength.STRONG,
                used=(
                    "gene_disease_validity",
                    "pathogenic_lof_variant_count",
                )
                + distribution_used
                + ("null_model_recaptures_disease",),
            )
        if prior_lof:
            return _MechanismAssessment(
                cap=PVS1Strength.MODERATE,
                used=(
                    "gene_disease_validity",
                    "pathogenic_lof_variant_count",
                )
                + distribution_used,
            )
        if mouse_positive:
            return _MechanismAssessment(
                cap=PVS1Strength.MODERATE,
                used=("gene_disease_validity", "null_model_recaptures_disease"),
            )

    # One adjudicated LoF-mechanism field is the smallest sufficient missing
    # item, even if a caller has supplied incomplete Table 1 ingredients.
    return _MechanismAssessment(
        candidate=_undetermined(
            "SVI2018.ELIGIBILITY.LOF_UNKNOWN",
            locator,
            "The gene–disease LoF mechanism is not established by a curated conclusion or a complete applicable ClinGen SVI Table 1 evidence combination.",
            ("lof_disease_mechanism",),
            used=used,
        )
    )


def _relevance(evidence: BiologicalFeatures) -> _RelevanceAssessment:
    """Evaluate transcript/exon relevance and transcript-level rescue.

    Population exon tolerance belongs only to the Figure 1 role-unknown
    percentage branch. It must not suppress an NMD or critical-region branch.
    """

    used: list[str] = []
    exclusions: list[Mapping[str, str]] = []
    transcript = evidence.biologically_relevant_transcript
    exon = evidence.exon_biological_relevance
    policy = evidence.transcript_selection_policy

    # A supplied known negative is decisive unless the two direct relevance
    # assertions disagree. It should not be hidden by a conflict in a
    # modifier such as alternative-transcript rescue.
    direct_relevance = (
        ("biologically_relevant_transcript", transcript),
        ("exon_biological_relevance", exon),
    )
    known_direct = tuple(
        (name, value)
        for name, value in direct_relevance
        if value.status == EvidenceStatus.KNOWN
    )
    direct_positive = any(_is_known_true(value) for _, value in known_direct)
    direct_negative = any(_is_known_false(value) for _, value in known_direct)
    if direct_positive and direct_negative:
        return _RelevanceAssessment(
            "UNDETERMINED",
            used=tuple(name for name, _ in known_direct),
            missing=("biologically_relevant_transcript",),
            exclusions=(_exclusion("transcript_or_exon_relevance", "CONFLICTING"),),
            reason="The supplied transcript and exon relevance assertions disagree.",
        )
    if direct_negative:
        name = next(name for name, value in known_direct if _is_known_false(value))
        return _RelevanceAssessment(
            "EXCLUDED",
            used=(name,),
            exclusions=(_exclusion(name, "EXCLUDES_PVS1"),),
            reason="The affected sequence is known not to be present in a biologically relevant transcript/exon.",
        )

    for name, value in direct_relevance:
        if value.status == EvidenceStatus.CONFLICTING:
            return _RelevanceAssessment(
                "UNDETERMINED",
                used=(name,),
                missing=(name,),
                exclusions=(_exclusion(name, "CONFLICTING"),),
                reason="Biological transcript/exon relevance evidence is conflicting.",
            )

    rescue = evidence.alternative_transcript_rescue
    if _is_known_true(rescue):
        return _RelevanceAssessment(
            "EXCLUDED",
            used=("alternative_transcript_rescue",),
            exclusions=(_exclusion("alternative_transcript_rescue", "EXCLUDES_PVS1"),),
            reason="Credible alternative-transcript rescue shows that the affected sequence is not required in all biologically relevant transcripts.",
        )
    if rescue.status == EvidenceStatus.CONFLICTING:
        return _RelevanceAssessment(
            "UNDETERMINED",
            used=("alternative_transcript_rescue",),
            missing=("alternative_transcript_rescue",),
            exclusions=(_exclusion("alternative_transcript_rescue", "CONFLICTING"),),
            reason="Alternative-transcript rescue evidence is conflicting.",
        )
    if rescue.status == EvidenceStatus.KNOWN:
        used.append("alternative_transcript_rescue")
        exclusions.append(_exclusion("alternative_transcript_rescue", "NOT_PRESENT"))

    if _is_known_true(transcript):
        used.append("biologically_relevant_transcript")
    elif _is_known_true(exon):
        used.append("exon_biological_relevance")
    elif policy.status == EvidenceStatus.KNOWN:
        used.append("transcript_selection_policy")
    elif policy.status == EvidenceStatus.CONFLICTING:
        return _RelevanceAssessment(
            "UNDETERMINED",
            used=tuple(used) + ("transcript_selection_policy",),
            missing=("transcript_selection_policy",),
            exclusions=tuple(exclusions) + (_exclusion("transcript_selection_policy", "CONFLICTING"),),
            reason="Transcript-selection policy evidence is conflicting.",
        )
    else:
        return _RelevanceAssessment(
            "UNDETERMINED",
            used=tuple(used),
            missing=("biologically_relevant_transcript",),
            exclusions=tuple(exclusions) + (_exclusion("transcript_or_exon_relevance", "UNKNOWN"),),
            reason="Biological relevance of the affected transcript/exon has not been established.",
        )
    exclusions.append(_exclusion("transcript_or_exon_relevance", "PRESENT"))
    return _RelevanceAssessment("OK", tuple(used), exclusions=tuple(exclusions))


def _relevance_candidate(
    relevance: _RelevanceAssessment,
    path: str,
    locator: str,
) -> Optional[_Candidate]:
    if relevance.kind == "OK":
        return None
    if relevance.kind == "EXCLUDED":
        return _not_met(path, locator, relevance.reason, used=relevance.used, exclusions=relevance.exclusions)
    return _undetermined(
        path,
        locator,
        relevance.reason,
        relevance.missing,
        used=relevance.used,
        exclusions=relevance.exclusions,
    )


def _non_nmd_ladder(
    evidence: BiologicalFeatures,
    path_prefix: str,
    label: str,
    *,
    guideline: GuidelineReference = GENERAL_SVI_GUIDELINE,
    evidence_mode: Optional[str] = None,
) -> _Candidate:
    """The SVI Figure 1 escape-NMD/in-frame percentage-strength ladder."""

    locator = "Figure 1 — critical region or percentage-loss branch"
    relevance = _relevance(evidence)
    blocked = _relevance_candidate(relevance, f"{path_prefix}.RELEVANCE", locator)
    if blocked:
        return replace(blocked, guideline=guideline, evidence_mode=evidence_mode)

    critical = evidence.critical_protein_region
    if critical.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{path_prefix}.CRITICAL_REGION_CONFLICT",
            locator,
            "Evidence about whether the altered region is critical to protein function is conflicting.",
            ("critical_protein_region",),
            used=relevance.used + ("critical_protein_region",),
            exclusions=relevance.exclusions,
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    if _is_known_true(critical):
        return _activated(
            PVS1Strength.STRONG,
            f"{path_prefix}.CRITICAL_REGION",
            "Figure 1 — truncated/altered region critical to protein function",
            f"{label} affects a region established as critical to protein function.",
            used=relevance.used + ("critical_protein_region",),
            exclusions=relevance.exclusions,
            guideline=guideline,
            evidence_mode=evidence_mode,
        )

    # The role-unknown route is a real independent Figure 1 path; lack of a
    # critical-region assertion does not become a negative value or block it.
    tolerance = evidence.population_evidence_for_exon_lof_tolerance
    if tolerance.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{path_prefix}.TOLERANCE_CONFLICT",
            locator,
            "Population evidence about LoF tolerance is conflicting for the role-unknown percentage branch.",
            ("population_evidence_for_exon_lof_tolerance",),
            used=relevance.used + ("population_evidence_for_exon_lof_tolerance",),
            exclusions=relevance.exclusions
            + (_exclusion("population_evidence_for_exon_lof_tolerance", "CONFLICTING"),),
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    if _is_known_true(tolerance):
        return _not_met(
            f"{path_prefix}.TOLERANT_EXON",
            locator,
            "The role-unknown percentage branch cannot apply because LoF variants in the affected exon are frequent in the general population.",
            used=relevance.used + ("population_evidence_for_exon_lof_tolerance",),
            exclusions=relevance.exclusions
            + (_exclusion("population_evidence_for_exon_lof_tolerance", "EXCLUDES_PVS1"),),
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    if tolerance.status == EvidenceStatus.UNKNOWN:
        return _undetermined(
            f"{path_prefix}.PERCENTAGE_LOSS",
            locator,
            "The role-unknown percentage branch requires exon LoF-tolerance evidence.",
            ("population_evidence_for_exon_lof_tolerance",),
            used=relevance.used,
            exclusions=relevance.exclusions,
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    percent = evidence.protein_truncation_percent
    if percent.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{path_prefix}.PERCENTAGE_LOSS_CONFLICT",
            locator,
            "Protein loss/alteration proportion evidence is conflicting.",
            ("protein_truncation_percent",),
            used=relevance.used + ("protein_truncation_percent",),
            exclusions=relevance.exclusions,
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    if percent.status == EvidenceStatus.UNKNOWN:
        return _undetermined(
            f"{path_prefix}.PERCENTAGE_LOSS",
            locator,
            "The role-unknown percentage branch requires the proportion of protein removed or altered.",
            ("protein_truncation_percent",),
            used=relevance.used,
            exclusions=relevance.exclusions,
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    if percent.value > 10:
        strength = PVS1Strength.STRONG
        branch = "GREATER_THAN_10_PERCENT"
    elif percent.value < 10:
        strength = PVS1Strength.MODERATE
        branch = "LESS_THAN_10_PERCENT"
    else:
        return _undetermined(
            f"{path_prefix}.EXACTLY_10_PERCENT",
            locator,
            "ClinGen SVI Figure 1 specifies >10% and <10% protein loss, but not exactly 10%; a scoped VCEP rule or expert adjudication is required.",
            ("protein_truncation_percent",),
            used=relevance.used + ("protein_truncation_percent",),
            exclusions=relevance.exclusions,
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    return _activated(
        strength,
        f"{path_prefix}.{branch}",
        f"Figure 1 — {branch.replace('_', ' ').lower()} protein altered",
        f"{label} escapes NMD or preserves the reading frame; the affected exon is relevant and non-tolerant, and the altered protein proportion meets the Figure 1 {branch.replace('_', ' ').lower()} branch.",
        used=relevance.used + ("population_evidence_for_exon_lof_tolerance", "protein_truncation_percent"),
        exclusions=relevance.exclusions,
        guideline=guideline,
        evidence_mode=evidence_mode,
    )


def _nmd_or_escape(
    evidence: BiologicalFeatures,
    path_prefix: str,
    label: str,
    *,
    guideline: GuidelineReference = GENERAL_SVI_GUIDELINE,
    evidence_mode: Optional[str] = None,
) -> _Candidate:
    """Evaluate null branches without allowing NMD uncertainty to hide a veto."""

    nmd = evidence.nmd_prediction
    locator = "Figure 1 — predicted NMD branch"
    # Known transcript rescue or lack of biological relevance is decisive
    # independently of the eventual NMD strength branch. Evaluate it before
    # an unknown NMD result so a stronger-looking early branch cannot mask a
    # known exclusion.
    relevance = _relevance(evidence)
    blocked = _relevance_candidate(relevance, f"{path_prefix}.RELEVANCE", locator)
    if blocked:
        return replace(blocked, guideline=guideline, evidence_mode=evidence_mode)
    if nmd.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{path_prefix}.NMD_CONFLICT",
            locator,
            "NMD prediction evidence is conflicting.",
            ("nmd_prediction",),
            used=("nmd_prediction",),
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    if nmd.status == EvidenceStatus.UNKNOWN:
        return _undetermined(
            f"{path_prefix}.NMD_UNKNOWN",
            locator,
            "NMD prediction is required to distinguish the Very Strong branch from the NMD-escape branches; it cannot be inferred from a premature stop alone.",
            ("nmd_prediction",),
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    if nmd.value is False:
        result = _non_nmd_ladder(
            evidence,
            f"{path_prefix}.NMD_ESCAPE",
            label,
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
        return replace(result, used=stable_unique(("nmd_prediction",) + result.used))

    return _activated(
        PVS1Strength.VERY_STRONG,
        f"{path_prefix}.NMD.RELEVANT",
        "Figure 1 — predicted NMD; exon present in biologically relevant transcript(s)",
        f"{label} is predicted to undergo NMD and affects a biologically relevant transcript/exon.",
        used=("nmd_prediction",) + relevance.used,
        exclusions=relevance.exclusions,
        guideline=guideline,
        evidence_mode=evidence_mode,
    )


def _evaluate_null_variant(evidence: BiologicalFeatures) -> _Candidate:
    return _nmd_or_escape(
        evidence,
        "SVI2018.NF",
        "The nonsense/frameshift variant",
    )


def _evaluate_in_frame_splice(
    evidence: BiologicalFeatures,
    *,
    guideline: GuidelineReference,
    evidence_mode: Optional[str],
    prefix: str,
) -> _Candidate:
    # Walker et al. 2023 describes a distinct PVS1(RNA) Very Strong route
    # when an RNA-confirmed in-frame event removes undisputed clinically
    # relevant residues. Do not infer that higher bar from the generic
    # critical-region feature alone.
    rna_critical = evidence.rna_in_frame_loss_of_clinically_relevant_residues
    if evidence_mode == "RNA" and rna_critical.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{prefix}.RNA_CLINICALLY_RELEVANT_RESIDUES_CONFLICT",
            "PVS1(RNA) in-frame clinically relevant residues branch",
            "Evidence about whether the RNA-confirmed in-frame loss removes clinically relevant residues is conflicting.",
            ("rna_in_frame_loss_of_clinically_relevant_residues",),
            used=("rna_in_frame_loss_of_clinically_relevant_residues",),
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    if evidence_mode == "RNA" and _is_known_true(rna_critical):
        relevance = _relevance(evidence)
        blocked = _relevance_candidate(
            relevance,
            f"{prefix}.RNA_RELEVANCE",
            "PVS1(RNA) in-frame clinically relevant residues branch",
        )
        if blocked is not None:
            return replace(blocked, guideline=guideline, evidence_mode=evidence_mode)
        return _activated(
            PVS1Strength.VERY_STRONG,
            f"{prefix}.RNA_CLINICALLY_RELEVANT_RESIDUES",
            "PVS1(RNA) — in-frame loss of undisputed clinically relevant residues",
            "A usable RNA assay establishes an in-frame splice outcome that removes undisputed clinically relevant residues.",
            used=relevance.used
            + ("rna_in_frame_loss_of_clinically_relevant_residues",),
            exclusions=relevance.exclusions,
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    return _non_nmd_ladder(
        evidence,
        prefix,
        "The in-frame splice outcome",
        guideline=guideline,
        evidence_mode=evidence_mode,
    )


def _resolve_splice_outcomes(
    evidence: BiologicalFeatures,
    *,
    require_observed: bool,
    prefix: str,
    locator: str,
    guideline: GuidelineReference,
    evidence_mode: Optional[str],
) -> _SpliceOutcomeAssessment:
    """Obtain a resolved splice profile without assuming a splice consequence.

    Generic canonical evaluation may derive a frame outcome only from an
    explicitly supplied exon-skipping/cryptic-site assertion *and* an
    explicitly supplied frame result. RNA evaluation never substitutes a
    prediction for an observed splice profile.
    """

    outcomes = evidence.splice_outcomes
    if outcomes.status == EvidenceStatus.CONFLICTING:
        return _SpliceOutcomeAssessment(
            candidate=_undetermined(
                f"{prefix}.OUTCOME_CONFLICT",
                locator,
                "The expected or observed splice outcome is conflicting.",
                ("splice_effect",),
                used=("splice_outcomes",),
                guideline=guideline,
                evidence_mode=evidence_mode,
            )
        )
    if outcomes.status == EvidenceStatus.KNOWN:
        if require_observed and not any(
            item.kind in {EvidenceKind.OBSERVED, EvidenceKind.CURATED}
            for item in outcomes.provenance
        ):
            return _SpliceOutcomeAssessment(
                candidate=_undetermined(
                    f"{prefix}.OUTCOME_NOT_OBSERVED",
                    locator,
                    "PVS1(RNA) requires observed or curated RNA splice-profile provenance, not a prediction-only or unprovenanced outcome.",
                    ("splice_effect",),
                    used=("splice_outcomes",),
                    guideline=guideline,
                    evidence_mode=evidence_mode,
                )
            )
        return _SpliceOutcomeAssessment(outcomes=outcomes.value, used=("splice_outcomes",))

    if require_observed:
        return _SpliceOutcomeAssessment(
            candidate=_undetermined(
                f"{prefix}.OUTCOME_UNKNOWN",
                locator,
                "PVS1(RNA) requires a resolved allele-relevant RNA splice profile.",
                ("splice_effect",),
                guideline=guideline,
                evidence_mode=evidence_mode,
            )
        )

    signals = (
        ("predicted_exon_skipping", evidence.predicted_exon_skipping),
        ("cryptic_splice_site_usage", evidence.cryptic_splice_site_usage),
    )
    positive_signals = tuple(name for name, value in signals if _is_known_true(value))
    if not positive_signals:
        conflicting = next(
            (name for name, value in signals if value.status == EvidenceStatus.CONFLICTING),
            None,
        )
        if conflicting is not None:
            return _SpliceOutcomeAssessment(
                candidate=_undetermined(
                    f"{prefix}.OUTCOME_PREDICTION_CONFLICT",
                    locator,
                    "Splice-effect prediction evidence is conflicting.",
                    (conflicting,),
                    used=(conflicting,),
                    guideline=guideline,
                    evidence_mode=evidence_mode,
                )
            )
        return _SpliceOutcomeAssessment(
            candidate=_undetermined(
                f"{prefix}.OUTCOME_UNKNOWN",
                locator,
                "Canonical splice position alone does not establish exon skipping, cryptic-site use, reading-frame effect, or NMD.",
                ("splice_effect",),
                guideline=guideline,
                evidence_mode=evidence_mode,
            )
        )

    frame = evidence.reading_frame_disruption
    if frame.status == EvidenceStatus.CONFLICTING:
        return _SpliceOutcomeAssessment(
            candidate=_undetermined(
                f"{prefix}.FRAME_CONFLICT",
                locator,
                "The explicitly supplied reading-frame consequence for the splice prediction is conflicting.",
                ("reading_frame_disruption",),
                used=positive_signals + ("reading_frame_disruption",),
                guideline=guideline,
                evidence_mode=evidence_mode,
            )
        )
    if frame.status == EvidenceStatus.UNKNOWN:
        return _SpliceOutcomeAssessment(
            candidate=_undetermined(
                f"{prefix}.FRAME_UNKNOWN",
                locator,
                "An explicit exon-skipping/cryptic-site prediction still requires a supplied reading-frame consequence.",
                ("reading_frame_disruption",),
                used=positive_signals,
                guideline=guideline,
                evidence_mode=evidence_mode,
            )
        )
    outcome = SpliceOutcome.OUT_OF_FRAME if frame.value else SpliceOutcome.IN_FRAME
    return _SpliceOutcomeAssessment(
        outcomes=(outcome,),
        used=positive_signals + ("reading_frame_disruption",),
    )


def _assess_rna_profile(evidence: BiologicalFeatures) -> _RNAAssessment:
    """Assess whether direct RNA evidence is usable for a PVS1(RNA) branch."""

    prefix = "SVI2023.RNA_SPLICE"
    locator = RNA_SPLICE_GUIDELINE.source_locator
    rna = evidence.rna_splicing_evidence
    if rna.status == EvidenceStatus.CONFLICTING:
        return _RNAAssessment(
            "UNUSABLE",
            candidate=_undetermined(
                f"{prefix}.EVIDENCE_CONFLICT",
                locator,
                "RNA splicing evidence is conflicting.",
                ("rna_splicing_evidence",),
                used=("rna_splicing_evidence",),
                guideline=RNA_SPLICE_GUIDELINE,
                evidence_mode="RNA",
            ),
        )
    if rna.status == EvidenceStatus.UNKNOWN:
        return _RNAAssessment("ABSENT", reason="No direct RNA assay conclusion was supplied.")
    if not any(
        item.kind in {EvidenceKind.OBSERVED, EvidenceKind.CURATED}
        for item in rna.provenance
    ):
        return _RNAAssessment(
            "UNUSABLE",
            used=("rna_splicing_evidence",),
            reason="The supplied RNA assertion lacks observed or curated RNA-assay provenance.",
        )

    interpretable = evidence.rna_assay_interpretable
    if interpretable.status == EvidenceStatus.CONFLICTING:
        return _RNAAssessment(
            "UNUSABLE",
            candidate=_undetermined(
                f"{prefix}.ASSAY_INTERPRETABILITY_CONFLICT",
                locator,
                "RNA assay interpretability evidence is conflicting.",
                ("rna_assay_interpretable",),
                used=("rna_splicing_evidence", "rna_assay_interpretable"),
                guideline=RNA_SPLICE_GUIDELINE,
                evidence_mode="RNA",
            ),
        )
    if interpretable.status == EvidenceStatus.UNKNOWN:
        return _RNAAssessment(
            "UNUSABLE",
            candidate=_undetermined(
                f"{prefix}.ASSAY_INTERPRETABILITY_UNKNOWN",
                locator,
                "A positive RNA assertion requires an assessment that the assay is interpretable for this gene-disease context.",
                ("rna_assay_interpretable",),
                used=("rna_splicing_evidence",),
                guideline=RNA_SPLICE_GUIDELINE,
                evidence_mode="RNA",
            ),
            canonical_generic_fallback=True,
        )
    if interpretable.value is False:
        return _RNAAssessment(
            "UNUSABLE",
            used=("rna_splicing_evidence", "rna_assay_interpretable"),
            reason="The supplied RNA assay is known not to be interpretable for PVS1(RNA).",
        )

    # The high-level RNA conclusion and a resolved, allele-relevant RNA
    # isoform profile are two assertions about the same experiment.  Neither
    # is allowed to silently override the other.  Only an observed/curated
    # profile participates here; a predictor-labelled splice outcome will be
    # rejected by `_resolve_splice_outcomes` below rather than manufacturing a
    # cross-evidence conflict.
    outcomes = evidence.splice_outcomes
    observed_profile = (
        outcomes.status == EvidenceStatus.KNOWN
        and any(
            item.kind in {EvidenceKind.OBSERVED, EvidenceKind.CURATED}
            for item in outcomes.provenance
        )
    )
    if observed_profile:
        assert outcomes.value is not None
        profile_has_aberrant_outcome = any(
            outcome != SpliceOutcome.NO_ABERRANT_SPLICING
            for outcome in outcomes.value
        )
        if rna.value != profile_has_aberrant_outcome:
            return _RNAAssessment(
                "UNUSABLE",
                candidate=_undetermined(
                    f"{prefix}.PROFILE_CONFLICT",
                    locator,
                    "The direct RNA assay conclusion conflicts with its resolved observed splice profile.",
                    ("rna_splicing_evidence", "splice_effect"),
                    used=(
                        "rna_splicing_evidence",
                        "rna_assay_interpretable",
                        "splice_outcomes",
                    ),
                    guideline=RNA_SPLICE_GUIDELINE,
                    evidence_mode="RNA",
                ),
            )
    if rna.value is False:
        return _RNAAssessment(
            "UNUSABLE",
            candidate=_not_met(
                f"{prefix}.NO_ABERRANT_SPLICING",
                locator,
                "An interpretable direct RNA conclusion indicates no aberrant splicing suitable for PVS1.",
                used=("rna_splicing_evidence", "rna_assay_interpretable"),
                guideline=RNA_SPLICE_GUIDELINE,
                evidence_mode="RNA",
            ),
        )
    sufficient = evidence.rna_aberrant_splicing_sufficient_for_lof
    used = (
        "rna_splicing_evidence",
        "rna_assay_interpretable",
        "rna_aberrant_splicing_sufficient_for_lof",
    )
    if sufficient.status == EvidenceStatus.CONFLICTING:
        return _RNAAssessment(
            "UNUSABLE",
            candidate=_undetermined(
                f"{prefix}.PROFILE_SUFFICIENCY_CONFLICT",
                locator,
                "The clinical sufficiency of the aberrant RNA profile is conflicting.",
                ("rna_aberrant_splicing_sufficient_for_lof",),
                used=used,
                guideline=RNA_SPLICE_GUIDELINE,
                evidence_mode="RNA",
            ),
        )
    if sufficient.status == EvidenceStatus.UNKNOWN:
        return _RNAAssessment(
            "UNUSABLE",
            candidate=_undetermined(
                f"{prefix}.PROFILE_SUFFICIENCY_UNKNOWN",
                locator,
                "A usable RNA assay still requires a curated assessment that the aberrant profile is sufficient for LoF in this gene-disease context.",
                ("rna_aberrant_splicing_sufficient_for_lof",),
                used=used[:-1],
                guideline=RNA_SPLICE_GUIDELINE,
                evidence_mode="RNA",
            ),
            canonical_generic_fallback=True,
        )
    if sufficient.value is False:
        return _RNAAssessment(
            "UNUSABLE",
            candidate=_not_met(
                f"{prefix}.PROFILE_NOT_SUFFICIENT",
                locator,
                "The observed aberrant RNA profile is known not to be sufficient for PVS1 in this gene-disease context.",
                used=used,
                guideline=RNA_SPLICE_GUIDELINE,
                evidence_mode="RNA",
            ),
        )

    # A positive direct-RNA assertion needs an allele-relevant observed
    # product before it can drive PVS1(RNA).  The absent/predicted-only product
    # is incomplete rather than contradictory, so a canonical +/-1/2 variant
    # may still use independent generic prediction evidence. A conflicting
    # observed product remains material and cannot be bypassed.
    if outcomes.status == EvidenceStatus.CONFLICTING:
        return _RNAAssessment(
            "UNUSABLE",
            candidate=_undetermined(
                f"{prefix}.OUTCOME_CONFLICT",
                locator,
                "The observed RNA splice-product evidence is conflicting.",
                ("splice_effect",),
                used=used + ("splice_outcomes",),
                guideline=RNA_SPLICE_GUIDELINE,
                evidence_mode="RNA",
            ),
        )
    if outcomes.status == EvidenceStatus.UNKNOWN:
        return _RNAAssessment(
            "UNUSABLE",
            candidate=_undetermined(
                f"{prefix}.OUTCOME_UNKNOWN",
                locator,
                "PVS1(RNA) requires a resolved allele-relevant RNA splice profile.",
                ("splice_effect",),
                used=used,
                guideline=RNA_SPLICE_GUIDELINE,
                evidence_mode="RNA",
            ),
            canonical_generic_fallback=True,
        )
    if not any(
        item.kind in {EvidenceKind.OBSERVED, EvidenceKind.CURATED}
        for item in outcomes.provenance
    ):
        return _RNAAssessment(
            "UNUSABLE",
            candidate=_undetermined(
                f"{prefix}.OUTCOME_NOT_OBSERVED",
                locator,
                "PVS1(RNA) requires observed or curated RNA splice-profile provenance, not a prediction-only or unprovenanced outcome.",
                ("splice_effect",),
                used=used + ("splice_outcomes",),
                guideline=RNA_SPLICE_GUIDELINE,
                evidence_mode="RNA",
            ),
            canonical_generic_fallback=True,
        )
    return _RNAAssessment("READY", used=used)


def _evaluate_splice_legacy(
    evidence: BiologicalFeatures,
    *,
    canonical: bool,
    guideline: GuidelineReference = GENERAL_SVI_GUIDELINE,
    evidence_mode: Optional[str] = None,
) -> _Candidate:
    """Evaluate resolved splice outcomes, choosing the lowest credible strength."""

    prefix = "SVI2018.CANONICAL_SPLICE" if canonical else "SVI2023.RNA_SPLICE"
    locator = "Figure 1 — GT–AG ±1,2 splice-site branch" if canonical else guideline.source_locator
    rescue = evidence.splice_rescue_evidence
    if rescue.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{prefix}.RESCUE_CONFLICT",
            locator,
            "Splice-rescue evidence is conflicting.",
            ("splice_rescue_evidence",),
            used=("splice_rescue_evidence",),
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    if _is_known_true(rescue):
        return _not_met(
            f"{prefix}.RESCUED",
            locator,
            "Credible splice-rescue evidence prevents application of PVS1.",
            used=("splice_rescue_evidence",),
            exclusions=(_exclusion("splice_rescue_evidence", "EXCLUDES_PVS1"),),
            guideline=guideline,
            evidence_mode=evidence_mode,
        )

    outcomes = evidence.splice_outcomes
    if outcomes.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{prefix}.OUTCOME_CONFLICT",
            locator,
            "The expected or observed splice outcome is conflicting.",
            ("splice_effect",),
            used=("splice_outcomes",),
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    if outcomes.status == EvidenceStatus.UNKNOWN:
        return _undetermined(
            f"{prefix}.OUTCOME_UNKNOWN",
            locator,
            "Canonical splice position alone does not establish exon skipping, cryptic-site use, reading-frame effect, or NMD.",
            ("splice_effect",),
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    if outcomes.value == (SpliceOutcome.NO_ABERRANT_SPLICING,):
        return _not_met(
            f"{prefix}.NO_ABERRANT_SPLICING",
            locator,
            "The resolved splice evidence indicates no aberrant splice consequence suitable for PVS1.",
            used=("splice_outcomes",),
            guideline=guideline,
            evidence_mode=evidence_mode,
        )
    damaging = tuple(
        outcome
        for outcome in outcomes.value
        if outcome != SpliceOutcome.NO_ABERRANT_SPLICING
    )
    if len(damaging) != len(outcomes.value) and not _is_known_false(rescue):
        return _undetermined(
            f"{prefix}.MIXED_NORMAL_AND_ABERRANT",
            locator,
            "Both normal and aberrant splice outcomes are credible; adequacy of normal-transcript rescue has not been resolved.",
            ("splice_rescue_evidence",),
            used=("splice_outcomes",),
            guideline=guideline,
            evidence_mode=evidence_mode,
        )

    candidates: list[_Candidate] = []
    for outcome in damaging:
        if outcome == SpliceOutcome.OUT_OF_FRAME:
            candidates.append(
                _nmd_or_escape(
                    evidence,
                    f"{prefix}.OUT_OF_FRAME",
                    "The out-of-frame splice outcome",
                    guideline=guideline,
                    evidence_mode=evidence_mode,
                )
            )
        elif outcome == SpliceOutcome.IN_FRAME:
            candidates.append(
                _evaluate_in_frame_splice(
                    evidence,
                    guideline=guideline,
                    evidence_mode=evidence_mode,
                    prefix=f"{prefix}.IN_FRAME",
                )
            )

    unresolved = next((candidate for candidate in candidates if candidate.status == DecisionStatus.UNDETERMINED), None)
    if unresolved is not None:
        return replace(
            unresolved,
            used=stable_unique(("splice_outcomes",) + unresolved.used),
        )
    not_met = next((candidate for candidate in candidates if candidate.status == DecisionStatus.NOT_MET), None)
    if not_met is not None:
        return replace(not_met, used=stable_unique(("splice_outcomes",) + not_met.used))
    weakest = min(candidates, key=lambda candidate: _STRENGTH_RANK[candidate.base_strength])
    if len(candidates) == 1:
        return replace(weakest, used=stable_unique(("splice_outcomes",) + weakest.used))
    return replace(
        weakest,
        path=f"{prefix}.MULTIPLE_OUTCOMES.LOWEST_STRENGTH",
        locator="Figure 1 — apply lowest strength across credible splice outcomes",
        reason="Multiple credible splice outcomes were evaluated; the lowest applicable PVS1 strength was selected conservatively.",
        used=stable_unique(
            feature for candidate in candidates for feature in ("splice_outcomes",) + candidate.used
        ),
        exclusions=tuple(
            exclusion for candidate in candidates for exclusion in candidate.exclusions
        ),
    )


def _evaluate_splice(
    evidence: BiologicalFeatures,
    *,
    canonical: bool,
    guideline: GuidelineReference = GENERAL_SVI_GUIDELINE,
    evidence_mode: Optional[str] = None,
) -> _Candidate:
    """Evaluate resolved splice outcomes, choosing the lowest credible strength."""

    is_rna = evidence_mode == "RNA"
    prefix = (
        "SVI2023.RNA_SPLICE"
        if is_rna
        else "SVI2018.CANONICAL_SPLICE"
        if canonical
        else "SVI2023.RNA_SPLICE"
    )
    locator = (
        guideline.source_locator
        if is_rna or not canonical
        else "Figure 1 — GT–AG ±1,2 splice-site branch"
    )
    rescue = evidence.splice_rescue_evidence
    if _is_known_true(rescue):
        return _not_met(
            f"{prefix}.RESCUED",
            locator,
            "A functional splice-rescue isoform is sufficient to prevent application of PVS1.",
            used=("splice_rescue_evidence",),
            exclusions=(_exclusion("splice_rescue_evidence", "EXCLUDES_PVS1"),),
            guideline=guideline,
            evidence_mode=evidence_mode,
        )

    assessment = _resolve_splice_outcomes(
        evidence,
        require_observed=is_rna,
        prefix=prefix,
        locator=locator,
        guideline=guideline,
        evidence_mode=evidence_mode,
    )
    if assessment.candidate is not None:
        return assessment.candidate
    assert assessment.outcomes is not None
    outcomes = assessment.outcomes
    outcome_used = assessment.used

    # A directly resolved normal profile is decisive. A conflicting rescue
    # modifier cannot turn it into a PVS1 result.
    if outcomes == (SpliceOutcome.NO_ABERRANT_SPLICING,):
        return _not_met(
            f"{prefix}.NO_ABERRANT_SPLICING",
            locator,
            "The resolved splice evidence indicates no aberrant splice consequence suitable for PVS1.",
            used=outcome_used,
            guideline=guideline,
            evidence_mode=evidence_mode,
        )

    damaging = tuple(
        outcome for outcome in outcomes if outcome != SpliceOutcome.NO_ABERRANT_SPLICING
    )
    if rescue.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{prefix}.RESCUE_CONFLICT",
            locator,
            "Functional splice-rescue evidence is conflicting and could change the PVS1 result.",
            ("splice_rescue_evidence",),
            used=outcome_used + ("splice_rescue_evidence",),
            guideline=guideline,
            evidence_mode=evidence_mode,
        )

    if len(damaging) != len(outcomes):
        residual = evidence.residual_normal_transcript_sufficient_for_function
        if _is_known_true(residual):
            return _not_met(
                f"{prefix}.RESIDUAL_NORMAL_TRANSCRIPT_SUFFICIENT",
                locator,
                "The residual normal transcript level is assessed as sufficient for function in this gene-disease context.",
                used=outcome_used
                + ("residual_normal_transcript_sufficient_for_function",),
                exclusions=(
                    _exclusion(
                        "residual_normal_transcript_sufficient_for_function",
                        "EXCLUDES_PVS1",
                    ),
                ),
                guideline=guideline,
                evidence_mode=evidence_mode,
            )
        if residual.status == EvidenceStatus.CONFLICTING:
            return _undetermined(
                f"{prefix}.RESIDUAL_NORMAL_TRANSCRIPT_CONFLICT",
                locator,
                "The functional sufficiency of the residual normal transcript is conflicting.",
                ("residual_normal_transcript_sufficient_for_function",),
                used=outcome_used
                + ("residual_normal_transcript_sufficient_for_function",),
                guideline=guideline,
                evidence_mode=evidence_mode,
            )
        if residual.status == EvidenceStatus.UNKNOWN:
            return _undetermined(
                f"{prefix}.MIXED_NORMAL_AND_ABERRANT",
                locator,
                "Both normal and aberrant splice outcomes are credible; a gene-disease-specific assessment of residual normal transcript sufficiency is required.",
                ("residual_normal_transcript_sufficient_for_function",),
                used=outcome_used,
                guideline=guideline,
                evidence_mode=evidence_mode,
            )
        # A known false residual-sufficiency conclusion permits evaluation of
        # the LoF outcomes; it is not a generic numerical rescue threshold.
        outcome_used = outcome_used + (
            "residual_normal_transcript_sufficient_for_function",
        )

    candidates: list[_Candidate] = []
    for outcome in damaging:
        if outcome == SpliceOutcome.OUT_OF_FRAME:
            candidates.append(
                _nmd_or_escape(
                    evidence,
                    f"{prefix}.OUT_OF_FRAME",
                    "The out-of-frame splice outcome",
                    guideline=guideline,
                    evidence_mode=evidence_mode,
                )
            )
        elif outcome == SpliceOutcome.IN_FRAME:
            candidates.append(
                _evaluate_in_frame_splice(
                    evidence,
                    guideline=guideline,
                    evidence_mode=evidence_mode,
                    prefix=f"{prefix}.IN_FRAME",
                )
            )

    unresolved = next(
        (candidate for candidate in candidates if candidate.status == DecisionStatus.UNDETERMINED),
        None,
    )
    if unresolved is not None:
        return replace(unresolved, used=stable_unique(outcome_used + unresolved.used))
    not_met = next(
        (candidate for candidate in candidates if candidate.status == DecisionStatus.NOT_MET),
        None,
    )
    if not_met is not None:
        return replace(not_met, used=stable_unique(outcome_used + not_met.used))
    weakest = min(candidates, key=lambda candidate: _STRENGTH_RANK[candidate.base_strength])
    if len(candidates) == 1:
        return replace(weakest, used=stable_unique(outcome_used + weakest.used))
    return replace(
        weakest,
        path=f"{prefix}.MULTIPLE_OUTCOMES.LOWEST_STRENGTH",
        locator="Figure 1 — apply lowest strength across credible splice outcomes",
        reason="Multiple credible splice outcomes were evaluated; the lowest applicable PVS1 strength was selected conservatively.",
        used=stable_unique(
            feature for candidate in candidates for feature in outcome_used + candidate.used
        ),
        exclusions=tuple(
            exclusion for candidate in candidates for exclusion in candidate.exclusions
        ),
    )


def _frame(evidence: BiologicalFeatures) -> tuple[Optional[ReadingFrame], Optional[_Candidate], tuple[str, ...]]:
    """Use either explicit CNV frame or the redundant boolean, never infer it."""

    explicit = evidence.cnv_reading_frame
    boolean = evidence.reading_frame_disruption
    path = "SVI2018.CNV.READING_FRAME"
    locator = "Figure 1 — reading-frame consequence"
    if explicit.status == EvidenceStatus.CONFLICTING:
        return None, _undetermined(path, locator, "CNV reading-frame evidence is conflicting.", ("cnv_reading_frame",), used=("cnv_reading_frame",)), ()
    if boolean.status == EvidenceStatus.CONFLICTING:
        return None, _undetermined(path, locator, "Reading-frame disruption evidence is conflicting.", ("reading_frame_disruption",), used=("reading_frame_disruption",)), ()
    explicit_value = explicit.value if explicit.status == EvidenceStatus.KNOWN else None
    boolean_value = boolean.value if boolean.status == EvidenceStatus.KNOWN else None
    mapped_boolean = None
    if boolean_value is not None:
        mapped_boolean = ReadingFrame.DISRUPTED if boolean_value else ReadingFrame.PRESERVED
    if explicit_value is not None and mapped_boolean is not None and explicit_value != mapped_boolean:
        return None, _undetermined(
            path,
            locator,
            "The explicit CNV frame and reading-frame disruption inputs disagree.",
            ("cnv_reading_frame",),
            used=("cnv_reading_frame", "reading_frame_disruption"),
        ), ()
    value = explicit_value or mapped_boolean
    if value is None:
        return None, _undetermined(
            path,
            locator,
            "The structural variant's reading-frame consequence is required and cannot be inferred from its type.",
            ("cnv_reading_frame",),
        ), ()
    used = ("cnv_reading_frame",) if explicit_value is not None else ("reading_frame_disruption",)
    return value, None, used


def _require_affected_exons(
    evidence: BiologicalFeatures,
    path: str,
    locator: str,
) -> tuple[Optional[_Candidate], tuple[str, ...]]:
    """Require the event to be demonstrably exon-level before generic CNV PVS1."""

    affected = evidence.affected_exons
    if affected.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            path,
            locator,
            "Affected-exon evidence is conflicting.",
            ("affected_exons",),
            used=("affected_exons",),
        ), ()
    if affected.status == EvidenceStatus.UNKNOWN:
        return _undetermined(
            path,
            locator,
            "The generic structural-variant pathway requires the affected exon(s) to be identified.",
            ("affected_exons",),
        ), ()
    return None, ("affected_exons",)


def _evaluate_exon_deletion(evidence: BiologicalFeatures) -> _Candidate:
    deletion_locator = "Figure 1 — single/multi-exon deletion"
    exon_issue, exon_used = _require_affected_exons(
        evidence,
        "SVI2018.EXON_DELETION.AFFECTED_EXONS",
        deletion_locator,
    )
    if exon_issue is not None:
        return exon_issue
    frame, issue, frame_used = _frame(evidence)
    if issue is not None:
        return replace(
            issue,
            path="SVI2018.EXON_DELETION.READING_FRAME",
            used=stable_unique(exon_used + issue.used),
        )
    assert frame is not None
    if frame == ReadingFrame.DISRUPTED:
        result = _nmd_or_escape(
            evidence,
            "SVI2018.EXON_DELETION.OUT_OF_FRAME",
            "The out-of-frame exon deletion",
        )
        return replace(result, used=stable_unique(exon_used + frame_used + result.used))

    locator = "Figure 1 — single/multi-exon deletion preserving reading frame"
    relevance = _relevance(evidence)
    blocked = _relevance_candidate(relevance, "SVI2018.EXON_DELETION.IN_FRAME.RELEVANCE", locator)
    if blocked:
        return replace(blocked, used=stable_unique(exon_used + frame_used + blocked.used))
    critical = evidence.critical_protein_region
    if critical.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            "SVI2018.EXON_DELETION.IN_FRAME.CRITICAL_CONFLICT",
            locator,
            "Critical-region evidence for the in-frame genomic exon deletion is conflicting.",
            ("critical_protein_region",),
            used=exon_used + frame_used + relevance.used + ("critical_protein_region",),
            exclusions=relevance.exclusions,
        )
    if _is_known_true(critical):
        return _activated(
            PVS1Strength.STRONG,
            "SVI2018.EXON_DELETION.IN_FRAME.CRITICAL_REGION",
            "Figure 1 — in-frame exon deletion affecting critical protein region",
            "The in-frame genomic exon deletion affects a region established as critical to protein function.",
            used=exon_used + frame_used + relevance.used + ("critical_protein_region",),
            exclusions=relevance.exclusions,
        )
    if _is_known_false(critical):
        return _not_met(
            "SVI2018.EXON_DELETION.IN_FRAME.NOT_CRITICAL",
            locator,
            "The generic SVI 2018 in-frame genomic exon-deletion branch supports PVS1 only when the removed region is critical to protein function.",
            used=exon_used + frame_used + relevance.used + ("critical_protein_region",),
            exclusions=relevance.exclusions,
        )
    return _undetermined(
        "SVI2018.EXON_DELETION.IN_FRAME.CRITICAL_UNKNOWN",
        locator,
        "The generic SVI 2018 in-frame genomic exon-deletion branch requires critical-region evidence; it has no generic percentage-loss branch.",
        ("critical_protein_region",),
        used=exon_used + frame_used + relevance.used,
        exclusions=relevance.exclusions,
    )


def _evaluate_whole_gene_deletion(evidence: BiologicalFeatures) -> _Candidate:
    return _activated(
        PVS1Strength.VERY_STRONG,
        "SVI2018.WHOLE_GENE_DELETION",
        "Figure 1 — full gene deletion",
        "A full-gene deletion is a null event in a gene–disease assertion with an eligible LoF mechanism.",
        # Routing may have been established by either a direct
        # ``WHOLE_GENE_DELETION`` consequence or an explicit whole-gene flag.
        # The final selector audit adds whichever known routing fact was
        # actually supplied; do not report an omitted flag as evidence used.
        used=(),
    )


def _evaluate_duplication(evidence: BiologicalFeatures) -> _Candidate:
    prefix = "SVI2018.DUPLICATION"
    # Known exclusions close this generic branch without demanding evidence
    # that would only be needed to prove an activation.
    whole = evidence.whole_gene_duplication
    if _is_known_true(whole):
        return _not_met(
            f"{prefix}.WHOLE_GENE",
            "Figure 1 - intragenic duplication branch",
            "A full-gene tandem duplication is not a generic PVS1 event.",
            used=("whole_gene_duplication",),
            exclusions=(_exclusion("whole_gene_duplication", "EXCLUDES_PVS1"),),
        )

    tandem = evidence.tandem_duplication_status
    if (
        tandem.status == EvidenceStatus.KNOWN
        and tandem.value == TandemDuplicationStatus.PROVEN_NOT_TANDEM
    ):
        return _not_met(
            f"{prefix}.NOT_TANDEM",
            "Figure 1 - intragenic duplication branch",
            "PVS1 does not apply to a duplication proven not to be tandem.",
            used=("tandem_duplication_status",),
            exclusions=(_exclusion("tandem_duplication_status", "EXCLUDES_PVS1"),),
        )

    explicit_frame = evidence.cnv_reading_frame
    boolean_frame = evidence.reading_frame_disruption
    frame_is_known_preserved = (
        _is_known_false(boolean_frame)
        and explicit_frame.status != EvidenceStatus.CONFLICTING
        and not (
            explicit_frame.status == EvidenceStatus.KNOWN
            and explicit_frame.value == ReadingFrame.DISRUPTED
        )
    ) or (
        explicit_frame.status == EvidenceStatus.KNOWN
        and explicit_frame.value == ReadingFrame.PRESERVED
        and boolean_frame.status != EvidenceStatus.CONFLICTING
        and not _is_known_true(boolean_frame)
    )
    if frame_is_known_preserved:
        frame_used = (
            ("cnv_reading_frame",)
            if explicit_frame.status == EvidenceStatus.KNOWN
            else ("reading_frame_disruption",)
        )
        return _not_met(
            f"{prefix}.IN_FRAME",
            "Figure 1 - intragenic duplication branch",
            "The duplication is known not to disrupt the reading frame; generic SVI PVS1 does not apply.",
            used=frame_used,
        )

    nmd = evidence.nmd_prediction
    if _is_known_false(nmd):
        return _not_met(
            f"{prefix}.NMD_NOT_PREDICTED",
            "Figure 1 - intragenic duplication branch",
            "PVS1 does not apply to this duplication because NMD is not predicted.",
            used=("nmd_prediction",),
        )
    locator = "Figure 1 — intragenic duplication branch"
    contained = evidence.duplication_contained_within_gene
    if contained.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{prefix}.CONTAINMENT_CONFLICT", locator,
            "Evidence about whether the duplication is completely contained within the gene is conflicting.",
            ("duplication_contained_within_gene",), used=("duplication_contained_within_gene",),
        )
    if contained.status == EvidenceStatus.UNKNOWN:
        return _undetermined(
            f"{prefix}.CONTAINMENT_UNKNOWN", locator,
            "Generic SVI PVS1 duplication evaluation requires evidence that the duplicated segment is completely contained within the relevant gene.",
            ("duplication_contained_within_gene",),
        )
    if contained.value is False:
        return _not_met(
            f"{prefix}.NOT_INTRAGENIC", locator,
            "The duplication is known not to be completely contained within the relevant gene.",
            used=("duplication_contained_within_gene",),
            exclusions=(_exclusion("duplication_contained_within_gene", "EXCLUDES_PVS1"),),
        )
    whole = evidence.whole_gene_duplication
    if whole.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{prefix}.WHOLE_GENE_CONFLICT", locator,
            "Whether the duplication includes the full gene is conflicting.",
            ("whole_gene_duplication",),
            used=("duplication_contained_within_gene", "whole_gene_duplication"),
        )
    if _is_known_true(whole):
        return _not_met(
            f"{prefix}.WHOLE_GENE", locator,
            "A full-gene tandem duplication is not a generic PVS1 event.",
            used=("duplication_contained_within_gene", "whole_gene_duplication"),
            exclusions=(_exclusion("whole_gene_duplication", "EXCLUDES_PVS1"),),
        )
    if whole.status == EvidenceStatus.UNKNOWN:
        return _undetermined(
            f"{prefix}.WHOLE_GENE_UNKNOWN",
            locator,
            "Generic PVS1 duplication evaluation requires confirmation that the event is not a full-gene duplication.",
            ("whole_gene_duplication",),
            used=("duplication_contained_within_gene",),
        )
    exon_issue, exon_used = _require_affected_exons(
        evidence,
        f"{prefix}.AFFECTED_EXONS",
        locator,
    )
    if exon_issue is not None:
        return replace(
            exon_issue,
            used=stable_unique(("duplication_contained_within_gene",) + exon_issue.used),
        )
    tandem = evidence.tandem_duplication_status
    if tandem.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{prefix}.TANDEM_CONFLICT", locator,
            "Tandem-insertion evidence is conflicting.",
            ("tandem_duplication_status",),
            used=("duplication_contained_within_gene",) + exon_used + ("tandem_duplication_status",),
        )
    if tandem.status == EvidenceStatus.UNKNOWN:
        return _undetermined(
            f"{prefix}.TANDEM_UNKNOWN", locator,
            "The duplication is neither proven nor explicitly assessed as presumed tandem.",
            ("tandem_duplication_status",),
            used=("duplication_contained_within_gene",) + exon_used,
        )
    if tandem.value == TandemDuplicationStatus.PROVEN_NOT_TANDEM:
        return _not_met(
            f"{prefix}.NOT_TANDEM", locator,
            "PVS1 does not apply to a duplication proven not to be tandem.",
            used=("duplication_contained_within_gene",) + exon_used + ("tandem_duplication_status",),
            exclusions=(_exclusion("tandem_duplication_status", "EXCLUDES_PVS1"),),
        )
    frame, issue, frame_used = _frame(evidence)
    if issue is not None:
        return replace(issue, path=f"{prefix}.READING_FRAME", used=stable_unique(("duplication_contained_within_gene",) + exon_used + ("tandem_duplication_status",) + issue.used))
    assert frame is not None
    if frame == ReadingFrame.PRESERVED:
        return _not_met(
            f"{prefix}.IN_FRAME", locator,
            "The duplication is known not to disrupt the reading frame; generic SVI PVS1 does not apply.",
            used=("duplication_contained_within_gene",) + exon_used + ("tandem_duplication_status",) + frame_used,
        )
    nmd = evidence.nmd_prediction
    if nmd.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{prefix}.NMD_CONFLICT", locator,
            "NMD prediction evidence is conflicting.", ("nmd_prediction",),
            used=("duplication_contained_within_gene",) + exon_used + ("tandem_duplication_status",) + frame_used + ("nmd_prediction",),
        )
    if nmd.status == EvidenceStatus.UNKNOWN:
        return _undetermined(
            f"{prefix}.NMD_UNKNOWN", locator,
            "The duplication branch cannot assign PVS1 without a supplied NMD prediction.", ("nmd_prediction",),
            used=("duplication_contained_within_gene",) + exon_used + ("tandem_duplication_status",) + frame_used,
        )
    if nmd.value is False:
        return _not_met(
            f"{prefix}.NMD_NOT_PREDICTED", locator,
            "PVS1 does not apply to this duplication because NMD is not predicted.",
            used=("duplication_contained_within_gene",) + exon_used + ("tandem_duplication_status", "nmd_prediction") + frame_used,
        )
    relevance = _relevance(evidence)
    blocked = _relevance_candidate(relevance, f"{prefix}.RELEVANCE", locator)
    shared = ("duplication_contained_within_gene",) + exon_used + ("tandem_duplication_status", "nmd_prediction") + frame_used
    if blocked:
        return replace(blocked, used=stable_unique(shared + blocked.used))
    strength = (
        PVS1Strength.VERY_STRONG
        if tandem.value == TandemDuplicationStatus.PROVEN_TANDEM
        else PVS1Strength.STRONG
    )
    placement = "proven" if tandem.value == TandemDuplicationStatus.PROVEN_TANDEM else "presumed"
    return _activated(
        strength,
        f"{prefix}.{placement.upper()}_TANDEM.OUT_OF_FRAME.NMD",
        f"Figure 1 — {placement} tandem duplication, disrupted frame, predicted NMD",
        f"The {placement}-tandem intragenic duplication disrupts the reading frame and is predicted to undergo NMD.",
        used=shared + relevance.used,
        exclusions=relevance.exclusions,
    )


def _evaluate_initiation_codon(evidence: BiologicalFeatures) -> _Candidate:
    prefix = "SVI2018.INITIATION_CODON"
    locator = "Figure 1 — initiation codon branch"
    alternative = evidence.alternative_start_codon
    # A functional alternative start is a decisive branch exclusion. It is
    # independent of the later relevance prerequisite for a positive result,
    # so a missing relevance assertion must not hide this known negative.
    if alternative.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{prefix}.ALTERNATIVE_START_CONFLICT", locator,
            "Evidence about functional alternative start codons is conflicting.",
            ("alternative_start_codon",), used=("alternative_start_codon",),
        )
    if _is_known_true(alternative):
        return _not_met(
            f"{prefix}.ALTERNATIVE_START_PRESENT", locator,
            "A different functional transcript uses an alternative start codon, so generic PVS1 does not apply.",
            used=("alternative_start_codon",),
            exclusions=(_exclusion("alternative_start_codon", "EXCLUDES_PVS1"),),
        )
    relevance = _relevance(evidence)
    blocked = _relevance_candidate(relevance, f"{prefix}.RELEVANCE", locator)
    if blocked:
        return blocked
    if alternative.status == EvidenceStatus.UNKNOWN:
        return _undetermined(
            f"{prefix}.ALTERNATIVE_START_UNKNOWN", locator,
            "Initiation-codon PVS1 evaluation requires assessment of functional alternative start codons.",
            ("alternative_start_codon",), used=relevance.used, exclusions=relevance.exclusions,
        )
    upstream = evidence.upstream_pathogenic_variant_evidence
    if upstream.status == EvidenceStatus.CONFLICTING:
        return _undetermined(
            f"{prefix}.UPSTREAM_PATHOGENIC_CONFLICT", locator,
            "Evidence for pathogenic variants upstream of the closest downstream in-frame start is conflicting.",
            ("upstream_pathogenic_variant_evidence",),
            used=relevance.used + ("alternative_start_codon", "upstream_pathogenic_variant_evidence"), exclusions=relevance.exclusions,
        )
    if upstream.status == EvidenceStatus.UNKNOWN:
        return _undetermined(
            f"{prefix}.UPSTREAM_PATHOGENIC_UNKNOWN", locator,
            "The initiation-codon branch requires assessment of pathogenic variants upstream of the closest downstream in-frame start.",
            ("upstream_pathogenic_variant_evidence",),
            used=relevance.used + ("alternative_start_codon",), exclusions=relevance.exclusions,
        )
    strength = PVS1Strength.MODERATE if upstream.value else PVS1Strength.SUPPORTING
    branch = "UPSTREAM_PATHOGENIC" if upstream.value else "NO_UPSTREAM_PATHOGENIC"
    return _activated(
        strength,
        f"{prefix}.{branch}",
        f"Figure 1 — no functional alternative start; {branch.replace('_', ' ').lower()}",
        "No functional alternative start codon is known; the documented upstream pathogenic-variant evidence supports the applicable generic initiation-codon strength.",
        used=relevance.used + ("alternative_start_codon", "upstream_pathogenic_variant_evidence"),
        exclusions=relevance.exclusions,
    )


def _identify_pathway_legacy(evidence: BiologicalFeatures) -> _PathwayAssessment:
    """Historical routing retained temporarily for source compatibility tests.

    `PVS1Engine` uses `_identify_pathway_checked`, which adds cross-field
    consistency handling while preserving the same generic branch outcomes.
    """

    event = evidence.cnv_event_type
    whole_deletion = evidence.whole_gene_deletion
    if event.status == EvidenceStatus.CONFLICTING:
        return _PathwayAssessment(candidate=_undetermined(
            "SVI2018.PATHWAY.CNV_EVENT_CONFLICT", "Figure 1 — variant class",
            "The structural-event type is conflicting.", ("cnv_event_type",), used=("cnv_event_type",),
        ))
    if whole_deletion.status == EvidenceStatus.CONFLICTING:
        return _PathwayAssessment(candidate=_undetermined(
            "SVI2018.PATHWAY.WHOLE_GENE_DELETION_CONFLICT", "Figure 1 — deletion branch",
            "Evidence about whether the deletion covers the full gene is conflicting.", ("whole_gene_deletion",), used=("whole_gene_deletion",),
        ))
    if _is_known_true(whole_deletion):
        if event.status == EvidenceStatus.KNOWN and event.value == StructuralEvent.DUPLICATION:
            return _PathwayAssessment(candidate=_undetermined(
                "SVI2018.PATHWAY.STRUCTURAL_CONFLICT", "Figure 1 — variant class",
                "A full-gene deletion conflicts with a duplication event type.", ("cnv_event_type",),
                used=("whole_gene_deletion", "cnv_event_type"),
            ))
        return _PathwayAssessment(pathway=Pathway.WHOLE_GENE_DELETION)
    if event.status == EvidenceStatus.KNOWN:
        if event.value == StructuralEvent.DELETION:
            if whole_deletion.status == EvidenceStatus.UNKNOWN:
                if evidence.affected_exons.status == EvidenceStatus.KNOWN:
                    return _PathwayAssessment(pathway=Pathway.EXON_DELETION)
                if evidence.affected_exons.status == EvidenceStatus.CONFLICTING:
                    return _PathwayAssessment(candidate=_undetermined(
                        "SVI2018.PATHWAY.DELETION_EXTENT_CONFLICT", "Figure 1 — deletion branch",
                        "Affected-exon evidence is conflicting and whole-gene deletion status is unknown.", ("affected_exons",),
                        used=("cnv_event_type", "affected_exons"),
                    ))
                return _PathwayAssessment(candidate=_undetermined(
                    "SVI2018.PATHWAY.DELETION_EXTENT_UNKNOWN", "Figure 1 — deletion branch",
                    "The deletion pathway requires knowing whether the event is whole-gene or exon-level.", ("whole_gene_deletion",),
                    used=("cnv_event_type",),
                ))
            return _PathwayAssessment(pathway=Pathway.EXON_DELETION)
        return _PathwayAssessment(pathway=Pathway.DUPLICATION)

    consequence = evidence.variant_consequence
    if consequence.status == EvidenceStatus.CONFLICTING:
        return _PathwayAssessment(candidate=_undetermined(
            "SVI2018.PATHWAY.CONSEQUENCE_CONFLICT", "Figure 1 — variant class",
            "Variant consequence evidence is conflicting.", ("variant_consequence",), used=("variant_consequence",),
        ))
    if consequence.status == EvidenceStatus.UNKNOWN:
        if _is_known_true(evidence.start_codon_affected):
            return _PathwayAssessment(pathway=Pathway.INITIATION_CODON)
        return _PathwayAssessment(candidate=_undetermined(
            "SVI2018.PATHWAY.CONSEQUENCE_UNKNOWN", "Figure 1 — variant class",
            "A supported PVS1 variant consequence has not been supplied.", ("variant_consequence",),
        ))
    if consequence.value in {VariantConsequence.STOP_GAIN, VariantConsequence.FRAMESHIFT}:
        return _PathwayAssessment(pathway=Pathway.NULL_VARIANT)
    if consequence.value == VariantConsequence.CANONICAL_SPLICE:
        return _PathwayAssessment(pathway=Pathway.CANONICAL_SPLICE)
    if consequence.value == VariantConsequence.SPLICE_SITE:
        position = evidence.splice_site_position
        if position.status == EvidenceStatus.CONFLICTING:
            return _PathwayAssessment(candidate=_undetermined(
                "SVI2018.PATHWAY.SPLICE_POSITION_CONFLICT", "Figure 1 — GT–AG ±1,2 splice-site branch",
                "Splice-site position evidence is conflicting.", ("splice_site_position",), used=("variant_consequence", "splice_site_position"),
            ))
        if position.status == EvidenceStatus.UNKNOWN:
            return _PathwayAssessment(candidate=_undetermined(
                "SVI2018.PATHWAY.SPLICE_POSITION_UNKNOWN", "Figure 1 — GT–AG ±1,2 splice-site branch",
                "A generic splice-site consequence requires canonical ±1/2 position or RNA evidence for a non-canonical pathway.", ("splice_site_position",), used=("variant_consequence",),
            ))
        if abs(position.value) in {1, 2}:
            return _PathwayAssessment(pathway=Pathway.CANONICAL_SPLICE)
        return _PathwayAssessment(pathway=Pathway.NONCANONICAL_SPLICE)
    if consequence.value == VariantConsequence.START_LOST:
        return _PathwayAssessment(pathway=Pathway.INITIATION_CODON)
    if consequence.value == VariantConsequence.EXON_DELETION:
        return _PathwayAssessment(pathway=Pathway.EXON_DELETION)
    if consequence.value == VariantConsequence.WHOLE_GENE_DELETION:
        return _PathwayAssessment(pathway=Pathway.WHOLE_GENE_DELETION)
    if consequence.value == VariantConsequence.DUPLICATION:
        return _PathwayAssessment(pathway=Pathway.DUPLICATION)
    return _PathwayAssessment(candidate=_not_applicable(
        "SVI2018.PATHWAY.UNSUPPORTED_CONSEQUENCE",
        "The supplied variant consequence is outside the generic PVS1-supported variant classes.",
        used=("variant_consequence",),
    ))


def _pathway_issue(
    path: str,
    reason: str,
    missing: Iterable[str],
    used: Iterable[str],
) -> _PathwayAssessment:
    return _PathwayAssessment(
        candidate=_undetermined(
            path,
            "Figure 1 - variant class",
            reason,
            missing,
            used=used,
        )
    )


def _identify_pathway_checked(evidence: BiologicalFeatures) -> _PathwayAssessment:
    """Route generic evidence after reconciling material class assertions.

    This intentionally does not derive a copy-number class, splice position,
    or start-loss effect.  It only selects a compatible supplied branch or
    returns the smallest class fact that must be reconciled.
    """

    consequence = evidence.variant_consequence
    event = evidence.cnv_event_type
    whole_deletion = evidence.whole_gene_deletion
    whole_duplication = evidence.whole_gene_duplication
    start = evidence.start_codon_affected

    if consequence.status == EvidenceStatus.CONFLICTING:
        return _pathway_issue(
            "SVI2018.PATHWAY.CONSEQUENCE_CONFLICT",
            "Variant consequence evidence is conflicting.",
            ("variant_consequence",),
            ("variant_consequence",),
        )
    if event.status == EvidenceStatus.CONFLICTING:
        return _pathway_issue(
            "SVI2018.PATHWAY.CNV_EVENT_CONFLICT",
            "The structural-event type is conflicting.",
            ("cnv_event_type",),
            ("cnv_event_type",),
        )

    consequence_value = (
        consequence.value if consequence.status == EvidenceStatus.KNOWN else None
    )

    def conflict(reason: str, fields: tuple[str, ...]) -> _PathwayAssessment:
        return _pathway_issue(
            "SVI2018.PATHWAY.STRUCTURAL_CONFLICT",
            reason,
            fields,
            fields,
        )

    # This evaluator receives evidence for one variant/event. An independently
    # supplied small-variant/CNV anchor cannot be reconciled with a conflicting
    # assertion that the same event is a full-gene deletion or duplication.
    # Treat this as a material preflight conflict before a VCEP rule can take
    # precedence. Position/exon anchors are included because a scoped VCEP
    # rule can validly use them when a generic consequence is absent.
    if (
        consequence_value in {
            VariantConsequence.STOP_GAIN,
            VariantConsequence.FRAMESHIFT,
            VariantConsequence.CANONICAL_SPLICE,
            VariantConsequence.SPLICE_SITE,
            VariantConsequence.START_LOST,
            VariantConsequence.OTHER,
        }
        # A supplied start-codon effect is also a direct small-variant route
        # when no consequence label is available.
        or _is_known_true(start)
        or evidence.splice_site_position.status == EvidenceStatus.KNOWN
        or evidence.affected_exons.status == EvidenceStatus.KNOWN
    ):
        conflicting_full_gene = tuple(
            name
            for name, value in (
                ("whole_gene_deletion", whole_deletion),
                ("whole_gene_duplication", whole_duplication),
            )
            if value.status == EvidenceStatus.CONFLICTING
        )
        if conflicting_full_gene:
            return conflict(
                "The small-variant consequence cannot be reconciled with conflicting full-gene structural evidence.",
                ("variant_consequence",) + conflicting_full_gene,
            )

    # Explicit CNV event and consequence must describe compatible classes.
    if event.status == EvidenceStatus.KNOWN:
        if event.value == StructuralEvent.DELETION:
            if consequence_value not in {
                None,
                VariantConsequence.EXON_DELETION,
                VariantConsequence.WHOLE_GENE_DELETION,
            }:
                return conflict(
                    "The deletion event type conflicts with the supplied variant consequence.",
                    ("cnv_event_type", "variant_consequence"),
                )
            if _is_known_true(whole_duplication) or whole_duplication.status == EvidenceStatus.CONFLICTING:
                return conflict(
                    "Deletion and whole-gene-duplication assertions cannot both describe the same event.",
                    ("cnv_event_type", "whole_gene_duplication"),
                )
            if whole_deletion.status == EvidenceStatus.CONFLICTING:
                return _pathway_issue(
                    "SVI2018.PATHWAY.WHOLE_GENE_DELETION_CONFLICT",
                    "Evidence about whether the deletion covers the full gene is conflicting.",
                    ("whole_gene_deletion",),
                    ("cnv_event_type", "whole_gene_deletion"),
                )
            if _is_known_true(whole_deletion):
                if consequence_value == VariantConsequence.EXON_DELETION:
                    return conflict(
                        "An exon-deletion consequence conflicts with an explicit full-gene deletion assertion.",
                        ("variant_consequence", "whole_gene_deletion"),
                    )
                return _PathwayAssessment(pathway=Pathway.WHOLE_GENE_DELETION)
            if consequence_value == VariantConsequence.WHOLE_GENE_DELETION:
                if _is_known_false(whole_deletion):
                    return conflict(
                        "A whole-gene deletion consequence conflicts with an explicit non-whole-gene assertion.",
                        ("variant_consequence", "whole_gene_deletion"),
                    )
                return _PathwayAssessment(pathway=Pathway.WHOLE_GENE_DELETION)
            if consequence_value == VariantConsequence.EXON_DELETION or _is_known_false(whole_deletion):
                return _PathwayAssessment(pathway=Pathway.EXON_DELETION)
            if evidence.affected_exons.status == EvidenceStatus.KNOWN:
                return _PathwayAssessment(pathway=Pathway.EXON_DELETION)
            if evidence.affected_exons.status == EvidenceStatus.CONFLICTING:
                return _pathway_issue(
                    "SVI2018.PATHWAY.DELETION_EXTENT_CONFLICT",
                    "Affected-exon evidence is conflicting while deletion extent is unresolved.",
                    ("affected_exons",),
                    ("cnv_event_type", "affected_exons"),
                )
            return _pathway_issue(
                "SVI2018.PATHWAY.DELETION_EXTENT_UNKNOWN",
                "The deletion pathway requires knowing whether the event is whole-gene or exon-level.",
                ("whole_gene_deletion",),
                ("cnv_event_type",),
            )

        if consequence_value not in {None, VariantConsequence.DUPLICATION}:
            return conflict(
                "The duplication event type conflicts with the supplied variant consequence.",
                ("cnv_event_type", "variant_consequence"),
            )
        if _is_known_true(whole_deletion) or whole_deletion.status == EvidenceStatus.CONFLICTING:
            return conflict(
                "Duplication and full-gene-deletion assertions cannot both describe the same event.",
                ("cnv_event_type", "whole_gene_deletion"),
            )
        if whole_duplication.status == EvidenceStatus.CONFLICTING:
            return _pathway_issue(
                "SVI2018.PATHWAY.WHOLE_GENE_DUPLICATION_CONFLICT",
                "Evidence about whether the duplication covers the full gene is conflicting.",
                ("whole_gene_duplication",),
                ("cnv_event_type", "whole_gene_duplication"),
            )
        return _PathwayAssessment(pathway=Pathway.DUPLICATION)

    # Direct structural consequences can route the generic branch without an
    # explicit `cnv_event_type`.
    if consequence_value in {
        VariantConsequence.EXON_DELETION,
        VariantConsequence.WHOLE_GENE_DELETION,
    }:
        if _is_known_true(whole_duplication) or whole_duplication.status == EvidenceStatus.CONFLICTING:
            return conflict(
                "Deletion and whole-gene-duplication assertions cannot both describe the same event.",
                ("variant_consequence", "whole_gene_duplication"),
            )
        if whole_deletion.status == EvidenceStatus.CONFLICTING:
            return _pathway_issue(
                "SVI2018.PATHWAY.WHOLE_GENE_DELETION_CONFLICT",
                "Evidence about whether the deletion covers the full gene is conflicting.",
                ("whole_gene_deletion",),
                ("variant_consequence", "whole_gene_deletion"),
            )
        if _is_known_true(whole_deletion):
            if consequence_value == VariantConsequence.EXON_DELETION:
                return conflict(
                    "An exon-deletion consequence conflicts with an explicit full-gene deletion assertion.",
                    ("variant_consequence", "whole_gene_deletion"),
                )
            return _PathwayAssessment(pathway=Pathway.WHOLE_GENE_DELETION)
        if consequence_value == VariantConsequence.WHOLE_GENE_DELETION:
            if _is_known_false(whole_deletion):
                return conflict(
                    "A whole-gene deletion consequence conflicts with an explicit non-whole-gene assertion.",
                    ("variant_consequence", "whole_gene_deletion"),
                )
            return _PathwayAssessment(pathway=Pathway.WHOLE_GENE_DELETION)
        return _PathwayAssessment(pathway=Pathway.EXON_DELETION)

    if consequence_value == VariantConsequence.DUPLICATION:
        if _is_known_true(whole_deletion) or whole_deletion.status == EvidenceStatus.CONFLICTING:
            return conflict(
                "Duplication and full-gene-deletion assertions cannot both describe the same event.",
                ("variant_consequence", "whole_gene_deletion"),
            )
        if whole_duplication.status == EvidenceStatus.CONFLICTING:
            return _pathway_issue(
                "SVI2018.PATHWAY.WHOLE_GENE_DUPLICATION_CONFLICT",
                "Evidence about whether the duplication covers the full gene is conflicting.",
                ("whole_gene_duplication",),
                ("variant_consequence", "whole_gene_duplication"),
            )
        return _PathwayAssessment(pathway=Pathway.DUPLICATION)

    # Full-gene assertions are themselves a structural route only if no
    # incompatible small-variant consequence was supplied.
    if _is_known_true(whole_deletion) or _is_known_true(whole_duplication):
        if consequence_value is not None:
            return conflict(
                "The full-gene structural assertion conflicts with the supplied variant consequence.",
                ("variant_consequence", "whole_gene_deletion", "whole_gene_duplication"),
            )
        if (
            _is_known_true(whole_deletion)
            and whole_duplication.status == EvidenceStatus.CONFLICTING
        ) or (
            _is_known_true(whole_duplication)
            and whole_deletion.status == EvidenceStatus.CONFLICTING
        ):
            return conflict(
                "A full-gene structural assertion cannot be reconciled with conflicting evidence for the opposite structural event.",
                ("whole_gene_deletion", "whole_gene_duplication"),
            )
        if _is_known_true(whole_deletion) and _is_known_true(whole_duplication):
            return conflict(
                "Full-gene deletion and duplication assertions cannot both describe the same event.",
                ("whole_gene_deletion", "whole_gene_duplication"),
            )
        return _PathwayAssessment(
            pathway=(
                Pathway.WHOLE_GENE_DELETION
                if _is_known_true(whole_deletion)
                else Pathway.DUPLICATION
            )
        )

    if consequence.status == EvidenceStatus.UNKNOWN:
        if start.status == EvidenceStatus.CONFLICTING:
            return _pathway_issue(
                "SVI2018.PATHWAY.START_CODON_CONFLICT",
                "Start-codon involvement is conflicting while no variant consequence is available.",
                ("start_codon_affected",),
                ("start_codon_affected",),
            )
        if _is_known_true(start):
            return _PathwayAssessment(pathway=Pathway.INITIATION_CODON)
        return _pathway_issue(
            "SVI2018.PATHWAY.CONSEQUENCE_UNKNOWN",
            "A supported PVS1 variant consequence has not been supplied.",
            ("variant_consequence",),
            (),
        )

    # Met1 stop-gain/frameshift annotations are routed to the initiation
    # branch when the supplied start-codon fact says that codon is affected.
    if start.status == EvidenceStatus.CONFLICTING:
        return _pathway_issue(
            "SVI2018.PATHWAY.START_CODON_CONFLICT",
            "Start-codon involvement is conflicting and could change the applicable PVS1 pathway.",
            ("start_codon_affected",),
            ("variant_consequence", "start_codon_affected"),
        )
    if _is_known_true(start):
        if consequence_value in {
            VariantConsequence.STOP_GAIN,
            VariantConsequence.FRAMESHIFT,
            VariantConsequence.START_LOST,
            VariantConsequence.OTHER,
        }:
            return _PathwayAssessment(pathway=Pathway.INITIATION_CODON)
        return conflict(
            "The supplied start-codon assertion conflicts with the non-initiation variant class.",
            ("variant_consequence", "start_codon_affected"),
        )
    if consequence_value == VariantConsequence.START_LOST:
        if _is_known_false(start):
            return conflict(
                "A start-loss consequence conflicts with an explicit negative start-codon assertion.",
                ("variant_consequence", "start_codon_affected"),
            )
        return _PathwayAssessment(pathway=Pathway.INITIATION_CODON)
    if consequence_value in {VariantConsequence.STOP_GAIN, VariantConsequence.FRAMESHIFT}:
        return _PathwayAssessment(pathway=Pathway.NULL_VARIANT)
    if consequence_value == VariantConsequence.CANONICAL_SPLICE:
        position = evidence.splice_site_position
        if position.status == EvidenceStatus.CONFLICTING:
            return _pathway_issue(
                "SVI2018.PATHWAY.SPLICE_POSITION_CONFLICT",
                "Canonical-splice and splice-position evidence is conflicting.",
                ("splice_site_position",),
                ("variant_consequence", "splice_site_position"),
            )
        if position.status == EvidenceStatus.KNOWN and abs(position.value) not in {1, 2}:
            return conflict(
                "A canonical splice consequence conflicts with a non-canonical splice-site position.",
                ("variant_consequence", "splice_site_position"),
            )
        return _PathwayAssessment(pathway=Pathway.CANONICAL_SPLICE)
    if consequence_value == VariantConsequence.SPLICE_SITE:
        position = evidence.splice_site_position
        if position.status == EvidenceStatus.CONFLICTING:
            return _pathway_issue(
                "SVI2018.PATHWAY.SPLICE_POSITION_CONFLICT",
                "Splice-site position evidence is conflicting.",
                ("splice_site_position",),
                ("variant_consequence", "splice_site_position"),
            )
        if position.status == EvidenceStatus.UNKNOWN:
            if _has_direct_rna_splice_assertion(evidence):
                # Canonical position is required to select the prediction
                # branch, but an observed allele-relevant RNA assay supplies
                # the splice outcome directly.
                return _PathwayAssessment(pathway=Pathway.NONCANONICAL_SPLICE)
            return _pathway_issue(
                "SVI2018.PATHWAY.SPLICE_POSITION_UNKNOWN",
                "A generic splice-site consequence requires canonical +/-1/2 position or RNA evidence for a non-canonical pathway.",
                ("splice_site_position",),
                ("variant_consequence",),
            )
        return _PathwayAssessment(
            pathway=(
                Pathway.CANONICAL_SPLICE
                if abs(position.value) in {1, 2}
                else Pathway.NONCANONICAL_SPLICE
            )
        )
    if (
        consequence_value == VariantConsequence.OTHER
        and _has_direct_rna_splice_assertion(evidence)
    ):
        # A synonymous/intronic/other sequence annotation can still have a
        # directly observed splice consequence. Do not let the annotation
        # label hide a qualified PVS1(RNA) route.
        return _PathwayAssessment(pathway=Pathway.NONCANONICAL_SPLICE)
    return _PathwayAssessment(candidate=_not_applicable(
        "SVI2018.PATHWAY.UNSUPPORTED_CONSEQUENCE",
        "The supplied variant consequence is outside the generic PVS1-supported variant classes.",
        used=("variant_consequence",),
    ))


def _pathway_hint(evidence: BiologicalFeatures) -> Optional[str]:
    """Infer only a broad route for VCEP filtering, never a generic result.

    Unlike `_identify_pathway`, this intentionally does not demand all facts
    required by the generic flowchart.  It keeps a path-specific VCEP rule
    from being applied to an obviously different variant class, while letting
    an explicit ``ANY`` rule handle an unclassifiable generic splice call.
    """

    consequence = evidence.variant_consequence
    event = evidence.cnv_event_type
    if consequence.status == EvidenceStatus.CONFLICTING or event.status == EvidenceStatus.CONFLICTING:
        return None

    if event.status == EvidenceStatus.KNOWN:
        if consequence.status == EvidenceStatus.KNOWN:
            if event.value == StructuralEvent.DELETION and consequence.value not in {
                VariantConsequence.EXON_DELETION,
                VariantConsequence.WHOLE_GENE_DELETION,
            }:
                return None
            if event.value == StructuralEvent.DUPLICATION and consequence.value != VariantConsequence.DUPLICATION:
                return None
        if event.value == StructuralEvent.DELETION:
            if _is_known_true(evidence.whole_gene_deletion) or (
                consequence.status == EvidenceStatus.KNOWN
                and consequence.value == VariantConsequence.WHOLE_GENE_DELETION
            ):
                return Pathway.WHOLE_GENE_DELETION
            return Pathway.EXON_DELETION
        return Pathway.DUPLICATION

    if consequence.status == EvidenceStatus.KNOWN:
        value = consequence.value
        if _is_known_true(evidence.start_codon_affected) and value in {
            VariantConsequence.STOP_GAIN,
            VariantConsequence.FRAMESHIFT,
            VariantConsequence.START_LOST,
            VariantConsequence.OTHER,
        }:
            return Pathway.INITIATION_CODON
        mapping = {
            VariantConsequence.STOP_GAIN: Pathway.NULL_VARIANT,
            VariantConsequence.FRAMESHIFT: Pathway.NULL_VARIANT,
            VariantConsequence.CANONICAL_SPLICE: Pathway.CANONICAL_SPLICE,
            VariantConsequence.START_LOST: Pathway.INITIATION_CODON,
            VariantConsequence.EXON_DELETION: Pathway.EXON_DELETION,
            VariantConsequence.WHOLE_GENE_DELETION: Pathway.WHOLE_GENE_DELETION,
            VariantConsequence.DUPLICATION: Pathway.DUPLICATION,
        }
        if value == VariantConsequence.SPLICE_SITE:
            position = evidence.splice_site_position
            if position.status == EvidenceStatus.KNOWN:
                return (
                    Pathway.CANONICAL_SPLICE
                    if abs(position.value) in {1, 2}
                    else Pathway.NONCANONICAL_SPLICE
                )
            if _has_direct_rna_splice_assertion(evidence):
                return Pathway.NONCANONICAL_SPLICE
            return None
        if value == VariantConsequence.OTHER and _has_direct_rna_splice_assertion(evidence):
            return Pathway.NONCANONICAL_SPLICE
        return mapping.get(value)

    if _is_known_true(evidence.start_codon_affected):
        return Pathway.INITIATION_CODON
    return None


def _vcep_candidate(evidence: BiologicalFeatures) -> Optional[_Candidate]:
    """Resolve a caller-pinned VCEP declaration before generic routing.

    A VCEP manifest can contain a disease-specific condition that is more
    informative than a generic annotation.  It therefore has to be consulted
    before a missing generic splice position (or a generic unsupported
    consequence) produces a terminal result.  The manifest is still required
    to be versioned, scoped, and explicit; unresolved or overlapping rules are
    returned as safe ``UNDETERMINED`` results.
    """

    hint = _pathway_hint(evidence)
    specific = select_specific_rule(
        evidence,
        hint,
        only_any=hint is None,
        permit_anchored_path_specific=hint is None,
    )
    if specific.conflict:
        candidate = _undetermined(
            "VCEP.SELECTION.CONFLICT",
            "Applicable VCEP specification selection",
            specific.conflict,
            (),
            used=("gene_specific_pvs1_rules",),
        )
        return replace(
            candidate,
            guideline_conditions={
                "framework": "caller-supplied VCEP specification selection",
                "selection": "CONFLICT",
                "detail": specific.conflict,
            },
        )
    if specific.missing:
        candidate = _undetermined(
            "VCEP.SELECTION.UNRESOLVED",
            "Applicable VCEP specification selection",
            "A higher-precedence VCEP rule may apply but its scope or conditions are unresolved.",
            specific.missing,
            used=("gene_specific_pvs1_rules",) + tuple(specific.missing),
        )
        return replace(
            candidate,
            guideline_conditions={
                "framework": "caller-supplied VCEP specification selection",
                "selection": "UNRESOLVED",
                "missing_scope_or_conditions": tuple(specific.missing),
            },
        )
    if specific.matched is None:
        return None

    spec, rule = specific.matched
    scope_conditions = {
        field: value
        for field, value in (
            ("gene_id", spec.scope.gene_id),
            ("disease_id", spec.scope.disease_id),
            ("transcript_id", spec.scope.transcript_id),
            ("transcript_version", spec.scope.transcript_version),
            ("mode_of_inheritance", spec.scope.mode_of_inheritance),
        )
        if value is not None
    }
    return _Candidate(
        status=rule.status,
        base_strength=rule.strength,
        path=rule.decision_path,
        locator=rule.source_locator,
        reason=(
            "A scoped, versioned VCEP PVS1 rule explicitly matches the "
            "supplied evidence and takes precedence over the generic SVI branch."
        ),
        used=stable_unique(
            ("gene_specific_pvs1_rules",)
            + tuple(scope_conditions)
            + tuple(rule.when)
        ),
        guideline=guideline_for_specific_rule(spec, rule),
        guideline_conditions={
            "scope": scope_conditions,
            "pathway": rule.pathway,
            "when": dict(rule.when),
        },
        evidence_mode=rule.evidence_mode,
    )


def _pathway_selector_fields(
    candidate: _Candidate, evidence: BiologicalFeatures
) -> tuple[str, ...]:
    """Return known input facts that selected the generic decision branch.

    Branch helpers intentionally focus on their local minimum evidence.  The
    final audit must additionally show the consequence/structural fact that
    selected that branch, without serializing unrelated supplied inputs.
    """

    if candidate.path.startswith("VCEP.SELECTION"):
        return ()
    path = str(candidate.guideline_conditions.get("pathway", candidate.path))
    candidates: tuple[str, ...]
    if "WHOLE_GENE_DELETION" in path:
        candidates = ("whole_gene_deletion", "cnv_event_type", "variant_consequence")
    elif "EXON_DELETION" in path:
        candidates = ("cnv_event_type", "variant_consequence", "whole_gene_deletion")
    elif "DUPLICATION" in path:
        candidates = ("cnv_event_type", "variant_consequence", "whole_gene_duplication")
    elif "INITIATION_CODON" in path:
        candidates = ("variant_consequence", "start_codon_affected")
    elif "SPLICE" in path:
        candidates = ("variant_consequence", "splice_site_position")
    elif ".NF" in path or "NONSENSE" in path:
        candidates = ("variant_consequence", "start_codon_affected")
    else:
        candidates = ("variant_consequence",)
    return tuple(
        field
        for field in candidates
        if evidence.evidence(field).status == EvidenceStatus.KNOWN
        and not (
            field == "start_codon_affected"
            and evidence.start_codon_affected.value is False
        )
    )


def _audit_feature_name(name: str) -> str:
    """Map public aliases back to the validated schema for audit snapshots."""

    return {"splice_effect": "splice_outcomes"}.get(name, name)


def _generic_guideline_conditions(
    candidate: _Candidate,
    evidence: BiologicalFeatures,
    used: Iterable[str],
) -> Mapping[str, Any]:
    """Create an auditable condition record for a generic SVI/RNA branch.

    The published flowcharts do not expose a machine-readable condition DSL.
    Recording the deterministic engine branch, source locator, and the exact
    normalized values used to satisfy it is therefore preferable to inventing
    guideline rule identifiers or rewriting the source as pseudo-rules.
    """

    evaluated: dict[str, Mapping[str, Any]] = {}
    for field in used:
        canonical = _audit_feature_name(field)
        if canonical in BiologicalFeatures.__dataclass_fields__:
            evaluated[field] = evidence.evidence(canonical).to_dict()
    if (
        candidate.evidence_mode == "RNA"
        and candidate.guideline == RNA_SPLICE_GUIDELINE
    ):
        result: dict[str, Any] = {
            "framework": (
                "ClinGen SVI PVS1(RNA) direct-splice evidence with "
                "ClinGen SVI 2018 strength mapping"
            ),
            "branch": candidate.path,
            "rna_evidence_source": {
                "authority": RNA_SPLICE_GUIDELINE.authority,
                "name": RNA_SPLICE_GUIDELINE.document,
                "version": RNA_SPLICE_GUIDELINE.version,
                "source_url": RNA_SPLICE_GUIDELINE.source_url,
                "source_locator": RNA_SPLICE_GUIDELINE.source_locator,
                "rule_id": RNA_SPLICE_GUIDELINE.rule_id,
            },
            "evaluated_conditions": evaluated,
        }
        # RNA guidance establishes that the assay/profile can be used for
        # PVS1. Where the resolved product enters the published Figure 1
        # NMD/in-frame strength ladder, retain that second source explicitly
        # instead of attaching its locator to the 2023 RNA document.
        if candidate.locator.startswith("Figure 1"):
            result["strength_resolution_source"] = {
                "authority": GENERAL_SVI_GUIDELINE.authority,
                "name": GENERAL_SVI_GUIDELINE.document,
                "version": GENERAL_SVI_GUIDELINE.version,
                "source_url": GENERAL_SVI_GUIDELINE.source_url,
                "source_locator": candidate.locator,
                "rule_id": GENERAL_SVI_GUIDELINE.rule_id,
            }
        return result
    return {
        "framework": "ClinGen SVI generic PVS1 decision tree",
        "branch": candidate.path,
        "source_locator": candidate.locator,
        "evaluated_conditions": evaluated,
    }


def _evaluate_generic_pathway(
    evidence: BiologicalFeatures, pathway: str
) -> _Candidate:
    """Evaluate a routed generic branch without applying LoF eligibility/caps.

    This separation lets a known pathway exclusion remain decisive even if the
    positive gene-disease LoF eligibility evidence is absent. Only an
    activation is subsequently gated and capped by Table 1.
    """

    if pathway == Pathway.NULL_VARIANT:
        return _evaluate_null_variant(evidence)
    if pathway == Pathway.CANONICAL_SPLICE:
        rna = _assess_rna_profile(evidence)
        if rna.candidate is not None and not rna.canonical_generic_fallback:
            return rna.candidate
        if rna.state == "READY":
            candidate = _evaluate_splice(
                evidence,
                canonical=True,
                guideline=RNA_SPLICE_GUIDELINE,
                evidence_mode="RNA",
            )
            return replace(
                candidate,
                used=stable_unique(rna.used + candidate.used),
            )
        # Canonical +/-1/2 prediction remains a valid generic route when no
        # usable direct RNA assay is supplied. An incomplete positive RNA
        # assertion remains terminal for the non-canonical RNA-only route,
        # but does not add a generic canonical gate.
        return _evaluate_splice(evidence, canonical=True)
    if pathway == Pathway.NONCANONICAL_SPLICE:
        rna = _assess_rna_profile(evidence)
        if rna.candidate is not None:
            return rna.candidate
        if rna.state == "READY":
            candidate = _evaluate_splice(
                evidence,
                canonical=False,
                guideline=RNA_SPLICE_GUIDELINE,
                evidence_mode="RNA",
            )
            return replace(candidate, used=stable_unique(rna.used + candidate.used))
        return _not_applicable(
            "SVI2023.RNA_SPLICE.NOT_ESTABLISHED",
            "A non-canonical splice variant requires a usable observed RNA assay with an aberrant LoF profile (or a scoped VCEP rule) before PVS1 can be evaluated. "
            + rna.reason,
            used=tuple(
                name
                for name in ("variant_consequence", "splice_site_position")
                if evidence.evidence(name).status == EvidenceStatus.KNOWN
            )
            + rna.used,
        )
    if pathway == Pathway.EXON_DELETION:
        return _evaluate_exon_deletion(evidence)
    if pathway == Pathway.WHOLE_GENE_DELETION:
        return _evaluate_whole_gene_deletion(evidence)
    if pathway == Pathway.DUPLICATION:
        return _evaluate_duplication(evidence)
    if pathway == Pathway.INITIATION_CODON:
        return _evaluate_initiation_codon(evidence)
    # Defensive: all identifiers above are exhaustive.
    return _not_applicable("PVS1.UNSUPPORTED", "No supported PVS1 pathway was identified.")


class PVS1Engine:
    """Deterministic evaluator for one supplied biological feature set."""

    def evaluate(self, biological_features: Union[BiologicalFeatures, Mapping[str, Any]]) -> PVS1Decision:
        evidence = (
            biological_features.normalized()
            if isinstance(biological_features, BiologicalFeatures)
            else BiologicalFeatures.from_mapping(biological_features)
        )

        # A supplied VCEP snapshot is part of the caller's clinical input, not
        # a lazy fallback.  Validate every declaration before any generic
        # short-circuit so malformed specifications never disappear behind an
        # otherwise decisive branch.
        validate_specific_rule_declarations(evidence)

        # A VCEP declaration can resolve missing generic-routing facts, but it
        # cannot make two supplied variant-class assertions compatible. Check
        # only material cross-field conflicts before considering an override;
        # ordinary missing-data candidates remain eligible for VCEP handling.
        preflight = _identify_pathway_checked(evidence)
        if (
            preflight.candidate is not None
            and preflight.candidate.path.endswith("_CONFLICT")
        ):
            return self._finalize(preflight.candidate, evidence)

        # Specific, pinned VCEP guidance takes precedence even when the
        # generic variant-class annotation is incomplete.  This is evaluated
        # before `_identify_pathway`, not merely before the LoF-mechanism cap.
        specific_candidate = _vcep_candidate(evidence)
        if specific_candidate is not None:
            return self._finalize(
                specific_candidate,
                evidence,
                apply_mechanism_cap=False,
            )

        pathway_assessment = preflight
        if pathway_assessment.candidate is not None:
            return self._finalize(pathway_assessment.candidate, evidence)
        assert pathway_assessment.pathway is not None
        pathway = pathway_assessment.pathway

        # A known non-LoF mechanism is a global PVS1 exclusion. Other
        # mechanism uncertainty must not hide a known local branch exclusion.
        mechanism = _mechanism_assessment(evidence)
        if (
            mechanism.candidate is not None
            and mechanism.candidate.status == DecisionStatus.NOT_MET
        ):
            return self._finalize(mechanism.candidate, evidence)

        # Known branch exclusions are decisive even when positive LoF
        # eligibility is unknown. Only an activation is later gated/capped.
        candidate = _evaluate_generic_pathway(evidence, pathway)
        if candidate.status in {DecisionStatus.NOT_MET, DecisionStatus.NOT_APPLICABLE}:
            return self._finalize(candidate, evidence)

        if mechanism.candidate is not None:
            return self._finalize(mechanism.candidate, evidence)
        assert mechanism.cap is not None

        return self._finalize(candidate, evidence, mechanism)

    @staticmethod
    def _finalize(
        candidate: _Candidate,
        evidence: BiologicalFeatures,
        mechanism: Optional[_MechanismAssessment] = None,
        *,
        apply_mechanism_cap: bool = True,
    ) -> PVS1Decision:
        used = stable_unique(
            _pathway_selector_fields(candidate, evidence) + candidate.used
        )
        reason = candidate.reason
        strength = PVS1Strength.NONE
        if mechanism is not None:
            used = stable_unique(mechanism.used + used)
        if candidate.status == DecisionStatus.ACTIVATED:
            strength = candidate.base_strength
            if apply_mechanism_cap and mechanism is not None and mechanism.cap is not None:
                capped = _weaker(strength, mechanism.cap)
                if capped != strength:
                    reason += (
                        " The final strength is capped by the applicable ClinGen SVI 2018 "
                        "Table 1 gene–disease LoF-mechanism evidence tier."
                    )
                    strength = capped

        # RNA-mode decisions have two authorities: the 2023 guidance
        # establishes usable direct-splice evidence, while some resolved
        # products map to a 2018 Figure 1 strength row. Keep the primary
        # guideline's own locator here; `_generic_guideline_conditions`
        # records the second source/row when it was used.
        guideline = (
            candidate.guideline
            if (
                candidate.evidence_mode == "RNA"
                and candidate.guideline == RNA_SPLICE_GUIDELINE
            )
            else replace(candidate.guideline, source_locator=candidate.locator)
        )
        provenance_fields = list(used) + list(candidate.missing)
        for exclusion in candidate.exclusions:
            feature = exclusion.get("feature")
            feature = _audit_feature_name(feature or "")
            if feature in BiologicalFeatures.__dataclass_fields__:
                provenance_fields.append(feature)
        provenance: list[Provenance] = []
        for field in stable_unique(provenance_fields):
            canonical = _audit_feature_name(field)
            if canonical in BiologicalFeatures.__dataclass_fields__:
                provenance.extend(evidence.evidence(canonical).provenance)

        evidence_used: dict[str, Mapping[str, Any]] = {}
        for field in used:
            canonical = _audit_feature_name(field)
            if canonical in BiologicalFeatures.__dataclass_fields__:
                # Preserve a public alias in the key if it was the field the
                # branch reports, while using the canonical parsed evidence.
                evidence_used[field] = evidence.evidence(canonical).to_dict()
        guideline_conditions = (
            dict(candidate.guideline_conditions)
            if candidate.guideline_conditions
            else _generic_guideline_conditions(candidate, evidence, used)
        )
        return PVS1Decision(
            status=candidate.status,
            strength=strength,
            guideline=guideline,
            decision_path=candidate.path,
            reason=reason,
            features_used=stable_unique(used),
            features_missing=stable_unique(candidate.missing),
            exclusions_evaluated=tuple(candidate.exclusions),
            evidence_provenance=_unique_provenance(provenance),
            evidence_used=evidence_used,
            guideline_conditions=guideline_conditions,
            evidence_mode=candidate.evidence_mode,
        )


def evaluate_pvs1(
    biological_features: Union[BiologicalFeatures, Mapping[str, Any]],
) -> PVS1Decision:
    """Evaluate a supplied biological feature set under PVS1 only.

    The function does no annotation, database querying, or ACMG evidence
    combination.  It is deterministic and does not mutate the input mapping.
    """

    return PVS1Engine().evaluate(biological_features)
