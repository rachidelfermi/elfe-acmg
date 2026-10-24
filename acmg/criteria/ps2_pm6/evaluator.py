"""Adaptive ClinGen SVI PS2/PM6 evaluator.

This implements the ClinGen SVI v1.1 point table as a *pooled evidence
framework*, not as a generic feature score.  Every point comes from an
eligible, independently documented de novo observation; all other fields
establish whether that observation is eligible or which published table row
applies.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Iterable, Mapping, Optional, Sequence, Union

from acmg.criteria.pvs1.models import DecisionStatus, EvidenceStatus, EvidenceValue, GuidelineReference, Provenance

from .models import (
    DeNovoEvidence,
    DeNovoObservation,
    DeNovoStrength,
    GenotypeState,
    InheritancePattern,
    ObservationDecision,
    ObservationStatus,
    PS2PM6Decision,
    ParentageStatus,
    PhenotypeCategory,
    TestAdequacy,
    known,
    stable_unique,
    unique_provenance,
)
from .specifications import (
    DEFAULT_SVI_STRENGTH_THRESHOLDS,
    ScoringPolicyOverride,
    VCEPDeNovoRule,
    VCEPDeNovoSpecification,
    guideline_for_vcep_rule,
    parse_specifications,
    select_vcep_rule,
)


GENERAL_SVI_GUIDELINE = GuidelineReference(
    authority="ClinGen Sequence Variant Interpretation Working Group",
    document="Recommendation for de novo Criteria (PS2/PM6)",
    version="v1.1 (approved 2018-03-18; updated 2021-05-05)",
    source_url="https://clinicalgenome.org/site/assets/files/3461/svi_proposal_for_de_novo_criteria_v1_1.pdf",
    source_locator="Table 1, Table 2, and Additional considerations",
)

ACMG_AMP_GUIDELINE = GuidelineReference(
    authority="ACMG/AMP",
    document="Standards and guidelines for the interpretation and reporting of sequence variants",
    version="2015 (PMID: 25741868)",
    source_url="https://pubmed.ncbi.nlm.nih.gov/25741868/",
    source_locator="PS2 and PM6 definitions",
)


_DEFAULT_POINTS: Mapping[PhenotypeCategory, Mapping[ParentageStatus, float]] = {
    PhenotypeCategory.HIGHLY_SPECIFIC: {
        ParentageStatus.CONFIRMED: 2.0,
        ParentageStatus.ASSUMED: 1.0,
    },
    PhenotypeCategory.CONSISTENT: {
        ParentageStatus.CONFIRMED: 1.0,
        ParentageStatus.ASSUMED: 0.5,
    },
    PhenotypeCategory.HETEROGENEOUS: {
        ParentageStatus.CONFIRMED: 0.5,
        ParentageStatus.ASSUMED: 0.25,
    },
    PhenotypeCategory.INCONSISTENT: {
        ParentageStatus.CONFIRMED: 0.0,
        ParentageStatus.ASSUMED: 0.0,
    },
}

_DEFAULT_THRESHOLDS: Mapping[DeNovoStrength, float] = DEFAULT_SVI_STRENGTH_THRESHOLDS

_STRENGTH_ORDER = (
    DeNovoStrength.NONE,
    DeNovoStrength.SUPPORTING,
    DeNovoStrength.MODERATE,
    DeNovoStrength.STRONG,
    DeNovoStrength.VERY_STRONG,
)


def _exclusion(feature: str, outcome: str, detail: str = "") -> Mapping[str, str]:
    result: dict[str, str] = {"feature": feature, "outcome": outcome}
    if detail:
        result["detail"] = detail
    return result


def _case_id(observation: DeNovoObservation) -> Optional[str]:
    if observation.case_id.status == EvidenceStatus.KNOWN:
        return observation.case_id.value
    if observation.publication_case_id.status == EvidenceStatus.KNOWN:
        return observation.publication_case_id.value
    return None


def _case_identity_keys(observation: DeNovoObservation) -> tuple[str, ...]:
    """Stable identifiers used to prevent patient/report double counting.

    Local patient identifiers and literature/source-case identifiers are
    separate namespaces.  Equal-looking tokens in those namespaces are not
    adequate proof that two people are the same, so a cross-field collision
    cannot silently discard an independent observation.  Callers can link a
    report safely by supplying the relevant identifier on both records.
    """

    keys: list[str] = []
    for namespace, value in (
        ("case", observation.case_id),
        ("publication", observation.publication_case_id),
    ):
        if value.status == EvidenceStatus.KNOWN and isinstance(value.value, str):
            token = value.value.strip().upper()
            if token:
                keys.append(f"{namespace}:{token}")
    return tuple(dict.fromkeys(keys))


def _canonical_provenance(values: Iterable[Provenance]) -> tuple[Provenance, ...]:
    """Deduplicate source records in a deterministic order for merged cases."""

    return tuple(sorted(
        unique_provenance(values),
        key=lambda item: (
            item.source.casefold(),
            item.kind.value,
            item.identifier or "",
            item.method or "",
            item.note or "",
        ),
    ))


def _merge_evidence_values(values: Sequence[EvidenceValue[Any]]) -> EvidenceValue[Any]:
    """Merge duplicate reports without converting disagreement into a choice."""

    provenance = _canonical_provenance(
        item for value in values for item in value.provenance
    )
    if any(value.status == EvidenceStatus.CONFLICTING for value in values):
        return EvidenceValue(status=EvidenceStatus.CONFLICTING, provenance=provenance)
    known_values = [value.value for value in values if value.status == EvidenceStatus.KNOWN]
    if not known_values:
        return EvidenceValue(status=EvidenceStatus.UNKNOWN, provenance=provenance)
    baseline = known_values[0]
    if any(value != baseline for value in known_values[1:]):
        return EvidenceValue(status=EvidenceStatus.CONFLICTING, provenance=provenance)
    return EvidenceValue(status=EvidenceStatus.KNOWN, value=baseline, provenance=provenance)


def _merge_identifier_values(values: Sequence[EvidenceValue[Any]]) -> EvidenceValue[Any]:
    """Different local aliases can identify one report; retain one deterministically."""

    provenance = _canonical_provenance(
        item for value in values for item in value.provenance
    )
    if any(value.status == EvidenceStatus.CONFLICTING for value in values):
        return EvidenceValue(status=EvidenceStatus.CONFLICTING, provenance=provenance)
    identifiers = [
        value.value for value in values
        if value.status == EvidenceStatus.KNOWN and isinstance(value.value, str)
    ]
    if not identifiers:
        return EvidenceValue(status=EvidenceStatus.UNKNOWN, provenance=provenance)
    selected = min(identifiers, key=lambda value: (value.casefold(), value))
    return EvidenceValue(status=EvidenceStatus.KNOWN, value=selected, provenance=provenance)


def _merge_case_independence(values: Sequence[EvidenceValue[Any]]) -> EvidenceValue[bool]:
    """Merge independence assertions without letting a duplicate report veto itself.

    A report marked non-independent can merely identify itself as another
    report of the already-linked patient.  At the representative-group level,
    any positive assertion documents independence from *other* groups; a
    negative assertion matters only when no report asserts independence.
    """

    provenance = _canonical_provenance(
        item for value in values for item in value.provenance
    )
    if any(value.status == EvidenceStatus.CONFLICTING for value in values):
        return EvidenceValue(status=EvidenceStatus.CONFLICTING, provenance=provenance)
    known_values = [value.value for value in values if value.status == EvidenceStatus.KNOWN]
    if any(value is True for value in known_values):
        return EvidenceValue(status=EvidenceStatus.KNOWN, value=True, provenance=provenance)
    if any(value is False for value in known_values):
        return EvidenceValue(status=EvidenceStatus.KNOWN, value=False, provenance=provenance)
    return EvidenceValue(status=EvidenceStatus.UNKNOWN, provenance=provenance)


def _merge_duplicate_observations(
    observations: Sequence[DeNovoObservation],
) -> tuple[Mapping[int, DeNovoObservation], Mapping[int, int]]:
    """Coalesce reports linked by a stable patient/publication identifier.

    A paper and a laboratory record frequently contain complementary fields for
    the same individual.  Selecting whichever arrived first would make a
    decision input-order dependent, while simply counting both would inflate
    points.  We therefore merge agreement/complementary evidence into the
    earliest report and retain later reports as auditable duplicates.  Any
    genuinely incompatible fact becomes ``CONFLICTING``.
    """

    parent = list(range(len(observations)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    seen: dict[str, int] = {}
    for index, observation in enumerate(observations):
        for identity in _case_identity_keys(observation):
            if identity in seen:
                union(index, seen[identity])
            else:
                seen[identity] = index

    groups: dict[int, list[int]] = {}
    for index in range(len(observations)):
        groups.setdefault(find(index), []).append(index)

    representatives: dict[int, DeNovoObservation] = {}
    duplicate_of: dict[int, int] = {}
    for indexes in groups.values():
        representative = min(indexes)
        group = [observations[index] for index in indexes]
        replacements: dict[str, EvidenceValue[Any]] = {}
        for name in DeNovoObservation.__dataclass_fields__:
            values = [getattr(observation, name) for observation in group]
            if name in {"case_id", "publication_case_id"}:
                replacements[name] = _merge_identifier_values(values)
            elif name == "case_independence":
                replacements[name] = _merge_case_independence(values)
            else:
                replacements[name] = _merge_evidence_values(values)
        representatives[representative] = replace(group[0], **replacements)
        for index in indexes:
            if index != representative:
                duplicate_of[index] = representative
    return representatives, duplicate_of


def _merged_duplicate_decision(
    index: int,
    original_observation: DeNovoObservation,
    merged_observation: DeNovoObservation,
    representative: int,
) -> ObservationDecision:
    identity_features = tuple(
        name
        for name in ("case_id", "publication_case_id")
        if original_observation.evidence(name).status == EvidenceStatus.KNOWN
    )
    return ObservationDecision(
        index=index,
        case_id=_case_id(original_observation) or _case_id(merged_observation),
        duplicate_of=representative,
        status=ObservationStatus.DUPLICATE,
        parentage_status=None,
        phenotype_category=None,
        raw_points=0.0,
        contributed_points=0.0,
        features_used=identity_features,
        exclusions_evaluated=(_exclusion("case_id", "DUPLICATE"),),
        reason="This report was merged with an earlier record for the same stable case/publication identifier and was not counted twice.",
    )


def _is_inadequate_raw_metric(value: Any, *, variant_support: bool = False) -> bool:
    """Recognize an explicit inadequate/no-coverage conclusion only.

    Sequencing depth and allele balance have assay- and laboratory-specific
    thresholds.  A numerical value such as 8 is therefore preserved for audit
    but never silently compared with an invented cut-off.  Zero coverage or an
    explicit qualitative ``LOW``/``INADEQUATE`` conclusion is decisive.
    """

    if isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return value <= 0
    if isinstance(value, str):
        return value.strip().upper().replace("-", "_").replace(" ", "_") in {
            "LOW",
            "POOR",
            "INADEQUATE",
            "INSUFFICIENT",
            "NO_COVERAGE",
            "ZERO_COVERAGE",
            "FAIL",
        }
    if isinstance(value, Mapping):
        coverage_keys = {
            "DEPTH", "READ_DEPTH", "TOTAL_DEPTH", "COVERAGE", "TOTAL_COVERAGE",
            "CALLABLE_DEPTH", "PARENTAL_READ_DEPTH",
        }
        support_keys = {
            "ALT_READS", "ALTERNATE_READS", "VARIANT_READS", "SUPPORTING_READS",
            "READ_SUPPORT", "ALT_COUNT", "VARIANT_COUNT",
        }
        role_keys = {"MOTHER", "MATERNAL", "FATHER", "PATERNAL", "PARENT_1", "PARENT_2"}
        for key, item in value.items():
            token = str(key).strip().upper().replace("-", "_").replace(" ", "_")
            # A qualitative laboratory conclusion is meaningful wherever it
            # appears in a structured metric.
            if _has_explicit_inadequate_qualifier(item):
                return True
            numeric_measurement = token in coverage_keys or (
                variant_support and token in support_keys
            ) or (not variant_support and token in role_keys)
            if numeric_measurement:
                if isinstance(item, (int, float)) and not isinstance(item, bool) and item <= 0:
                    return True
                if isinstance(item, Mapping) and _is_inadequate_raw_metric(
                    item, variant_support=variant_support
                ):
                    return True
            elif isinstance(item, Mapping) and _is_inadequate_raw_metric(
                item, variant_support=variant_support
            ):
                return True
        return False
    if isinstance(value, (tuple, list)):
        # An unlabeled vector is not enough to distinguish alternate reads
        # from coverage.  Nested named measurements and explicit qualitative
        # assessments above remain decisive.
        return any(
            _is_inadequate_raw_metric(item, variant_support=variant_support)
            for item in value
            if isinstance(item, Mapping) or isinstance(item, str)
        )
    return False


def _has_explicit_inadequate_qualifier(value: Any) -> bool:
    """Find a qualitative laboratory conclusion without treating a number as one."""

    if isinstance(value, str):
        return value.strip().upper().replace("-", "_").replace(" ", "_") in {
            "LOW", "POOR", "INADEQUATE", "INSUFFICIENT", "NO_COVERAGE", "ZERO_COVERAGE", "FAIL",
        }
    if isinstance(value, Mapping):
        return any(_has_explicit_inadequate_qualifier(item) for item in value.values())
    if isinstance(value, (tuple, list)):
        return any(_has_explicit_inadequate_qualifier(item) for item in value)
    return False


def _is_relation_disproven(observation: DeNovoObservation) -> bool:
    relationship = observation.family_relationship_structure
    if relationship.status != EvidenceStatus.KNOWN or not isinstance(relationship.value, str):
        return False
    token = "_".join(relationship.value.upper().replace("-", " ").split())
    # These are deliberately finite, explicit normalized phrases rather than
    # a fuzzy free-text search.  A documented donor/non-biological parent is
    # decisive; an ambiguous narrative such as "possible donor" remains audit
    # context and must not be silently interpreted as a biological conclusion.
    return token in {
        "NON_BIOLOGICAL",
        "NONBIOLOGICAL",
        "NON_BIOLOGICAL_PARENT",
        "NONBIOLOGICAL_PARENT",
        "NON_BIOLOGICAL_MOTHER",
        "NONBIOLOGICAL_MOTHER",
        "MOTHER_NON_BIOLOGICAL",
        "MOTHER_NONBIOLOGICAL",
        "NON_BIOLOGICAL_FATHER",
        "NONBIOLOGICAL_FATHER",
        "FATHER_NON_BIOLOGICAL",
        "FATHER_NONBIOLOGICAL",
        "DONOR",
        "DONOR_EGG",
        "EGG_DONOR",
        "DONOR_SPERM",
        "SPERM_DONOR",
        "SURROGACY_CONFLICT",
    }


def _phenotype(
    observation: DeNovoObservation,
) -> tuple[Optional[PhenotypeCategory], tuple[str, ...], tuple[str, ...], Optional[str]]:
    """Resolve the published SVI phenotype row from direct or component data.

    ``phenotype_category`` is a convenience normalization field.  Upstream
    pipelines are also allowed to supply the task's native components:
    consistency, specificity, and genetic heterogeneity.  Those components
    describe the same Table 1 row; they are not an extra evidence score.
    """

    category = observation.phenotype_category
    specificity = observation.phenotype_specificity
    consistency = observation.phenotype_consistency
    heterogeneity = observation.phenotype_heterogeneity
    fields = {
        "phenotype_category": category,
        "phenotype_specificity": specificity,
        "phenotype_consistency": consistency,
        "phenotype_heterogeneity": heterogeneity,
    }
    conflicting = tuple(name for name, value in fields.items() if value.status == EvidenceStatus.CONFLICTING)
    if conflicting:
        return None, (), conflicting, "Phenotype evidence is conflicting."

    def used_fields(*names: str) -> tuple[str, ...]:
        return tuple(name for name in names if fields[name].status == EvidenceStatus.KNOWN)

    def conflict(names: tuple[str, ...], detail: str) -> tuple[None, tuple[str, ...], tuple[str, ...], str]:
        return None, used_fields(*names), names, detail

    def compatible_with_components(selected: PhenotypeCategory, names: tuple[str, ...]):
        if consistency.status == EvidenceStatus.KNOWN:
            expected_consistency = selected != PhenotypeCategory.INCONSISTENT
            if consistency.value != expected_consistency:
                return conflict(
                    names + ("phenotype_consistency",),
                    "Phenotype category conflicts with phenotype-consistency evidence.",
                )
        if heterogeneity.status == EvidenceStatus.KNOWN:
            expected_heterogeneity = selected == PhenotypeCategory.HETEROGENEOUS
            if heterogeneity.value != expected_heterogeneity:
                return conflict(
                    names + ("phenotype_heterogeneity",),
                    "Phenotype category conflicts with phenotype-heterogeneity evidence.",
                )
        return None

    if category.status == EvidenceStatus.KNOWN:
        selected = category.value
        allowed_specificity = {
            PhenotypeCategory.HIGHLY_SPECIFIC: {PhenotypeCategory.HIGHLY_SPECIFIC},
            PhenotypeCategory.CONSISTENT: {PhenotypeCategory.CONSISTENT},
            # A pipeline can state "consistent" and separately state high
            # genetic heterogeneity; together those are the heterogeneous row.
            PhenotypeCategory.HETEROGENEOUS: {
                PhenotypeCategory.CONSISTENT,
                PhenotypeCategory.HETEROGENEOUS,
            },
            PhenotypeCategory.INCONSISTENT: {PhenotypeCategory.INCONSISTENT},
        }
        if (
            specificity.status == EvidenceStatus.KNOWN
            and specificity.value not in allowed_specificity[selected]
        ):
            return conflict(
                ("phenotype_category", "phenotype_specificity"),
                "Direct phenotype category conflicts with phenotype-specificity evidence.",
            )
        component_problem = compatible_with_components(selected, ("phenotype_category",))
        if component_problem is not None:
            return component_problem
        return selected, used_fields(
            "phenotype_category", "phenotype_specificity", "phenotype_consistency", "phenotype_heterogeneity"
        ), (), None

    # An explicitly incompatible phenotype receives zero points.  It cannot
    # coexist with a positive-consistency or high-heterogeneity assertion.
    if consistency.status == EvidenceStatus.KNOWN and consistency.value is False:
        if heterogeneity.status == EvidenceStatus.KNOWN and heterogeneity.value:
            return conflict(
                ("phenotype_consistency", "phenotype_heterogeneity"),
                "An inconsistent phenotype cannot simultaneously be highly genetically heterogeneous for this disorder.",
            )
        if specificity.status == EvidenceStatus.KNOWN and specificity.value != PhenotypeCategory.INCONSISTENT:
            return conflict(
                ("phenotype_consistency", "phenotype_specificity"),
                "Phenotype-specificity and phenotype-consistency evidence conflict.",
            )
        return PhenotypeCategory.INCONSISTENT, used_fields(
            "phenotype_consistency", "phenotype_specificity", "phenotype_heterogeneity"
        ), (), None

    # Heterogeneity takes precedence over a non-highly-specific consistent
    # description.  It is the SVI "consistent but highly genetically
    # heterogeneous" row, not contradictory duplicate evidence.
    if heterogeneity.status == EvidenceStatus.KNOWN and heterogeneity.value:
        if specificity.status == EvidenceStatus.KNOWN and specificity.value in {
            PhenotypeCategory.HIGHLY_SPECIFIC,
            PhenotypeCategory.INCONSISTENT,
        }:
            return conflict(
                ("phenotype_specificity", "phenotype_heterogeneity"),
                "Highly specific or inconsistent phenotype evidence conflicts with high genetic heterogeneity.",
            )
        if (
            consistency.status == EvidenceStatus.KNOWN and consistency.value
        ) or (
            specificity.status == EvidenceStatus.KNOWN
            and specificity.value in {PhenotypeCategory.CONSISTENT, PhenotypeCategory.HETEROGENEOUS}
        ):
            return PhenotypeCategory.HETEROGENEOUS, used_fields(
                "phenotype_consistency", "phenotype_specificity", "phenotype_heterogeneity"
            ), (), None
        return None, used_fields("phenotype_heterogeneity"), ("phenotype_consistency",), None

    if specificity.status == EvidenceStatus.KNOWN:
        selected = specificity.value
        if selected == PhenotypeCategory.HETEROGENEOUS:
            if heterogeneity.status == EvidenceStatus.KNOWN and not heterogeneity.value:
                return conflict(
                    ("phenotype_specificity", "phenotype_heterogeneity"),
                    "Phenotype-specificity and phenotype-heterogeneity evidence conflict.",
                )
            return PhenotypeCategory.HETEROGENEOUS, used_fields(
                "phenotype_specificity", "phenotype_consistency", "phenotype_heterogeneity"
            ), (), None
        component_problem = compatible_with_components(selected, ("phenotype_specificity",))
        if component_problem is not None:
            return component_problem
        return selected, used_fields(
            "phenotype_specificity", "phenotype_consistency", "phenotype_heterogeneity"
        ), (), None

    if consistency.status == EvidenceStatus.KNOWN and consistency.value:
        # Positive consistency establishes the conservative (not highly
        # specific) row unless a separate specificity assertion says more.
        return PhenotypeCategory.CONSISTENT, used_fields(
            "phenotype_consistency", "phenotype_heterogeneity"
        ), (), None
    return None, (), ("phenotype_consistency",), None


def _proband_present(
    observation: DeNovoObservation,
) -> tuple[Optional[bool], tuple[str, ...], tuple[str, ...], Optional[str]]:
    present = observation.proband_variant_present
    genotype = observation.proband_genotype
    if present.status == EvidenceStatus.CONFLICTING or genotype.status == EvidenceStatus.CONFLICTING:
        return None, (), tuple(
            name
            for name, value in (("proband_variant_present", present), ("proband_genotype", genotype))
            if value.status == EvidenceStatus.CONFLICTING
        ), "Proband variant evidence is conflicting."
    from_boolean: Optional[bool] = present.value if present.status == EvidenceStatus.KNOWN else None
    from_genotype: Optional[bool]
    if genotype.status != EvidenceStatus.KNOWN:
        from_genotype = None
    elif genotype.value in {GenotypeState.VARIANT_PRESENT, GenotypeState.MOSAIC}:
        from_genotype = True
    elif genotype.value == GenotypeState.VARIANT_ABSENT:
        from_genotype = False
    else:
        # OTHER and NO_CALL are not a curated reference call.  Treating an
        # unfamiliar genotype label as absence could incorrectly suppress a
        # valid observation instead of reporting its missing activation fact.
        from_genotype = None
    if from_boolean is not None and from_genotype is not None and from_boolean != from_genotype:
        return None, ("proband_variant_present", "proband_genotype"), (
            "proband_variant_present",
            "proband_genotype",
        ), "Proband Boolean and genotype evidence conflict."
    result = from_boolean if from_boolean is not None else from_genotype
    if result is None:
        return None, (), ("proband_variant_present",), None
    used = tuple(
        name
        for name, available in (("proband_variant_present", from_boolean), ("proband_genotype", from_genotype))
        if available is not None
    )
    return result, used, (), None


def _parent_absence(
    observation: DeNovoObservation,
    role: str,
) -> tuple[Optional[bool], tuple[str, ...], tuple[str, ...], Optional[str], tuple[Mapping[str, str], ...]]:
    genotype = observation.maternal_genotype if role == "maternal" else observation.paternal_genotype
    field = f"{role}_genotype"
    if genotype.status == EvidenceStatus.CONFLICTING:
        return None, (), (field,), f"{role.capitalize()} genotype evidence is conflicting.", ()
    if genotype.status == EvidenceStatus.UNKNOWN or genotype.value == GenotypeState.NO_CALL:
        return None, (), (field,), None, ()
    if genotype.value in {GenotypeState.VARIANT_PRESENT, GenotypeState.MOSAIC}:
        detail = "mosaic parental call" if genotype.value == GenotypeState.MOSAIC else "parental variant call"
        return False, (field,), (), None, (_exclusion(field, "PRESENT", detail),)
    if genotype.value == GenotypeState.VARIANT_ABSENT:
        return True, (field,), (), None, ()
    return None, (field,), (field,), None, ()


def _parentage(
    observation: DeNovoObservation,
) -> tuple[Optional[ParentageStatus], tuple[str, ...], tuple[str, ...], Optional[str]]:
    maternity = observation.maternity_confirmed
    paternity = observation.paternity_confirmed
    disproven = observation.parentage_disproven
    method = observation.parentage_confirmation_method
    quality = observation.parentage_confirmation_quality
    asserted_status = observation.independent_case_parentage_status
    fields = {
        "maternity_confirmed": maternity,
        "paternity_confirmed": paternity,
        "parentage_disproven": disproven,
        "parentage_confirmation_method": method,
        "parentage_confirmation_quality": quality,
        "independent_case_parentage_status": asserted_status,
    }
    conflicts = tuple(name for name, value in fields.items() if value.status == EvidenceStatus.CONFLICTING)
    if conflicts:
        return None, (), conflicts, "Parental-relationship evidence is conflicting."
    if known(disproven):
        return ParentageStatus.INVALID, ("parentage_disproven",), (), "A biological parental relationship is disproven."
    if quality.status == EvidenceStatus.KNOWN and quality.value == TestAdequacy.INADEQUATE:
        # Inadequate parentage testing means the relationships are not
        # *confirmed* for PS2.  It does not undo otherwise tested-negative
        # parental genotypes or turn an apparent de novo event into a known
        # non-de-novo event, so the generic SVI pathway is PM6/assumed unless
        # a relationship has actually been disproven above.
        used = tuple(
            name
            for name, value in (
                ("maternity_confirmed", maternity),
                ("paternity_confirmed", paternity),
                ("independent_case_parentage_status", asserted_status),
            )
            if value.status == EvidenceStatus.KNOWN
        ) + ("parentage_confirmation_quality",)
        return ParentageStatus.ASSUMED, used, (), None
    if known(maternity) and known(paternity):
        if (
            asserted_status.status == EvidenceStatus.KNOWN
            and asserted_status.value == ParentageStatus.ASSUMED
        ):
            return None, (), ("independent_case_parentage_status",), (
                "The asserted assumed-parentage status conflicts with confirmed maternity and paternity."
            )
        used = ["maternity_confirmed", "paternity_confirmed"]
        if method.status == EvidenceStatus.KNOWN:
            used.append("parentage_confirmation_method")
        if quality.status == EvidenceStatus.KNOWN:
            used.append("parentage_confirmation_quality")
        if asserted_status.status == EvidenceStatus.KNOWN:
            used.append("independent_case_parentage_status")
        return ParentageStatus.CONFIRMED, tuple(used), (), None
    # The SVI distinction is confirmed versus unconfirmed parental
    # relationships.  A single confirmed parent is therefore an assumed-
    # parentage observation, not a partially confirmed PS2 observation.
    used = tuple(
        name
        for name, value in (
            ("maternity_confirmed", maternity),
            ("paternity_confirmed", paternity),
            ("independent_case_parentage_status", asserted_status),
        )
        if value.status == EvidenceStatus.KNOWN
    )
    if (
        asserted_status.status == EvidenceStatus.KNOWN
        and asserted_status.value == ParentageStatus.CONFIRMED
    ):
        return None, used, ("independent_case_parentage_status",), (
            "The asserted confirmed-parentage status is not supported by confirmed maternity and paternity."
        )
    return ParentageStatus.ASSUMED, used, (), None


def _strength_for(points: float, thresholds: Mapping[DeNovoStrength, float]) -> DeNovoStrength:
    if points <= 0:
        return DeNovoStrength.NONE
    result = DeNovoStrength.NONE
    for strength in (
        DeNovoStrength.SUPPORTING,
        DeNovoStrength.MODERATE,
        DeNovoStrength.STRONG,
        DeNovoStrength.VERY_STRONG,
    ):
        if points >= thresholds[strength]:
            result = strength
    return result


def _downgrade(strength: DeNovoStrength) -> DeNovoStrength:
    position = _STRENGTH_ORDER.index(strength)
    return _STRENGTH_ORDER[max(0, position - 1)]


def _point_value(
    category: PhenotypeCategory,
    parentage: ParentageStatus,
    override: Optional[ScoringPolicyOverride],
) -> float:
    if override is not None:
        value = override.points.get(category, {}).get(parentage)
        if value is not None:
            return value
    return _DEFAULT_POINTS[category][parentage]


@dataclass(frozen=True)
class _CaseAssessment:
    decision: ObservationDecision
    observation: DeNovoObservation
    ar_downgrade: bool = False
    rule: Optional[tuple[VCEPDeNovoSpecification, VCEPDeNovoRule]] = None


def _known_decisive_exclusion(
    index: int,
    observation: DeNovoObservation,
    case_id: Optional[str],
) -> Optional[_CaseAssessment]:
    """Return a conclusive generic exclusion before missing facts can mask it.

    This intentionally examines only explicit, self-contained negative facts.
    It never promotes missing evidence, and conflicts on a field remain for the
    main pathway to report as unresolved.  For example, a known paternal
    transmission or refuted gene–disease relationship stays decisive even if
    the maternal call is unavailable.
    """

    relation_disproven = _is_relation_disproven(observation)
    parentage_disproven = observation.parentage_disproven
    if relation_disproven and (
        parentage_disproven.status == EvidenceStatus.CONFLICTING
        or (
            parentage_disproven.status == EvidenceStatus.KNOWN
            and parentage_disproven.value is False
        )
    ):
        return _CaseAssessment(
            _undetermined_case(
                index,
                case_id,
                ("parentage_disproven", "family_relationship_structure"),
                "Structured parentage-disproval and family-relationship evidence conflict.",
            ),
            observation,
        )
    if known(parentage_disproven) or relation_disproven:
        field = "parentage_disproven" if known(observation.parentage_disproven) else "family_relationship_structure"
        return _CaseAssessment(
            _not_met_case(
                index,
                case_id,
                None,
                None,
                (field,),
                "The supplied family structure disproves a biological parental relationship.",
                (_exclusion(field, "NON_BIOLOGICAL"),),
            ),
            observation,
        )

    present, proband_used, _, proband_conflict = _proband_present(observation)
    if proband_conflict is None and present is False:
        return _CaseAssessment(
            _not_met_case(
                index,
                case_id,
                None,
                None,
                proband_used,
                "The variant is not present in the proband.",
                (_exclusion("proband_variant_present", "ABSENT"),),
            ),
            observation,
        )

    parents = (
        ("maternal", _parent_absence(observation, "maternal")),
        ("paternal", _parent_absence(observation, "paternal")),
    )
    for role, (absent, used, _, _, exclusions) in parents:
        if absent is False:
            return _CaseAssessment(
                _not_met_case(
                    index,
                    case_id,
                    None,
                    None,
                    used,
                    f"The variant is present in the {role}; this observation is not de novo.",
                    exclusions,
                ),
                observation,
            )

    quality = observation.variant_call_quality
    if quality.status == EvidenceStatus.KNOWN and quality.value == TestAdequacy.INADEQUATE:
        return _CaseAssessment(
            _not_met_case(
                index,
                case_id,
                None,
                None,
                ("variant_call_quality",),
                "Variant-call quality is inadequate for a de novo observation.",
                (_exclusion("variant_call_quality", "INADEQUATE"),),
            ),
            observation,
        )
    support = observation.proband_read_support
    if (
        support.status == EvidenceStatus.KNOWN
        and _is_inadequate_raw_metric(support.value, variant_support=True)
    ):
        return _CaseAssessment(
            _not_met_case(
                index,
                case_id,
                None,
                None,
                ("proband_read_support",),
                "Supplied proband read support explicitly indicates an inadequate variant call.",
                (_exclusion("proband_read_support", "INADEQUATE"),),
            ),
            observation,
        )
    for name in (
        "maternal_variant_absence_confidence",
        "paternal_variant_absence_confidence",
        "parental_testing_adequacy",
    ):
        assessment = observation.evidence(name)
        if assessment.status == EvidenceStatus.KNOWN and assessment.value == TestAdequacy.INADEQUATE:
            return _CaseAssessment(
                _not_met_case(
                    index,
                    case_id,
                    None,
                    None,
                    (name,),
                    "Parental testing is inadequate to establish variant absence.",
                    (_exclusion(name, "INADEQUATE"),),
                ),
                observation,
            )
    for name in ("parental_read_depth", "parental_allele_balance"):
        metric = observation.evidence(name)
        inadequate = (
            metric.status == EvidenceStatus.KNOWN
            and (
                (name == "parental_read_depth" and _is_inadequate_raw_metric(metric.value))
                or (
                    name == "parental_allele_balance"
                    and _has_explicit_inadequate_qualifier(metric.value)
                )
            )
        )
        if inadequate:
            return _CaseAssessment(
                _not_met_case(
                    index,
                    case_id,
                    None,
                    None,
                    (name,),
                    "Supplied parental test metrics explicitly indicate inadequate testing.",
                    (_exclusion(name, "INADEQUATE"),),
                ),
                observation,
            )
    mosaic = observation.parental_mosaicism_evidence
    if mosaic.status == EvidenceStatus.KNOWN and mosaic.value:
        return _CaseAssessment(
            _not_met_case(
                index,
                case_id,
                None,
                None,
                ("parental_mosaicism_evidence",),
                "Credible parental mosaicism prevents treating this observation as de novo.",
                (_exclusion("parental_mosaicism_evidence", "PRESENT"),),
            ),
            observation,
        )

    validity = observation.gene_disease_validity
    if validity.status == EvidenceStatus.KNOWN:
        if validity.value in {validity.value.REFUTED, validity.value.DISPUTED}:
            return _CaseAssessment(
                _not_applicable_case(
                    index,
                    case_id,
                    ("gene_disease_validity",),
                    "The supplied gene–disease association is refuted or disputed.",
                    (_exclusion("gene_disease_validity", "NOT_ESTABLISHED"),),
                ),
                observation,
            )
        if validity.value == validity.value.LIMITED:
            return _CaseAssessment(
                _not_met_case(
                    index,
                    case_id,
                    None,
                    None,
                    ("gene_disease_validity",),
                    "Limited gene–disease validity does not establish the disease context for generic PS2/PM6.",
                    (_exclusion("gene_disease_validity", "LIMITED"),),
                ),
                observation,
            )
    alternative = observation.alternative_molecular_diagnosis
    if alternative.status == EvidenceStatus.KNOWN and alternative.value:
        return _CaseAssessment(
            _not_met_case(
                index,
                case_id,
                None,
                None,
                ("alternative_molecular_diagnosis",),
                "A credible alternative molecular diagnosis explains the phenotype for this observation.",
                (_exclusion("alternative_molecular_diagnosis", "ALTERNATIVE_EXPLANATION"),),
            ),
            observation,
        )
    family_history = observation.family_history_consistent
    if family_history.status == EvidenceStatus.KNOWN and family_history.value is False:
        return _CaseAssessment(
            _not_met_case(
                index,
                case_id,
                None,
                None,
                ("family_history_consistent",),
                "The family history is not consistent with a de novo occurrence for this case.",
                (_exclusion("family_history_consistent", "INCONSISTENT"),),
            ),
            observation,
        )
    inheritance = observation.inheritance_pattern
    if (
        inheritance.status == EvidenceStatus.KNOWN
        and inheritance.value in {InheritancePattern.OTHER, InheritancePattern.MITOCHONDRIAL}
    ):
        return _CaseAssessment(
            _not_applicable_case(
                index,
                case_id,
                ("inheritance_pattern",),
                "This inheritance mechanism is outside the generic biparental PS2/PM6 pathway.",
                (_exclusion("inheritance_pattern", "OUTSIDE_GENERIC_SCOPE"),),
            ),
            observation,
        )
    category, phenotype_used, _, phenotype_conflict = _phenotype(observation)
    if phenotype_conflict is None and category == PhenotypeCategory.INCONSISTENT:
        return _CaseAssessment(
            _not_met_case(
                index,
                case_id,
                None,
                category,
                phenotype_used,
                "The phenotype is not consistent with the gene/disease association and receives zero SVI points.",
                (_exclusion("phenotype_category", "INCONSISTENT"),),
            ),
            observation,
        )
    return None


def _known_vcep_exclusion(
    index: int,
    observation: DeNovoObservation,
    specifications: Sequence[VCEPDeNovoSpecification],
    case_id: Optional[str],
) -> Optional[_CaseAssessment]:
    """Apply a resolvable, explicit VCEP exclusion before activation gates.

    A VCEP rule that only says an observation is out of scope does not need a
    complete de novo call to be decisive.  Only an actually matched exclusion
    is returned here; unresolved or scoring-only VCEP rows continue through
    the normal minimum-sufficient-evidence pathway.
    """

    if not specifications:
        return None
    parentage, parentage_used, _, parentage_problem = _parentage(observation)
    category, phenotype_used, _, phenotype_problem = _phenotype(observation)
    # Conflicting derived facts must not be coerced into a VCEP condition.  A
    # rule that does not reference them can still match safely.
    derived = {
        "parentage_status": None if parentage_problem else parentage,
        "phenotype_category": None if phenotype_problem else category,
    }
    selection = select_vcep_rule(specifications, observation, derived)
    if selection.matched is None:
        return None
    specification, rule = selection.matched
    if not rule.exclude:
        return None
    used = list(parentage_used) + list(phenotype_used)
    for field in ("gene_id", "disease_id", "inheritance_pattern"):
        if observation.evidence(field).status == EvidenceStatus.KNOWN:
            used.append(field)
    reason = rule.exclusion_reason or "A scoped VCEP rule excludes this observation."
    exclusions = (
        _exclusion("applicable_vcep_de_novo_rules", "VCEP_EXCLUSION", rule.source_locator),
    )
    if rule.exclusion_status == ObservationStatus.NOT_APPLICABLE:
        decision = _not_applicable_case(
            index,
            case_id,
            used,
            reason,
            exclusions,
            parentage=parentage,
            category=category,
            guideline_locator=rule.source_locator,
        )
    else:
        decision = _not_met_case(
            index,
            case_id,
            parentage,
            category,
            used,
            reason,
            exclusions,
            guideline_locator=rule.source_locator,
        )
    return _CaseAssessment(decision, observation, rule=selection.matched)


def _case(
    index: int,
    observation: DeNovoObservation,
    specifications: Sequence[VCEPDeNovoSpecification],
    *,
    multiple_observations: bool,
    duplicate: bool,
) -> _CaseAssessment:
    case_id = _case_id(observation)
    if duplicate:
        return _CaseAssessment(
            ObservationDecision(
                index=index,
                case_id=case_id,
                status=ObservationStatus.DUPLICATE,
                parentage_status=None,
                phenotype_category=None,
                raw_points=0.0,
                contributed_points=0.0,
                reason="This observation has the same stable case/publication identifier as an earlier observation and was not counted twice.",
                exclusions_evaluated=(_exclusion("case_id", "DUPLICATE"),),
            ),
            observation,
        )
    independence = observation.case_independence
    if multiple_observations:
        if independence.status == EvidenceStatus.CONFLICTING:
            return _CaseAssessment(_undetermined_case(index, case_id, ("case_independence",), "Case-independence evidence is conflicting."), observation)
        if independence.status == EvidenceStatus.KNOWN and independence.value is False:
            return _CaseAssessment(
                ObservationDecision(
                    index=index, case_id=case_id, status=ObservationStatus.DUPLICATE,
                    parentage_status=None, phenotype_category=None, raw_points=0.0, contributed_points=0.0,
                    reason="This observation is explicitly not independent and was not pooled.",
                    features_used=("case_independence",),
                    exclusions_evaluated=(_exclusion("case_independence", "NOT_INDEPENDENT"),),
                ), observation,
            )
        if independence.status == EvidenceStatus.UNKNOWN:
            return _CaseAssessment(_undetermined_case(index, case_id, ("case_independence",), "Each distinct observation must document independence before it can be pooled."), observation)

    vcep_exclusion = _known_vcep_exclusion(index, observation, specifications, case_id)
    if vcep_exclusion is not None:
        return vcep_exclusion

    # Sample identity is a prerequisite for trusting any purported trio
    # transmission.  Preserve its established unresolved outcome instead of
    # letting a parent call from a potentially misidentified sample decide the
    # case in the generic decisive-exclusion pass below.
    sample = observation.sample_identity_confirmation
    if sample.status == EvidenceStatus.CONFLICTING:
        return _CaseAssessment(_undetermined_case(index, case_id, ("sample_identity_confirmation",), "Sample-identity evidence is conflicting."), observation)
    if sample.status == EvidenceStatus.KNOWN and sample.value is False:
        return _CaseAssessment(_undetermined_case(index, case_id, ("sample_identity_confirmation",), "Sample identity is not confirmed and may invalidate the trio observation."), observation)
    relationship = observation.family_relationship_structure
    if relationship.status == EvidenceStatus.CONFLICTING:
        return _CaseAssessment(_undetermined_case(index, case_id, ("family_relationship_structure",), "Family-relationship structure evidence is conflicting."), observation)

    decisive = _known_decisive_exclusion(index, observation, case_id)
    if decisive is not None:
        return decisive

    # Credibility and identity exclusions are checked before table scoring.
    asserted_parentage = observation.independent_case_parentage_status
    if asserted_parentage.status == EvidenceStatus.CONFLICTING:
        return _CaseAssessment(
            _undetermined_case(
                index,
                case_id,
                ("independent_case_parentage_status",),
                "The independent-case parentage-status evidence is conflicting.",
            ),
            observation,
        )
    if known(observation.parentage_disproven) or _is_relation_disproven(observation):
        field = "parentage_disproven" if known(observation.parentage_disproven) else "family_relationship_structure"
        return _CaseAssessment(_not_met_case(index, case_id, None, None, (field,), "The supplied family structure disproves a biological parental relationship.", (_exclusion(field, "NON_BIOLOGICAL"),)), observation)

    present, used, missing, conflict = _proband_present(observation)
    if conflict:
        return _CaseAssessment(_undetermined_case(index, case_id, missing, conflict, used), observation)
    if present is None:
        return _CaseAssessment(_undetermined_case(index, case_id, missing, "The proband variant has not been established.", used), observation)
    if not present:
        return _CaseAssessment(_not_met_case(index, case_id, None, None, used, "The variant is not present in the proband.", (_exclusion("proband_variant_present", "ABSENT"),)), observation)

    proband_mosaic = observation.proband_mosaicism_evidence
    if proband_mosaic.status == EvidenceStatus.CONFLICTING:
        return _CaseAssessment(
            _undetermined_case(
                index,
                case_id,
                ("proband_mosaicism_evidence",),
                "Proband-mosaicism evidence is conflicting.",
                used,
            ),
            observation,
        )
    if proband_mosaic.status == EvidenceStatus.KNOWN:
        used += ("proband_mosaicism_evidence",)
        genotype = observation.proband_genotype
        if (
            genotype.status == EvidenceStatus.KNOWN
            and genotype.value == GenotypeState.MOSAIC
            and proband_mosaic.value is False
        ):
            return _CaseAssessment(
                _undetermined_case(
                    index,
                    case_id,
                    ("proband_genotype", "proband_mosaicism_evidence"),
                    "Proband genotype and mosaicism evidence conflict.",
                    used,
                ),
                observation,
            )

    mother_absent, mother_used, mother_missing, mother_conflict, mother_exclusions = _parent_absence(observation, "maternal")
    if mother_conflict:
        return _CaseAssessment(_undetermined_case(index, case_id, mother_missing, mother_conflict, used + mother_used), observation)
    if mother_absent is None:
        return _CaseAssessment(_undetermined_case(index, case_id, mother_missing, "Maternal testing does not establish variant absence.", used + mother_used), observation)
    if not mother_absent:
        return _CaseAssessment(_not_met_case(index, case_id, None, None, used + mother_used, "The variant is present in the mother; this observation is not de novo.", mother_exclusions), observation)

    father_absent, father_used, father_missing, father_conflict, father_exclusions = _parent_absence(observation, "paternal")
    if father_conflict:
        return _CaseAssessment(_undetermined_case(index, case_id, father_missing, father_conflict, used + mother_used + father_used), observation)
    if father_absent is None:
        return _CaseAssessment(_undetermined_case(index, case_id, father_missing, "Paternal testing does not establish variant absence.", used + mother_used + father_used), observation)
    if not father_absent:
        return _CaseAssessment(_not_met_case(index, case_id, None, None, used + mother_used + father_used, "The variant is present in the father; this observation is not de novo.", father_exclusions), observation)
    used += mother_used + father_used

    quality = observation.variant_call_quality
    if quality.status == EvidenceStatus.CONFLICTING:
        return _CaseAssessment(_undetermined_case(index, case_id, ("variant_call_quality",), "Variant-call quality evidence is conflicting.", used), observation)
    if quality.status == EvidenceStatus.UNKNOWN:
        return _CaseAssessment(_undetermined_case(index, case_id, ("variant_call_quality",), "Technical credibility of the variant observation is not established.", used), observation)
    used += ("variant_call_quality",)
    if quality.value == TestAdequacy.INADEQUATE:
        return _CaseAssessment(_not_met_case(index, case_id, None, None, used, "Variant-call quality is inadequate for a de novo observation.", (_exclusion("variant_call_quality", "INADEQUATE"),)), observation)
    proband_support = observation.proband_read_support
    if proband_support.status == EvidenceStatus.CONFLICTING:
        return _CaseAssessment(_undetermined_case(index, case_id, ("proband_read_support",), "Proband read-support evidence is conflicting.", used), observation)
    if proband_support.status == EvidenceStatus.KNOWN:
        used += ("proband_read_support",)
        if _is_inadequate_raw_metric(proband_support.value, variant_support=True):
            return _CaseAssessment(_not_met_case(index, case_id, None, None, used, "Supplied proband read support explicitly indicates an inadequate variant call.", (_exclusion("proband_read_support", "INADEQUATE"),)), observation)
    for name in (
        "maternal_variant_absence_confidence",
        "paternal_variant_absence_confidence",
        "parental_testing_adequacy",
    ):
        assessment = observation.evidence(name)
        if assessment.status == EvidenceStatus.CONFLICTING:
            return _CaseAssessment(_undetermined_case(index, case_id, (name,), f"{name} is conflicting.", used), observation)
        if assessment.status == EvidenceStatus.KNOWN:
            used += (name,)
            if assessment.value == TestAdequacy.INADEQUATE:
                return _CaseAssessment(_not_met_case(index, case_id, None, None, used, "Parental testing is inadequate to establish variant absence.", (_exclusion(name, "INADEQUATE"),)), observation)
    for name in ("parental_read_depth", "parental_allele_balance"):
        metric = observation.evidence(name)
        if metric.status == EvidenceStatus.CONFLICTING:
            return _CaseAssessment(_undetermined_case(index, case_id, (name,), f"{name} is conflicting.", used), observation)
        if metric.status == EvidenceStatus.KNOWN:
            used += (name,)
            # An allele balance is context, not a universally calibrated
            # adequacy threshold: a zero parental alternate fraction can be
            # exactly the expected negative call.  It can still carry an
            # explicit laboratory conclusion such as ``INADEQUATE``.
            inadequate = (
                name == "parental_read_depth" and _is_inadequate_raw_metric(metric.value)
            ) or (
                name == "parental_allele_balance" and _has_explicit_inadequate_qualifier(metric.value)
            )
            if inadequate:
                return _CaseAssessment(_not_met_case(index, case_id, None, None, used, "Supplied parental test metrics explicitly indicate inadequate testing.", (_exclusion(name, "INADEQUATE"),)), observation)
    mosaic = observation.parental_mosaicism_evidence
    if mosaic.status == EvidenceStatus.CONFLICTING:
        return _CaseAssessment(_undetermined_case(index, case_id, ("parental_mosaicism_evidence",), "Parental-mosaicism evidence is conflicting.", used), observation)
    if mosaic.status == EvidenceStatus.KNOWN:
        used += ("parental_mosaicism_evidence",)
        if mosaic.value:
            return _CaseAssessment(_not_met_case(index, case_id, None, None, used, "Credible parental mosaicism prevents treating this observation as de novo.", (_exclusion("parental_mosaicism_evidence", "PRESENT"),)), observation)

    validity = observation.gene_disease_validity
    if validity.status == EvidenceStatus.CONFLICTING:
        return _CaseAssessment(_undetermined_case(index, case_id, ("gene_disease_validity",), "Gene–disease validity evidence is conflicting.", used), observation)
    if validity.status == EvidenceStatus.UNKNOWN:
        return _CaseAssessment(_undetermined_case(index, case_id, ("gene_disease_validity",), "The gene–disease association is not established.", used), observation)
    used += ("gene_disease_validity",)
    if validity.value in {validity.value.REFUTED, validity.value.DISPUTED}:
        return _CaseAssessment(_not_applicable_case(index, case_id, used, "The supplied gene–disease association is refuted or disputed.", (_exclusion("gene_disease_validity", "NOT_ESTABLISHED"),)), observation)
    if validity.value == validity.value.LIMITED:
        return _CaseAssessment(_not_met_case(index, case_id, None, None, used, "Limited gene–disease validity does not establish the disease context for generic PS2/PM6.", (_exclusion("gene_disease_validity", "LIMITED"),)), observation)

    alternative = observation.alternative_molecular_diagnosis
    if alternative.status == EvidenceStatus.CONFLICTING:
        return _CaseAssessment(_undetermined_case(index, case_id, ("alternative_molecular_diagnosis",), "Alternative-diagnosis evidence is conflicting.", used), observation)
    if alternative.status == EvidenceStatus.KNOWN:
        used += ("alternative_molecular_diagnosis",)
        if alternative.value:
            return _CaseAssessment(_not_met_case(index, case_id, None, None, used, "A credible alternative molecular diagnosis explains the phenotype for this observation.", (_exclusion("alternative_molecular_diagnosis", "ALTERNATIVE_EXPLANATION"),)), observation)

    family_history = observation.family_history_consistent
    if family_history.status == EvidenceStatus.CONFLICTING:
        return _CaseAssessment(_undetermined_case(index, case_id, ("family_history_consistent",), "Family-history evidence is conflicting."), observation)
    if family_history.status == EvidenceStatus.KNOWN:
        used += ("family_history_consistent",)
        if family_history.value is False:
            return _CaseAssessment(_not_met_case(index, case_id, None, None, used, "The family history is not consistent with a de novo occurrence for this case.", (_exclusion("family_history_consistent", "INCONSISTENT"),)), observation)

    inheritance = observation.inheritance_pattern
    if inheritance.status == EvidenceStatus.CONFLICTING:
        return _CaseAssessment(_undetermined_case(index, case_id, ("inheritance_pattern",), "Inheritance-mechanism evidence is conflicting.", used), observation)
    # SVI v1.1 makes inheritance-specific adjustments only when that context
    # is established (notably the AR one-level reduction).  It is not a
    # universal prerequisite for the generic point table, so its omission
    # cannot block an otherwise minimum-sufficient PS2/PM6 observation.
    if inheritance.status == EvidenceStatus.KNOWN:
        used += ("inheritance_pattern",)
    if inheritance.status == EvidenceStatus.KNOWN and inheritance.value in {InheritancePattern.OTHER, InheritancePattern.MITOCHONDRIAL}:
        return _CaseAssessment(_not_applicable_case(index, case_id, used, "This inheritance mechanism is outside the generic biparental PS2/PM6 pathway.", (_exclusion("inheritance_pattern", "OUTSIDE_GENERIC_SCOPE"),)), observation)

    category, phenotype_used, phenotype_missing, phenotype_conflict = _phenotype(observation)
    if phenotype_conflict:
        return _CaseAssessment(_undetermined_case(index, case_id, phenotype_missing, phenotype_conflict, used + phenotype_used), observation)
    if category is None:
        return _CaseAssessment(_undetermined_case(index, case_id, phenotype_missing, "Phenotype specificity cannot be assigned to a published SVI table row.", used + phenotype_used), observation)
    used += phenotype_used
    if category == PhenotypeCategory.INCONSISTENT:
        return _CaseAssessment(_not_met_case(index, case_id, None, category, used, "The phenotype is not consistent with the gene/disease association and receives zero SVI points.", (_exclusion("phenotype_category", "INCONSISTENT"),)), observation)

    parentage, parentage_used, parentage_missing, parentage_problem = _parentage(observation)
    if parentage is None:
        return _CaseAssessment(_undetermined_case(index, case_id, parentage_missing, parentage_problem or "Parentage status is unresolved.", used + parentage_used), observation)
    used += parentage_used
    if parentage == ParentageStatus.INVALID:
        return _CaseAssessment(_not_met_case(index, case_id, None, category, used, parentage_problem or "Parentage is invalid.", (_exclusion("parentage_confirmation_quality", "INADEQUATE"),)), observation)
    germline = observation.apparent_germline_mosaicism
    if germline.status == EvidenceStatus.CONFLICTING:
        return _CaseAssessment(_undetermined_case(index, case_id, ("apparent_germline_mosaicism",), "Apparent germline-mosaicism evidence is conflicting.", used), observation)
    if germline.status == EvidenceStatus.KNOWN:
        used += ("apparent_germline_mosaicism",)
        if germline.value and parentage != ParentageStatus.CONFIRMED:
            return _CaseAssessment(_not_met_case(index, case_id, parentage, category, used, "Apparent germline mosaicism requires confirmed parental relationships before de novo criteria apply.", (_exclusion("apparent_germline_mosaicism", "PARENTAGE_NOT_CONFIRMED"),)), observation)

    selection = select_vcep_rule(
        specifications,
        observation,
        {"parentage_status": parentage, "phenotype_category": category},
    )
    if selection.conflict:
        return _CaseAssessment(_undetermined_case(index, case_id, (), selection.conflict, used), observation)
    if selection.missing:
        return _CaseAssessment(_undetermined_case(index, case_id, selection.missing, "A potentially applicable VCEP rule cannot be resolved safely.", used), observation)
    rule_pair = selection.matched
    if rule_pair is not None:
        _, rule = rule_pair
        if rule.exclude:
            reason = rule.exclusion_reason or "A scoped VCEP rule excludes this observation."
            exclusions = (_exclusion("applicable_vcep_de_novo_rules", "VCEP_EXCLUSION", rule.source_locator),)
            if rule.exclusion_status == ObservationStatus.NOT_APPLICABLE:
                decision = _not_applicable_case(
                    index,
                    case_id,
                    used,
                    reason,
                    exclusions,
                    parentage=parentage,
                    category=category,
                    guideline_locator=rule.source_locator,
                )
            else:
                decision = _not_met_case(
                    index,
                    case_id,
                    parentage,
                    category,
                    used,
                    reason,
                    exclusions,
                    guideline_locator=rule.source_locator,
                )
            return _CaseAssessment(decision, observation, rule=rule_pair)

    ar_downgrade = False
    if inheritance.status == EvidenceStatus.KNOWN and inheritance.value == InheritancePattern.AUTOSOMAL_RECESSIVE:
        second = observation.additional_pathogenic_variant_present
        if second.status == EvidenceStatus.CONFLICTING:
            return _CaseAssessment(_undetermined_case(index, case_id, ("additional_pathogenic_variant_present",), "Additional-allele evidence is conflicting for an autosomal-recessive disorder.", used), observation, rule=rule_pair)
        if second.status == EvidenceStatus.UNKNOWN:
            return _CaseAssessment(_undetermined_case(index, case_id, ("additional_pathogenic_variant_present",), "The autosomal-recessive one-level adjustment cannot be determined without the additional-allele assessment.", used), observation, rule=rule_pair)
        used += ("additional_pathogenic_variant_present",)
        ar_downgrade = second.value is False

    override = rule_pair[1].scoring if rule_pair is not None else None
    points = _point_value(category, parentage, override)
    locator = rule_pair[1].source_locator if rule_pair is not None else GENERAL_SVI_GUIDELINE.source_locator
    return _CaseAssessment(
        ObservationDecision(
            index=index,
            case_id=case_id,
            status=ObservationStatus.ELIGIBLE,
            parentage_status=parentage,
            phenotype_category=category,
            raw_points=points,
            contributed_points=points,
            features_used=stable_unique(used),
            guideline_locator=locator,
            reason="Eligible de novo observation scored with the applicable ClinGen SVI/VCEP point-table row.",
        ),
        observation,
        ar_downgrade=ar_downgrade,
        rule=rule_pair,
    )


def _undetermined_case(
    index: int,
    case_id: Optional[str],
    missing: Iterable[str],
    reason: str,
    used: Iterable[str] = (),
) -> ObservationDecision:
    return ObservationDecision(
        index=index, case_id=case_id, status=ObservationStatus.UNDETERMINED,
        parentage_status=None, phenotype_category=None, raw_points=0.0, contributed_points=0.0,
        features_used=stable_unique(used), missing_critical_features=stable_unique(missing), reason=reason,
    )


def _not_met_case(
    index: int,
    case_id: Optional[str],
    parentage: Optional[ParentageStatus],
    category: Optional[PhenotypeCategory],
    used: Iterable[str],
    reason: str,
    exclusions: Iterable[Mapping[str, str]] = (),
    *,
    guideline_locator: Optional[str] = None,
) -> ObservationDecision:
    return ObservationDecision(
        index=index, case_id=case_id, status=ObservationStatus.NOT_MET,
        parentage_status=parentage, phenotype_category=category, raw_points=0.0, contributed_points=0.0,
        features_used=stable_unique(used), exclusions_evaluated=tuple(exclusions),
        guideline_locator=guideline_locator, reason=reason,
    )


def _not_applicable_case(
    index: int,
    case_id: Optional[str],
    used: Iterable[str],
    reason: str,
    exclusions: Iterable[Mapping[str, str]] = (),
    *,
    parentage: Optional[ParentageStatus] = None,
    category: Optional[PhenotypeCategory] = None,
    guideline_locator: Optional[str] = None,
) -> ObservationDecision:
    return ObservationDecision(
        index=index, case_id=case_id, status=ObservationStatus.NOT_APPLICABLE,
        parentage_status=parentage, phenotype_category=category, raw_points=0.0, contributed_points=0.0,
        features_used=stable_unique(used), exclusions_evaluated=tuple(exclusions),
        guideline_locator=guideline_locator, reason=reason,
    )


def _final_policy(
    assessments: Sequence[_CaseAssessment],
) -> tuple[Optional[Mapping[DeNovoStrength, float]], Optional[float], Optional[str], Optional[str]]:
    """Return thresholds, heterogeneous cap, label override, or a conflict reason."""

    threshold_sets: list[Mapping[DeNovoStrength, float]] = []
    caps: list[float] = []
    labels: list[str] = []
    for assessment in assessments:
        if assessment.rule is None or assessment.decision.status != ObservationStatus.ELIGIBLE:
            continue
        override = assessment.rule[1].scoring
        if override.thresholds:
            threshold_sets.append(override.thresholds)
        if override.heterogeneous_points_cap is not None:
            caps.append(override.heterogeneous_points_cap)
        if override.assigned_criterion is not None:
            labels.append(override.assigned_criterion)

    def one_value(values: Sequence[Any], label: str) -> tuple[Optional[Any], Optional[str]]:
        if not values:
            return None, None
        baseline = values[0]
        if any(value != baseline for value in values[1:]):
            return None, f"Multiple VCEP rules specify conflicting {label}."
        return baseline, None

    threshold_override, conflict = one_value(threshold_sets, "strength thresholds")
    if conflict:
        return None, None, None, conflict
    cap_override, conflict = one_value(caps, "heterogeneous-phenotype caps")
    if conflict:
        return None, None, None, conflict
    label_override, conflict = one_value(labels, "criterion labels")
    if conflict:
        return None, None, None, conflict
    thresholds = dict(_DEFAULT_THRESHOLDS)
    if threshold_override:
        thresholds.update(threshold_override)
    return thresholds, cap_override if cap_override is not None else 1.0, label_override, None


def _qualified_feature(index: int, count: int, field: str) -> str:
    return field if count == 1 else f"observations[{index}].{field}"


def _pooled_context_problem(
    assessments: Sequence[_CaseAssessment],
) -> tuple[tuple[str, ...], Optional[str]]:
    """Ensure a generic pooled assertion stays within one disease context.

    The SVI point table pools independent observations for the same
    gene/disease/inheritance setting, not arbitrary observations that happen
    to be present in one API request.  Completely omitted context remains
    nonblocking for backward-compatible generic use; once a context field is
    supplied for any eligible representative, every pooled representative must
    agree on it.
    """

    eligible = [
        item for item in assessments
        if item.decision.status == ObservationStatus.ELIGIBLE
    ]
    if len(eligible) < 2:
        return (), None
    unresolved: list[str] = []
    for field in ("gene_id", "disease_id", "inheritance_pattern"):
        values = [item.observation.evidence(field) for item in eligible]
        if any(value.status == EvidenceStatus.CONFLICTING for value in values):
            unresolved.append(field)
            continue
        known_values = [value.value for value in values if value.status == EvidenceStatus.KNOWN]
        if not known_values:
            continue
        if len(known_values) != len(values):
            unresolved.append(field)
            continue
        baseline = known_values[0]
        if field in {"gene_id", "disease_id"}:
            equal = lambda value: (
                isinstance(value, str)
                and isinstance(baseline, str)
                and value.strip().casefold() == baseline.strip().casefold()
            )
        else:
            equal = lambda value: value == baseline
        if any(not equal(value) for value in known_values[1:]):
            unresolved.append(field)
    if unresolved:
        return (
            tuple(dict.fromkeys(unresolved)),
            "Eligible observations do not establish one consistent gene, disease, and inheritance context for pooled PS2/PM6 evidence.",
        )
    return (), None


class PS2PM6Engine:
    """Deterministic shared evaluator for the mutually exclusive PS2/PM6 pair."""

    def evaluate(self, biological_features: Union[DeNovoEvidence, Mapping[str, Any]]) -> PS2PM6Decision:
        evidence = (
            biological_features
            if isinstance(biological_features, DeNovoEvidence)
            else DeNovoEvidence.from_mapping(biological_features)
        )
        vcep_raw = evidence.applicable_vcep_de_novo_rules
        if vcep_raw.status == EvidenceStatus.CONFLICTING:
            return self._global_undetermined(evidence, ("applicable_vcep_de_novo_rules",), "VCEP specification evidence is conflicting.")
        specifications: tuple[VCEPDeNovoSpecification, ...] = ()
        if vcep_raw.status == EvidenceStatus.KNOWN:
            specifications = parse_specifications(vcep_raw.value or ())

        merged_observations, duplicate_of = _merge_duplicate_observations(evidence.observations)
        multiple_representative_groups = len(merged_observations) > 1
        assessments: list[_CaseAssessment] = []
        for index, original_observation in enumerate(evidence.observations):
            if index in duplicate_of:
                representative = duplicate_of[index]
                merged_observation = merged_observations[representative]
                assessments.append(
                    _CaseAssessment(
                        _merged_duplicate_decision(
                            index,
                            original_observation,
                            merged_observation,
                            representative,
                        ),
                        # Use merged evidence for deterministic provenance;
                        # the decision itself retains the original record ID
                        # and representative index for audit.
                        merged_observation,
                    )
                )
                continue
            observation = merged_observations[index]
            assessment = _case(
                index,
                observation,
                specifications,
                multiple_observations=multiple_representative_groups,
                duplicate=False,
            )
            # Score from the merged evidence group, but retain the submitted
            # report's own display identifier on its audit row.  Canonical
            # merged identity is available through the linked duplicate rows
            # and merged provenance without rewriting a laboratory's record.
            original_case_id = _case_id(original_observation)
            if original_case_id is not None and assessment.decision.case_id != original_case_id:
                assessment = replace(
                    assessment,
                    decision=replace(assessment.decision, case_id=original_case_id),
                )
            assessments.append(assessment)

        # A pooled PS2/PM6 assertion describes one variant in one
        # gene/disease context.  Do not let an override for one VCEP scope
        # change thresholds or provenance for a generic/different-scope case.
        eligible_assessments = [
            item for item in assessments
            if item.decision.status == ObservationStatus.ELIGIBLE
        ]
        context_missing, context_problem = _pooled_context_problem(assessments)
        if context_problem is not None:
            return self._global_undetermined(
                evidence,
                context_missing,
                context_problem,
                assessments,
            )
        scoped_specs = [item.rule[0] for item in eligible_assessments if item.rule is not None]
        if scoped_specs:
            first_specification = scoped_specs[0]
            has_generic_eligible_case = any(item.rule is None for item in eligible_assessments)
            if has_generic_eligible_case or any(
                specification != first_specification for specification in scoped_specs[1:]
            ):
                return self._global_undetermined(
                    evidence,
                    (),
                    "Eligible observations span generic or multiple VCEP gene/disease scopes and cannot be pooled under one PS2/PM6 policy.",
                    assessments,
                )

        thresholds, heterogeneous_cap, label_override, policy_conflict = _final_policy(assessments)
        if policy_conflict:
            return self._global_undetermined(evidence, (), policy_conflict, assessments)
        assert thresholds is not None and heterogeneous_cap is not None

        # Apply SVI's one-point aggregate cap only to the high-heterogeneity
        # row.  The partial final contribution is retained per case for audit.
        heterogeneous_used = 0.0
        finalized: list[_CaseAssessment] = []
        for assessment in assessments:
            decision = assessment.decision
            contributed = decision.raw_points
            if decision.status == ObservationStatus.ELIGIBLE and decision.phenotype_category == PhenotypeCategory.HETEROGENEOUS:
                remaining = max(0.0, heterogeneous_cap - heterogeneous_used)
                contributed = min(contributed, remaining)
                heterogeneous_used += contributed
            if contributed != decision.contributed_points:
                decision = replace(
                    decision,
                    contributed_points=contributed,
                    reason=decision.reason + " The applicable heterogeneous-phenotype aggregate cap limits this case's contribution.",
                )
            finalized.append(replace(assessment, decision=decision))
        assessments = finalized
        # A case remains a qualifying confirmed/assumed observation even when
        # the published heterogeneous-phenotype *aggregate* cap trims its
        # final numeric contribution.  Use raw eligibility for the one-code
        # PS2 precedence rule and reported case counts, and capped values only
        # for the pooled strength calculation.
        eligible_cases = [item for item in assessments if item.decision.status == ObservationStatus.ELIGIBLE and item.decision.raw_points > 0]
        contributing_cases = [item for item in eligible_cases if item.decision.contributed_points > 0]
        # Do not round before threshold comparison: a VCEP can publish a
        # point allocation with more than four decimal places, and rounding a
        # value below 0.5 up to 0.5 would create unsupported Supporting data.
        total_points = sum(item.decision.contributed_points for item in contributing_cases)
        confirmed = sum(1 for item in eligible_cases if item.decision.parentage_status == ParentageStatus.CONFIRMED)
        assumed = sum(1 for item in eligible_cases if item.decision.parentage_status == ParentageStatus.ASSUMED)
        strength = _strength_for(total_points, thresholds)

        # SVI instructs a one-level reduction for an AR disease without a
        # second pathogenic/likely pathogenic allele.  It is a strength
        # adjustment after pooling, not a fabricated change to point values.
        if strength != DeNovoStrength.NONE and any(item.ar_downgrade for item in eligible_cases):
            strength = _downgrade(strength)

        # A supplied case-count summary cannot create cases.  If it contradicts
        # explicitly independent input it makes the evidence set unresolved.
        count = evidence.number_of_independent_de_novo_cases
        if count.status == EvidenceStatus.CONFLICTING:
            return self._global_undetermined(evidence, ("number_of_independent_de_novo_cases",), "The independent-case count is conflicting.", assessments)
        if count.status == EvidenceStatus.KNOWN and count.value is not None and count.value < len(eligible_cases):
            return self._global_undetermined(evidence, ("number_of_independent_de_novo_cases",), "The supplied independent-case count conflicts with pooled eligible observations.", assessments)

        inferred_criterion = "PS2" if confirmed else "PM6"
        if label_override is not None and label_override != inferred_criterion:
            return self._global_undetermined(
                evidence,
                (),
                f"The applicable VCEP criterion label ({label_override}) conflicts with the qualifying {inferred_criterion} parentage evidence.",
                assessments,
            )

        if strength != DeNovoStrength.NONE:
            assigned = inferred_criterion
            guideline = self._guideline_for(assessments)
            ar_note = " The final strength was reduced one level for autosomal-recessive disease without an additional pathogenic/likely pathogenic allele." if any(item.ar_downgrade for item in eligible_cases) else ""
            return self._decision(
                evidence=evidence,
                assessments=assessments,
                status=DecisionStatus.ACTIVATED,
                assigned=assigned,
                strength=strength,
                points=total_points,
                confirmed=confirmed,
                assumed=assumed,
                guideline=guideline,
                path="SVI2021.PS2_PM6.POOLED",
                reason=(
                    f"{total_points:g} eligible ClinGen SVI de novo point(s) meet the {strength.value} threshold; "
                    f"{assigned} is the sole assigned criterion."
                    + ar_note
                ),
            )

        unresolved = [item.decision for item in assessments if item.decision.status == ObservationStatus.UNDETERMINED]
        if unresolved:
            missing = stable_unique(field for item in unresolved for field in item.missing_critical_features)
            return self._decision(
                evidence=evidence,
                assessments=assessments,
                status=DecisionStatus.UNDETERMINED,
                assigned=None,
                strength=DeNovoStrength.NONE,
                points=total_points,
                confirmed=confirmed,
                assumed=assumed,
                guideline=self._guideline_for(assessments),
                path="SVI2021.PS2_PM6.UNRESOLVED",
                reason="No qualifying point threshold is established because one or more potentially relevant observations are unresolved.",
                extra_missing=missing,
            )

        # DUPLICATE is audit metadata rather than a substantive clinical
        # outcome.  It must not turn an otherwise out-of-scope VCEP decision
        # into NOT_MET merely because the same report was supplied twice.
        statuses = {
            item.decision.status
            for item in assessments
            if item.decision.status != ObservationStatus.DUPLICATE
        }
        status = DecisionStatus.NOT_APPLICABLE if statuses == {ObservationStatus.NOT_APPLICABLE} else DecisionStatus.NOT_MET
        return self._decision(
            evidence=evidence,
            assessments=assessments,
            status=status,
            assigned=None,
            strength=DeNovoStrength.NONE,
            points=total_points,
            confirmed=confirmed,
            assumed=assumed,
            guideline=self._guideline_for(assessments),
            path="SVI2021.PS2_PM6.NO_THRESHOLD",
            reason="No eligible de novo observations reach the ClinGen SVI Supporting threshold.",
        )

    def _global_undetermined(
        self,
        evidence: DeNovoEvidence,
        missing: Iterable[str],
        reason: str,
        assessments: Sequence[_CaseAssessment] = (),
    ) -> PS2PM6Decision:
        return self._decision(
            evidence=evidence,
            assessments=assessments,
            status=DecisionStatus.UNDETERMINED,
            assigned=None,
            strength=DeNovoStrength.NONE,
            points=0.0,
            confirmed=0,
            assumed=0,
            guideline=self._guideline_for(assessments) if assessments else GENERAL_SVI_GUIDELINE,
            path="SVI2021.PS2_PM6.UNRESOLVED",
            reason=reason,
            extra_missing=tuple(missing),
        )

    @staticmethod
    def _guideline_for(assessments: Sequence[_CaseAssessment]) -> GuidelineReference:
        # A scoped VCEP rule remains the governing source when it excludes a
        # case or declares it out of scope; filtering only eligible scoring
        # rows would falsely report the generic SVI document for that outcome.
        eligible_rules = [
            item.rule for item in assessments
            if item.rule is not None and item.decision.status == ObservationStatus.ELIGIBLE
        ]
        rules = eligible_rules or [item.rule for item in assessments if item.rule is not None]
        if rules:
            first = rules[0]
            assert first is not None
            if all(rule is not None and rule[0] == first[0] for rule in rules):
                specification = first[0]
                locators = tuple(dict.fromkeys(rule[1].source_locator for rule in rules if rule is not None))
                # A single VCEP document may have separate phenotype rows.
                # Keep the document identity at the top level and make every
                # applied source row visible rather than silently falling back
                # to the generic SVI reference.
                if len(locators) > 1:
                    return GuidelineReference(
                        authority=specification.authority,
                        document=specification.name,
                        version=specification.version,
                        source_url=specification.source_url,
                        source_locator="; ".join(locators),
                        rule_id=specification.specification_id,
                    )
                return guideline_for_vcep_rule(*first)
        return GENERAL_SVI_GUIDELINE

    @staticmethod
    def _decision(
        *,
        evidence: DeNovoEvidence,
        assessments: Sequence[_CaseAssessment],
        status: DecisionStatus,
        assigned: Optional[str],
        strength: DeNovoStrength,
        points: float,
        confirmed: int,
        assumed: int,
        guideline: GuidelineReference,
        path: str,
        reason: str,
        extra_missing: Iterable[str] = (),
    ) -> PS2PM6Decision:
        count = len(evidence.observations)
        observations = tuple(item.decision for item in assessments)
        used: list[str] = []
        missing: list[str] = list(extra_missing)
        exclusions: list[Mapping[str, str]] = []
        provenance: list[Provenance] = list(evidence.evidence_provenance)
        for assessment in assessments:
            decision = assessment.decision
            index = decision.index
            local_fields = set(decision.features_used) | set(decision.missing_critical_features)
            for exclusion in decision.exclusions_evaluated:
                field = exclusion.get("feature")
                if field:
                    local_fields.add(field)
            used.extend(_qualified_feature(index, count, field) for field in decision.features_used)
            missing.extend(_qualified_feature(index, count, field) for field in decision.missing_critical_features)
            exclusions.extend(decision.exclusions_evaluated)
            for field in local_fields:
                if field in DeNovoObservation.__dataclass_fields__:
                    provenance.extend(assessment.observation.evidence(field).provenance)
        if evidence.applicable_vcep_de_novo_rules.status == EvidenceStatus.KNOWN:
            provenance.extend(evidence.applicable_vcep_de_novo_rules.provenance)
        if evidence.number_of_independent_de_novo_cases.status == EvidenceStatus.KNOWN:
            provenance.extend(evidence.number_of_independent_de_novo_cases.provenance)
        return PS2PM6Decision(
            status=status,
            assigned_criterion=assigned,
            strength=strength,
            total_points=points,
            confirmed_de_novo_cases=confirmed,
            assumed_de_novo_cases=assumed,
            observations=observations,
            guideline=guideline,
            decision_path=path,
            reason=reason,
            features_used=stable_unique(used),
            missing_critical_features=stable_unique(missing),
            exclusions_evaluated=tuple(exclusions),
            evidence_provenance=_canonical_provenance(provenance),
        )


def evaluate_ps2_pm6(
    biological_features: Union[DeNovoEvidence, Mapping[str, Any]],
) -> PS2PM6Decision:
    """Evaluate pre-extracted evidence under the shared PS2/PM6 framework.

    The function is deterministic, does not mutate its input, and never
    performs annotation, trio calling, phenotype extraction, or full ACMG
    classification.
    """

    return PS2PM6Engine().evaluate(biological_features)


# A compact alias is useful for integrations that name the criterion pair as a
# single concept.  It intentionally returns the same mutually-exclusive result.
evaluate_de_novo_evidence = evaluate_ps2_pm6


__all__ = [
    "GENERAL_SVI_GUIDELINE",
    "PS2PM6Engine",
    "evaluate_de_novo_evidence",
    "evaluate_ps2_pm6",
]
