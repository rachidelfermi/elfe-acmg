"""Pinned, declarative VCEP PS1 rule support.

This module intentionally ships no mutable catalogue of gene-specific rules.
VCEP documents change independently and have distinct scope, evidence, and
strength requirements.  A deployment supplies a reviewed, versioned snapshot
of the exact rule it is authorized to apply.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Mapping, Optional

from ..evidence import Evidence, SchemaValidationError
from ..models import BiologicalFeatures, DecisionStatus, PS1Strength, ReferenceVariant, VariantConsequence
from ..output import GuidelineReference, RuleEvaluation
from .common import activated, exclusion, label, not_applicable, not_met, pairwise_evidence, undetermined
from .missense import _protein_equivalence, _splice_and_mechanism_compatibility
from .splicing import (
    _candidate_has_splice_context,
    _events_and_predictions_equivalent,
    _exonic_protein_effect_consistency,
    _functional_consequence_consistency,
    _motif_relationship,
    _reference_splice_mechanism,
    _rna_consistency,
    _transcript_compatibility,
)


_SUPPORTED_PATHWAYS = {
    "SAME_AMINO_ACID_CHANGE",
    "SPLICING_SAME_PREDICTED_EVENT",
}
_CLASSIFICATION_CONDITION_FIELDS = {
    "classification",
    "reference_classification",
    "classification_evidence",
    "reference_classification_evidence",
    "classification_review_status",
    "reference_classification_review_status",
    "classification_source",
    "reference_classification_source",
}
_VCEP_METADATA_FIELDS = (
    "vcep_gene_scope",
    "vcep_disease_scope",
    "vcep_specification_version",
    "vcep_ps1_eligibility_rules",
    "vcep_ps1_strength_rules",
    "vcep_splicing_rules",
    "vcep_reference_classification_requirements",
    "vcep_exclusion_rules",
)
_PAIRWISE_CANDIDATE_FIELDS = frozenset(
    {
        "protein_substitution_equivalence",
        "underlying_allele_equivalence",
        "protein_effect_mechanism_compatibility",
        "splice_event_equivalence",
        "splice_motif_equivalence",
        "splice_prediction_strength_comparable",
        "splice_transcript_equivalence",
    }
)
_NUMERIC_COMPARISON_ALIASES = {
    ">=": ">=", "gte": ">=", "ge": ">=", "min": ">=",
    ">": ">", "gt": ">",
    "<=": "<=", "lte": "<=", "le": "<=", "max": "<=",
    "<": "<", "lt": "<",
}


class MatchState(str, Enum):
    MATCH = "MATCH"
    NO_MATCH = "NO_MATCH"
    UNRESOLVED = "UNRESOLVED"


def _token(value: Any) -> str:
    return "_".join(
        piece
        for piece in "".join(ch if ch.isalnum() else " " for ch in str(value)).upper().split()
        if piece
    )


def _is_numeric_threshold_mapping(value: Any) -> bool:
    """Return whether a condition is an explicit numeric bound or range."""

    if not isinstance(value, Mapping) or not value or set(value) == {"in"}:
        return False
    for operator, threshold in value.items():
        if _NUMERIC_COMPARISON_ALIASES.get(str(operator).strip().lower()) is None:
            return False
        if isinstance(threshold, bool) or not isinstance(threshold, (int, float)):
            return False
    return True


def _validate_condition_value(value: Any, context: str) -> None:
    """Validate the deliberately small condition DSL at config-load time."""

    if not isinstance(value, Mapping):
        return
    if not value:
        raise SchemaValidationError(f"{context} cannot be an empty condition object")
    if set(value) == {"in"}:
        if not isinstance(value["in"], (list, tuple, set)):
            raise SchemaValidationError(f"{context}.in must be a sequence")
        return
    if not _is_numeric_threshold_mapping(value):
        raise SchemaValidationError(
            f"{context} supports only scalar equality, {{'in': [...]}}, or numeric comparison operators >=, >, <=, <"
        )


def _required_string(raw: Mapping[str, Any], key: str, context: str) -> str:
    value = raw.get(key)
    if not isinstance(value, str) or not value.strip():
        raise SchemaValidationError(f"{context}.{key} must be a non-empty string")
    return value.strip()


def _optional_string(raw: Mapping[str, Any], key: str, context: str) -> Optional[str]:
    value = raw.get(key)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise SchemaValidationError(f"{context}.{key} must be a non-empty string when supplied")
    return value.strip()


def _parse_status(raw: Any, context: str) -> DecisionStatus:
    try:
        return DecisionStatus(_token(raw))
    except ValueError as error:
        raise SchemaValidationError(f"{context}.status is not a valid decision status") from error


def _parse_strength(raw: Any, context: str) -> PS1Strength:
    try:
        return PS1Strength(_token(raw))
    except ValueError as error:
        raise SchemaValidationError(f"{context}.strength is not a valid PS1 strength") from error


@dataclass(frozen=True)
class VCEPScope:
    gene_id: str
    disease_id: str
    inheritance_context: Optional[str] = None
    candidate_transcript_id_version: Optional[str] = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any], context: str) -> "VCEPScope":
        return cls(
            gene_id=_required_string(raw, "gene_id", context),
            disease_id=_required_string(raw, "disease_id", context),
            inheritance_context=_optional_string(raw, "inheritance_context", context),
            candidate_transcript_id_version=_optional_string(raw, "candidate_transcript_id_version", context),
        )


@dataclass(frozen=True)
class VCEPDecisionRule:
    source_locator: str
    pathway: str
    when: Mapping[str, Any]
    reference_when: Mapping[str, Any]
    unless: Mapping[str, Any]
    reference_unless: Mapping[str, Any]
    status: DecisionStatus
    strength: PS1Strength
    decision_path: str

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any], context: str) -> "VCEPDecisionRule":
        if not isinstance(raw, Mapping):
            raise SchemaValidationError(f"{context} must be an object")
        source_locator = _required_string(raw, "source_locator", context)
        pathway = _required_string(raw, "pathway", context)
        when = raw.get("when", {})
        reference_when = raw.get("reference_when", {})
        unless = raw.get("unless", {})
        reference_unless = raw.get("reference_unless", {})
        for name, value in (
            ("when", when), ("reference_when", reference_when),
            ("unless", unless), ("reference_unless", reference_unless),
        ):
            if not isinstance(value, Mapping):
                raise SchemaValidationError(f"{context}.{name} must be an object")
            for field, expected in value.items():
                if not isinstance(field, str) or not field.strip():
                    raise SchemaValidationError(f"{context}.{name} has an invalid feature name")
                _validate_condition_value(expected, f"{context}.{name}.{field}")
        if not when and not reference_when:
            raise SchemaValidationError(f"{context} requires at least one explicit eligibility condition")
        decision = raw.get("decision")
        if not isinstance(decision, Mapping):
            raise SchemaValidationError(f"{context}.decision must be an object")
        status = _parse_status(decision.get("status"), f"{context}.decision")
        strength = _parse_strength(decision.get("strength", "NONE"), f"{context}.decision")
        if status == DecisionStatus.ACTIVATED and strength == PS1Strength.NONE:
            raise SchemaValidationError(f"{context} activates PS1 without an allowed strength")
        if status != DecisionStatus.ACTIVATED and strength != PS1Strength.NONE:
            raise SchemaValidationError(f"{context} gives non-NONE strength for non-activation")
        pathway_token = _token(pathway)
        if pathway_token not in _SUPPORTED_PATHWAYS:
            accepted = ", ".join(sorted(_SUPPORTED_PATHWAYS))
            raise SchemaValidationError(
                f"{context}.pathway must be one of {accepted}; put a VCEP-specific name in decision.decision_path"
            )
        if status == DecisionStatus.ACTIVATED and not (
            set(reference_when) & _CLASSIFICATION_CONDITION_FIELDS
        ):
            raise SchemaValidationError(
                f"{context} activates PS1 without an explicit reference-classification requirement"
            )
        return cls(
            source_locator=source_locator,
            pathway=pathway_token,
            when=dict(when),
            reference_when=dict(reference_when),
            unless=dict(unless),
            reference_unless=dict(reference_unless),
            status=status,
            strength=strength,
            decision_path=_required_string(decision, "decision_path", f"{context}.decision"),
        )


@dataclass(frozen=True)
class VCEPSpecification:
    authority: str
    name: str
    version: str
    source_url: str
    scope: VCEPScope
    rules: tuple[VCEPDecisionRule, ...]
    specification_id: Optional[str] = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any], index: int) -> "VCEPSpecification":
        if not isinstance(raw, Mapping):
            raise SchemaValidationError("vcep_specifications entries must be objects")
        context = f"vcep_specifications[{index}]"
        scope_raw = raw.get("scope", raw)
        if not isinstance(scope_raw, Mapping):
            raise SchemaValidationError(f"{context}.scope must be an object")
        source_url = _required_string(raw, "source_url", context)
        if not source_url.startswith(("https://", "http://")):
            raise SchemaValidationError(f"{context}.source_url must be an HTTP(S) URL")
        rules_raw = raw.get("rules")
        if not isinstance(rules_raw, (list, tuple)) or not rules_raw:
            raise SchemaValidationError(f"{context}.rules must be a non-empty list")
        authority = raw.get("authority", "ClinGen Variant Curation Expert Panel")
        if not isinstance(authority, str) or not authority.strip():
            raise SchemaValidationError(f"{context}.authority must be a non-empty string")
        return cls(
            authority=authority.strip(),
            name=_required_string(raw, "name", context),
            version=_required_string(raw, "version", context),
            source_url=source_url,
            scope=VCEPScope.from_mapping(scope_raw, f"{context}.scope"),
            rules=tuple(VCEPDecisionRule.from_mapping(item, f"{context}.rules[{position}]") for position, item in enumerate(rules_raw)),
            specification_id=_optional_string(raw, "specification_id", context),
        )


@dataclass(frozen=True)
class VCEPSelection:
    matched: Optional[tuple[VCEPSpecification, VCEPDecisionRule]] = None
    excluded: Optional[tuple[VCEPSpecification, VCEPDecisionRule, str]] = None
    missing: tuple[str, ...] = ()
    conflict: Optional[str] = None
    metadata_used: tuple[str, ...] = ()


def _reference_label(reference_index: int, field: str) -> str:
    """Return the public audit label for a reference-condition feature.

    VCEP rule files intentionally accept both concise fields (for example,
    ``classification``) and ``reference_classification``.  The result must
    still identify *which* comparison record supplied or lacked the evidence;
    otherwise provenance collection could accidentally look for a candidate
    feature with the same name.
    """

    suffix = field[len("reference_"):] if field.startswith("reference_") else field
    return f"references[{reference_index}].{suffix}"


def _matches_value(actual: Any, expected: Any) -> bool:
    if isinstance(expected, Mapping):
        _validate_condition_value(expected, "VCEP condition")
        if set(expected) == {"in"}:
            return any(_matches_value(actual, option) for option in expected["in"])

        # A VCEP may publish a gene-specific predictor threshold.  The engine
        # applies only the caller-pinned value; it never derives a generic
        # threshold from raw scores.  Several bounds can be combined to express
        # a closed or open range.
        normalized_ops: list[tuple[str, Any]] = []
        for operator, threshold in expected.items():
            token = str(operator).strip().lower()
            normalized = _NUMERIC_COMPARISON_ALIASES[token]
            normalized_ops.append((normalized, threshold))
        if isinstance(actual, bool) or not isinstance(actual, (int, float)):
            raise SchemaValidationError("a numeric VCEP threshold was used with a non-numeric feature")
        for operator, threshold in normalized_ops:
            if operator == ">=" and not actual >= threshold:
                return False
            if operator == ">" and not actual > threshold:
                return False
            if operator == "<=" and not actual <= threshold:
                return False
            if operator == "<" and not actual < threshold:
                return False
        return True
    if isinstance(expected, (list, tuple, set)):
        return any(_matches_value(actual, option) for option in expected)
    if hasattr(actual, "value"):
        actual = actual.value
    if hasattr(expected, "value"):
        expected = expected.value
    if isinstance(actual, str) and isinstance(expected, str):
        return _token(actual) == _token(expected)
    return actual == expected


def _candidate_evidence(evidence: BiologicalFeatures, field: str) -> Optional[Evidence[Any]]:
    if field in {
        "gene_id", "disease_id", "inheritance_context", "protein_substitution_equivalence",
        "underlying_allele_equivalence", "protein_effect_mechanism_compatibility",
        "splice_event_equivalence", "splice_motif_equivalence", "splice_prediction_strength_comparable",
        "splice_transcript_equivalence",
        "applicable_splicing_pvs1_strength", "splice_product_functional_consequence",
        "splice_rna_assay_reliability", "independent_reference_classifications",
    }:
        return getattr(evidence, field)
    if not field.startswith("candidate_"):
        return None
    suffix = field[len("candidate_"):]
    candidate = evidence.candidate
    direct = {
        "variant_id": candidate.variant_id,
        "genomic_allele": candidate.genomic_allele,
        "coding_hgvs": candidate.coding_hgvs,
        "protein_hgvs": candidate.protein_hgvs,
        "transcript_id_version": candidate.transcript_id_version,
        "protein_accession_version": candidate.protein_accession_version,
        "amino_acid_position": candidate.amino_acid_position,
        "reference_amino_acid": candidate.reference_amino_acid,
        "alternate_amino_acid": candidate.alternate_amino_acid,
        "variant_consequence": candidate.consequence,
    }
    if suffix in direct:
        return direct[suffix]
    splice_direct = {
        "splice_impact": candidate.splice.impact,
        "splice_site_type": candidate.splice.site_type,
        "splice_motif_position": candidate.splice.motif_position,
        "splice_motif_id": candidate.splice.motif_id,
        "splice_nucleotide_position": candidate.splice.nucleotide_position,
        "splice_prediction": candidate.splice.prediction,
        "spliceai_event_type": candidate.splice.spliceai_event_type,
        "spliceai_event_position": candidate.splice.spliceai_event_position,
        "splice_prediction_magnitude": candidate.splice.prediction_magnitude,
        "rna_splicing_evidence": candidate.splice.rna_splicing_evidence,
        "observed_splice_product": candidate.splice.observed_splice_product,
        "residual_normal_splicing": candidate.splice.residual_normal_splicing,
    }
    return splice_direct.get(suffix)


def _candidate_evidence_for_reference(
    evidence: BiologicalFeatures,
    reference: ReferenceVariant,
    reference_index: int,
    field: str,
) -> Optional[Evidence[Any]]:
    """Resolve pairwise VCEP conditions against the comparison in progress."""

    if field in _PAIRWISE_CANDIDATE_FIELDS:
        value, _ = pairwise_evidence(evidence, reference, reference_index, field)
        return value
    return _candidate_evidence(evidence, field)


def _candidate_condition_label(
    evidence: BiologicalFeatures,
    reference: ReferenceVariant,
    reference_index: int,
    field: str,
) -> str:
    """Return the audit label for the actual source of a candidate condition."""

    if field in _PAIRWISE_CANDIDATE_FIELDS:
        _, field_label = pairwise_evidence(evidence, reference, reference_index, field)
        return field_label
    return field


def _reference_evidence(reference: ReferenceVariant, field: str) -> Optional[Evidence[Any]]:
    # Accept both concise reference rule fields and reference_* spelling.
    suffix = field[len("reference_"):] if field.startswith("reference_") else field
    ann = reference.annotations
    direct = {
        "variant_id": ann.variant_id,
        "gene_id": reference.gene_id,
        "disease_id": reference.disease_id,
        "genomic_allele": ann.genomic_allele,
        "coding_hgvs": ann.coding_hgvs,
        "protein_hgvs": ann.protein_hgvs,
        "transcript_id_version": ann.transcript_id_version,
        "protein_accession_version": ann.protein_accession_version,
        "amino_acid_position": ann.amino_acid_position,
        "reference_amino_acid": ann.reference_amino_acid,
        "alternate_amino_acid": ann.alternate_amino_acid,
        "variant_consequence": ann.consequence,
        "classification": reference.classification,
        "classification_source": reference.classification_source,
        "classification_date": reference.classification_date,
        "classification_review_status": reference.review_status,
        "classification_evidence": reference.classification_evidence,
        "classification_independence": reference.classification_independence,
        "disease_association": reference.disease_association,
        "molecular_mechanism": reference.molecular_mechanism,
        "classification_conflicts": reference.classification_conflicts,
        "protein_substitution_equivalence": reference.protein_substitution_equivalence,
        "underlying_allele_equivalence": reference.underlying_allele_equivalence,
        "protein_effect_mechanism_compatibility": reference.protein_effect_mechanism_compatibility,
        "splice_event_equivalence": reference.splice_event_equivalence,
        "splice_motif_equivalence": reference.splice_motif_equivalence,
        "splice_prediction_strength_comparable": reference.splice_prediction_strength_comparable,
        "splice_transcript_equivalence": reference.splice_transcript_equivalence,
        "splice_impact": ann.splice.impact,
        "splice_site_type": ann.splice.site_type,
        "splice_motif_position": ann.splice.motif_position,
        "splice_motif_id": ann.splice.motif_id,
        "splice_nucleotide_position": ann.splice.nucleotide_position,
        "splice_prediction": ann.splice.prediction,
        "spliceai_event_type": ann.splice.spliceai_event_type,
        "spliceai_event_position": ann.splice.spliceai_event_position,
        "splice_prediction_magnitude": ann.splice.prediction_magnitude,
        "rna_splicing_evidence": ann.splice.rna_splicing_evidence,
        "observed_splice_product": ann.splice.observed_splice_product,
        "residual_normal_splicing": ann.splice.residual_normal_splicing,
    }
    return direct.get(suffix)


def _condition_match(
    conditions: Mapping[str, Any], resolver: Any
) -> tuple[MatchState, tuple[str, ...]]:
    missing: list[str] = []
    for field, expected in conditions.items():
        evidence = resolver(field)
        if evidence is None:
            raise SchemaValidationError(f"VCEP condition refers to unsupported feature {field!r}")
        if evidence.is_conflicting or evidence.is_unknown:
            missing.append(field)
            continue
        if not _matches_value(evidence.value, expected):
            return MatchState.NO_MATCH, ()
    return (MatchState.UNRESOLVED, tuple(dict.fromkeys(missing))) if missing else (MatchState.MATCH, ())


def _scope_match(scope: VCEPScope, evidence: BiologicalFeatures) -> tuple[MatchState, tuple[str, ...]]:
    conditions: dict[str, Any] = {"gene_id": scope.gene_id, "disease_id": scope.disease_id}
    if scope.inheritance_context is not None:
        conditions["inheritance_context"] = scope.inheritance_context
    if scope.candidate_transcript_id_version is not None:
        conditions["candidate_transcript_id_version"] = scope.candidate_transcript_id_version
    return _condition_match(conditions, lambda field: _candidate_evidence(evidence, field))


def _vcep_metadata_state(evidence: BiologicalFeatures) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Return supplied and conflicting legacy VCEP metadata declarations.

    The individual inventory fields are retained for typed intake, but are not
    independently executable rules.  A full immutable `vcep_specifications`
    snapshot is required to activate PS1.  Recognising these declarations here
    prevents them from being silently ignored and from permitting a generic
    fallback when an applicable VCEP rule may be incomplete.
    """

    supplied: list[str] = []
    conflicting: list[str] = []
    for field in _VCEP_METADATA_FIELDS:
        value = getattr(evidence, field)
        if value.is_conflicting:
            conflicting.append(field)
        elif value.is_known:
            supplied.append(field)
    return tuple(supplied), tuple(conflicting)


def _metadata_matches_specification(
    evidence: BiologicalFeatures, specification: VCEPSpecification
) -> tuple[bool, tuple[str, ...]]:
    """Validate redundant scope/version inventory fields against a snapshot."""

    checks = (
        ("vcep_gene_scope", specification.scope.gene_id),
        ("vcep_disease_scope", specification.scope.disease_id),
        ("vcep_specification_version", specification.version),
    )
    mismatches: list[str] = []
    for field, expected in checks:
        declared = getattr(evidence, field)
        if declared.is_known and not _matches_value(declared.value, expected):
            mismatches.append(field)
    return not mismatches, tuple(mismatches)


def _legacy_metadata_scope_state(
    evidence: BiologicalFeatures, metadata_used: tuple[str, ...]
) -> tuple[MatchState, tuple[str, ...]]:
    """Determine whether legacy VCEP metadata could describe this candidate.

    Inventory metadata is not executable, but it must not be ignored when it
    declares an apparently applicable VCEP and the pinned configuration does
    not contain that specification.
    """

    if not metadata_used:
        return MatchState.NO_MATCH, ()
    checks = (
        ("vcep_gene_scope", evidence.vcep_gene_scope, "gene_id", evidence.gene_id),
        ("vcep_disease_scope", evidence.vcep_disease_scope, "disease_id", evidence.disease_id),
    )
    missing: list[str] = []
    for declared_name, declared, candidate_name, candidate in checks:
        if declared.is_unknown:
            missing.append(declared_name)
            continue
        if candidate.is_unknown or candidate.is_conflicting:
            missing.append(candidate_name)
            continue
        if not _matches_value(declared.value, candidate.value):
            return MatchState.NO_MATCH, ()
    if missing:
        return MatchState.UNRESOLVED, tuple(dict.fromkeys(missing))
    return MatchState.MATCH, ()


def select_vcep_rule(
    evidence: BiologicalFeatures, reference: ReferenceVariant, reference_index: int = 0
) -> VCEPSelection:
    """Select one exact scoped VCEP rule for this reference, if any.

    A supplied but potentially applicable specification with missing scope or
    rule conditions blocks fallback to a generic pathway.  That prevents a
    generic rule from silently overriding a higher-precedence VCEP document.
    """

    metadata_used, metadata_conflicts = _vcep_metadata_state(evidence)
    if metadata_conflicts:
        return VCEPSelection(
            conflict="VCEP metadata declarations are conflicting: " + ", ".join(metadata_conflicts),
            metadata_used=metadata_used + metadata_conflicts,
        )
    supplied = evidence.vcep_specifications
    if supplied.is_unknown:
        if metadata_used:
            return VCEPSelection(
                missing=("vcep_specifications",),
                metadata_used=metadata_used,
            )
        return VCEPSelection()
    if supplied.is_conflicting:
        return VCEPSelection(
            conflict="VCEP specification configuration is conflicting",
            metadata_used=metadata_used,
        )
    specifications = tuple(
        VCEPSpecification.from_mapping(raw, index)
        for index, raw in enumerate(supplied.value or ())
    )
    matched_specs: list[VCEPSpecification] = []
    unresolved_scope: list[str] = []
    for specification in specifications:
        state, missing = _scope_match(specification.scope, evidence)
        if state == MatchState.MATCH:
            matched_specs.append(specification)
        elif state == MatchState.UNRESOLVED:
            unresolved_scope.extend(missing)
    if unresolved_scope:
        return VCEPSelection(
            missing=tuple(dict.fromkeys(unresolved_scope)), metadata_used=metadata_used
        )

    if not matched_specs:
        legacy_state, legacy_missing = _legacy_metadata_scope_state(evidence, metadata_used)
        if legacy_state != MatchState.NO_MATCH:
            return VCEPSelection(
                missing=tuple(dict.fromkeys(("vcep_specifications",) + legacy_missing)),
                metadata_used=metadata_used,
            )

    if matched_specs:
        compatible_specs: list[VCEPSpecification] = []
        metadata_mismatches: list[str] = []
        for specification in matched_specs:
            compatible, mismatches = _metadata_matches_specification(evidence, specification)
            if compatible:
                compatible_specs.append(specification)
            else:
                metadata_mismatches.extend(mismatches)
        if not compatible_specs:
            return VCEPSelection(
                conflict=(
                    "Supplied VCEP scope/version metadata conflicts with the pinned specification: "
                    + ", ".join(dict.fromkeys(metadata_mismatches))
                ),
                metadata_used=metadata_used,
            )
        matched_specs = compatible_specs

    applicable_versions = {specification.version for specification in matched_specs}
    if len(applicable_versions) > 1:
        return VCEPSelection(
            conflict=(
                "Multiple VCEP specification versions match the candidate scope; "
                "supply vcep_specification_version to pin one version."
            ),
            metadata_used=metadata_used,
        )

    matching: list[tuple[VCEPSpecification, VCEPDecisionRule]] = []
    exclusions: list[tuple[VCEPSpecification, VCEPDecisionRule, str]] = []
    unresolved: list[str] = []
    for specification in matched_specs:
        for rule in specification.rules:
            candidate_state, candidate_missing = _condition_match(
                rule.when,
                lambda field: _candidate_evidence_for_reference(
                    evidence, reference, reference_index, field
                ),
            )
            candidate_missing = tuple(
                _candidate_condition_label(evidence, reference, reference_index, field)
                for field in candidate_missing
            )
            if candidate_state == MatchState.NO_MATCH:
                continue
            reference_state, reference_missing = _condition_match(
                rule.reference_when, lambda field: _reference_evidence(reference, field)
            )
            reference_missing = tuple(
                _reference_label(reference_index, field) for field in reference_missing
            )
            if reference_state == MatchState.NO_MATCH:
                continue
            if rule.unless:
                unless_state, unless_missing = _condition_match(
                    rule.unless,
                    lambda field: _candidate_evidence_for_reference(
                        evidence, reference, reference_index, field
                    ),
                )
                unless_missing = tuple(
                    _candidate_condition_label(evidence, reference, reference_index, field)
                    for field in unless_missing
                )
                if unless_state == MatchState.MATCH:
                    exclusions.append((specification, rule, "unless"))
                    continue
            else:
                unless_state, unless_missing = MatchState.NO_MATCH, ()
            if rule.reference_unless:
                reference_unless_state, reference_unless_missing = _condition_match(
                    rule.reference_unless, lambda field: _reference_evidence(reference, field)
                )
                reference_unless_missing = tuple(
                    _reference_label(reference_index, field)
                    for field in reference_unless_missing
                )
                if reference_unless_state == MatchState.MATCH:
                    exclusions.append((specification, rule, "reference_unless"))
                    continue
            else:
                reference_unless_state, reference_unless_missing = MatchState.NO_MATCH, ()
            missing = candidate_missing + reference_missing + unless_missing + reference_unless_missing
            if candidate_state == MatchState.UNRESOLVED or reference_state == MatchState.UNRESOLVED or unless_state == MatchState.UNRESOLVED or reference_unless_state == MatchState.UNRESOLVED:
                unresolved.extend(missing)
                continue
            matching.append((specification, rule))

    if len(matching) > 1:
        return VCEPSelection(
            conflict="multiple equally applicable VCEP PS1 rules match the supplied evidence",
            metadata_used=metadata_used,
        )
    if unresolved:
        # There is no implicit rule priority in the compact VCEP schema.  An
        # unresolved scoped rule could be an exclusion or a different strength
        # path, so it must block a favorable rule rather than be ignored.
        return VCEPSelection(
            missing=tuple(dict.fromkeys(unresolved)), metadata_used=metadata_used
        )
    if exclusions:
        # `unless` conditions are documented VCEP exclusion gates.  Letting a
        # generic rule activate after one matches would silently override the
        # scoped specification.
        return VCEPSelection(excluded=exclusions[0], metadata_used=metadata_used)
    if len(matching) == 1:
        return VCEPSelection(matched=matching[0], metadata_used=metadata_used)
    return VCEPSelection(metadata_used=metadata_used)


def _has_pinned_splice_prediction_threshold(rule: VCEPDecisionRule) -> bool:
    """Whether a matched VCEP rule explicitly supplies a predictor threshold.

    This is intentionally narrower than a generic numeric condition: only a
    threshold on the *candidate* splice-prediction magnitude may replace
    Walker's generic pairwise strength assertion. A reference-only bound says
    nothing about whether the candidate prediction is sufficient.
    """

    for field, expected in rule.when.items():
        if field == "candidate_splice_prediction_magnitude" and _is_numeric_threshold_mapping(expected):
            return True
    return False


def _vcep_standard_safety(
    evidence: BiologicalFeatures,
    reference: ReferenceVariant,
    reference_index: int,
    guideline: GuidelineReference,
    decision_path: str,
) -> Optional[RuleEvaluation]:
    """Apply generic non-overridable safety checks to a VCEP AA rule.

    A VCEP can change a reference-classification requirement or code strength,
    but it cannot make different protein effects, a known splice mechanism, or
    an unresolved transcript/protein mapping into a valid same-AA PS1 match.
    """

    consequence = evidence.candidate.consequence
    if consequence.is_known and consequence.value in {
        VariantConsequence.SYNONYMOUS,
        VariantConsequence.INTRONIC,
    }:
        return not_applicable(
            decision_path, guideline,
            "The supplied candidate consequence is not a protein-altering same-amino-acid comparison.",
            reference_index=reference_index, used=("candidate_variant_consequence",),
        )
    protein = _protein_equivalence(evidence, reference, reference_index)
    if protein.state == "NOT_MET":
        return not_met(
            decision_path, guideline, protein.reason, reference_index=reference_index,
            used=protein.used, exclusions=protein.exclusions,
        )
    if protein.state == "UNDETERMINED":
        return undetermined(
            decision_path, guideline, protein.reason, protein.missing,
            reference_index=reference_index, used=protein.used,
            exclusions=protein.exclusions,
        )
    mechanism = _splice_and_mechanism_compatibility(evidence, reference, reference_index)
    used = protein.used + mechanism.used
    if mechanism.state == "NOT_MET":
        return not_met(
            decision_path, guideline, mechanism.reason, reference_index=reference_index,
            used=used, exclusions=protein.exclusions + mechanism.exclusions,
        )
    if mechanism.state == "UNDETERMINED":
        return undetermined(
            decision_path, guideline, mechanism.reason, mechanism.missing,
            reference_index=reference_index, used=used,
            exclusions=protein.exclusions + mechanism.exclusions,
        )
    return None


def _vcep_splice_safety(
    evidence: BiologicalFeatures,
    reference: ReferenceVariant,
    reference_index: int,
    guideline: GuidelineReference,
    decision_path: str,
    *,
    has_pinned_prediction_threshold: bool,
) -> Optional[RuleEvaluation]:
    """Apply generic non-overridable event-comparison checks to a VCEP rule.

    The VCEP configuration supplies its own verified strength and any
    gene-specific threshold.  It cannot waive exact event/motif/transcript
    comparability or known RNA/protein-effect contradictions.
    """

    if not _candidate_has_splice_context(evidence):
        return not_applicable(
            decision_path, guideline,
            "No candidate splice context was supplied for the VCEP splice comparison.",
            reference_index=reference_index,
        )

    mechanism = _reference_splice_mechanism(reference, reference_index)
    if mechanism.state == "NOT_MET":
        return not_met(decision_path, guideline, mechanism.reason, reference_index=reference_index,
            used=mechanism.used, exclusions=mechanism.exclusions)
    if mechanism.state == "UNDETERMINED":
        return undetermined(decision_path, guideline, mechanism.reason, mechanism.missing,
            reference_index=reference_index, used=mechanism.used,
            exclusions=mechanism.exclusions)

    transcript = _transcript_compatibility(evidence, reference, reference_index)
    if transcript.state == "NOT_MET":
        return not_met(decision_path, guideline, transcript.reason, reference_index=reference_index,
            used=transcript.used, exclusions=transcript.exclusions)
    if transcript.state == "UNDETERMINED":
        return undetermined(decision_path, guideline, transcript.reason, transcript.missing,
            reference_index=reference_index, used=transcript.used, exclusions=transcript.exclusions)

    events = _events_and_predictions_equivalent(
        evidence,
        reference,
        reference_index,
        require_prediction_strength_comparable=not has_pinned_prediction_threshold,
    )
    if events.state == "NOT_MET":
        return not_met(decision_path, guideline, events.reason, reference_index=reference_index,
            used=transcript.used + events.used, exclusions=events.exclusions)
    if events.state == "UNDETERMINED":
        return undetermined(decision_path, guideline, events.reason, events.missing,
            reference_index=reference_index, used=transcript.used + events.used,
            exclusions=events.exclusions)

    motif = _motif_relationship(evidence, reference, reference_index)
    if motif.state == "NOT_MET":
        return not_met(decision_path, guideline, motif.reason, reference_index=reference_index,
            used=transcript.used + events.used + motif.used, exclusions=motif.exclusions)
    if motif.state == "UNDETERMINED":
        return undetermined(decision_path, guideline, motif.reason, motif.missing,
            reference_index=reference_index, used=transcript.used + events.used + motif.used,
            exclusions=motif.exclusions)
    if motif.state == "NOT_APPLICABLE":
        return not_applicable(decision_path, guideline, motif.reason, reference_index=reference_index,
            used=transcript.used + events.used + motif.used)

    protein = _exonic_protein_effect_consistency(evidence, reference, reference_index)
    if protein.state == "NOT_MET":
        return not_met(decision_path, guideline, protein.reason, reference_index=reference_index,
            used=transcript.used + events.used + motif.used + protein.used,
            exclusions=protein.exclusions)
    if protein.state == "UNDETERMINED":
        return undetermined(decision_path, guideline, protein.reason, protein.missing,
            reference_index=reference_index, used=transcript.used + events.used + motif.used + protein.used,
            exclusions=protein.exclusions)

    rna = _rna_consistency(evidence, reference, reference_index)
    if rna.state == "NOT_MET":
        return not_met(decision_path, guideline, rna.reason, reference_index=reference_index,
            used=transcript.used + events.used + motif.used + protein.used + rna.used,
            exclusions=rna.exclusions)
    if rna.state == "UNDETERMINED":
        return undetermined(decision_path, guideline, rna.reason, rna.missing,
            reference_index=reference_index, used=transcript.used + events.used + motif.used + protein.used + rna.used,
            exclusions=rna.exclusions)

    functional = _functional_consequence_consistency(evidence)
    if functional.state == "UNDETERMINED":
        return undetermined(decision_path, guideline, functional.reason, functional.missing,
            reference_index=reference_index,
            used=transcript.used + events.used + motif.used + protein.used + rna.used + functional.used)
    return None


def evaluate_selected_vcep_rule(
    selection: VCEPSelection,
    evidence: BiologicalFeatures,
    reference: ReferenceVariant,
    reference_index: int,
) -> Optional[RuleEvaluation]:
    """Convert a selected rule (or safe selection failure) to an outcome."""

    generic = GuidelineReference(
        authority="ClinGen Variant Curation Expert Panel",
        name="VCEP PS1 specification selection",
        version="caller-supplied pinned configuration",
        source_url="https://cspec.genome.network/cspec/ui/svi/summary",
        source_locator="VCEP scope and rule selection",
        rule_id="PS1",
    )
    if selection.conflict:
        return undetermined(
            "VCEP_PS1_SELECTION", generic, selection.conflict, (), reference_index=reference_index,
            used=("vcep_specifications",) + selection.metadata_used,
        )
    if selection.missing:
        return undetermined(
            "VCEP_PS1_SELECTION", generic,
            "A supplied VCEP specification may apply, but its scope or exact rule conditions are unresolved.",
            selection.missing, reference_index=reference_index,
            used=("vcep_specifications",) + selection.metadata_used,
        )
    if selection.excluded is not None:
        specification, rule, exclusion_kind = selection.excluded
        guideline = GuidelineReference(
            authority=specification.authority,
            name=specification.name,
            version=specification.version,
            source_url=specification.source_url,
            source_locator=rule.source_locator,
            rule_id=specification.specification_id or "PS1",
        )
        exclusion_fields = rule.unless if exclusion_kind == "unless" else rule.reference_unless
        used = (
            ("vcep_specifications",)
            + selection.metadata_used
            + tuple(
                _candidate_condition_label(evidence, reference, reference_index, field)
                for field in rule.when
            )
            + tuple(_reference_label(reference_index, field) for field in rule.reference_when)
            + (
                tuple(
                    _candidate_condition_label(evidence, reference, reference_index, field)
                    for field in exclusion_fields
                )
                if exclusion_kind == "unless"
                else tuple(_reference_label(reference_index, field) for field in exclusion_fields)
            )
        )
        return not_met(
            rule.decision_path,
            guideline,
            "A pinned, disease-scoped VCEP rule explicitly excludes PS1 for the supplied evidence.",
            reference_index=reference_index,
            used=used,
        )
    if selection.matched is None:
        return None
    specification, rule = selection.matched
    guideline = GuidelineReference(
        authority=specification.authority,
        name=specification.name,
        version=specification.version,
        source_url=specification.source_url,
        source_locator=rule.source_locator,
        rule_id=specification.specification_id or "PS1",
    )
    used = (
        ("vcep_specifications",)
        + selection.metadata_used
        + tuple(
            _candidate_condition_label(evidence, reference, reference_index, field)
            for field in rule.when
        )
        + tuple(_reference_label(reference_index, field) for field in rule.reference_when)
        + tuple(
            _candidate_condition_label(evidence, reference, reference_index, field)
            for field in rule.unless
        )
        + tuple(_reference_label(reference_index, field) for field in rule.reference_unless)
    )
    if rule.status == DecisionStatus.ACTIVATED:
        safety = (
            _vcep_standard_safety(
                evidence, reference, reference_index, guideline, rule.decision_path
            )
            if rule.pathway == "SAME_AMINO_ACID_CHANGE"
            else _vcep_splice_safety(
                evidence,
                reference,
                reference_index,
                guideline,
                rule.decision_path,
                has_pinned_prediction_threshold=_has_pinned_splice_prediction_threshold(rule),
            )
        )
        if safety is not None:
            return RuleEvaluation(
                status=safety.status,
                strength=safety.strength,
                decision_path=safety.decision_path,
                guideline=safety.guideline,
                reason=safety.reason,
                reference_index=safety.reference_index,
                features_used=tuple(dict.fromkeys(used + safety.features_used)),
                missing_critical_features=safety.missing_critical_features,
                exclusions_evaluated=safety.exclusions_evaluated,
            )
        return activated(
            rule.strength, rule.decision_path, guideline,
            "A pinned, disease-scoped VCEP PS1 rule explicitly matches the supplied evidence.",
            reference_index=reference_index,
            used=used + ((label(reference_index, "molecular_mechanism"),)
                         if rule.pathway == "SPLICING_SAME_PREDICTED_EVENT" else ()),
        )
    if rule.status == DecisionStatus.NOT_MET:
        return not_met(
            rule.decision_path, guideline,
            "A pinned, disease-scoped VCEP rule explicitly excludes PS1 for the supplied evidence.",
            reference_index=reference_index, used=used,
        )
    if rule.status == DecisionStatus.NOT_APPLICABLE:
        return not_applicable(
            rule.decision_path, guideline,
            "A pinned, disease-scoped VCEP rule declares PS1 not applicable for the supplied evidence.",
            reference_index=reference_index, used=used,
        )
    return undetermined(
        rule.decision_path, guideline,
        "A pinned, disease-scoped VCEP rule declares PS1 unresolved for the supplied evidence.",
        (), reference_index=reference_index, used=used,
    )

