"""Versioned, scoped VCEP PVS1 specification support.

This module deliberately does not ship a partial catalogue of live VCEP rules.
Those documents are versioned, disease-scoped, and change independently of the
engine.  Instead, a caller supplies a pinned declarative snapshot.  The engine
validates its provenance and refuses to choose between conflicting matching
rules.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import re
from typing import Any, Mapping, Optional, Sequence

from .models import (
    BiologicalFeatures,
    DecisionStatus,
    EvidenceStatus,
    GuidelineReference,
    PVS1Strength,
    SchemaValidationError,
)


class MatchState(str, Enum):
    MATCH = "MATCH"
    NO_MATCH = "NO_MATCH"
    UNRESOLVED = "UNRESOLVED"


# These are engine pathway selectors, not ClinGen rule identifiers.  Keeping
# this finite prevents a misspelled VCEP pathway from being silently ignored.
_SUPPORTED_VCEP_PATHWAYS = frozenset(
    {
        "ANY",
        "NONSENSE_OR_FRAMESHIFT",
        "CANONICAL_SPLICE",
        "NONCANONICAL_SPLICE",
        "EXON_DELETION",
        "WHOLE_GENE_DELETION",
        "DUPLICATION",
        "INITIATION_CODON",
    }
)

_COLLECTION_CONDITION_FIELDS = frozenset({"affected_exons", "splice_outcomes"})
_CONDITION_OPERATORS = frozenset(
    {
        "in",
        "lt",
        "lte",
        "gt",
        "gte",
        "between",
        "contains_any",
        "contains_all",
        "contains_exact",
    }
)


def _token(value: Any) -> str:
    if isinstance(value, Enum):
        value = value.value
    return re.sub(r"[^A-Z0-9]+", "_", str(value).upper()).strip("_")


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
    if isinstance(raw, DecisionStatus):
        return raw
    try:
        return DecisionStatus(_token(raw))
    except ValueError as error:
        raise SchemaValidationError(f"{context}.status is not a valid decision status") from error


def _parse_strength(raw: Any, context: str) -> PVS1Strength:
    if isinstance(raw, PVS1Strength):
        return raw
    try:
        return PVS1Strength(_token(raw))
    except ValueError as error:
        raise SchemaValidationError(f"{context}.strength is not a valid PVS1 strength") from error


def _known_condition_value(field: str, raw: Any, context: str) -> Any:
    """Validate one declarative condition atom using the public feature schema."""

    try:
        parsed = BiologicalFeatures.from_mapping({field: raw}).evidence(field)
    except SchemaValidationError as error:
        raise SchemaValidationError(
            f"{context} has an invalid value for {field!r}: {error}"
        ) from error
    if parsed.status != EvidenceStatus.KNOWN:
        raise SchemaValidationError(
            f"{context} must use a known value, not an evidence state, for {field!r}"
        )
    return parsed.value


def _known_collection_member(field: str, raw: Any, context: str) -> Any:
    """Validate a member used by an explicit collection-condition operator."""

    if isinstance(raw, (Mapping, list, tuple, set)):
        raise SchemaValidationError(
            f"{context} collection operators require scalar members"
        )
    if field == "affected_exons":
        value = _known_condition_value(field, [raw], context)
        return value[0]
    if field == "splice_outcomes":
        value = _known_condition_value(field, raw, context)
        # A scalar public splice alias (for example ``frameshift_nmd``)
        # normalizes to a one-member tuple in the feature schema.
        return value[0]
    return _known_condition_value(field, raw, context)


def _validate_numeric(value: Any, field: str, context: str) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SchemaValidationError(
            f"{context} requires a numeric feature for its {field!r} condition"
        )


def _validate_condition(field: str, expected: Any, context: str) -> None:
    """Validate the small, non-executable VCEP condition DSL at load time.

    A condition is an exact value (or backwards-compatible list of exact
    alternatives), or one explicit operator.  This is intentionally narrower
    than Python expressions: a pinned clinical specification must not execute
    caller-provided code.
    """

    if field == "gene_specific_pvs1_rules":
        raise SchemaValidationError(
            f"{context} cannot recursively condition on gene_specific_pvs1_rules"
        )

    if isinstance(expected, Mapping):
        if len(expected) != 1:
            raise SchemaValidationError(
                f"{context} must contain exactly one supported condition operator"
            )
        operator, operand = next(iter(expected.items()))
        if operator not in _CONDITION_OPERATORS:
            choices = ", ".join(sorted(_CONDITION_OPERATORS))
            raise SchemaValidationError(
                f"{context} uses unsupported condition operator {operator!r}; expected one of {choices}"
            )
        if operator == "in":
            if not isinstance(operand, (list, tuple, set)) or not operand:
                raise SchemaValidationError(f"{context}.in must be a non-empty list")
            if field in _COLLECTION_CONDITION_FIELDS:
                raise SchemaValidationError(
                    f"{context}.in is ambiguous for collection feature {field!r}; "
                    "use contains_any, contains_all, or contains_exact"
                )
            for item in operand:
                _known_condition_value(field, item, context)
            return
        if operator in {"lt", "lte", "gt", "gte"}:
            value = _known_condition_value(field, operand, context)
            _validate_numeric(value, field, context)
            return
        if operator == "between":
            if not isinstance(operand, (list, tuple)) or len(operand) != 2:
                raise SchemaValidationError(f"{context}.between must be a two-item inclusive range")
            lower = _known_condition_value(field, operand[0], context)
            upper = _known_condition_value(field, operand[1], context)
            _validate_numeric(lower, field, context)
            _validate_numeric(upper, field, context)
            if lower > upper:
                raise SchemaValidationError(f"{context}.between lower bound exceeds upper bound")
            return
        if field not in _COLLECTION_CONDITION_FIELDS:
            raise SchemaValidationError(
                f"{context}.{operator} is only valid for a collection-valued feature"
            )
        if not isinstance(operand, (list, tuple, set)) or not operand:
            raise SchemaValidationError(f"{context}.{operator} must be a non-empty list")
        members = tuple(
            _known_collection_member(field, item, context) for item in operand
        )
        for index, member in enumerate(members):
            if any(
                _scalar_matches(member, earlier) for earlier in members[:index]
            ):
                raise SchemaValidationError(
                    f"{context}.{operator} must not contain duplicate members"
                )
        return

    if isinstance(expected, (list, tuple, set)):
        if not expected:
            raise SchemaValidationError(f"{context} must not use an empty condition list")
        if field in _COLLECTION_CONDITION_FIELDS:
            raise SchemaValidationError(
                f"{context} is ambiguous for collection feature {field!r}; "
                "use contains_any, contains_all, or contains_exact"
            )
        for item in expected:
            _known_condition_value(field, item, context)
        return
    _known_condition_value(field, expected, context)


def _normalize_condition(field: str, expected: Any, context: str) -> Any:
    """Validate and store conditions in the same canonical value space as input.

    Validation alone is insufficient: public feature aliases such as
    ``frameshift_variant`` or ``nonsense`` must also match the normalized enum
    stored on ``BiologicalFeatures`` at evaluation time. The returned DSL is
    still declarative and non-executable; only its operands are canonicalized.
    """

    _validate_condition(field, expected, context)

    def normalized_items(raw_items: Sequence[Any] | set[Any], parser: Any) -> tuple[Any, ...]:
        """Preserve declared sequence order, but canonicalize unordered sets.

        Set operands are accepted for convenient programmatic manifest
        construction. They cannot leak their hash-randomized iteration order
        into an auditable decision payload, however.
        """

        parsed = tuple(parser(item) for item in raw_items)
        if isinstance(raw_items, set):
            return tuple(
                sorted(
                    parsed,
                    key=lambda value: (
                        type(value).__name__,
                        _token(value) if isinstance(value, (str, Enum)) else repr(value),
                    ),
                )
            )
        return parsed

    if isinstance(expected, Mapping):
        operator, operand = next(iter(expected.items()))
        if operator == "in":
            return {
                operator: normalized_items(
                    operand,
                    lambda item: _known_condition_value(field, item, context),
                )
            }
        if operator in {"lt", "lte", "gt", "gte"}:
            return {operator: _known_condition_value(field, operand, context)}
        if operator == "between":
            return {
                operator: tuple(
                    _known_condition_value(field, item, context) for item in operand
                )
            }
        return {
            operator: normalized_items(
                operand,
                lambda item: _known_collection_member(field, item, context),
            )
        }
    if isinstance(expected, (list, tuple, set)):
        return normalized_items(
            expected,
            lambda item: _known_condition_value(field, item, context),
        )
    return _known_condition_value(field, expected, context)


@dataclass(frozen=True)
class VCEPScope:
    gene_id: str
    disease_id: str
    transcript_id: Optional[str] = None
    transcript_version: Optional[str] = None
    mode_of_inheritance: Optional[str] = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any], context: str) -> "VCEPScope":
        # The compact form keeps gene/disease at the specification top level.
        return cls(
            gene_id=_required_string(raw, "gene_id", context),
            disease_id=_required_string(raw, "disease_id", context),
            transcript_id=_optional_string(raw, "transcript_id", context),
            transcript_version=_optional_string(raw, "transcript_version", context),
            mode_of_inheritance=_optional_string(raw, "mode_of_inheritance", context),
        )


@dataclass(frozen=True)
class VCEPDecisionRule:
    source_locator: str
    pathway: str
    when: Mapping[str, Any]
    status: DecisionStatus
    strength: PVS1Strength
    decision_path: str
    evidence_mode: Optional[str] = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any], context: str) -> "VCEPDecisionRule":
        if not isinstance(raw, Mapping):
            raise SchemaValidationError(f"{context} must be an object")
        source_locator = _required_string(raw, "source_locator", context)
        pathway = _required_string(raw, "pathway", context)
        when_raw = raw.get("when")
        if not isinstance(when_raw, Mapping) or not when_raw:
            raise SchemaValidationError(
                f"{context}.when must be a non-empty mapping of explicit conditions"
            )
        aliases = {"splice_effect": "splice_outcomes"}
        when: dict[str, Any] = {}
        for field, expected in when_raw.items():
            normalized = aliases.get(field, field)
            if normalized in when:
                raise SchemaValidationError(
                    f"{context}.when supplies both {field!r} and {normalized!r}"
                )
            when[normalized] = expected
        unknown_fields = set(when) - set(BiologicalFeatures.__dataclass_fields__)
        if unknown_fields:
            raise SchemaValidationError(
                f"{context}.when refers to unsupported feature(s): "
                + ", ".join(sorted(unknown_fields))
            )
        for field, expected in tuple(when.items()):
            when[field] = _normalize_condition(
                field, expected, f"{context}.when.{field}"
            )
        pathway_token = _token(pathway)
        if pathway_token not in _SUPPORTED_VCEP_PATHWAYS:
            choices = ", ".join(sorted(_SUPPORTED_VCEP_PATHWAYS))
            raise SchemaValidationError(
                f"{context}.pathway has unsupported value {pathway!r}; expected one of {choices}"
            )
        decision = raw.get("decision")
        if not isinstance(decision, Mapping):
            raise SchemaValidationError(f"{context}.decision must be an object")
        status = _parse_status(decision.get("status"), f"{context}.decision")
        strength = _parse_strength(decision.get("strength", "NONE"), f"{context}.decision")
        if status == DecisionStatus.ACTIVATED and strength == PVS1Strength.NONE:
            raise SchemaValidationError(f"{context} activates PVS1 without a strength")
        if status != DecisionStatus.ACTIVATED and strength != PVS1Strength.NONE:
            raise SchemaValidationError(f"{context} has non-NONE strength for non-activation")
        decision_path = _required_string(decision, "decision_path", f"{context}.decision")
        evidence_mode = _optional_string(decision, "evidence_mode", f"{context}.decision")
        return cls(
            source_locator=source_locator,
            pathway=pathway_token,
            when=dict(when),
            status=status,
            strength=strength,
            decision_path=decision_path,
            evidence_mode=evidence_mode,
        )


@dataclass(frozen=True)
class VCEPSpecification:
    """A caller-pinned VCEP snapshot with enough identity to be auditable."""

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
            raise SchemaValidationError("gene_specific_pvs1_rules entries must be objects")
        context = f"gene_specific_pvs1_rules[{index}]"
        name = _required_string(raw, "name", context)
        version = _required_string(raw, "version", context)
        source_url = _required_string(raw, "source_url", context)
        if not source_url.startswith(("https://", "http://")):
            raise SchemaValidationError(f"{context}.source_url must be an HTTP(S) URL")
        scope_raw = raw.get("scope", raw)
        if not isinstance(scope_raw, Mapping):
            raise SchemaValidationError(f"{context}.scope must be an object")
        scope = VCEPScope.from_mapping(scope_raw, context)
        rules_raw = raw.get("rules")
        if not isinstance(rules_raw, (list, tuple)) or not rules_raw:
            raise SchemaValidationError(f"{context}.rules must be a non-empty list")
        rules = tuple(
            VCEPDecisionRule.from_mapping(rule, f"{context}.rules[{rule_index}]")
            for rule_index, rule in enumerate(rules_raw)
        )
        authority = raw.get("authority", "ClinGen VCEP")
        if not isinstance(authority, str) or not authority.strip():
            raise SchemaValidationError(f"{context}.authority must be a non-empty string")
        specification_id = _optional_string(raw, "specification_id", context)
        return cls(
            authority=authority.strip(),
            name=name,
            version=version,
            source_url=source_url,
            scope=scope,
            rules=rules,
            specification_id=specification_id,
        )


@dataclass(frozen=True)
class SpecificRuleSelection:
    """Result of VCEP rule selection, before a generic rule is considered."""

    matched: Optional[tuple[VCEPSpecification, VCEPDecisionRule]] = None
    missing: tuple[str, ...] = ()
    conflict: Optional[str] = None


def validate_specific_rule_declarations(evidence: BiologicalFeatures) -> None:
    """Eagerly validate every supplied VCEP snapshot.

    Schema validity must not depend on which generic branch happens to return
    first. Selection reparses the same immutable caller payload later so it
    can remain a small standalone function.
    """

    supplied = evidence.gene_specific_pvs1_rules
    if supplied.status != EvidenceStatus.KNOWN:
        return
    for index, raw in enumerate(supplied.value or ()):
        VCEPSpecification.from_mapping(raw, index)


def _scalar_matches(actual: Any, expected: Any) -> bool:
    actual_was_enum = isinstance(actual, Enum)
    expected_was_enum = isinstance(expected, Enum)
    if isinstance(actual, Enum):
        actual = actual.value
    if isinstance(expected, Enum):
        expected = expected.value
    if (actual_was_enum or expected_was_enum) and isinstance(actual, str) and isinstance(expected, str):
        return _token(actual) == _token(expected)
    return actual == expected


def _collection_members_match(actual: Any, expected: Sequence[Any], mode: str) -> bool:
    if not isinstance(actual, (tuple, list, set)):
        return False
    actual_items = tuple(actual)
    expected_items = tuple(expected)
    if mode == "contains_any":
        return any(
            _scalar_matches(actual_item, expected_item)
            for actual_item in actual_items
            for expected_item in expected_items
        )
    if mode == "contains_all":
        return all(
            any(_scalar_matches(actual_item, expected_item) for actual_item in actual_items)
            for expected_item in expected_items
        )
    if mode == "contains_exact":
        return len(actual_items) == len(expected_items) and _collection_members_match(
            actual_items, expected_items, "contains_all"
        )
    raise SchemaValidationError(f"unsupported collection condition operator {mode!r}")


def _value_matches(actual: Any, expected: Any) -> bool:
    """Compare a validated declarative condition without executing code."""

    if isinstance(expected, Mapping):
        operator, operand = next(iter(expected.items()))
        if operator == "in":
            return any(_value_matches(actual, item) for item in operand)
        if operator == "lt":
            return actual < operand
        if operator == "lte":
            return actual <= operand
        if operator == "gt":
            return actual > operand
        if operator == "gte":
            return actual >= operand
        if operator == "between":
            lower, upper = operand
            return lower <= actual <= upper
        if operator in {"contains_any", "contains_all", "contains_exact"}:
            return _collection_members_match(actual, operand, operator)
        raise SchemaValidationError(f"unsupported VCEP rule condition operator {operator!r}")
    if isinstance(expected, (list, tuple, set)):
        return any(_value_matches(actual, item) for item in expected)
    if isinstance(actual, (tuple, list, set)):
        # Scalar splice-outcome conditions preserve the historic "any credible
        # outcome" semantics. Collection-level intent must use contains_*.
        return any(_value_matches(item, expected) for item in actual)
    return _scalar_matches(actual, expected)


def _scope_match(scope: VCEPScope, evidence: BiologicalFeatures) -> tuple[MatchState, tuple[str, ...]]:
    required = {
        "gene_id": scope.gene_id,
        "disease_id": scope.disease_id,
        "transcript_id": scope.transcript_id,
        "transcript_version": scope.transcript_version,
        "mode_of_inheritance": scope.mode_of_inheritance,
    }
    unresolved: list[str] = []
    for field, expected in required.items():
        if expected is None:
            continue
        value = evidence.evidence(field)
        if value.status == EvidenceStatus.KNOWN and not _value_matches(value.value, expected):
            return MatchState.NO_MATCH, ()
        if value.status != EvidenceStatus.KNOWN:
            unresolved.append(field)
    if unresolved:
        return MatchState.UNRESOLVED, tuple(unresolved)
    return MatchState.MATCH, ()


def _rule_match(rule: VCEPDecisionRule, evidence: BiologicalFeatures) -> tuple[MatchState, tuple[str, ...]]:
    unresolved: list[str] = []
    for field, expected in rule.when.items():
        value = evidence.evidence(field)
        if value.status == EvidenceStatus.KNOWN and not _value_matches(value.value, expected):
            return MatchState.NO_MATCH, ()
        if value.status != EvidenceStatus.KNOWN:
            unresolved.append(field)
    if unresolved:
        return MatchState.UNRESOLVED, tuple(unresolved)
    return MatchState.MATCH, ()


def _known_rule_anchor_matches(
    rule: VCEPDecisionRule,
    evidence: BiologicalFeatures,
    field: str,
    allowed_values: Sequence[Any],
) -> bool:
    """Whether a rule's own condition matches an observed compatible anchor."""

    expected = rule.when.get(field)
    value = evidence.evidence(field)
    return (
        expected is not None
        and value.status == EvidenceStatus.KNOWN
        and _value_matches(value.value, expected)
        and any(_value_matches(value.value, allowed) for allowed in allowed_values)
    )


def _has_actual_supported_anchor(
    rule: VCEPDecisionRule, evidence: BiologicalFeatures
) -> bool:
    """Require an early ANY rule to match a real supported PVS1 class."""

    return any(
        (
            _known_rule_anchor_matches(
                rule,
                evidence,
                "variant_consequence",
                (
                    "STOP_GAIN",
                    "FRAMESHIFT",
                    "CANONICAL_SPLICE",
                    "SPLICE_SITE",
                    "START_LOST",
                    "EXON_DELETION",
                    "WHOLE_GENE_DELETION",
                    "DUPLICATION",
                ),
            ),
            _known_rule_anchor_matches(rule, evidence, "cnv_event_type", ("DELETION", "DUPLICATION")),
            _known_rule_anchor_matches(rule, evidence, "start_codon_affected", (True,)),
            _known_rule_anchor_matches(rule, evidence, "whole_gene_deletion", (True,)),
            _known_rule_anchor_matches(rule, evidence, "whole_gene_duplication", (True,)),
        )
    )


def _has_actual_compatible_pathway_anchor(
    rule: VCEPDecisionRule, evidence: BiologicalFeatures
) -> bool:
    """Ensure an early path-specific rule matches its actual biological class."""

    consequence = lambda *values: _known_rule_anchor_matches(
        rule, evidence, "variant_consequence", values
    )
    event = lambda *values: _known_rule_anchor_matches(
        rule, evidence, "cnv_event_type", values
    )
    position = lambda: _known_rule_anchor_matches(
        rule, evidence, "splice_site_position", (evidence.splice_site_position.value,)
    )
    affected_exons = lambda: (
        "affected_exons" in rule.when
        and evidence.affected_exons.status == EvidenceStatus.KNOWN
        and _value_matches(
            evidence.affected_exons.value,
            rule.when["affected_exons"],
        )
    )
    if rule.pathway == "NONSENSE_OR_FRAMESHIFT":
        return consequence("STOP_GAIN", "FRAMESHIFT")
    if rule.pathway == "CANONICAL_SPLICE":
        if consequence("CANONICAL_SPLICE"):
            return True
        splice_position = evidence.splice_site_position
        return (
            consequence("SPLICE_SITE")
            and splice_position.status == EvidenceStatus.KNOWN
            and abs(splice_position.value) in {1, 2}
        ) or (
            splice_position.status == EvidenceStatus.KNOWN
            and abs(splice_position.value) in {1, 2}
            and position()
        )
    if rule.pathway == "NONCANONICAL_SPLICE":
        splice_position = evidence.splice_site_position
        return (
            consequence("SPLICE_SITE")
            and splice_position.status == EvidenceStatus.KNOWN
            and abs(splice_position.value) not in {1, 2}
        ) or (
            splice_position.status == EvidenceStatus.KNOWN
            and abs(splice_position.value) not in {1, 2}
            and position()
        )
    if rule.pathway == "INITIATION_CODON":
        return consequence("START_LOST") or _known_rule_anchor_matches(
            rule, evidence, "start_codon_affected", (True,)
        )
    if rule.pathway == "WHOLE_GENE_DELETION":
        return consequence("WHOLE_GENE_DELETION") or _known_rule_anchor_matches(
            rule, evidence, "whole_gene_deletion", (True,)
        )
    if rule.pathway == "EXON_DELETION":
        return consequence("EXON_DELETION") or event("DELETION") or affected_exons()
    if rule.pathway == "DUPLICATION":
        return (
            consequence("DUPLICATION")
            or event("DUPLICATION")
            or affected_exons()
            or _known_rule_anchor_matches(
                rule, evidence, "whole_gene_duplication", (True,)
            )
        )
    return False


def select_specific_rule(
    evidence: BiologicalFeatures,
    pathway: Optional[str] = None,
    *,
    only_any: bool = False,
    permit_anchored_path_specific: bool = False,
) -> SpecificRuleSelection:
    """Select one exact, scoped VCEP rule or return a safe unresolved result.

    An unknown condition in a matching VCEP scope prevents falling back to the
    generic rule, because doing so could bypass the higher-precedence document.
    ``only_any`` is used before generic pathway identification when its branch
    cannot yet be inferred. It permits a deliberately broad ``ANY`` rule to
    resolve incomplete generic annotation, while preventing an unrelated
    path-specific rule from being applied solely because its local ``when``
    conditions happen to match. ``permit_anchored_path_specific`` additionally
    admits a path-specific rule only when its `when` expression declares and
    matches a pathway-compatible biological anchor: a variant/CNV/start fact,
    a canonical/non-canonical splice position, or affected exons for an
    exon-level structural path.

    Within one versioned specification, a path-specific rule is more specific
    than that document's ``ANY`` fallback. Rules from separate caller-supplied
    specifications have no implicit authority, version, or scope precedence:
    overlapping matches (or an unresolved competing rule) are returned rather
    than silently choosing one clinical declaration.
    """

    supplied = evidence.gene_specific_pvs1_rules
    if supplied.status == EvidenceStatus.UNKNOWN:
        return SpecificRuleSelection()
    if supplied.status == EvidenceStatus.CONFLICTING:
        return SpecificRuleSelection(
            conflict="gene-specific PVS1 specification evidence is conflicting"
        )
    specs = tuple(
        VCEPSpecification.from_mapping(raw, index)
        for index, raw in enumerate(supplied.value or ())
    )
    pathway_token = _token(pathway) if pathway is not None else None

    def eligible(rule: VCEPDecisionRule) -> bool:
        if only_any:
            if rule.pathway == "ANY":
                return _has_actual_supported_anchor(rule, evidence)
            return bool(
                permit_anchored_path_specific
                and _has_actual_compatible_pathway_anchor(rule, evidence)
            )
        return pathway_token is None or rule.pathway in {pathway_token, "ANY"}

    matched_specs: list[
        tuple[
            VCEPSpecification,
            tuple[tuple[VCEPDecisionRule, MatchState, tuple[str, ...]], ...],
        ]
    ] = []
    unresolved_scope: list[str] = []
    for spec in specs:
        eligible_rules = tuple(rule for rule in spec.rules if eligible(rule))
        if not eligible_rules:
            continue
        # Check condition feasibility before letting an incomplete scope block
        # generic evaluation. A known condition mismatch proves this rule
        # cannot apply regardless of missing scope fields.
        rule_assessments = tuple(
            (rule, *_rule_match(rule, evidence)) for rule in eligible_rules
        )
        if all(state == MatchState.NO_MATCH for _, state, _ in rule_assessments):
            continue
        scope_state, scope_missing = _scope_match(spec.scope, evidence)
        if scope_state == MatchState.MATCH:
            matched_specs.append((spec, rule_assessments))
        elif scope_state == MatchState.UNRESOLVED:
            unresolved_scope.extend(scope_missing)

    # An unresolved scope matters only when the declaration has a rule that
    # could apply to the route being evaluated.
    if unresolved_scope:
        return SpecificRuleSelection(missing=tuple(dict.fromkeys(unresolved_scope)))

    chosen_rules: list[tuple[VCEPSpecification, VCEPDecisionRule]] = []
    unresolved_rules: list[str] = []
    for spec, rule_assessments in matched_specs:
        specific_matches: list[VCEPDecisionRule] = []
        any_matches: list[VCEPDecisionRule] = []
        specific_missing: list[str] = []
        any_missing: list[str] = []
        for rule, state, missing in rule_assessments:
            if state == MatchState.MATCH:
                (any_matches if rule.pathway == "ANY" else specific_matches).append(rule)
            elif state == MatchState.UNRESOLVED:
                (any_missing if rule.pathway == "ANY" else specific_missing).extend(
                    missing
                )

        # A path-specific VCEP rule is more specific than an ANY fallback only
        # inside this one source document. An unresolved path-specific rule
        # prevents its own fallback from being selected.
        if specific_missing:
            unresolved_rules.extend(specific_missing)
            continue
        if len(specific_matches) > 1:
            return SpecificRuleSelection(
                conflict="multiple equally applicable VCEP PVS1 rules match within one specification"
            )
        if len(specific_matches) == 1:
            chosen_rules.append((spec, specific_matches[0]))
            continue
        if any_missing:
            unresolved_rules.extend(any_missing)
            continue
        if len(any_matches) > 1:
            return SpecificRuleSelection(
                conflict="multiple equally applicable ANY VCEP PVS1 rules match within one specification"
            )
        if len(any_matches) == 1:
            chosen_rules.append((spec, any_matches[0]))

    # A potentially applicable condition in a separate specification cannot
    # be bypassed by an exact match elsewhere: callers have supplied no
    # machine-readable authority or cross-document precedence relationship.
    if unresolved_rules:
        return SpecificRuleSelection(missing=tuple(dict.fromkeys(unresolved_rules)))
    if len(chosen_rules) > 1:
        return SpecificRuleSelection(
            conflict="multiple separately supplied VCEP PVS1 specifications match this evidence"
        )
    if len(chosen_rules) == 1:
        return SpecificRuleSelection(matched=chosen_rules[0])
    return SpecificRuleSelection()


def guideline_for_specific_rule(
    spec: VCEPSpecification, rule: VCEPDecisionRule
) -> GuidelineReference:
    return GuidelineReference(
        authority=spec.authority,
        document=spec.name,
        version=spec.version,
        source_url=spec.source_url,
        source_locator=rule.source_locator,
        rule_id=spec.specification_id,
    )
