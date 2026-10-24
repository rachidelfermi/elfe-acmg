"""Versioned, scoped ClinGen VCEP overrides for PS2/PM6.

The engine intentionally ships no partial catalogue of living VCEP
specifications.  A caller supplies the exact, versioned document snapshot that
is applicable to the gene/disease being interpreted.  Rules are declarative;
they can adjust phenotype point allocations, pool thresholds, restrictions, or
validate the final code label without executing caller-provided code.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
from typing import Any, Callable, Mapping, Optional, Sequence

from acmg.criteria.pvs1.models import EvidenceStatus, GuidelineReference, SchemaValidationError

from .models import (
    DeNovoObservation,
    DeNovoStrength,
    InheritancePattern,
    ObservationStatus,
    ParentageStatus,
    PhenotypeCategory,
    _OBSERVATION_ALIASES,
    _parse_adequacy,
    _parse_bool,
    _parse_genotype,
    _parse_phenotype,
    _parse_validity,
    _token,
)


DEFAULT_SVI_STRENGTH_THRESHOLDS: Mapping[DeNovoStrength, float] = {
    DeNovoStrength.SUPPORTING: 0.5,
    DeNovoStrength.MODERATE: 1.0,
    DeNovoStrength.STRONG: 2.0,
    DeNovoStrength.VERY_STRONG: 4.0,
}


# Keep condition parsing and condition matching on the same vocabulary as the
# biological input schema.  A VCEP author can therefore use a documented input
# alias such as ``NEGATIVE`` or ``GOOD`` without creating a rule that validates
# successfully but can never match a normalized observation.
_BOOLEAN_CONDITION_FIELDS = frozenset({
    "case_independence", "proband_variant_present", "maternity_confirmed",
    "paternity_confirmed", "parentage_disproven", "sample_identity_confirmation",
    "parental_mosaicism_evidence", "proband_mosaicism_evidence",
    "apparent_germline_mosaicism", "phenotype_consistency",
    "phenotype_heterogeneity", "alternative_molecular_diagnosis",
    "additional_pathogenic_variant_present", "family_history_consistent",
})
_GENOTYPE_CONDITION_FIELDS = frozenset({
    "proband_genotype", "maternal_genotype", "paternal_genotype",
})
_ADEQUACY_CONDITION_FIELDS = frozenset({
    "maternal_variant_absence_confidence", "paternal_variant_absence_confidence",
    "parental_testing_adequacy", "parentage_confirmation_quality", "variant_call_quality",
})
_STRING_CONDITION_FIELDS = frozenset({
    "case_id", "publication_case_id", "parentage_confirmation_method",
    "family_relationship_structure", "gene_id", "disease_id",
})


class MatchState(str, Enum):
    MATCH = "MATCH"
    NO_MATCH = "NO_MATCH"
    UNRESOLVED = "UNRESOLVED"


def _reject_unknown_keys(raw: Mapping[str, Any], allowed: set[str], context: str) -> None:
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise SchemaValidationError(f"{context} has unsupported field(s): " + ", ".join(unknown))


def _aliased_value(
    raw: Mapping[str, Any],
    canonical: str,
    alias: str,
    context: str,
    *,
    normalize: Optional[Callable[[Any], Any]] = None,
) -> Any:
    """Read a canonical/legacy pair without silently choosing one value."""

    canonical_present = canonical in raw
    alias_present = alias in raw
    if canonical_present and alias_present:
        left = raw[canonical]
        right = raw[alias]
        normalized_left = normalize(left) if normalize is not None else left
        normalized_right = normalize(right) if normalize is not None else right
        if normalized_left != normalized_right:
            raise SchemaValidationError(f"{context}.{canonical} conflicts with {alias}")
    if canonical_present:
        return raw[canonical]
    return raw.get(alias)


def _canonical_parentage_status(raw: Any, context: str) -> ParentageStatus:
    token = _token(raw)
    aliases = {
        "UNCONFIRMED": ParentageStatus.ASSUMED,
        "CONFIRMED_PARENTAGE": ParentageStatus.CONFIRMED,
    }
    if token in aliases:
        return aliases[token]
    try:
        result = ParentageStatus(token)
    except ValueError as error:
        raise SchemaValidationError(f"{context} is not a supported parentage status") from error
    if result not in {ParentageStatus.CONFIRMED, ParentageStatus.ASSUMED}:
        raise SchemaValidationError(f"{context} is not a supported parentage status")
    return result


def _normalize_condition_value(raw: Any, context: str, *, field: Optional[str] = None) -> Any:
    """Validate and canonicalize declarative VCEP condition grammar.

    Condition values are stored in their schema-normalized form.  This avoids
    the unsafe situation where a manifest accepts an upstream alias at parse
    time but compares that raw alias against a normalized enum at evaluation
    time.
    """

    if isinstance(raw, Mapping):
        if set(raw) != {"in"}:
            raise SchemaValidationError(
                f"{context} only supports a scalar, list, or {{'in': [...]}} condition"
            )
        values = raw["in"]
        if not isinstance(values, (list, tuple, set)) or not values:
            raise SchemaValidationError(f"{context}.in must be a non-empty list")
        return {
            "in": tuple(
                _normalize_condition_value(value, f"{context}.in[{index}]", field=field)
                for index, value in enumerate(values)
            )
        }
    if isinstance(raw, (list, tuple, set)):
        if not raw:
            raise SchemaValidationError(f"{context} cannot be an empty list")
        return tuple(
            _normalize_condition_value(value, f"{context}[{index}]", field=field)
            for index, value in enumerate(raw)
        )
    if raw is None or isinstance(raw, (Mapping, list, tuple, set)):
        raise SchemaValidationError(f"{context} must be a scalar condition value")
    if isinstance(raw, float) and not math.isfinite(raw):
        raise SchemaValidationError(f"{context} must not be NaN or Infinity")
    if field in _BOOLEAN_CONDITION_FIELDS:
        return _parse_bool(raw, context)
    if field in _GENOTYPE_CONDITION_FIELDS:
        return _parse_genotype(raw, context)
    if field in _ADEQUACY_CONDITION_FIELDS:
        return _parse_adequacy(raw, context)
    if field == "gene_disease_validity":
        return _parse_validity(raw, context)
    if field in _STRING_CONDITION_FIELDS:
        if not isinstance(raw, str) or not raw.strip():
            raise SchemaValidationError(f"{context} must be a non-empty string")
        return raw.strip()
    if field == "number_of_independent_de_novo_cases":
        if isinstance(raw, bool) or not isinstance(raw, int) or raw < 0:
            raise SchemaValidationError(f"{context} must be a non-negative integer")
        return raw
    if field == "inheritance_pattern":
        return _canonical_inheritance(raw, context)
    if field in {"parentage_status", "independent_case_parentage_status"}:
        return _canonical_parentage_status(raw, context)
    if field in {"phenotype_category", "phenotype_specificity"}:
        return _parse_phenotype(raw, context)
    return raw


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


def _number(raw: Any, context: str) -> float:
    if isinstance(raw, bool) or not isinstance(raw, (int, float)):
        raise SchemaValidationError(f"{context} must be a non-negative number")
    result = float(raw)
    if not math.isfinite(result) or result < 0:
        raise SchemaValidationError(f"{context} must be a finite non-negative number")
    return result


def _strength(raw: Any, context: str) -> DeNovoStrength:
    token = _token(raw)
    try:
        return DeNovoStrength(token)
    except ValueError as error:
        raise SchemaValidationError(f"{context} is not a PS2/PM6 strength") from error


def _canonical_inheritance(raw: Any, context: str) -> str:
    """Normalize VCEP scope abbreviations to the observation enum spelling."""

    if not isinstance(raw, str) or not raw.strip():
        raise SchemaValidationError(f"{context} must be a non-empty inheritance string")
    token = _token(raw)
    aliases = {
        "AD": InheritancePattern.AUTOSOMAL_DOMINANT,
        "AR": InheritancePattern.AUTOSOMAL_RECESSIVE,
        "XL": InheritancePattern.X_LINKED,
        "XLD": InheritancePattern.X_LINKED_DOMINANT,
        "XLR": InheritancePattern.X_LINKED_RECESSIVE,
    }
    if token in aliases:
        return aliases[token].value
    for member in InheritancePattern:
        if token in {_token(member.name), _token(member.value)}:
            return member.value
    raise SchemaValidationError(f"{context} is not a supported inheritance pattern")


@dataclass(frozen=True)
class VCEPScope:
    gene_id: str
    disease_id: str
    mode_of_inheritance: Optional[str] = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any], context: str) -> "VCEPScope":
        _reject_unknown_keys(
            raw,
            {"gene_id", "disease_id", "mode_of_inheritance", "inheritance_pattern"},
            context,
        )
        mode_of_inheritance = _aliased_value(
            raw,
            "mode_of_inheritance",
            "inheritance_pattern",
            context,
            normalize=lambda value: _canonical_inheritance(
                value, f"{context}.mode_of_inheritance"
            ),
        )
        if mode_of_inheritance is not None:
            mode_of_inheritance = _canonical_inheritance(
                mode_of_inheritance, f"{context}.mode_of_inheritance"
            )
        return cls(
            gene_id=_required_string(raw, "gene_id", context),
            disease_id=_required_string(raw, "disease_id", context),
            mode_of_inheritance=mode_of_inheritance,
        )


@dataclass(frozen=True)
class ScoringPolicyOverride:
    """Only documented VCEP changes and label assertions are represented."""

    points: Mapping[PhenotypeCategory, Mapping[ParentageStatus, float]]
    thresholds: Mapping[DeNovoStrength, float]
    heterogeneous_points_cap: Optional[float] = None
    assigned_criterion: Optional[str] = None

    @classmethod
    def empty(cls) -> "ScoringPolicyOverride":
        return cls(points={}, thresholds={})

    @property
    def changes_final_policy(self) -> bool:
        return bool(self.thresholds or self.heterogeneous_points_cap is not None)


def _parse_points(raw: Any, context: str) -> Mapping[PhenotypeCategory, Mapping[ParentageStatus, float]]:
    if raw is None:
        return {}
    if not isinstance(raw, Mapping):
        raise SchemaValidationError(f"{context} must be an object")
    categories: dict[PhenotypeCategory, Mapping[ParentageStatus, float]] = {}
    for raw_category, values in raw.items():
        try:
            category = _parse_phenotype(raw_category, f"{context}.{raw_category}")
        except SchemaValidationError as error:
            raise SchemaValidationError(f"{context} has invalid phenotype category {raw_category!r}") from error
        if category in categories:
            raise SchemaValidationError(
                f"{context} repeats the normalized phenotype category {category.value!r}"
            )
        if not isinstance(values, Mapping):
            raise SchemaValidationError(f"{context}.{raw_category} must map parentage statuses to points")
        point_values: dict[ParentageStatus, float] = {}
        for raw_parentage, raw_value in values.items():
            token = _token(raw_parentage)
            aliases = {"UNCONFIRMED": ParentageStatus.ASSUMED, "CONFIRMED_PARENTAGE": ParentageStatus.CONFIRMED}
            try:
                parentage = aliases[token] if token in aliases else ParentageStatus(token)
            except ValueError as error:
                raise SchemaValidationError(
                    f"{context}.{raw_category} has invalid parentage status {raw_parentage!r}"
                ) from error
            if parentage not in {ParentageStatus.CONFIRMED, ParentageStatus.ASSUMED}:
                raise SchemaValidationError(f"{context}.{raw_category} may only set CONFIRMED or ASSUMED points")
            if parentage in point_values:
                raise SchemaValidationError(
                    f"{context}.{raw_category} repeats the normalized parentage status {parentage.value!r}"
                )
            point_values[parentage] = _number(raw_value, f"{context}.{raw_category}.{raw_parentage}")
        categories[category] = point_values
    return categories


def _parse_thresholds(raw: Any, context: str) -> Mapping[DeNovoStrength, float]:
    if raw is None:
        return {}
    if not isinstance(raw, Mapping):
        raise SchemaValidationError(f"{context} must be an object")
    values: dict[DeNovoStrength, float] = {}
    allowed = {
        DeNovoStrength.SUPPORTING,
        DeNovoStrength.MODERATE,
        DeNovoStrength.STRONG,
        DeNovoStrength.VERY_STRONG,
    }
    for raw_strength, raw_value in raw.items():
        strength = _strength(raw_strength, f"{context}.{raw_strength}")
        if strength not in allowed:
            raise SchemaValidationError(f"{context} cannot set a NONE threshold")
        if strength in values:
            raise SchemaValidationError(
                f"{context} repeats the normalized strength {strength.value!r}"
            )
        value = _number(raw_value, f"{context}.{raw_strength}")
        if value <= 0:
            raise SchemaValidationError(f"{context}.{raw_strength} must be greater than zero")
        values[strength] = value
    ordered = [
        DeNovoStrength.SUPPORTING,
        DeNovoStrength.MODERATE,
        DeNovoStrength.STRONG,
        DeNovoStrength.VERY_STRONG,
    ]
    supplied = [values[item] for item in ordered if item in values]
    if supplied and supplied != sorted(supplied):
        raise SchemaValidationError(f"{context} must increase from Supporting to Very Strong")
    return values


def _validate_effective_thresholds(
    overrides: Mapping[DeNovoStrength, float],
    context: str,
) -> None:
    """Ensure a partial VCEP override remains valid after generic defaults merge."""

    effective = dict(DEFAULT_SVI_STRENGTH_THRESHOLDS)
    effective.update(overrides)
    ordered = (
        DeNovoStrength.SUPPORTING,
        DeNovoStrength.MODERATE,
        DeNovoStrength.STRONG,
        DeNovoStrength.VERY_STRONG,
    )
    values = [effective[strength] for strength in ordered]
    if any(value <= 0 for value in values) or any(
        lower >= higher for lower, higher in zip(values, values[1:])
    ):
        raise SchemaValidationError(
            f"{context} produces an invalid effective Supporting-to-Very Strong threshold ladder"
        )


@dataclass(frozen=True)
class VCEPDeNovoRule:
    source_locator: str
    when: Mapping[str, Any]
    scoring: ScoringPolicyOverride
    exclude: bool = False
    exclusion_status: ObservationStatus = ObservationStatus.NOT_MET
    exclusion_reason: Optional[str] = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any], context: str) -> "VCEPDeNovoRule":
        if not isinstance(raw, Mapping):
            raise SchemaValidationError(f"{context} must be an object")
        _reject_unknown_keys(
            raw,
            {
                "source_locator", "when", "point_overrides", "points",
                "strength_thresholds", "thresholds",
                "heterogeneous_points_cap", "heterogeneous_point_cap",
                "assigned_criterion", "criterion", "exclude", "exclusion_reason", "decision",
            },
            context,
        )
        source_locator = _required_string(raw, "source_locator", context)
        when_raw = raw.get("when", {})
        if not isinstance(when_raw, Mapping):
            raise SchemaValidationError(f"{context}.when must be an object")
        allowed_fields = set(DeNovoObservation.__dataclass_fields__) | {
            "parentage_status",
            "phenotype_category",
        }
        normalized_when: dict[str, Any] = {}
        for raw_field, value in when_raw.items():
            field = _OBSERVATION_ALIASES.get(raw_field, raw_field)
            if field not in allowed_fields:
                raise SchemaValidationError(
                    f"{context}.when refers to unsupported feature {raw_field!r}"
                )
            normalized_value = _normalize_condition_value(
                value, f"{context}.when.{raw_field}", field=field
            )
            if field in normalized_when and normalized_when[field] != normalized_value:
                raise SchemaValidationError(
                    f"{context}.when.{raw_field} conflicts with its canonical {field!r} condition"
                )
            normalized_when[field] = normalized_value
        points_raw = _aliased_value(raw, "point_overrides", "points", context)
        thresholds_raw = _aliased_value(raw, "strength_thresholds", "thresholds", context)
        cap_raw = _aliased_value(
            raw, "heterogeneous_points_cap", "heterogeneous_point_cap", context
        )
        label = _aliased_value(
            raw,
            "assigned_criterion",
            "criterion",
            context,
            normalize=lambda value: _token(value) if value is not None else None,
        )
        if label is not None:
            label = _token(label)
            if label not in {"PS2", "PM6"}:
                raise SchemaValidationError(f"{context}.assigned_criterion must be PS2 or PM6")
        decision = raw.get("decision")
        exclude_raw = raw.get("exclude", False)
        if not isinstance(exclude_raw, bool):
            raise SchemaValidationError(f"{context}.exclude must be a Boolean")
        exclude = exclude_raw
        exclusion_status = ObservationStatus.NOT_MET
        exclusion_reason = raw.get("exclusion_reason")
        if decision is not None:
            if not isinstance(decision, Mapping):
                raise SchemaValidationError(f"{context}.decision must be an object")
            _reject_unknown_keys(decision, {"status"}, f"{context}.decision")
            if "status" not in decision:
                raise SchemaValidationError(f"{context}.decision.status is required")
            status = _token(decision["status"])
            if status in {"NOT_MET", "NOT_APPLICABLE"}:
                if "exclude" in raw and exclude_raw is False:
                    raise SchemaValidationError(f"{context}.exclude conflicts with {context}.decision")
                exclude = True
                exclusion_status = ObservationStatus(status)
                exclusion_reason = exclusion_reason or f"VCEP rule returns {status}"
            elif status:
                raise SchemaValidationError(
                    f"{context}.decision supports only NOT_MET or NOT_APPLICABLE restrictions; "
                    "use point_overrides/thresholds for scoring changes"
                )
        if exclusion_reason is not None and (not isinstance(exclusion_reason, str) or not exclusion_reason.strip()):
            raise SchemaValidationError(f"{context}.exclusion_reason must be a non-empty string")
        if exclusion_reason is not None and not exclude:
            raise SchemaValidationError(f"{context}.exclusion_reason requires an exclusion decision")
        points = _parse_points(points_raw, f"{context}.point_overrides")
        thresholds = _parse_thresholds(thresholds_raw, f"{context}.strength_thresholds")
        _validate_effective_thresholds(thresholds, f"{context}.strength_thresholds")
        return cls(
            source_locator=source_locator,
            when=normalized_when,
            scoring=ScoringPolicyOverride(
                points=points,
                thresholds=thresholds,
                heterogeneous_points_cap=_number(cap_raw, f"{context}.heterogeneous_points_cap")
                if cap_raw is not None
                else None,
                assigned_criterion=label,
            ),
            exclude=exclude,
            exclusion_status=exclusion_status,
            exclusion_reason=exclusion_reason.strip() if isinstance(exclusion_reason, str) else None,
        )


@dataclass(frozen=True)
class VCEPDeNovoSpecification:
    authority: str
    name: str
    version: str
    source_url: str
    scope: VCEPScope
    rules: tuple[VCEPDeNovoRule, ...]
    specification_id: Optional[str] = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any], index: int) -> "VCEPDeNovoSpecification":
        if not isinstance(raw, Mapping):
            raise SchemaValidationError("applicable_vcep_de_novo_rules entries must be objects")
        context = f"applicable_vcep_de_novo_rules[{index}]"
        scope_fields = {
            "gene_id", "disease_id", "mode_of_inheritance", "inheritance_pattern",
        }
        _reject_unknown_keys(
            raw,
            {
                "authority", "name", "version", "source_url", "scope", "specification_id", "rules",
                *scope_fields,
            },
            context,
        )
        name = _required_string(raw, "name", context)
        version = _required_string(raw, "version", context)
        source_url = _required_string(raw, "source_url", context)
        if not source_url.startswith(("https://", "http://")):
            raise SchemaValidationError(f"{context}.source_url must be an HTTP(S) URL")
        nested_scope = raw.get("scope")
        if nested_scope is not None and not isinstance(nested_scope, Mapping):
            raise SchemaValidationError(f"{context}.scope must be an object")
        scope_raw: dict[str, Any] = dict(nested_scope or {})
        for field in scope_fields:
            if field not in raw:
                continue
            if field in scope_raw:
                left = scope_raw[field]
                right = raw[field]
                if field in {"mode_of_inheritance", "inheritance_pattern"}:
                    equal = _canonical_inheritance(left, f"{context}.scope.{field}") == _canonical_inheritance(
                        right, f"{context}.{field}"
                    )
                elif field in {"gene_id", "disease_id"}:
                    equal = _identifier_matches(left, right)
                else:
                    equal = left == right
                if not equal:
                    raise SchemaValidationError(f"{context}.scope.{field} conflicts with {context}.{field}")
            else:
                scope_raw[field] = raw[field]
        rules_raw = raw.get("rules")
        if not isinstance(rules_raw, (list, tuple)) or not rules_raw:
            raise SchemaValidationError(f"{context}.rules must be a non-empty list")
        authority = raw.get("authority", "ClinGen VCEP")
        if not isinstance(authority, str) or not authority.strip():
            raise SchemaValidationError(f"{context}.authority must be a non-empty string")
        specification_id = _optional_string(raw, "specification_id", context)
        return cls(
            authority=authority.strip(),
            name=name,
            version=version,
            source_url=source_url,
            scope=VCEPScope.from_mapping(scope_raw, f"{context}.scope"),
            rules=tuple(
                VCEPDeNovoRule.from_mapping(item, f"{context}.rules[{rule_index}]")
                for rule_index, item in enumerate(rules_raw)
            ),
            specification_id=specification_id,
        )


def parse_specifications(raw_specs: Sequence[Mapping[str, Any]]) -> tuple[VCEPDeNovoSpecification, ...]:
    return tuple(VCEPDeNovoSpecification.from_mapping(raw, index) for index, raw in enumerate(raw_specs))


def _value_matches(actual: Any, expected: Any) -> bool:
    if isinstance(expected, Mapping):
        if set(expected) != {"in"}:
            raise SchemaValidationError(
                "VCEP rule conditions only support a scalar, list, or {'in': [...]}"
            )
        expected = expected["in"]
    if isinstance(expected, (list, tuple, set)):
        return any(_value_matches(actual, item) for item in expected)
    actual_value = actual.value if isinstance(actual, Enum) else actual
    expected_value = expected.value if isinstance(expected, Enum) else expected
    if isinstance(actual_value, str) and isinstance(expected_value, str):
        return _token(actual_value) == _token(expected_value)
    return actual_value == expected_value


def _identifier_matches(actual: Any, expected: Any) -> bool:
    """Match stable identifiers without collapsing meaningful punctuation."""

    if not isinstance(actual, str) or not isinstance(expected, str):
        return actual == expected
    return actual.strip().casefold() == expected.strip().casefold()


def _field_value_matches(field: str, actual: Any, expected: Any) -> bool:
    if isinstance(expected, Mapping):
        if set(expected) != {"in"}:
            raise SchemaValidationError(
                "VCEP rule conditions only support a scalar, list, or {'in': [...]}"
            )
        return any(_field_value_matches(field, actual, item) for item in expected["in"])
    if isinstance(expected, (list, tuple, set)):
        return any(_field_value_matches(field, actual, item) for item in expected)
    if field in {"gene_id", "disease_id"}:
        return _identifier_matches(actual, expected)
    if field == "inheritance_pattern":
        actual_value = actual.value if isinstance(actual, Enum) else actual
        expected_value = expected.value if isinstance(expected, Enum) else expected
        if not isinstance(actual_value, str) or not isinstance(expected_value, str):
            return False
        return _canonical_inheritance(actual_value, "observed inheritance pattern") == _canonical_inheritance(
            expected_value, "VCEP inheritance condition"
        )
    if field == "parentage_status":
        actual_value = actual.value if isinstance(actual, Enum) else actual
        expected_value = expected.value if isinstance(expected, Enum) else expected
        aliases = {"UNCONFIRMED": "ASSUMED", "CONFIRMED_PARENTAGE": "CONFIRMED"}
        if not isinstance(actual_value, str) or not isinstance(expected_value, str):
            return False
        return _token(actual_value) == aliases.get(_token(expected_value), _token(expected_value))
    if isinstance(actual, bool) and isinstance(expected, str):
        token = expected.strip().upper()
        if token in {"TRUE", "YES"}:
            return actual is True
        if token in {"FALSE", "NO"}:
            return actual is False
    return _value_matches(actual, expected)


def _scope_match(scope: VCEPScope, observation: DeNovoObservation) -> tuple[MatchState, tuple[str, ...]]:
    expected = {
        "gene_id": scope.gene_id,
        "disease_id": scope.disease_id,
        "inheritance_pattern": scope.mode_of_inheritance,
    }
    unresolved: list[str] = []
    for field, target in expected.items():
        if target is None:
            continue
        value = observation.evidence(field)
        if value.status == EvidenceStatus.CONFLICTING:
            unresolved.append(field)
            continue
        if value.status == EvidenceStatus.UNKNOWN:
            unresolved.append(field)
            continue
        if not _field_value_matches(field, value.value, target):
            return MatchState.NO_MATCH, ()
    return (MatchState.UNRESOLVED, tuple(unresolved)) if unresolved else (MatchState.MATCH, ())


def _rule_match(
    rule: VCEPDeNovoRule,
    observation: DeNovoObservation,
    derived: Mapping[str, Any],
) -> tuple[MatchState, tuple[str, ...]]:
    unresolved: list[str] = []
    for field, target in rule.when.items():
        if field in derived:
            actual = derived[field]
            if actual is None:
                unresolved.append(field)
                continue
        else:
            evidence = observation.evidence(field)
            if evidence.status == EvidenceStatus.CONFLICTING:
                unresolved.append(field)
                continue
            if evidence.status == EvidenceStatus.UNKNOWN:
                unresolved.append(field)
                continue
            actual = evidence.value
        if not _field_value_matches(field, actual, target):
            return MatchState.NO_MATCH, ()
    return (MatchState.UNRESOLVED, tuple(unresolved)) if unresolved else (MatchState.MATCH, ())


@dataclass(frozen=True)
class VCEPSelection:
    matched: Optional[tuple[VCEPDeNovoSpecification, VCEPDeNovoRule]] = None
    missing: tuple[str, ...] = ()
    conflict: Optional[str] = None


def select_vcep_rule(
    specifications: Sequence[VCEPDeNovoSpecification],
    observation: DeNovoObservation,
    derived: Mapping[str, Any],
) -> VCEPSelection:
    """Select one matching VCEP rule before falling back to generic SVI.

    An unresolved potentially matching scope or rule deliberately blocks generic
    scoring for that observation: falling back could bypass a higher-precedence
    gene/disease specification.
    """

    scoped: list[VCEPDeNovoSpecification] = []
    unresolved_scope: list[str] = []
    for specification in specifications:
        state, missing = _scope_match(specification.scope, observation)
        if state == MatchState.MATCH:
            scoped.append(specification)
        elif state == MatchState.UNRESOLVED:
            unresolved_scope.extend(missing)
    if unresolved_scope:
        return VCEPSelection(missing=tuple(dict.fromkeys(unresolved_scope)))

    matching: list[tuple[VCEPDeNovoSpecification, VCEPDeNovoRule]] = []
    unresolved_rules: list[str] = []
    for specification in scoped:
        for rule in specification.rules:
            state, missing = _rule_match(rule, observation, derived)
            if state == MatchState.MATCH:
                matching.append((specification, rule))
            elif state == MatchState.UNRESOLVED:
                unresolved_rules.extend(missing)
    # An unresolved condition in a scoped VCEP declaration may be a more
    # specific restriction than a broad matching row.  Prefer an auditable
    # pause to silently bypassing that higher-precedence restriction.
    if unresolved_rules:
        return VCEPSelection(missing=tuple(dict.fromkeys(unresolved_rules)))
    if len(matching) > 1:
        return VCEPSelection(conflict="multiple equally applicable VCEP PS2/PM6 rules match this observation")
    if len(matching) == 1:
        return VCEPSelection(matched=matching[0])
    return VCEPSelection()


def guideline_for_vcep_rule(
    specification: VCEPDeNovoSpecification,
    rule: VCEPDeNovoRule,
) -> GuidelineReference:
    return GuidelineReference(
        authority=specification.authority,
        document=specification.name,
        version=specification.version,
        source_url=specification.source_url,
        source_locator=rule.source_locator,
        rule_id=specification.specification_id,
    )


__all__ = [
    "ScoringPolicyOverride",
    "VCEPDeNovoRule",
    "VCEPDeNovoSpecification",
    "VCEPSelection",
    "guideline_for_vcep_rule",
    "parse_specifications",
    "select_vcep_rule",
]
