"""Typed, provenance-bearing input and output models for PS2/PM6.

The generic PS2/PM6 recommendation is deliberately evaluated from already
curated biological observations.  This module does *not* call variants, infer
parentage from sequencing, or turn a missing call into a reference call.

The common evidence envelope is shared with :mod:`pvs1`: a value is either
known, unknown, or conflicting.  Boolean values retain their truth value
inside a known envelope, so ``False`` is never confused with ``UNKNOWN``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import Enum
import math
import re
from typing import Any, Callable, Generic, Iterable, Mapping, Optional, Sequence, TypeVar

from acmg.criteria.pvs1.models import (
    DecisionStatus,
    EvidenceKind,
    EvidenceStatus,
    EvidenceValue,
    GeneDiseaseValidity,
    GuidelineReference,
    Provenance,
    SchemaValidationError,
    TextEnum,
    stable_unique,
)


T = TypeVar("T")
_MISSING = object()


_OBSERVATION_ALIASES: Mapping[str, str] = {
    "observation_id": "case_id",
    "patient_id": "case_id",
    "literature_case_id": "publication_case_id",
    "independent_case_phenotype_specificity": "phenotype_specificity",
    "mode_of_inheritance": "inheritance_pattern",
    "parental_test_adequacy": "parental_testing_adequacy",
    "technical_quality": "variant_call_quality",
    "second_pathogenic_variant_present": "additional_pathogenic_variant_present",
}


class DeNovoStrength(TextEnum):
    """Allowed ACMG/AMP strength modifiers for PS2 and PM6."""

    VERY_STRONG = "VERY_STRONG"
    STRONG = "STRONG"
    MODERATE = "MODERATE"
    SUPPORTING = "SUPPORTING"
    NONE = "NONE"


# A descriptive alias makes integrations written before the package name was
# settled readable without introducing a second enum type.
PS2PM6Strength = DeNovoStrength


class GenotypeState(TextEnum):
    """A call with respect to the variant being evaluated, not a zygosity API."""

    VARIANT_PRESENT = "VARIANT_PRESENT"
    VARIANT_ABSENT = "VARIANT_ABSENT"
    MOSAIC = "MOSAIC"
    NO_CALL = "NO_CALL"
    OTHER = "OTHER"


class TestAdequacy(TextEnum):
    """Curated adequacy conclusion; raw depth alone has no universal cutoff."""

    ADEQUATE = "ADEQUATE"
    INADEQUATE = "INADEQUATE"


class PhenotypeCategory(TextEnum):
    """The four ClinGen SVI Table 1 phenotype categories."""

    HIGHLY_SPECIFIC = "HIGHLY_SPECIFIC"
    CONSISTENT = "CONSISTENT"
    HETEROGENEOUS = "HETEROGENEOUS"
    INCONSISTENT = "INCONSISTENT"


class ParentageStatus(TextEnum):
    CONFIRMED = "CONFIRMED"
    ASSUMED = "ASSUMED"
    UNRESOLVED = "UNRESOLVED"
    INVALID = "INVALID"


def _parse_parentage_status(raw: Any, field_name: str) -> ParentageStatus:
    """Parse the per-case SVI confirmed-versus-assumed parentage assertion."""

    aliases = {
        "UNCONFIRMED": ParentageStatus.ASSUMED,
        "CONFIRMED_PARENTAGE": ParentageStatus.CONFIRMED,
    }
    result = _enum_value(raw, ParentageStatus, field_name, aliases)  # type: ignore[assignment]
    if result not in {ParentageStatus.CONFIRMED, ParentageStatus.ASSUMED}:
        raise SchemaValidationError(
            f"{field_name} must be CONFIRMED or ASSUMED when known"
        )
    return result


class InheritancePattern(TextEnum):
    AUTOSOMAL_DOMINANT = "AUTOSOMAL_DOMINANT"
    AUTOSOMAL_RECESSIVE = "AUTOSOMAL_RECESSIVE"
    X_LINKED = "X_LINKED"
    X_LINKED_DOMINANT = "X_LINKED_DOMINANT"
    X_LINKED_RECESSIVE = "X_LINKED_RECESSIVE"
    MITOCHONDRIAL = "MITOCHONDRIAL"
    OTHER = "OTHER"


class ObservationStatus(TextEnum):
    ELIGIBLE = "ELIGIBLE"
    NOT_MET = "NOT_MET"
    UNDETERMINED = "UNDETERMINED"
    NOT_APPLICABLE = "NOT_APPLICABLE"
    DUPLICATE = "DUPLICATE"


def _token(value: Any) -> str:
    if isinstance(value, Enum):
        value = value.value
    return re.sub(r"[^A-Z0-9]+", "_", str(value).upper()).strip("_")


def _enum_value(
    raw: Any,
    enum_type: type[TextEnum],
    field_name: str,
    aliases: Optional[Mapping[str, TextEnum]] = None,
) -> TextEnum:
    if isinstance(raw, enum_type):
        return raw
    if not isinstance(raw, str):
        raise SchemaValidationError(
            f"{field_name} must be a {enum_type.__name__} string, not {raw!r}"
        )
    token = _token(raw)
    if aliases and token in aliases:
        return aliases[token]
    for member in enum_type:
        if token in {_token(member.name), _token(member.value)}:
            return member
    choices = ", ".join(member.value for member in enum_type)
    raise SchemaValidationError(
        f"{field_name} has unsupported value {raw!r}; expected one of {choices}"
    )


def _parse_bool(raw: Any, field_name: str) -> bool:
    if isinstance(raw, bool):
        return raw
    if isinstance(raw, str):
        token = raw.strip().lower()
        if token in {"true", "yes", "present", "confirmed", "adequate", "pass", "consistent"}:
            return True
        if token in {"false", "no", "absent", "not_confirmed", "fail", "inconsistent"}:
            return False
    raise SchemaValidationError(f"{field_name} must be true or false, not {raw!r}")


def _parse_string(raw: Any, field_name: str) -> str:
    if not isinstance(raw, str) or not raw.strip():
        raise SchemaValidationError(f"{field_name} must be a non-empty string")
    return raw.strip()


def _parse_int(raw: Any, field_name: str) -> int:
    if isinstance(raw, bool) or not isinstance(raw, int) or raw < 0:
        raise SchemaValidationError(f"{field_name} must be a non-negative integer, not {raw!r}")
    return raw


def _parse_float(raw: Any, field_name: str) -> float:
    if isinstance(raw, bool) or not isinstance(raw, (int, float)):
        raise SchemaValidationError(f"{field_name} must be a number, not {raw!r}")
    result = float(raw)
    if not math.isfinite(result):
        raise SchemaValidationError(f"{field_name} must be a finite number, not {raw!r}")
    return result


def _parse_metric(raw: Any, field_name: str) -> Any:
    """Retain assay-specific metrics while rejecting non-finite numbers.

    The evaluator intentionally does not impose universal depth or allele
    balance cut-offs.  It still must not accept IEEE NaN/Infinity, which can
    otherwise bypass every comparison and make audit output non-JSON-safe.
    """

    def validate(value: Any, path: str) -> None:
        if isinstance(value, bool):
            return
        if isinstance(value, int):
            return
        if isinstance(value, float):
            if not math.isfinite(value):
                raise SchemaValidationError(f"{path} must not contain NaN or Infinity")
            return
        if isinstance(value, Mapping):
            for key, item in value.items():
                validate(item, f"{path}.{key}")
            return
        if isinstance(value, (list, tuple, set)):
            for index, item in enumerate(value):
                validate(item, f"{path}[{index}]")

    validate(raw, field_name)
    return raw


def _parse_genotype(raw: Any, field_name: str) -> GenotypeState:
    if isinstance(raw, bool):
        return GenotypeState.VARIANT_PRESENT if raw else GenotypeState.VARIANT_ABSENT
    aliases = {
        "PRESENT": GenotypeState.VARIANT_PRESENT,
        "POSITIVE": GenotypeState.VARIANT_PRESENT,
        "POS": GenotypeState.VARIANT_PRESENT,
        "HET": GenotypeState.VARIANT_PRESENT,
        "HETEROZYGOUS": GenotypeState.VARIANT_PRESENT,
        "HOMOZYGOUS_ALTERNATE": GenotypeState.VARIANT_PRESENT,
        "HOMOZYGOUS_ALT": GenotypeState.VARIANT_PRESENT,
        "0_1": GenotypeState.VARIANT_PRESENT,
        "1_0": GenotypeState.VARIANT_PRESENT,
        "1_1": GenotypeState.VARIANT_PRESENT,
        "0_1_": GenotypeState.VARIANT_PRESENT,
        "0_1_ALT": GenotypeState.VARIANT_PRESENT,
        "ABSENT": GenotypeState.VARIANT_ABSENT,
        "NEGATIVE": GenotypeState.VARIANT_ABSENT,
        "NEG": GenotypeState.VARIANT_ABSENT,
        "REFERENCE": GenotypeState.VARIANT_ABSENT,
        "HOMOZYGOUS_REFERENCE": GenotypeState.VARIANT_ABSENT,
        "HOM_REF": GenotypeState.VARIANT_ABSENT,
        "0_0": GenotypeState.VARIANT_ABSENT,
        "MOSAIC": GenotypeState.MOSAIC,
        "LOW_LEVEL_MOSAIC": GenotypeState.MOSAIC,
        "NO_CALL": GenotypeState.NO_CALL,
        "NOCALL": GenotypeState.NO_CALL,
        "MISSING": GenotypeState.NO_CALL,
        "DOT": GenotypeState.NO_CALL,
        "_": GenotypeState.NO_CALL,
    }
    # VCF phase separators are irrelevant to presence/absence.
    if isinstance(raw, str) and raw.strip() in {"./.", ".|."}:
        return GenotypeState.NO_CALL
    return _enum_value(raw, GenotypeState, field_name, aliases)  # type: ignore[return-value]


def _parse_adequacy(raw: Any, field_name: str) -> TestAdequacy:
    aliases = {
        "HIGH": TestAdequacy.ADEQUATE,
        "SUFFICIENT": TestAdequacy.ADEQUATE,
        "PASS": TestAdequacy.ADEQUATE,
        "VALIDATED": TestAdequacy.ADEQUATE,
        "GOOD": TestAdequacy.ADEQUATE,
        "LOW": TestAdequacy.INADEQUATE,
        "POOR": TestAdequacy.INADEQUATE,
        "INSUFFICIENT": TestAdequacy.INADEQUATE,
        "FAIL": TestAdequacy.INADEQUATE,
    }
    return _enum_value(raw, TestAdequacy, field_name, aliases)  # type: ignore[return-value]


def _parse_phenotype(raw: Any, field_name: str) -> PhenotypeCategory:
    aliases = {
        "HIGHLY_SPECIFIC_FOR_GENE": PhenotypeCategory.HIGHLY_SPECIFIC,
        "HIGHLY_CONSISTENT": PhenotypeCategory.HIGHLY_SPECIFIC,
        "SPECIFIC": PhenotypeCategory.HIGHLY_SPECIFIC,
        "HIGH": PhenotypeCategory.HIGHLY_SPECIFIC,
        "CONSISTENT_NOT_HIGHLY_SPECIFIC": PhenotypeCategory.CONSISTENT,
        "NOT_HIGHLY_SPECIFIC": PhenotypeCategory.CONSISTENT,
        "NOT_SPECIFIC": PhenotypeCategory.CONSISTENT,
        "GENETICALLY_HETEROGENEOUS": PhenotypeCategory.HETEROGENEOUS,
        "CONSISTENT_HIGH_GENETIC_HETEROGENEITY": PhenotypeCategory.HETEROGENEOUS,
        "INCOMPATIBLE": PhenotypeCategory.INCONSISTENT,
        "NOT_CONSISTENT": PhenotypeCategory.INCONSISTENT,
    }
    return _enum_value(raw, PhenotypeCategory, field_name, aliases)  # type: ignore[return-value]


def _parse_inheritance(raw: Any, field_name: str) -> InheritancePattern:
    aliases = {
        "AD": InheritancePattern.AUTOSOMAL_DOMINANT,
        "AUTOSOMAL_DOMINANT": InheritancePattern.AUTOSOMAL_DOMINANT,
        "AR": InheritancePattern.AUTOSOMAL_RECESSIVE,
        "AUTOSOMAL_RECESSIVE": InheritancePattern.AUTOSOMAL_RECESSIVE,
        "XL": InheritancePattern.X_LINKED,
        "X_LINKED": InheritancePattern.X_LINKED,
        "XLD": InheritancePattern.X_LINKED_DOMINANT,
        "X_LINKED_DOMINANT": InheritancePattern.X_LINKED_DOMINANT,
        "XLR": InheritancePattern.X_LINKED_RECESSIVE,
        "X_LINKED_RECESSIVE": InheritancePattern.X_LINKED_RECESSIVE,
        "MT": InheritancePattern.MITOCHONDRIAL,
        "MITO": InheritancePattern.MITOCHONDRIAL,
        "INCOMPATIBLE": InheritancePattern.OTHER,
        "UNSUPPORTED": InheritancePattern.OTHER,
    }
    return _enum_value(raw, InheritancePattern, field_name, aliases)  # type: ignore[return-value]


def _parse_validity(raw: Any, field_name: str) -> GeneDiseaseValidity:
    aliases = {
        "ESTABLISHED": GeneDiseaseValidity.DEFINITIVE,
        "VALID": GeneDiseaseValidity.DEFINITIVE,
    }
    return _enum_value(raw, GeneDiseaseValidity, field_name, aliases)  # type: ignore[return-value]


def _parse_mapping_list(raw: Any, field_name: str) -> tuple[Mapping[str, Any], ...]:
    values = raw if isinstance(raw, (list, tuple)) else (raw,)
    if not values:
        raise SchemaValidationError(f"{field_name} cannot be empty when supplied")
    result: list[Mapping[str, Any]] = []
    for value in values:
        if not isinstance(value, Mapping):
            raise SchemaValidationError(f"{field_name} must contain objects")
        result.append(dict(value))
    return tuple(result)


def _parse_provenance_list(raw: Any, field_name: str) -> tuple[Provenance, ...]:
    values = raw if isinstance(raw, (list, tuple)) else (raw,)
    if not values:
        raise SchemaValidationError(f"{field_name} cannot be empty when supplied")
    return tuple(Provenance.from_mapping(value) for value in values)


def _parse_evidence(
    raw: Any,
    parser: Callable[[Any, str], T],
    field_name: str,
    *,
    boolean: bool = False,
) -> EvidenceValue[T]:
    """Parse compact values and provenance envelopes without collapsing state."""

    if raw is _MISSING or raw is None:
        return EvidenceValue()
    if isinstance(raw, EvidenceValue):
        return raw

    state_raw: Any = "KNOWN"
    provenance_raw: Any = ()
    value_raw = raw
    # A Boolean envelope may use the task's ``TRUE``/``FALSE`` spelling in
    # either ``state`` or ``status``.  Keep that assertion even when the
    # other spelling is the generic ``KNOWN`` state, so an envelope such as
    # ``{"state": "KNOWN", "status": "FALSE", "value": True}`` cannot
    # quietly lose its negative assertion.
    asserted_boolean: Optional[bool] = None
    wrapped = isinstance(raw, Mapping) and any(
        key in raw for key in ("value", "state", "status", "provenance")
    )
    if wrapped:
        state_raw = raw.get("state", raw.get("status", "KNOWN"))
        provenance_raw = raw.get("provenance", ())
        value_raw = raw.get("value", _MISSING)
        if "state" in raw and "status" in raw:
            def state_signature(candidate: Any) -> tuple[str, Optional[bool]]:
                token = _token(candidate)
                if boolean and token in {"TRUE", "FALSE"}:
                    return EvidenceStatus.KNOWN.value, token == "TRUE"
                aliases = {
                    "UNAVAILABLE": EvidenceStatus.UNKNOWN.value,
                    "NOT_ASSESSED": EvidenceStatus.UNKNOWN.value,
                    "CONFLICT": EvidenceStatus.CONFLICTING.value,
                    "TRUE": EvidenceStatus.KNOWN.value,
                    "FALSE": EvidenceStatus.KNOWN.value,
                }
                if token in aliases:
                    return aliases[token], None
                return token, None

            state_signature_value = state_signature(raw["state"])
            status_signature_value = state_signature(raw["status"])
            if (
                state_signature_value[0] != status_signature_value[0]
                or (
                    state_signature_value[1] is not None
                    and status_signature_value[1] is not None
                    and state_signature_value[1] != status_signature_value[1]
                )
            ):
                raise SchemaValidationError(
                    f"{field_name}.state conflicts with {field_name}.status"
                )
            boolean_assertions = tuple(
                value
                for _, value in (state_signature_value, status_signature_value)
                if value is not None
            )
            if boolean_assertions:
                asserted_boolean = boolean_assertions[0]

    state_token = _token(state_raw)
    # The task's TRUE/FALSE evidence-state spelling is accepted for Booleans
    # while the serialized shared model remains KNOWN + a Boolean value.
    if boolean and state_token in {"TRUE", "FALSE"}:
        state = EvidenceStatus.KNOWN
        state_assertion = state_token == "TRUE"
        if asserted_boolean is not None and asserted_boolean != state_assertion:
            raise SchemaValidationError(
                f"{field_name}.state conflicts with {field_name}.status"
            )
        asserted_boolean = state_assertion
    else:
        aliases = {
            "TRUE": EvidenceStatus.KNOWN,
            "FALSE": EvidenceStatus.KNOWN,
            "UNAVAILABLE": EvidenceStatus.UNKNOWN,
            "NOT_ASSESSED": EvidenceStatus.UNKNOWN,
            "CONFLICT": EvidenceStatus.CONFLICTING,
        }
        state = _enum_value(state_raw, EvidenceStatus, f"{field_name}.state", aliases)  # type: ignore[assignment]

    if not wrapped and isinstance(raw, str):
        shorthand = raw.strip().lower()
        if shorthand in {"unknown", "unavailable", "not_assessed", "not assessed"}:
            state = EvidenceStatus.UNKNOWN
        elif shorthand in {"conflicting", "conflict"}:
            state = EvidenceStatus.CONFLICTING

    if boolean and asserted_boolean is not None and value_raw is _MISSING:
        value_raw = asserted_boolean

    if provenance_raw is None:
        provenance_raw = ()
    if isinstance(provenance_raw, Mapping):
        provenance_raw = (provenance_raw,)
    if not isinstance(provenance_raw, (list, tuple)):
        raise SchemaValidationError(f"{field_name}.provenance must be a list of objects")
    provenance = tuple(Provenance.from_mapping(item) for item in provenance_raw)

    if state != EvidenceStatus.KNOWN:
        return EvidenceValue(status=state, provenance=provenance)
    if value_raw is _MISSING or value_raw is None:
        raise SchemaValidationError(f"{field_name} is known but has no value")
    parsed = parser(value_raw, field_name)
    if asserted_boolean is not None and parsed != asserted_boolean:
        raise SchemaValidationError(
            f"{field_name} has {state_token} state but contradictory value {value_raw!r}"
        )
    return EvidenceValue(status=state, value=parsed, provenance=provenance)


@dataclass(frozen=True)
class DeNovoObservation:
    """Evidence for one proband/observation of the variant.

    Each field is independently provenance-bearing.  The evaluator only uses
    the subset that changes the selected pathway, so optional audit data never
    becomes an accidental feature-count requirement.
    """

    case_id: EvidenceValue[str] = field(default_factory=EvidenceValue.unknown)
    publication_case_id: EvidenceValue[str] = field(default_factory=EvidenceValue.unknown)
    case_independence: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)

    proband_variant_present: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    proband_genotype: EvidenceValue[GenotypeState] = field(default_factory=EvidenceValue.unknown)
    maternal_genotype: EvidenceValue[GenotypeState] = field(default_factory=EvidenceValue.unknown)
    paternal_genotype: EvidenceValue[GenotypeState] = field(default_factory=EvidenceValue.unknown)
    maternal_variant_absence_confidence: EvidenceValue[TestAdequacy] = field(default_factory=EvidenceValue.unknown)
    paternal_variant_absence_confidence: EvidenceValue[TestAdequacy] = field(default_factory=EvidenceValue.unknown)
    parental_testing_adequacy: EvidenceValue[TestAdequacy] = field(default_factory=EvidenceValue.unknown)
    parental_read_depth: EvidenceValue[Any] = field(default_factory=EvidenceValue.unknown)
    parental_allele_balance: EvidenceValue[Any] = field(default_factory=EvidenceValue.unknown)

    maternity_confirmed: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    paternity_confirmed: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    parentage_confirmation_method: EvidenceValue[str] = field(default_factory=EvidenceValue.unknown)
    parentage_confirmation_quality: EvidenceValue[TestAdequacy] = field(default_factory=EvidenceValue.unknown)
    parentage_disproven: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    family_relationship_structure: EvidenceValue[str] = field(default_factory=EvidenceValue.unknown)
    sample_identity_confirmation: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    independent_case_parentage_status: EvidenceValue[ParentageStatus] = field(
        default_factory=EvidenceValue.unknown
    )

    variant_call_quality: EvidenceValue[TestAdequacy] = field(default_factory=EvidenceValue.unknown)
    proband_read_support: EvidenceValue[Any] = field(default_factory=EvidenceValue.unknown)
    parental_mosaicism_evidence: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    proband_mosaicism_evidence: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    apparent_germline_mosaicism: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)

    inheritance_pattern: EvidenceValue[InheritancePattern] = field(default_factory=EvidenceValue.unknown)
    gene_disease_validity: EvidenceValue[GeneDiseaseValidity] = field(default_factory=EvidenceValue.unknown)
    phenotype_category: EvidenceValue[PhenotypeCategory] = field(default_factory=EvidenceValue.unknown)
    phenotype_consistency: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    phenotype_specificity: EvidenceValue[PhenotypeCategory] = field(default_factory=EvidenceValue.unknown)
    phenotype_heterogeneity: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    alternative_molecular_diagnosis: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    additional_pathogenic_variant_present: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)
    family_history_consistent: EvidenceValue[bool] = field(default_factory=EvidenceValue.unknown)

    gene_id: EvidenceValue[str] = field(default_factory=EvidenceValue.unknown)
    disease_id: EvidenceValue[str] = field(default_factory=EvidenceValue.unknown)

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "DeNovoObservation":
        if not isinstance(raw, Mapping):
            raise SchemaValidationError("each de novo observation must be an object")
        normalized: dict[str, Any] = {}
        for name, value in raw.items():
            target = _OBSERVATION_ALIASES.get(name, name)
            if target in normalized:
                raise SchemaValidationError(f"{name} and {target} both supply the same feature")
            normalized[target] = value

        # Parentage status is allowed as a consistency check, never as a
        # substitute for tested parents.  It is translated to the minimum
        # facts it can safely establish.
        parentage_hint_raw = normalized.get("independent_case_parentage_status", _MISSING)
        if parentage_hint_raw is not _MISSING:
            parentage_hint = _parse_evidence(
                parentage_hint_raw,
                _parse_parentage_status,
                "independent_case_parentage_status",
            )
            if (
                parentage_hint.status == EvidenceStatus.KNOWN
                and parentage_hint.value == ParentageStatus.CONFIRMED
            ):
                for field_name in ("maternity_confirmed", "paternity_confirmed"):
                    normalized.setdefault(field_name, True)
            elif parentage_hint.status in {EvidenceStatus.UNKNOWN} or (
                parentage_hint.status == EvidenceStatus.KNOWN
                and parentage_hint.value == ParentageStatus.ASSUMED
            ):
                # This is an explicit per-case assertion that full parental
                # relationship confirmation is unavailable.  Do not let a
                # shared top-level confirmed-parentage default silently turn
                # it into PS2.
                normalized.setdefault("maternity_confirmed", False)
                normalized.setdefault("paternity_confirmed", False)

        parsers: Mapping[str, tuple[Callable[[Any, str], Any], bool]] = {
            "case_id": (_parse_string, False),
            "publication_case_id": (_parse_string, False),
            "case_independence": (_parse_bool, True),
            "proband_variant_present": (_parse_bool, True),
            "proband_genotype": (_parse_genotype, False),
            "maternal_genotype": (_parse_genotype, False),
            "paternal_genotype": (_parse_genotype, False),
            "maternal_variant_absence_confidence": (_parse_adequacy, False),
            "paternal_variant_absence_confidence": (_parse_adequacy, False),
            "parental_testing_adequacy": (_parse_adequacy, False),
            "parental_read_depth": (_parse_metric, False),
            "parental_allele_balance": (_parse_metric, False),
            "maternity_confirmed": (_parse_bool, True),
            "paternity_confirmed": (_parse_bool, True),
            "parentage_confirmation_method": (_parse_string, False),
            "parentage_confirmation_quality": (_parse_adequacy, False),
            "parentage_disproven": (_parse_bool, True),
            "family_relationship_structure": (_parse_string, False),
            "sample_identity_confirmation": (_parse_bool, True),
            "independent_case_parentage_status": (_parse_parentage_status, False),
            "variant_call_quality": (_parse_adequacy, False),
            "proband_read_support": (_parse_metric, False),
            "parental_mosaicism_evidence": (_parse_bool, True),
            "proband_mosaicism_evidence": (_parse_bool, True),
            "apparent_germline_mosaicism": (_parse_bool, True),
            "inheritance_pattern": (_parse_inheritance, False),
            "gene_disease_validity": (_parse_validity, False),
            "phenotype_category": (_parse_phenotype, False),
            "phenotype_consistency": (_parse_bool, True),
            "phenotype_specificity": (_parse_phenotype, False),
            "phenotype_heterogeneity": (_parse_bool, True),
            "alternative_molecular_diagnosis": (_parse_bool, True),
            "additional_pathogenic_variant_present": (_parse_bool, True),
            "family_history_consistent": (_parse_bool, True),
            "gene_id": (_parse_string, False),
            "disease_id": (_parse_string, False),
        }
        unknown = sorted(set(normalized) - set(parsers))
        if unknown:
            raise SchemaValidationError("unsupported de novo feature(s): " + ", ".join(unknown))
        values = {
            name: _parse_evidence(normalized.get(name, _MISSING), parser, name, boolean=boolean)
            for name, (parser, boolean) in parsers.items()
        }
        if parentage_hint_raw is not _MISSING:
            hint = values["independent_case_parentage_status"]
            maternity = values["maternity_confirmed"]
            paternity = values["paternity_confirmed"]
            both_confirmed = known(maternity) and known(paternity)
            if hint.status == EvidenceStatus.KNOWN and hint.value == ParentageStatus.CONFIRMED and not both_confirmed:
                raise SchemaValidationError(
                    "independent_case_parentage_status=CONFIRMED conflicts with maternity/paternity evidence"
                )
            if hint.status == EvidenceStatus.KNOWN and hint.value == ParentageStatus.ASSUMED and both_confirmed:
                raise SchemaValidationError(
                    "independent_case_parentage_status=ASSUMED conflicts with confirmed parentage"
                )
        return cls(**values)

    def evidence(self, field_name: str) -> EvidenceValue[Any]:
        try:
            return getattr(self, field_name)
        except AttributeError as error:
            raise SchemaValidationError(f"unknown de novo feature {field_name!r}") from error

    def inherit_unknowns_from(self, shared: "DeNovoObservation") -> "DeNovoObservation":
        """Backward-compatible helper for all-unknown local observations.

        Prefer :meth:`inherit_omitted_from` in parsers: it retains the crucial
        distinction between an omitted field and an explicitly supplied
        ``UNKNOWN`` envelope.
        """

        replacements = {
            name: getattr(shared, name)
            for name in self.__dataclass_fields__
            if getattr(self, name).status == EvidenceStatus.UNKNOWN
            and getattr(shared, name).status != EvidenceStatus.UNKNOWN
        }
        return replace(self, **replacements)

    def inherit_omitted_from(
        self,
        shared: "DeNovoObservation",
        supplied_fields: Iterable[str],
    ) -> "DeNovoObservation":
        """Inherit shared context only for raw fields absent from this case."""

        supplied = set(supplied_fields)
        replacements = {
            name: getattr(shared, name)
            for name in self.__dataclass_fields__
            if name not in supplied
            and getattr(self, name).status == EvidenceStatus.UNKNOWN
            and getattr(shared, name).status != EvidenceStatus.UNKNOWN
        }
        return replace(self, **replacements)


@dataclass(frozen=True)
class DeNovoEvidence:
    """A validated de novo evidence set containing one or more observations."""

    observations: tuple[DeNovoObservation, ...]
    applicable_vcep_de_novo_rules: EvidenceValue[tuple[Mapping[str, Any], ...]] = field(
        default_factory=EvidenceValue.unknown
    )
    number_of_independent_de_novo_cases: EvidenceValue[int] = field(
        default_factory=EvidenceValue.unknown
    )
    evidence_provenance: tuple[Provenance, ...] = ()

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "DeNovoEvidence":
        if not isinstance(raw, Mapping):
            raise SchemaValidationError("biological_features must be a mapping")
        aliases = {"independent_cases": "observations", "cases": "observations", "gene_specific_de_novo_rules": "applicable_vcep_de_novo_rules"}
        normalized: dict[str, Any] = {}
        for name, value in raw.items():
            target = aliases.get(name)
            if target is None:
                observation_target = _OBSERVATION_ALIASES.get(name)
                target = observation_target or name
            if target in normalized:
                raise SchemaValidationError(f"{name} and {target} both supply the same feature")
            normalized[target] = value

        top_level = {
            "observations",
            "applicable_vcep_de_novo_rules",
            "number_of_independent_de_novo_cases",
            "evidence_provenance",
        }
        observation_keys = set(DeNovoObservation.__dataclass_fields__)
        accepted_observation_inputs = observation_keys | set(_OBSERVATION_ALIASES)
        unknown = sorted(set(normalized) - top_level - accepted_observation_inputs)
        if unknown:
            raise SchemaValidationError("unsupported de novo feature(s): " + ", ".join(unknown))

        shared_raw = {
            name: value
            for name, value in normalized.items()
            if name in observation_keys or name in _OBSERVATION_ALIASES
        }
        shared = DeNovoObservation.from_mapping(shared_raw)
        observations_raw = normalized.get("observations", _MISSING)
        if observations_raw is _MISSING:
            observations = (shared,)
        else:
            if not isinstance(observations_raw, (list, tuple)) or not observations_raw:
                raise SchemaValidationError("observations must be a non-empty list of objects")
            parsed_observations: list[DeNovoObservation] = []
            for item in observations_raw:
                if not isinstance(item, Mapping):
                    raise SchemaValidationError("each de novo observation must be an object")
                supplied: set[str] = set()
                for raw_name in item:
                    target = _OBSERVATION_ALIASES.get(raw_name, raw_name)
                    if target == "independent_case_parentage_status":
                        # The per-case assertion controls the translated
                        # parentage facts as well; never let shared confirmed
                        # parentage overwrite an explicit assumed/unknown
                        # local assertion.
                        supplied.update({
                            "independent_case_parentage_status",
                            "maternity_confirmed",
                            "paternity_confirmed",
                        })
                    elif target in DeNovoObservation.__dataclass_fields__:
                        supplied.add(target)
                parsed_observations.append(
                    DeNovoObservation.from_mapping(item).inherit_omitted_from(shared, supplied)
                )
            observations = tuple(parsed_observations)

        vcep = _parse_evidence(
            normalized.get("applicable_vcep_de_novo_rules", _MISSING),
            _parse_mapping_list,
            "applicable_vcep_de_novo_rules",
        )
        count = _parse_evidence(
            normalized.get("number_of_independent_de_novo_cases", _MISSING),
            _parse_int,
            "number_of_independent_de_novo_cases",
        )
        provenance_raw = normalized.get("evidence_provenance", ())
        if provenance_raw is None:
            provenance_raw = ()
        provenance = _parse_provenance_list(provenance_raw, "evidence_provenance") if provenance_raw else ()
        return cls(
            observations=observations,
            applicable_vcep_de_novo_rules=vcep,
            number_of_independent_de_novo_cases=count,
            evidence_provenance=provenance,
        )


@dataclass(frozen=True)
class ObservationDecision:
    """Auditable result for exactly one supplied patient/case observation."""

    index: int
    status: ObservationStatus
    parentage_status: Optional[ParentageStatus]
    phenotype_category: Optional[PhenotypeCategory]
    raw_points: float
    contributed_points: float
    reason: str
    case_id: Optional[str] = None
    duplicate_of: Optional[int] = None
    features_used: tuple[str, ...] = ()
    missing_critical_features: tuple[str, ...] = ()
    exclusions_evaluated: tuple[Mapping[str, str], ...] = ()
    guideline_locator: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "index": self.index,
            "case_id": self.case_id,
            "duplicate_of": self.duplicate_of,
            "status": self.status.value,
            "parentage_status": self.parentage_status.value if self.parentage_status else None,
            "phenotype_category": self.phenotype_category.value if self.phenotype_category else None,
            "raw_points": self.raw_points,
            "contributed_points": self.contributed_points,
            "features_used": list(self.features_used),
            "missing_critical_features": list(self.missing_critical_features),
            "exclusions_evaluated": [dict(item) for item in self.exclusions_evaluated],
            "guideline_locator": self.guideline_locator,
            "reason": self.reason,
        }


@dataclass(frozen=True)
class PS2PM6Decision:
    """One pooled, mutually exclusive PS2 *or* PM6 decision."""

    status: DecisionStatus
    assigned_criterion: Optional[str]
    strength: DeNovoStrength
    total_points: float
    confirmed_de_novo_cases: int
    assumed_de_novo_cases: int
    observations: tuple[ObservationDecision, ...]
    guideline: GuidelineReference
    decision_path: str
    reason: str
    features_used: tuple[str, ...] = ()
    missing_critical_features: tuple[str, ...] = ()
    exclusions_evaluated: tuple[Mapping[str, str], ...] = ()
    evidence_provenance: tuple[Provenance, ...] = ()

    def __post_init__(self) -> None:
        if self.status == DecisionStatus.ACTIVATED:
            if self.assigned_criterion not in {"PS2", "PM6"}:
                raise ValueError("an activated PS2/PM6 decision requires exactly one criterion")
            if self.strength == DeNovoStrength.NONE:
                raise ValueError("an activated PS2/PM6 decision requires a non-NONE strength")
        elif self.assigned_criterion is not None or self.strength != DeNovoStrength.NONE:
            raise ValueError("a non-activated PS2/PM6 decision cannot assign a criterion or strength")

    @property
    def criteria(self) -> tuple[str, str]:
        return ("PS2", "PM6")

    @property
    def criterion(self) -> Optional[str]:
        """Compatibility alias for consumers of the PVS1-shaped result API."""

        return self.assigned_criterion

    @property
    def features_missing(self) -> tuple[str, ...]:
        return self.missing_critical_features

    def to_dict(self) -> dict[str, Any]:
        missing = list(self.missing_critical_features)
        return {
            "criteria": list(self.criteria),
            "criterion": self.assigned_criterion,
            "assigned_criterion": self.assigned_criterion,
            "status": self.status.value,
            "strength": self.strength.value,
            "total_points": self.total_points,
            "confirmed_de_novo_cases": self.confirmed_de_novo_cases,
            "assumed_de_novo_cases": self.assumed_de_novo_cases,
            "observations": [item.to_dict() for item in self.observations],
            "guideline": self.guideline.to_dict(),
            "decision_path": self.decision_path,
            "features_used": list(self.features_used),
            "missing_critical_features": missing,
            "features_missing": missing,
            "exclusions_evaluated": [dict(item) for item in self.exclusions_evaluated],
            "evidence_provenance": [item.to_dict() for item in self.evidence_provenance],
            "reason": self.reason,
        }


def known(value: EvidenceValue[Any], expected: Any = True) -> bool:
    return value.status == EvidenceStatus.KNOWN and value.value == expected


def unique_provenance(values: Iterable[Provenance]) -> tuple[Provenance, ...]:
    result: list[Provenance] = []
    seen: set[tuple[Any, ...]] = set()
    for value in values:
        key = (value.source, value.kind.value, value.identifier, value.method, value.note)
        if key not in seen:
            seen.add(key)
            result.append(value)
    return tuple(result)


__all__ = [
    "DecisionStatus",
    "DeNovoEvidence",
    "DeNovoObservation",
    "DeNovoStrength",
    "EvidenceKind",
    "EvidenceStatus",
    "EvidenceValue",
    "GeneDiseaseValidity",
    "GenotypeState",
    "GuidelineReference",
    "InheritancePattern",
    "ObservationDecision",
    "ObservationStatus",
    "PS2PM6Decision",
    "PS2PM6Strength",
    "ParentageStatus",
    "PhenotypeCategory",
    "Provenance",
    "SchemaValidationError",
    "TestAdequacy",
    "known",
    "stable_unique",
]
