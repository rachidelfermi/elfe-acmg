"""Published PS3 strength rules and versioned VCEP configuration parsing."""

from __future__ import annotations

from typing import Any, Mapping

from .models import (
    DecisionStatus,
    GuidelineReference,
    SchemaValidationError,
    Strength,
    VCEPRule,
    VCEPSpecification,
    token,
)


ACMG_2015 = GuidelineReference(
    name="ACMG/AMP sequence variant interpretation guidelines",
    version="2015",
    rule_id="PS3",
    source_url="https://pubmed.ncbi.nlm.nih.gov/25741868/",
)

CLINGEN_SVI_2019 = GuidelineReference(
    name="ClinGen SVI functional evidence recommendations",
    version="2019",
    rule_id="Brnich2019-PS3/BS3",
    source_url="https://pubmed.ncbi.nlm.nih.gov/31892348/",
)

# Brnich et al. 2019, Table 3. Comparisons intentionally preserve the
# publication's strict greater-than symbols. No threshold is inferred from
# sensitivity, specificity, effect size, or statistical significance alone.
def strength_from_oddspath(value: float) -> Strength:
    if value > 350:
        return Strength.VERY_STRONG
    if value > 18.7:
        return Strength.STRONG
    if value > 4.3:
        return Strength.MODERATE
    if value > 2.1:
        return Strength.SUPPORTING
    return Strength.NONE


def benign_strength_from_oddspath(value: float) -> Strength:
    """Return the Table 3 BS3-equivalent tier for audit only.

    The PS3 evaluator never activates the returned benign criterion.  Strong
    is the maximum benign tier recommended by Brnich et al.
    """
    if value < 0.053:
        return Strength.STRONG
    if value < 0.23:
        return Strength.MODERATE
    if value < 0.48:
        return Strength.SUPPORTING
    return Strength.NONE


def parse_strength(raw: Any, field: str = "strength") -> Strength:
    try:
        return raw if isinstance(raw, Strength) else Strength(token(raw))
    except ValueError as exc:
        raise SchemaValidationError(f"{field} is invalid: {raw!r}") from exc


def parse_status(raw: Any, field: str = "status") -> DecisionStatus:
    try:
        return raw if isinstance(raw, DecisionStatus) else DecisionStatus(token(raw))
    except ValueError as exc:
        raise SchemaValidationError(f"{field} is invalid: {raw!r}") from exc


def parse_vcep_specification(raw: Mapping[str, Any]) -> VCEPSpecification:
    if not isinstance(raw, Mapping):
        raise SchemaValidationError("vcep_ps3_specification entries must be objects")
    required = ("name", "version", "source_url", "gene_id")
    missing = [name for name in required if not raw.get(name)]
    specification_id = raw.get("specification_id", raw.get("rule_id"))
    if not specification_id:
        missing.append("specification_id")
    if missing:
        raise SchemaValidationError(
            "VCEP specifications require pinned provenance: " + ", ".join(missing)
        )
    rules_raw = raw.get("rules")
    if not isinstance(rules_raw, (list, tuple)) or not rules_raw:
        raise SchemaValidationError("VCEP specification.rules must be a non-empty list")
    rules: list[VCEPRule] = []
    for index, item in enumerate(rules_raw):
        if not isinstance(item, Mapping):
            raise SchemaValidationError(f"VCEP rule {index} must be an object")
        decision = item.get("decision", item)
        if not isinstance(decision, Mapping):
            raise SchemaValidationError(f"VCEP rule {index}.decision must be an object")
        rule_id = str(item.get("rule_id", item.get("source_locator", f"rule-{index + 1}")))
        when = item.get("when", {})
        if not isinstance(when, Mapping):
            raise SchemaValidationError(f"VCEP rule {rule_id}.when must be an object")
        status = parse_status(decision.get("status", "ACTIVATED"), f"VCEP {rule_id}.status")
        strength = parse_strength(decision.get("strength", "NONE"), f"VCEP {rule_id}.strength")
        if status == DecisionStatus.ACTIVATED and strength == Strength.NONE:
            raise SchemaValidationError("activated VCEP rules require a non-NONE strength")
        if status != DecisionStatus.ACTIVATED and strength != Strength.NONE:
            raise SchemaValidationError("non-activated VCEP rules must use strength NONE")
        rules.append(VCEPRule(
            rule_id=rule_id,
            when=dict(when),
            status=status,
            strength=strength,
            reason=str(decision.get("reason", item.get("reason", "The pinned VCEP rule matched."))),
            decision_path=str(decision.get("decision_path", f"VCEP.{rule_id}")),
        ))
    return VCEPSpecification(
        name=str(raw["name"]),
        version=str(raw["version"]),
        specification_id=str(specification_id),
        source_url=str(raw["source_url"]),
        gene_id=str(raw["gene_id"]),
        disease_id=str(raw["disease_id"]) if raw.get("disease_id") else None,
        rules=tuple(rules),
        exclusive=bool(raw.get("exclusive", True)),
    )


def vcep_reference(spec: VCEPSpecification, rule: VCEPRule) -> GuidelineReference:
    return GuidelineReference(
        name=spec.name,
        version=spec.version,
        rule_id=f"{spec.specification_id}:{rule.rule_id}",
        source_url=spec.source_url,
    )
