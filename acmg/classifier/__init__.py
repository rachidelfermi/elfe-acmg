"""Unified evaluator orchestration with conservative classification."""

from __future__ import annotations

from dataclasses import fields, is_dataclass
from enum import Enum
from importlib import import_module
from types import UnionType
from typing import Any, Union, get_args, get_origin, get_type_hints

from ..classification_rules import STRENGTHS, classify_criteria
from ..criterion_registry import ALL_CRITERIA, DEPRECATED, EVALUATORS, evaluator_for
from ..evidence_dependencies import resolve_dependencies


def _plain(value):
    if isinstance(value, Enum):
        return value.value
    if is_dataclass(value):
        return {f.name: _plain(getattr(value, f.name)) for f in fields(value)}
    if isinstance(value, dict):
        return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    return value


def _construct(hint, value):
    if value is None:
        return None
    origin, args = get_origin(hint), get_args(hint)
    if origin in (Union, UnionType):
        choices = [arg for arg in args if arg is not type(None)]
        if isinstance(value, dict):
            for choice in choices:
                candidate = get_origin(choice) or choice
                if isinstance(candidate, type) and is_dataclass(candidate):
                    return _construct(choice, value)
        return _construct(choices[0], value)
    if origin in (tuple, list) and isinstance(value, (tuple, list)):
        converted = [_construct(args[0], v) for v in value]
        return tuple(converted) if origin is tuple else converted
    if isinstance(origin, type) and is_dataclass(origin):
        hint = origin
    elif origin is not None and args:
        return _construct(next(a for a in args if a is not type(None)), value)
    if isinstance(hint, type) and issubclass(hint, Enum):
        return hint(value)
    if isinstance(hint, type) and is_dataclass(hint) and isinstance(value, dict):
        names = {f.name for f in fields(hint)}
        hints = get_type_hints(hint)
        return hint(**{k: _construct(hints[k], v) for k, v in value.items() if k in names})
    return value


_TYPED = {
    "BA1": ("acmg.criteria.ba1.models", "BA1Evidence"),
    "BS1": ("acmg.criteria.bs1.models", "BS1Features"),
    "BP2": ("acmg.criteria.bp2.models", "BiologicalFeatures"),
    "BP3": ("acmg.criteria.bp3.models", "BiologicalFeatures"),
    "BP4": ("acmg.criteria.bp4.models", "BiologicalFeatures"),
    "PP4": ("acmg.criteria.pp4.models", "PP4Features"),
}


def _evaluate(code, payload):
    if code in _TYPED:
        module, name = _TYPED[code]
        payload = _construct(getattr(import_module(module), name), payload)
    decision = evaluator_for(code)(payload)
    row = _plain(decision.to_dict() if hasattr(decision, "to_dict") else decision)
    if not isinstance(row, dict):
        raise ValueError("Evaluator returned no structured decision")
    row["criterion"] = (row.get("assigned_criterion") or "PS2_PM6") if code == "PS2_PM6" else code
    row["status"] = str(row.get("status", "UNDETERMINED"))
    row["strength"] = str(row.get("strength") or "NONE")
    if row["status"] not in {"ACTIVATED", "NOT_MET", "UNDETERMINED", "NOT_APPLICABLE"}:
        raise ValueError(f"Invalid status: {row['status']}")
    if row["strength"] not in STRENGTHS:
        raise ValueError(f"Invalid strength: {row['strength']}")
    if row["status"] == "ACTIVATED" and row["strength"] == "NONE":
        raise ValueError("Activated criterion has no strength")
    if row["status"] != "ACTIVATED":
        row["strength"] = "NONE"
    row.setdefault("features_used", [])
    row.setdefault("missing_critical_features", row.get("features_missing", []))
    row.setdefault("decision_path", row.get("path"))
    row.setdefault("exclusions_evaluated", row.get("exclusions", []))
    row.setdefault("evidence_provenance", row.get("provenance", []))
    row.setdefault("guideline", None)
    row.setdefault("reason", "")
    row["counted"] = row["status"] == "ACTIVATED"
    return row


def classify_variant(biological_features: dict[str, Any]) -> dict[str, Any]:
    """Accept ``{"variant_id": ..., "criteria": {"PVS1": {...}, ...}}``."""
    if not isinstance(biological_features, dict):
        raise TypeError("biological_features must be a mapping")
    inputs = biological_features.get("criteria", {})
    if not isinstance(inputs, dict):
        raise TypeError("criteria must be a mapping")
    unknown = set(inputs) - ALL_CRITERIA
    if unknown:
        raise ValueError(f"Unknown criteria: {sorted(unknown)}")
    if "PS2_PM6" in inputs and ("PS2" in inputs or "PM6" in inputs):
        raise ValueError("Supply de novo evidence once under PS2_PM6")
    if "PS2" in inputs and "PM6" in inputs and inputs["PS2"] != inputs["PM6"]:
        raise ValueError("PS2 and PM6 cannot use separate assessments")
    if biological_features.get("classification_method", "ACMG_AMP_2015") != "ACMG_AMP_2015":
        raise ValueError("Only ACMG_AMP_2015 combination rules are implemented")

    rows, problems = [], []
    # Criterion payloads are independent; reject an explicit cross-criterion
    # identity mismatch before their evidence can be combined.
    identities = []
    for code, payload in inputs.items():
        if code in DEPRECATED or not isinstance(payload, dict):
            continue
        direct = payload.get("variant_id", payload.get("candidate_variant_id"))
        nested = payload.get("variant")
        if direct:
            identities.append((code, direct))
        if isinstance(nested, dict) and nested.get("variant_id"):
            identities.append((code, nested["variant_id"]))
    if biological_features.get("variant_id"):
        identities.append(("variant", biological_features["variant_id"]))
    if len({value for _, value in identities}) > 1:
        problems.append({"criterion": "variant_identity", "reason":
                         f"Conflicting variant IDs in {', '.join(code for code, _ in identities)}"})
    for field in ("gene_id", "disease_id"):
        contexts = [(code, payload[field]) for code, payload in inputs.items()
                    if code not in DEPRECATED and isinstance(payload, dict)
                    and isinstance(payload.get(field), str) and payload[field]]
        if isinstance(biological_features.get(field), str) and biological_features[field]:
            contexts.append(("variant", biological_features[field]))
        if len({value for _, value in contexts}) > 1:
            problems.append({"criterion": field, "reason":
                             f"Conflicting {field} values in {', '.join(code for code, _ in contexts)}"})
    for code in EVALUATORS:
        payload = inputs.get(code)
        if code == "PS2_PM6":
            payload = inputs.get(code, inputs.get("PS2", inputs.get("PM6")))
        if payload is None:
            continue
        try:
            rows.append(_evaluate(code, payload))
        except (ImportError, AttributeError, TypeError, ValueError, KeyError) as error:
            problems.append({"criterion": code, "reason": f"Evaluator unavailable or input unsupported: {error}"})
    for code in sorted(DEPRECATED & set(inputs)):
        rows.append({"criterion": code, "status": "NOT_APPLICABLE", "strength": "NONE",
                     "decision_path": "DISCONTINUED", "features_used": [],
                     "missing_critical_features": [], "exclusions_evaluated": ["deprecated_criterion"],
                     "evidence_provenance": [],
                     "guideline": {"name": "ClinGen SVI reputable-source recommendation", "version": "2018"},
                     "reason": "Deprecated reputable-source criterion is excluded.", "counted": False})
    counted, conflicts, resolved = resolve_dependencies(rows, inputs)
    unsupported = [r["criterion"] for r in counted if r["strength"] == "VERY_STRONG" and r["criterion"] != "PVS1"]
    if unsupported:
        problems.append({"criterion": ", ".join(unsupported), "reason": "No verified Table 5 rule for non-PVS1 Very Strong evidence."})
    if biological_features.get("vcep_specification"):
        problems.append({"criterion": "VCEP", "reason": "Complete VCEP scope, prohibitions and combination rules are not registered."})
    classification, reason = ("VUS", "Required guideline or evaluator is unsupported.") if problems else classify_criteria(counted)
    # The SVI PM2 adjustment permits PVS1 + PM2_Supporting, specifically.
    # Do not generalize it to an arbitrary single supporting criterion.
    if classification == "LIKELY_PATHOGENIC" and len(counted) == 2 and not any(
        r["criterion"] == "PM2" and r["strength"] == "SUPPORTING" for r in counted
    ) and any(r["criterion"] == "PVS1" and r["strength"] == "VERY_STRONG" for r in counted):
        other = next(r for r in counted if r["criterion"] != "PVS1")
        if other["strength"] == "SUPPORTING":
            classification, reason = "VUS", "A single non-PM2 supporting criterion does not meet a verified combination."
    if conflicts:
        classification, reason = "VUS", "Conflicting evidence requires review."
    suppressed = {r["criterion"] for r in rows if "suppression_reason" in r}
    if any(d.get("unresolved", True) and set(d["criteria"]) & suppressed for d in resolved):
        classification, reason = "VUS", "Unresolved evidence dependency requires review."
    by_status = lambda status: [r for r in rows if r["status"] == status]
    supplied = {r["criterion"] for r in rows}
    return {
        "variant_id": biological_features.get("variant_id"),
        "classification": classification, "classification_method": "ACMG_AMP_2015",
        "criteria_activated": by_status("ACTIVATED"),
        "criteria_not_met": by_status("NOT_MET"),
        "criteria_undetermined": by_status("UNDETERMINED") + [
            {"criterion": c, "status": "UNDETERMINED", "reason": "No evidence supplied"}
            for c in sorted(set(EVALUATORS) - {"PS2_PM6"} - supplied) if c not in inputs
        ] + problems,
        "criteria_not_applicable": by_status("NOT_APPLICABLE"),
        "evidence_conflicts": conflicts, "evidence_dependencies_resolved": resolved,
        "guideline": {"name": "ACMG/AMP Table 5 with ClinGen SVI PM2 combination update",
                      "version": "2015 + SVI PM2 v1.0 (2020)",
                      "sources": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/",
                                  "https://www.clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf"]},
        "reason": reason,
    }
