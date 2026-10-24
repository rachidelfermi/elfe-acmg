"""Offline BP6 policy selection; external assertions are context only."""

from collections.abc import Mapping

from .models import BP6Result, ExternalAssertion


_ASSERTION_FIELDS = {
    "external_classification": "classification",
    "classification_source": "source",
    "classification_date": "date",
    "classification_review_status": "review_status",
    "classification_evidence_available": "evidence_available",
    "classification_evidence_independently_assessable": "evidence_independently_assessable",
    "classification_conflicts": "conflicts",
    "classification_provenance": "provenance",
}

# CSpec GN080 v2.1, F9 / hemophilia B, explicitly marks BP6 Not Applicable.
_VERIFIED_VCEP = {
    ("GN080", "2.1"): {
        "id": "GN080",
        "version": "2.1",
        "name": "ClinGen Coagulation Factor Deficiency VCEP F9",
        "gene_id": "HGNC:3551",
        "disease_id": "MONDO:0010604",
        "bp6": "EXCLUDED",
        "source": "https://cspec.genome.network/cspec/ui/svi/doc/GN080?version=2.1",
    }
}


def _assertions(data: Mapping) -> list[ExternalAssertion]:
    raw = data.get("external_assertions")
    if raw is None:
        raw = [data] if any(key in data for key in _ASSERTION_FIELDS) else []
    if isinstance(raw, (str, bytes, Mapping)):
        raise TypeError("external_assertions must be a sequence of assertion records")
    result = []
    for item in raw:
        if isinstance(item, ExternalAssertion):
            result.append(item)
        elif isinstance(item, Mapping):
            values = {target: item[source] for source, target in _ASSERTION_FIELDS.items() if source in item}
            values.update({key: item[key] for key in ExternalAssertion.__dataclass_fields__ if key in item})
            result.append(ExternalAssertion(**values))
        else:
            raise TypeError("external_assertions must contain mappings or ExternalAssertion records")
    return result


def evaluate_bp6(biological_features: Mapping | None = None) -> BP6Result:
    """Evaluate BP6 using current ClinGen guidance unless history is explicit.

    Historical selection documents the 2015 rule but cannot verify the
    reputation, recency, independence or provenance of an external assertion.
    It therefore never awards BP6 automatically.
    """
    if biological_features is None:
        data = {}
    elif isinstance(biological_features, Mapping):
        data = biological_features
    else:
        raise TypeError("biological_features must be a mapping")

    policy = data.get("bp6_usage_policy", "CURRENT_CLINGEN")
    if policy not in ("CURRENT_CLINGEN", "HISTORICAL_ACMG_2015"):
        raise ValueError("bp6_usage_policy must be CURRENT_CLINGEN or HISTORICAL_ACMG_2015")
    name = data.get("guideline_name")
    expected_names = ({"ClinGen SVI", "ClinGen SVI PP5/BP6 recommendation"}
                      if policy == "CURRENT_CLINGEN" else {"ACMG/AMP", "ACMG/AMP 2015"})
    if name is not None and name not in expected_names:
        raise ValueError("guideline_name does not match the selected policy")
    if data.get("guideline_version") is not None:
        raise ValueError("Neither selected guideline has a verified numbered version")

    assertions = _assertions(data)
    variant = {key: data[key] for key in
               ("variant_id", "gene_id", "disease_id", "normalized_variant_representation")
               if key in data}
    features_used = ["bp6_usage_policy"]
    policy_conflicts = []
    vcep = None
    requested = data.get("applicable_vcep_specification")
    if requested is not None:
        if not isinstance(requested, Mapping):
            raise TypeError("applicable_vcep_specification must be a mapping")
        candidate = _VERIFIED_VCEP.get((requested.get("id"), requested.get("version")))
        if (candidate is not None and data.get("gene_id") in (candidate["gene_id"], "F9")
                and data.get("disease_id") == candidate["disease_id"]):
            vcep = dict(candidate)
            features_used.extend(["applicable_vcep_specification", "gene_id", "disease_id"])
        else:
            policy_conflicts.append(
                "The supplied VCEP specification is unknown, unversioned, or does not match the supplied gene and disease; its BP6 rule was not applied."
            )
        if requested.get("bp6") == "ALLOWED" and policy == "CURRENT_CLINGEN":
            policy_conflicts.append(
                "The supplied VCEP claims BP6 is allowed, which conflicts with the ClinGen SVI recommendation; verify the specification and its source."
            )

    recommendations = []
    if any(a.classification in ("BENIGN", "LIKELY_BENIGN") for a in assertions):
        recommendations.append(
            "Treat external benign classifications as context; assess their underlying primary evidence under appropriate ACMG benign criteria without counting it twice."
        )
    elif any(a.evidence_independently_assessable is True for a in assertions):
        recommendations.append(
            "Assess independently available primary evidence under appropriate ACMG criteria without counting it twice."
        )
    if any(a.conflicts is True or a.classification == "CONFLICTING" for a in assertions):
        recommendations.append("Review conflicting assertions, their provenance, and underlying evidence.")

    if policy == "CURRENT_CLINGEN":
        return BP6Result(
            criterion="BP6", status="NOT_APPLICABLE", strength="NONE",
            guideline={"name": "ClinGen SVI PP5/BP6 recommendation", "version": None,
                       "source": "https://pubmed.ncbi.nlm.nih.gov/29543229/"},
            reason="BP6 is not recommended under the selected ClinGen policy.",
            variant=variant, external_assertions=assertions, features_used=features_used,
            recommendations=recommendations, applicable_vcep_specification=vcep,
            policy_conflicts=policy_conflicts,
        )

    if vcep is not None:
        status = "NOT_APPLICABLE"
        reason = "The applicable, versioned VCEP specification explicitly excludes BP6."
        missing = []
    else:
        status = "UNDETERMINED"
        reason = (
            "Historical ACMG/AMP 2015 BP6 described a recent benign report from a reputable "
            "source whose evidence was unavailable for independent evaluation. Source reputation, "
            "recency and evidence independence require manual verification; no BP6 evidence is assigned."
        )
        missing = ["verified_historical_source_conditions"]
    return BP6Result(
        criterion="BP6", status=status, strength="NONE",
        guideline={"name": "ACMG/AMP 2015", "version": None,
                   "source": "https://pubmed.ncbi.nlm.nih.gov/25741868/"},
        reason=reason, variant=variant, external_assertions=assertions,
        features_used=features_used, missing_critical_features=missing,
        recommendations=recommendations, applicable_vcep_specification=vcep,
        policy_conflicts=policy_conflicts,
    )
