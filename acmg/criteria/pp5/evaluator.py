"""PP5 policy evaluation. No assertion is converted into criterion evidence."""

from collections.abc import Mapping

from .models import ExternalAssertion, PP5Result


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

# Verified CSpec example; only exact ID, version, gene and disease match applies.
_VERIFIED_VCEP = {
    ("GN039", "2.2"): {
        "id": "GN039",
        "version": "2.2",
        "name": "ClinGen RASopathy VCEP NRAS",
        "gene_id": "HGNC:7989",
        "disease_id": "MONDO:0021060",
        "pp5": "EXCLUDED",
        "source": "https://cspec.genome.network/cspec/ui/svi/doc/GN039?version=2.2",
    }
}


def _assertions(data: Mapping) -> list[ExternalAssertion]:
    raw = data.get("external_assertions")
    if raw is None:
        raw = [data] if any(key in data for key in _ASSERTION_FIELDS) else []
    result = []
    for item in raw:
        if isinstance(item, ExternalAssertion):
            result.append(item)
        elif isinstance(item, Mapping):
            values = {target: item[source] for source, target in _ASSERTION_FIELDS.items() if source in item}
            # A list may also contain records using the concise dataclass field names.
            values.update({key: item[key] for key in ExternalAssertion.__dataclass_fields__ if key in item})
            result.append(ExternalAssertion(**values))
        else:
            raise TypeError("external_assertions must contain mappings or ExternalAssertion records")
    return result


def evaluate_pp5(biological_features: Mapping | None = None) -> PP5Result:
    """Evaluate supplied evidence under current ClinGen policy by default.

    Historical selection is audit-only: source reputation and independence cannot
    be established from a classification label or ClinVar review status alone.
    """
    data = biological_features or {}
    if not isinstance(data, Mapping):
        raise TypeError("biological_features must be a mapping")
    policy = data.get("pp5_usage_policy", "CURRENT_CLINGEN")
    if policy not in ("CURRENT_CLINGEN", "HISTORICAL_ACMG_2015"):
        raise ValueError("pp5_usage_policy must be CURRENT_CLINGEN or HISTORICAL_ACMG_2015")
    name = data.get("guideline_name")
    expected_names = ({"ClinGen SVI", "ClinGen SVI PP5/BP6 recommendation"}
                      if policy == "CURRENT_CLINGEN" else {"ACMG/AMP", "ACMG/AMP 2015"})
    if name is not None and name not in expected_names:
        raise ValueError("guideline_name does not match the selected policy")
    version = data.get("guideline_version")
    if version is not None:
        raise ValueError("guideline_version does not match the selected verified policy")

    assertions = _assertions(data)
    variant = {key: data[key] for key in
               ("variant_id", "gene_id", "disease_id", "variant_normalization")
               if key in data}
    features_used = ["pp5_usage_policy"]
    conflicts = []
    vcep = None
    requested = data.get("applicable_vcep_specification")
    if requested is not None:
        if not isinstance(requested, Mapping):
            raise TypeError("applicable_vcep_specification must be a mapping")
        candidate = _VERIFIED_VCEP.get((requested.get("id"), requested.get("version")))
        if candidate and data.get("gene_id") in (candidate["gene_id"], "NRAS") and data.get("disease_id") == candidate["disease_id"]:
            vcep = dict(candidate)
            features_used.extend(["applicable_vcep_specification", "gene_id", "disease_id"])
        else:
            conflicts.append("The supplied VCEP specification is unknown, unversioned, or does not match the supplied gene and disease; its PP5 rule was not applied.")
    if requested and requested.get("pp5") == "ALLOWED":
        conflicts.append("The supplied VCEP claims PP5 is allowed, which conflicts with the ClinGen SVI recommendation; verify the specification and its source.")

    recommendations = []
    if any(a.evidence_available is True or a.evidence_independently_assessable is True for a in assertions):
        recommendations.append("Review underlying evidence independently under the appropriate ACMG/AMP criteria; do not count the classification itself or reuse evidence twice.")
    if any(a.conflicts is True or a.classification == "CONFLICTING" for a in assertions):
        recommendations.append("Resolve conflicting external assertions by reviewing their underlying evidence and provenance.")

    if policy == "CURRENT_CLINGEN":
        return PP5Result(
            criterion="PP5", status="NOT_APPLICABLE", strength="NONE",
            guideline={"name": "ClinGen SVI PP5/BP6 recommendation", "version": None,
                       "source": "https://pubmed.ncbi.nlm.nih.gov/29543229/"},
            reason="PP5 is not recommended under the selected ClinGen guidance.",
            variant=variant, external_assertions=assertions, features_used=features_used,
            recommendations=recommendations, applicable_vcep_specification=vcep,
            policy_conflicts=conflicts,
        )

    if vcep is not None:
        return PP5Result(
            criterion="PP5", status="NOT_APPLICABLE", strength="NONE",
            guideline={"name": "ACMG/AMP 2015", "version": None,
                       "source": "https://pubmed.ncbi.nlm.nih.gov/25741868/"},
            reason="The applicable, versioned VCEP specification explicitly excludes PP5.",
            variant=variant, external_assertions=assertions, features_used=features_used,
            recommendations=recommendations, applicable_vcep_specification=vcep,
            policy_conflicts=conflicts,
        )
    return PP5Result(
        criterion="PP5", status="UNDETERMINED", strength="NONE",
        guideline={"name": "ACMG/AMP 2015", "version": None,
                   "source": "https://pubmed.ncbi.nlm.nih.gov/25741868/"},
        reason="Historical PP5 was defined for a pathogenic report from a reputable source when its evidence was unavailable for independent evaluation. This outdated policy requires verified source-specific conditions and manual assessment; no PP5 evidence is assigned.",
        variant=variant, external_assertions=assertions, features_used=features_used,
        missing_critical_features=["verified_source_specific_conditions"],
        recommendations=recommendations, policy_conflicts=conflicts,
    )
