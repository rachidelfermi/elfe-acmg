"""Deterministic BS1 decision using only complete, sourced rule pathways."""

from .guidelines import matching_rule
from .models import BS1Features, BS1Result, EvidenceState, Status, Strength


def evaluate_bs1(features: BS1Features) -> BS1Result:
    """Evaluate BS1; missing data never becomes a zero frequency or a cutoff."""
    if not isinstance(features, BS1Features):
        raise TypeError("evaluate_bs1 expects BS1Features")
    rule = features.guideline or matching_rule(features.gene, features.disease, features.inheritance)
    path = "VCEP_SPECIFIC_FREQUENCY_THRESHOLD" if rule and rule.rule_id in {
        "GN100:BS1", "1529230095:BS1"} else "DISEASE_SPECIFIC_FREQUENCY_THRESHOLD"
    citation = ({"name": rule.name, "version": rule.version, "rule_id": rule.rule_id,
                 "source_url": rule.source_url} if rule else None)

    def result(status: Status, reason: str, *, strength: Strength = Strength.NONE,
               frequency: float | None = None, missing: tuple[str, ...] = (),
               used: tuple[str, ...] = (), exclusions: tuple[str, ...] = ()) -> BS1Result:
        return BS1Result("BS1", status, strength, path if rule else None, frequency,
                         rule.threshold if rule else None, used, missing, exclusions,
                         citation, reason)

    if features.ba1_established:
        return result(Status.NOT_APPLICABLE, "BA1 is established; BS1 is not counted separately.",
                      exclusions=("BA1 established",))
    if rule is None:
        return result(Status.UNDETERMINED, "No verified BS1 threshold for this gene, disease and inheritance.",
                      missing=("applicable verified BS1 threshold",))
    if not rule.verified:
        return result(Status.UNDETERMINED, "The supplied guideline rule is not attested as verified.",
                      missing=("verified guideline rule",))
    if (features.gene, features.disease, features.inheritance) != (rule.gene, rule.disease, rule.inheritance):
        return result(Status.NOT_APPLICABLE, "Guideline gene, disease or inheritance scope does not match.",
                      exclusions=("guideline scope mismatch",))
    if features.variant.variant_type in rule.excluded_variant_types:
        return result(Status.NOT_APPLICABLE, "This variant class is outside the guideline's BS1 scope.",
                      exclusions=("variant class",))
    if features.variant_mapping != EvidenceState.KNOWN or features.gene_disease_validity != EvidenceState.KNOWN:
        return result(Status.UNDETERMINED, "Variant mapping or gene-disease applicability is unresolved.",
                      missing=("variant mapping and gene-disease validity",))
    if not features.variant.normalized and features.variant.variant_type == "INDEL":
        return result(Status.UNDETERMINED, "Indel normalization is unverified.",
                      missing=("normalized indel identity",))
    if features.exception_status == EvidenceState.KNOWN:
        if not features.exception_source:
            return result(Status.UNDETERMINED, "A reported exception has no supporting source.",
                          missing=("exception source",))
        if not features.exception_disease_scope:
            return result(Status.UNDETERMINED, "A reported exception has no disease scope.",
                          missing=("exception disease scope",))
        if features.exception_disease_scope != features.disease:
            return result(Status.UNDETERMINED, "The reported exception concerns another disease; review applicability.",
                          missing=("applicable exception review",))
        return result(Status.NOT_APPLICABLE, "A sourced, variant-specific BS1 exception applies.",
                      exclusions=("known BS1 exception",))
    if features.exception_status == EvidenceState.CONFLICTING:
        return result(Status.UNDETERMINED, "Conflicting evidence about a BS1 exception.",
                      missing=("resolved exception status",))
    if features.founder_variant or features.risk_allele or features.hypomorphic_allele:
        if features.exception_status != EvidenceState.NOT_APPLICABLE:
            return result(Status.UNDETERMINED, "Special variant context needs an explicit exception review.",
                          missing=("reviewed variant-specific exception status",))
    if features.conflicting_evidence:
        return result(Status.UNDETERMINED, "Conflicting evidence prevents BS1 activation.",
                      missing=("resolved evidence conflict",))

    eligible = []
    issues = []
    for record in features.frequencies:
        if record.variant != features.variant:
            issues.append("exact variant, assembly or normalization mismatch")
            continue
        if record.metric != rule.metric:
            issues.append("required frequency metric")
            continue
        if rule.database and record.database.lower() != rule.database.lower():
            issues.append("required population database")
            continue
        if record.quality != EvidenceState.KNOWN or record.filter_status != "PASS" or record.technical_limitation:
            issues.append("validated population record quality")
            continue
        if record.frequency is None:
            issues.append("population frequency")
            continue
        if rule.minimum_allele_number is not None and (
            record.allele_number is None or record.allele_number < rule.minimum_allele_number
        ):
            issues.append("minimum allele number")
            continue
        if rule.minimum_exome_depth is not None:
            if record.platform not in {"exome", "genome"}:
                issues.append("sequencing platform for depth rule")
                continue
            if record.platform == "exome" and (record.depth is None or record.depth < rule.minimum_exome_depth):
                issues.append("minimum exome read depth")
                continue
        if rule.eligible_ancestries is not None and record.ancestry not in rule.eligible_ancestries:
            issues.append("eligible ancestry")
            continue
        if record.metric.startswith("ANCESTRY_") and (
            not record.ancestry or record.sampling != EvidenceState.KNOWN
        ):
            issues.append("reliable ancestry sampling")
            continue
        if rule.requires_non_founder and record.founder_population:
            issues.append("non-founder population evidence")
            continue
        eligible.append(record)

    if not eligible:
        return result(Status.UNDETERMINED, "No complete, reliable record has the guideline's required metric.",
                      missing=tuple(dict.fromkeys(issues)) or ("population frequency",))

    above = [r for r in eligible if (r.frequency > rule.threshold if rule.operator == ">"
                                    else r.frequency >= rule.threshold)]
    below = [r for r in eligible if r not in above]
    if above and below:
        return result(Status.UNDETERMINED, "Reliable population records straddle the BS1 threshold; reconcile sources.",
                      missing=("resolved population source conflict",))
    selected = max(eligible, key=lambda r: r.frequency)  # deterministic if input order is fixed
    used = ("exact variant identity", "population record quality", rule.metric,
            "gene/disease/inheritance scope", "guideline threshold")
    if above:
        if features.pm2_established:
            return result(Status.UNDETERMINED, "PM2 and BS1 conflict on population frequency assumptions.",
                          frequency=selected.frequency, missing=("resolved PM2/BS1 conflict",), used=used)
        return result(Status.ACTIVATED,
                      f"{rule.metric} {selected.frequency:g} {rule.operator} {rule.threshold:g} under {rule.rule_id}.",
                      strength=rule.strength, frequency=selected.frequency, used=used,
                      exclusions=("BA1 not established", "no known BS1 exception", "record passes quality checks"))
    return result(Status.NOT_MET,
                  f"Reliable {rule.metric} {selected.frequency:g} does not satisfy {rule.operator} {rule.threshold:g}.",
                  frequency=selected.frequency, used=used)
