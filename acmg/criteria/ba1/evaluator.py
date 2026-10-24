"""Conservative, deterministic BA1 evaluation from supplied population evidence."""

from .guidelines import GENERAL, SVI_EXCEPTION_ALLELES, select_guideline
from .models import BA1Evidence, BA1Result, EvidenceState, PopulationRecord


def _result(status, reason, rule=None, record=None, frequency=None, features=(), exclusions=()):
    return BA1Result(
        "BA1", status, "STAND_ALONE" if status == "ACTIVATED" else None,
        ("GENERAL_POPULATION_FREQUENCY" if rule == GENERAL else "VCEP_POPULATION_FREQUENCY")
        if rule else None,
        frequency, rule.threshold if rule else None, tuple(features), tuple(exclusions),
        {"name": rule.name, "version": rule.version, "rule_id": rule.rule_id}
        if rule else None, reason,
    )


def _frequency(record: PopulationRecord, metric: str) -> tuple[float | None, str | None]:
    if record.allele_count is not None and record.allele_number is not None:
        if (record.allele_number <= 0 or record.allele_count < 0
                or record.allele_count > record.allele_number):
            return None, "Invalid allele count or allele number"
    elif record.allele_count is not None or record.allele_number is not None:
        return None, "Allele count and allele number must be supplied together"
    value = getattr(record, metric)
    if value is not None and not 0 <= value <= 1:
        return None, "Frequency outside [0, 1]"
    if metric == "allele_frequency" and record.allele_count is not None:
        calculated = record.allele_count / record.allele_number
        if value is not None and abs(calculated - value) > max(0.00001, 1 / record.allele_number):
            return None, "Allele frequency conflicts with allele count/number"
        value = calculated
    return value, None


def _identity_ok(variant, record):
    if not all((variant.assembly, variant.chromosome, variant.position,
                variant.reference, variant.alternate)):
        return False
    return (variant.assembly == record.assembly
            and variant.chromosome.removeprefix("chr") == record.chromosome.removeprefix("chr")
            and variant.position == record.position
            and variant.reference.upper() == record.reference.upper()
            and variant.alternate.upper() == record.alternate.upper()
            and variant.normalized is True and record.normalized is True)


def evaluate_ba1(evidence: BA1Evidence) -> BA1Result:
    """Evaluate only BA1; inputs are validated evidence, never raw annotations."""
    if not isinstance(evidence, BA1Evidence):
        raise TypeError("evaluate_ba1 expects BA1Evidence")
    variant = evidence.variant
    rule = select_guideline(evidence.gene_id, evidence.disease_id)
    if rule is None:
        return _result("UNDETERMINED", "Disease-specific BA1 rule cannot be selected without the matching disease scope")
    if rule == GENERAL:
        if evidence.guideline_modification_check != EvidenceState.KNOWN or evidence.known_modification is None:
            return _result("UNDETERMINED", "Gene/disease BA1 modification status is unknown", rule)
        if evidence.known_modification:
            return _result("UNDETERMINED", "A BA1 modification exists but no verified bundled specification matches", rule)
    elif (evidence.guideline_modification_check == EvidenceState.KNOWN
          and evidence.known_modification is False):
        return _result("UNDETERMINED", "Selected VCEP rule conflicts with modification check", rule)
    if variant.variant_type not in rule.allowed_variant_types:
        return _result("NOT_APPLICABLE", "No verified BA1 pathway for this variant class", rule)
    chromosome = (variant.chromosome or "").removeprefix("chr").upper()
    if chromosome in ("M", "MT"):
        return _result("NOT_APPLICABLE", "Mitochondrial frequency requires a separate verified rule", rule)
    if chromosome in ("X", "Y") and variant.pseudoautosomal_region is not True:
        return _result("UNDETERMINED", "Sex chromosome population frequency requires a validated ploidy-specific pathway", rule)
    if variant.mapping_confidence != EvidenceState.KNOWN:
        return _result("UNDETERMINED", "Exact variant mapping is unverified", rule)
    if variant.normalized is not True:
        return _result("UNDETERMINED", "Variant normalization is unverified", rule)
    if evidence.frequency_conflict:
        return _result("UNDETERMINED", "Population frequency sources conflict", rule)
    if evidence.pathogenic_evidence_conflict:
        return _result("UNDETERMINED", "Pathogenic evidence conflicts with BA1; review required", rule)
    if (evidence.reduced_penetrance or evidence.late_onset or evidence.risk_allele
            or evidence.hypomorphic_allele or evidence.founder_pathogenic_context):
        return _result("UNDETERMINED", "Special disease or allele context requires variant-specific review", rule)

    checked = ["variant_class", "variant_mapping", "guideline_scope"]
    if variant.clingen_allele_id in SVI_EXCEPTION_ALLELES and rule == GENERAL:
        return _result("NOT_APPLICABLE", "Variant is on the ClinGen SVI BA1 exception list", rule,
                       exclusions=checked + ["SVI_exception_list_2018"])
    exception = evidence.exception_check
    if exception.state != EvidenceState.KNOWN or exception.is_exception is None:
        return _result("UNDETERMINED", "Applicable BA1 exception status is unverified", rule,
                       exclusions=checked)
    if exception.disease_id not in (None, evidence.disease_id):
        return _result("UNDETERMINED", "Exception check has a different disease scope", rule)
    if not exception.source or not exception.version:
        return _result("UNDETERMINED", "Exception check lacks source or version", rule)
    checked.append("verified_exception_check")
    if exception.is_exception:
        return _result("NOT_APPLICABLE", "Verified BA1 exception applies", rule,
                       exclusions=checked)

    eligible: list[tuple[float, PopulationRecord]] = []
    unresolved = []
    for record in evidence.populations:
        if not _identity_ok(variant, record):
            unresolved.append("population variant/build/normalization mismatch")
            continue
        if rule.database and record.database.lower() != rule.database.lower():
            unresolved.append("population database does not meet VCEP requirement")
            continue
        if rule.continental_only and (record.continental is not True or record.founder_or_bottlenecked is True):
            unresolved.append("population is not a validated general continental group")
            continue
        if rule.minimum_allele_number and (record.allele_number is None
                                           or record.allele_number < rule.minimum_allele_number):
            unresolved.append("qualifying allele number is missing or too small")
            continue
        if (record.filter_pass is not True or record.quality != EvidenceState.KNOWN
                or record.mapping_ambiguous is not False
                or record.callability_adequate is not True
                or record.sampling_bias is True or record.technical_limitation is True):
            unresolved.append("population quality or callability is not established")
            continue
        if variant.variant_type in ("SV", "CNV") and (
                record.population_dataset != "gnomAD SV"
                or record.detection_sensitivity_adequate is not True):
            unresolved.append("structural variant frequency requires validated gnomAD SV detection")
            continue
        value, error = _frequency(record, rule.metric)
        if error or value is None:
            unresolved.append(error or f"{rule.metric} is missing")
            continue
        eligible.append((value, record))

    if not eligible:
        return _result("UNDETERMINED", "; ".join(dict.fromkeys(unresolved)) or "No population evidence supplied", rule,
                       exclusions=checked)
    value, record = max(eligible, key=lambda pair: pair[0])
    features = ("variant_identity", "population_database", "population_quality",
                "allele_number", rule.metric, "exception_check")
    met = value > rule.threshold if rule.operator == ">" else value >= rule.threshold
    if met:
        # A second apparently qualifying but unresolved record needs review if it
        # explicitly signals an artifact or conflicting frequency. A clean record
        # remains sufficient when other records simply lack optional data.
        return _result("ACTIVATED", f"{record.ancestry} {rule.metric} {value:.6g} {rule.operator} {rule.threshold:.6g}",
                       rule, record, value, features, checked + ["population_reliability"])
    if unresolved:
        return _result("UNDETERMINED", "Qualified frequency is below threshold, but other records are unresolved: "
                       + "; ".join(dict.fromkeys(unresolved)), rule, record, value, features, checked)
    return _result("NOT_MET", f"Highest qualified {rule.metric} {value:.6g} does not meet {rule.operator} {rule.threshold:.6g}",
                   rule, record, value, features, checked + ["population_reliability"])
