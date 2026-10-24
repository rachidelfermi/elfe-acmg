"""Conservative BS4 decisions from already verified family evidence."""

from .guidelines import BRCA2, brca2_lr_strength, select_guideline
from .models import BS4Input, Decision, Evidence, Observation, State


def _known(value: Evidence, expected: object) -> bool:
    return value.state == State.KNOWN and value.value == expected


def _missing(value: Evidence) -> bool:
    return value.state != State.KNOWN


def _result(context: BS4Input, guideline, status: str, reason: str,
            *, strength: str = "NONE", path: str = "NONE", events=(),
            missing=(), exclusions=(), features=(), conflicts=()) -> Decision:
    return Decision("BS4", status, strength, path, len(events),
                    tuple(sorted({o.family_id for o in events})),
                    tuple(sorted(set(features))), tuple(sorted(set(missing))),
                    tuple(sorted(set(exclusions))), guideline.citation(), reason,
                    tuple(sorted(set(conflicts))), tuple(context.other_criteria_overlap))


def _observation(observation: Observation, guideline) -> tuple[str, list[str], list[str]]:
    """Return informative, excluded, or incomplete, plus reasons and assessed exclusions."""
    o = observation
    if _known(o.affected, False) or _known(o.variant_present, True):
        return "not_informative", [], []
    exclusions = []
    for name in ("alternative_diagnosis", "phenocopy_evidence",
                 "genetic_heterogeneity", "locus_heterogeneity",
                 "phenotype_misclassification", "mosaicism_explanation",
                 "genotyping_artifact", "copy_number_conflict", "phasing_conflict",
                 "recombination_explanation", "sample_identity_conflict",
                 "pedigree_conflict"):
        value = getattr(o, name)
        if _known(value, True):
            return "excluded", [name], [name]
        if _known(value, False):
            exclusions.append(name)
    required = ("affected", "variant_present", "genotype_reliable",
                "phenotype_reliable", "phenotype_consistent",
                "biological_relationship_confirmed", "phenocopy_excluded",
                "heterogeneity_excluded")
    missing = []
    for name in required:
        value = getattr(o, name)
        if _missing(value):
            missing.append(f"{o.family_id}/{o.individual_id}.{name}")
        elif name not in ("variant_present",) and value.value is not True:
            return "excluded", [name], exclusions
        elif name == "variant_present" and value.value is not False:
            return "excluded", [name], exclusions
    if guideline.phenotype_score_minimum is not None:
        if _missing(o.phenotype_score):
            missing.append(f"{o.family_id}/{o.individual_id}.phenotype_score")
        elif not isinstance(o.phenotype_score.value, (int, float)) or isinstance(o.phenotype_score.value, bool):
            return "excluded", ["invalid_phenotype_score"], exclusions
        elif o.phenotype_score.value < guideline.phenotype_score_minimum:
            return "excluded", ["phenotype_score_below_vcep_threshold"], exclusions
    if o.evidence_quality.state == State.CONFLICTING:
        missing.append(f"{o.family_id}/{o.individual_id}.evidence_quality")
    elif _known(o.evidence_quality, False):
        return "excluded", ["evidence_quality"], exclusions
    if missing:
        return "incomplete", missing, exclusions
    return "informative", [], exclusions


def _quantitative(context: BS4Input, guideline) -> Decision:
    if not guideline.quantitative:
        return _result(context, guideline, "NOT_MET",
                       "This guideline has no verified quantitative BS4 threshold.")
    required = ("quantitative_likelihood_ratio", "quantitative_validated",
                "quantitative_model_appropriate", "quantitative_independent")
    missing = [name for name in required if _missing(getattr(context, name))]
    if missing:
        return _result(context, guideline, "UNDETERMINED",
                       "Validated, independent likelihood evidence and model assumptions are incomplete.",
                       missing=missing)
    if not all(_known(getattr(context, name), True) for name in required[1:]):
        return _result(context, guideline, "NOT_MET",
                       "Likelihood evidence failed validation, independence, or model checks.")
    lr = context.quantitative_likelihood_ratio.value
    if not isinstance(lr, (int, float)) or isinstance(lr, bool) or not 0 <= lr < float("inf"):
        return _result(context, guideline, "UNDETERMINED",
                       "Likelihood ratio must be a finite nonnegative number.",
                       missing=("quantitative_likelihood_ratio",))
    strength = brca2_lr_strength(lr)
    if strength is None:
        return _result(context, guideline, "NOT_MET",
                       "The supplied likelihood ratio does not reach the BRCA2 BS4 Supporting threshold.")
    return _result(context, guideline, "ACTIVATED",
                   f"Validated BRCA2 segregation LR {lr:g} meets the {strength} threshold in GN097 v1.2.",
                   strength=strength, path="QUANTITATIVE_NONSEGREGATION",
                   features=required)


def evaluate_bs4(biological_features: BS4Input | dict) -> Decision:
    """Evaluate one variant/disease pair; never derive genotypes or clinical status."""
    context = (biological_features if isinstance(biological_features, BS4Input)
               else BS4Input.from_mapping(biological_features))
    guideline = select_guideline(context.gene_id, context.disease_id)
    if (guideline.gene_id is not None and context.vcep_specification_version is not None
            and context.vcep_specification_version != guideline.version):
        return _result(context, guideline, "UNDETERMINED",
                       "The supplied VCEP version differs from the pinned specification.",
                       missing=("vcep_specification_version",))
    identity_missing = [name for name in ("variant_id", "gene_id", "disease_id")
                        if not getattr(context, name)]
    if identity_missing:
        return _result(context, guideline, "UNDETERMINED",
                       "Variant, gene, and disease scope must be identified.", missing=identity_missing)
    if _known(context.gene_disease_validity, False):
        return _result(context, guideline, "NOT_APPLICABLE",
                       "The supplied gene–disease relationship is invalid.")
    if _known(context.inheritance_compatible, False):
        return _result(context, guideline, "NOT_APPLICABLE",
                       "The inheritance mechanism does not predict this segregation.")
    if _known(context.variant_reliable, False):
        return _result(context, guideline, "NOT_MET", "The candidate variant call is unreliable.")
    for name in ("variant_mapping_confidence", "variant_call_quality"):
        value = getattr(context, name)
        if _known(value, False):
            return _result(context, guideline, "NOT_MET",
                           f"The supplied {name} assessment is unreliable.")
        if value.state == State.CONFLICTING:
            return _result(context, guideline, "UNDETERMINED",
                           f"The supplied {name} assessment conflicts.", missing=(name,))
    context_missing = [name for name in ("gene_disease_validity", "inheritance_compatible",
                                       "variant_reliable") if _missing(getattr(context, name))]
    if context_missing:
        return _result(context, guideline, "UNDETERMINED",
                       "Essential variant and disease context is missing or conflicting.",
                       missing=context_missing)
    if _known(context.pp1_conflict, True) or context.pp1_conflict.state == State.CONFLICTING:
        return _result(context, guideline, "UNDETERMINED",
                       "Credible PP1/BS4 family evidence conflicts and needs resolution.",
                       conflicts=("PP1_BS4",))
    if _known(context.family_overlap, True) or context.family_overlap.state == State.CONFLICTING:
        return _result(context, guideline, "UNDETERMINED",
                       "Family overlap or duplicate reports are unresolved.",
                       conflicts=("FAMILY_OVERLAP",))
    if guideline == BRCA2:
        return _quantitative(context, guideline)

    if not context.observations:
        return _result(context, guideline, "UNDETERMINED",
                       "No individual-level segregation observations were supplied.",
                       missing=("observations",))
    unique = {}
    duplicate_conflict = False
    for o in context.observations:
        key = (o.family_id, o.individual_id)
        if key in unique and o != unique[key]:
            duplicate_conflict = True
        unique.setdefault(key, o)
    if duplicate_conflict:
        return _result(context, guideline, "UNDETERMINED",
                       "Duplicate individual records disagree.", conflicts=("DUPLICATE_RECORD",))
    informative, missing, exclusions = [], [], []
    incomplete = False
    for o in unique.values():
        state, details, assessed = _observation(o, guideline)
        exclusions.extend(assessed)
        if state == "informative":
            informative.append(o)
        elif state == "incomplete":
            incomplete = True
            missing.extend(details)
        elif state == "excluded":
            exclusions.extend(details)
    if informative:
        strength = (guideline.multiple_event_strength if len(informative) >= 2
                    else guideline.one_event_strength)
        path = ("MULTIPLE_NONSEGREGATION" if len(informative) >= 2
                else "CONFIRMED_NONSEGREGATION")
        return _result(context, guideline, "ACTIVATED",
                       f"{len(informative)} distinct affected noncarrier(s) satisfy {guideline.rule_id}; "
                       f"applicable strength is {strength}.", strength=strength, path=path,
                       events=informative, exclusions=exclusions,
                       features=("affected", "variant_present", "genotype_reliable",
                                 "phenotype_reliable", "phenotype_consistent",
                                 "biological_relationship_confirmed", "phenocopy_excluded",
                                 "heterogeneity_excluded") +
                                (("phenotype_score",) if guideline.phenotype_score_minimum else ()))
    if incomplete:
        return _result(context, guideline, "UNDETERMINED",
                       "Potential nonsegregation cannot be established from the supplied evidence.",
                       missing=missing, exclusions=exclusions)
    return _result(context, guideline, "NOT_MET",
                   "No reliable, eligible affected noncarrier observation was found.",
                   exclusions=exclusions)
