"""Conservative, adaptive evaluation of BP1 only."""

from typing import Any, Mapping

from .guidelines import ACMG_2015, Specification, find_specification
from .models import BP1Features, BP1Result, EvidenceState


def evaluate_bp1(features: BP1Features | Mapping[str, Any]) -> BP1Result:
    if isinstance(features, Mapping):
        features = BP1Features.from_mapping(features)
    if not isinstance(features, BP1Features):
        raise TypeError("features must be BP1Features or a mapping")

    used: list[str] = []
    missing: list[str] = []
    exclusions: list[str] = []

    def value(name: str) -> Any:
        evidence = features.evidence(name)
        used.append(name)
        if evidence.state == EvidenceState.KNOWN:
            return evidence.value
        missing.append(name)
        return None

    def result(status: str, reason: str, path: str | None = None,
               spec: Specification | None = None, strength: str = "NONE") -> BP1Result:
        return BP1Result("BP1", status, strength, path, tuple(dict.fromkeys(used)),
                         tuple(dict.fromkeys(missing)), tuple(dict.fromkeys(exclusions)),
                         spec.guideline if spec else ACMG_2015.copy(), reason)

    gene = value("gene_id")
    disease = value("disease_id")
    if gene is None or disease is None:
        return result("UNDETERMINED", "A gene and a disease are required to choose a BP1 rule.")
    gene, disease = str(gene).upper(), str(disease).upper()
    spec = find_specification(gene, disease)
    requested = features.evidence("vcep_specification_id")
    if requested.state == EvidenceState.KNOWN and requested.value != (spec.spec_id if spec else None):
        return result("UNDETERMINED", "Requested VCEP specification does not match the pinned gene–disease scope.")
    if spec:
        version = features.evidence("vcep_specification_version")
        if version.state == EvidenceState.KNOWN and str(version.value) != spec.version:
            return result("UNDETERMINED", "Requested VCEP version is not the pinned verified version.", spec=spec)
        if spec.mode == "DISABLED":
            used.append("vcep_specification_id")
            exclusions.append("VCEP_BP1_NOT_APPLICABLE")
            return result("NOT_APPLICABLE", "This VCEP specification marks BP1 not applicable for this gene–disease pair.", spec=spec)

    validity = value("gene_disease_validity")
    if validity is False:
        return result("NOT_MET", "The supplied gene–disease relationship is not valid.", spec=spec)
    if validity is None:
        return result("UNDETERMINED", "Gene–disease validity is unresolved.", spec=spec)

    transcript = value("biologically_relevant_transcript")
    if transcript is False:
        return result("NOT_MET", "The consequence is not on a biologically relevant transcript.", spec=spec)
    if transcript is None:
        return result("UNDETERMINED", "Relevant transcript status is unresolved.", spec=spec)

    consequence = value("variant_consequence")
    if consequence is None:
        return result("UNDETERMINED", "Variant consequence is unresolved.", spec=spec)
    consequence = str(consequence).upper()

    if spec and spec.mode == "RASOPATHY_LOF":
        return _evaluate_rasopathy(features, consequence, spec, value, result)

    if consequence != "MISSENSE":
        return result("NOT_MET", "General BP1 requires a missense consequence.")

    prior = features.evidence("prior_criteria")
    if prior.state == EvidenceState.KNOWN:
        used.append("prior_criteria")
        if not isinstance(prior.value, (list, tuple, set, frozenset)) or not all(isinstance(x, str) for x in prior.value):
            raise ValueError("prior_criteria must be a sequence of criterion names")
        if "PP2" in {x.upper() for x in prior.value}:
            exclusions.append("PP2_PREMISE_CONFLICT")
            return result("UNDETERMINED", "A prior PP2 assignment conflicts with BP1's disease-mechanism premise; resolve the evidence before applying BP1.")

    confidence = features.evidence("missense_consequence_confidence")
    if confidence.state == EvidenceState.KNOWN:
        used.append("missense_consequence_confidence")
        if str(confidence.value).upper() == "LOW":
            return result("UNDETERMINED", "The missense consequence has low confidence.")
    elif confidence.state == EvidenceState.CONFLICTING:
        missing.append("missense_consequence_confidence")
        return result("UNDETERMINED", "The missense consequence is conflicting.")

    # Observed decisive exceptions override general gene-level assessments.
    for field, label in (
        ("candidate_is_known_pathogenic_missense", "known pathogenic missense variant"),
        ("candidate_region_missense_pathogenicity", "pathogenic missense region"),
        ("candidate_residue_missense_pathogenicity", "pathogenic missense residue"),
        ("missense_specific_disease_subtype", "missense-specific disease subtype"),
    ):
        e = features.evidence(field)
        if e.state == EvidenceState.KNOWN:
            used.append(field)
            exclusions.append(field)
            if e.value is True:
                return result("NOT_MET", f"Candidate overlaps a {label}.")
        elif e.state == EvidenceState.CONFLICTING:
            missing.append(field)
            return result("UNDETERMINED", f"Evidence for {label} is conflicting.")

    for field in ("missense_gain_of_function_mechanism", "missense_dominant_negative_mechanism"):
        e = features.evidence(field)
        if e.state == EvidenceState.KNOWN and e.value is True:
            used.append(field)
            exclusions.append(field)
            return result("NOT_MET", "An established missense disease mechanism is present.")
        if e.state == EvidenceState.CONFLICTING:
            missing.append(field)
            return result("UNDETERMINED", "The missense disease mechanism is conflicting.")

    splice = features.evidence("candidate_splice_effect")
    if splice.state == EvidenceState.KNOWN:
        used.append("candidate_splice_effect")
        exclusions.append("candidate_splice_effect")
        if str(splice.value).upper() in {"DAMAGING", "PATHOGENIC"}:
            return result("NOT_MET", "Established damaging splicing supersedes a purely missense BP1 assessment.")
    elif splice.state == EvidenceState.CONFLICTING:
        missing.append("candidate_splice_effect")
        return result("UNDETERMINED", "Candidate splice evidence is conflicting.")
    prediction = features.evidence("candidate_splice_prediction")
    if prediction.state == EvidenceState.KNOWN and str(prediction.value).upper() in {"HIGH_RISK", "DAMAGING"}:
        used.append("candidate_splice_prediction")
        missing.append("candidate_splice_effect")
        return result("UNDETERMINED", "A credible splice prediction requires resolution of the actual effect.")
    alternative = features.evidence("alternative_transcript_consequence")
    if alternative.state == EvidenceState.KNOWN and str(alternative.value).upper() in {"TRUNCATING", "CANONICAL_SPLICE"}:
        used.append("alternative_transcript_consequence")
        return result("UNDETERMINED", "A relevant alternative transcript has a potentially truncating consequence.")
    alternative_splice = features.evidence("alternative_splice_consequence")
    if alternative_splice.state == EvidenceState.KNOWN and str(alternative_splice.value).upper() in {"DAMAGING", "PATHOGENIC"}:
        used.append("alternative_splice_consequence")
        return result("UNDETERMINED", "A relevant alternative splice consequence requires mechanism assessment.")

    assessed = value("exceptions_assessed")
    if assessed is False:
        return result("UNDETERMINED", "Missense mechanism exceptions have not been assessed.")
    if assessed is None:
        return result("UNDETERMINED", "A documented exception assessment is required.")

    missense = value("missense_common_mechanism")
    if missense is True:
        return result("NOT_MET", "Missense variants are an established common disease mechanism.")
    if missense is None:
        return result("UNDETERMINED", "Whether missense is a common disease mechanism is unresolved.")

    predominant = features.evidence("truncating_predominant")
    spectrum = features.evidence("spectrum_truncating_predominant")
    if predominant.state == EvidenceState.CONFLICTING or spectrum.state == EvidenceState.CONFLICTING:
        missing.append("truncating_predominant" if predominant.state == EvidenceState.CONFLICTING else "spectrum_truncating_predominant")
        return result("UNDETERMINED", "Evidence about truncating predominance is conflicting.")
    if predominant.state == EvidenceState.KNOWN and predominant.value is False and spectrum.state == EvidenceState.KNOWN and spectrum.value is False:
        used.extend(("truncating_predominant", "spectrum_truncating_predominant"))
        return result("NOT_MET", "Both supplied assessments reject truncating predominance.")
    if predominant.state == EvidenceState.KNOWN and predominant.value is True:
        used.append("truncating_predominant")
        quality = value("mechanism_evidence_quality")
        source = value("mechanism_evidence_source")
        if quality not in {"CURATED", "HIGH"} or not source:
            missing.extend(("mechanism_evidence_quality", "mechanism_evidence_source"))
            return result("UNDETERMINED", "Truncating predominance needs a cited, curated disease-specific assessment.")
        if spectrum.state == EvidenceState.KNOWN and spectrum.value is False:
            used.append("spectrum_truncating_predominant")
            return result("UNDETERMINED", "Mechanism and spectrum assessments conflict.")
        return result("ACTIVATED", "A cited, curated assessment establishes truncating predominance for this disease; missense is not a common mechanism and assessed exceptions do not apply.",
                      "TRUNCATING_DISEASE_MECHANISM", strength="SUPPORTING")
    if predominant.state == EvidenceState.KNOWN and predominant.value is False:
        used.append("truncating_predominant")
        if spectrum.state != EvidenceState.KNOWN or spectrum.value is not True:
            return result("NOT_MET", "Truncating variants are not the predominant pathogenic mechanism.")
        return result("UNDETERMINED", "Mechanism and spectrum assessments conflict.")

    if spectrum.state == EvidenceState.KNOWN and spectrum.value is True:
        used.append("spectrum_truncating_predominant")
        source = value("mechanism_evidence_source")
        if not source:
            missing.append("mechanism_evidence_source")
            return result("UNDETERMINED", "The curated spectrum assessment needs a cited source.")
        checks = ("spectrum_evidence_quality", "spectrum_classification_independent", "spectrum_unique_variants",
                  "spectrum_disease_match", "spectrum_ascertainment_bias_addressed")
        vals = {name: value(name) for name in checks}
        if vals["spectrum_classification_independent"] is False:
            return result("UNDETERMINED", "Reference classifications use circular BP1 evidence.")
        if vals["spectrum_unique_variants"] is False:
            return result("UNDETERMINED", "Duplicated submissions are not independent pathogenic variants.")
        if vals["spectrum_disease_match"] is False:
            return result("NOT_MET", "The pathogenic spectrum is for another disease.")
        if vals["spectrum_ascertainment_bias_addressed"] is False:
            return result("UNDETERMINED", "Spectrum ascertainment bias is unresolved.")
        if vals["spectrum_evidence_quality"] not in {"CURATED", "HIGH"} or any(vals[x] is None for x in checks[1:]):
            return result("UNDETERMINED", "Curated, independent, disease-matched spectrum evidence is incomplete.")
        return result("ACTIVATED", "A curated disease-matched pathogenic spectrum establishes truncating predominance without circular classifications; missense is not a common mechanism.",
                      "CURATED_PATHOGENIC_SPECTRUM", strength="SUPPORTING")
    if spectrum.state == EvidenceState.KNOWN and spectrum.value is False:
        used.append("spectrum_truncating_predominant")
        return result("NOT_MET", "Curated pathogenic spectrum does not show truncating predominance.")
    missing.extend(("truncating_predominant", "spectrum_truncating_predominant"))
    return result("UNDETERMINED", "No verified disease-specific assessment establishes truncating predominance; raw counts alone do not suffice.")


def _evaluate_rasopathy(features: BP1Features, consequence: str, spec: Specification,
                        value: Any, result: Any) -> BP1Result:
    allowed = {"STOP_GAINED", "NONSENSE", "FRAMESHIFT", "CANONICAL_SPLICE", "START_LOST",
               "INITIATION_CODON", "WHOLE_GENE_DELETION", "MULTI_EXON_DELETION"}
    if consequence not in allowed:
        return result("NOT_MET", "This RASopathy VCEP BP1 rule requires a specified truncating/LOF variant class.", spec=spec)
    lof = value("loss_of_function_disease_mechanism")
    if lof is True:
        return result("NOT_MET", "An established LOF correlation to this disease excludes this VCEP BP1 rule.", spec=spec)
    if lof is None:
        return result("UNDETERMINED", "LOF correlation to this disease must be resolved.", spec=spec)
    gof = value("missense_gain_of_function_mechanism")
    if gof is False:
        return result("NOT_MET", "The VCEP rule requires a missense gain-of-function disease mechanism.", spec=spec)
    if gof is None:
        return result("UNDETERMINED", "The RASopathy gain-of-function mechanism is unresolved.", spec=spec)
    return result("ACTIVATED", "MAP2K1 RASopathy VCEP BP1 applies to this truncating class with a gain-of-function disease mechanism and no established LOF correlation.",
                  "VCEP_RASOPATHY_LOF", spec=spec, strength=spec.strength)
