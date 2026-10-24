"""Evaluate only BP7; input evidence is never fetched or inferred."""
from math import isfinite
from .guidelines import GUIDELINES
from .models import BP7Features, BP7Result, EvidenceState, RNAEvidence, SplicePrediction


def evaluate_bp7(features: BP7Features | dict) -> BP7Result:
    """Return a deterministic criterion decision for a pinned guideline.

    Dict input uses the same keys as BP7Features, with nested `prediction`/`rna`
    dictionaries. Unknown keys and invalid values raise ValueError/TypeError.
    """
    if isinstance(features, dict):
        data = dict(features)
        if isinstance(data.get("prediction"), dict):
            data["prediction"] = SplicePrediction(**data["prediction"])
        if isinstance(data.get("rna"), dict):
            data["rna"] = RNAEvidence(**data["rna"])
        features = BP7Features(**data)
    if not isinstance(features, BP7Features):
        raise TypeError("features must be BP7Features or a matching dictionary")
    if features.guideline not in GUIDELINES:
        raise ValueError("Unknown guideline; select a pinned rule")
    rule = GUIDELINES[features.guideline]
    citation = {"name": rule.name, "version": rule.version, "rule_id": rule.rule_id, "source": rule.source}
    used: list[str] = []
    exclusions: list[str] = []

    def out(status, reason, *, strength="NONE", path=None, missing=(), note=None):
        return BP7Result("BP7", status, strength, path, tuple(dict.fromkeys(used)), tuple(missing), tuple(dict.fromkeys(exclusions)), citation, reason, note)

    if rule.name == "ClinGen TP53 VCEP":
        if features.gene_id not in ("TP53", "HGNC:11998") or features.disease_id not in ("MONDO:0018875", "Li-Fraumeni syndrome"):
            return out("NOT_APPLICABLE", "TP53 VCEP rule requires TP53 and Li-Fraumeni syndrome scope.")
        used.extend(("gene_id", "disease_id"))

    missing = []
    for field in ("variant_consequence", "consequence_confident", "biologically_relevant_transcript", "variant_mapping_valid", "transcript_id_version"):
        if getattr(features, field) is None:
            missing.append(field)
    if missing:
        return out("UNDETERMINED", "Variant consequence, mapping, and transcript relevance must be established.", missing=missing)
    used.extend(("variant_consequence", "consequence_confident", "biologically_relevant_transcript", "variant_mapping_valid", "transcript_id_version"))
    if features.variant_mapping_valid is False or features.consequence_confident is False:
        return out("UNDETERMINED", "Variant mapping or transcript-specific consequence is unreliable.")
    if features.biologically_relevant_transcript is False:
        return out("NOT_APPLICABLE", "The assessed transcript is not biologically relevant.")
    if features.variant_consequence not in (("synonymous", "intronic") if rule.permits_intronic else ("synonymous",)):
        return out("NOT_APPLICABLE", "Variant class is not eligible under the selected rule.")

    for field in ("established_regulatory_splice_effect", "alternative_transcript_effect", "known_pathogenic_synonymous_mechanism", "established_alternative_pathogenic_effect"):
        value = getattr(features, field)
        if value is not None:
            exclusions.append(field)
            if value:
                used.append(field)
                return out("NOT_MET", f"Credible contradictory evidence: {field}.")
    rna = features.rna
    if rna is not None:
        exclusions.append("rna_splicing_evidence")
        if rna.conflict or (rna.normal_splicing is True and rna.abnormal_splicing is True):
            return out("UNDETERMINED", "RNA evidence is conflicting.")
        if rna.abnormal_splicing is True:
            if (rna.validated is True and rna.controls is True and rna.relevant_tissue is True and
                    rna.variant_match is True and rna.transcript_id_version == features.transcript_id_version):
                used.append("rna.abnormal_splicing")
                return out("NOT_MET", "Validated, matching RNA evidence shows abnormal splicing.")
            return out("UNDETERMINED", "Possible abnormal RNA splicing requires resolution before BP7 can be assigned.")

    # RNA is a distinct experimental path. No arbitrary fraction or assay-quality weight is inferred.
    rna_ready = bool(rna and rna.normal_splicing is True and rna.abnormal_splicing is False and
                     rna.validated is True and rna.controls is True and rna.relevant_tissue is True and
                     rna.adequate_sensitivity is True and rna.variant_match is True and
                     rna.transcript_id_version == features.transcript_id_version and not rna.already_used_for_bs3)
    if rna and rna.normal_splicing is True and not rna_ready:
        exclusions.append("rna_assay_reliability_or_overlap")

    if features.variant_consequence == "synonymous":
        if rule.requires_conservation:
            pass  # ACMG 2015 does not define a universal splice-position cutoff.
        elif features.exon_core_splice_motif is None and not rna_ready:
            return out("UNDETERMINED", "Exon position is needed for the computational BP7 path.", missing=("exon_core_splice_motif",))
        elif features.exon_core_splice_motif is True and not rna_ready:
            return out("NOT_MET", "Synonymous variant lies in the first or last three exon bases.")
    else:
        pos = features.intronic_position
        if pos is None and not rna_ready:
            return out("UNDETERMINED", "Intronic donor/acceptor position is needed.", missing=("intronic_position",))
        if pos is not None and not (pos >= 7 or pos <= -21) and not rna_ready:
            return out("NOT_MET", "Intronic variant is within the designated splice region.")
    if features.canonical_splice_site_overlap is True and not rna_ready:
        return out("NOT_MET", "Canonical splice-site overlap excludes computational BP7.")
    if features.regulatory_overlap is True and features.established_regulatory_splice_effect is None and not rna_ready:
        return out("UNDETERMINED", "A known regulatory overlap needs assessment of its splice effect.", missing=("established_regulatory_splice_effect",))

    prediction = features.prediction
    prediction_valid = False
    if prediction:
        if prediction.conflict or prediction.quality == EvidenceState.CONFLICTING:
            return out("UNDETERMINED", "Splicing predictions conflict.")
        prediction_valid = (prediction.variant_match is True and prediction.applicable is True and
                            prediction.calibrated is True and prediction.quality == EvidenceState.KNOWN and
                            prediction.transcript_id_version == features.transcript_id_version and
                            prediction.genome_assembly is not None and prediction.genome_assembly == features.genome_assembly)
    if prediction and not prediction_valid and not rna_ready:
        return out("UNDETERMINED", "Splice prediction lacks valid calibration, mapping, build, transcript, or quality.", missing=("valid_splice_prediction",))

    benign_prediction = False
    if prediction_valid:
        used.append("prediction")
        if rule.spliceai_cutoff is None:
            if prediction.no_splice_effect is None:
                return out("UNDETERMINED", "Original ACMG requires a validated no-splice-effect assessment.", missing=("prediction.no_splice_effect",))
            benign_prediction = prediction.no_splice_effect is True
        else:
            if prediction.predictor != "SpliceAI" or prediction.score is None:
                if not rna_ready:
                    return out("UNDETERMINED", "This pinned rule requires a complete SpliceAI delta score.", missing=("prediction.score",))
            else:
                if not isfinite(prediction.score) or not 0 <= prediction.score <= 1:
                    raise ValueError("SpliceAI delta score must be finite and between 0 and 1")
                benign_prediction = prediction.score <= rule.spliceai_cutoff
                used.append("prediction.score")

    if rna_ready and prediction_valid and not benign_prediction and rule.spliceai_cutoff is not None:
        return out("UNDETERMINED", "Normal RNA and computational splice evidence disagree; resolve the conflict before assigning BP7.")

    if rna_ready and rule.permits_rna_strong and (features.guideline != "TP53_2.1.0" or benign_prediction):
        used.append("rna.normal_splicing")
        return out("ACTIVATED", "Validated RNA assay with controls shows no variant-associated abnormal splicing.", strength="STRONG", path="NO_ABNORMAL_SPLICING_RNA", note="Use this assay under BP7_Strong (RNA), not BS3; do not independently count the same assay.")
    if not prediction:
        return out("UNDETERMINED", "No computational splicing assessment was supplied.", missing=("prediction",))
    if not prediction_valid:
        return out("UNDETERMINED", "A valid splice prediction is required.", missing=("valid_splice_prediction",))
    if not benign_prediction:
        return out("NOT_MET", "Computational evidence does not meet the pinned benign splicing rule.")

    if rule.requires_conservation:
        if features.conservation in (None, "unknown", "conflicting") or features.conservation_reliable is not True:
            return out("UNDETERMINED", "Reliable nucleotide conservation assessment is required by ACMG 2015.", missing=("reliable_conservation",))
        used.append("conservation")
        if features.conservation == "highly_conserved":
            return out("NOT_MET", "Nucleotide is highly conserved under the supplied reliable assessment.")
        if features.conservation != "not_highly_conserved":
            raise ValueError("Invalid conservation category")

    note = "BP4 is a prerequisite in this rule and may accompany BP7; this evaluator does not emit BP4."
    return out("ACTIVATED", "Eligible variant has a validated benign splice assessment and meets the selected BP7 rule.", strength="SUPPORTING", path="SYNONYMOUS_NO_SPLICE_EFFECT" if features.variant_consequence == "synonymous" else "INTRONIC_NO_SPLICE_EFFECT", note=note if not rule.requires_conservation else None)
