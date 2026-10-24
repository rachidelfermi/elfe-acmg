"""Deterministic PM1 rules; unknown evidence never means false."""

from .guidelines import ACMG, VCEPS
from .models import BiologicalFeatures, Decision, Evidence, Region, State


def _value(item: Evidence):
    return item.value if item.state == State.KNOWN else None


def _decision(status, reason, guideline=ACMG, path=None, region=None,
              used=(), missing=(), exclusions=(), strength="NONE"):
    return Decision("PM1", status, strength, path, region, tuple(used),
                    tuple(dict.fromkeys(missing)), tuple(exclusions), guideline, reason)


def _vcep(features: BiologicalFeatures, rule: dict) -> Decision:
    if features.transcript_id_version != rule["transcript"]:
        status = "UNDETERMINED" if features.transcript_id_version is None else "NOT_MET"
        return _decision(status, "The VCEP protein coordinates require its specified transcript.",
                         rule, missing=("transcript_id_version",) if status == "UNDETERMINED" else ())
    pos = features.protein_position
    if pos is None:
        return _decision("UNDETERMINED", "Protein position is missing.", rule,
                         missing=("protein_position",))
    if features.variant_consequence is None:
        return _decision("UNDETERMINED", "Variant consequence is missing.", rule,
                         missing=("variant_consequence",))
    if features.variant_consequence.lower() != "missense":
        return _decision("NOT_APPLICABLE", "This implemented VCEP PM1 pathway is for missense variants.", rule)

    gene = features.gene_id
    if gene == "TP53":
        if pos in rule["codons"]:
            return _decision("ACTIVATED", "Variant affects a TP53 VCEP hotspot codon.", rule,
                             "VCEP_HOTSPOT_CODON", f"TP53:{pos}",
                             ("protein_position", "transcript_id_version", "variant_consequence"), strength="MODERATE")
        count = _value(features.somatic_same_amino_acid_observations)
        if count is None or not features.amino_acid_change:
            return _decision("UNDETERMINED", "TP53 somatic hotspot observations are not established.",
                             rule, missing=("somatic_same_amino_acid_observations", "amino_acid_change"))
        if features.somatic_source_version != "cancerhotspots.org:v2":
            return _decision("UNDETERMINED", "The specified somatic hotspot source/version is required.",
                             rule, missing=("somatic_source_version",))
        if count >= 10:
            strength = "MODERATE"
        elif count >= 2:
            strength = "SUPPORTING"
        else:
            return _decision("NOT_MET", "Somatic observation count is below the VCEP minimum.", rule)
        return _decision("ACTIVATED", "Same amino-acid change meets the TP53 VCEP somatic occurrence rule.",
                         rule, "VCEP_SOMATIC_HOTSPOT", f"TP53:{pos}",
                         ("somatic_same_amino_acid_observations", "somatic_source_version",
                          "amino_acid_change"), strength=strength)
    if gene == "DICER1":
        if features.pm5_strength and features.pm5_strength != "NONE":
            return _decision("NOT_MET", "DICER1 VCEP prohibits combining PM1 with PM5.",
                             rule, exclusions=("PM5",))
        if pos in rule["metal_codons"]:
            strength, path = "MODERATE", "VCEP_METAL_BINDING_RESIDUE"
        elif rule["domain"][0] <= pos <= rule["domain"][1]:
            strength, path = "SUPPORTING", "VCEP_RNASE_IIIB_DOMAIN"
        else:
            return _decision("NOT_MET", "Position is outside the DICER1 VCEP PM1 regions.", rule)
        return _decision("ACTIVATED", "Position meets the DICER1 VCEP PM1 rule.", rule,
                         path, f"DICER1:{pos}", ("protein_position", "variant_consequence"),
                         exclusions=("PM5",), strength=strength)
    if gene == "F9":
        if pos in rule["critical_codons"]:
            strength, path = "STRONG", "VCEP_CATALYTIC_OR_ACTIVATION_RESIDUE"
        elif features.exon_number in rule["exons"]:
            strength, path = "MODERATE", "VCEP_EXON_REGION"
        elif features.exon_number is None:
            return _decision("UNDETERMINED", "Exon number is needed for the F9 regional rule.", rule,
                             missing=("exon_number",))
        else:
            return _decision("NOT_MET", "Position is outside the F9 VCEP PM1 rules.", rule)
        if strength == "STRONG" and features.pm5_strength in ("SUPPORTING", "MODERATE", "STRONG"):
            return _decision("UNDETERMINED", "F9 PM1 and PM5 combined weight is capped at Strong; coordinate codes externally.",
                             rule, exclusions=("PM5_COMBINATION_CAP",))
        return _decision("ACTIVATED", "Position meets the F9 VCEP PM1 rule.", rule, path,
                         f"F9:{pos}", ("protein_position", "variant_consequence") +
                         (("exon_number",) if path == "VCEP_EXON_REGION" else ()),
                         exclusions=("PM5_COMBINATION_CAP",), strength=strength)
    raise AssertionError("Unexpected VCEP")


def _region(features: BiologicalFeatures, region: Region) -> Decision:
    if region.start < 1 or region.end < region.start:
        return _decision("UNDETERMINED", "Region protein coordinates are invalid.",
                         region=region.region_id, missing=("valid_region_boundaries",))
    if region.kind == "GENERIC_DOMAIN":
        return _decision("NOT_MET", "A named domain alone does not establish PM1.", region=region.region_id)
    if region.kind not in ("HOTSPOT", "CRITICAL_DOMAIN", "ENRICHMENT"):
        return _decision("NOT_MET", "No verified PM1 pathway for this region type.", region=region.region_id)
    if region.transcript_id_version != features.transcript_id_version or (
        region.protein_accession_version and features.protein_accession_version and
        region.protein_accession_version != features.protein_accession_version
    ):
        return _decision("NOT_MET", "Region and variant protein coordinates are incompatible.",
                         region=region.region_id)
    if region.protein_accession_version and not features.protein_accession_version:
        return _decision("UNDETERMINED", "Protein accession is needed to match region coordinates.",
                         region=region.region_id, missing=("protein_accession_version",))
    if not region.start <= features.protein_position <= region.end:
        return _decision("NOT_MET", "Variant is outside this region.", region=region.region_id)
    if region.kind == "ENRICHMENT":
        return _decision("UNDETERMINED", "ACMG/AMP 2015 supplies no universal enrichment threshold; a scoped VCEP rule is required.",
                         region=region.region_id, missing=("applicable_published_enrichment_rule",))
    if _value(region.excluded) is True:
        return _decision("NOT_MET", "A documented region exclusion applies.", region=region.region_id,
                         exclusions=("region.excluded",))
    if _value(region.lacks_benign_variation) is False:
        return _decision("NOT_MET", "Credible benign variation is present in the region.",
                         region=region.region_id, exclusions=("benign_variation",))
    if _value(region.mechanism_compatible) is False or _value(features.variant_mechanism_compatibility) is False:
        return _decision("NOT_MET", "Region or variant mechanism is incompatible with the disease.",
                         region=region.region_id, exclusions=("mechanism_mismatch",))
    if _value(region.classification_independent) is False:
        return _decision("UNDETERMINED", "Regional pathogenic evidence depends circularly on PM1 or the index variant.",
                         region=region.region_id, missing=("independent_regional_classifications",),
                         exclusions=("classification_circularity",))
    if _value(region.classification_conflicts) is True:
        return _decision("UNDETERMINED", "Reference variant classifications conflict.",
                         region=region.region_id, missing=("resolved_classification_conflicts",))
    if _value(region.ascertainment_bias_resolved) is False:
        return _decision("UNDETERMINED", "Regional ascertainment bias undermines the evidence.",
                         region=region.region_id, missing=("resolved_ascertainment_bias",))
    required = ["established", "disease_relevant", "mechanism_compatible",
                "benign_variation_assessed", "lacks_benign_variation",
                "classification_independent", "evidence_quality"]
    required.append("pathogenic_clustering" if region.kind == "HOTSPOT" else "critical_function")
    missing = [f"region.{name}" for name in required if _value(getattr(region, name)) is None]
    if region.definition_source is None:
        missing.append("region.definition_source")
    if _value(region.population_coverage_adequate) is False:
        missing.append("adequate_population_coverage")
    if any(_value(getattr(region, name)) is False for name in required):
        return _decision("NOT_MET", "A required regional PM1 condition is disproven.",
                         region=region.region_id)
    if missing:
        return _decision("UNDETERMINED", "Essential regional PM1 evidence is missing or conflicting.",
                         region=region.region_id, missing=missing)
    path = "VALIDATED_MUTATIONAL_HOTSPOT" if region.kind == "HOTSPOT" else "ESTABLISHED_CRITICAL_DOMAIN"
    return _decision("ACTIVATED", "Curated region satisfies the ACMG/AMP PM1 definition.",
                     path=path, region=region.region_id,
                     used=("protein_position", "transcript_id_version") + tuple(f"region.{x}" for x in required),
                     exclusions=("region.excluded", "benign_variation", "classification_circularity"),
                     strength="MODERATE")


def evaluate_pm1(biological_features: BiologicalFeatures | dict) -> Decision:
    """Evaluate only PM1 from pre-extracted evidence; return one deterministic decision."""
    f = BiologicalFeatures.from_dict(biological_features) if isinstance(biological_features, dict) else biological_features
    if not f.gene_id or not f.disease_id:
        return _decision("UNDETERMINED", "Gene and disease scope are required.",
                         missing=tuple(x for x in ("gene_id", "disease_id") if not getattr(f, x)))
    if _value(f.gene_disease_validity) is False:
        return _decision("NOT_APPLICABLE", "The gene-disease relationship is not established.")
    if _value(f.mapping_confidence) is False:
        return _decision("UNDETERMINED", "Protein position mapping is unreliable.", missing=("mapping_confidence",))
    rule = VCEPS.get((f.gene_id, f.disease_id))
    if rule:
        return _vcep(f, rule)
    if f.gene_id in {gene for gene, _ in VCEPS}:
        return _decision("UNDETERMINED", "This gene has a VCEP PM1 specification, but the supplied disease is outside its implemented scope.",
                         missing=("applicable_disease_specification",))
    if f.protein_position is None or f.transcript_id_version is None:
        return _decision("UNDETERMINED", "Protein location and transcript are required.",
                         missing=tuple(x for x in ("protein_position", "transcript_id_version") if getattr(f, x) is None))
    if f.variant_consequence is None:
        return _decision("UNDETERMINED", "Variant consequence is required for the implemented protein PM1 pathways.",
                         missing=("variant_consequence",))
    if f.variant_consequence.lower() != "missense":
        return _decision("NOT_APPLICABLE", "The implemented ACMG/AMP protein PM1 pathways are for missense variants.")
    if not f.regions:
        return _decision("UNDETERMINED", "No curated regional evidence was supplied.",
                         missing=("regions",))
    by_id = {}
    for region in f.regions:
        previous = by_id.setdefault(region.region_id, region)
        if previous != region:
            return _decision("UNDETERMINED", "Conflicting annotations were supplied for one region ID.",
                             region=region.region_id, missing=("reconciled_region_annotation",))
    results = [_region(f, r) for r in f.regions]
    for result in results:
        if result.status == "ACTIVATED":
            return result
    unresolved = [r for r in results if r.status == "UNDETERMINED"]
    if unresolved:
        return _decision("UNDETERMINED", "At least one potentially applicable region remains unresolved.",
                         missing=tuple(x for r in unresolved for x in r.missing_critical_features))
    return _decision("NOT_MET", "No supplied region satisfies a verified PM1 pathway.")
