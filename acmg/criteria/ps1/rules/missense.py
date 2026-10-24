"""ACMG/AMP 2015 standard (amino-acid equivalence) PS1 pathway."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from ..evidence import Evidence
from ..guidelines import ACMG_PS1
from ..models import (
    BiologicalFeatures,
    DecisionStatus,
    MolecularMechanism,
    PS1Strength,
    ReferenceClassification,
    ReferenceVariant,
    SpliceImpact,
    SplicePrediction,
    VariantAnnotations,
    VariantConsequence,
)
from ..output import RuleEvaluation
from .common import (
    ReferencePreconditions,
    activated,
    common_reference_preconditions,
    exclusion,
    label,
    normalized,
    not_applicable,
    not_met,
    pairwise_evidence,
    require_reference_evidence_quality,
    undetermined,
)


PATH = "SAME_AMINO_ACID_CHANGE"


@dataclass(frozen=True)
class _Assessment:
    state: str  # OK, NOT_MET, UNDETERMINED
    used: tuple[str, ...] = ()
    missing: tuple[str, ...] = ()
    reason: str = ""
    exclusions: tuple[dict[str, str], ...] = ()


def _known_strings_match(left: Evidence[str], right: Evidence[str]) -> Optional[bool]:
    if left.is_conflicting or right.is_conflicting:
        return None
    if not left.is_known or not right.is_known:
        return None
    return normalized(left.value or "") == normalized(right.value or "")


def _protein_equivalence(
    evidence: BiologicalFeatures, reference: ReferenceVariant, reference_index: int
) -> _Assessment:
    """Validate protein substitution equivalence without trusting HGVS alone."""

    candidate = evidence.candidate
    ref = reference.annotations
    used: list[str] = []
    missing: list[str] = []

    explicit, explicit_label = pairwise_evidence(
        evidence, reference, reference_index, "protein_substitution_equivalence"
    )
    if explicit.is_conflicting:
        return _Assessment(
            "UNDETERMINED", (explicit_label,),
            (explicit_label,),
            "The supplied protein-substitution-equivalence assertion is conflicting.",
        )

    comparison_pairs = (
        ("candidate_protein_hgvs", candidate.protein_hgvs, label(reference_index, "protein_hgvs"), ref.protein_hgvs),
        ("candidate_amino_acid_position", candidate.amino_acid_position, label(reference_index, "amino_acid_position"), ref.amino_acid_position),
        ("candidate_reference_amino_acid", candidate.reference_amino_acid, label(reference_index, "reference_amino_acid"), ref.reference_amino_acid),
        ("candidate_alternate_amino_acid", candidate.alternate_amino_acid, label(reference_index, "alternate_amino_acid"), ref.alternate_amino_acid),
    )
    known_mismatch = False
    protein_effect_known_match = False
    components_known = True
    components_equal = True
    for candidate_name, candidate_value, reference_name, reference_value in comparison_pairs:
        if candidate_value.is_conflicting or reference_value.is_conflicting:
            return _Assessment(
                "UNDETERMINED", tuple(used), (candidate_name if candidate_value.is_conflicting else reference_name,),
                "Protein-effect annotations are conflicting.",
            )
        if candidate_value.is_known and reference_value.is_known:
            used.extend((candidate_name, reference_name))
            left = candidate_value.value
            right = reference_value.value
            match = normalized(left) == normalized(right) if isinstance(left, str) and isinstance(right, str) else left == right
            if not match:
                known_mismatch = True
            if candidate_name == "candidate_protein_hgvs" and match:
                protein_effect_known_match = True
        elif candidate_name != "candidate_protein_hgvs":
            components_known = False
            components_equal = False

    if components_known and not known_mismatch:
        protein_effect_known_match = True
    if known_mismatch:
        if explicit.is_true:
            return _Assessment(
                "UNDETERMINED", tuple(used) + (explicit_label,),
                (explicit_label,),
                "The declared protein equivalence conflicts with normalized amino-acid annotations.",
            )
        return _Assessment(
            "NOT_MET", tuple(used), (),
            "Candidate and reference do not produce the same amino-acid substitution.",
            (exclusion(explicit_label, "DIFFERENT_PROTEIN_EFFECT"),),
        )

    if explicit.is_false:
        return _Assessment(
            "NOT_MET", tuple(used) + (explicit_label,), (),
            "Validated annotation states that the protein substitutions are not equivalent.",
            (exclusion(explicit_label, "NOT_EQUIVALENT"),),
        )

    transcript_match = _known_strings_match(candidate.transcript_id_version, ref.transcript_id_version)
    accession_match = _known_strings_match(candidate.protein_accession_version, ref.protein_accession_version)
    if candidate.transcript_id_version.is_conflicting or ref.transcript_id_version.is_conflicting:
        return _Assessment(
            "UNDETERMINED", tuple(used), ("candidate_transcript_id_version", label(reference_index, "transcript_id_version")),
            "Transcript mapping evidence is conflicting.",
        )
    if candidate.protein_accession_version.is_conflicting or ref.protein_accession_version.is_conflicting:
        return _Assessment(
            "UNDETERMINED", tuple(used), ("candidate_protein_accession_version", label(reference_index, "protein_accession_version")),
            "Protein reference-sequence evidence is conflicting.",
        )
    if candidate.transcript_id_version.is_known and ref.transcript_id_version.is_known:
        used.extend(("candidate_transcript_id_version", label(reference_index, "transcript_id_version")))
    if candidate.protein_accession_version.is_known and ref.protein_accession_version.is_known:
        used.extend(("candidate_protein_accession_version", label(reference_index, "protein_accession_version")))

    if explicit.is_true:
        # `protein_substitution_equivalence` is a derived assertion, not a
        # replacement for the normalized protein observations from which it
        # must have been derived.  It can resolve reference-sequence mapping
        # differences only after a same protein effect is independently shown.
        if not protein_effect_known_match:
            return _Assessment(
                "UNDETERMINED", tuple(used) + (explicit_label,),
                ("candidate_protein_hgvs", label(reference_index, "protein_hgvs")),
                "A derived protein-equivalence assertion cannot be accepted without normalized amino-acid effect annotations.",
            )
        return _Assessment("OK", tuple(used) + (explicit_label,))

    if transcript_match is False or accession_match is False:
        return _Assessment(
            "UNDETERMINED", tuple(used),
            (explicit_label,),
            "Protein HGVS equality cannot establish equivalence across different transcript or protein reference sequences.",
        )
    if protein_effect_known_match and transcript_match is True and accession_match is True:
        return _Assessment("OK", tuple(used))

    if not protein_effect_known_match:
        missing.extend(("candidate_protein_hgvs", label(reference_index, "protein_hgvs")))
    if transcript_match is None:
        missing.extend(("candidate_transcript_id_version", label(reference_index, "transcript_id_version")))
    if accession_match is None:
        missing.extend(("candidate_protein_accession_version", label(reference_index, "protein_accession_version")))
    missing.append(explicit_label)
    return _Assessment(
        "UNDETERMINED", tuple(used), tuple(dict.fromkeys(missing)),
        "Validated protein-substitution equivalence has not been established.",
    )


def _variant_has_no_splice_mechanism(
    annotations: VariantAnnotations,
    prefix: str,
    *,
    reference_index: int,
) -> _Assessment:
    """Require a documented absence of a predicted/observed splice mechanism."""

    splice = annotations.splice
    used: list[str] = []
    label_for = lambda name: f"{prefix}{name}" if prefix == "candidate_" else label(reference_index, name)
    impact_name = label_for("splice_impact")
    prediction_name = label_for("splice_prediction")
    rna_name = label_for("rna_splicing_evidence")
    product_name = label_for("observed_splice_product")
    residual_name = label_for("residual_normal_splicing")

    if (
        splice.impact.is_conflicting
        or splice.prediction.is_conflicting
        or splice.rna_splicing_evidence.is_conflicting
        or splice.residual_normal_splicing.is_conflicting
    ):
        field = (
            impact_name
            if splice.impact.is_conflicting
            else prediction_name
            if splice.prediction.is_conflicting
            else rna_name
            if splice.rna_splicing_evidence.is_conflicting
            else residual_name
        )
        return _Assessment("UNDETERMINED", tuple(used), (field,), "Splice-mechanism evidence is conflicting.")
    if splice.observed_splice_product.is_conflicting:
        return _Assessment(
            "UNDETERMINED", tuple(used), (product_name,),
            "Observed RNA product evidence is conflicting.",
        )
    if splice.rna_splicing_evidence.is_false and (
        splice.observed_splice_product.is_known or splice.residual_normal_splicing.is_known
    ):
        supplied_name = (
            product_name if splice.observed_splice_product.is_known else residual_name
        )
        return _Assessment(
            "UNDETERMINED", tuple(used) + (rna_name, supplied_name),
            (rna_name, supplied_name),
            "RNA-splicing evidence is marked false despite a supplied RNA result.",
        )

    abnormal = False
    if splice.impact.is_known:
        used.append(impact_name)
        if splice.impact.value in {SpliceImpact.POSSIBLE, SpliceImpact.ABERRANT}:
            abnormal = True
    if splice.prediction.is_known:
        used.append(prediction_name)
        if splice.prediction.value == SplicePrediction.SPLICE_ALTERING:
            abnormal = True
        if splice.prediction.value == SplicePrediction.UNCERTAIN and not splice.impact.is_known:
            return _Assessment("UNDETERMINED", tuple(used), (prediction_name,), "Splice prediction is explicitly uncertain.")
    if splice.rna_splicing_evidence.is_known:
        used.append(rna_name)
    if splice.rna_splicing_evidence.is_true and splice.observed_splice_product.is_unknown:
        return _Assessment("UNDETERMINED", tuple(used), (product_name,), "RNA splicing evidence was supplied without its observed product.")
    if splice.observed_splice_product.is_known:
        used.append(product_name)
        product = normalized(splice.observed_splice_product.value or "")
        if product not in {"NORMAL", "NO_ABERRANT_SPLICING", "NO_IMPACT"}:
            abnormal = True
    if splice.residual_normal_splicing.is_known:
        used.append(residual_name)
        if splice.residual_normal_splicing.value is False:
            abnormal = True
        elif (
            splice.rna_splicing_evidence.is_unknown
            and splice.observed_splice_product.is_unknown
        ):
            return _Assessment(
                "UNDETERMINED", tuple(used), (product_name,),
                "Residual-normal-splicing evidence was supplied without a complete RNA splice-product result.",
            )
    if abnormal:
        return _Assessment(
            "NOT_MET", tuple(used), (),
            "A predicted or observed splice-mediated effect prevents use of the standard protein-equivalence PS1 path.",
            (exclusion(impact_name, "SPLICE_MEDIATED_OR_POSSIBLE"),),
        )

    if (splice.impact.is_known and splice.impact.value == SpliceImpact.NO_IMPACT) or (
        splice.prediction.is_known and splice.prediction.value == SplicePrediction.NO_IMPACT
    ):
        return _Assessment("OK", tuple(used))
    return _Assessment(
        "UNDETERMINED", tuple(used), (impact_name, prediction_name),
        "A no-splice-impact conclusion is required before standard amino-acid PS1 can be applied.",
    )


def _splice_and_mechanism_compatibility(
    evidence: BiologicalFeatures, reference: ReferenceVariant, reference_index: int
) -> _Assessment:
    explicit, explicit_label = pairwise_evidence(
        evidence, reference, reference_index, "protein_effect_mechanism_compatibility"
    )
    if explicit.is_conflicting:
        return _Assessment(
            "UNDETERMINED", (explicit_label,),
            (explicit_label,), "Protein-mechanism compatibility is conflicting."
        )
    if explicit.is_false:
        return _Assessment(
            "NOT_MET", (explicit_label,), (),
            "The supplied evidence establishes incompatible molecular mechanisms.",
            (exclusion(explicit_label, "INCOMPATIBLE"),),
        )

    reference_mechanism = reference.molecular_mechanism
    reference_mechanism_label = label(reference_index, "molecular_mechanism")
    if reference_mechanism.is_conflicting:
        return _Assessment("UNDETERMINED", (), (reference_mechanism_label,), "Reference molecular mechanism is conflicting.")
    if reference_mechanism.is_known and reference_mechanism.value != MolecularMechanism.PROTEIN_ALTERATION:
        return _Assessment(
            "NOT_MET", (reference_mechanism_label,), (),
            "The reference is established through a non-protein-alteration mechanism.",
            (exclusion(reference_mechanism_label, "INCOMPATIBLE_MECHANISM"),),
        )

    if reference_mechanism.is_unknown and not explicit.is_true:
        return _Assessment(
            "UNDETERMINED", (),
            (reference_mechanism_label, explicit_label),
            "Protein-effect mechanism compatibility has not been established for the reference.",
        )

    candidate_assessment = _variant_has_no_splice_mechanism(
        evidence.candidate, "candidate_", reference_index=reference_index
    )
    reference_assessment = _variant_has_no_splice_mechanism(
        reference.annotations, "reference_", reference_index=reference_index
    )
    used = candidate_assessment.used + reference_assessment.used
    if candidate_assessment.state == "NOT_MET" or reference_assessment.state == "NOT_MET":
        culprit = candidate_assessment if candidate_assessment.state == "NOT_MET" else reference_assessment
        return _Assessment("NOT_MET", used, (), culprit.reason, culprit.exclusions)
    if candidate_assessment.state == "UNDETERMINED" or reference_assessment.state == "UNDETERMINED":
        missing = candidate_assessment.missing + reference_assessment.missing
        # A curated direct compatibility conclusion can satisfy this check only
        # when no contrary splice observation was found.
        if explicit.is_true:
            return _Assessment("OK", used + (explicit_label,))
        return _Assessment(
            "UNDETERMINED", used, tuple(dict.fromkeys(missing)),
            "Splice-related mechanism compatibility has not been established.",
        )
    return _Assessment("OK", used + ((explicit_label,) if explicit.is_true else ()))


def evaluate_standard_amino_acid_path(
    evidence: BiologicalFeatures, reference: ReferenceVariant, reference_index: int
) -> RuleEvaluation:
    """Evaluate the non-splicing ACMG/AMP 2015 PS1 pathway for one reference."""

    candidate_consequence = evidence.candidate.consequence
    if candidate_consequence.is_known and candidate_consequence.value in {
        VariantConsequence.SYNONYMOUS,
        VariantConsequence.INTRONIC,
    }:
        return not_applicable(
            PATH, ACMG_PS1,
            "The supplied candidate consequence is not a protein-altering missense pathway.",
            reference_index=reference_index, used=("candidate_variant_consequence",),
        )

    preconditions: ReferencePreconditions = common_reference_preconditions(
        evidence, reference, reference_index, path=PATH, guideline=ACMG_PS1
    )
    if preconditions.terminal is not None:
        return preconditions.terminal

    quality = require_reference_evidence_quality(
        evidence, reference, reference_index, path=PATH, guideline=ACMG_PS1,
        used=preconditions.used, exclusions=preconditions.exclusions,
    )
    if quality is not None:
        return quality

    classification = reference.classification
    classification_label = label(reference_index, "classification")
    # The generic path explicitly requires a curated established reference
    # record.  Preserve that feature in the final audit/provenance trail on
    # the successful branch as well as on quality failures.
    used = preconditions.used + (
        label(reference_index, "classification_evidence"),
        classification_label,
    )
    if classification.value != ReferenceClassification.PATHOGENIC:
        return not_met(
            PATH, ACMG_PS1,
            "Generic ACMG/AMP PS1 requires a previously established Pathogenic, not merely Likely Pathogenic, reference.",
            reference_index=reference_index, used=used,
            exclusions=preconditions.exclusions + (exclusion(classification_label, "NOT_PATHOGENIC_FOR_GENERIC_PS1"),),
        )

    protein = _protein_equivalence(evidence, reference, reference_index)
    used += protein.used
    if protein.state == "NOT_MET":
        return not_met(PATH, ACMG_PS1, protein.reason, reference_index=reference_index, used=used,
            exclusions=preconditions.exclusions + protein.exclusions)
    if protein.state == "UNDETERMINED":
        return undetermined(PATH, ACMG_PS1, protein.reason, protein.missing, reference_index=reference_index,
            used=used, exclusions=preconditions.exclusions + protein.exclusions)

    mechanism = _splice_and_mechanism_compatibility(evidence, reference, reference_index)
    used += mechanism.used
    if mechanism.state == "NOT_MET":
        return not_met(PATH, ACMG_PS1, mechanism.reason, reference_index=reference_index, used=used,
            exclusions=preconditions.exclusions + mechanism.exclusions)
    if mechanism.state == "UNDETERMINED":
        return undetermined(PATH, ACMG_PS1, mechanism.reason, mechanism.missing, reference_index=reference_index,
            used=used, exclusions=preconditions.exclusions + mechanism.exclusions)

    return activated(
        PS1Strength.STRONG, PATH, ACMG_PS1,
        "The candidate and an independently established Pathogenic reference have the same validated amino-acid substitution from different alleles, with no predicted or observed splice-mediated mechanism.",
        reference_index=reference_index, used=used, exclusions=preconditions.exclusions,
    )

