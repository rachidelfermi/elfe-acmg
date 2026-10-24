"""ClinGen SVI (Walker et al. 2023) splicing-related PS1 pathway."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from ..guidelines import CLINGEN_SPLICE_PS1
from ..models import (
    BiologicalFeatures,
    DecisionStatus,
    MolecularMechanism,
    PS1Strength,
    ReferenceClassification,
    ReferenceVariant,
    RnaAssayReliability,
    SpliceBaselineCode,
    SpliceImpact,
    SpliceMotifPosition,
    SpliceProductFunctionalConsequence,
    SplicePrediction,
    SpliceSiteType,
    VariantConsequence,
)
from ..output import RuleEvaluation
from .common import (
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


PATH = "SPLICING_SAME_PREDICTED_EVENT"


@dataclass(frozen=True)
class _Assessment:
    state: str  # OK, NOT_MET, UNDETERMINED, NOT_APPLICABLE
    used: tuple[str, ...] = ()
    missing: tuple[str, ...] = ()
    reason: str = ""
    exclusions: tuple[dict[str, str], ...] = ()
    same_nucleotide: Optional[bool] = None


def _candidate_has_splice_context(evidence: BiologicalFeatures) -> bool:
    candidate = evidence.candidate
    consequence = candidate.consequence
    if consequence.is_known and consequence.value in {
        VariantConsequence.SPLICE_REGION,
        VariantConsequence.SPLICE_DONOR,
        VariantConsequence.SPLICE_ACCEPTOR,
    }:
        return True
    splice = candidate.splice
    return any((
        (
            splice.impact.is_known and splice.impact.value is not None and splice.impact.value.value != "NO_IMPACT"
        ),
        (
            splice.prediction.is_known and splice.prediction.value is not None and splice.prediction.value.value != "NO_IMPACT"
        ),
        splice.site_type.is_known,
        splice.motif_position.is_known,
        splice.motif_id.is_known,
        splice.nucleotide_position.is_known,
        splice.spliceai_event_type.is_known,
        splice.spliceai_event_position.is_known,
        splice.rna_splicing_evidence.is_true,
        splice.observed_splice_product.is_known,
    ))


def _reference_splice_mechanism(
    reference: ReferenceVariant, reference_index: int
) -> _Assessment:
    """Require the reference pathogenicity to be attributable to splicing."""

    field = label(reference_index, "molecular_mechanism")
    mechanism = reference.molecular_mechanism
    if mechanism.is_conflicting:
        return _Assessment(
            "UNDETERMINED", (), (field,),
            "The reference molecular mechanism is conflicting.",
        )
    if not mechanism.is_known:
        return _Assessment(
            "UNDETERMINED", (), (field,),
            "The reference molecular mechanism is not established for a splice PS1 comparison.",
        )
    if mechanism.value != MolecularMechanism.SPLICING:
        return _Assessment(
            "NOT_MET", (field,), (),
            "The reference pathogenic mechanism is not splicing.",
            (exclusion(field, "INCOMPATIBLE_SPLICE_MECHANISM"),),
        )
    return _Assessment("OK", (field,))


def _transcript_compatibility(
    evidence: BiologicalFeatures, reference: ReferenceVariant, reference_index: int
) -> _Assessment:
    """Require a common curation transcript for a splice-event comparison.

    Walker et al. position and motif comparisons are relative to the curation
    transcript.  A caller may provide a curated mapping assertion when two
    transcript labels legitimately describe the same compared event.
    """

    explicit, explicit_label = pairwise_evidence(
        evidence, reference, reference_index, "splice_transcript_equivalence"
    )
    candidate = evidence.candidate.transcript_id_version
    ref = reference.annotations.transcript_id_version
    candidate_label = "candidate_transcript_id_version"
    reference_label = label(reference_index, "transcript_id_version")
    used: list[str] = []
    if explicit.is_conflicting:
        return _Assessment(
            "UNDETERMINED", (explicit_label,), (explicit_label,),
            "Splice-transcript equivalence is conflicting.",
        )
    if candidate.is_conflicting or ref.is_conflicting:
        field = candidate_label if candidate.is_conflicting else reference_label
        return _Assessment(
            "UNDETERMINED", (), (field,),
            "Splice transcript mapping evidence is conflicting.",
        )
    transcript_match: Optional[bool] = None
    if candidate.is_known and ref.is_known:
        used.extend((candidate_label, reference_label))
        transcript_match = normalized(candidate.value or "") == normalized(ref.value or "")
    if transcript_match is True:
        if explicit.is_false:
            return _Assessment(
                "UNDETERMINED", tuple(used) + (explicit_label,), (explicit_label,),
                "The splice-transcript equivalence assertion contradicts identical transcript labels.",
            )
        return _Assessment(
            "OK", tuple(used) + ((explicit_label,) if explicit.is_true else ())
        )
    if explicit.is_false:
        return _Assessment(
            "NOT_MET", tuple(used) + (explicit_label,), (),
            "The supplied evidence establishes that the splice comparisons are not on equivalent transcript mappings.",
            (exclusion(explicit_label, "DIFFERENT_TRANSCRIPT_MAPPING"),),
        )
    if explicit.is_true:
        return _Assessment("OK", tuple(used) + (explicit_label,))
    return _Assessment(
        "UNDETERMINED", tuple(used),
        (candidate_label, reference_label, explicit_label),
        "The transcript mapping for the splice-event comparison has not been established.",
    )


def _exonic_protein_effect_consistency(
    evidence: BiologicalFeatures, reference: ReferenceVariant, reference_index: int
) -> _Assessment:
    """Account for encoded missense effects when a splice comparison is exonic.

    The generic Table 2 route can be driven by the same predicted splice
    event, but Walker et al. directs curators to consider an exonic variant's
    possible encoded amino-acid effect.  This gate does not invent a new
    strength: it prevents activation while a known additional protein effect
    is incompatible or has not been reconciled.
    """

    candidate = evidence.candidate
    ref = reference.annotations
    compatibility, compatibility_label = pairwise_evidence(
        evidence, reference, reference_index, "protein_effect_mechanism_compatibility"
    )
    aa_equivalence, aa_equivalence_label = pairwise_evidence(
        evidence, reference, reference_index, "protein_substitution_equivalence"
    )
    if compatibility.is_conflicting or aa_equivalence.is_conflicting:
        field = compatibility_label if compatibility.is_conflicting else aa_equivalence_label
        return _Assessment(
            "UNDETERMINED", (field,), (field,),
            "Exonic protein-effect compatibility evidence is conflicting.",
        )
    if compatibility.is_false or aa_equivalence.is_false:
        field = compatibility_label if compatibility.is_false else aa_equivalence_label
        return _Assessment(
            "NOT_MET", (field,), (),
            "The supplied evidence establishes an incompatible encoded protein effect for this splice comparison.",
            (exclusion(field, "INCOMPATIBLE_EXONIC_PROTEIN_EFFECT"),),
        )

    candidate_consequence = candidate.consequence
    reference_consequence = ref.consequence
    if candidate_consequence.is_conflicting or reference_consequence.is_conflicting:
        field = "candidate_variant_consequence" if candidate_consequence.is_conflicting else label(reference_index, "variant_consequence")
        return _Assessment(
            "UNDETERMINED", (), (field,),
            "Exonic consequence annotation is conflicting.",
        )
    exonic = (
        (candidate_consequence.is_known and candidate_consequence.value == VariantConsequence.MISSENSE)
        or (reference_consequence.is_known and reference_consequence.value == VariantConsequence.MISSENSE)
    )
    pairs = (
        ("candidate_protein_hgvs", candidate.protein_hgvs, label(reference_index, "protein_hgvs"), ref.protein_hgvs),
        ("candidate_amino_acid_position", candidate.amino_acid_position, label(reference_index, "amino_acid_position"), ref.amino_acid_position),
        ("candidate_reference_amino_acid", candidate.reference_amino_acid, label(reference_index, "reference_amino_acid"), ref.reference_amino_acid),
        ("candidate_alternate_amino_acid", candidate.alternate_amino_acid, label(reference_index, "alternate_amino_acid"), ref.alternate_amino_acid),
    )
    used: list[str] = []
    known_mismatch = False
    protein_match = aa_equivalence.is_true
    component_pairs_known = True
    any_protein_annotation = False
    for candidate_name, candidate_value, reference_name, reference_value in pairs:
        if candidate_value.is_conflicting or reference_value.is_conflicting:
            field = candidate_name if candidate_value.is_conflicting else reference_name
            return _Assessment(
                "UNDETERMINED", tuple(used), (field,),
                "Encoded protein-effect annotation is conflicting.",
            )
        if candidate_value.is_known or reference_value.is_known:
            any_protein_annotation = True
        if candidate_value.is_known and reference_value.is_known:
            used.extend((candidate_name, reference_name))
            left, right = candidate_value.value, reference_value.value
            equal = normalized(left) == normalized(right) if isinstance(left, str) and isinstance(right, str) else left == right
            if not equal:
                known_mismatch = True
            if candidate_name == "candidate_protein_hgvs" and equal:
                protein_match = True
        elif candidate_name != "candidate_protein_hgvs":
            component_pairs_known = False
    if component_pairs_known and any_protein_annotation and not known_mismatch:
        protein_match = True

    if not exonic and not any_protein_annotation:
        return _Assessment("OK")
    if known_mismatch:
        if compatibility.is_true or aa_equivalence.is_true:
            return _Assessment(
                "UNDETERMINED", tuple(used) + (compatibility_label,), (compatibility_label,),
                "A declared protein-effect compatibility conflicts with distinct encoded amino-acid annotations.",
            )
        return _Assessment(
            "UNDETERMINED", tuple(used), (compatibility_label,),
            "Distinct encoded amino-acid effects require a curated compatibility conclusion before generic splice PS1 can be applied.",
        )
    if protein_match or compatibility.is_true:
        extra = ()
        if aa_equivalence.is_true:
            extra += (aa_equivalence_label,)
        if compatibility.is_true:
            extra += (compatibility_label,)
        return _Assessment("OK", tuple(used) + extra)
    return _Assessment(
        "UNDETERMINED", tuple(used), (compatibility_label,),
        "A possible encoded missense effect has not been reconciled for this splice comparison.",
    )


def _functional_consequence_consistency(evidence: BiologicalFeatures) -> _Assessment:
    """Detect an explicit conflict between an LoF baseline and splice product."""

    consequence = evidence.splice_product_functional_consequence
    baseline = evidence.applicable_splicing_pvs1_strength
    if consequence.is_conflicting:
        return _Assessment(
            "UNDETERMINED", (), ("splice_product_functional_consequence",),
            "Splice-product functional consequence is conflicting.",
        )
    if (
        consequence.is_known
        and consequence.value == SpliceProductFunctionalConsequence.NOT_LOSS_OF_FUNCTION
        and baseline.is_known
        and baseline.value
        in {
            SpliceBaselineCode.PVS1,
            SpliceBaselineCode.PVS1_STRONG,
            SpliceBaselineCode.PVS1_MODERATE,
            SpliceBaselineCode.PVS1_SUPPORTING,
        }
    ):
        return _Assessment(
            "UNDETERMINED", ("splice_product_functional_consequence", "applicable_splicing_pvs1_strength"),
            ("splice_product_functional_consequence",),
            "The observed or curated non-loss-of-function product conflicts with the supplied PVS1 baseline.",
        )
    return _Assessment(
        "OK", ("splice_product_functional_consequence",) if consequence.is_known else ()
    )


def _events_and_predictions_equivalent(
    evidence: BiologicalFeatures,
    reference: ReferenceVariant,
    reference_index: int,
    *,
    require_prediction_strength_comparable: bool = True,
) -> _Assessment:
    candidate = evidence.candidate.splice
    ref = reference.annotations.splice
    used: list[str] = []
    missing: list[str] = []

    candidate_prediction_name = "candidate_splice_prediction"
    reference_prediction_name = label(reference_index, "splice_prediction")
    for value, field in ((candidate.prediction, candidate_prediction_name), (ref.prediction, reference_prediction_name)):
        if value.is_conflicting:
            return _Assessment("UNDETERMINED", tuple(used), (field,), "Splice prediction evidence is conflicting.")
        if value.is_unknown:
            missing.append(field)
        else:
            used.append(field)
            if value.value != SplicePrediction.SPLICE_ALTERING:
                if value.value == SplicePrediction.UNCERTAIN:
                    return _Assessment(
                        "UNDETERMINED", tuple(used), (field,),
                        "A splice-PS1 comparison requires an established splice-altering prediction; the supplied prediction is uncertain.",
                    )
                return _Assessment(
                    "NOT_MET", tuple(used), (),
                    "Walker et al. Table 2 requires a splice-altering predicted event for both variants.",
                    (exclusion(field, "NO_COMPARABLE_SPLICE_ALTERING_PREDICTION"),),
                )
    if missing:
        return _Assessment(
            "UNDETERMINED", tuple(used), tuple(missing),
            "A splice-PS1 comparison requires predicted splice events for both variants.",
        )

    # A curated/observed no-impact conclusion and a splice-altering prediction
    # are material disagreement, not evidence that may be silently ignored.
    for value, field in (
        (candidate.impact, "candidate_splice_impact"),
        (ref.impact, label(reference_index, "splice_impact")),
    ):
        if value.is_conflicting:
            return _Assessment(
                "UNDETERMINED", tuple(used), (field,),
                "Splice-impact evidence is conflicting.",
            )
        if value.is_known:
            used.append(field)
            if value.value == SpliceImpact.NO_IMPACT:
                return _Assessment(
                    "UNDETERMINED", tuple(used), (field,),
                    "A no-impact splice conclusion conflicts with the supplied splice-altering prediction.",
                )

    explicit, explicit_label = pairwise_evidence(
        evidence, reference, reference_index, "splice_event_equivalence"
    )
    if explicit.is_conflicting:
        return _Assessment("UNDETERMINED", tuple(used), (explicit_label,), "Splice-event equivalence is conflicting.")
    if explicit.is_false:
        return _Assessment(
            "NOT_MET", tuple(used) + (explicit_label,), (),
            "The supplied evidence establishes non-equivalent splice events.",
            (exclusion(explicit_label, "DIFFERENT_EVENTS"),),
        )

    event_pairs = (
        ("candidate_spliceai_event_type", candidate.spliceai_event_type, label(reference_index, "spliceai_event_type"), ref.spliceai_event_type),
        ("candidate_spliceai_event_position", candidate.spliceai_event_position, label(reference_index, "spliceai_event_position"), ref.spliceai_event_position),
    )
    raw_match: Optional[bool] = True
    for candidate_name, candidate_value, reference_name, reference_value in event_pairs:
        if candidate_value.is_conflicting or reference_value.is_conflicting:
            field = candidate_name if candidate_value.is_conflicting else reference_name
            return _Assessment("UNDETERMINED", tuple(used), (field,), "Predicted splice-event annotation is conflicting.")
        if not candidate_value.is_known or not reference_value.is_known:
            raw_match = None
            continue
        used.extend((candidate_name, reference_name))
        if normalized(candidate_value.value or "") != normalized(reference_value.value or ""):
            raw_match = False

    if raw_match is False:
        if explicit.is_true:
            return _Assessment(
                "UNDETERMINED", tuple(used) + (explicit_label,),
                (explicit_label,),
                "The declared splice-event equivalence conflicts with event type or position annotations.",
            )
        return _Assessment(
            "NOT_MET", tuple(used), (), "The predicted splice events are not precisely equivalent.",
            (exclusion(explicit_label, "DIFFERENT_EVENTS"),),
        )
    if raw_match is None and not explicit.is_true:
        return _Assessment(
            "UNDETERMINED", tuple(used),
            ("candidate_spliceai_event_type", label(reference_index, "spliceai_event_type"),
             "candidate_spliceai_event_position", label(reference_index, "spliceai_event_position"),
             explicit_label),
            "Precise equivalence of the predicted splice event has not been established.",
        )
    if explicit.is_true:
        used.append(explicit_label)

    if not require_prediction_strength_comparable:
        # A pinned VCEP may document and apply a gene-specific numerical
        # predictor threshold.  The caller of this helper has already verified
        # that such a condition matched; generic Walker comparison strength is
        # therefore not silently converted into a universal cutoff here.
        return _Assessment("OK", tuple(used))

    strength, strength_label = pairwise_evidence(
        evidence, reference, reference_index, "splice_prediction_strength_comparable"
    )
    if strength.is_conflicting:
        return _Assessment("UNDETERMINED", tuple(used), (strength_label,), "Prediction-strength comparison is conflicting.")
    if strength.is_unknown:
        return _Assessment(
            "UNDETERMINED", tuple(used), (strength_label,),
            "Walker et al. Table 2 requires a curated finding that the candidate prediction is similar to or stronger than the reference prediction; raw scores are not assigned an invented universal threshold.",
        )
    used.append(strength_label)
    if strength.value is False:
        return _Assessment(
            "NOT_MET", tuple(used), (), "Candidate splice prediction is weaker than the comparison prediction.",
            (exclusion(strength_label, "WEAKER_THAN_REFERENCE"),),
        )
    return _Assessment("OK", tuple(used))


def _motif_relationship(
    evidence: BiologicalFeatures, reference: ReferenceVariant, reference_index: int
) -> _Assessment:
    candidate = evidence.candidate.splice
    ref = reference.annotations.splice
    used: list[str] = []
    required = (
        ("candidate_splice_site_type", candidate.site_type),
        ("candidate_splice_motif_position", candidate.motif_position),
        (label(reference_index, "splice_site_type"), ref.site_type),
        (label(reference_index, "splice_motif_position"), ref.motif_position),
    )
    for field, value in required:
        if value.is_conflicting:
            return _Assessment("UNDETERMINED", tuple(used), (field,), "Splice-motif location evidence is conflicting.")
        if value.is_unknown:
            return _Assessment("UNDETERMINED", tuple(used), (field,), "Splice-motif location is required by Walker et al. Table 2.")
        used.append(field)
    if candidate.site_type.value not in {SpliceSiteType.DONOR, SpliceSiteType.ACCEPTOR}:
        return _Assessment("NOT_APPLICABLE", tuple(used), (), "The candidate is not located in a donor or acceptor motif covered by the generic Walker Table 2 pathway.")
    if ref.site_type.value != candidate.site_type.value:
        return _Assessment(
            "NOT_MET", tuple(used), (), "Candidate and reference do not affect the same type of splice motif.",
            (exclusion(label(reference_index, "splice_site_type"), "DIFFERENT_MOTIF"),),
        )
    if candidate.motif_position.value == SpliceMotifPosition.OUTSIDE_MOTIF or ref.motif_position.value == SpliceMotifPosition.OUTSIDE_MOTIF:
        return _Assessment(
            "NOT_APPLICABLE", tuple(used), (), "At least one variant is outside the donor/acceptor motif range covered by generic splice-PS1 guidance.",
        )

    candidate_position = candidate.nucleotide_position
    reference_position = ref.nucleotide_position
    if candidate_position.is_conflicting or reference_position.is_conflicting:
        field = "candidate_splice_nucleotide_position" if candidate_position.is_conflicting else label(reference_index, "splice_nucleotide_position")
        return _Assessment("UNDETERMINED", tuple(used), (field,), "Splice nucleotide position is conflicting.")
    if candidate_position.is_known and reference_position.is_known:
        used.extend(("candidate_splice_nucleotide_position", label(reference_index, "splice_nucleotide_position")))
        same_nucleotide = normalized(candidate_position.value or "") == normalized(reference_position.value or "")
    else:
        same_nucleotide = None

    explicit, explicit_label = pairwise_evidence(
        evidence, reference, reference_index, "splice_motif_equivalence"
    )
    if explicit.is_conflicting:
        return _Assessment("UNDETERMINED", tuple(used), (explicit_label,), "Splice-motif equivalence is conflicting.")
    if explicit.is_false:
        if same_nucleotide is True:
            return _Assessment(
                "UNDETERMINED", tuple(used) + (explicit_label,), (explicit_label,),
                "A different-motif assertion conflicts with identical splice-nucleotide positions.",
            )
        return _Assessment(
            "NOT_MET", tuple(used) + (explicit_label,), (), "The variants are established to affect different splice motifs.",
            (exclusion(explicit_label, "DIFFERENT_MOTIF"),),
        )

    candidate_id = candidate.motif_id
    reference_id = ref.motif_id
    raw_match: Optional[bool] = None
    if candidate_id.is_conflicting or reference_id.is_conflicting:
        field = "candidate_splice_motif_id" if candidate_id.is_conflicting else label(reference_index, "splice_motif_id")
        return _Assessment("UNDETERMINED", tuple(used), (field,), "Splice motif identifier is conflicting.")
    if candidate_id.is_known and reference_id.is_known:
        used.extend(("candidate_splice_motif_id", label(reference_index, "splice_motif_id")))
        raw_match = normalized(candidate_id.value or "") == normalized(reference_id.value or "")
    if raw_match is False:
        if explicit.is_true or same_nucleotide is True:
            return _Assessment(
                "UNDETERMINED", tuple(used) + ((explicit_label,) if explicit.is_true else ()),
                (explicit_label,),
                "Splice-motif identifiers conflict with the supplied motif relationship.",
            )
        return _Assessment(
            "NOT_MET", tuple(used), (), "Candidate and reference motif identifiers differ.",
            (exclusion(explicit_label, "DIFFERENT_MOTIF"),),
        )

    # For the Table 2 same-nucleotide row, a common curation transcript,
    # identical nucleotide position, and matching donor/acceptor context are
    # already sufficient to establish the motif.  Arbitrary motif identifiers
    # are required only for comparisons at different or unknown coordinates.
    if same_nucleotide is True:
        if explicit.is_true:
            used.append(explicit_label)
        return _Assessment("OK", tuple(used), same_nucleotide=True)

    if raw_match is None and not explicit.is_true:
        return _Assessment(
            "UNDETERMINED", tuple(used), ("candidate_splice_motif_id", label(reference_index, "splice_motif_id"), explicit_label),
            "The same donor or acceptor motif has not been established.",
        )
    if explicit.is_true:
        used.append(explicit_label)
    return _Assessment("OK", tuple(used), same_nucleotide=same_nucleotide)


def _rna_consistency(
    evidence: BiologicalFeatures, reference: ReferenceVariant, reference_index: int
) -> _Assessment:
    """Treat available reliable RNA data as a potential material contradiction.

    RNA experiments are not required for the generic predicted-event route.
    Once an RNA result is supplied as reliable, however, an incompatible result
    cannot be silently ignored.
    """

    reliability = evidence.splice_rna_assay_reliability
    candidate = evidence.candidate.splice
    ref = reference.annotations.splice
    _, event_label = pairwise_evidence(
        evidence, reference, reference_index, "splice_event_equivalence"
    )
    used: list[str] = []
    if reliability.is_conflicting:
        return _Assessment("UNDETERMINED", (), ("splice_rna_assay_reliability",), "RNA assay reliability is conflicting.")
    if reliability.is_known:
        used.append("splice_rna_assay_reliability")

    rna_rows = (
        (
            "candidate_rna_splicing_evidence",
            candidate.rna_splicing_evidence,
            "candidate_observed_splice_product",
            candidate.observed_splice_product,
            "candidate_residual_normal_splicing",
            candidate.residual_normal_splicing,
        ),
        (
            label(reference_index, "rna_splicing_evidence"),
            ref.rna_splicing_evidence,
            label(reference_index, "observed_splice_product"),
            ref.observed_splice_product,
            label(reference_index, "residual_normal_splicing"),
            ref.residual_normal_splicing,
        ),
    )
    for rna_name, rna_evidence, product_name, product, residual_name, residual in rna_rows:
        if rna_evidence.is_conflicting:
            return _Assessment("UNDETERMINED", tuple(used), (rna_name,), "RNA splicing evidence is conflicting.")
        if product.is_conflicting:
            return _Assessment("UNDETERMINED", tuple(used), (product_name,), "Observed splice-product evidence is conflicting.")
        if residual.is_conflicting:
            return _Assessment("UNDETERMINED", tuple(used), (residual_name,), "Residual-normal-splicing evidence is conflicting.")
        if rna_evidence.is_false and (product.is_known or residual.is_known):
            supplied = (product_name,) if product.is_known else (residual_name,)
            return _Assessment(
                "UNDETERMINED", tuple(used) + (rna_name,) + supplied,
                (rna_name,) + supplied,
                "RNA-splicing evidence is marked false despite a supplied RNA result.",
            )

    if reliability.is_known and reliability.value == RnaAssayReliability.INADEQUATE:
        # An explicitly inadequate assay cannot overrule a complete
        # prediction-based Table 2 path, but contradictory RNA input was
        # rejected above rather than being silently ignored.
        return _Assessment("OK", tuple(used))

    candidate_residual = candidate.residual_normal_splicing
    reference_residual = ref.residual_normal_splicing
    if candidate_residual.is_known or reference_residual.is_known:
        if candidate_residual.is_known:
            used.append("candidate_residual_normal_splicing")
        if reference_residual.is_known:
            used.append(label(reference_index, "residual_normal_splicing"))
        if not candidate_residual.is_known or not reference_residual.is_known:
            missing = (
                "candidate_residual_normal_splicing"
                if not candidate_residual.is_known
                else label(reference_index, "residual_normal_splicing")
            )
            return _Assessment(
                "UNDETERMINED", tuple(used), (missing,),
                "A residual-normal-splicing profile was supplied for only one comparison variant.",
            )
        if candidate_residual.value != reference_residual.value:
            if reliability.is_known and reliability.value == RnaAssayReliability.ESTABLISHED:
                return _Assessment(
                    "NOT_MET", tuple(used), (), "Reliable RNA evidence shows different residual normal-splicing profiles.",
                    (exclusion(event_label, "RESIDUAL_SPLICING_DIFFERENCE"),),
                )
            return _Assessment(
                "UNDETERMINED", tuple(used), ("splice_rna_assay_reliability",),
                "Residual normal-splicing profiles differ, but their reliability is not established.",
            )

    candidate_has_rna = (
        candidate.rna_splicing_evidence.is_true or candidate.observed_splice_product.is_known
    )
    reference_has_rna = (
        ref.rna_splicing_evidence.is_true or ref.observed_splice_product.is_known
    )
    if not candidate_has_rna and not reference_has_rna:
        return _Assessment("OK", tuple(used))

    candidate_product = candidate.observed_splice_product
    reference_product = ref.observed_splice_product
    if candidate.rna_splicing_evidence.is_true:
        used.append("candidate_rna_splicing_evidence")
    if ref.rna_splicing_evidence.is_true:
        used.append(label(reference_index, "rna_splicing_evidence"))
    product_rows = (
        ("candidate_observed_splice_product", candidate_product, candidate_has_rna),
        (label(reference_index, "observed_splice_product"), reference_product, reference_has_rna),
    )
    for field, value, has_rna in product_rows:
        if not has_rna:
            continue
        if value.is_unknown:
            return _Assessment("UNDETERMINED", tuple(used), (field,), "RNA evidence was supplied without its observed splice product.")
        used.append(field)

    # Both variants were required above to have splice-altering predictions.
    # A reliable normal RNA product for *either* variant directly contradicts
    # that premise, even if the other variant has no RNA experiment.
    normal_products = {"NORMAL", "NO_ABERRANT_SPLICING", "NO_IMPACT"}
    for field, value, has_rna in product_rows:
        if has_rna and normalized(value.value or "") in normal_products:
            if reliability.is_known and reliability.value == RnaAssayReliability.ESTABLISHED:
                return _Assessment(
                    "NOT_MET", tuple(used), (),
                    "Reliable RNA evidence shows normal splicing despite the splice-altering prediction.",
                    (exclusion(field, "RNA_CONTRADICTS_SPLICE_PREDICTION"),),
                )
            return _Assessment(
                "UNDETERMINED", tuple(used), ("splice_rna_assay_reliability",),
                "An RNA result potentially contradicts the predicted splice event, but assay reliability is not established.",
            )

    if candidate_has_rna and reference_has_rna:
        if normalized(candidate_product.value or "") != normalized(reference_product.value or ""):
            if reliability.is_known and reliability.value == RnaAssayReliability.ESTABLISHED:
                return _Assessment(
                    "NOT_MET", tuple(used), (), "Reliable RNA data show different splice products.",
                    (exclusion(event_label, "RNA_PRODUCTS_DIFFER"),),
                )
            return _Assessment(
                "UNDETERMINED", tuple(used), ("splice_rna_assay_reliability",),
                "RNA assays report different splice products, but their reliability is not established.",
            )

    return _Assessment("OK", tuple(used))


def _matrix_strength(
    evidence: BiologicalFeatures,
    reference: ReferenceVariant,
    reference_index: int,
    motif: _Assessment,
) -> RuleEvaluation | _Assessment:
    """Apply the six published rows of Walker et al. 2023 Table 2 exactly."""

    baseline = evidence.applicable_splicing_pvs1_strength
    classification = reference.classification
    used = motif.used
    if baseline.is_conflicting:
        return _Assessment("UNDETERMINED", used, ("applicable_splicing_pvs1_strength",), "Applicable baseline PP3/PVS1 code is conflicting.")
    if baseline.is_unknown:
        return _Assessment(
            "UNDETERMINED", used, ("applicable_splicing_pvs1_strength",),
            "Walker et al. Table 2 strength cannot be selected without the independently curated applicable PP3/PVS1 baseline code.",
        )
    used = used + ("applicable_splicing_pvs1_strength", label(reference_index, "classification"))
    candidate_position = evidence.candidate.splice.motif_position.value
    reference_position = reference.annotations.splice.motif_position.value
    reference_classification = classification.value
    same_nucleotide = motif.same_nucleotide

    if same_nucleotide is True and candidate_position != reference_position:
        return _Assessment(
            "UNDETERMINED", used,
            ("candidate_splice_motif_position", label(reference_index, "splice_motif_position")),
            "The same supplied nucleotide position has incompatible canonical versus non-canonical motif annotations.",
        )

    if candidate_position == SpliceMotifPosition.WITHIN_MOTIF:
        if baseline.value != SpliceBaselineCode.PP3:
            return _Assessment("NOT_APPLICABLE", used, (), "The Table 2 non-canonical splice row requires PP3 as the candidate baseline code.")
        if same_nucleotide is True:
            strength = PS1Strength.STRONG if reference_classification == ReferenceClassification.PATHOGENIC else PS1Strength.MODERATE
            return activated(
                strength, PATH, CLINGEN_SPLICE_PS1,
                "Walker et al. 2023 Table 2: non-canonical candidate with PP3 at the same nucleotide as a reference with the same predicted splice event.",
                reference_index=reference_index, used=used,
            )
        if same_nucleotide is None:
            return _Assessment("UNDETERMINED", used, ("candidate_splice_nucleotide_position", label(reference_index, "splice_nucleotide_position")), "Nucleotide-position relation is strength-critical for the non-canonical Table 2 rows.")
        # Different nucleotides in the same motif: Table 2 row 2.
        strength = PS1Strength.MODERATE if reference_classification == ReferenceClassification.PATHOGENIC else PS1Strength.SUPPORTING
        return activated(
            strength, PATH, CLINGEN_SPLICE_PS1,
            "Walker et al. 2023 Table 2: non-canonical candidate with PP3 and a reference in the same splice motif with the same predicted event.",
            reference_index=reference_index, used=used,
        )

    if candidate_position != SpliceMotifPosition.CANONICAL_DINUCLEOTIDE:
        return _Assessment("NOT_APPLICABLE", used, (), "Candidate splice position is outside the generic Table 2 categories.")

    if baseline.value not in {
        SpliceBaselineCode.PVS1,
        SpliceBaselineCode.PVS1_STRONG,
        SpliceBaselineCode.PVS1_MODERATE,
        SpliceBaselineCode.PVS1_SUPPORTING,
    }:
        return _Assessment("NOT_APPLICABLE", used, (), "Canonical splice-site Table 2 rows require an applicable PVS1 baseline code.")

    if reference_position == SpliceMotifPosition.WITHIN_MOTIF:
        if baseline.value == SpliceBaselineCode.PVS1:
            strength = PS1Strength.SUPPORTING
        else:
            strength = PS1Strength.MODERATE if reference_classification == ReferenceClassification.PATHOGENIC else PS1Strength.SUPPORTING
        return activated(
            strength, PATH, CLINGEN_SPLICE_PS1,
            "Walker et al. 2023 Table 2: canonical candidate and comparison variant in the same splice motif outside the canonical dinucleotide with a matched predicted event.",
            reference_index=reference_index, used=used,
        )

    if reference_position != SpliceMotifPosition.CANONICAL_DINUCLEOTIDE:
        return _Assessment("NOT_APPLICABLE", used, (), "Comparison position is outside the generic canonical Table 2 categories.")

    if reference_classification == ReferenceClassification.LIKELY_PATHOGENIC:
        return not_applicable(
            PATH, CLINGEN_SPLICE_PS1,
            "Walker et al. 2023 Table 2 marks a comparison within the same canonical donor/acceptor dinucleotide to a Likely Pathogenic reference as not applicable.",
            reference_index=reference_index, used=used,
        )
    strength = PS1Strength.SUPPORTING if baseline.value == SpliceBaselineCode.PVS1 else PS1Strength.STRONG
    return activated(
        strength, PATH, CLINGEN_SPLICE_PS1,
        "Walker et al. 2023 Table 2: canonical candidate and Pathogenic reference within the same donor/acceptor +/-1,2 dinucleotide with the same predicted event.",
        reference_index=reference_index, used=used,
    )


def evaluate_splicing_path(
    evidence: BiologicalFeatures, reference: ReferenceVariant, reference_index: int
) -> RuleEvaluation:
    """Evaluate the generic ClinGen SVI splice-PS1 path for one reference."""

    if not _candidate_has_splice_context(evidence):
        return not_applicable(
            PATH, CLINGEN_SPLICE_PS1,
            "No candidate splice context was supplied for the generic splicing pathway.",
            reference_index=reference_index,
        )

    preconditions = common_reference_preconditions(
        evidence, reference, reference_index, path=PATH, guideline=CLINGEN_SPLICE_PS1
    )
    if preconditions.terminal is not None:
        return preconditions.terminal
    quality = require_reference_evidence_quality(
        evidence, reference, reference_index, path=PATH, guideline=CLINGEN_SPLICE_PS1,
        used=preconditions.used, exclusions=preconditions.exclusions,
    )
    if quality is not None:
        return quality

    # `require_reference_evidence_quality` has confirmed this assertion.  It
    # remains part of the minimum sufficient evidence and therefore needs to
    # appear in an activated pathway's audit/provenance output.
    precondition_used = preconditions.used + (
        label(reference_index, "classification_evidence"),
    )

    mechanism = _reference_splice_mechanism(reference, reference_index)
    if mechanism.state == "NOT_MET":
        return not_met(PATH, CLINGEN_SPLICE_PS1, mechanism.reason, reference_index=reference_index,
            used=precondition_used + mechanism.used,
            exclusions=preconditions.exclusions + mechanism.exclusions)
    if mechanism.state == "UNDETERMINED":
        return undetermined(PATH, CLINGEN_SPLICE_PS1, mechanism.reason, mechanism.missing,
            reference_index=reference_index, used=precondition_used + mechanism.used,
            exclusions=preconditions.exclusions + mechanism.exclusions)
    precondition_used += mechanism.used

    transcript = _transcript_compatibility(evidence, reference, reference_index)
    if transcript.state == "NOT_MET":
        return not_met(PATH, CLINGEN_SPLICE_PS1, transcript.reason, reference_index=reference_index,
            used=precondition_used + transcript.used,
            exclusions=preconditions.exclusions + transcript.exclusions)
    if transcript.state == "UNDETERMINED":
        return undetermined(PATH, CLINGEN_SPLICE_PS1, transcript.reason, transcript.missing, reference_index=reference_index,
            used=precondition_used + transcript.used,
            exclusions=preconditions.exclusions + transcript.exclusions)
    core_used = precondition_used + transcript.used

    events = _events_and_predictions_equivalent(evidence, reference, reference_index)
    if events.state == "NOT_MET":
        return not_met(PATH, CLINGEN_SPLICE_PS1, events.reason, reference_index=reference_index,
            used=core_used + events.used, exclusions=preconditions.exclusions + events.exclusions)
    if events.state == "UNDETERMINED":
        return undetermined(PATH, CLINGEN_SPLICE_PS1, events.reason, events.missing, reference_index=reference_index,
            used=core_used + events.used, exclusions=preconditions.exclusions + events.exclusions)

    motif = _motif_relationship(evidence, reference, reference_index)
    if motif.state == "NOT_MET":
        return not_met(PATH, CLINGEN_SPLICE_PS1, motif.reason, reference_index=reference_index,
            used=core_used + events.used + motif.used, exclusions=preconditions.exclusions + motif.exclusions)
    if motif.state == "UNDETERMINED":
        return undetermined(PATH, CLINGEN_SPLICE_PS1, motif.reason, motif.missing, reference_index=reference_index,
            used=core_used + events.used + motif.used, exclusions=preconditions.exclusions + motif.exclusions)
    if motif.state == "NOT_APPLICABLE":
        return not_applicable(PATH, CLINGEN_SPLICE_PS1, motif.reason, reference_index=reference_index,
            used=core_used + events.used + motif.used)

    protein = _exonic_protein_effect_consistency(evidence, reference, reference_index)
    if protein.state == "NOT_MET":
        return not_met(PATH, CLINGEN_SPLICE_PS1, protein.reason, reference_index=reference_index,
            used=core_used + events.used + motif.used + protein.used,
            exclusions=preconditions.exclusions + events.exclusions + motif.exclusions + protein.exclusions)
    if protein.state == "UNDETERMINED":
        return undetermined(PATH, CLINGEN_SPLICE_PS1, protein.reason, protein.missing, reference_index=reference_index,
            used=core_used + events.used + motif.used + protein.used,
            exclusions=preconditions.exclusions + events.exclusions + motif.exclusions + protein.exclusions)

    rna = _rna_consistency(evidence, reference, reference_index)
    if rna.state == "NOT_MET":
        return not_met(PATH, CLINGEN_SPLICE_PS1, rna.reason, reference_index=reference_index,
            used=core_used + events.used + motif.used + protein.used + rna.used,
            exclusions=preconditions.exclusions + events.exclusions + motif.exclusions + protein.exclusions + rna.exclusions)
    if rna.state == "UNDETERMINED":
        return undetermined(PATH, CLINGEN_SPLICE_PS1, rna.reason, rna.missing, reference_index=reference_index,
            used=core_used + events.used + motif.used + protein.used + rna.used,
            exclusions=preconditions.exclusions + events.exclusions + motif.exclusions + protein.exclusions + rna.exclusions)

    functional = _functional_consequence_consistency(evidence)
    if functional.state == "UNDETERMINED":
        return undetermined(PATH, CLINGEN_SPLICE_PS1, functional.reason, functional.missing, reference_index=reference_index,
            used=core_used + events.used + motif.used + protein.used + rna.used + functional.used,
            exclusions=preconditions.exclusions + events.exclusions + motif.exclusions + protein.exclusions + rna.exclusions)

    matrix = _matrix_strength(evidence, reference, reference_index, motif)
    if isinstance(matrix, RuleEvaluation):
        # Matrix output contains its own source and core fields; preserve all
        # prerequisite audit evidence accumulated by this pathway.
        return RuleEvaluation(
            status=matrix.status,
            strength=matrix.strength,
            decision_path=matrix.decision_path,
            guideline=matrix.guideline,
            reason=matrix.reason,
            reference_index=reference_index,
            features_used=tuple(dict.fromkeys(core_used + events.used + motif.used + protein.used + rna.used + functional.used + matrix.features_used)),
            missing_critical_features=matrix.missing_critical_features,
            exclusions_evaluated=preconditions.exclusions + events.exclusions + motif.exclusions + protein.exclusions + rna.exclusions + matrix.exclusions_evaluated,
        )
    if matrix.state == "UNDETERMINED":
        return undetermined(PATH, CLINGEN_SPLICE_PS1, matrix.reason, matrix.missing, reference_index=reference_index,
            used=core_used + events.used + motif.used + protein.used + rna.used + functional.used + matrix.used,
            exclusions=preconditions.exclusions + events.exclusions + motif.exclusions + protein.exclusions + rna.exclusions + matrix.exclusions)
    if matrix.state == "NOT_APPLICABLE":
        return not_applicable(PATH, CLINGEN_SPLICE_PS1, matrix.reason, reference_index=reference_index,
            used=core_used + events.used + motif.used + protein.used + rna.used + functional.used + matrix.used)
    return not_met(PATH, CLINGEN_SPLICE_PS1, matrix.reason, reference_index=reference_index,
        used=core_used + events.used + motif.used + protein.used + rna.used + functional.used + matrix.used,
        exclusions=preconditions.exclusions + events.exclusions + motif.exclusions + protein.exclusions + rna.exclusions + matrix.exclusions)

