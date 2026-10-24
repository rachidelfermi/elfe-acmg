"""Deterministic, non-scoring PS1 decision engine."""

from __future__ import annotations

from dataclasses import replace
from typing import Any, Iterable, Mapping, Optional, Union

from .evidence import Evidence, Provenance, stable_unique
from .guidelines import ACMG_PS1, CLINGEN_SPLICE_PS1
from .models import (
    BiologicalFeatures,
    DecisionStatus,
    PS1Strength,
    ReferenceVariant,
)
from .output import GuidelineReference, PS1Decision, RuleEvaluation
from .rules.common import common_reference_preconditions, label, normalized, not_applicable, undetermined
from .rules.missense import evaluate_standard_amino_acid_path
from .rules.splicing import evaluate_splicing_path
from .rules.vcep import evaluate_selected_vcep_rule, select_vcep_rule


_STRENGTH_RANK = {
    PS1Strength.NONE: -1,
    PS1Strength.SUPPORTING: 0,
    PS1Strength.MODERATE: 1,
    PS1Strength.STRONG: 2,
}


def _merge(base: RuleEvaluation, addition: RuleEvaluation) -> RuleEvaluation:
    """Attach universal-precondition audit records to a chosen VCEP result."""

    return replace(
        addition,
        features_used=stable_unique(base.features_used + addition.features_used),
        missing_critical_features=stable_unique(base.missing_critical_features + addition.missing_critical_features),
        exclusions_evaluated=base.exclusions_evaluated + addition.exclusions_evaluated,
    )


def _resolve_paths(paths: tuple[RuleEvaluation, ...], reference_index: int) -> RuleEvaluation:
    """Resolve standard and splice paths for one reference without scoring."""

    activated_paths = [item for item in paths if item.status == DecisionStatus.ACTIVATED]
    if activated_paths:
        # A strength is selected only because one complete published path
        # independently permits it.  Multiple paths are never added together.
        return max(
            activated_paths,
            key=lambda item: (_STRENGTH_RANK[item.strength], -paths.index(item)),
        )
    unresolved = [item for item in paths if item.status == DecisionStatus.UNDETERMINED]
    if unresolved:
        # Combine only missing requirements from independently applicable,
        # unresolved paths; no missing optional fields are introduced.
        first = unresolved[0]
        return replace(
            first,
            missing_critical_features=stable_unique(
                tuple(name for item in unresolved for name in item.missing_critical_features)
            ),
            exclusions_evaluated=tuple(
                exclusion for item in unresolved for exclusion in item.exclusions_evaluated
            ),
        )
    not_met_paths = [item for item in paths if item.status == DecisionStatus.NOT_MET]
    if not_met_paths:
        return not_met_paths[0]
    return paths[0] if paths else not_applicable(
        "PS1_NO_PATH", ACMG_PS1, "No supported PS1 pathway was identified.", reference_index=reference_index
    )


def _reference_identifier(reference: ReferenceVariant, index: int) -> str:
    value = reference.annotations.variant_id
    if value.is_known and value.value:
        return value.value
    return f"reference_{index + 1}"


def _reference_summary(reference: ReferenceVariant, index: int, evaluation: RuleEvaluation) -> dict[str, str]:
    classification = reference.classification.value.value if reference.classification.is_known and reference.classification.value else "UNKNOWN"
    return {
        "identifier": _reference_identifier(reference, index),
        "classification": classification,
        "status": evaluation.status.value,
        "strength": evaluation.strength.value,
        "decision_path": evaluation.decision_path,
    }


def _candidate_field(evidence: BiologicalFeatures, suffix: str) -> Optional[Evidence[Any]]:
    direct: dict[str, Evidence[Any]] = {
        "gene_id": evidence.gene_id,
        "disease_id": evidence.disease_id,
        "inheritance_context": evidence.inheritance_context,
        "protein_substitution_equivalence": evidence.protein_substitution_equivalence,
        "underlying_allele_equivalence": evidence.underlying_allele_equivalence,
        "protein_effect_mechanism_compatibility": evidence.protein_effect_mechanism_compatibility,
        "splice_event_equivalence": evidence.splice_event_equivalence,
        "splice_motif_equivalence": evidence.splice_motif_equivalence,
        "splice_prediction_strength_comparable": evidence.splice_prediction_strength_comparable,
        "splice_transcript_equivalence": evidence.splice_transcript_equivalence,
        "applicable_splicing_pvs1_strength": evidence.applicable_splicing_pvs1_strength,
        "splice_product_functional_consequence": evidence.splice_product_functional_consequence,
        "splice_rna_assay_reliability": evidence.splice_rna_assay_reliability,
        "independent_reference_classifications": evidence.independent_reference_classifications,
        "vcep_specifications": evidence.vcep_specifications,
        "vcep_gene_scope": evidence.vcep_gene_scope,
        "vcep_disease_scope": evidence.vcep_disease_scope,
        "vcep_specification_version": evidence.vcep_specification_version,
        "vcep_ps1_eligibility_rules": evidence.vcep_ps1_eligibility_rules,
        "vcep_ps1_strength_rules": evidence.vcep_ps1_strength_rules,
        "vcep_splicing_rules": evidence.vcep_splicing_rules,
        "vcep_reference_classification_requirements": evidence.vcep_reference_classification_requirements,
        "vcep_exclusion_rules": evidence.vcep_exclusion_rules,
        "feature_conflict_status": evidence.feature_conflict_status,
    }
    if suffix in direct:
        return direct[suffix]
    candidate = evidence.candidate
    direct_candidate: dict[str, Evidence[Any]] = {
        "variant_id": candidate.variant_id,
        "genomic_allele": candidate.genomic_allele,
        "coding_hgvs": candidate.coding_hgvs,
        "protein_hgvs": candidate.protein_hgvs,
        "transcript_id_version": candidate.transcript_id_version,
        "protein_accession_version": candidate.protein_accession_version,
        "amino_acid_position": candidate.amino_acid_position,
        "reference_amino_acid": candidate.reference_amino_acid,
        "alternate_amino_acid": candidate.alternate_amino_acid,
        "variant_consequence": candidate.consequence,
        "splice_impact": candidate.splice.impact,
        "splice_site_type": candidate.splice.site_type,
        "splice_motif_position": candidate.splice.motif_position,
        "splice_motif_id": candidate.splice.motif_id,
        "splice_nucleotide_position": candidate.splice.nucleotide_position,
        "splice_prediction": candidate.splice.prediction,
        "spliceai_event_type": candidate.splice.spliceai_event_type,
        "spliceai_event_position": candidate.splice.spliceai_event_position,
        "splice_prediction_magnitude": candidate.splice.prediction_magnitude,
        "rna_splicing_evidence": candidate.splice.rna_splicing_evidence,
        "observed_splice_product": candidate.splice.observed_splice_product,
        "residual_normal_splicing": candidate.splice.residual_normal_splicing,
    }
    return direct_candidate.get(suffix)


def _reference_field(reference: ReferenceVariant, suffix: str) -> Optional[Evidence[Any]]:
    ann = reference.annotations
    direct: dict[str, Evidence[Any]] = {
        "variant_id": ann.variant_id,
        "gene_id": reference.gene_id,
        "disease_id": reference.disease_id,
        "genomic_allele": ann.genomic_allele,
        "coding_hgvs": ann.coding_hgvs,
        "protein_hgvs": ann.protein_hgvs,
        "transcript_id_version": ann.transcript_id_version,
        "protein_accession_version": ann.protein_accession_version,
        "amino_acid_position": ann.amino_acid_position,
        "reference_amino_acid": ann.reference_amino_acid,
        "alternate_amino_acid": ann.alternate_amino_acid,
        "variant_consequence": ann.consequence,
        "classification": reference.classification,
        "classification_source": reference.classification_source,
        "classification_date": reference.classification_date,
        "classification_review_status": reference.review_status,
        "classification_evidence": reference.classification_evidence,
        "classification_independence": reference.classification_independence,
        "disease_association": reference.disease_association,
        "molecular_mechanism": reference.molecular_mechanism,
        "classification_conflicts": reference.classification_conflicts,
        "protein_substitution_equivalence": reference.protein_substitution_equivalence,
        "underlying_allele_equivalence": reference.underlying_allele_equivalence,
        "protein_effect_mechanism_compatibility": reference.protein_effect_mechanism_compatibility,
        "splice_event_equivalence": reference.splice_event_equivalence,
        "splice_motif_equivalence": reference.splice_motif_equivalence,
        "splice_prediction_strength_comparable": reference.splice_prediction_strength_comparable,
        "splice_transcript_equivalence": reference.splice_transcript_equivalence,
        "splice_impact": ann.splice.impact,
        "splice_site_type": ann.splice.site_type,
        "splice_motif_position": ann.splice.motif_position,
        "splice_motif_id": ann.splice.motif_id,
        "splice_nucleotide_position": ann.splice.nucleotide_position,
        "splice_prediction": ann.splice.prediction,
        "spliceai_event_type": ann.splice.spliceai_event_type,
        "spliceai_event_position": ann.splice.spliceai_event_position,
        "splice_prediction_magnitude": ann.splice.prediction_magnitude,
        "rna_splicing_evidence": ann.splice.rna_splicing_evidence,
        "observed_splice_product": ann.splice.observed_splice_product,
        "residual_normal_splicing": ann.splice.residual_normal_splicing,
    }
    return direct.get(suffix)


def _provenance_for_features(
    evidence: BiologicalFeatures, fields: Iterable[str]
) -> tuple[Provenance, ...]:
    provenance: list[Provenance] = []
    for field in stable_unique(tuple(fields)):
        value: Optional[Evidence[Any]] = None
        if field.startswith("candidate_"):
            value = _candidate_field(evidence, field[len("candidate_"):])
        elif field.startswith("references["):
            close = field.find("]")
            if close != -1 and len(field) > close + 2:
                try:
                    index = int(field[len("references["):close])
                except ValueError:
                    index = -1
                if 0 <= index < len(evidence.references):
                    value = _reference_field(evidence.references[index], field[close + 2:])
        else:
            value = _candidate_field(evidence, field)
        if value is not None:
            provenance.extend(value.provenance)
    return tuple(dict.fromkeys(provenance))


class PS1Engine:
    """Evaluate supplied biological features under PS1 only.

    The engine does not combine criteria, download classifications, or run
    splice predictors.  It is deterministic and does not mutate its input.
    """

    def evaluate(self, biological_features: Union[BiologicalFeatures, Mapping[str, Any]]) -> PS1Decision:
        evidence = (
            biological_features
            if isinstance(biological_features, BiologicalFeatures)
            else BiologicalFeatures.from_mapping(biological_features)
        )
        if not evidence.references:
            candidate = undetermined(
                "PS1_REFERENCE_REQUIRED", ACMG_PS1,
                "PS1 requires at least one independently assessed comparison reference variant.",
                ("references",),
            )
            return self._finalize(candidate, evidence, (), ())

        for index, reference in enumerate(evidence.references):
            for other_index, other in enumerate(evidence.references[:index]):
                first_id = reference.annotations.variant_id
                other_id = other.annotations.variant_id
                if not (first_id.is_known and other_id.is_known and
                        normalized(first_id.value or "") == normalized(other_id.value or "")):
                    continue
                if (reference.disease_id.is_known and other.disease_id.is_known and
                        normalized(reference.disease_id.value or "") != normalized(other.disease_id.value or "")):
                    continue
                first_class = reference.classification
                other_class = other.classification
                if (first_class.is_known and other_class.is_known and
                        first_class.value != other_class.value):
                    fields = (label(other_index, "classification"), label(index, "classification"))
                    conflict = undetermined(
                        "PS1_DUPLICATE_REFERENCE_CONFLICT", ACMG_PS1,
                        "The same reference variant has conflicting classifications for the disease context.",
                        fields, used=(label(other_index, "variant_id"), label(index, "variant_id")),
                    )
                    return self._finalize(conflict, evidence, (), ())

        reference_results: list[RuleEvaluation] = []
        reference_paths: list[tuple[RuleEvaluation, ...]] = []
        for index, reference in enumerate(evidence.references):
            # Universal safety constraints (disease context, independent P/LP
            # classification, and a distinct allele) also apply to VCEP rules.
            universal = common_reference_preconditions(
                evidence, reference, index, path="PS1_REFERENCE_ELIGIBILITY", guideline=ACMG_PS1
            )
            if universal.terminal is not None:
                reference_results.append(universal.terminal)
                reference_paths.append((universal.terminal,))
                continue

            vcep = evaluate_selected_vcep_rule(
                select_vcep_rule(evidence, reference, index), evidence, reference, index
            )
            if vcep is not None:
                result = _merge(
                    RuleEvaluation(
                        status=DecisionStatus.NOT_APPLICABLE,
                        strength=PS1Strength.NONE,
                        decision_path="PS1_REFERENCE_ELIGIBILITY",
                        guideline=ACMG_PS1,
                        reason="Universal PS1 comparison preconditions are satisfied.",
                        reference_index=index,
                        features_used=universal.used,
                        exclusions_evaluated=universal.exclusions,
                    ),
                    vcep,
                )
                reference_results.append(result)
                reference_paths.append((result,))
                continue

            standard = evaluate_standard_amino_acid_path(evidence, reference, index)
            splicing = evaluate_splicing_path(evidence, reference, index)
            resolved = _resolve_paths((standard, splicing), index)
            reference_results.append(resolved)
            reference_paths.append((standard, splicing))

        final = self._resolve_references(reference_results)
        return self._finalize(final, evidence, tuple(reference_results), tuple(reference_paths))

    @staticmethod
    def _resolve_references(results: list[RuleEvaluation]) -> RuleEvaluation:
        activated = [result for result in results if result.status == DecisionStatus.ACTIVATED]
        if activated:
            return max(
                activated,
                key=lambda result: (_STRENGTH_RANK[result.strength], -results.index(result)),
            )
        unresolved = [result for result in results if result.status == DecisionStatus.UNDETERMINED]
        if unresolved:
            first = unresolved[0]
            return replace(
                first,
                missing_critical_features=stable_unique(
                    tuple(field for result in unresolved for field in result.missing_critical_features)
                ),
            )
        not_met_results = [result for result in results if result.status == DecisionStatus.NOT_MET]
        if not_met_results:
            return not_met_results[0]
        return results[0] if results else not_applicable(
            "PS1_NO_REFERENCE", ACMG_PS1, "No comparison reference was available."
        )

    @staticmethod
    def _finalize(
        final: RuleEvaluation,
        evidence: BiologicalFeatures,
        reference_results: tuple[RuleEvaluation, ...],
        reference_paths: tuple[tuple[RuleEvaluation, ...], ...],
    ) -> PS1Decision:
        activated_results = [item for item in reference_results if item.status == DecisionStatus.ACTIVATED]
        reference_variants_used: tuple[Mapping[str, str], ...] = ()
        if activated_results:
            reference_variants_used = tuple(
                _reference_summary(evidence.references[item.reference_index or 0], item.reference_index or 0, item)
                for item in activated_results
            )
        # When several references independently activate PS1, every one is a
        # reported supporting comparison.  Preserve the features/provenance for
        # all of them rather than only the strongest tie-breaking result.
        audit_results = activated_results if final.status == DecisionStatus.ACTIVATED else [final]
        final_features_used = stable_unique(
            tuple(field for result in audit_results for field in result.features_used)
        )
        final_missing = stable_unique(
            tuple(field for result in audit_results for field in result.missing_critical_features)
        )
        final_exclusions = tuple(
            exclusion for result in audit_results for exclusion in result.exclusions_evaluated
        )
        feature_fields = final_features_used + final_missing
        for exclusion in final_exclusions:
            feature = exclusion.get("feature")
            if feature:
                feature_fields += (feature,)
        provenance = _provenance_for_features(evidence, feature_fields)

        audit_rows: list[Mapping[str, Any]] = []
        for index, (result, paths) in enumerate(zip(reference_results, reference_paths)):
            row: dict[str, Any] = _reference_summary(evidence.references[index], index, result)
            row["reason"] = result.reason
            row["pathways"] = [
                {
                    "decision_path": path.decision_path,
                    "status": path.status.value,
                    "strength": path.strength.value,
                    "guideline": path.guideline.to_dict(),
                    "reason": path.reason,
                    "features_used": list(path.features_used),
                    "missing_critical_features": list(path.missing_critical_features),
                    "exclusions_evaluated": [dict(exclusion) for exclusion in path.exclusions_evaluated],
                }
                for path in paths
            ]
            audit_rows.append(row)
        reason = final.reason
        if len(activated_results) > 1:
            reason += " Multiple references were evaluated independently; the final strength is the strongest single complete pathway and was not summed across references."
        return PS1Decision(
            status=final.status,
            strength=final.strength,
            decision_path=final.decision_path,
            guideline=final.guideline,
            reason=reason,
            reference_variants_used=reference_variants_used,
            features_used=final_features_used,
            missing_critical_features=final_missing,
            exclusions_evaluated=final_exclusions,
            evidence_provenance=provenance,
            reference_evaluations=tuple(audit_rows),
        )


def evaluate_ps1(
    biological_features: Union[BiologicalFeatures, Mapping[str, Any]],
) -> PS1Decision:
    """Evaluate a pre-extracted biological feature set under PS1 only."""

    return PS1Engine().evaluate(biological_features)

