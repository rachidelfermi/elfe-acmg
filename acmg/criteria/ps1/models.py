"""Typed biological feature schema for standalone ACMG/ClinGen PS1 evaluation.

The schema intentionally models candidate and reference observations
separately.  It does not annotate variants, resolve HGVS, look up ClinVar, or
infer a gene--disease relationship; those jobs belong to the upstream pipeline
and its curation process.  All fields are optional evidence assertions.  A
rule asks only for the assertions that can change its own decision.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Mapping, Optional

from .evidence import (
    Evidence,
    EvidenceState,
    SchemaValidationError,
    _MISSING,
    parse_bool,
    parse_enum,
    parse_evidence,
    parse_number,
    parse_string,
)


class TextEnum(str, Enum):
    """Enum which serializes directly to a JSON string."""


class DecisionStatus(TextEnum):
    ACTIVATED = "ACTIVATED"
    NOT_MET = "NOT_MET"
    UNDETERMINED = "UNDETERMINED"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class PS1Strength(TextEnum):
    STRONG = "STRONG"
    MODERATE = "MODERATE"
    SUPPORTING = "SUPPORTING"
    NONE = "NONE"


class ReferenceClassification(TextEnum):
    PATHOGENIC = "PATHOGENIC"
    LIKELY_PATHOGENIC = "LIKELY_PATHOGENIC"
    VUS = "VUS"
    LIKELY_BENIGN = "LIKELY_BENIGN"
    BENIGN = "BENIGN"


class ClassificationEvidence(TextEnum):
    """Curated adequacy of the reference classification's evidence record."""

    ESTABLISHED = "ESTABLISHED"
    INSUFFICIENT = "INSUFFICIENT"


class ReviewStatus(TextEnum):
    EXPERT_PANEL = "EXPERT_PANEL"
    PRACTICE_GUIDELINE = "PRACTICE_GUIDELINE"
    MULTIPLE_SUBMITTERS_NO_CONFLICT = "MULTIPLE_SUBMITTERS_NO_CONFLICT"
    SINGLE_SUBMITTER = "SINGLE_SUBMITTER"
    OTHER = "OTHER"


class MolecularMechanism(TextEnum):
    PROTEIN_ALTERATION = "PROTEIN_ALTERATION"
    SPLICING = "SPLICING"
    OTHER = "OTHER"


class VariantConsequence(TextEnum):
    MISSENSE = "MISSENSE"
    SPLICE_REGION = "SPLICE_REGION"
    SPLICE_DONOR = "SPLICE_DONOR"
    SPLICE_ACCEPTOR = "SPLICE_ACCEPTOR"
    SYNONYMOUS = "SYNONYMOUS"
    INTRONIC = "INTRONIC"
    OTHER = "OTHER"


class SpliceImpact(TextEnum):
    NO_IMPACT = "NO_IMPACT"
    POSSIBLE = "POSSIBLE"
    ABERRANT = "ABERRANT"


class SplicePrediction(TextEnum):
    NO_IMPACT = "NO_IMPACT"
    SPLICE_ALTERING = "SPLICE_ALTERING"
    UNCERTAIN = "UNCERTAIN"


class SpliceSiteType(TextEnum):
    DONOR = "DONOR"
    ACCEPTOR = "ACCEPTOR"
    NONE = "NONE"
    OTHER = "OTHER"


class SpliceMotifPosition(TextEnum):
    CANONICAL_DINUCLEOTIDE = "CANONICAL_DINUCLEOTIDE"
    WITHIN_MOTIF = "WITHIN_MOTIF"
    OUTSIDE_MOTIF = "OUTSIDE_MOTIF"


class SpliceProductFunctionalConsequence(TextEnum):
    LOSS_OF_FUNCTION = "LOSS_OF_FUNCTION"
    NOT_LOSS_OF_FUNCTION = "NOT_LOSS_OF_FUNCTION"
    UNKNOWN = "UNKNOWN"


class RnaAssayReliability(TextEnum):
    ESTABLISHED = "ESTABLISHED"
    LIMITED = "LIMITED"
    INADEQUATE = "INADEQUATE"


class SpliceBaselineCode(TextEnum):
    """Already-curated baseline code required by Walker et al. Table 2.

    This is an input, not a call to a PVS1 evaluator.  A PS1-only module must
    not recreate gene/disease LoF logic to determine it.
    """

    PP3 = "PP3"
    PVS1 = "PVS1"
    PVS1_STRONG = "PVS1_STRONG"
    PVS1_MODERATE = "PVS1_MODERATE"
    PVS1_SUPPORTING = "PVS1_SUPPORTING"


def _enum_parser(
    enum_type: type[Enum], aliases: Optional[Mapping[str, Enum]] = None
) -> Callable[[Any, str], Enum]:
    return lambda raw, field_name: parse_enum(enum_type, raw, field_name, aliases)


def _positive_int(raw: Any, field_name: str) -> int:
    if isinstance(raw, bool) or not isinstance(raw, int) or raw < 1:
        raise SchemaValidationError(f"{field_name} must be a positive integer")
    return raw


def _score(raw: Any, field_name: str) -> float:
    value = parse_number(raw, field_name)
    if not 0 <= value <= 1:
        raise SchemaValidationError(f"{field_name} must be between 0 and 1")
    return value


_CLASSIFICATION_ALIASES: Mapping[str, ReferenceClassification] = {
    "P": ReferenceClassification.PATHOGENIC,
    "LP": ReferenceClassification.LIKELY_PATHOGENIC,
    "UNCERTAIN_SIGNIFICANCE": ReferenceClassification.VUS,
    "LIKELY_B": ReferenceClassification.LIKELY_BENIGN,
    "B": ReferenceClassification.BENIGN,
}

_CONSEQUENCE_ALIASES: Mapping[str, VariantConsequence] = {
    "MISSENSE_VARIANT": VariantConsequence.MISSENSE,
    "SPLICE_VARIANT": VariantConsequence.SPLICE_REGION,
    "SPLICE_SITE": VariantConsequence.SPLICE_REGION,
    "SPLICE_DONOR_VARIANT": VariantConsequence.SPLICE_DONOR,
    "SPLICE_ACCEPTOR_VARIANT": VariantConsequence.SPLICE_ACCEPTOR,
    "SYNONYMOUS_VARIANT": VariantConsequence.SYNONYMOUS,
}

_SPLICE_IMPACT_ALIASES: Mapping[str, SpliceImpact] = {
    "NORMAL": SpliceImpact.NO_IMPACT,
    "NO_SPLICE_EFFECT": SpliceImpact.NO_IMPACT,
    "NO_PREDICTED_OR_OBSERVED_IMPACT": SpliceImpact.NO_IMPACT,
    "SPLICE_ALTERING": SpliceImpact.ABERRANT,
    "ABNORMAL": SpliceImpact.ABERRANT,
}

_SPLICE_PREDICTION_ALIASES: Mapping[str, SplicePrediction] = {
    "NORMAL": SplicePrediction.NO_IMPACT,
    "NO_SPLICE_EFFECT": SplicePrediction.NO_IMPACT,
    "ABERRANT": SplicePrediction.SPLICE_ALTERING,
    "SPLICE_IMPACT": SplicePrediction.SPLICE_ALTERING,
}

_MOTIF_ALIASES: Mapping[str, SpliceMotifPosition] = {
    "CANONICAL": SpliceMotifPosition.CANONICAL_DINUCLEOTIDE,
    "PLUS_MINUS_1_2": SpliceMotifPosition.CANONICAL_DINUCLEOTIDE,
    "PM_1_2": SpliceMotifPosition.CANONICAL_DINUCLEOTIDE,
    "SPLICE_MOTIF": SpliceMotifPosition.WITHIN_MOTIF,
    "OUTSIDE_CANONICAL": SpliceMotifPosition.WITHIN_MOTIF,
}

_BASELINE_ALIASES: Mapping[str, SpliceBaselineCode] = {
    "PVS1_VERY_STRONG": SpliceBaselineCode.PVS1,
    "PVS1_VS": SpliceBaselineCode.PVS1,
    "PVS1_S": SpliceBaselineCode.PVS1_STRONG,
    "PVS1_M": SpliceBaselineCode.PVS1_MODERATE,
    "PVS1_P": SpliceBaselineCode.PVS1_SUPPORTING,
}


@dataclass(frozen=True)
class SpliceEvidence:
    """Candidate- or reference-specific splicing observations."""

    impact: Evidence[SpliceImpact] = field(default_factory=Evidence.unknown)
    site_type: Evidence[SpliceSiteType] = field(default_factory=Evidence.unknown)
    motif_position: Evidence[SpliceMotifPosition] = field(default_factory=Evidence.unknown)
    motif_id: Evidence[str] = field(default_factory=Evidence.unknown)
    nucleotide_position: Evidence[str] = field(default_factory=Evidence.unknown)
    prediction: Evidence[SplicePrediction] = field(default_factory=Evidence.unknown)
    spliceai_event_type: Evidence[str] = field(default_factory=Evidence.unknown)
    spliceai_event_position: Evidence[str] = field(default_factory=Evidence.unknown)
    prediction_magnitude: Evidence[float] = field(default_factory=Evidence.unknown)
    rna_splicing_evidence: Evidence[bool] = field(default_factory=Evidence.unknown)
    observed_splice_product: Evidence[str] = field(default_factory=Evidence.unknown)
    residual_normal_splicing: Evidence[bool] = field(default_factory=Evidence.unknown)

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any], prefix: str) -> "SpliceEvidence":
        """Read either prefixed flat input or an unprefixed reference object."""

        def take(name: str, *aliases: str) -> Any:
            for candidate in (f"{prefix}{name}", name, *aliases):
                if candidate in raw:
                    return raw[candidate]
            return _MISSING

        return cls(
            impact=parse_evidence(
                take("splice_impact"),
                _enum_parser(SpliceImpact, _SPLICE_IMPACT_ALIASES),
                f"{prefix}splice_impact",
            ),
            site_type=parse_evidence(
                take("splice_site_type"),
                _enum_parser(SpliceSiteType),
                f"{prefix}splice_site_type",
            ),
            motif_position=parse_evidence(
                take("splice_motif_position"),
                _enum_parser(SpliceMotifPosition, _MOTIF_ALIASES),
                f"{prefix}splice_motif_position",
            ),
            motif_id=parse_evidence(
                take("splice_motif_id"), parse_string, f"{prefix}splice_motif_id"
            ),
            nucleotide_position=parse_evidence(
                take("splice_nucleotide_position"), parse_string, f"{prefix}splice_nucleotide_position"
            ),
            prediction=parse_evidence(
                take("splice_prediction"),
                _enum_parser(SplicePrediction, _SPLICE_PREDICTION_ALIASES),
                f"{prefix}splice_prediction",
            ),
            spliceai_event_type=parse_evidence(
                take("spliceai_event_type", f"{prefix}splice_event"),
                parse_string,
                f"{prefix}spliceai_event_type",
            ),
            spliceai_event_position=parse_evidence(
                take("spliceai_event_position", f"{prefix}splice_event_position"),
                parse_string,
                f"{prefix}spliceai_event_position",
            ),
            prediction_magnitude=parse_evidence(
                take("splice_prediction_magnitude"), _score, f"{prefix}splice_prediction_magnitude"
            ),
            rna_splicing_evidence=parse_evidence(
                take("rna_splicing_evidence"),
                parse_bool,
                f"{prefix}rna_splicing_evidence",
                boolean=True,
            ),
            observed_splice_product=parse_evidence(
                take("observed_splice_product"), parse_string, f"{prefix}observed_splice_product"
            ),
            residual_normal_splicing=parse_evidence(
                take("residual_normal_splicing"),
                parse_bool,
                f"{prefix}residual_normal_splicing",
                boolean=True,
            ),
        )


@dataclass(frozen=True)
class VariantAnnotations:
    """Normalized annotations for one candidate or comparison variant."""

    variant_id: Evidence[str] = field(default_factory=Evidence.unknown)
    genomic_allele: Evidence[str] = field(default_factory=Evidence.unknown)
    coding_hgvs: Evidence[str] = field(default_factory=Evidence.unknown)
    protein_hgvs: Evidence[str] = field(default_factory=Evidence.unknown)
    transcript_id_version: Evidence[str] = field(default_factory=Evidence.unknown)
    protein_accession_version: Evidence[str] = field(default_factory=Evidence.unknown)
    amino_acid_position: Evidence[int] = field(default_factory=Evidence.unknown)
    reference_amino_acid: Evidence[str] = field(default_factory=Evidence.unknown)
    alternate_amino_acid: Evidence[str] = field(default_factory=Evidence.unknown)
    consequence: Evidence[VariantConsequence] = field(default_factory=Evidence.unknown)
    splice: SpliceEvidence = field(default_factory=SpliceEvidence)

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any], prefix: str) -> "VariantAnnotations":
        def take(name: str, *aliases: str) -> Any:
            for candidate in (f"{prefix}{name}", name, *aliases):
                if candidate in raw:
                    return raw[candidate]
            return _MISSING

        return cls(
            variant_id=parse_evidence(take("variant_id", "id", "identifier"), parse_string, f"{prefix}variant_id"),
            genomic_allele=parse_evidence(take("genomic_allele"), parse_string, f"{prefix}genomic_allele"),
            coding_hgvs=parse_evidence(take("coding_hgvs", "nucleotide_change"), parse_string, f"{prefix}coding_hgvs"),
            protein_hgvs=parse_evidence(take("protein_hgvs", "protein_change"), parse_string, f"{prefix}protein_hgvs"),
            transcript_id_version=parse_evidence(take("transcript_id_version"), parse_string, f"{prefix}transcript_id_version"),
            protein_accession_version=parse_evidence(take("protein_accession_version"), parse_string, f"{prefix}protein_accession_version"),
            amino_acid_position=parse_evidence(take("amino_acid_position"), _positive_int, f"{prefix}amino_acid_position"),
            reference_amino_acid=parse_evidence(take("reference_amino_acid"), parse_string, f"{prefix}reference_amino_acid"),
            alternate_amino_acid=parse_evidence(take("alternate_amino_acid"), parse_string, f"{prefix}alternate_amino_acid"),
            consequence=parse_evidence(
                take("variant_consequence", "consequence"),
                _enum_parser(VariantConsequence, _CONSEQUENCE_ALIASES),
                f"{prefix}variant_consequence",
            ),
            splice=SpliceEvidence.from_mapping(raw, prefix),
        )


@dataclass(frozen=True)
class ReferenceVariant:
    """A comparison variant and its pathogenicity assertions."""

    annotations: VariantAnnotations
    gene_id: Evidence[str] = field(default_factory=Evidence.unknown)
    disease_id: Evidence[str] = field(default_factory=Evidence.unknown)
    classification: Evidence[ReferenceClassification] = field(default_factory=Evidence.unknown)
    classification_source: Evidence[str] = field(default_factory=Evidence.unknown)
    classification_date: Evidence[str] = field(default_factory=Evidence.unknown)
    review_status: Evidence[ReviewStatus] = field(default_factory=Evidence.unknown)
    classification_evidence: Evidence[ClassificationEvidence] = field(default_factory=Evidence.unknown)
    classification_independence: Evidence[bool] = field(default_factory=Evidence.unknown)
    disease_association: Evidence[bool] = field(default_factory=Evidence.unknown)
    molecular_mechanism: Evidence[MolecularMechanism] = field(default_factory=Evidence.unknown)
    classification_conflicts: Evidence[bool] = field(default_factory=Evidence.unknown)
    # Pairwise assertions may differ for each candidate/reference comparison.
    # A top-level assertion remains a backwards-compatible default when a
    # reference-specific assertion is absent.
    protein_substitution_equivalence: Evidence[bool] = field(default_factory=Evidence.unknown)
    underlying_allele_equivalence: Evidence[bool] = field(default_factory=Evidence.unknown)
    protein_effect_mechanism_compatibility: Evidence[bool] = field(default_factory=Evidence.unknown)
    splice_event_equivalence: Evidence[bool] = field(default_factory=Evidence.unknown)
    splice_motif_equivalence: Evidence[bool] = field(default_factory=Evidence.unknown)
    splice_prediction_strength_comparable: Evidence[bool] = field(default_factory=Evidence.unknown)
    splice_transcript_equivalence: Evidence[bool] = field(default_factory=Evidence.unknown)
    # Preserve whether a nested assertion was supplied even when its explicit
    # state is UNKNOWN.  That lets a reference-specific unknown override a
    # backwards-compatible top-level default instead of inheriting it.
    pairwise_assertions_supplied: frozenset[str] = field(default_factory=frozenset)

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any], *, flat: bool = False) -> "ReferenceVariant":
        prefix = "reference_" if flat else ""

        def take(name: str, *aliases: str) -> Any:
            for candidate in (f"{prefix}{name}", name, *aliases):
                if candidate in raw:
                    return raw[candidate]
            return _MISSING

        def supplied(name: str, *aliases: str) -> bool:
            return any(
                candidate in raw
                for candidate in (f"{prefix}{name}", name, *aliases)
            )

        pairwise_inputs = {
            "protein_substitution_equivalence": (),
            "underlying_allele_equivalence": (),
            "protein_effect_mechanism_compatibility": (),
            "splice_event_equivalence": (),
            "splice_motif_equivalence": (),
            "splice_prediction_strength_comparable": (),
            "splice_transcript_equivalence": ("splice_transcript_mapping_equivalence",),
        }

        return cls(
            annotations=VariantAnnotations.from_mapping(raw, prefix),
            gene_id=parse_evidence(take("gene_id"), parse_string, f"{prefix}gene_id"),
            disease_id=parse_evidence(take("disease_id"), parse_string, f"{prefix}disease_id"),
            classification=parse_evidence(
                take("classification"),
                _enum_parser(ReferenceClassification, _CLASSIFICATION_ALIASES),
                f"{prefix}classification",
            ),
            classification_source=parse_evidence(take("classification_source"), parse_string, f"{prefix}classification_source"),
            classification_date=parse_evidence(take("classification_date"), parse_string, f"{prefix}classification_date"),
            review_status=parse_evidence(take("classification_review_status", "review_status"), _enum_parser(ReviewStatus), f"{prefix}classification_review_status"),
            classification_evidence=parse_evidence(
                take("classification_evidence"),
                _enum_parser(
                    ClassificationEvidence,
                    {"CREDIBLE": ClassificationEvidence.ESTABLISHED, "ADEQUATE": ClassificationEvidence.ESTABLISHED},
                ),
                f"{prefix}classification_evidence",
            ),
            classification_independence=parse_evidence(
                take("classification_independence"), parse_bool, f"{prefix}classification_independence", boolean=True
            ),
            disease_association=parse_evidence(
                take("disease_association", "disease_match"), parse_bool, f"{prefix}disease_association", boolean=True
            ),
            molecular_mechanism=parse_evidence(
                take("molecular_mechanism"), _enum_parser(MolecularMechanism), f"{prefix}molecular_mechanism"
            ),
            classification_conflicts=parse_evidence(
                take("classification_conflicts"), parse_bool, f"{prefix}classification_conflicts", boolean=True
            ),
            protein_substitution_equivalence=parse_evidence(
                take("protein_substitution_equivalence"), parse_bool,
                f"{prefix}protein_substitution_equivalence", boolean=True,
            ),
            underlying_allele_equivalence=parse_evidence(
                take("underlying_allele_equivalence"), parse_bool,
                f"{prefix}underlying_allele_equivalence", boolean=True,
            ),
            protein_effect_mechanism_compatibility=parse_evidence(
                take("protein_effect_mechanism_compatibility"), parse_bool,
                f"{prefix}protein_effect_mechanism_compatibility", boolean=True,
            ),
            splice_event_equivalence=parse_evidence(
                take("splice_event_equivalence"), parse_bool,
                f"{prefix}splice_event_equivalence", boolean=True,
            ),
            splice_motif_equivalence=parse_evidence(
                take("splice_motif_equivalence"), parse_bool,
                f"{prefix}splice_motif_equivalence", boolean=True,
            ),
            splice_prediction_strength_comparable=parse_evidence(
                take("splice_prediction_strength_comparable"), parse_bool,
                f"{prefix}splice_prediction_strength_comparable", boolean=True,
            ),
            splice_transcript_equivalence=parse_evidence(
                take("splice_transcript_equivalence", "splice_transcript_mapping_equivalence"),
                parse_bool, f"{prefix}splice_transcript_equivalence", boolean=True,
            ),
            pairwise_assertions_supplied=frozenset(
                name
                for name, aliases in pairwise_inputs.items()
                if supplied(name, *aliases)
            ),
        )


def _raw_rule_list(raw: Any, field_name: str) -> tuple[Mapping[str, Any], ...]:
    values = raw if isinstance(raw, (list, tuple)) else (raw,)
    if not values:
        raise SchemaValidationError(f"{field_name} cannot be empty when supplied")
    result: list[Mapping[str, Any]] = []
    for value in values:
        if not isinstance(value, Mapping):
            raise SchemaValidationError(f"{field_name} must contain configuration objects")
        result.append(dict(value))
    return tuple(result)


@dataclass(frozen=True)
class BiologicalFeatures:
    """All PS1-relevant evidence for one candidate and zero or more references."""

    gene_id: Evidence[str] = field(default_factory=Evidence.unknown)
    disease_id: Evidence[str] = field(default_factory=Evidence.unknown)
    inheritance_context: Evidence[str] = field(default_factory=Evidence.unknown)
    candidate: VariantAnnotations = field(default_factory=VariantAnnotations)
    references: tuple[ReferenceVariant, ...] = ()
    protein_substitution_equivalence: Evidence[bool] = field(default_factory=Evidence.unknown)
    underlying_allele_equivalence: Evidence[bool] = field(default_factory=Evidence.unknown)
    protein_effect_mechanism_compatibility: Evidence[bool] = field(default_factory=Evidence.unknown)
    splice_event_equivalence: Evidence[bool] = field(default_factory=Evidence.unknown)
    splice_motif_equivalence: Evidence[bool] = field(default_factory=Evidence.unknown)
    splice_prediction_strength_comparable: Evidence[bool] = field(default_factory=Evidence.unknown)
    splice_transcript_equivalence: Evidence[bool] = field(default_factory=Evidence.unknown)
    splice_product_functional_consequence: Evidence[SpliceProductFunctionalConsequence] = field(default_factory=Evidence.unknown)
    splice_rna_assay_reliability: Evidence[RnaAssayReliability] = field(default_factory=Evidence.unknown)
    applicable_splicing_pvs1_strength: Evidence[SpliceBaselineCode] = field(default_factory=Evidence.unknown)
    independent_reference_classifications: Evidence[int] = field(default_factory=Evidence.unknown)
    vcep_specifications: Evidence[tuple[Mapping[str, Any], ...]] = field(default_factory=Evidence.unknown)
    vcep_gene_scope: Evidence[str] = field(default_factory=Evidence.unknown)
    vcep_disease_scope: Evidence[str] = field(default_factory=Evidence.unknown)
    vcep_specification_version: Evidence[str] = field(default_factory=Evidence.unknown)
    vcep_ps1_eligibility_rules: Evidence[tuple[Mapping[str, Any], ...]] = field(default_factory=Evidence.unknown)
    vcep_ps1_strength_rules: Evidence[tuple[Mapping[str, Any], ...]] = field(default_factory=Evidence.unknown)
    vcep_splicing_rules: Evidence[tuple[Mapping[str, Any], ...]] = field(default_factory=Evidence.unknown)
    vcep_reference_classification_requirements: Evidence[tuple[Mapping[str, Any], ...]] = field(default_factory=Evidence.unknown)
    vcep_exclusion_rules: Evidence[tuple[Mapping[str, Any], ...]] = field(default_factory=Evidence.unknown)
    feature_conflict_status: Evidence[bool] = field(default_factory=Evidence.unknown)

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "BiologicalFeatures":
        if not isinstance(raw, Mapping):
            raise SchemaValidationError("biological_features must be a mapping")

        aliases = {
            "reference_variants": "references",
            "applicable_vcep_ps1_specification": "vcep_specifications",
            "gene_specific_ps1_rules": "vcep_specifications",
            "candidate_protein_change": "candidate_protein_hgvs",
            "reference_protein_change": "reference_protein_hgvs",
            "candidate_nucleotide_change": "candidate_coding_hgvs",
            "reference_nucleotide_change": "reference_coding_hgvs",
            "protein_change_match": "protein_substitution_equivalence",
            "splice_event_comparability": "splice_event_equivalence",
            "splice_prediction_comparability": "splice_prediction_strength_comparable",
            "candidate_prediction_at_least_reference": "splice_prediction_strength_comparable",
            "splice_transcript_mapping_equivalence": "splice_transcript_equivalence",
        }
        normalized: dict[str, Any] = {}
        for name, value in raw.items():
            target = aliases.get(name, name)
            if target in normalized:
                raise SchemaValidationError(f"{name} and {target} both supply the same feature")
            normalized[target] = value

        reference_prefixes = (
            "reference_variant_id",
            "reference_gene_id",
            "reference_disease_id",
            "reference_genomic_allele",
            "reference_coding_hgvs",
            "reference_protein_hgvs",
            "reference_transcript_id_version",
            "reference_protein_accession_version",
            "reference_amino_acid_position",
            "reference_reference_amino_acid",
            "reference_alternate_amino_acid",
            "reference_variant_consequence",
            "reference_splice_",
            "reference_spliceai_",
            "reference_rna_",
            "reference_observed_",
            "reference_residual_",
            "reference_classification",
            "reference_disease_association",
            "reference_molecular_mechanism",
            "reference_protein_substitution_equivalence",
            "reference_underlying_allele_equivalence",
            "reference_protein_effect_mechanism_compatibility",
            "reference_splice_event_equivalence",
            "reference_splice_motif_equivalence",
            "reference_splice_prediction_strength_comparable",
            "reference_splice_transcript_equivalence",
        )
        flat_reference_keys = {
            key
            for key in normalized
            if any(key.startswith(prefix) for prefix in reference_prefixes)
        }
        references_raw = normalized.get("references", _MISSING)
        if references_raw is not _MISSING and flat_reference_keys:
            raise SchemaValidationError("use either references or flat reference_* fields, not both")
        if references_raw is _MISSING:
            references = (
                (ReferenceVariant.from_mapping(normalized, flat=True),) if flat_reference_keys else ()
            )
        else:
            if not isinstance(references_raw, (list, tuple)):
                raise SchemaValidationError("references must be a list of reference objects")
            references = tuple(
                ReferenceVariant.from_mapping(item)
                for item in references_raw
                if isinstance(item, Mapping)
            )
            if len(references) != len(references_raw):
                raise SchemaValidationError("references entries must be objects")

        parsers: Mapping[str, tuple[Callable[[Any, str], Any], bool]] = {
            "gene_id": (parse_string, False),
            "disease_id": (parse_string, False),
            "inheritance_context": (parse_string, False),
            "protein_substitution_equivalence": (parse_bool, True),
            "underlying_allele_equivalence": (parse_bool, True),
            "protein_effect_mechanism_compatibility": (parse_bool, True),
            "splice_event_equivalence": (parse_bool, True),
            "splice_motif_equivalence": (parse_bool, True),
            "splice_prediction_strength_comparable": (parse_bool, True),
            "splice_transcript_equivalence": (parse_bool, True),
            "splice_product_functional_consequence": (_enum_parser(SpliceProductFunctionalConsequence), False),
            "splice_rna_assay_reliability": (_enum_parser(RnaAssayReliability), False),
            "applicable_splicing_pvs1_strength": (_enum_parser(SpliceBaselineCode, _BASELINE_ALIASES), False),
            "independent_reference_classifications": (_positive_int, False),
            "vcep_specifications": (_raw_rule_list, False),
            "vcep_gene_scope": (parse_string, False),
            "vcep_disease_scope": (parse_string, False),
            "vcep_specification_version": (parse_string, False),
            "vcep_ps1_eligibility_rules": (_raw_rule_list, False),
            "vcep_ps1_strength_rules": (_raw_rule_list, False),
            "vcep_splicing_rules": (_raw_rule_list, False),
            "vcep_reference_classification_requirements": (_raw_rule_list, False),
            "vcep_exclusion_rules": (_raw_rule_list, False),
            "feature_conflict_status": (parse_bool, True),
        }

        recognised = set(parsers) | {"references"}
        candidate_stems = {
            "candidate_variant_id", "candidate_genomic_allele", "candidate_coding_hgvs",
            "candidate_protein_hgvs", "candidate_transcript_id_version",
            "candidate_protein_accession_version", "candidate_amino_acid_position",
            "candidate_reference_amino_acid", "candidate_alternate_amino_acid",
            "candidate_variant_consequence", "candidate_splice_impact", "candidate_splice_site_type",
            "candidate_splice_motif_position", "candidate_splice_prediction",
            "candidate_splice_motif_id", "candidate_splice_nucleotide_position",
            "candidate_spliceai_event_type", "candidate_spliceai_event_position",
            "candidate_splice_event", "candidate_splice_event_position",
            "candidate_splice_prediction_magnitude", "candidate_rna_splicing_evidence",
            "candidate_observed_splice_product", "candidate_residual_normal_splicing",
        }
        recognised |= candidate_stems | flat_reference_keys
        unknown = sorted(set(normalized) - recognised)
        if unknown:
            raise SchemaValidationError("unsupported biological feature(s): " + ", ".join(unknown))

        values = {
            name: parse_evidence(normalized.get(name, _MISSING), parser, name, boolean=boolean)
            for name, (parser, boolean) in parsers.items()
        }
        return cls(
            gene_id=values["gene_id"],
            disease_id=values["disease_id"],
            inheritance_context=values["inheritance_context"],
            candidate=VariantAnnotations.from_mapping(normalized, "candidate_"),
            references=references,
            protein_substitution_equivalence=values["protein_substitution_equivalence"],
            underlying_allele_equivalence=values["underlying_allele_equivalence"],
            protein_effect_mechanism_compatibility=values["protein_effect_mechanism_compatibility"],
            splice_event_equivalence=values["splice_event_equivalence"],
            splice_motif_equivalence=values["splice_motif_equivalence"],
            splice_prediction_strength_comparable=values["splice_prediction_strength_comparable"],
            splice_transcript_equivalence=values["splice_transcript_equivalence"],
            splice_product_functional_consequence=values["splice_product_functional_consequence"],
            splice_rna_assay_reliability=values["splice_rna_assay_reliability"],
            applicable_splicing_pvs1_strength=values["applicable_splicing_pvs1_strength"],
            independent_reference_classifications=values["independent_reference_classifications"],
            vcep_specifications=values["vcep_specifications"],
            vcep_gene_scope=values["vcep_gene_scope"],
            vcep_disease_scope=values["vcep_disease_scope"],
            vcep_specification_version=values["vcep_specification_version"],
            vcep_ps1_eligibility_rules=values["vcep_ps1_eligibility_rules"],
            vcep_ps1_strength_rules=values["vcep_ps1_strength_rules"],
            vcep_splicing_rules=values["vcep_splicing_rules"],
            vcep_reference_classification_requirements=values["vcep_reference_classification_requirements"],
            vcep_exclusion_rules=values["vcep_exclusion_rules"],
            feature_conflict_status=values["feature_conflict_status"],
        )

