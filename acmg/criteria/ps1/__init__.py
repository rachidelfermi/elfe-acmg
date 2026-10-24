"""Standalone adaptive ACMG/ClinGen PS1 evaluation package.

Only PS1 is evaluated.  The package does not classify a variant, retrieve
external records, annotate sequence data, or combine ACMG/AMP criteria.
"""

from .evaluator import PS1Engine, evaluate_ps1
from .evidence import Evidence, EvidenceKind, EvidenceState, Provenance, SchemaValidationError
from .guidelines import GENERIC_RULE_TABLE, rule_table_as_dicts
from .models import (
    BiologicalFeatures,
    ClassificationEvidence,
    DecisionStatus,
    MolecularMechanism,
    PS1Strength,
    ReferenceClassification,
    ReferenceVariant,
    ReviewStatus,
    RnaAssayReliability,
    SpliceBaselineCode,
    SpliceImpact,
    SpliceMotifPosition,
    SplicePrediction,
    SpliceProductFunctionalConsequence,
    SpliceEvidence,
    SpliceSiteType,
    VariantAnnotations,
    VariantConsequence,
)
from .output import GuidelineReference, PS1Decision

__all__ = [
    "BiologicalFeatures",
    "ClassificationEvidence",
    "DecisionStatus",
    "Evidence",
    "EvidenceKind",
    "EvidenceState",
    "GENERIC_RULE_TABLE",
    "GuidelineReference",
    "MolecularMechanism",
    "PS1Decision",
    "PS1Engine",
    "PS1Strength",
    "Provenance",
    "ReferenceClassification",
    "ReferenceVariant",
    "ReviewStatus",
    "RnaAssayReliability",
    "SchemaValidationError",
    "SpliceBaselineCode",
    "SpliceImpact",
    "SpliceMotifPosition",
    "SplicePrediction",
    "SpliceProductFunctionalConsequence",
    "SpliceEvidence",
    "SpliceSiteType",
    "VariantAnnotations",
    "VariantConsequence",
    "evaluate_ps1",
    "rule_table_as_dicts",
]

