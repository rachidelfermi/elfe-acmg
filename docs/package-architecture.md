# Package architecture

`acmg.classifier.classify_variant` accepts already extracted evidence in a
criterion keyed mapping. `acmg.criterion_registry` selects canonical evaluators
from `acmg.criteria`. `acmg.evidence_dependencies` resolves overlap before
`acmg.classification_rules` applies the existing ACMG/AMP combination rules.

The canonical criterion packages are `pvs1`, `ps1`, `ps2_pm6`, `ps3`, `ps4`,
`pm1` through `pm5`, `pp1` through `pp5`, `ba1`, `bs1` through `bs4`, and
`bp1` through `bp7`, all within `acmg.criteria`. The original 28 criteria are
represented. PS2 and PM6 use the one shared `ps2_pm6` path in the classifier;
PP5 and BP6 are represented but inactive.

The old `pvs1`, `ps2pm6`, `ps2_pm6`, and `adaptive_*` imports remain for
compatibility. Their Python modules delegate to `acmg.criteria` modules and
contain no scoring rules. The older standalone PM6 evaluator is retained at
`acmg.criteria.pm6` because its public result schema and verified HNF1A and
SLC9A6 VCEP pathways differ from the shared evaluator. It is not called by
the classifier. Converging the two PM6 interfaces requires a separate
scientific equivalence review; substituting one for the other here could
change clinical decisions.

The classifier supports ACMG/AMP 2015 combinations with the existing ClinGen
SVI PM2 update. A complete VCEP-wide or Bayesian final classification remains
outside its verified combination rules. Criterion-specific VCEP pathways stay
within their respective evaluator packages.
