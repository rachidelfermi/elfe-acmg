# ACMG/ClinGen engine architecture

This guide describes the repository as implemented. All examples are synthetic decision-engine inputs, not clinical assertions. The engine takes **already extracted biological evidence**, evaluates only supplied criteria, reconciles documented evidence reuse, and applies its supported ACMG/AMP classification combinations. It does not annotate variants, query databases, analyze patients, or establish clinical validity.

## 1. Project at a glance

The installed distribution is `adaptive-acmg-evidence` (Python 3.12 or newer); its primary import is `acmg`. All 28 original ACMG/AMP criteria are represented by 28 criterion packages. The registry runs 25 evaluator entries: 24 single active codes plus one shared `PS2_PM6` entry, which assigns **either** PS2 **or** PM6. PP5 and BP6 retain independent audit APIs but are inactive in unified classification under the [ClinGen SVI recommendation to discontinue reputable-source criteria](https://pubmed.ncbi.nlm.nih.gov/29543229/). The standalone `acmg.criteria.pm6` API is retained for its distinct, verified VCEP behavior; the unified classifier uses `acmg.criteria.ps2_pm6`.

Criterion evaluation returns a decision and evidence strength for **one criterion**. `acmg.classify_variant` combines eligible decisions into `PATHOGENIC`, `LIKELY_PATHOGENIC`, `VUS`, `LIKELY_BENIGN`, or `BENIGN`. It uses [ACMG/AMP 2015 Table 5](https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/) plus the implemented [ClinGen SVI PM2 Supporting update](https://www.clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf). It does not implement a general Bayesian or VCEP-wide final classifier.

## 2. Repository map and file responsibilities

The tree and inventories below cover **125 canonical Python source files**, **122 compatibility Python source files**, **70 test/support Python files**, the Markdown documentation, and root configuration. Braces in the compact tree list the actual filenames in a directory; they are not wildcard omissions. Generated caches and build outputs are listed separately after the source inventory.

### Actual source tree
```text
acmg/  (repository root)
|-- pyproject.toml, .gitignore, README.md, acmg_feature_to_rule_mapping.csv
|-- acmg/
|   |-- __init__.py, classification_rules.py, criterion_registry.py, evidence_dependencies.py
|   |-- classifier/__init__.py
|   `-- criteria/__init__.py
|       |-- ba1/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- bp1/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- bp2/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- bp3/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- bp4/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- bp5/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- bp6/{__init__.py, evaluator.py, models.py}
|       |-- bp7/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- bs1/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- bs2/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- bs3/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- bs4/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- pm1/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- pm2/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- pm3/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- pm4/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- pm5/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- pm6/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- pp1/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- pp2/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- pp3/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- pp4/{__init__.py, evaluator.py, guidelines.py, models.py}
|       |-- pp5/{__init__.py, evaluator.py, models.py}
|       |-- ps1/{__init__.py, evaluator.py, evidence.py, guidelines.py, models.py, output.py}
|       |   `-- rules/{__init__.py, common.py, missense.py, splicing.py, vcep.py}
|       |-- ps2_pm6/{__init__.py, evaluator.py, models.py, specifications.py}
|       |-- ps3/{__init__.py, evaluator.py, functional_assays.py, guidelines.py, models.py}
|       |-- ps4/{__init__.py, evaluator.py, guidelines.py, models.py, statistics.py}
|       |-- pvs1/{__init__.py, evaluator.py, models.py, specifications.py}
|-- adaptive_ba1/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_bp1/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_bp2/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_bp3/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_bp4/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_bp5/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_bp6/{__init__.py, evaluator.py, models.py}
|-- adaptive_bp7/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_bs1/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_bs2/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_bs3/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_bs4/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_pm1/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_pm2/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_pm3/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_pm4/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_pm5/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_pm6/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_pp1/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_pp2/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_pp3/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_pp4/{__init__.py, evaluator.py, guidelines.py, models.py}
|-- adaptive_pp5/{__init__.py, evaluator.py, models.py}
|-- adaptive_ps1/{__init__.py, evaluator.py, evidence.py, guidelines.py, models.py, output.py, rules/__init__.py, rules/common.py, rules/missense.py, rules/splicing.py, rules/vcep.py}
|-- adaptive_ps2_pm6/{__init__.py}
|-- adaptive_ps2pm6/{__init__.py}
|-- adaptive_ps3/{__init__.py, evaluator.py, functional_assays.py, guidelines.py, models.py}
|-- adaptive_ps4/{__init__.py, evaluator.py, guidelines.py, models.py, statistics.py}
|-- ps2_pm6/{__init__.py}
|-- ps2pm6/{__init__.py, evaluator.py, models.py, specifications.py}
|-- pvs1/{__init__.py, evaluator.py, models.py, specifications.py}
|-- tests/
|   |-- __init__.py
|   |-- pm3_fixtures.py
|   |-- pp3_fixtures.py
|   |-- ps3_fixtures.py
|   |-- test_adaptive_pm6.py
|   |-- test_adaptive_ps1_engine.py
|   |-- test_assay_validation.py
|   |-- test_ba1.py
|   |-- test_bp1_eligibility.py
|   |-- test_bp1_exceptions.py
|   |-- test_bp1_mechanism.py
|   |-- test_bp1_missing_evidence.py
|   |-- test_bp1_vcep.py
|   |-- test_bp2.py
|   |-- test_bp3_functional_regions.py
|   |-- test_bp3_missing_evidence.py
|   |-- test_bp3_repeat_indels.py
|   |-- test_bp3_vcep.py
|   |-- test_bp4_calibration.py
|   |-- test_bp4_missense.py
|   |-- test_bp4_splicing.py
|   |-- test_bp4_vcep.py
|   |-- test_bp5_engine.py
|   |-- test_bp6.py
|   |-- test_bp7_engine.py
|   |-- test_bs1.py
|   |-- test_bs2.py
|   |-- test_bs3_calibration.py
|   |-- test_bs3_functional_assays.py
|   |-- test_bs3_missing_evidence.py
|   |-- test_bs3_splicing.py
|   |-- test_bs3_vcep.py
|   |-- test_bs4.py
|   |-- test_calibration.py
|   |-- test_classification.py
|   |-- test_domains.py
|   |-- test_evidence_dependencies.py
|   |-- test_homozygous.py
|   |-- test_hotspots.py
|   |-- test_inframe_indels.py
|   |-- test_missense.py
|   |-- test_missing_evidence.py
|   |-- test_multiple_cases.py
|   |-- test_pm1_missing_evidence.py
|   |-- test_pm1_vcep.py
|   |-- test_pm2.py
|   |-- test_pm3_missing_evidence.py
|   |-- test_pm4_missing_evidence.py
|   |-- test_pm5_engine.py
|   |-- test_pp1.py
|   |-- test_pp2.py
|   |-- test_pp3_calibration.py
|   |-- test_pp3_missense.py
|   |-- test_pp3_missing_evidence.py
|   |-- test_pp3_splicing.py
|   |-- test_pp3_vcep.py
|   |-- test_pp4_engine.py
|   |-- test_pp5.py
|   |-- test_ps2pm6_engine.py
|   |-- test_ps3.py
|   |-- test_ps4.py
|   |-- test_pvs1_engine.py
|   |-- test_repeat_regions.py
|   |-- test_splicing.py
|   |-- test_stop_loss.py
|   |-- test_trans.py
|   |-- test_unified_classification.py
|   |-- test_unphased.py
|   |-- test_vcep.py
|   |-- test_vcep_pp3.py
`-- docs/
    |-- bp1-guideline-matrix.md
    |-- BP4_SOURCE_RULE_MATRIX.md
    |-- bp7-guideline-matrix.md
    |-- bs2-guideline-matrix.md
    |-- bs3-guideline-matrix.md
    |-- bs4-guideline-matrix.md
    |-- coverage-matrix.md
    |-- decision-table.md
    |-- feature-contract.md
    |-- package-architecture.md
    |-- pm1-guideline-matrix.md
    |-- pm2-guideline-matrix.md
    |-- pm4-guideline-matrix.md
    |-- pm5-guideline-matrix.md
    |-- pp1-guideline-matrix.md
    |-- pp2-guideline-matrix.md
    |-- pp3-guideline-matrix.md
    |-- pp5-policy.md
    |-- PROJECT_ARCHITECTURE.md
    |-- ps1-decision-table.md
    |-- ps2-pm6-decision-table.md
    |-- ps3-decision-table.md
    |-- ps4-guideline-matrix.md
    |-- unified-classification-source-matrix.md
    |-- unified-scientific-validation.md
|-- adaptive_ba1/README.md  (legacy package documentation)
|-- adaptive_bp2/README.md  (legacy package documentation)
|-- adaptive_bp3/README.md  (legacy package documentation)
|-- adaptive_bp5/README.md  (legacy package documentation)
|-- adaptive_bp6/README.md  (legacy package documentation)
|-- adaptive_bs1/README.md  (legacy package documentation)
|-- adaptive_bs4/README.md  (legacy package documentation)
|-- adaptive_pm3/SOURCES.md  (legacy package documentation)
|-- adaptive_pm6/README.md  (legacy package documentation)
|-- adaptive_pm6/SOURCES.md  (legacy package documentation)
|-- adaptive_pp4/README.md  (legacy package documentation)
```

### Root and integration files
| File | Responsibility | Used by |
| --- | --- | --- |
| `acmg/__init__.py` | Exports `classify_variant`. | External callers |
| `acmg/classifier/__init__.py` | Validates explicit IDs/context; dispatches supplied criteria; normalizes rows; resolves dependencies; combines evidence; returns a dictionary. | Public API |
| `acmg/criterion_registry.py` | Maps 25 evaluator entries to canonical package functions; lists disabled codes and accepted de novo aliases. | Classifier |
| `acmg/evidence_dependencies.py` | Marks documented overlapping activated rows uncounted and records contradictions. | Classifier |
| `acmg/classification_rules.py` | ACMG/AMP Table 5 combinations and implemented SVI PM2 update. | Classifier |
| `acmg/criteria/__init__.py` | Marks the canonical criterion namespace. | Criterion imports |
### Canonical criterion source files
| File | Responsibility | Used by |
| --- | --- | --- |
| `acmg/criteria/ba1/__init__.py` | Exports the public BA1 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/ba1/evaluator.py` | BA1 scientific decision pathways and strength assignment. | BA1 public package |
| `acmg/criteria/ba1/guidelines.py` | BA1 rule thresholds, references, VCEP scope and guideline selection. | BA1 evaluator |
| `acmg/criteria/ba1/models.py` | BA1 input, evidence-state and decision/result models. | BA1 evaluator and callers |
| `acmg/criteria/bp1/__init__.py` | Exports the public BP1 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/bp1/evaluator.py` | BP1 scientific decision pathways and strength assignment. | BP1 public package |
| `acmg/criteria/bp1/guidelines.py` | BP1 rule thresholds, references, VCEP scope and guideline selection. | BP1 evaluator |
| `acmg/criteria/bp1/models.py` | BP1 input, evidence-state and decision/result models. | BP1 evaluator and callers |
| `acmg/criteria/bp2/__init__.py` | Exports the public BP2 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/bp2/evaluator.py` | BP2 scientific decision pathways and strength assignment. | BP2 public package |
| `acmg/criteria/bp2/guidelines.py` | BP2 rule thresholds, references, VCEP scope and guideline selection. | BP2 evaluator |
| `acmg/criteria/bp2/models.py` | BP2 input, evidence-state and decision/result models. | BP2 evaluator and callers |
| `acmg/criteria/bp3/__init__.py` | Exports the public BP3 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/bp3/evaluator.py` | BP3 scientific decision pathways and strength assignment. | BP3 public package |
| `acmg/criteria/bp3/guidelines.py` | BP3 rule thresholds, references, VCEP scope and guideline selection. | BP3 evaluator |
| `acmg/criteria/bp3/models.py` | BP3 input, evidence-state and decision/result models. | BP3 evaluator and callers |
| `acmg/criteria/bp4/__init__.py` | Exports the public BP4 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/bp4/evaluator.py` | BP4 scientific decision pathways and strength assignment. | BP4 public package |
| `acmg/criteria/bp4/guidelines.py` | BP4 rule thresholds, references, VCEP scope and guideline selection. | BP4 evaluator |
| `acmg/criteria/bp4/models.py` | BP4 input, evidence-state and decision/result models. | BP4 evaluator and callers |
| `acmg/criteria/bp5/__init__.py` | Exports the public BP5 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/bp5/evaluator.py` | BP5 scientific decision pathways and strength assignment. | BP5 public package |
| `acmg/criteria/bp5/guidelines.py` | BP5 rule thresholds, references, VCEP scope and guideline selection. | BP5 evaluator |
| `acmg/criteria/bp5/models.py` | BP5 input, evidence-state and decision/result models. | BP5 evaluator and callers |
| `acmg/criteria/bp6/__init__.py` | Exports the public BP6 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/bp6/evaluator.py` | BP6 scientific decision pathways and strength assignment. | BP6 public package |
| `acmg/criteria/bp6/models.py` | BP6 input, evidence-state and decision/result models. | BP6 evaluator and callers |
| `acmg/criteria/bp7/__init__.py` | Exports the public BP7 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/bp7/evaluator.py` | BP7 scientific decision pathways and strength assignment. | BP7 public package |
| `acmg/criteria/bp7/guidelines.py` | BP7 rule thresholds, references, VCEP scope and guideline selection. | BP7 evaluator |
| `acmg/criteria/bp7/models.py` | BP7 input, evidence-state and decision/result models. | BP7 evaluator and callers |
| `acmg/criteria/bs1/__init__.py` | Exports the public BS1 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/bs1/evaluator.py` | BS1 scientific decision pathways and strength assignment. | BS1 public package |
| `acmg/criteria/bs1/guidelines.py` | BS1 rule thresholds, references, VCEP scope and guideline selection. | BS1 evaluator |
| `acmg/criteria/bs1/models.py` | BS1 input, evidence-state and decision/result models. | BS1 evaluator and callers |
| `acmg/criteria/bs2/__init__.py` | Exports the public BS2 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/bs2/evaluator.py` | BS2 scientific decision pathways and strength assignment. | BS2 public package |
| `acmg/criteria/bs2/guidelines.py` | BS2 rule thresholds, references, VCEP scope and guideline selection. | BS2 evaluator |
| `acmg/criteria/bs2/models.py` | BS2 input, evidence-state and decision/result models. | BS2 evaluator and callers |
| `acmg/criteria/bs3/__init__.py` | Exports the public BS3 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/bs3/evaluator.py` | BS3 scientific decision pathways and strength assignment. | BS3 public package |
| `acmg/criteria/bs3/guidelines.py` | BS3 rule thresholds, references, VCEP scope and guideline selection. | BS3 evaluator |
| `acmg/criteria/bs3/models.py` | BS3 input, evidence-state and decision/result models. | BS3 evaluator and callers |
| `acmg/criteria/bs4/__init__.py` | Exports the public BS4 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/bs4/evaluator.py` | BS4 scientific decision pathways and strength assignment. | BS4 public package |
| `acmg/criteria/bs4/guidelines.py` | BS4 rule thresholds, references, VCEP scope and guideline selection. | BS4 evaluator |
| `acmg/criteria/bs4/models.py` | BS4 input, evidence-state and decision/result models. | BS4 evaluator and callers |
| `acmg/criteria/pm1/__init__.py` | Exports the public PM1 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/pm1/evaluator.py` | PM1 scientific decision pathways and strength assignment. | PM1 public package |
| `acmg/criteria/pm1/guidelines.py` | PM1 rule thresholds, references, VCEP scope and guideline selection. | PM1 evaluator |
| `acmg/criteria/pm1/models.py` | PM1 input, evidence-state and decision/result models. | PM1 evaluator and callers |
| `acmg/criteria/pm2/__init__.py` | Exports the public PM2 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/pm2/evaluator.py` | PM2 scientific decision pathways and strength assignment. | PM2 public package |
| `acmg/criteria/pm2/guidelines.py` | PM2 rule thresholds, references, VCEP scope and guideline selection. | PM2 evaluator |
| `acmg/criteria/pm2/models.py` | PM2 input, evidence-state and decision/result models. | PM2 evaluator and callers |
| `acmg/criteria/pm3/__init__.py` | Exports the public PM3 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/pm3/evaluator.py` | PM3 scientific decision pathways and strength assignment. | PM3 public package |
| `acmg/criteria/pm3/guidelines.py` | PM3 rule thresholds, references, VCEP scope and guideline selection. | PM3 evaluator |
| `acmg/criteria/pm3/models.py` | PM3 input, evidence-state and decision/result models. | PM3 evaluator and callers |
| `acmg/criteria/pm4/__init__.py` | Exports the public PM4 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/pm4/evaluator.py` | PM4 scientific decision pathways and strength assignment. | PM4 public package |
| `acmg/criteria/pm4/guidelines.py` | PM4 rule thresholds, references, VCEP scope and guideline selection. | PM4 evaluator |
| `acmg/criteria/pm4/models.py` | PM4 input, evidence-state and decision/result models. | PM4 evaluator and callers |
| `acmg/criteria/pm5/__init__.py` | Exports the public PM5 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/pm5/evaluator.py` | PM5 scientific decision pathways and strength assignment. | PM5 public package |
| `acmg/criteria/pm5/guidelines.py` | PM5 rule thresholds, references, VCEP scope and guideline selection. | PM5 evaluator |
| `acmg/criteria/pm5/models.py` | PM5 input, evidence-state and decision/result models. | PM5 evaluator and callers |
| `acmg/criteria/pm6/__init__.py` | Exports the public PM6 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/pm6/evaluator.py` | PM6 scientific decision pathways and strength assignment. | PM6 public package |
| `acmg/criteria/pm6/guidelines.py` | PM6 rule thresholds, references, VCEP scope and guideline selection. | PM6 evaluator |
| `acmg/criteria/pm6/models.py` | PM6 input, evidence-state and decision/result models. | PM6 evaluator and callers |
| `acmg/criteria/pp1/__init__.py` | Exports the public PP1 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/pp1/evaluator.py` | PP1 scientific decision pathways and strength assignment. | PP1 public package |
| `acmg/criteria/pp1/guidelines.py` | PP1 rule thresholds, references, VCEP scope and guideline selection. | PP1 evaluator |
| `acmg/criteria/pp1/models.py` | PP1 input, evidence-state and decision/result models. | PP1 evaluator and callers |
| `acmg/criteria/pp2/__init__.py` | Exports the public PP2 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/pp2/evaluator.py` | PP2 scientific decision pathways and strength assignment. | PP2 public package |
| `acmg/criteria/pp2/guidelines.py` | PP2 rule thresholds, references, VCEP scope and guideline selection. | PP2 evaluator |
| `acmg/criteria/pp2/models.py` | PP2 input, evidence-state and decision/result models. | PP2 evaluator and callers |
| `acmg/criteria/pp3/__init__.py` | Exports the public PP3 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/pp3/evaluator.py` | PP3 scientific decision pathways and strength assignment. | PP3 public package |
| `acmg/criteria/pp3/guidelines.py` | PP3 rule thresholds, references, VCEP scope and guideline selection. | PP3 evaluator |
| `acmg/criteria/pp3/models.py` | PP3 input, evidence-state and decision/result models. | PP3 evaluator and callers |
| `acmg/criteria/pp4/__init__.py` | Exports the public PP4 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/pp4/evaluator.py` | PP4 scientific decision pathways and strength assignment. | PP4 public package |
| `acmg/criteria/pp4/guidelines.py` | PP4 rule thresholds, references, VCEP scope and guideline selection. | PP4 evaluator |
| `acmg/criteria/pp4/models.py` | PP4 input, evidence-state and decision/result models. | PP4 evaluator and callers |
| `acmg/criteria/pp5/__init__.py` | Exports the public PP5 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/pp5/evaluator.py` | PP5 scientific decision pathways and strength assignment. | PP5 public package |
| `acmg/criteria/pp5/models.py` | PP5 input, evidence-state and decision/result models. | PP5 evaluator and callers |
| `acmg/criteria/ps1/__init__.py` | Exports the public PS1 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/ps1/evaluator.py` | PS1 scientific decision pathways and strength assignment. | PS1 public package |
| `acmg/criteria/ps1/evidence.py` | Typed evidence states, provenance and parsing helpers for PS1. | PS1 public package |
| `acmg/criteria/ps1/guidelines.py` | PS1 rule thresholds, references, VCEP scope and guideline selection. | PS1 evaluator |
| `acmg/criteria/ps1/models.py` | PS1 input, evidence-state and decision/result models. | PS1 evaluator and callers |
| `acmg/criteria/ps1/output.py` | PS1 decision and guideline output models. | PS1 public package |
| `acmg/criteria/ps1/rules/__init__.py` | PS1 rule subpackage marker. | Registry or independent caller |
| `acmg/criteria/ps1/rules/common.py` | Shared PS1 preconditions and rule-result helpers. | PS1 public package |
| `acmg/criteria/ps1/rules/missense.py` | Same-amino-acid protein comparison pathway. | PS1 public package |
| `acmg/criteria/ps1/rules/splicing.py` | Splice-event comparison pathway. | PS1 public package |
| `acmg/criteria/ps1/rules/vcep.py` | Scoped, versioned PS1 VCEP rule selection and evaluation. | PS1 public package |
| `acmg/criteria/ps2_pm6/__init__.py` | Exports the public PS2_PM6 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/ps2_pm6/evaluator.py` | PS2_PM6 scientific decision pathways and strength assignment. | PS2_PM6 public package |
| `acmg/criteria/ps2_pm6/models.py` | PS2_PM6 input, evidence-state and decision/result models. | PS2_PM6 evaluator and callers |
| `acmg/criteria/ps2_pm6/specifications.py` | ClinGen SVI and scoped de novo specification tables. | PS2_PM6 evaluator |
| `acmg/criteria/ps3/__init__.py` | Exports the public PS3 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/ps3/evaluator.py` | PS3 scientific decision pathways and strength assignment. | PS3 public package |
| `acmg/criteria/ps3/functional_assays.py` | PS3 assay calibration, relevance and individual assay assessment. | PS3 evaluator |
| `acmg/criteria/ps3/guidelines.py` | PS3 rule thresholds, references, VCEP scope and guideline selection. | PS3 evaluator |
| `acmg/criteria/ps3/models.py` | PS3 input, evidence-state and decision/result models. | PS3 evaluator and callers |
| `acmg/criteria/ps4/__init__.py` | Exports the public PS4 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/ps4/evaluator.py` | PS4 scientific decision pathways and strength assignment. | PS4 public package |
| `acmg/criteria/ps4/guidelines.py` | PS4 rule thresholds, references, VCEP scope and guideline selection. | PS4 evaluator |
| `acmg/criteria/ps4/models.py` | PS4 input, evidence-state and decision/result models. | PS4 evaluator and callers |
| `acmg/criteria/ps4/statistics.py` | Case-control statistics used by PS4. | PS4 evaluator |
| `acmg/criteria/pvs1/__init__.py` | Exports the public PVS1 evaluator and selected native types. | Registry or independent caller |
| `acmg/criteria/pvs1/evaluator.py` | PVS1 scientific decision pathways and strength assignment. | PVS1 public package |
| `acmg/criteria/pvs1/models.py` | PVS1 input, evidence-state and decision/result models. | PVS1 evaluator and callers |
| `acmg/criteria/pvs1/specifications.py` | PVS1 guideline decision trees and pinned specification handling. | PVS1 evaluator |
### Compatibility Python source files
Each listed Python file is a thin import delegate to the matching canonical module; no legacy directory retains criterion scoring code. The exceptions in route are `adaptive_pm6` ? `acmg.criteria.pm6` and the three extra shared de novo alias packages ? `acmg.criteria.ps2_pm6`.
| File | Responsibility | Used by |
| --- | --- | --- |
| `adaptive_ba1/__init__.py` | Re-exports `acmg.criteria.ba1`. | Legacy imports |
| `adaptive_ba1/evaluator.py` | Re-exports `acmg.criteria.ba1.evaluator`. | Legacy imports |
| `adaptive_ba1/guidelines.py` | Re-exports `acmg.criteria.ba1.guidelines`. | Legacy imports |
| `adaptive_ba1/models.py` | Re-exports `acmg.criteria.ba1.models`. | Legacy imports |
| `adaptive_bp1/__init__.py` | Re-exports `acmg.criteria.bp1`. | Legacy imports |
| `adaptive_bp1/evaluator.py` | Re-exports `acmg.criteria.bp1.evaluator`. | Legacy imports |
| `adaptive_bp1/guidelines.py` | Re-exports `acmg.criteria.bp1.guidelines`. | Legacy imports |
| `adaptive_bp1/models.py` | Re-exports `acmg.criteria.bp1.models`. | Legacy imports |
| `adaptive_bp2/__init__.py` | Re-exports `acmg.criteria.bp2`. | Legacy imports |
| `adaptive_bp2/evaluator.py` | Re-exports `acmg.criteria.bp2.evaluator`. | Legacy imports |
| `adaptive_bp2/guidelines.py` | Re-exports `acmg.criteria.bp2.guidelines`. | Legacy imports |
| `adaptive_bp2/models.py` | Re-exports `acmg.criteria.bp2.models`. | Legacy imports |
| `adaptive_bp3/__init__.py` | Re-exports `acmg.criteria.bp3`. | Legacy imports |
| `adaptive_bp3/evaluator.py` | Re-exports `acmg.criteria.bp3.evaluator`. | Legacy imports |
| `adaptive_bp3/guidelines.py` | Re-exports `acmg.criteria.bp3.guidelines`. | Legacy imports |
| `adaptive_bp3/models.py` | Re-exports `acmg.criteria.bp3.models`. | Legacy imports |
| `adaptive_bp4/__init__.py` | Re-exports `acmg.criteria.bp4`. | Legacy imports |
| `adaptive_bp4/evaluator.py` | Re-exports `acmg.criteria.bp4.evaluator`. | Legacy imports |
| `adaptive_bp4/guidelines.py` | Re-exports `acmg.criteria.bp4.guidelines`. | Legacy imports |
| `adaptive_bp4/models.py` | Re-exports `acmg.criteria.bp4.models`. | Legacy imports |
| `adaptive_bp5/__init__.py` | Re-exports `acmg.criteria.bp5`. | Legacy imports |
| `adaptive_bp5/evaluator.py` | Re-exports `acmg.criteria.bp5.evaluator`. | Legacy imports |
| `adaptive_bp5/guidelines.py` | Re-exports `acmg.criteria.bp5.guidelines`. | Legacy imports |
| `adaptive_bp5/models.py` | Re-exports `acmg.criteria.bp5.models`. | Legacy imports |
| `adaptive_bp6/__init__.py` | Re-exports `acmg.criteria.bp6`. | Legacy imports |
| `adaptive_bp6/evaluator.py` | Re-exports `acmg.criteria.bp6.evaluator`. | Legacy imports |
| `adaptive_bp6/models.py` | Re-exports `acmg.criteria.bp6.models`. | Legacy imports |
| `adaptive_bp7/__init__.py` | Re-exports `acmg.criteria.bp7`. | Legacy imports |
| `adaptive_bp7/evaluator.py` | Re-exports `acmg.criteria.bp7.evaluator`. | Legacy imports |
| `adaptive_bp7/guidelines.py` | Re-exports `acmg.criteria.bp7.guidelines`. | Legacy imports |
| `adaptive_bp7/models.py` | Re-exports `acmg.criteria.bp7.models`. | Legacy imports |
| `adaptive_bs1/__init__.py` | Re-exports `acmg.criteria.bs1`. | Legacy imports |
| `adaptive_bs1/evaluator.py` | Re-exports `acmg.criteria.bs1.evaluator`. | Legacy imports |
| `adaptive_bs1/guidelines.py` | Re-exports `acmg.criteria.bs1.guidelines`. | Legacy imports |
| `adaptive_bs1/models.py` | Re-exports `acmg.criteria.bs1.models`. | Legacy imports |
| `adaptive_bs2/__init__.py` | Re-exports `acmg.criteria.bs2`. | Legacy imports |
| `adaptive_bs2/evaluator.py` | Re-exports `acmg.criteria.bs2.evaluator`. | Legacy imports |
| `adaptive_bs2/guidelines.py` | Re-exports `acmg.criteria.bs2.guidelines`. | Legacy imports |
| `adaptive_bs2/models.py` | Re-exports `acmg.criteria.bs2.models`. | Legacy imports |
| `adaptive_bs3/__init__.py` | Re-exports `acmg.criteria.bs3`. | Legacy imports |
| `adaptive_bs3/evaluator.py` | Re-exports `acmg.criteria.bs3.evaluator`. | Legacy imports |
| `adaptive_bs3/guidelines.py` | Re-exports `acmg.criteria.bs3.guidelines`. | Legacy imports |
| `adaptive_bs3/models.py` | Re-exports `acmg.criteria.bs3.models`. | Legacy imports |
| `adaptive_bs4/__init__.py` | Re-exports `acmg.criteria.bs4`. | Legacy imports |
| `adaptive_bs4/evaluator.py` | Re-exports `acmg.criteria.bs4.evaluator`. | Legacy imports |
| `adaptive_bs4/guidelines.py` | Re-exports `acmg.criteria.bs4.guidelines`. | Legacy imports |
| `adaptive_bs4/models.py` | Re-exports `acmg.criteria.bs4.models`. | Legacy imports |
| `adaptive_pm1/__init__.py` | Re-exports `acmg.criteria.pm1`. | Legacy imports |
| `adaptive_pm1/evaluator.py` | Re-exports `acmg.criteria.pm1.evaluator`. | Legacy imports |
| `adaptive_pm1/guidelines.py` | Re-exports `acmg.criteria.pm1.guidelines`. | Legacy imports |
| `adaptive_pm1/models.py` | Re-exports `acmg.criteria.pm1.models`. | Legacy imports |
| `adaptive_pm2/__init__.py` | Re-exports `acmg.criteria.pm2`. | Legacy imports |
| `adaptive_pm2/evaluator.py` | Re-exports `acmg.criteria.pm2.evaluator`. | Legacy imports |
| `adaptive_pm2/guidelines.py` | Re-exports `acmg.criteria.pm2.guidelines`. | Legacy imports |
| `adaptive_pm2/models.py` | Re-exports `acmg.criteria.pm2.models`. | Legacy imports |
| `adaptive_pm3/__init__.py` | Re-exports `acmg.criteria.pm3`. | Legacy imports |
| `adaptive_pm3/evaluator.py` | Re-exports `acmg.criteria.pm3.evaluator`. | Legacy imports |
| `adaptive_pm3/guidelines.py` | Re-exports `acmg.criteria.pm3.guidelines`. | Legacy imports |
| `adaptive_pm3/models.py` | Re-exports `acmg.criteria.pm3.models`. | Legacy imports |
| `adaptive_pm4/__init__.py` | Re-exports `acmg.criteria.pm4`. | Legacy imports |
| `adaptive_pm4/evaluator.py` | Re-exports `acmg.criteria.pm4.evaluator`. | Legacy imports |
| `adaptive_pm4/guidelines.py` | Re-exports `acmg.criteria.pm4.guidelines`. | Legacy imports |
| `adaptive_pm4/models.py` | Re-exports `acmg.criteria.pm4.models`. | Legacy imports |
| `adaptive_pm5/__init__.py` | Re-exports `acmg.criteria.pm5`. | Legacy imports |
| `adaptive_pm5/evaluator.py` | Re-exports `acmg.criteria.pm5.evaluator`. | Legacy imports |
| `adaptive_pm5/guidelines.py` | Re-exports `acmg.criteria.pm5.guidelines`. | Legacy imports |
| `adaptive_pm5/models.py` | Re-exports `acmg.criteria.pm5.models`. | Legacy imports |
| `adaptive_pm6/__init__.py` | Re-exports `acmg.criteria.pm6`. | Legacy imports |
| `adaptive_pm6/evaluator.py` | Re-exports `acmg.criteria.pm6.evaluator`. | Legacy imports |
| `adaptive_pm6/guidelines.py` | Re-exports `acmg.criteria.pm6.guidelines`. | Legacy imports |
| `adaptive_pm6/models.py` | Re-exports `acmg.criteria.pm6.models`. | Legacy imports |
| `adaptive_pp1/__init__.py` | Re-exports `acmg.criteria.pp1`. | Legacy imports |
| `adaptive_pp1/evaluator.py` | Re-exports `acmg.criteria.pp1.evaluator`. | Legacy imports |
| `adaptive_pp1/guidelines.py` | Re-exports `acmg.criteria.pp1.guidelines`. | Legacy imports |
| `adaptive_pp1/models.py` | Re-exports `acmg.criteria.pp1.models`. | Legacy imports |
| `adaptive_pp2/__init__.py` | Re-exports `acmg.criteria.pp2`. | Legacy imports |
| `adaptive_pp2/evaluator.py` | Re-exports `acmg.criteria.pp2.evaluator`. | Legacy imports |
| `adaptive_pp2/guidelines.py` | Re-exports `acmg.criteria.pp2.guidelines`. | Legacy imports |
| `adaptive_pp2/models.py` | Re-exports `acmg.criteria.pp2.models`. | Legacy imports |
| `adaptive_pp3/__init__.py` | Re-exports `acmg.criteria.pp3`. | Legacy imports |
| `adaptive_pp3/evaluator.py` | Re-exports `acmg.criteria.pp3.evaluator`. | Legacy imports |
| `adaptive_pp3/guidelines.py` | Re-exports `acmg.criteria.pp3.guidelines`. | Legacy imports |
| `adaptive_pp3/models.py` | Re-exports `acmg.criteria.pp3.models`. | Legacy imports |
| `adaptive_pp4/__init__.py` | Re-exports `acmg.criteria.pp4`. | Legacy imports |
| `adaptive_pp4/evaluator.py` | Re-exports `acmg.criteria.pp4.evaluator`. | Legacy imports |
| `adaptive_pp4/guidelines.py` | Re-exports `acmg.criteria.pp4.guidelines`. | Legacy imports |
| `adaptive_pp4/models.py` | Re-exports `acmg.criteria.pp4.models`. | Legacy imports |
| `adaptive_pp5/__init__.py` | Re-exports `acmg.criteria.pp5`. | Legacy imports |
| `adaptive_pp5/evaluator.py` | Re-exports `acmg.criteria.pp5.evaluator`. | Legacy imports |
| `adaptive_pp5/models.py` | Re-exports `acmg.criteria.pp5.models`. | Legacy imports |
| `adaptive_ps1/__init__.py` | Re-exports `acmg.criteria.ps1`. | Legacy imports |
| `adaptive_ps1/evaluator.py` | Re-exports `acmg.criteria.ps1.evaluator`. | Legacy imports |
| `adaptive_ps1/evidence.py` | Re-exports `acmg.criteria.ps1.evidence`. | Legacy imports |
| `adaptive_ps1/guidelines.py` | Re-exports `acmg.criteria.ps1.guidelines`. | Legacy imports |
| `adaptive_ps1/models.py` | Re-exports `acmg.criteria.ps1.models`. | Legacy imports |
| `adaptive_ps1/output.py` | Re-exports `acmg.criteria.ps1.output`. | Legacy imports |
| `adaptive_ps1/rules/__init__.py` | Re-exports `acmg.criteria.ps1.rules`. | Legacy imports |
| `adaptive_ps1/rules/common.py` | Re-exports `acmg.criteria.ps1.rules.common`. | Legacy imports |
| `adaptive_ps1/rules/missense.py` | Re-exports `acmg.criteria.ps1.rules.missense`. | Legacy imports |
| `adaptive_ps1/rules/splicing.py` | Re-exports `acmg.criteria.ps1.rules.splicing`. | Legacy imports |
| `adaptive_ps1/rules/vcep.py` | Re-exports `acmg.criteria.ps1.rules.vcep`. | Legacy imports |
| `adaptive_ps2_pm6/__init__.py` | Re-exports `acmg.criteria.ps2_pm6`. | Legacy imports |
| `adaptive_ps2pm6/__init__.py` | Re-exports `acmg.criteria.ps2_pm6`. | Legacy imports |
| `adaptive_ps3/__init__.py` | Re-exports `acmg.criteria.ps3`. | Legacy imports |
| `adaptive_ps3/evaluator.py` | Re-exports `acmg.criteria.ps3.evaluator`. | Legacy imports |
| `adaptive_ps3/functional_assays.py` | Re-exports `acmg.criteria.ps3.functional_assays`. | Legacy imports |
| `adaptive_ps3/guidelines.py` | Re-exports `acmg.criteria.ps3.guidelines`. | Legacy imports |
| `adaptive_ps3/models.py` | Re-exports `acmg.criteria.ps3.models`. | Legacy imports |
| `adaptive_ps4/__init__.py` | Re-exports `acmg.criteria.ps4`. | Legacy imports |
| `adaptive_ps4/evaluator.py` | Re-exports `acmg.criteria.ps4.evaluator`. | Legacy imports |
| `adaptive_ps4/guidelines.py` | Re-exports `acmg.criteria.ps4.guidelines`. | Legacy imports |
| `adaptive_ps4/models.py` | Re-exports `acmg.criteria.ps4.models`. | Legacy imports |
| `adaptive_ps4/statistics.py` | Re-exports `acmg.criteria.ps4.statistics`. | Legacy imports |
| `ps2_pm6/__init__.py` | Re-exports `acmg.criteria.ps2_pm6`. | Legacy imports |
| `ps2pm6/__init__.py` | Re-exports `acmg.criteria.ps2_pm6`. | Legacy imports |
| `ps2pm6/evaluator.py` | Re-exports `acmg.criteria.ps2_pm6.evaluator`. | Legacy imports |
| `ps2pm6/models.py` | Re-exports `acmg.criteria.ps2_pm6.models`. | Legacy imports |
| `ps2pm6/specifications.py` | Re-exports `acmg.criteria.ps2_pm6.specifications`. | Legacy imports |
| `pvs1/__init__.py` | Re-exports `acmg.criteria.pvs1`. | Legacy imports |
| `pvs1/evaluator.py` | Re-exports `acmg.criteria.pvs1.evaluator`. | Legacy imports |
| `pvs1/models.py` | Re-exports `acmg.criteria.pvs1.models`. | Legacy imports |
| `pvs1/specifications.py` | Re-exports `acmg.criteria.pvs1.specifications`. | Legacy imports |
### Tests and support modules
| File | Purpose |
| --- | --- |
| `tests/__init__.py` | Test package marker. |
| `tests/pm3_fixtures.py` | Shared synthetic fixture builders for PM3. |
| `tests/pp3_fixtures.py` | Shared synthetic fixture builders for PP3. |
| `tests/ps3_fixtures.py` | Shared synthetic fixture builders for PS3. |
| `tests/test_adaptive_pm6.py` | Synthetic checks for adaptive pm6 behavior and boundaries. |
| `tests/test_adaptive_ps1_engine.py` | Synthetic checks for adaptive ps1 engine behavior and boundaries. |
| `tests/test_assay_validation.py` | Synthetic checks for assay validation behavior and boundaries. |
| `tests/test_ba1.py` | Synthetic checks for ba1 behavior and boundaries. |
| `tests/test_bp1_eligibility.py` | Synthetic checks for bp1 eligibility behavior and boundaries. |
| `tests/test_bp1_exceptions.py` | Synthetic checks for bp1 exceptions behavior and boundaries. |
| `tests/test_bp1_mechanism.py` | Synthetic checks for bp1 mechanism behavior and boundaries. |
| `tests/test_bp1_missing_evidence.py` | Synthetic checks for bp1 missing evidence behavior and boundaries. |
| `tests/test_bp1_vcep.py` | Synthetic checks for bp1 vcep behavior and boundaries. |
| `tests/test_bp2.py` | Synthetic checks for bp2 behavior and boundaries. |
| `tests/test_bp3_functional_regions.py` | Synthetic checks for bp3 functional regions behavior and boundaries. |
| `tests/test_bp3_missing_evidence.py` | Synthetic checks for bp3 missing evidence behavior and boundaries. |
| `tests/test_bp3_repeat_indels.py` | Synthetic checks for bp3 repeat indels behavior and boundaries. |
| `tests/test_bp3_vcep.py` | Synthetic checks for bp3 vcep behavior and boundaries. |
| `tests/test_bp4_calibration.py` | Synthetic checks for bp4 calibration behavior and boundaries. |
| `tests/test_bp4_missense.py` | Synthetic checks for bp4 missense behavior and boundaries. |
| `tests/test_bp4_splicing.py` | Synthetic checks for bp4 splicing behavior and boundaries. |
| `tests/test_bp4_vcep.py` | Synthetic checks for bp4 vcep behavior and boundaries. |
| `tests/test_bp5_engine.py` | Synthetic checks for bp5 engine behavior and boundaries. |
| `tests/test_bp6.py` | Synthetic checks for bp6 behavior and boundaries. |
| `tests/test_bp7_engine.py` | Synthetic checks for bp7 engine behavior and boundaries. |
| `tests/test_bs1.py` | Synthetic checks for bs1 behavior and boundaries. |
| `tests/test_bs2.py` | Synthetic checks for bs2 behavior and boundaries. |
| `tests/test_bs3_calibration.py` | Synthetic checks for bs3 calibration behavior and boundaries. |
| `tests/test_bs3_functional_assays.py` | Synthetic checks for bs3 functional assays behavior and boundaries. |
| `tests/test_bs3_missing_evidence.py` | Synthetic checks for bs3 missing evidence behavior and boundaries. |
| `tests/test_bs3_splicing.py` | Synthetic checks for bs3 splicing behavior and boundaries. |
| `tests/test_bs3_vcep.py` | Synthetic checks for bs3 vcep behavior and boundaries. |
| `tests/test_bs4.py` | Synthetic checks for bs4 behavior and boundaries. |
| `tests/test_calibration.py` | Synthetic checks for calibration behavior and boundaries. |
| `tests/test_classification.py` | Synthetic checks for classification behavior and boundaries. |
| `tests/test_domains.py` | Synthetic checks for domains behavior and boundaries. |
| `tests/test_evidence_dependencies.py` | Synthetic checks for evidence dependencies behavior and boundaries. |
| `tests/test_homozygous.py` | Synthetic checks for homozygous behavior and boundaries. |
| `tests/test_hotspots.py` | Synthetic checks for hotspots behavior and boundaries. |
| `tests/test_inframe_indels.py` | Synthetic checks for inframe indels behavior and boundaries. |
| `tests/test_missense.py` | Synthetic checks for missense behavior and boundaries. |
| `tests/test_missing_evidence.py` | Synthetic checks for missing evidence behavior and boundaries. |
| `tests/test_multiple_cases.py` | Synthetic checks for multiple cases behavior and boundaries. |
| `tests/test_pm1_missing_evidence.py` | Synthetic checks for pm1 missing evidence behavior and boundaries. |
| `tests/test_pm1_vcep.py` | Synthetic checks for pm1 vcep behavior and boundaries. |
| `tests/test_pm2.py` | Synthetic checks for pm2 behavior and boundaries. |
| `tests/test_pm3_missing_evidence.py` | Synthetic checks for pm3 missing evidence behavior and boundaries. |
| `tests/test_pm4_missing_evidence.py` | Synthetic checks for pm4 missing evidence behavior and boundaries. |
| `tests/test_pm5_engine.py` | Synthetic checks for pm5 engine behavior and boundaries. |
| `tests/test_pp1.py` | Synthetic checks for pp1 behavior and boundaries. |
| `tests/test_pp2.py` | Synthetic checks for pp2 behavior and boundaries. |
| `tests/test_pp3_calibration.py` | Synthetic checks for pp3 calibration behavior and boundaries. |
| `tests/test_pp3_missense.py` | Synthetic checks for pp3 missense behavior and boundaries. |
| `tests/test_pp3_missing_evidence.py` | Synthetic checks for pp3 missing evidence behavior and boundaries. |
| `tests/test_pp3_splicing.py` | Synthetic checks for pp3 splicing behavior and boundaries. |
| `tests/test_pp3_vcep.py` | Synthetic checks for pp3 vcep behavior and boundaries. |
| `tests/test_pp4_engine.py` | Synthetic checks for pp4 engine behavior and boundaries. |
| `tests/test_pp5.py` | Synthetic checks for pp5 behavior and boundaries. |
| `tests/test_ps2pm6_engine.py` | Synthetic checks for ps2pm6 engine behavior and boundaries. |
| `tests/test_ps3.py` | Synthetic checks for ps3 behavior and boundaries. |
| `tests/test_ps4.py` | Synthetic checks for ps4 behavior and boundaries. |
| `tests/test_pvs1_engine.py` | Synthetic checks for pvs1 engine behavior and boundaries. |
| `tests/test_repeat_regions.py` | Synthetic checks for repeat regions behavior and boundaries. |
| `tests/test_splicing.py` | Synthetic checks for splicing behavior and boundaries. |
| `tests/test_stop_loss.py` | Synthetic checks for stop loss behavior and boundaries. |
| `tests/test_trans.py` | Synthetic checks for trans behavior and boundaries. |
| `tests/test_unified_classification.py` | Synthetic checks for unified classification behavior and boundaries. |
| `tests/test_unphased.py` | Synthetic checks for unphased behavior and boundaries. |
| `tests/test_vcep.py` | Synthetic checks for vcep behavior and boundaries. |
| `tests/test_vcep_pp3.py` | Synthetic checks for vcep pp3 behavior and boundaries. |
### Documentation and configuration
| File | Responsibility |
| --- | --- |
| `adaptive_ba1/README.md` | Preserved standalone-package usage, source, or policy notes. |
| `adaptive_bp2/README.md` | Preserved standalone-package usage, source, or policy notes. |
| `adaptive_bp3/README.md` | Preserved standalone-package usage, source, or policy notes. |
| `adaptive_bp5/README.md` | Preserved standalone-package usage, source, or policy notes. |
| `adaptive_bp6/README.md` | Preserved standalone-package usage, source, or policy notes. |
| `adaptive_bs1/README.md` | Preserved standalone-package usage, source, or policy notes. |
| `adaptive_bs4/README.md` | Preserved standalone-package usage, source, or policy notes. |
| `adaptive_pm3/SOURCES.md` | Preserved standalone-package usage, source, or policy notes. |
| `adaptive_pm6/README.md` | Preserved standalone-package usage, source, or policy notes. |
| `adaptive_pm6/SOURCES.md` | Preserved standalone-package usage, source, or policy notes. |
| `adaptive_pp4/README.md` | Preserved standalone-package usage, source, or policy notes. |
| `docs/bp1-guideline-matrix.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/BP4_SOURCE_RULE_MATRIX.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/bp7-guideline-matrix.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/bs2-guideline-matrix.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/bs3-guideline-matrix.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/bs4-guideline-matrix.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/coverage-matrix.md` | PVS1 branch and feature coverage. |
| `docs/decision-table.md` | Criterion or engine decision-path table. |
| `docs/feature-contract.md` | Cross-criterion feature names, local aliases and dependency summary. |
| `docs/package-architecture.md` | Consolidation layout and compatibility summary. |
| `docs/pm1-guideline-matrix.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/pm2-guideline-matrix.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/pm4-guideline-matrix.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/pm5-guideline-matrix.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/pp1-guideline-matrix.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/pp2-guideline-matrix.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/pp3-guideline-matrix.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/pp5-policy.md` | Inactive PP5 policy and sources. |
| `docs/PROJECT_ARCHITECTURE.md` | This current end-to-end architecture guide. |
| `docs/ps1-decision-table.md` | Criterion or engine decision-path table. |
| `docs/ps2-pm6-decision-table.md` | Criterion or engine decision-path table. |
| `docs/ps3-decision-table.md` | Criterion or engine decision-path table. |
| `docs/ps4-guideline-matrix.md` | Criterion-specific guideline-to-rule mapping and scope. |
| `docs/unified-classification-source-matrix.md` | Final-combination source-to-code mapping. |
| `docs/unified-scientific-validation.md` | Criterion coverage and remaining validation limits. |
| `README.md` | Main installation, usage and criterion overview. |
| `pyproject.toml` | Build backend, package discovery, Python requirement, pytest configuration and test extra. |
| `.gitignore` | Excludes generated caches, build output and local virtual environments. |
| `acmg_feature_to_rule_mapping.csv` | Reference feature-to-rule index at repository root; no runtime module reads it. |

**Generated, not scientific source:** `__pycache__/`, `.pytest_cache/`, `build/`, and `adaptive_acmg_evidence.egg-info/` may be present. They are excluded from the source tree above and by `.gitignore`; do not treat `build/lib` copies as canonical code.

**Inventory count:** 125 canonical Python files, 122 compatibility Python files, 70 test/support Python files, and 37 Markdown files (including this report).

## 3. Execution flow

```mermaid
flowchart LR
    A[Already extracted biological features] --> B[acmg.classify_variant]
    B --> C[Validate supplied codes and explicit identity/context]
    C --> D[criterion_registry: select evaluator]
    D --> E[acmg.criteria: independent criterion decisions]
    E --> F[Normalize result rows]
    F --> G[evidence_dependencies: count, suppress, or flag]
    G --> H[classification_rules: supported combinations]
    H --> I[Result dictionary with decisions and audit details]
```

`acmg/__init__.py` exports `classify_variant` from `acmg/classifier/__init__.py`. That function checks the input mapping, criterion codes, mutually exclusive de novo input, the requested classification method, and explicit variant/gene/disease identifiers. It iterates `criterion_registry.EVALUATORS` in insertion order and calls each supplied evaluator. Six typed evaluators are constructed from nested dictionaries at this boundary; the others parse their own native input. Each evaluator applies its **own** guideline scope, eligibility, exclusions, and strength rules. The classifier converts decisions to dictionaries, validates status and strength, and supplies common audit fields. `resolve_dependencies` marks overlapping rows `counted=False` or records conflicts. `classify_criteria` applies the supported combinations; unresolved evaluator, scope, or evidence conflicts keep the final result at `VUS`. The caller receives a dictionary, not a dataclass.

Missing criterion payloads are not evaluated and appear as no-evidence `UNDETERMINED` placeholders for registered single codes. The classifier does not require a giant feature record or run biological feature extraction.

## 4. Where biological features enter

The public input shape is `{"variant_id": optional_id, "criteria": {"CODE": criterion_specific_payload}}`. `classification_method` may be omitted (defaults to `ACMG_AMP_2015`); any other value is rejected. Optional top-level `gene_id` and `disease_id`, when supplied, must agree with explicit criterion payload IDs. Explicitly different variant, gene, or disease IDs produce a `VUS` with an unresolved problem row. The check compares supplied identifiers literally; it does not reconcile ontology synonyms or equivalently named transcripts.

This complete synthetic example uses a supported PVS1 frameshift/NMD branch and the generic PM2 confirmed-absence branch. Every field shown inside those two payloads supplies eligibility, source/quality, or a branch condition; no other criterion input is required.

```python
from acmg import classify_variant

pvs1 = {
    "gene_id": "GENE-TEST", "disease_id": "DISEASE-TEST",
    "transcript_id": "TX-TEST.1", "variant_consequence": "frameshift_variant",
    "lof_disease_mechanism": "established", "gene_disease_validity": "definitive",
    "biologically_relevant_transcript": True, "nmd_prediction": True,
}
pm2 = {
    "variant_id": "synthetic-1", "gene_id": "GENE-TEST",
    "disease_id": "DISEASE-TEST", "inheritance_mode": "AUTOSOMAL_DOMINANT",
    "genome_assembly": "GRCh38",
    "normalized_variant_representation": "NC_000001.11:g.100A>G",
    "population_records": [{
        "name": "gnomAD", "version": "4.1", "genome_assembly": "GRCh38",
        "normalized_variant_representation": "NC_000001.11:g.100A>G",
        "allele_count": 0, "allele_number": 100000,
        "identity_confirmed": True, "filter_pass": True,
        "site_callable": True, "technically_reliable": True,
    }],
}
result = classify_variant({
    "variant_id": "synthetic-1", "criteria": {"PVS1": pvs1, "PM2": pm2}
})
summary = {
    "classification": result["classification"],
    "criteria_activated": [
        (row["criterion"], row["strength"], row["counted"])
        for row in result["criteria_activated"]
    ],
    "undetermined_count": len(result["criteria_undetermined"]),
    "conflicts": result["evidence_conflicts"],
    "reason": result["reason"],
}
print(summary)
```

The verified printed summary is:

```python
{'classification': 'LIKELY_PATHOGENIC',
 'criteria_activated': [('PVS1', 'VERY_STRONG', True), ('PM2', 'SUPPORTING', True)],
 'undetermined_count': 22, 'conflicts': [],
 'reason': 'PVS1 plus PM2_Supporting: ClinGen SVI PM2 recommendation v1.0 (2020).'}
```

Use `{"criteria": {"PVS1": pvs1}}` to evaluate only PVS1; the other criteria may be omitted. PVS1 alone produces `VUS` despite an activated `VERY_STRONG` criterion because it does not meet a final combination. For an explicitly unknown essential fact, `{"criteria": {"PVS1": {**pvs1, "nmd_prediction": {"state": "UNKNOWN"}}}}` produces an `UNDETERMINED` PVS1 row with `missing_critical_features == ["nmd_prediction"]`. This is different from `False`. Adding an unrelated field to a typed criterion payload is ignored at the classifier boundary; it does not supply evidence. Criterion parsers may additionally accept their own `KNOWN`, `UNKNOWN`, `CONFLICTING`, and `NOT_APPLICABLE` wrappers with source/version/provenance.

For de novo evidence, supply one payload under `PS2_PM6` (or one of its accepted `PS2`/`PM6` aliases). Supplying `PS2_PM6` together with either alias is rejected; supplying both `PS2` and `PM6` with different assessments is also rejected. The shared evaluator assigns at most one code.

## 5. Adaptive criterion activation

The table names **one verified sufficient pathway per criterion**, not every supported pathway or a universal required schema. Its outcomes were executed against the current canonical functions using the synthetic fixtures linked in the last column. Fields outside a chosen path may be omitted. The source files and criterion-specific decision tables govern other pathways.

| Criterion and verified path | Sufficient supplied facts for that path; optional examples omitted | Known exclusion and verified outcome | Activated outcome; fixture |
| --- | --- | --- | --- |
| PVS1, predicted-NMD frameshift | Gene/disease, transcript, frameshift consequence, established LoF mechanism, definitive gene-disease validity, relevant transcript, predicted NMD. No CNV, RNA, or functional assay fields. | Irrelevant transcript → `NOT_MET/NONE`. | `ACTIVATED/VERY_STRONG`; `test_pvs1_engine.py::eligible_null` |
| PS1, same protein change | Candidate allele/coding/protein/transcript/protein accession and missense/no-splice facts, plus **one** independently classified pathogenic reference with matching gene/disease, protein and transcript context. Other references are optional. | Reference classification not independent → `NOT_MET/NONE`. | `ACTIVATED/STRONG`; `test_adaptive_ps1_engine.py::standard_features` |
| PS2/PM6, one de novo case | Proband variant present, both parental genotypes absent, adequate call, gene-disease validity, compatible inheritance, phenotype category, and parentage status. Other cases and unrelated clinical facts are optional. | Poor variant call → `NOT_MET/NONE`. | Confirmed highly specific → `PS2/STRONG`; assumed parentage with the same phenotype → `PM6/MODERATE`; `test_ps2pm6_engine.py::eligible_case` |
| PM2, confirmed absence | Variant/build/normalized representation, gene/disease/inheritance, and one matching, callable, technically reliable population record with version, zero allele count and adequate allele number. Prevalence is optional. | Nonzero allele count under this generic path → `NOT_MET/NONE`. | `ACTIVATED/SUPPORTING`; `test_pm2.py::base` |
| PP3, calibrated missense | Missense consequence, relevant transcript, selected REVEL predictor with score, calibration, version and matching transcript. Other predictors are optional. | Conflicting computational evidence → `UNDETERMINED/NONE`. | REVEL 0.644 → `ACTIVATED/SUPPORTING`; `test_pp3_missense.py::case` |
| BA1, general population rule | Typed variant identity with verified mapping and normalization; matched quality/callability-checked population record with frequency and count/number; verified exception and guideline-modification checks. No phenotype is required. | Verified exception → `NOT_APPLICABLE`; no stand-alone strength. | Frequency 0.06 with 600/10000 → `ACTIVATED/STAND_ALONE`; `test_ba1.py::base` |
| BS3, calibrated normal assay | Disease mechanism and one mechanism-compatible, validated assay with reliable normal interpretation, accepted calibration and odds path. Other assays are optional. | Damaging assay result → `NOT_MET/NONE`. | `ACTIVATED/STRONG`; `test_bs3_functional_assays.py::case` |
| BP7, calibrated synonymous/no-splice path | Verified variant/transcript/build context, confident synonymous consequence, relevant transcript, no core splice motif, and one matching calibrated benign splice prediction. Conservation is optional on this path. | Core splice motif present → `NOT_MET/NONE`. | `ACTIVATED/SUPPORTING`; `test_bp7_engine.py::base` |

## 6. Independent criterion APIs

Every row below names an importable public function. Input model names are scoped to the listed package; similarly named classes in different packages are **not** interchangeable. A `dict` in the input column means a mapping is accepted directly; typed-only evaluators can be called with their exported model or supplied via the unified classifier's typed boundary. `ACMG` means its implemented general rule; `SVI` and `VCEP` mean only the versioned/scoped rules present in that package, not every published specification.

| Criterion | Canonical package / public function | Main input | Purpose | Implemented scope | Output |
| --- | --- | --- | --- | --- | --- |
| PVS1 | `acmg.criteria.pvs1.evaluate_pvs1` | BiologicalFeatures or mapping | LoF/null-variant impact | ACMG; SVI LoF/RNA; pinned VCEP | PVS1Decision |
| PS1 | `acmg.criteria.ps1.evaluate_ps1` | BiologicalFeatures or mapping | Same protein or splice effect as pathogenic reference | ACMG; SVI splice; pinned VCEP | PS1Decision |
| PS2 | `acmg.criteria.ps2_pm6.evaluate_ps2_pm6 / evaluate_de_novo` | DeNovoEvidence or mapping | Confirmed de novo | ClinGen SVI v1.1; scoped VCEP | PS2PM6Decision |
| PS3 | `acmg.criteria.ps3.evaluate_ps3` | PS3Evidence or mapping | Damaging functional assay | ACMG; SVI; scoped VCEP | PS3Decision |
| PS4 | `acmg.criteria.ps4.evaluate_ps4` | PS4Input or dict | Affected-case enrichment | ACMG; selected VCEP | PS4Result |
| PM1 | `acmg.criteria.pm1.evaluate_pm1` | BiologicalFeatures or dict | Hotspot/critical domain | ACMG; selected VCEP | Decision |
| PM2 | `acmg.criteria.pm2.evaluate_pm2` | PM2Input or dict | Population absence/rarity | ACMG; SVI Supporting; selected VCEP | PM2Result |
| PM3 | `acmg.criteria.pm3.evaluate_pm3` | PM3Input or dict | Recessive trans/co-occurrence | ClinGen SVI; selected VCEP | dict |
| PM4 | `acmg.criteria.pm4.evaluate_pm4` | BiologicalFeatures or mapping | Qualifying protein-length change | ACMG; selected VCEP | PM4Result |
| PM5 | `acmg.criteria.pm5.evaluate_pm5` | BiologicalFeatures or mapping | Different pathogenic missense at same residue | ACMG; selected VCEP | PM5Result |
| PM6 | `acmg.criteria.ps2_pm6.evaluate_ps2_pm6 / evaluate_de_novo` | DeNovoEvidence or mapping | Assumed de novo, shared with PS2 | ClinGen SVI v1.1; scoped VCEP | PS2PM6Decision |
| PP1 | `acmg.criteria.pp1.evaluate_pp1` | PP1Input or dict | Disease segregation | ACMG; selected VCEP | PP1Result |
| PP2 | `acmg.criteria.pp2.evaluate_pp2` | BiologicalFeatures or mapping | Disease-specific missense mechanism | ACMG; selected VCEP | Decision |
| PP3 | `acmg.criteria.pp3.evaluate_pp3` | Features or mapping | Calibrated damaging prediction | ClinGen SVI; selected VCEP | Result |
| PP4 | `acmg.criteria.pp4.evaluate_pp4` | PP4Features (typed) | Specific phenotype/family history | ACMG; selected VCEP | PP4Result |
| PP5 | `acmg.criteria.pp5.evaluate_pp5` | Mapping or None | External pathogenic assertion, inactive | ClinGen discontinuation policy | PP5Result |
| BA1 | `acmg.criteria.ba1.evaluate_ba1` | BA1Evidence (typed) | Stand-alone benign frequency | ACMG; SVI exceptions; selected VCEP | BA1Result |
| BS1 | `acmg.criteria.bs1.evaluate_bs1` | BS1Features (typed) | Disease-specific high frequency | Selected versioned VCEP | BS1Result |
| BS2 | `acmg.criteria.bs2.evaluate_bs2` | BS2Input or dict | Healthy-person genotype | ACMG; selected VCEP | BS2Result |
| BS3 | `acmg.criteria.bs3.evaluate_bs3` | BS3Input or mapping | Normal functional assay | ACMG; SVI; selected VCEP | BS3Result |
| BS4 | `acmg.criteria.bs4.evaluate_bs4` | BS4Input or dict | Nonsegregation | ACMG; selected VCEP | Decision |
| BP1 | `acmg.criteria.bp1.evaluate_bp1` | BP1Features or mapping | Missense in truncating-mechanism context | ACMG; selected VCEP | BP1Result |
| BP2 | `acmg.criteria.bp2.evaluate_bp2` | BiologicalFeatures (typed) | Benign phase/co-occurrence | ACMG; selected VCEP | BP2Result |
| BP3 | `acmg.criteria.bp3.evaluate_bp3` | BiologicalFeatures (typed) | In-frame repeat-region indel | ACMG; selected VCEP | BP3Result |
| BP4 | `acmg.criteria.bp4.evaluate_bp4` | BiologicalFeatures (typed) | Calibrated benign prediction | ClinGen SVI; selected VCEP | BP4Result |
| BP5 | `acmg.criteria.bp5.evaluate_bp5` | BP5Features or mapping | Alternative molecular diagnosis | ACMG; selected VCEP | BP5Result |
| BP6 | `acmg.criteria.bp6.evaluate_bp6` | Mapping or None | External benign assertion, inactive | ClinGen discontinuation policy | BP6Result |
| BP7 | `acmg.criteria.bp7.evaluate_bp7` | BP7Features or dict | Synonymous/no splice effect | ACMG; SVI; selected VCEP | BP7Result |

The separate compatibility API `acmg.criteria.pm6.evaluate_pm6(PM6Input | Mapping) -> PM6Result` preserves its own HNF1A/SLC9A6 behavior and is not registered for unified classification. PP5/BP6 functions return audit results but their codes are deliberately absent from `EVALUATORS`.

For example, `from acmg.criteria.ps1 import evaluate_ps1`, `from acmg.criteria.pm3 import evaluate_pm3`, and `from acmg.criteria.bp4 import evaluate_bp4` import the real functions. The independent functions return their native result types. `evaluate_pm3` returns a dictionary; most other functions return a typed decision object, many with `.to_dict()`. Inspect each result model before assuming its field names. The unified classifier provides one consistent dictionary boundary.

### Selecting a pinned guideline

General rules are selected from each evaluator's input and bundled guidance. A scoped VCEP is selected only when that evaluator's gene/disease and, where required, transcript/version or explicit profile match. For example, the following **independent** BP7 call selects the implemented TP53 2.1.0 rule; all values are synthetic:

```python
from acmg.criteria.bp7 import BP7Features, EvidenceState, SplicePrediction, evaluate_bp7

tp53 = BP7Features(
    variant_id="synthetic", gene_id="TP53", disease_id="MONDO:0018875",
    genome_assembly="GRCh38", transcript_id_version="NM_TEST.1",
    variant_consequence="synonymous", consequence_confident=True,
    biologically_relevant_transcript=True, variant_mapping_valid=True,
    exon_core_splice_motif=False, guideline="TP53_2.1.0",
    prediction=SplicePrediction(
        score=0.05, transcript_id_version="NM_TEST.1",
        genome_assembly="GRCh38", variant_match=True,
        applicable=True, calibrated=True, quality=EvidenceState.KNOWN,
    ),
)
decision = evaluate_bp7(tp53)
print(decision.status, decision.strength, decision.guideline["version"])
# ACTIVATED SUPPORTING 2.1.0
```

This illustrates selection behavior, not a real TP53 variant assertion. An unsupported or mismatched VCEP declaration is not treated as a generic positive result. PVS1 and PS1 can also take caller-pinned specifications, with their own provenance/scope checks. At the classifier level, a nonempty top-level `vcep_specification` is reported as unsupported for final VCEP-wide combining and the classification stays `VUS`, even if individual criteria can use scoped rules.

## 7. How to read the classifier response

`classify_variant(features)` returns a plain `dict`; call `result["classification"]`, not `result.to_dict()`. Its top-level keys are `variant_id`, `classification`, `classification_method`, `criteria_activated`, `criteria_not_met`, `criteria_undetermined`, `criteria_not_applicable`, `evidence_conflicts`, `evidence_dependencies_resolved`, `guideline`, and `reason`. `guideline` identifies the final combination framework; each evaluated criterion row has its own guideline when available.

```python
label = result["classification"]
for row in result["criteria_activated"]:
    print(row["criterion"], row["strength"], row["decision_path"])
for row in result["criteria_undetermined"]:
    print(row["criterion"], row.get("missing_critical_features"), row["reason"])
for row in result["criteria_not_applicable"]:
    print(row["criterion"], row["reason"])
for conflict in result["evidence_conflicts"]:
    print(conflict["criteria"], conflict["reason"])
for row in result["criteria_activated"]:
    print(row["evidence_provenance"], row["reason"])
```

Evaluated rows carry `criterion`, `status`, `strength`, `decision_path`, `features_used`, `missing_critical_features`, `exclusions_evaluated`, `evidence_provenance`, `guideline`, `reason`, and `counted`; native evaluators can add more fields. Empty provenance means none was supplied or emitted, **not** that a source was independently verified. A suppressed activated row stays in `criteria_activated` with `counted=False` and `suppression_reason`. No-input placeholders and integration problem rows in `criteria_undetermined` have fewer fields; use `.get()` on those rows.

`ACTIVATED` means a criterion's supported pathway met its evidence conditions. `NOT_MET` means supplied facts refute that path. `UNDETERMINED` means critical evidence, scope, quality, or a conflict remains unresolved. `NOT_APPLICABLE` means the rule cannot be applied in the supplied context, including inactive PP5/BP6. `PATHOGENIC`, `LIKELY_PATHOGENIC`, `LIKELY_BENIGN`, and `BENIGN` are final evidence combinations; `VUS` covers insufficient, contradictory, or unsupported combining. An activated criterion alone need not produce a non-VUS final classification.

## 8. Shared facts and evidence dependencies

The repository has **criterion-specific** evidence sections, not one shared biological feature object. These are actual overlapping concepts and local field names:

| Concept | Consumers and current representation | Normalization / reuse |
| --- | --- | --- |
| Variant consequence | `variant_consequence`, `candidate_variant_consequence`, `consequence` in PVS1, PS1, PM4/PM5, PP2, BP1/BP3/BP7 and others | Each parser normalizes only its own accepted labels; the classifier checks explicit variant IDs, not semantic equality of consequence calls. |
| Population frequency | PM2 `population_records`, BA1 `populations`, BS1 `frequencies`; allele frequency versus filtering allele frequency are distinct | Separate quality and scope rules. BA1 supersedes BS1 when both activate; contradictory PM2 plus benign frequency is flagged. |
| SpliceAI / REVEL | PP3 `predictions` and selected predictor; BP4 `predictions`; BP7 `prediction`; PVS1/PS1 splice evidence | Local calibration and transcript checks. Opposing PP3/BP4 rows sharing `predictor_used` are suppressed; PVS1/PP3 reuse needs an explicit overlap flag. |
| Functional assays | PS3 `functional_assays`, BS3 `assays` and assay IDs | Independent parsers; shared `assays_used` IDs suppress opposing PS3/BS3 rows, separate opposing assays flag a conflict. |
| Family/de novo | PS2/PM6 `observations`, PP1 family data, BS4 observations | Shared de novo engine pools cases once; PP1/BS4 shared `families_used` IDs suppress opposing rows. Other cross-criterion case identity is not generally inferred. |
| Phase | PM3 case `phase` and method; BP2 observations with phase/reliability | Each criterion checks its own phase. No universal PM3/BP2 cross-row suppression exists. |
| Protein length / repeat | PM4 `repeat_region_status`, BP3 `repeat_region_status` plus functional context | If both activate for one variant, both are suppressed as contradictory interpretations. |
| Disease mechanism / phenotype | PVS1 LoF, PP2/BP1 mechanism, PS2/PM6 phenotype, PP4 phenotype and related fields | Criterion-specific meanings remain distinct. No general PP4/de novo phenotype identity resolver exists. |

`acmg/evidence_dependencies.py` runs after evaluation. It suppresses the **same** opposing PS3/BS3 assay, PP3/BP4 predictor label, or PP1/BS4 family ID; if both opposing rows activate without a shared identifier, it records an unresolved conflict and the final class becomes `VUS`. PM4/BP3 are both suppressed if active. A PVS1/PP3 splice overlap suppresses PP3 only when the PP3 payload explicitly says `pvs1_splice_overlap=True` or `rna_splice_evidence_used=True`. The de novo package assigns PS2 or PM6 once and records that resolution. BA1 supersedes BS1; PM2 opposed to activated BA1 or BS1 is flagged and the active population rows are suppressed. Suppressed rows are auditable through `counted=False`; the resolver does not silently assume independence from shared source publications.

Current boundaries matter: PM3/BP2 phase overlap and PP4/de novo phenotype reuse are checked only to the extent of their native evaluator flags, not by a general cross-criterion identity registry. Predictor labels and family IDs are coarse identifiers, and absent IDs cannot establish independence. The classifier does not infer that two gene/disease ontology identifiers are equivalent. Review unresolved overlap before clinical use.

## 9. Guideline selection and final combining

Criterion-specific `guidelines.py` or `specifications.py` files hold the implemented ACMG/AMP, SVI, and selected versioned VCEP rules. Evaluators choose by their own accepted gene, disease, transcript, inheritance, input version, or caller-pinned specification; there is no universal VCEP catalogue. Incomplete or conflicting declarations generally yield `UNDETERMINED` or `NOT_APPLICABLE` at criterion level. The exact pathway and version are returned in each available criterion `guideline` field. The standalone PM6 HNF1A/SLC9A6 behavior is preserved separately from the classifier's shared de novo path; no silent merge occurs.

The combination module counts activated, unsuppressed rows by reported strength under ACMG/AMP 2015 Table 5. The classifier also handles the implemented SVI `PVS1 + PM2_Supporting` likely-pathogenic combination. Credible opposing pathogenic and benign evidence, evidence-dependency conflicts, unsupported non-PVS1 `VERY_STRONG` generic combining, explicit cross-payload identity/context mismatches, and explicit top-level VCEP-wide requests stay `VUS`. No arbitrary point scoring or Bayesian conversion is added at this layer; de novo points remain inside the shared PS2/PM6 evaluator.

## 10. Installation, Windows execution, and compatibility

`pyproject.toml` declares Python `>=3.12`, setuptools as the build backend, **no required third-party runtime dependencies**, and `pytest>=7` only in the `test` extra. There is no CLI or service. From the repository root in PowerShell:

```powershell
$python312 = 'C:\Users\PC\AppData\Local\Programs\Python\Python312\python.exe'
& $python312 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
.\.venv\Scripts\python.exe -c "from acmg import classify_variant; print(classify_variant({})['classification'])"
.\.venv\Scripts\python.exe -m pytest -q
```

The import prints `VUS`; use the input example above for an evaluated criterion. The first line is the installed Python 3.12 path on the inspected machine; adjust that path on another machine. `python -m pip install .` also builds the one distribution. `pytest` discovers `tests/` from `pyproject.toml`. All evaluated criterion code runs offline from supplied evidence.

Former top-level `pvs1`, `ps2pm6`, `ps2_pm6`, and `adaptive_*` Python modules delegate to their matching `acmg.criteria` implementations, preserving legacy imports. They remain packaged by the setuptools discovery pattern. `adaptive_ps2pm6`, `adaptive_ps2_pm6`, and `ps2_pm6` are aliases of the shared de novo API. `adaptive_pm6` instead delegates to `acmg.criteria.pm6`, an intentionally distinct standalone PM6 result schema with scoped HNF1A and SLC9A6 pathways; it is **not** the classifier's PS2/PM6 scorer. Keep those two input/result contracts separate. Wrapper READMEs and source matrices are historical usage documentation, not duplicate scientific implementations.

## 11. Scientific and technical limits

- The engine accepts already-extracted biological facts. There is no ANNOVAR, VEP, VCF, gnomAD/ClinVar extraction, raw prediction execution, patient annotation, or external database downloader.
- The tests use synthetic fixtures. Passing them verifies implemented decisions and integration, not diagnostic accuracy or validation on patient cohorts.
- VCEP coverage is limited to bundled or accepted pinned/versioned scopes. A complete VCEP-wide combination and Bayesian classification are unsupported.
- Exact identity/context strings are checked where explicitly supplied; ontology equivalence, all transcript/allele reconciliation, and all evidence-independence relationships are not resolved automatically.
- Missing essential evidence is unresolved, not positive. Extra irrelevant information is not an activation requirement. Unsupported biological situations should not be treated as verified decisions.
- PP5 and BP6 remain inactive. The standalone PM6 interface remains intentionally distinct pending a separate scientific equivalence review.

For scientific source details, follow the existing [validation summary](unified-scientific-validation.md), [feature contract](feature-contract.md), [classification source matrix](unified-classification-source-matrix.md), and criterion-specific matrices listed in the inventory below.
