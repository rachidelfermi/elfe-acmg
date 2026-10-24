"""Synthetic cases for PM3 unit tests; none are published validation cases."""

from adaptive_pm3 import Case, PM3Input, SecondVariant


def second(**changes):
    values = dict(variant_id="v2", gene_id="PAH", genotype="HET",
                  classification="PATHOGENIC", classification_independent=True,
                  classification_credible=True,
                  disease_match=True, rare_for_disease=True, call_valid=True)
    values.update(changes)
    return SecondVariant(**values)


def case(**changes):
    values = dict(case_id="case1", affected=True, clinical_disease_match=True,
                  candidate_genotype="HET", candidate_call_valid=True,
                  second_variant=second(), phase="TRANS", phase_method="molecular")
    values.update(changes)
    return Case(**values)


def data(*cases, **changes):
    values = dict(candidate_variant_id="v1", gene_id="PAH", disease_id="PKU",
                  inheritance_mode="AUTOSOMAL_RECESSIVE", candidate_rare_for_disease=True,
                  cases=list(cases) or [case()])
    values.update(changes)
    return PM3Input(**values)
