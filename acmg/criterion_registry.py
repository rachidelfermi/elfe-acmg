"""Criterion entry points into the canonical scientific packages."""

from importlib import import_module


EVALUATORS = {
    "PVS1": ("acmg.criteria.pvs1", "evaluate_pvs1"),
    "PS1": ("acmg.criteria.ps1", "evaluate_ps1"),
    "PS2_PM6": ("acmg.criteria.ps2_pm6", "evaluate_ps2_pm6"),
    "PS3": ("acmg.criteria.ps3", "evaluate_ps3"),
    "PS4": ("acmg.criteria.ps4", "evaluate_ps4"),
    "PM1": ("acmg.criteria.pm1", "evaluate_pm1"),
    "PM2": ("acmg.criteria.pm2", "evaluate_pm2"),
    "PM3": ("acmg.criteria.pm3", "evaluate_pm3"),
    "PM4": ("acmg.criteria.pm4", "evaluate_pm4"),
    "PM5": ("acmg.criteria.pm5", "evaluate_pm5"),
    "PP1": ("acmg.criteria.pp1", "evaluate_pp1"),
    "PP2": ("acmg.criteria.pp2", "evaluate_pp2"),
    "PP3": ("acmg.criteria.pp3", "evaluate_pp3"),
    "PP4": ("acmg.criteria.pp4", "evaluate_pp4"),
    "BA1": ("acmg.criteria.ba1", "evaluate_ba1"),
    "BS1": ("acmg.criteria.bs1", "evaluate_bs1"),
    "BS2": ("acmg.criteria.bs2", "evaluate_bs2"),
    "BS3": ("acmg.criteria.bs3", "evaluate_bs3"),
    "BS4": ("acmg.criteria.bs4", "evaluate_bs4"),
    "BP1": ("acmg.criteria.bp1", "evaluate_bp1"),
    "BP2": ("acmg.criteria.bp2", "evaluate_bp2"),
    "BP3": ("acmg.criteria.bp3", "evaluate_bp3"),
    "BP4": ("acmg.criteria.bp4", "evaluate_bp4"),
    "BP5": ("acmg.criteria.bp5", "evaluate_bp5"),
    "BP7": ("acmg.criteria.bp7", "evaluate_bp7"),
}

DEPRECATED = frozenset({"PP5", "BP6"})
ALL_CRITERIA = frozenset(EVALUATORS) | DEPRECATED | {"PS2", "PM6"}


def evaluator_for(code: str):
    module, function = EVALUATORS[code]
    return getattr(import_module(module), function)
