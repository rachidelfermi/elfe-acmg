"""Compatibility import path for :mod:`acmg.criteria.pm5.evaluator`."""

from acmg.criteria.pm5.evaluator import *  # noqa: F401,F403
from acmg.criteria.pm5.evaluator import __dict__ as _canonical_namespace

def __getattr__(name):
    try:
        return _canonical_namespace[name]
    except KeyError:
        raise AttributeError(name) from None
