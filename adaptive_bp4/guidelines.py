"""Compatibility import path for :mod:`acmg.criteria.bp4.guidelines`."""

from acmg.criteria.bp4.guidelines import *  # noqa: F401,F403
from acmg.criteria.bp4.guidelines import __dict__ as _canonical_namespace

def __getattr__(name):
    try:
        return _canonical_namespace[name]
    except KeyError:
        raise AttributeError(name) from None
