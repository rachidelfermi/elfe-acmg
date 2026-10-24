"""Compatibility import path for :mod:`acmg.criteria.pvs1.specifications`."""

from acmg.criteria.pvs1.specifications import *  # noqa: F401,F403
from acmg.criteria.pvs1.specifications import __dict__ as _canonical_namespace

def __getattr__(name):
    try:
        return _canonical_namespace[name]
    except KeyError:
        raise AttributeError(name) from None
