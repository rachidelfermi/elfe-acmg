"""Compatibility import path for :mod:`acmg.criteria.ps2_pm6.specifications`."""

from acmg.criteria.ps2_pm6.specifications import *  # noqa: F401,F403
from acmg.criteria.ps2_pm6.specifications import __dict__ as _canonical_namespace

def __getattr__(name):
    try:
        return _canonical_namespace[name]
    except KeyError:
        raise AttributeError(name) from None
