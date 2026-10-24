"""Compatibility import path for :mod:`acmg.criteria.ps1.output`."""

from acmg.criteria.ps1.output import *  # noqa: F401,F403
from acmg.criteria.ps1.output import __dict__ as _canonical_namespace

def __getattr__(name):
    try:
        return _canonical_namespace[name]
    except KeyError:
        raise AttributeError(name) from None
