"""Compatibility import path for :mod:`acmg.criteria.ps1.rules.missense`."""

from acmg.criteria.ps1.rules.missense import *  # noqa: F401,F403
from acmg.criteria.ps1.rules.missense import __dict__ as _canonical_namespace

def __getattr__(name):
    try:
        return _canonical_namespace[name]
    except KeyError:
        raise AttributeError(name) from None
