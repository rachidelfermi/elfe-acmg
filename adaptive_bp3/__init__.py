"""Compatibility import path for :mod:`acmg.criteria.bp3`."""

from acmg.criteria.bp3 import *  # noqa: F401,F403
from acmg.criteria.bp3 import __dict__ as _canonical_namespace

def __getattr__(name):
    try:
        return _canonical_namespace[name]
    except KeyError:
        raise AttributeError(name) from None
