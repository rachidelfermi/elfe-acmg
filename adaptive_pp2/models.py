"""Compatibility import path for :mod:`acmg.criteria.pp2.models`."""

from acmg.criteria.pp2.models import *  # noqa: F401,F403
from acmg.criteria.pp2.models import __dict__ as _canonical_namespace

def __getattr__(name):
    try:
        return _canonical_namespace[name]
    except KeyError:
        raise AttributeError(name) from None
