"""Offline, conservative evaluation of ACMG/AMP criterion PP5."""

from .evaluator import evaluate_pp5
from .models import PP5Result

__all__ = ["evaluate_pp5", "PP5Result"]
