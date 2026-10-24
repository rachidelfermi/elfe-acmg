"""Integration of the existing ACMG/ClinGen criterion evaluators."""

from .classifier import classify_variant

__all__ = ["classify_variant"]
