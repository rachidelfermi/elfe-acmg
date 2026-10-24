"""Versioned local source-to-rule mapping for supported PS1 pathways.

The evaluator intentionally has no network dependency.  These immutable source
records are the verified generic rules implemented by this release.  VCEP
rules are supplied as pinned configuration objects by the caller and are never
silently guessed from a gene name.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from .output import GuidelineReference


ACMG_2015_URL = "https://www.nature.com/articles/gim201530.pdf"
WALKER_2023_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC10357475/"
CSPEC_URL = "https://cspec.genome.network/cspec/ui/svi/summary"


ACMG_PS1 = GuidelineReference(
    authority="ACMG/AMP",
    name="Standards and guidelines for the interpretation of sequence variants",
    version="2015 (PMID: 25741868; DOI: 10.1038/gim.2015.30)",
    source_url=ACMG_2015_URL,
    source_locator="Table 3 - PS1; PS1 Same amino acid change",
    rule_id="PS1",
)

CLINGEN_SPLICE_PS1 = GuidelineReference(
    authority="ClinGen Sequence Variant Interpretation Splicing Subgroup",
    name="Using the ACMG/AMP framework to capture evidence related to predicted and observed impact on splicing",
    version="2023 (PMID: 37352859; DOI: 10.1016/j.ajhg.2023.06.002)",
    source_url=WALKER_2023_URL,
    source_locator="Table 2 - PS1 code weights for same predicted splicing event",
    rule_id="PS1",
)


@dataclass(frozen=True)
class RuleTableRow:
    """Human-readable rule coverage metadata exposed to audit tools."""

    path: str
    source: GuidelineReference
    eligibility: str
    exclusions: str
    strengths: str


GENERIC_RULE_TABLE: tuple[RuleTableRow, ...] = (
    RuleTableRow(
        path="SAME_AMINO_ACID_CHANGE",
        source=ACMG_PS1,
        eligibility=(
            "Same validated amino-acid substitution, distinct underlying allele, "
            "independently established Pathogenic reference, compatible disease "
            "context, and no predicted or observed splice-mediated mechanism."
        ),
        exclusions=(
            "Same allele, different protein effect, invalid/unresolved protein "
            "reference equivalence, non-Pathogenic/circular/conflicting reference, "
            "or a splice-mediated/incompatible mechanism."
        ),
        strengths="Strong only.",
    ),
    RuleTableRow(
        path="SPLICING_SAME_PREDICTED_EVENT",
        source=CLINGEN_SPLICE_PS1,
        eligibility=(
            "Precisely matched predicted event, candidate prediction similar to or "
            "stronger than the reference prediction, same donor/acceptor motif or "
            "same nucleotide as applicable, and an independently classified P/LP reference."
        ),
        exclusions=(
            "Different allele requirement unmet, event or prediction-strength "
            "comparison absent/conflicting, positions outside the Table 2 cases, "
            "or a circular/conflicting reference classification."
        ),
        strengths="Strong, Moderate, or Supporting only as listed in Walker et al. Table 2.",
    ),
)


def rule_table_as_dicts() -> tuple[Mapping[str, str], ...]:
    """Return local rule coverage metadata without exposing mutable globals."""

    return tuple(
        {
            "path": row.path,
            "source": row.source.name,
            "version": row.source.version,
            "source_url": row.source.source_url,
            "source_locator": row.source.source_locator,
            "eligibility": row.eligibility,
            "exclusions": row.exclusions,
            "strengths": row.strengths,
        }
        for row in GENERIC_RULE_TABLE
    )

