"""Pinned, verified BP1 specifications; no caller-defined strength rules."""

from dataclasses import dataclass


ACMG_2015 = {
    "name": "ACMG/AMP",
    "version": "2015",
    "rule_id": "BP1",
}


@dataclass(frozen=True)
class Specification:
    spec_id: str
    gene_id: str
    disease_id: str
    version: str
    mode: str
    url: str
    strength: str = "NONE"

    @property
    def guideline(self) -> dict[str, str]:
        return {"name": "ClinGen CSpec", "version": self.version, "rule_id": f"{self.spec_id}:BP1"}


# Each scope and BP1 disposition is stated in the linked CSpec document.
SPECIFICATIONS = (
    Specification("GN024", "DICER1", "MONDO:0100216", "1.4", "DISABLED",
                  "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN024"),
    Specification("GN012", "RYR1", "MONDO:0018493", "2.0", "DISABLED",
                  "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN012"),
    Specification("GN045", "MAP2K1", "MONDO:0021060", "2.3", "RASOPATHY_LOF",
                  "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN045", "SUPPORTING"),
)


def find_specification(gene_id: str | None, disease_id: str | None) -> Specification | None:
    return next((s for s in SPECIFICATIONS if s.gene_id == gene_id and s.disease_id == disease_id), None)
