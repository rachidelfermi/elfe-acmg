"""Small, pinned source-to-rule table; caller-supplied specifications are not trusted."""

from dataclasses import dataclass


ACMG_2015 = {"name": "ACMG/AMP", "version": "2015", "rule_id": "BP5"}


@dataclass(frozen=True)
class Specification:
    spec_id: str
    gene_ids: tuple[str, ...]
    disease_ids: tuple[str, ...]
    version: str
    mode: str
    url: str

    @property
    def guideline(self) -> dict[str, str]:
        return {"name": "ClinGen CSpec", "version": self.version,
                "rule_id": f"{self.spec_id}:BP5"}


# GN011: BP5 explicitly prohibited for Glanzmann thrombasthenia.
# GN092: BRCA1 clinical LR only; ordinary co-observation is not BP5 evidence.
# GN127: BP5/BP2 uses a supplemental case-point table not implemented here.
SPECIFICATIONS = (
    Specification("GN011", ("ITGA2B", "ITGB3", "HGNC:6138", "HGNC:6156"),
                  ("MONDO:0100326",), "2.1", "DISABLED",
                  "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN011"),
    Specification("GN092", ("BRCA1", "HGNC:1100"), ("MONDO:0700268",),
                  "1.2", "CLINICAL_LR",
                  "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN092"),
    Specification("GN127", ("RRAS2", "HGNC:17271"), ("MONDO:0021060",),
                  "1.3", "REQUIRES_CASE_SCORING",
                  "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN127"),
)


def find_specification(gene_id: str, disease_id: str) -> Specification | None:
    return next((s for s in SPECIFICATIONS if gene_id in s.gene_ids and
                 disease_id in s.disease_ids), None)
