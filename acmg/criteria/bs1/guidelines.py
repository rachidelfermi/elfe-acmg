"""Small, version-pinned rule registry transcribed from ClinGen CSpec."""

from .models import GuidelineRule


BUILTIN_RULES = (
    GuidelineRule(
        name="ClinGen Cardiomyopathy VCEP TPM1", version="1.0", rule_id="GN100:BS1",
        source_url="https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN100",
        gene="TPM1", disease="MONDO:0005045", inheritance="AD",
        threshold=0.0001, metric="POPMAX_FAF", operator=">=", database="gnomAD",
        verified=True,
    ),
    GuidelineRule(
        name="ClinGen SCID VCEP RMRP", version="1.3", rule_id="1529230095:BS1",
        source_url="https://erepo.clinicalgenome.org/cspec/ui/svi/doc/1529230095",
        gene="RMRP", disease="MONDO:0009595", inheritance="AR",
        threshold=0.00089, metric="POPMAX_FAF", operator=">", database="gnomAD",
        minimum_exome_depth=20, verified=True,
    ),
)


def matching_rule(gene: str, disease: str, inheritance: str) -> GuidelineRule | None:
    return next((r for r in BUILTIN_RULES if (r.gene, r.disease, r.inheritance) ==
                 (gene, disease, inheritance)), None)
