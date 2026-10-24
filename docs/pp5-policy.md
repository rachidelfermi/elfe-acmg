# PP5 policy and provenance

- **Original rule:** [Richards et al., ACMG/AMP 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/), PP5: a reputable source reports a variant pathogenic while its supporting evidence is unavailable for independent assessment. Its evidence strength was Supporting. This is retained for historical documentation only. “2015” is its publication year, not a numbered version; the API reports `version: null`.
- **Default rule:** [Biesecker and Harrison, ClinGen SVI, 2018](https://pubmed.ncbi.nlm.nih.gov/29543229/) recommends discontinuing PP5 and BP6. The publication is dated 2018; it does not assign a numbered guideline version, so the API reports `version: null`. The [ClinGen guidance index](https://clinicalgenome.org/tools/clingen-variant-classification-guidance/) still lists that recommendation (index last updated July 2025). A classification is context, not independent variant evidence; reusing the evidence underlying it can double count.
- **Verified gene-specific example:** [ClinGen RASopathy VCEP NRAS, CSpec GN039 version 2.2](https://cspec.genome.network/cspec/ui/svi/doc/GN039?version=2.2), scoped to HGNC:7989 and MONDO:0021060, explicitly marks PP5 Not Applicable. The evaluator applies this specification only when its ID, version, gene and disease all match. Other or unknown VCEP specifications are reported as unverified and do not change the current ClinGen decision. The evaluator never infers a VCEP exception from a user-supplied `pp5` flag.

`pp5_usage_policy` accepts `CURRENT_CLINGEN` (default) or `HISTORICAL_ACMG_2015`.
Historical selection returns `UNDETERMINED` and no strength unless a verified
applicable VCEP excludes PP5, in which case it returns `NOT_APPLICABLE`.
Source reputation, identity matching, assertion currency, and evidence
independence cannot be established from a classification or review status alone.
The historical path therefore never awards PP5 automatically.

The input can contain the flat `external_classification` and related
`classification_*` fields, or an `external_assertions` list of those mappings
or `ExternalAssertion` records. Missing values remain `None`; classifications
can explicitly be `UNKNOWN` or `CONFLICTING`. The current policy needs no
variant metadata or external assertion. `features_used` lists only inputs that
control the PP5 policy decision; assertions and supplied variant identity are
retained separately for audit. An explicit `guideline_version` must match the
selected verified policy. Neither the ACMG/AMP 2015 publication nor the ClinGen
2018 recommendation has a numbered version in this implementation.
