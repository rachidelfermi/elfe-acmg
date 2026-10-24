# PM2 rule matrix

This package evaluates **only PM2** from supplied, already-extracted evidence. The rule set was checked against the sources below on 2026-10-23. All test variants are synthetic rule-boundary examples, not independently curated clinical classifications.

| Scope | Minimum sufficient evidence | Positive result | Known failure or unresolved evidence | Source |
| --- | --- | --- | --- | --- |
| General ACMG/AMP 2015 as revised by ClinGen SVI 2020 | Exact matched variant and assembly, reliable callable source, AC=0, AN>0; no established BA1/BS1 or unresolved BS2 | PM2 Supporting | Nonzero count does not satisfy the general absence branch. For recessive disease without a pinned rarity cutoff, return `UNDETERMINED`. Conflicting calls or unavailable counts remain unresolved. | [ACMG/AMP 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/); [ClinGen SVI PM2 v1.0](https://www.clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf) |
| ACADVL / VLCAD deficiency, MONDO:0008723, autosomal recessive, VCEP v2.2 | Reliable gnomAD record; complete continental ancestry assessment; at least one eligible ancestry group with AN>2000; maximum eligible ancestry MAF <0.001 | PM2 Supporting | MAF ≥0.001 in an eligible ancestry group is `NOT_MET`; unknown ancestry completeness or frequency is `UNDETERMINED`. Observed homozygotes require expert discussion. | [ClinGen ACADVL VCEP v2.2](https://erepo.clinicalgenome.org/cspec/ui/svi/doc/639509375) |
| F9 / hemophilia B, MONDO:0010604, X-linked, VCEP v2.1 | Reliable matched callable population record, zero male hemizygotes, positive number of male alleles | PM2 Supporting | Any observed male hemizygote is `NOT_MET`; missing male count or male allele number is `UNDETERMINED`. Female/global absence alone is insufficient. | [ClinGen Coagulation Factor Deficiency VCEP v2.1](https://erepo.clinicalgenome.org/cspec/ui/svi/doc/1527859052) |

## Evidence handling

- A population record must have matching assembly and normalized allele, confirmed identity, passing filter, callable site, and explicit technical reliability. Explicitly adverse mapping, detection, sample, or sampling flags prevent that record from supporting PM2.
- An absent database row is not an allele count of zero. A zero count requires a positive allele number. Missing ancestry frequencies never imply rarity.
- BA1 or BS1 already established for the same variant/disease prevents PM2 activation. BS2 and explicit frequency conflicts require reconciliation. Multiple sources are not multiplied into extra pathogenic evidence.
- General SVI gives a strength adjustment, not a universal numerical rarity threshold. The evaluator therefore uses no invented general frequency cutoff, and does not derive one from prevalence or penetrance.
- The two VCEP rules above are pinned. Other VCEP specifications and their thresholds are **not implemented**; users should check the applicable current specification before relying on the general pathway for such a gene/disease pair.
- Structural, multiallelic, indel, and mitochondrial input is held unresolved until the relevant normalization/detection context is explicitly confirmed. This package does not normalize variants, estimate coverage, query databases, assess genotype clinical meaning, or classify pathogenicity.
