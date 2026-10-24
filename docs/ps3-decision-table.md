# PS3 guideline-to-code decision table

This table describes the generic rules implemented by `adaptive_ps3`. VCEP
rules take precedence only when a supplied rule is versioned, has a source,
matches gene/disease scope, and matches every declared condition.

| Question | Minimum sufficient evidence | Outcome |
|---|---|---|
| Is the assay applicable? | Known gene, disease, disease mechanism; variant-level assay; known relevance/compatibility to that mechanism | Continue. Known irrelevance is `NOT_APPLICABLE`; missing or conflicting critical facts are `UNDETERMINED`. |
| Is the specific assay usable? | Normal/wild-type and abnormal/null controls plus replicates; or exceptionally, documented dynamic range, separation, and calibration; or documented broad historical acceptance/validated kit | Continue. Known absent controls, failed reproducibility, artifacts, batch effects, or saturation prevent PS3. |
| Is formal pathogenic OddsPath available? | Rigorous OddsPath calibration and an interpretable variant result | `>2.1` Supporting; `>4.3` Moderate; `>18.7` Strong; `>350` Very Strong. Boundaries are strict. |
| No formal statistics: are clinical validation controls available? | Mixed LP/P and LB/B controls classified independently of the assay/PS3/BS3, with acceptable classification quality | At most Supporting for 1–10 total controls; Moderate for at least 11. No generic Strong fallback exists. |
| Is the variant result pathogenic functional evidence? | Result in the assay's validated abnormal range and compatible with the disease mechanism | Activate PS3 at the assay-validation strength. Effect magnitude alone is insufficient. |
| Is the result normal? | Result in a validated normal range | PS3 `NOT_MET`; flag for separate BS3 consideration without activating BS3. |
| Is the result intermediate? | Intermediate/hypomorphic/indeterminate range | `UNDETERMINED`; generic code assigns no strength. |
| Do multiple assays agree? | Two or more applicable results with the same direction | Use the best-validated, most mechanism-relevant assay. Never add strengths. |
| Do multiple assays conflict? | Abnormal and normal credible results | A more relevant and at least as validated assay may override; otherwise `UNDETERMINED`. |

Primary sources:

- Richards et al. 2015, ACMG/AMP standards, PMID 25741868:
  <https://pubmed.ncbi.nlm.nih.gov/25741868/>
- Brnich et al. 2019, ClinGen SVI PS3/BS3 recommendations, PMID 31892348:
  <https://pubmed.ncbi.nlm.nih.gov/31892348/>
- ClinGen Criteria Specification Registry (current VCEP specifications):
  <https://cspec.genome.network/cspec/ui/svi/summary>

The tests use synthetic fixtures so their provenance is not confused with
published expert-curated variants. Production users must pin the applicable
VCEP source/version and independently validate this software in their quality
system.
