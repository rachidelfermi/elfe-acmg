# BS4 source-to-rule matrix

| Source | Scope | Implemented rule |
|---|---|---|
| [Richards et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/) | General BS4 | Lack of segregation in affected family members is Strong benign evidence. Evaluate phenocopies, clinical status, and biological relationships. No universal observation-count or LOD threshold is supplied. |
| [ClinGen variant curation SOP v1.0, 2019](https://www.clinicalgenome.org/site/assets/files/3677/clingen_variant-curation_sopv1.pdf) | General interpretation | BS4 concerns phenotype-positive/genotype-negative affected members. Unaffected carriers ordinarily concern BS2. VCEP specifications take precedence when applicable. |
| [ClinGen PIK3R1 GN160 v1.0](https://erepo.clinicalgenome.org/cspec/ui/svi/doc/GN160) | PIK3R1, MONDO:1060136 | Each affected noncarrier needs at least 6 PS4 phenotype points; one qualifying member gives Supporting, at least two give Strong. |
| [ClinGen CTLA4 specification v1.0](https://cspec.genome.network/cspec/ui/svi/doc/1567102241?version=1.0.0) | CTLA4, MONDO:0014493 | Each affected noncarrier needs at least 6 PS4 phenotype points; one qualifying member gives Supporting, more than one gives Strong. |
| [ClinGen ENIGMA BRCA2 GN097 v1.2](https://cspec.genome.network/cspec/ui/svi/doc/GN097?version=1.2) | BRCA2, MONDO:0700269 | Validated co-segregation likelihood ratio (LR) ≤0.48: Supporting; ≤0.23: Moderate; ≤0.05: Strong. This API caps the published Very Strong tier at Strong and discloses that limit. |

## Evidence policy

The qualitative path requires a confirmed affected relative, reliable variant absence, reliable phenotype, confirmed biological relationship, compatible inheritance, and explicit assessment excluding phenocopy and heterogeneous cause for that observation. The same relative is never counted twice. Known alternative diagnoses or sample, genotype, mosaicism, or pedigree errors exclude that observation. Unknown essential facts cause `UNDETERMINED` unless another complete, independent observation already meets BS4. An unaffected carrier alone does not establish BS4 under the general rule.

The quantitative path is restricted to the pinned BRCA2 specification. The caller supplies a validated LR and confirms that the disease model and evidence independence are appropriate. Raw pedigrees and LOD scores are not converted to LRs. An applicable VCEP outside the pinned set must be reviewed manually; this offline package cannot discover new specifications. No case fixtures in this package are externally curated.
