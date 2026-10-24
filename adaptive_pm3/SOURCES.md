# PM3 source to rule map

The implementation uses only these pinned specifications. Clinical curators must
check whether a newer, applicable VCEP specification supersedes them.

| Implemented rule | Primary source |
| --- | --- |
| Recessive, affected patient; pathogenic or likely pathogenic partner | [ACMG/AMP 2015 PM3](https://pubmed.ncbi.nlm.nih.gov/25741868/); [SVI v1.0, p. 1](https://www.clinicalgenome.org/site/assets/files/3717/svi_proposal_for_pm3_criterion_-_version_1.pdf) |
| Trans P/LP 1.0; unknown phase P 0.5, LP 0.25; homozygous 0.5; trans rare VUS 0.25 | [SVI v1.0, Table 1](https://www.clinicalgenome.org/site/assets/files/3717/svi_proposal_for_pm3_criterion_-_version_1.pdf) |
| Supporting 0.5; Moderate 1; Strong 2; Very Strong 4 | [SVI v1.0, Table 2](https://www.clinicalgenome.org/site/assets/files/3717/svi_proposal_for_pm3_criterion_-_version_1.pdf) |
| Both variants sufficiently rare, no circular partner classification; homozygous cap 1.0, VUS cap 0.5; single informative parent can establish trans | [SVI v1.0, considerations](https://www.clinicalgenome.org/site/assets/files/3717/svi_proposal_for_pm3_criterion_-_version_1.pdf) |
| GAA/Pompe diagnosis, pseudodeficiency context, repeated unphased pair cap of two cases, partner classified under VCEP rules; no VUS row in GAA table | [LSD VCEP GAA v2, pp. 10–11](https://www.clinicalgenome.org/site/assets/files/6993/clingen_lsd_acmg_specifications_v2-1.pdf) |

The VCEP does **not** change the general PM3 strength thresholds. Its repeated
unphased pair limit and VUS omission are the implemented GAA-specific changes.
No other VCEP is implemented; selecting an unknown specification fails closed.
