"""Only the scoped rules below are encoded; registry rules may be revised."""
from .models import Guideline

ACMG = Guideline("ACMG/AMP", "2015", "PP4", "https://pubmed.ncbi.nlm.nih.gov/25741868/")
ADA = Guideline("ClinGen SCID VCEP ADA", "2.2", "GN114:PP4", "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN114")
ADA_GENE = {"ADA", "HGNC:186"}
ADA_DISEASE = {"MONDO:0007064"}
GT = Guideline("ClinGen Platelet Disorders VCEP", "2.1", "GN011:PP4", "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN011")
GT_GENES = {"ITGA2B", "HGNC:6138", "ITGB3", "HGNC:6156"}
GT_DISEASE = {"MONDO:0100326"}
CDH1 = Guideline("ClinGen CDH1 VCEP", "3.1", "GN007:PP4", "https://www.clinicalgenome.org/site/assets/files/7580/clingen_cdh1_acmg_specifications_v3_1.pdf")
CDH1_GENE = {"CDH1", "HGNC:999"}
CDH1_DISEASE = {"MONDO:0007648"}

