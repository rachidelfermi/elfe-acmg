"""Table 5 combinations and end-to-end evaluator integration."""

import pytest

from acmg.classification_rules import classify_criteria
from acmg.classifier import classify_variant


def rows(*items):
    return [{"criterion": code, "strength": strength} for code, strength in items]


@pytest.mark.parametrize("items,expected", [
    ([('PVS1','VERY_STRONG'),('PS1','STRONG')], 'PATHOGENIC'),
    ([('PVS1','VERY_STRONG'),('PM1','MODERATE'),('PM3','MODERATE')], 'PATHOGENIC'),
    ([('PVS1','VERY_STRONG'),('PM1','MODERATE'),('PP1','SUPPORTING')], 'PATHOGENIC'),
    ([('PVS1','VERY_STRONG'),('PP1','SUPPORTING'),('PP2','SUPPORTING')], 'PATHOGENIC'),
    ([('PS1','STRONG'),('PS4','STRONG')], 'PATHOGENIC'),
    ([('PS1','STRONG'),('PM1','MODERATE'),('PM3','MODERATE'),('PM4','MODERATE')], 'PATHOGENIC'),
    ([('PS1','STRONG'),('PM1','MODERATE'),('PM3','MODERATE'),('PP1','SUPPORTING'),('PP2','SUPPORTING')], 'PATHOGENIC'),
    ([('PS1','STRONG'),('PM1','MODERATE'),('PP1','SUPPORTING'),('PP2','SUPPORTING'),('PP3','SUPPORTING'),('PP4','SUPPORTING')], 'PATHOGENIC'),
    ([('PVS1','VERY_STRONG'),('PM2','MODERATE')], 'LIKELY_PATHOGENIC'),
    ([('PS1','STRONG'),('PM2','MODERATE')], 'LIKELY_PATHOGENIC'),
    ([('PS1','STRONG'),('PP1','SUPPORTING'),('PP2','SUPPORTING')], 'LIKELY_PATHOGENIC'),
    ([('PM1','MODERATE'),('PM3','MODERATE'),('PM4','MODERATE')], 'LIKELY_PATHOGENIC'),
    ([('PM1','MODERATE'),('PM3','MODERATE'),('PP1','SUPPORTING'),('PP2','SUPPORTING')], 'LIKELY_PATHOGENIC'),
    ([('PM1','MODERATE'),('PP1','SUPPORTING'),('PP2','SUPPORTING'),('PP3','SUPPORTING'),('PP4','SUPPORTING')], 'LIKELY_PATHOGENIC'),
    ([('BA1','STAND_ALONE')], 'BENIGN'),
    ([('BS1','STRONG'),('BS2','STRONG')], 'BENIGN'),
    ([('BS1','STRONG'),('BP1','SUPPORTING')], 'LIKELY_BENIGN'),
    ([('BP1','SUPPORTING'),('BP2','SUPPORTING')], 'LIKELY_BENIGN'),
    ([('PVS1','VERY_STRONG')], 'VUS'),
    ([('PVS1','VERY_STRONG'),('PM2','SUPPORTING')], 'LIKELY_PATHOGENIC'),
    ([('BA1','STAND_ALONE'),('PVS1','VERY_STRONG')], 'VUS'),
])
def test_table_5(items, expected):
    assert classify_criteria(rows(*items))[0] == expected


PVS1 = {
    'gene_id':'GENE-TEST', 'disease_id':'DISEASE-TEST', 'transcript_id':'TX-TEST.1',
    'variant_consequence':'frameshift_variant', 'lof_disease_mechanism':'established',
    'gene_disease_validity':'definitive', 'biologically_relevant_transcript':True,
    'nmd_prediction':True,
}

DE_NOVO = {
    'proband_variant_present':True, 'maternal_genotype':'ABSENT',
    'paternal_genotype':'ABSENT', 'variant_call_quality':'ADEQUATE',
    'maternity_confirmed':True, 'paternity_confirmed':True,
    'gene_disease_validity':'DEFINITIVE', 'inheritance_pattern':'AUTOSOMAL_DOMINANT',
    'phenotype_category':'HIGHLY_SPECIFIC',
}

BA1 = {
    'variant': {'variant_id':'v1', 'assembly':'GRCh38','chromosome':'1','position':100,
                'reference':'A','alternate':'G','normalized':True,'mapping_confidence':'KNOWN'},
    'populations': [{'database':'gnomAD','version':'4.1','assembly':'GRCh38',
                     'chromosome':'1','position':100,'reference':'A','alternate':'G',
                     'normalized':True,'ancestry':'AFR','continental':True,
                     'allele_frequency':0.06,'allele_count':600,'allele_number':10000,
                     'filter_pass':True,'quality':'KNOWN','mapping_ambiguous':False,
                     'callability_adequate':True}],
    'exception_check': {'state':'KNOWN','is_exception':False,'source':'ClinGen SVI list','version':'2018'},
    'guideline_modification_check':'KNOWN','known_modification':False,
}


def test_real_pvs1_and_de_novo_evaluators():
    result = classify_variant({'variant_id':'v1', 'criteria':{'PVS1':PVS1,'PS2_PM6':DE_NOVO}})
    assert result['classification'] == 'PATHOGENIC'
    assert {x['criterion'] for x in result['criteria_activated']} == {'PVS1','PS2'}
    assert result['evidence_dependencies_resolved']


def test_real_ba1_dataclass_adapter():
    result = classify_variant({'variant_id':'v1','criteria':{'BA1':BA1}})
    assert result['classification'] == 'BENIGN'
    assert result['criteria_activated'][0]['strength'] == 'STAND_ALONE'


def test_missing_critical_and_optional_evidence():
    incomplete = dict(PVS1)
    incomplete.pop('nmd_prediction')
    result = classify_variant({'criteria':{'PVS1':incomplete}})
    assert result['classification'] == 'VUS'
    assert result['criteria_undetermined'][0]['criterion'] == 'PVS1'


def test_deprecated_and_unsupported_vcep():
    result = classify_variant({'criteria':{'PP5':{},'BP6':{},'BA1':BA1},
                               'vcep_specification':{'name':'unverified'}})
    assert result['classification'] == 'VUS'
    assert {x['criterion'] for x in result['criteria_not_applicable']} == {'PP5','BP6'}


def test_one_de_novo_assessment_required():
    with pytest.raises(ValueError):
        classify_variant({'criteria':{'PS2_PM6':DE_NOVO,'PM6':DE_NOVO}})
