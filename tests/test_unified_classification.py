"""ACMG/AMP 2015 Table 5 boundaries and integration smoke tests."""

import pytest

from acmg import classify_variant
from acmg.classification_rules import classify_criteria


def rows(*codes):
    return [{"criterion": code, "strength": strength} for code, strength in codes]


@pytest.mark.parametrize("criteria,expected", [
    ([], "VUS"),
    ([('PVS1', 'VERY_STRONG'), ('PS1', 'STRONG')], 'PATHOGENIC'),
    ([('PVS1', 'VERY_STRONG'), ('PM1', 'MODERATE'), ('PM2', 'MODERATE')], 'PATHOGENIC'),
    ([('PVS1', 'VERY_STRONG'), ('PM1', 'MODERATE'), ('PP1', 'SUPPORTING')], 'PATHOGENIC'),
    ([('PVS1', 'VERY_STRONG'), ('PP1', 'SUPPORTING'), ('PP2', 'SUPPORTING')], 'PATHOGENIC'),
    ([('PS1', 'STRONG'), ('PS3', 'STRONG')], 'PATHOGENIC'),
    ([('PS1', 'STRONG'), ('PM1', 'MODERATE'), ('PM2', 'MODERATE'), ('PM3', 'MODERATE')], 'PATHOGENIC'),
    ([('PS1', 'STRONG'), ('PM1', 'MODERATE'), ('PM2', 'MODERATE'), ('PP1', 'SUPPORTING'), ('PP2', 'SUPPORTING')], 'PATHOGENIC'),
    ([('PS1', 'STRONG'), ('PM1', 'MODERATE'), ('PP1', 'SUPPORTING'), ('PP2', 'SUPPORTING'), ('PP3', 'SUPPORTING'), ('PP4', 'SUPPORTING')], 'PATHOGENIC'),
    ([('PVS1', 'VERY_STRONG'), ('PM1', 'MODERATE')], 'LIKELY_PATHOGENIC'),
    ([('PVS1', 'VERY_STRONG'), ('PM2', 'SUPPORTING')], 'LIKELY_PATHOGENIC'),
    ([('PS1', 'STRONG'), ('PM1', 'MODERATE')], 'LIKELY_PATHOGENIC'),
    ([('PS1', 'STRONG'), ('PP1', 'SUPPORTING'), ('PP2', 'SUPPORTING')], 'LIKELY_PATHOGENIC'),
    ([('PM1', 'MODERATE'), ('PM2', 'MODERATE'), ('PM3', 'MODERATE')], 'LIKELY_PATHOGENIC'),
    ([('PM1', 'MODERATE'), ('PM2', 'MODERATE'), ('PP1', 'SUPPORTING'), ('PP2', 'SUPPORTING')], 'LIKELY_PATHOGENIC'),
    ([('PM1', 'MODERATE'), ('PP1', 'SUPPORTING'), ('PP2', 'SUPPORTING'), ('PP3', 'SUPPORTING'), ('PP4', 'SUPPORTING')], 'LIKELY_PATHOGENIC'),
    ([('BA1', 'STAND_ALONE')], 'BENIGN'),
    ([('BS1', 'STRONG'), ('BS2', 'STRONG')], 'BENIGN'),
    ([('BS1', 'STRONG'), ('BP1', 'SUPPORTING')], 'LIKELY_BENIGN'),
    ([('BP1', 'SUPPORTING'), ('BP2', 'SUPPORTING')], 'LIKELY_BENIGN'),
    ([('BA1', 'STAND_ALONE'), ('PP1', 'SUPPORTING')], 'VUS'),
    ([('PS2', 'VERY_STRONG'), ('PM1', 'MODERATE')], 'VUS'),
])
def test_table_5(criteria, expected):
    assert classify_criteria(rows(*criteria))[0] == expected


def test_missing_evidence_is_not_negative():
    result = classify_variant({})
    assert result['classification'] == 'VUS'
    assert result['criteria_activated'] == []
    assert any(row['criterion'] == 'PVS1' for row in result['criteria_undetermined'])


def test_deprecated_codes_cannot_activate():
    result = classify_variant({'criteria': {
        'PP5': {'external_classification': 'PATHOGENIC'},
        'BP6': {'external_classification': 'BENIGN'}}})
    assert result['classification'] == 'VUS'
    assert all(row['status'] == 'NOT_APPLICABLE' for row in result['criteria_not_applicable'])


def test_real_pp3_evaluator_preserves_decision_provenance():
    evidence = {'PP3': {
        'variant_id': 'synthetic', 'variant_consequence': 'MISSENSE',
        'transcript_id_version': 'NM_TEST.1', 'biologically_relevant_transcript': True,
        'selected_missense_predictor': 'REVEL',
        'predictions': [{'name': 'REVEL', 'score': 0.80, 'version': '1',
                         'calibration': 'PEJAVER_2022', 'variant_id': 'synthetic',
                         'transcript_id_version': 'NM_TEST.1'}]}}
    result = classify_variant({'criteria': evidence})
    assert result['classification'] == 'VUS'
    assert result['criteria_activated'][0]['criterion'] == 'PP3'
    assert result['criteria_activated'][0]['guideline']
    assert result['criteria_activated'][0]['features_used']
    assert all(key in result['criteria_activated'][0] for key in (
        'decision_path', 'missing_critical_features', 'exclusions_evaluated',
        'evidence_provenance', 'reason'))


def test_unsupported_criterion_is_rejected():
    with pytest.raises(ValueError, match='Unknown criteria'):
        classify_variant({'criteria': {'PS9': {}}})


def test_explicit_variant_identity_mismatch_prevents_combination():
    result = classify_variant({'variant_id': 'v1', 'criteria': {
        'PVS1': {'variant_id': 'v2'}}})
    assert result['classification'] == 'VUS'
    assert any(row['criterion'] == 'variant_identity' for row in result['criteria_undetermined'])


@pytest.mark.parametrize('field', ['gene_id', 'disease_id'])
def test_explicit_disease_context_mismatch_prevents_combination(field):
    pvs1 = {
        'gene_id': 'GENE-TEST', 'disease_id': 'DISEASE-TEST',
        'transcript_id': 'TX-TEST.1', 'variant_consequence': 'frameshift_variant',
        'lof_disease_mechanism': 'established',
        'gene_disease_validity': 'definitive',
        'biologically_relevant_transcript': True, 'nmd_prediction': True,
    }
    pm2 = {
        'variant_id': 'synthetic-1', 'genome_assembly': 'GRCh38',
        'normalized_variant_representation': 'NC_000001.11:g.100A>G',
        'gene_id': 'GENE-TEST', 'disease_id': 'DISEASE-TEST',
        'inheritance_mode': 'AUTOSOMAL_DOMINANT',
        'population_records': [{
            'name': 'gnomAD', 'version': '4.1', 'genome_assembly': 'GRCh38',
            'normalized_variant_representation': 'NC_000001.11:g.100A>G',
            'allele_count': 0, 'allele_number': 100000,
            'identity_confirmed': True, 'filter_pass': True,
            'site_callable': True, 'technically_reliable': True,
        }],
    }
    consistent = classify_variant({'variant_id': 'synthetic-1',
                                   'criteria': {'PVS1': pvs1, 'PM2': pm2}})
    assert consistent['classification'] == 'LIKELY_PATHOGENIC'
    inconsistent = classify_variant({'variant_id': 'synthetic-1',
                                     'criteria': {'PVS1': pvs1,
                                                  'PM2': {**pm2, field: 'OTHER'}}})
    assert inconsistent['classification'] == 'VUS'
    assert any(row['criterion'] == field for row in inconsistent['criteria_undetermined'])


@pytest.mark.parametrize('state,value,expected', [
    ('KNOWN', True, 'ACTIVATED'),
    ('KNOWN', False, 'NOT_MET'),
    ('UNKNOWN', None, 'UNDETERMINED'),
    ('CONFLICTING', None, 'UNDETERMINED'),
])
def test_unified_bp3_respects_evidence_wrapper_state(state, value, expected):
    bp3 = {
        'gene_id': 'GENE1', 'disease_id': 'DISEASE1',
        'variant_consequence': 'INFRAME_DELETION',
        'reading_frame_preserved': True,
        'biologically_relevant_transcript': True,
        'protein_consequence_confidence': True,
        'repeat_region_status': {'state': state, 'value': value, 'source': 'curator'},
        'repeat_region_functional_importance': 'NO_KNOWN_FUNCTION',
        'repeat_length_disease_mechanism': False,
    }
    result = classify_variant({'criteria': {'BP3': bp3}})
    matches = [row for key in ('criteria_activated', 'criteria_not_met',
                              'criteria_undetermined') for row in result[key]
               if row['criterion'] == 'BP3']
    assert len(matches) == 1
    assert matches[0]['status'] == expected


def test_unrelated_extracted_feature_does_not_block_bp3_activation():
    bp3 = {
        'gene_id': 'GENE1', 'disease_id': 'DISEASE1',
        'variant_consequence': 'INFRAME_DELETION',
        'reading_frame_preserved': True,
        'biologically_relevant_transcript': True,
        'protein_consequence_confidence': True,
        'repeat_region_status': True,
        'repeat_region_functional_importance': 'NO_KNOWN_FUNCTION',
        'repeat_length_disease_mechanism': False,
    }
    without_extra = classify_variant({'criteria': {'BP3': bp3}})
    with_extra = classify_variant({'criteria': {'BP3': {**bp3, 'unrelated_annotation': 1}}})
    assert with_extra['criteria_activated'] == without_extra['criteria_activated']
