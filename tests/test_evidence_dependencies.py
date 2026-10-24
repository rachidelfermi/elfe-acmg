from acmg.evidence_dependencies import resolve_dependencies


def row(code, **extra):
    return {'criterion':code,'status':'ACTIVATED','strength':'SUPPORTING',**extra}


def test_shared_prediction_cannot_count_both_directions():
    rows = [row('PP3', predictor_used='REVEL'), row('BP4', predictor_used='REVEL')]
    counted, conflicts, resolved = resolve_dependencies(rows, {})
    assert counted == [] and not conflicts and resolved


def test_independent_functional_conflict_is_preserved():
    rows = [row('PS3', assays_used=['a']), row('BS3', assays_used=['b'])]
    counted, conflicts, resolved = resolve_dependencies(rows, {})
    assert len(counted) == 2 and conflicts and not resolved


def test_splice_overlap_suppresses_only_reused_pp3():
    rows = [row('PVS1'),row('PP3')]
    counted, _, _ = resolve_dependencies(rows, {'PP3':{'pvs1_splice_overlap':True}})
    assert [r['criterion'] for r in counted] == ['PVS1']


def test_ba1_supersedes_bs1():
    rows = [row('BA1'),row('BS1')]
    counted, _, resolved = resolve_dependencies(rows, {})
    assert [r['criterion'] for r in counted] == ['BA1']
    assert resolved[0]['unresolved'] is False
