from adaptive_pm3 import evaluate_pm3
from tests.pm3_fixtures import case, data


def homo(**changes):
    values = dict(candidate_genotype="HOM_ALT", second_variant=None,
                  homozygous_observation=True, homozygous_genotype_valid=True)
    values.update(changes)
    return case(**values)


def test_homozygous_supporting_and_deletion_exclusion():
    result = evaluate_pm3(data(homo()))
    assert result["total_points"] == 0.5 and result["strength"] == "SUPPORTING"
    assert evaluate_pm3(data(homo(deletion_mimicking_homozygosity=True)))["status"] == "NOT_MET"


def test_homozygous_validation_and_cap():
    assert evaluate_pm3(data(homo(homozygous_genotype_valid=None)))["status"] == "UNDETERMINED"
    cases = [homo(case_id=f"c{i}", family_id=f"f{i}", independence_confirmed=True) for i in range(3)]
    assert evaluate_pm3(data(*cases))["total_points"] == 1.0
