from copy import deepcopy

from adaptive_ps3 import evaluate_ps3

from .ps3_fixtures import assay, evidence


def test_missing_assay_controls_is_undetermined():
    raw = assay(positive_control_count=None, known_pathogenic_control_results=None)
    result = evaluate_ps3(evidence(raw))
    assert result.status.value == "UNDETERMINED"
    assert "positive_control_count" in result.missing_critical_features


def test_missing_essential_disease_mechanism_is_undetermined():
    raw = evidence()
    raw.pop("disease_mechanism")
    result = evaluate_ps3(raw)
    assert result.status.value == "UNDETERMINED"
    assert "disease_mechanism" in result.missing_critical_features


def test_missing_optional_features_do_not_block_activation():
    # No localization, binding, expression, blinding, sensitivity, or specificity.
    result = evaluate_ps3(evidence())
    assert result.status.value == "ACTIVATED"


def test_minimum_complete_pathway_does_not_require_inventory_exhaustion():
    minimal = {
        "gene_id": "GENE1",
        "disease_id": "DISEASE1",
        "disease_mechanism": "loss_of_function",
        "functional_assays": [{
            "assay_measures_variant_effect": True,
            "assay_measures_disease_mechanism": True,
            "variant_result_classification": "abnormal",
            "assay_validation_status": "historically_accepted",
        }],
    }
    result = evaluate_ps3(minimal)
    assert (result.status.value, result.strength.value) == ("ACTIVATED", "SUPPORTING")


def test_irrelevant_features_do_not_change_decision():
    baseline = evidence()
    enriched = deepcopy(baseline)
    enriched["functional_assays"][0].update({
        "altered_protein_localization": "unknown",
        "experimental_blinding": True,
        "protein_expression_level": 0.47,
        "assay_protocol": "synthetic protocol metadata",
    })
    first = evaluate_ps3(baseline)
    second = evaluate_ps3(enriched)
    assert (first.status, first.strength) == (second.status, second.strength)


def test_no_assay_is_undetermined_not_false_negative():
    result = evaluate_ps3(evidence(functional_assays=[]))
    assert result.status.value == "UNDETERMINED"
    assert result.missing_critical_features == ("functional_assays",)
