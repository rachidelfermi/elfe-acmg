from copy import deepcopy


def assay(**overrides):
    value = {
        "assay_id": "assay-1",
        "assay_type": "quantitative",
        "assay_measures_variant_effect": True,
        "assay_measures_disease_mechanism": True,
        "assay_validation_status": "validated",
        "variant_result_classification": "abnormal",
        "variant_result_within_assay_range": True,
        "positive_control_count": 1,
        "negative_control_count": 1,
        "technical_replicates": 3,
        "known_pathogenic_control_results": ["P1", "P2", "P3", "P4", "P5"],
        "known_benign_control_results": ["B1", "B2", "B3", "B4", "B5", "B6"],
        "assay_validation_sample_size": 11,
        "control_variant_classification_quality": True,
        "control_classification_independence": True,
    }
    value.update(overrides)
    return value


def evidence(*assays, **overrides):
    value = {
        "variant_id": "v1",
        "gene_id": "GENE1",
        "disease_id": "DISEASE1",
        "variant_consequence": "missense_variant",
        "disease_mechanism": "loss_of_function",
        "functional_assays": list(assays or (assay(),)),
    }
    value.update(overrides)
    return value


def formal_assay(oddspath=20, **overrides):
    value = assay(
        assay_validation_status="formally_calibrated",
        assay_odds_pathogenicity=oddspath,
        assay_validation_method="formal OddsPath analysis",
    )
    value.update(overrides)
    return value


def supporting_assay(**overrides):
    value = assay(
        known_pathogenic_control_results=["P1", "P2", "P3"],
        known_benign_control_results=["B1", "B2", "B3"],
        assay_validation_sample_size=6,
        assay_validation_status="limited",
    )
    value.update(overrides)
    return value


def clone(value):
    return deepcopy(value)


def vcep_spec():
    return {
        "name": "Synthetic pinned VCEP PS3 specification",
        "version": "1.0",
        "specification_id": "SYN-PS3-1",
        "source_url": "https://example.invalid/vcep/syn-ps3-1",
        "gene_id": "GENE1",
        "disease_id": "DISEASE1",
        "rules": [
            {
                "rule_id": "accepted-functional-assay",
                "when": {"assay_name": "VCEP accepted", "variant_result_classification": "abnormal"},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "STRONG",
                    "decision_path": "VCEP_ACCEPTED_ASSAY",
                    "reason": "The pinned VCEP specification accepts this assay at Strong.",
                },
            },
            {
                "rule_id": "rejected-functional-assay",
                "when": {"assay_name": "VCEP rejected"},
                "decision": {
                    "status": "NOT_MET",
                    "strength": "NONE",
                    "decision_path": "VCEP_REJECTED_ASSAY",
                    "reason": "The pinned VCEP specification rejects this assay.",
                },
            },
            {
                "rule_id": "moderate-adjustment",
                "when": {"assay_name": "VCEP moderate", "variant_result_normalized": {"operator": "<=", "value": 0.25}},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "MODERATE",
                    "decision_path": "VCEP_THRESHOLD_MODERATE",
                },
            },
        ],
    }
