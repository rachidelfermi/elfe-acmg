"""Synthetic acceptance tests for the shared ClinGen SVI PS2/PM6 engine."""

from __future__ import annotations

from copy import deepcopy
import unittest

from ps2pm6 import SchemaValidationError, evaluate_ps2_pm6


def eligible_case(**overrides):
    value = {
        "proband_variant_present": True,
        "maternal_genotype": "ABSENT",
        "paternal_genotype": "ABSENT",
        "variant_call_quality": "ADEQUATE",
        "maternity_confirmed": True,
        "paternity_confirmed": True,
        "gene_disease_validity": "DEFINITIVE",
        "inheritance_pattern": "AUTOSOMAL_DOMINANT",
        "phenotype_category": "HIGHLY_SPECIFIC",
    }
    value.update(overrides)
    return value


class PS2PM6EngineTests(unittest.TestCase):
    def decision(self, **overrides):
        return evaluate_ps2_pm6(eligible_case(**overrides))

    def assert_outcome(self, result, status, strength, criterion=None, points=None):
        self.assertEqual(result.status.value, status)
        self.assertEqual(result.strength.value, strength)
        self.assertEqual(result.assigned_criterion, criterion)
        if points is not None:
            self.assertEqual(result.total_points, points)

    def test_confirmed_highly_specific_is_ps2_strong(self):
        self.assert_outcome(self.decision(), "ACTIVATED", "STRONG", "PS2", 2.0)

    def test_confirmed_consistent_is_ps2_moderate(self):
        self.assert_outcome(
            self.decision(phenotype_category="CONSISTENT"),
            "ACTIVATED", "MODERATE", "PS2", 1.0,
        )

    def test_confirmed_heterogeneous_is_ps2_supporting(self):
        self.assert_outcome(
            self.decision(phenotype_category="HETEROGENEOUS"),
            "ACTIVATED", "SUPPORTING", "PS2", 0.5,
        )

    def test_assumed_highly_specific_is_pm6_moderate(self):
        self.assert_outcome(
            self.decision(maternity_confirmed=False, paternity_confirmed=False),
            "ACTIVATED", "MODERATE", "PM6", 1.0,
        )

    def test_assumed_consistent_is_pm6_supporting(self):
        self.assert_outcome(
            self.decision(
                maternity_confirmed=False,
                paternity_confirmed=False,
                phenotype_category="CONSISTENT",
            ),
            "ACTIVATED", "SUPPORTING", "PM6", 0.5,
        )

    def test_assumed_heterogeneous_has_no_threshold(self):
        self.assert_outcome(
            self.decision(
                maternity_confirmed=False,
                paternity_confirmed=False,
                phenotype_category="HETEROGENEOUS",
            ),
            "NOT_MET", "NONE", None, 0.25,
        )

    def test_missing_maternal_or_paternal_test_is_never_assumed_de_novo(self):
        for field in ("maternal_genotype", "paternal_genotype"):
            raw = eligible_case()
            raw.pop(field)
            result = evaluate_ps2_pm6(raw)
            with self.subTest(field=field):
                self.assert_outcome(result, "UNDETERMINED", "NONE")
                self.assertIn(field, result.missing_critical_features)

    def test_other_proband_genotype_is_not_misread_as_a_reference_call(self):
        raw = eligible_case(proband_genotype="OTHER")
        raw.pop("proband_variant_present")
        result = evaluate_ps2_pm6(raw)
        self.assert_outcome(result, "UNDETERMINED", "NONE")
        self.assertIn("proband_variant_present", result.missing_critical_features)

    def test_explicit_case_unknown_is_not_overwritten_by_shared_context(self):
        shared = eligible_case()
        result = evaluate_ps2_pm6({
            **shared,
            "observations": [
                {"case_id": "one", "case_independence": True},
                {
                    "case_id": "two", "case_independence": True,
                    "maternal_genotype": {"state": "UNKNOWN"},
                    "paternal_genotype": {"state": "UNKNOWN"},
                },
            ],
        })
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PS2", 2.0)
        self.assertEqual(result.observations[1].status.value, "UNDETERMINED")

    def test_explicit_assumed_parentage_is_not_overwritten_by_shared_confirmed_parentage(self):
        shared = eligible_case()
        result = evaluate_ps2_pm6({
            **shared,
            "observations": [{
                "case_id": "one", "case_independence": True,
                "independent_case_parentage_status": "ASSUMED",
            }],
        })
        self.assert_outcome(result, "ACTIVATED", "MODERATE", "PM6", 1.0)

    def test_parentage_status_envelopes_preserve_unknown_and_conflicting_states(self):
        shared = eligible_case()
        unknown = evaluate_ps2_pm6({
            **shared,
            "observations": [{
                "case_id": "unknown-parentage",
                "independent_case_parentage_status": {"state": "UNKNOWN"},
            }],
        })
        self.assert_outcome(unknown, "ACTIVATED", "MODERATE", "PM6", 1.0)
        conflicting = evaluate_ps2_pm6({
            **shared,
            "observations": [{
                "case_id": "conflicting-parentage",
                "independent_case_parentage_status": {"state": "CONFLICTING"},
            }],
        })
        self.assert_outcome(conflicting, "UNDETERMINED", "NONE")

    def test_one_confirmed_parent_is_assumed_pm6(self):
        result = self.decision(maternity_confirmed=True, paternity_confirmed=False)
        self.assert_outcome(result, "ACTIVATED", "MODERATE", "PM6", 1.0)

    def test_inadequate_parentage_confirmation_downgrades_to_assumed_pm6(self):
        result = self.decision(parentage_confirmation_quality="INADEQUATE")
        self.assert_outcome(result, "ACTIVATED", "MODERATE", "PM6", 1.0)

    def test_flat_input_accepts_independent_case_aliases(self):
        raw = eligible_case(
            technical_quality="ADEQUATE",
            independent_case_parentage_status="CONFIRMED",
            independent_case_phenotype_specificity="HIGHLY_SPECIFIC",
        )
        raw.pop("variant_call_quality")
        raw.pop("phenotype_category")
        self.assert_outcome(evaluate_ps2_pm6(raw), "ACTIVATED", "STRONG", "PS2", 2.0)

    def test_unconfirmed_parentage_with_negative_parents_is_pm6(self):
        result = self.decision(maternity_confirmed="UNKNOWN", paternity_confirmed="UNKNOWN")
        self.assert_outcome(result, "ACTIVATED", "MODERATE", "PM6", 1.0)

    def test_inherited_from_either_parent_is_not_met(self):
        for field in ("maternal_genotype", "paternal_genotype"):
            with self.subTest(field=field):
                result = self.decision(**{field: "PRESENT"})
                self.assert_outcome(result, "NOT_MET", "NONE")

    def test_decisive_parental_transmission_is_not_masked_by_the_other_missing_parent(self):
        raw = eligible_case(paternal_genotype="PRESENT")
        raw.pop("maternal_genotype")
        self.assert_outcome(evaluate_ps2_pm6(raw), "NOT_MET", "NONE")

    def test_decisive_context_and_quality_exclusions_are_not_masked_by_missing_parental_calls(self):
        cases = (
            ("REFUTED", "AUTOSOMAL_DOMINANT", "ADEQUATE", "NOT_APPLICABLE"),
            ("DEFINITIVE", "OTHER", "ADEQUATE", "NOT_APPLICABLE"),
            ("DEFINITIVE", "AUTOSOMAL_DOMINANT", "POOR", "NOT_MET"),
        )
        for validity, inheritance, quality, expected in cases:
            raw = eligible_case(
                gene_disease_validity=validity,
                inheritance_pattern=inheritance,
                variant_call_quality=quality,
            )
            raw.pop("maternal_genotype")
            with self.subTest(validity=validity, inheritance=inheritance, quality=quality):
                self.assert_outcome(evaluate_ps2_pm6(raw), expected, "NONE")

    def test_parental_mosaicism_prevents_de_novo_activation(self):
        result = self.decision(parental_mosaicism_evidence=True)
        self.assert_outcome(result, "NOT_MET", "NONE")

    def test_proband_mosaicism_is_not_an_automatic_generic_exclusion(self):
        result = self.decision(proband_genotype="MOSAIC", proband_mosaicism_evidence=True)
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PS2", 2.0)

    def test_conflicting_proband_mosaicism_is_not_silently_ignored(self):
        result = self.decision(proband_mosaicism_evidence={"state": "CONFLICTING"})
        self.assert_outcome(result, "UNDETERMINED", "NONE")

    def test_explicitly_low_parental_depth_blocks_activation_without_inventing_a_numeric_cutoff(self):
        result = self.decision(parental_read_depth="LOW")
        self.assert_outcome(result, "NOT_MET", "NONE")
        result = self.decision(parental_read_depth={"mother": 7, "father": 8})
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PS2", 2.0)

    def test_explicitly_inadequate_parental_allele_balance_blocks_but_zero_does_not(self):
        for value in ("INADEQUATE", "LOW", {"quality": "FAIL"}):
            with self.subTest(value=value):
                self.assert_outcome(self.decision(parental_allele_balance=value), "NOT_MET", "NONE")
        self.assert_outcome(self.decision(parental_allele_balance={"mother": 0, "father": 0}), "ACTIVATED", "STRONG", "PS2", 2.0)

    def test_poor_variant_call_quality_blocks_activation(self):
        self.assert_outcome(self.decision(variant_call_quality="POOR"), "NOT_MET", "NONE")

    def test_explicitly_inadequate_proband_read_support_blocks_activation(self):
        for value in (0, "LOW", {"alt_reads": 0}):
            with self.subTest(value=value):
                self.assert_outcome(self.decision(proband_read_support=value), "NOT_MET", "NONE")

    def test_named_reference_or_parental_alternate_counts_do_not_masquerade_as_zero_coverage(self):
        self.assert_outcome(
            self.decision(proband_read_support={"reference_reads": 0, "alternate_reads": 10}),
            "ACTIVATED", "STRONG", "PS2", 2.0,
        )

    def test_raw_metrics_reject_nonfinite_numbers(self):
        for field, value in (
            ("proband_read_support", float("nan")),
            ("parental_read_depth", {"mother": float("inf")}),
            ("parental_allele_balance", [float("-inf")]),
        ):
            with self.subTest(field=field):
                with self.assertRaises(SchemaValidationError):
                    self.decision(**{field: value})
        self.assert_outcome(
            self.decision(parental_read_depth={
                "mother": {"total_depth": 30, "alternate_reads": 0},
                "father": {"total_depth": 30, "alternate_reads": 0},
            }),
            "ACTIVATED", "STRONG", "PS2", 2.0,
        )

    def test_sample_identity_conflict_is_undetermined(self):
        result = self.decision(sample_identity_confirmation=False)
        self.assert_outcome(result, "UNDETERMINED", "NONE")
        self.assertIn("sample_identity_confirmation", result.missing_critical_features)

    def test_explicit_nonbiological_family_relationship_phrases_block_activation(self):
        for structure in (
            "non-biological father",
            "non biological mother",
            "DONOR_EGG",
        ):
            with self.subTest(structure=structure):
                self.assert_outcome(
                    self.decision(family_relationship_structure=structure),
                    "NOT_MET",
                    "NONE",
                )

    def test_multiple_independent_confirmed_cases_reach_very_strong(self):
        first = eligible_case(case_id="one", case_independence=True)
        second = eligible_case(case_id="two", case_independence=True)
        result = evaluate_ps2_pm6({"observations": [first, second]})
        self.assert_outcome(result, "ACTIVATED", "VERY_STRONG", "PS2", 4.0)

    def test_multiple_independent_assumed_cases_are_pooled_as_pm6(self):
        first = eligible_case(case_id="one", case_independence=True, maternity_confirmed=False, paternity_confirmed=False)
        second = eligible_case(case_id="two", case_independence=True, maternity_confirmed=False, paternity_confirmed=False)
        result = evaluate_ps2_pm6({"observations": [first, second]})
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PM6", 2.0)

    def test_mixed_confirmed_and_assumed_cases_use_one_ps2_code(self):
        first = eligible_case(case_id="one", case_independence=True)
        second = eligible_case(case_id="two", case_independence=True, maternity_confirmed=False, paternity_confirmed=False)
        result = evaluate_ps2_pm6({"observations": [first, second]})
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PS2", 3.0)
        self.assertEqual(result.to_dict()["criteria"], ["PS2", "PM6"])
        self.assertNotIn("PM6", [result.assigned_criterion])

    def test_duplicate_published_or_patient_case_is_not_counted_twice(self):
        first = eligible_case(case_id="same", case_independence=True)
        duplicate = eligible_case(case_id="same", case_independence=True)
        result = evaluate_ps2_pm6({"observations": [first, duplicate]})
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PS2", 2.0)
        self.assertEqual(result.observations[1].status.value, "DUPLICATE")

    def test_matching_publication_case_identifier_prevents_double_counting_even_when_local_ids_differ(self):
        first = eligible_case(case_id="local-a", publication_case_id="paper-1-case-7", case_independence=True)
        duplicate = eligible_case(case_id="local-b", publication_case_id="paper-1-case-7", case_independence=True)
        result = evaluate_ps2_pm6({"observations": [first, duplicate]})
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PS2", 2.0)
        self.assertEqual(result.observations[1].status.value, "DUPLICATE")

    def test_duplicate_reports_merge_complementary_evidence_independent_of_input_order(self):
        incomplete = {"case_id": "same", "case_independence": True}
        complete = eligible_case(case_id="same", case_independence=True)
        first = evaluate_ps2_pm6({"observations": [incomplete, complete]})
        second = evaluate_ps2_pm6({"observations": [complete, incomplete]})
        self.assert_outcome(first, "ACTIVATED", "STRONG", "PS2", 2.0)
        self.assert_outcome(second, "ACTIVATED", "STRONG", "PS2", 2.0)
        self.assertEqual(first.observations[1].status.value, "DUPLICATE")
        self.assertEqual(second.observations[1].status.value, "DUPLICATE")

    def test_duplicate_reports_with_opposing_biological_calls_remain_undetermined(self):
        first = eligible_case(case_id="same", case_independence=True)
        second = eligible_case(case_id="same", case_independence=True, maternal_genotype="PRESENT")
        result = evaluate_ps2_pm6({"observations": [first, second]})
        self.assert_outcome(result, "UNDETERMINED", "NONE")

    def test_duplicate_report_independence_flag_does_not_disqualify_its_representative_case(self):
        first = eligible_case(case_id="same", case_independence=True)
        duplicate = eligible_case(case_id="same", case_independence=False)
        result = evaluate_ps2_pm6({"observations": [first, duplicate]})
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PS2", 2.0)
        self.assertEqual(result.observations[1].status.value, "DUPLICATE")
        self.assertEqual(result.observations[1].duplicate_of, 0)

    def test_undocumented_distinct_case_is_not_pooled_regardless_of_order(self):
        unknown = eligible_case(case_id="unknown")
        documented = eligible_case(case_id="documented", case_independence=True)
        for observations in ((unknown, documented), (documented, unknown)):
            with self.subTest(order=[item["case_id"] for item in observations]):
                self.assert_outcome(
                    evaluate_ps2_pm6({"observations": list(observations)}),
                    "ACTIVATED",
                    "STRONG",
                    "PS2",
                    2.0,
                )

    def test_generic_pool_rejects_mixed_known_gene_disease_or_inheritance_context(self):
        cases = (
            (
                eligible_case(case_id="one", case_independence=True, gene_id="G1", disease_id="D1"),
                eligible_case(case_id="two", case_independence=True, gene_id="G2", disease_id="D1"),
            ),
            (
                eligible_case(
                    case_id="one",
                    case_independence=True,
                    gene_id="G1",
                    disease_id="D1",
                    inheritance_pattern="AUTOSOMAL_DOMINANT",
                ),
                eligible_case(
                    case_id="two",
                    case_independence=True,
                    gene_id="G1",
                    disease_id="D1",
                    inheritance_pattern="AUTOSOMAL_RECESSIVE",
                    additional_pathogenic_variant_present=True,
                ),
            ),
        )
        for observations in cases:
            with self.subTest(context=[item.get("gene_id") for item in observations]):
                self.assert_outcome(
                    evaluate_ps2_pm6({"observations": list(observations)}),
                    "UNDETERMINED",
                    "NONE",
                )

    def test_generic_pool_compares_identifier_context_case_insensitively_without_collapsing_punctuation(self):
        first = eligible_case(
            case_id="one", case_independence=True,
            gene_id="HGNC:5", disease_id="MONDO:0000001",
        )
        second = eligible_case(
            case_id="two", case_independence=True,
            gene_id="hgnc:5", disease_id="mondo:0000001",
        )
        self.assert_outcome(
            evaluate_ps2_pm6({"observations": [first, second]}),
            "ACTIVATED",
            "VERY_STRONG",
            "PS2",
            4.0,
        )

    def test_local_and_publication_identifiers_do_not_merge_on_an_unqualified_collision(self):
        local = eligible_case(case_id="1", case_independence=True)
        publication = eligible_case(publication_case_id="1", case_independence=True)
        result = evaluate_ps2_pm6({"observations": [local, publication]})
        self.assert_outcome(result, "ACTIVATED", "VERY_STRONG", "PS2", 4.0)

    def test_merged_duplicate_provenance_is_deterministic_and_preserves_all_sources(self):
        first = eligible_case(
            case_id={
                "value": "local-a",
                "provenance": {"source": "laboratory-a", "identifier": "A"},
            },
            publication_case_id={
                "value": "paper-1",
                "provenance": {"source": "paper", "identifier": "P"},
            },
            case_independence=True,
        )
        second = eligible_case(
            case_id={
                "value": "local-b",
                "provenance": {"source": "laboratory-b", "identifier": "B"},
            },
            publication_case_id={
                "value": "paper-1",
                "provenance": {"source": "paper", "identifier": "P"},
            },
            case_independence=True,
        )
        forward = evaluate_ps2_pm6({"observations": [first, second]})
        reverse = evaluate_ps2_pm6({"observations": [second, first]})
        sources = {(item.source, item.identifier) for item in forward.evidence_provenance}
        self.assertTrue({("laboratory-a", "A"), ("laboratory-b", "B"), ("paper", "P")} <= sources)
        self.assertEqual(
            [item.to_dict() for item in forward.evidence_provenance],
            [item.to_dict() for item in reverse.evidence_provenance],
        )
        self.assertEqual(forward.observations[0].case_id, "local-a")
        self.assertEqual(reverse.observations[0].case_id, "local-b")

    def test_non_independent_case_is_not_counted_twice(self):
        first = eligible_case(case_id="one", case_independence=True)
        second = eligible_case(case_id="two", case_independence=False)
        result = evaluate_ps2_pm6({"observations": [first, second]})
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PS2", 2.0)

    def test_first_non_independent_case_is_not_order_dependently_counted(self):
        first = eligible_case(case_id="one", case_independence=False)
        second = eligible_case(case_id="two", case_independence=True)
        result = evaluate_ps2_pm6({"observations": [first, second]})
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PS2", 2.0)
        self.assertEqual(result.observations[0].status.value, "DUPLICATE")

    def test_conflicting_phenotype_evidence_is_undetermined(self):
        result = self.decision(phenotype_category="HIGHLY_SPECIFIC", phenotype_consistency=False)
        self.assert_outcome(result, "UNDETERMINED", "NONE")

    def test_native_phenotype_components_select_the_heterogeneous_svi_row(self):
        raw = eligible_case(
            phenotype_consistency=True,
            phenotype_specificity="CONSISTENT",
            phenotype_heterogeneity=True,
        )
        raw.pop("phenotype_category")
        result = evaluate_ps2_pm6(raw)
        self.assert_outcome(result, "ACTIVATED", "SUPPORTING", "PS2", 0.5)
        self.assertEqual(result.observations[0].phenotype_category.value, "HETEROGENEOUS")

    def test_positive_consistency_without_a_specificity_claim_uses_the_conservative_row(self):
        raw = eligible_case(phenotype_consistency=True)
        raw.pop("phenotype_category")
        result = evaluate_ps2_pm6(raw)
        self.assert_outcome(result, "ACTIVATED", "MODERATE", "PS2", 1.0)

    def test_conflicting_family_structure_is_undetermined(self):
        result = self.decision(family_relationship_structure={"state": "CONFLICTING"})
        self.assert_outcome(result, "UNDETERMINED", "NONE")

    def test_conflicting_parentage_method_is_undetermined(self):
        result = self.decision(parentage_confirmation_method={"state": "CONFLICTING"})
        self.assert_outcome(result, "UNDETERMINED", "NONE")
        self.assertIn("parentage_confirmation_method", result.missing_critical_features)

    def test_inconsistent_family_history_is_not_met(self):
        self.assert_outcome(self.decision(family_history_consistent=False), "NOT_MET", "NONE")

    def test_incompatible_phenotype_is_not_met(self):
        result = self.decision(phenotype_category="INCONSISTENT")
        self.assert_outcome(result, "NOT_MET", "NONE")

    def test_incompatible_generic_inheritance_is_not_applicable(self):
        result = self.decision(inheritance_pattern="OTHER")
        self.assert_outcome(result, "NOT_APPLICABLE", "NONE")

    def test_ar_without_second_pathogenic_allele_downgrades_one_level(self):
        result = self.decision(
            inheritance_pattern="AUTOSOMAL_RECESSIVE",
            additional_pathogenic_variant_present=False,
        )
        self.assert_outcome(result, "ACTIVATED", "MODERATE", "PS2", 2.0)

    def test_ar_requires_an_explicit_additional_allele_assessment(self):
        result = self.decision(inheritance_pattern="AUTOSOMAL_RECESSIVE")
        self.assert_outcome(result, "UNDETERMINED", "NONE")
        self.assertIn("additional_pathogenic_variant_present", result.missing_critical_features)

    def test_heterogeneous_phenotype_contributions_are_capped_at_one_point(self):
        observations = [
            eligible_case(case_id=str(index), case_independence=True, phenotype_category="HETEROGENEOUS")
            for index in range(3)
        ]
        result = evaluate_ps2_pm6({"observations": observations})
        self.assert_outcome(result, "ACTIVATED", "MODERATE", "PS2", 1.0)
        self.assertEqual(result.observations[-1].contributed_points, 0.0)

    def test_confirmed_case_still_selects_ps2_when_the_heterogeneous_cap_trims_its_contribution(self):
        observations = [
            eligible_case(
                case_id=f"assumed-{index}", case_independence=True,
                maternity_confirmed=False, paternity_confirmed=False,
                phenotype_category="HETEROGENEOUS",
            )
            for index in range(4)
        ]
        observations.append(eligible_case(case_id="confirmed", case_independence=True, phenotype_category="HETEROGENEOUS"))
        result = evaluate_ps2_pm6({"observations": observations})
        self.assert_outcome(result, "ACTIVATED", "MODERATE", "PS2", 1.0)
        self.assertEqual(result.confirmed_de_novo_cases, 1)

    def test_apparent_germline_mosaicism_requires_confirmed_parentage(self):
        result = self.decision(
            apparent_germline_mosaicism=True,
            maternity_confirmed=False,
            paternity_confirmed=False,
        )
        self.assert_outcome(result, "NOT_MET", "NONE")

    def test_vcep_point_modification_is_scoped_versioned_and_traceable(self):
        specification = {
            "name": "Synthetic PS2/PM6 VCEP rule",
            "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "specification_id": "SYNTHETIC-PS2PM6",
            "gene_id": "GENE1",
            "disease_id": "DISEASE1",
            "rules": [{
                "source_locator": "Synthetic published table",
                "when": {"phenotype_category": "HIGHLY_SPECIFIC"},
                "point_overrides": {"HIGHLY_SPECIFIC": {"CONFIRMED": 1, "ASSUMED": 0.5}},
            }],
        }
        result = evaluate_ps2_pm6(eligible_case(
            gene_id="GENE1", disease_id="DISEASE1",
            applicable_vcep_de_novo_rules=[specification],
        ))
        self.assert_outcome(result, "ACTIVATED", "MODERATE", "PS2", 1.0)
        self.assertEqual(result.guideline.version, "1.0.0")
        self.assertEqual(result.guideline.rule_id, "SYNTHETIC-PS2PM6")

    def test_vcep_scope_mismatch_does_not_override_generic_rule(self):
        specification = {
            "name": "Synthetic PS2/PM6 VCEP rule",
            "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1", "disease_id": "OTHER",
            "rules": [{"source_locator": "rule", "when": {}, "point_overrides": {"HIGHLY_SPECIFIC": {"CONFIRMED": 1}}}],
        }
        result = evaluate_ps2_pm6(eligible_case(gene_id="GENE1", disease_id="DISEASE1", applicable_vcep_de_novo_rules=[specification]))
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PS2", 2.0)

    def test_vcep_scope_accepts_inheritance_pattern_alias(self):
        specification = {
            "name": "Synthetic inheritance-scoped VCEP rule",
            "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1", "disease_id": "DISEASE1",
            "inheritance_pattern": "AD",
            "rules": [{"source_locator": "rule", "when": {}, "point_overrides": {"HIGHLY_SPECIFIC": {"CONFIRMED": 1}}}],
        }
        result = evaluate_ps2_pm6(eligible_case(gene_id="GENE1", disease_id="DISEASE1", applicable_vcep_de_novo_rules=[specification]))
        self.assert_outcome(result, "ACTIVATED", "MODERATE", "PS2", 1.0)

    def test_vcep_conditions_normalize_documented_observation_aliases_and_values(self):
        specification = {
            "name": "Alias-aware VCEP rule", "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1", "disease_id": "DISEASE1",
            "rules": [{
                "source_locator": "alias row",
                "when": {
                    "maternal_genotype": "NEGATIVE",
                    "technical_quality": "GOOD",
                    "mode_of_inheritance": "AD",
                },
                "point_overrides": {"HIGHLY_SPECIFIC": {"CONFIRMED": 1}},
            }],
        }
        result = evaluate_ps2_pm6(eligible_case(
            gene_id="GENE1",
            disease_id="DISEASE1",
            applicable_vcep_de_novo_rules=[specification],
        ))
        self.assert_outcome(result, "ACTIVATED", "MODERATE", "PS2", 1.0)

    def test_vcep_rejects_repeated_normalized_table_keys(self):
        base = {
            "name": "Invalid VCEP rule", "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1", "disease_id": "DISEASE1",
        }
        repeated_category = {
            **base,
            "rules": [{
                "source_locator": "repeated category", "when": {},
                "point_overrides": {
                    "HIGHLY_SPECIFIC": {"CONFIRMED": 1},
                    "highly-specific": {"CONFIRMED": 3},
                },
            }],
        }
        repeated_threshold = {
            **base,
            "rules": [{
                "source_locator": "repeated threshold", "when": {},
                "strength_thresholds": {"SUPPORTING": 0.5, "supporting": 0.75},
            }],
        }
        for specification in (repeated_category, repeated_threshold):
            with self.subTest(rule=specification["rules"][0]["source_locator"]):
                with self.assertRaises(SchemaValidationError):
                    evaluate_ps2_pm6(eligible_case(
                        gene_id="GENE1",
                        disease_id="DISEASE1",
                        applicable_vcep_de_novo_rules=[specification],
                    ))

    def test_vcep_rejects_invalid_effective_thresholds_and_zero_thresholds(self):
        base = {
            "name": "Invalid threshold VCEP rule", "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1", "disease_id": "DISEASE1",
        }
        non_monotonic = {
            **base,
            "rules": [{
                "source_locator": "invalid partial ladder", "when": {},
                "strength_thresholds": {"MODERATE": 0.25},
            }],
        }
        zero = {
            **base,
            "rules": [{
                "source_locator": "zero threshold", "when": {},
                "strength_thresholds": {"SUPPORTING": 0},
            }],
        }
        for specification in (non_monotonic, zero):
            with self.subTest(specification=specification["rules"][0]["source_locator"]):
                with self.assertRaises(SchemaValidationError):
                    evaluate_ps2_pm6(eligible_case(
                        gene_id="GENE1", disease_id="DISEASE1",
                        applicable_vcep_de_novo_rules=[specification],
                    ))

    def test_vcep_rejects_malformed_condition_even_when_its_scope_does_not_match(self):
        specification = {
            "name": "Invalid VCEP rule", "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1", "disease_id": "OTHER",
            "rules": [{
                "source_locator": "invalid condition", "when": {"phenotype_category": "TYPO"},
            }],
        }
        with self.assertRaises(SchemaValidationError):
            evaluate_ps2_pm6(eligible_case(applicable_vcep_de_novo_rules=[specification]))

    def test_vcep_rejects_conflicting_canonical_and_alias_rule_fields(self):
        specification = {
            "name": "Invalid VCEP rule", "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1", "disease_id": "DISEASE1",
            "rules": [{
                "source_locator": "contradictory aliases", "when": {},
                "point_overrides": {"HIGHLY_SPECIFIC": {"CONFIRMED": 1}},
                "points": {"HIGHLY_SPECIFIC": {"CONFIRMED": 2}},
            }],
        }
        with self.assertRaises(SchemaValidationError):
            evaluate_ps2_pm6(eligible_case(
                gene_id="GENE1", disease_id="DISEASE1",
                applicable_vcep_de_novo_rules=[specification],
            ))

    def test_vcep_unknown_condition_does_not_win_over_a_later_known_mismatch(self):
        specification = {
            "name": "Synthetic VCEP rule", "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1", "disease_id": "DISEASE1",
            "rules": [{
                "source_locator": "nonmatching ordered condition",
                "when": {
                    "family_history_consistent": True,
                    "phenotype_category": "CONSISTENT",
                },
                "point_overrides": {"HIGHLY_SPECIFIC": {"CONFIRMED": 1}},
            }],
        }
        result = evaluate_ps2_pm6(eligible_case(
            gene_id="GENE1", disease_id="DISEASE1",
            applicable_vcep_de_novo_rules=[specification],
        ))
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PS2", 2.0)

    def test_vcep_cannot_relabel_assumed_parentage_as_ps2(self):
        specification = {
            "name": "Synthetic VCEP rule", "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1", "disease_id": "DISEASE1",
            "rules": [{"source_locator": "invalid label", "when": {}, "assigned_criterion": "PS2"}],
        }
        result = evaluate_ps2_pm6(eligible_case(
            gene_id="GENE1", disease_id="DISEASE1",
            maternity_confirmed=False, paternity_confirmed=False,
            applicable_vcep_de_novo_rules=[specification],
        ))
        self.assert_outcome(result, "UNDETERMINED", "NONE")

    def test_vcep_policies_from_distinct_scopes_cannot_be_pooled(self):
        def specification(gene, disease):
            return {
                "name": f"{gene} VCEP", "version": "1.0.0",
                "source_url": "https://example.invalid/ps2pm6/v1",
                "gene_id": gene, "disease_id": disease,
                "rules": [{"source_locator": "row", "when": {}}],
            }

        first = eligible_case(case_id="one", case_independence=True, gene_id="G1", disease_id="D1")
        second = eligible_case(case_id="two", case_independence=True, gene_id="G2", disease_id="D2")
        result = evaluate_ps2_pm6({
            "observations": [first, second],
            "applicable_vcep_de_novo_rules": [specification("G1", "D1"), specification("G2", "D2")],
        })
        self.assert_outcome(result, "UNDETERMINED", "NONE")

    def test_vcep_not_applicable_decision_preserves_its_status_and_source(self):
        specification = {
            "name": "Synthetic PS2/PM6 VCEP rule",
            "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "specification_id": "SYNTHETIC-PS2PM6",
            "gene_id": "GENE1",
            "disease_id": "DISEASE1",
            "rules": [{
                "source_locator": "VCEP out-of-scope row",
                "when": {},
                "decision": {"status": "NOT_APPLICABLE"},
            }],
        }
        result = evaluate_ps2_pm6(eligible_case(
            gene_id="GENE1", disease_id="DISEASE1",
            applicable_vcep_de_novo_rules=[specification],
        ))
        self.assert_outcome(result, "NOT_APPLICABLE", "NONE")
        self.assertEqual(result.guideline.version, "1.0.0")
        self.assertEqual(result.observations[0].guideline_locator, "VCEP out-of-scope row")

    def test_resolvable_vcep_exclusion_is_not_masked_by_missing_generic_activation_facts(self):
        specification = {
            "name": "Synthetic PS2/PM6 VCEP rule",
            "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1",
            "disease_id": "DISEASE1",
            "rules": [{
                "source_locator": "VCEP out-of-scope row",
                "when": {},
                "decision": {"status": "NOT_APPLICABLE"},
            }],
        }
        raw = eligible_case(
            gene_id="GENE1",
            disease_id="DISEASE1",
            applicable_vcep_de_novo_rules=[specification],
        )
        raw.pop("maternal_genotype")
        result = evaluate_ps2_pm6(raw)
        self.assert_outcome(result, "NOT_APPLICABLE", "NONE")
        self.assertEqual(result.guideline.version, "1.0.0")

    def test_duplicate_audit_record_does_not_change_not_applicable_terminal_status(self):
        specification = {
            "name": "Synthetic PS2/PM6 VCEP rule",
            "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1",
            "disease_id": "DISEASE1",
            "rules": [{
                "source_locator": "VCEP out-of-scope row",
                "when": {},
                "decision": {"status": "NOT_APPLICABLE"},
            }],
        }
        first = eligible_case(
            case_id="same", case_independence=True,
            gene_id="GENE1", disease_id="DISEASE1",
        )
        duplicate = eligible_case(
            case_id="same", case_independence=True,
            gene_id="GENE1", disease_id="DISEASE1",
        )
        result = evaluate_ps2_pm6({
            "observations": [first, duplicate],
            "applicable_vcep_de_novo_rules": [specification],
        })
        self.assert_outcome(result, "NOT_APPLICABLE", "NONE")

    def test_vcep_point_value_below_threshold_is_not_rounded_up(self):
        specification = {
            "name": "Synthetic precise-point VCEP rule",
            "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1", "disease_id": "DISEASE1",
            "rules": [{
                "source_locator": "precise row", "when": {},
                "point_overrides": {"HIGHLY_SPECIFIC": {"CONFIRMED": 0.49996}},
            }],
        }
        result = evaluate_ps2_pm6(eligible_case(
            gene_id="GENE1", disease_id="DISEASE1", applicable_vcep_de_novo_rules=[specification],
        ))
        self.assert_outcome(result, "NOT_MET", "NONE", None, 0.49996)

    def test_vcep_rejects_nonfinite_point_values(self):
        specification = {
            "name": "Invalid VCEP rule", "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1", "disease_id": "DISEASE1",
            "rules": [{
                "source_locator": "invalid", "when": {},
                "point_overrides": {"HIGHLY_SPECIFIC": {"CONFIRMED": float("inf")}},
            }],
        }
        with self.assertRaises(SchemaValidationError):
            evaluate_ps2_pm6(eligible_case(gene_id="GENE1", disease_id="DISEASE1", applicable_vcep_de_novo_rules=[specification]))

    def test_potentially_applicable_vcep_with_unknown_scope_blocks_generic_fallback(self):
        specification = {
            "name": "Synthetic PS2/PM6 VCEP rule",
            "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1", "disease_id": "DISEASE1",
            "rules": [{"source_locator": "rule", "when": {}, "point_overrides": {"HIGHLY_SPECIFIC": {"CONFIRMED": 1}}}],
        }
        raw = eligible_case(applicable_vcep_de_novo_rules=[specification])
        result = evaluate_ps2_pm6(raw)
        self.assert_outcome(result, "UNDETERMINED", "NONE")
        self.assertIn("gene_id", result.missing_critical_features)

    def test_unresolved_more_specific_vcep_condition_cannot_be_bypassed_by_broad_match(self):
        specification = {
            "name": "Synthetic PS2/PM6 VCEP rule",
            "version": "1.0.0",
            "source_url": "https://example.invalid/ps2pm6/v1",
            "gene_id": "GENE1", "disease_id": "DISEASE1",
            "rules": [
                {"source_locator": "broad", "when": {}, "point_overrides": {"HIGHLY_SPECIFIC": {"CONFIRMED": 1}}},
                {"source_locator": "restricted", "when": {"family_history_consistent": True}, "exclude": True},
            ],
        }
        result = evaluate_ps2_pm6(eligible_case(
            gene_id="GENE1", disease_id="DISEASE1", applicable_vcep_de_novo_rules=[specification],
        ))
        self.assert_outcome(result, "UNDETERMINED", "NONE")
        self.assertIn("family_history_consistent", result.missing_critical_features)

    def test_missing_optional_fields_do_not_block_minimum_sufficient_pathway(self):
        raw = eligible_case()
        before = deepcopy(raw)
        result = evaluate_ps2_pm6(raw)
        self.assertEqual(raw, before)
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PS2", 2.0)
        self.assertNotIn("parental_read_depth", result.missing_critical_features)

    def test_missing_inheritance_context_does_not_block_the_generic_svi_table(self):
        raw = eligible_case()
        raw.pop("inheritance_pattern")
        self.assert_outcome(evaluate_ps2_pm6(raw), "ACTIVATED", "STRONG", "PS2", 2.0)

    def test_missing_essential_quality_does_not_award_points(self):
        raw = eligible_case()
        raw.pop("variant_call_quality")
        result = evaluate_ps2_pm6(raw)
        self.assert_outcome(result, "UNDETERMINED", "NONE", None, 0.0)
        self.assertIn("variant_call_quality", result.missing_critical_features)

    def test_explicit_unknown_and_conflicting_states_are_not_coerced_to_false(self):
        result = self.decision(maternal_genotype={"state": "UNKNOWN"})
        self.assert_outcome(result, "UNDETERMINED", "NONE")

    def test_contradictory_boolean_envelope_is_rejected(self):
        with self.assertRaises(SchemaValidationError):
            self.decision(maternity_confirmed={"state": "FALSE", "value": True})
        with self.assertRaises(SchemaValidationError):
            self.decision(maternity_confirmed={"state": "KNOWN", "status": "FALSE", "value": True})
        result = self.decision(paternal_genotype={"state": "CONFLICTING"})
        self.assert_outcome(result, "UNDETERMINED", "NONE")

    def test_provenance_is_preserved_for_a_decision_critical_feature(self):
        result = self.decision(proband_variant_present={
            "state": "TRUE",
            "provenance": {"source": "validated trio", "kind": "OBSERVED", "identifier": "TRIO-1"},
        })
        self.assert_outcome(result, "ACTIVATED", "STRONG", "PS2")
        self.assertTrue(any(item.identifier == "TRIO-1" for item in result.evidence_provenance))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
