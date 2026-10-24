"""Acceptance tests for the evidence-adaptive PVS1 evaluator.

The fixtures in this module are deliberately synthetic.  They exercise the
published generic ClinGen SVI branches and the public API without requiring a
network connection, a transcript database, or a built-in VCEP catalogue.
"""

from __future__ import annotations

from copy import deepcopy
import unittest

from pvs1 import evaluate_pvs1


class PVS1TestCase(unittest.TestCase):
    """Shared synthetic evidence and assertions for public-API tests."""

    def eligible_null(self, **overrides):
        """Return the smallest ordinary null-variant fixture that can activate.

        Deliberately do not populate CNV, splicing, start-loss, or protein
        fields here.  Their absence is part of the minimum-evidence test.
        """

        evidence = {
            "gene_id": "GENE-TEST",
            "disease_id": "DISEASE-TEST",
            "transcript_id": "TX-TEST.1",
            "variant_consequence": "frameshift_variant",
            "lof_disease_mechanism": "established",
            "gene_disease_validity": "definitive",
            "biologically_relevant_transcript": True,
            "nmd_prediction": True,
        }
        evidence.update(overrides)
        return evidence

    def cnv_evidence(self, **overrides):
        """Return an eligible fixture without a small-variant consequence."""

        evidence = self.eligible_null()
        evidence.pop("variant_consequence")
        evidence.pop("nmd_prediction")
        evidence.update(overrides)
        return evidence

    def start_loss_evidence(self, **overrides):
        evidence = self.eligible_null(
            variant_consequence="start_lost",
            alternative_start_codon=False,
            upstream_pathogenic_variant_evidence=False,
        )
        evidence.pop("nmd_prediction")
        evidence.update(overrides)
        return evidence

    def evaluate(self, evidence):
        try:
            return evaluate_pvs1(evidence)
        except Exception as error:  # pragma: no cover - gives useful test output
            self.fail(f"evaluate_pvs1 unexpectedly raised {error!r} for {evidence!r}")

    def assert_decision(
        self,
        evidence,
        status,
        strength,
        *,
        missing=(),
        require_audit=False,
    ):
        decision = self.evaluate(evidence)
        self.assertEqual(status, decision.status.value)
        self.assertEqual(strength, decision.strength.value)

        if missing:
            # Each vector is built so that there is one smallest sufficient
            # request for additional evidence.  Do not let unrelated unknown
            # fields leak into the adaptive missing-evidence result.
            self.assertEqual(tuple(missing), tuple(decision.features_missing))

        if status == "ACTIVATED":
            self.assertEqual([], list(decision.features_missing))

        payload = decision.to_dict()
        self.assertEqual("PVS1", payload["criterion"])
        self.assertEqual(status, payload["status"])
        self.assertEqual(strength, payload["strength"])
        self.assertIn("features_used", payload)
        self.assertIn("features_missing", payload)
        self.assertIn("exclusions_evaluated", payload)
        self.assertIn("evidence_provenance", payload)
        self.assertIn("reason", payload)

        if require_audit:
            self.assertTrue(payload["decision_path"])
            guideline = payload["guideline"]
            self.assertTrue(guideline["name"])
            self.assertTrue(guideline["version"])
            # A generic published flowchart can identify a branch by its
            # source locator rather than an invented rule ID.  VCEP fixtures
            # may additionally carry a pinned specification ID.
            self.assertIn("rule_id", guideline)
            self.assertTrue(guideline["rule_id"])
            self.assertTrue(guideline["source_url"])
            self.assertTrue(guideline["source_locator"])
            self.assertTrue(payload["features_used"])

        return decision


class GenericNullVariantTests(PVS1TestCase):
    def test_frameshift_with_predicted_nmd_uses_very_strong_minimum_evidence(self):
        decision = self.assert_decision(
            self.eligible_null(),
            "ACTIVATED",
            "VERY_STRONG",
            require_audit=True,
        )
        self.assertEqual([], list(decision.features_missing))

    def test_curated_lof_conclusion_does_not_require_raw_mechanism_table_fields(self):
        evidence = self.eligible_null()
        evidence.pop("gene_disease_validity")
        self.assert_decision(evidence, "ACTIVATED", "VERY_STRONG")

    def test_single_exon_gene_can_satisfy_the_svi_table_1_fallback(self):
        evidence = self.eligible_null(
            pathogenic_lof_variant_count=3,
            pathogenic_lof_variant_fraction_percent=11,
            single_exon_gene=True,
        )
        evidence.pop("lof_disease_mechanism")
        self.assert_decision(evidence, "ACTIVATED", "VERY_STRONG")

    def test_negative_prior_pathogenic_lof_count_is_rejected(self):
        with self.assertRaisesRegex(ValueError, r"greater than or equal to zero"):
            evaluate_pvs1({"pathogenic_lof_variant_count": -1})

    def test_transcript_selection_policy_can_establish_relevance(self):
        evidence = self.eligible_null(
            transcript_selection_policy="disease_curated",
        )
        evidence.pop("biologically_relevant_transcript")
        self.assert_decision(evidence, "ACTIVATED", "VERY_STRONG")

    def test_curated_lof_mechanism_strength_caps_are_applied_without_scoring(self):
        with self.subTest("one-level cap"):
            self.assert_decision(
                self.eligible_null(lof_disease_mechanism="cap_strong"),
                "ACTIVATED",
                "STRONG",
            )
        with self.subTest("two-level cap"):
            self.assert_decision(
                self.eligible_null(lof_disease_mechanism="cap_moderate"),
                "ACTIVATED",
                "MODERATE",
            )

    def test_stop_gain_in_final_exon_with_small_truncation_is_moderate(self):
        evidence = self.eligible_null(
            variant_consequence="stop_gained",
            nmd_prediction=False,
            exon_location="final_exon",
            critical_protein_region=False,
            population_evidence_for_exon_lof_tolerance=False,
            protein_truncation_percent=5,
        )
        self.assert_decision(evidence, "ACTIVATED", "MODERATE", require_audit=True)

    def test_frameshift_without_nmd_uses_protein_loss_ladder(self):
        with self.subTest("greater than ten percent"):
            self.assert_decision(
                self.eligible_null(
                    nmd_prediction=False,
                    critical_protein_region=False,
                    population_evidence_for_exon_lof_tolerance=False,
                    protein_truncation_percent=15,
                ),
                "ACTIVATED",
                "STRONG",
            )

        with self.subTest("less than ten percent"):
            self.assert_decision(
                self.eligible_null(
                    nmd_prediction=False,
                    critical_protein_region=False,
                    population_evidence_for_exon_lof_tolerance=False,
                    protein_truncation_percent=5,
                ),
                "ACTIVATED",
                "MODERATE",
            )

    def test_nmd_escape_critical_region_is_strong(self):
        self.assert_decision(
            self.eligible_null(
                nmd_prediction=False,
                critical_protein_region=True,
            ),
            "ACTIVATED",
            "STRONG",
            require_audit=True,
        )

    def test_exact_ten_percent_boundary_is_not_silently_invented(self):
        """SVI 2018 states >10% and <10%, not an unqualified <= rule."""

        self.assert_decision(
            self.eligible_null(
                nmd_prediction=False,
                critical_protein_region=False,
                population_evidence_for_exon_lof_tolerance=False,
                protein_truncation_percent=10,
            ),
            "UNDETERMINED",
            "NONE",
        )

    def test_missing_lof_mechanism_is_undetermined_not_negative(self):
        evidence = self.eligible_null()
        evidence.pop("lof_disease_mechanism")
        self.assert_decision(
            evidence,
            "UNDETERMINED",
            "NONE",
            missing=("lof_disease_mechanism",),
        )

    def test_known_absent_lof_mechanism_is_not_met(self):
        self.assert_decision(
            self.eligible_null(lof_disease_mechanism="not_established"),
            "NOT_MET",
            "NONE",
        )

    def test_direct_established_lof_and_refuted_validity_are_undetermined(self):
        """Incompatible curated eligibility conclusions must not be silently chosen."""

        self.assert_decision(
            self.eligible_null(
                lof_disease_mechanism="established",
                gene_disease_validity="refuted",
            ),
            "UNDETERMINED",
            "NONE",
            missing=("gene_disease_validity",),
        )

    def test_missing_relevant_transcript_is_undetermined(self):
        evidence = self.eligible_null()
        evidence.pop("biologically_relevant_transcript")
        self.assert_decision(
            evidence,
            "UNDETERMINED",
            "NONE",
            missing=("biologically_relevant_transcript",),
        )

    def test_known_irrelevant_transcript_excludes_pvs1(self):
        self.assert_decision(
            self.eligible_null(biologically_relevant_transcript=False),
            "NOT_MET",
            "NONE",
        )

    def test_missing_nmd_does_not_silently_assign_a_lower_strength(self):
        evidence = self.eligible_null()
        evidence.pop("nmd_prediction")
        self.assert_decision(
            evidence,
            "UNDETERMINED",
            "NONE",
            missing=("nmd_prediction",),
        )

    def test_population_lof_tolerance_only_vetoes_the_role_unknown_percentage_branch(self):
        with self.subTest("predicted NMD remains very strong"):
            self.assert_decision(
                self.eligible_null(population_evidence_for_exon_lof_tolerance=True),
                "ACTIVATED",
                "VERY_STRONG",
            )

        with self.subTest("critical region remains strong when NMD escapes"):
            self.assert_decision(
                self.eligible_null(
                    nmd_prediction=False,
                    critical_protein_region=True,
                    population_evidence_for_exon_lof_tolerance=True,
                ),
                "ACTIVATED",
                "STRONG",
            )

        with self.subTest("role-unknown percentage branch is excluded"):
            self.assert_decision(
                self.eligible_null(
                    nmd_prediction=False,
                    critical_protein_region=False,
                    population_evidence_for_exon_lof_tolerance=True,
                    protein_truncation_percent=15,
                ),
                "NOT_MET",
                "NONE",
            )

    def test_alternative_transcript_rescue_is_an_exclusion(self):
        self.assert_decision(
            self.eligible_null(alternative_transcript_rescue=True),
            "NOT_MET",
            "NONE",
        )

    def test_irrelevant_fields_do_not_change_an_activated_result(self):
        minimal = self.eligible_null()
        enriched = deepcopy(minimal)
        enriched.update(
            {
                "variant_type": "small_variant",
                "coding_position": 42,
                "splice_site_position": 2,
                "start_codon_affected": False,
                "alternative_start_codon": False,
                "upstream_pathogenic_variant_evidence": False,
            }
        )

        baseline = self.assert_decision(minimal, "ACTIVATED", "VERY_STRONG")
        with_extra = self.assert_decision(enriched, "ACTIVATED", "VERY_STRONG")
        self.assertEqual(baseline.status.value, with_extra.status.value)
        self.assertEqual(baseline.strength.value, with_extra.strength.value)
        self.assertEqual(baseline.to_dict(), with_extra.to_dict())

    def test_conflicting_decision_critical_evidence_is_undetermined(self):
        evidence = self.eligible_null(
            nmd_prediction={
                "value": True,
                "state": "conflicting",
                "provenance": [
                    {
                        "source": "synthetic conflicting NMD reports",
                        "kind": "CURATED",
                        "identifier": "TEST-NMD-CONFLICT",
                    }
                ],
            }
        )
        decision = self.assert_decision(
            evidence,
            "UNDETERMINED",
            "NONE",
        )
        self.assertTrue(
            "nmd_prediction" in decision.features_missing
            or "conflict" in decision.to_dict()["reason"].lower()
        )

    def test_cross_field_transcript_and_exon_relevance_contradiction_is_undetermined(self):
        self.assert_decision(
            self.eligible_null(exon_biological_relevance=False),
            "UNDETERMINED",
            "NONE",
            missing=("biologically_relevant_transcript",),
        )

    def test_result_is_deterministic_and_does_not_mutate_input(self):
        evidence = self.eligible_null()
        original = deepcopy(evidence)
        first = self.evaluate(evidence).to_dict()
        second = self.evaluate(evidence).to_dict()
        reordered = dict(reversed(list(evidence.items())))
        reordered_result = self.evaluate(reordered).to_dict()
        self.assertEqual(original, evidence)
        self.assertEqual(first, second)
        self.assertEqual(first, reordered_result)


class SplicePathwayTests(PVS1TestCase):
    @staticmethod
    def observed_rna_evidence():
        return {
            "value": True,
            "state": "known",
            "provenance": [
                {
                    "source": "synthetic allele-specific RNA assay",
                    "kind": "OBSERVED",
                    "identifier": "TEST-RNA-OBSERVED-001",
                }
            ],
        }

    @staticmethod
    def observed_splice_effect(outcome="frameshift_nmd"):
        return {
            "value": outcome,
            "state": "known",
            "provenance": [
                {
                    "source": "synthetic RNA isoform interpretation",
                    "kind": "OBSERVED",
                    "identifier": "TEST-SPLICE-OUTCOME-001",
                }
            ],
        }

    def canonical_splice(self, **overrides):
        evidence = self.eligible_null(
            variant_consequence="canonical_splice_variant",
            splice_effect="frameshift_nmd",
            nmd_prediction=True,
        )
        evidence.update(overrides)
        return evidence

    def test_zero_splice_site_offset_is_rejected(self):
        with self.assertRaisesRegex(ValueError, r"non-zero signed splice-site offset"):
            evaluate_pvs1({
                "variant_consequence": "splice_site",
                "splice_site_position": 0,
            })

    def test_canonical_splice_with_known_frameshift_nmd_is_very_strong(self):
        decision = self.assert_decision(
            self.canonical_splice(),
            "ACTIVATED",
            "VERY_STRONG",
            require_audit=True,
        )
        self.assertNotIn("evidence_mode", decision.to_dict())

    def test_cryptic_splice_prediction_requires_and_uses_an_explicit_frame(self):
        evidence = self.eligible_null(
            variant_consequence="canonical_splice_variant",
            cryptic_splice_site_usage=True,
            reading_frame_disruption=True,
            nmd_prediction=True,
        )
        self.assert_decision(evidence, "ACTIVATED", "VERY_STRONG")

    def test_canonical_splice_with_unknown_outcome_is_undetermined(self):
        evidence = self.eligible_null(variant_consequence="canonical_splice_variant")
        evidence.pop("nmd_prediction")
        self.assert_decision(
            evidence,
            "UNDETERMINED",
            "NONE",
            missing=("splice_effect",),
        )

    def test_canonical_splice_with_known_normal_outcome_is_not_met(self):
        self.assert_decision(
            self.canonical_splice(splice_effect="no_aberrant_splicing"),
            "NOT_MET",
            "NONE",
        )

    def test_in_frame_exon_skipping_with_critical_region_is_strong(self):
        self.assert_decision(
            self.canonical_splice(
                splice_effect="in_frame",
                nmd_prediction=False,
                critical_protein_region=True,
            ),
            "ACTIVATED",
            "STRONG",
        )

    def test_in_frame_exon_skipping_uses_percent_ladder_when_not_critical(self):
        common = {
            "splice_effect": "in_frame",
            "nmd_prediction": False,
            "critical_protein_region": False,
            "population_evidence_for_exon_lof_tolerance": False,
        }
        with self.subTest("greater than ten percent"):
            self.assert_decision(
                self.canonical_splice(**common, protein_truncation_percent=15),
                "ACTIVATED",
                "STRONG",
            )
        with self.subTest("less than ten percent"):
            self.assert_decision(
                self.canonical_splice(**common, protein_truncation_percent=5),
                "ACTIVATED",
                "MODERATE",
            )

    def test_splice_rescue_prevents_early_return_from_nmd_branch(self):
        self.assert_decision(
            self.canonical_splice(splice_rescue_evidence=True),
            "NOT_MET",
            "NONE",
        )

    def test_multiple_credible_splice_outcomes_use_the_lowest_strength(self):
        self.assert_decision(
            self.canonical_splice(
                splice_effect=["frameshift_nmd", "in_frame"],
                nmd_prediction=True,
                critical_protein_region=False,
                population_evidence_for_exon_lof_tolerance=False,
                protein_truncation_percent=5,
            ),
            "ACTIVATED",
            "MODERATE",
        )

    def test_noncanonical_splice_without_rna_evidence_is_not_applicable(self):
        evidence = self.eligible_null(
            variant_consequence="splice_site",
            splice_site_position=5,
            splice_effect="frameshift_nmd",
        )
        self.assert_decision(evidence, "NOT_APPLICABLE", "NONE")

    def test_noncanonical_splice_with_known_no_rna_evidence_is_not_met(self):
        evidence = self.eligible_null(
            variant_consequence="splice_site",
            splice_site_position=5,
            rna_splicing_evidence={
                "value": False,
                "state": "known",
                "provenance": [
                    {
                        "source": "synthetic observed RNA assay",
                        "kind": "OBSERVED",
                        "identifier": "TEST-RNA-NO-ABERRANT-001",
                    }
                ],
            },
            rna_assay_interpretable=True,
        )
        self.assert_decision(evidence, "NOT_MET", "NONE")

    def test_noncanonical_splice_with_observed_interpretable_lof_rna_is_auditable(self):
        evidence = self.eligible_null(
            variant_consequence="splice_site",
            splice_site_position=5,
            rna_splicing_evidence=self.observed_rna_evidence(),
            rna_assay_interpretable=True,
            rna_aberrant_splicing_sufficient_for_lof=True,
            splice_effect=self.observed_splice_effect(),
            nmd_prediction=True,
        )
        decision = self.assert_decision(
            evidence,
            "ACTIVATED",
            "VERY_STRONG",
            require_audit=True,
        )
        payload = decision.to_dict()
        self.assertEqual("RNA", payload["evidence_mode"])
        self.assertTrue(
            any(
                item.get("identifier") == "TEST-RNA-OBSERVED-001"
                for item in payload["evidence_provenance"]
            )
        )
        self.assertEqual("PVS1(RNA) direct-splice guidance", payload["guideline"]["source_locator"])
        self.assertEqual(
            "ClinGen SVI PVS1(RNA) direct-splice evidence with ClinGen SVI 2018 strength mapping",
            payload["guideline_conditions"]["framework"],
        )
        self.assertEqual(
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC6185798/",
            payload["guideline_conditions"]["strength_resolution_source"]["source_url"],
        )

    def test_direct_rna_can_route_unpositioned_or_other_spliceogenic_annotation(self):
        """Observed RNA can establish the splice branch without a predictor label."""

        for label, route in (
            ("unpositioned splice-site annotation", {"variant_consequence": "splice_site"}),
            ("synonymous annotation", {"variant_consequence": "synonymous_variant"}),
        ):
            with self.subTest(label=label):
                decision = self.assert_decision(
                    self.eligible_null(
                        **route,
                        rna_splicing_evidence=self.observed_rna_evidence(),
                        rna_assay_interpretable=True,
                        rna_aberrant_splicing_sufficient_for_lof=True,
                        splice_effect=self.observed_splice_effect(),
                    ),
                    "ACTIVATED",
                    "VERY_STRONG",
                )
                self.assertEqual("RNA", decision.evidence_mode)
                self.assertNotIn("splice_site_position", decision.features_used)

    def test_noninterpretable_rna_does_not_hide_missing_splice_position(self):
        decision = self.assert_decision(
            self.eligible_null(
                variant_consequence="splice_site",
                rna_splicing_evidence=self.observed_rna_evidence(),
                rna_assay_interpretable=False,
                rna_aberrant_splicing_sufficient_for_lof=True,
                splice_effect=self.observed_splice_effect(),
            ),
            "UNDETERMINED",
            "NONE",
            missing=("splice_site_position",),
        )
        self.assertEqual("SVI2018.PATHWAY.SPLICE_POSITION_UNKNOWN", decision.decision_path)

    def test_noncanonical_positive_rna_claim_requires_assay_and_lof_sufficiency(self):
        common = {
            "variant_consequence": "splice_site",
            "splice_site_position": 5,
            "rna_splicing_evidence": self.observed_rna_evidence(),
            "splice_effect": self.observed_splice_effect(),
            "nmd_prediction": True,
        }
        with self.subTest("assay interpretability missing"):
            self.assert_decision(
                self.eligible_null(
                    **common,
                    rna_aberrant_splicing_sufficient_for_lof=True,
                ),
                "UNDETERMINED",
                "NONE",
                missing=("rna_assay_interpretable",),
            )
        with self.subTest("LoF sufficiency missing"):
            self.assert_decision(
                self.eligible_null(
                    **common,
                    rna_assay_interpretable=True,
                ),
                "UNDETERMINED",
                "NONE",
                missing=("rna_aberrant_splicing_sufficient_for_lof",),
            )
        with self.subTest("LoF insufficiency known"):
            self.assert_decision(
                self.eligible_null(
                    **common,
                    rna_assay_interpretable=True,
                    rna_aberrant_splicing_sufficient_for_lof=False,
                ),
                "NOT_MET",
                "NONE",
            )

    def test_predicted_only_rna_provenance_never_activates_rna_mode(self):
        evidence = self.eligible_null(
            variant_consequence="splice_site",
            splice_site_position=5,
            rna_splicing_evidence={
                "value": True,
                "state": "known",
                "provenance": [
                    {
                        "source": "synthetic splice predictor",
                        "kind": "PREDICTED",
                        "identifier": "TEST-RNA-PREDICTED-ONLY",
                    }
                ],
            },
            rna_assay_interpretable=True,
            rna_aberrant_splicing_sufficient_for_lof=True,
            splice_effect=self.observed_splice_effect(),
            nmd_prediction=True,
        )
        decision = self.evaluate(evidence)
        self.assertNotEqual("ACTIVATED", decision.status.value)
        self.assertEqual("NONE", decision.strength.value)
        self.assertNotEqual("RNA", decision.to_dict().get("evidence_mode"))

    def test_usable_canonical_rna_evidence_is_labeled_as_rna_mode(self):
        evidence = self.canonical_splice(
            rna_splicing_evidence=self.observed_rna_evidence(),
            rna_assay_interpretable=True,
            rna_aberrant_splicing_sufficient_for_lof=True,
            splice_effect=self.observed_splice_effect(),
        )
        decision = self.assert_decision(
            evidence,
            "ACTIVATED",
            "VERY_STRONG",
            require_audit=True,
        )
        self.assertEqual("RNA", decision.to_dict()["evidence_mode"])

    def test_incomplete_positive_rna_does_not_gate_canonical_generic_evaluation(self):
        """RNA-only requirements cannot become extra canonical requirements."""

        with self.subTest("assay interpretability missing"):
            decision = self.assert_decision(
                self.canonical_splice(
                    rna_splicing_evidence=self.observed_rna_evidence(),
                    rna_aberrant_splicing_sufficient_for_lof=True,
                ),
                "ACTIVATED",
                "VERY_STRONG",
            )
            self.assertNotEqual("RNA", decision.to_dict().get("evidence_mode"))

        with self.subTest("RNA LoF sufficiency missing"):
            decision = self.assert_decision(
                self.canonical_splice(
                    rna_splicing_evidence=self.observed_rna_evidence(),
                    rna_assay_interpretable=True,
                ),
                "ACTIVATED",
                "VERY_STRONG",
            )
            self.assertNotEqual("RNA", decision.to_dict().get("evidence_mode"))

        with self.subTest("observed RNA product missing"):
            decision = self.assert_decision(
                self.canonical_splice(
                    rna_splicing_evidence=self.observed_rna_evidence(),
                    rna_assay_interpretable=True,
                    rna_aberrant_splicing_sufficient_for_lof=True,
                    splice_effect=None,
                    predicted_exon_skipping=True,
                    reading_frame_disruption=True,
                ),
                "ACTIVATED",
                "VERY_STRONG",
            )
            self.assertNotEqual("RNA", decision.to_dict().get("evidence_mode"))

    def test_conflicting_rna_in_frame_criticality_blocks_a_lower_generic_branch(self):
        evidence = self.eligible_null(
            variant_consequence="splice_site",
            splice_site_position=5,
            rna_splicing_evidence=self.observed_rna_evidence(),
            rna_assay_interpretable=True,
            rna_aberrant_splicing_sufficient_for_lof=True,
            splice_effect=self.observed_splice_effect("in_frame"),
            critical_protein_region=True,
            rna_in_frame_loss_of_clinically_relevant_residues={
                "value": True,
                "state": "conflicting",
                "provenance": [{"source": "synthetic RNA review", "kind": "CURATED"}],
            },
        )
        self.assert_decision(
            evidence,
            "UNDETERMINED",
            "NONE",
            missing=("rna_in_frame_loss_of_clinically_relevant_residues",),
        )

    def test_mixed_rna_lof_and_normal_outcomes_require_residual_function_assessment(self):
        """A normal isoform only changes the result when its function is known."""

        common = {
            "variant_consequence": "splice_site",
            "splice_site_position": 5,
            "rna_splicing_evidence": self.observed_rna_evidence(),
            "rna_assay_interpretable": True,
            "rna_aberrant_splicing_sufficient_for_lof": True,
            "splice_effect": self.observed_splice_effect(
                ["frameshift_nmd", "no_aberrant_splicing"]
            ),
            "nmd_prediction": True,
        }
        with self.subTest("residual normal transcript is sufficient for function"):
            self.assert_decision(
                self.eligible_null(
                    **common,
                    residual_normal_transcript_sufficient_for_function=True,
                ),
                "NOT_MET",
                "NONE",
            )
        with self.subTest("residual normal transcript is insufficient for function"):
            decision = self.assert_decision(
                self.eligible_null(
                    **common,
                    residual_normal_transcript_sufficient_for_function=False,
                ),
                "ACTIVATED",
                "VERY_STRONG",
            )
            self.assertEqual("RNA", decision.to_dict()["evidence_mode"])
        with self.subTest("residual normal transcript function is unresolved"):
            self.assert_decision(
                self.eligible_null(**common),
                "UNDETERMINED",
                "NONE",
                missing=("residual_normal_transcript_sufficient_for_function",),
            )


class StructuralVariantTests(PVS1TestCase):
    def test_whole_gene_deletion_is_very_strong_without_transcript_details(self):
        evidence = {
            "gene_id": "GENE-DEL",
            "disease_id": "DISEASE-DEL",
            "lof_disease_mechanism": "established",
            "gene_disease_validity": "definitive",
            "cnv_event_type": "deletion",
            "whole_gene_deletion": True,
        }
        self.assert_decision(evidence, "ACTIVATED", "VERY_STRONG", require_audit=True)

    def test_direct_whole_gene_deletion_consequence_does_not_audit_an_unknown_flag(self):
        decision = self.assert_decision(
            self.cnv_evidence(
                variant_consequence="whole_gene_deletion",
                lof_disease_mechanism="established",
            ),
            "ACTIVATED",
            "VERY_STRONG",
        )
        self.assertIn("variant_consequence", decision.features_used)
        self.assertNotIn("whole_gene_deletion", decision.features_used)
        self.assertNotIn("whole_gene_deletion", decision.evidence_used)

    def test_exon_level_deletion_with_frameshift_and_nmd_is_very_strong(self):
        self.assert_decision(
            self.cnv_evidence(
                cnv_event_type="deletion",
                whole_gene_deletion=False,
                affected_exons=["3"],
                cnv_reading_frame="frameshift",
                nmd_prediction=True,
            ),
            "ACTIVATED",
            "VERY_STRONG",
        )

    def test_in_frame_exon_deletion_requires_critical_region_for_generic_strong(self):
        self.assert_decision(
            self.cnv_evidence(
                cnv_event_type="deletion",
                whole_gene_deletion=False,
                affected_exons=["3"],
                cnv_reading_frame="in_frame",
                critical_protein_region=True,
            ),
            "ACTIVATED",
            "STRONG",
        )

    def test_noncritical_in_frame_exon_deletion_does_not_borrow_splice_percent_ladder(self):
        self.assert_decision(
            self.cnv_evidence(
                cnv_event_type="deletion",
                whole_gene_deletion=False,
                affected_exons=["3"],
                cnv_reading_frame="in_frame",
                critical_protein_region=False,
                population_evidence_for_exon_lof_tolerance=False,
                protein_truncation_percent=15,
            ),
            "NOT_MET",
            "NONE",
        )

    def test_in_frame_exon_deletion_with_unknown_critical_region_is_undetermined(self):
        self.assert_decision(
            self.cnv_evidence(
                cnv_event_type="deletion",
                whole_gene_deletion=False,
                affected_exons=["3"],
                cnv_reading_frame="in_frame",
            ),
            "UNDETERMINED",
            "NONE",
            missing=("critical_protein_region",),
        )

    def test_cross_field_cnv_frame_contradiction_is_undetermined(self):
        self.assert_decision(
            self.cnv_evidence(
                cnv_event_type="deletion",
                whole_gene_deletion=False,
                affected_exons=["3"],
                cnv_reading_frame="frameshift",
                reading_frame_disruption=False,
            ),
            "UNDETERMINED",
            "NONE",
            missing=("cnv_reading_frame",),
        )

    def test_proven_and_presumed_tandem_duplication_strengths(self):
        common = {
            "cnv_event_type": "duplication",
            "affected_exons": ["3"],
            "whole_gene_duplication": False,
            "duplication_contained_within_gene": True,
            "cnv_reading_frame": "frameshift",
            "nmd_prediction": True,
        }
        with self.subTest("proven tandem"):
            self.assert_decision(
                self.cnv_evidence(**common, tandem_duplication_status="proven_tandem"),
                "ACTIVATED",
                "VERY_STRONG",
            )
        with self.subTest("presumed tandem"):
            self.assert_decision(
                self.cnv_evidence(**common, tandem_duplication_status="presumed_tandem"),
                "ACTIVATED",
                "STRONG",
            )

    def test_duplication_with_unknown_decision_critical_structure_is_undetermined(self):
        evidence = self.cnv_evidence(
            cnv_event_type="duplication",
            affected_exons=["3"],
            whole_gene_duplication=False,
            duplication_contained_within_gene=True,
            tandem_duplication_status="proven_tandem",
        )
        self.assert_decision(
            evidence,
            "UNDETERMINED",
            "NONE",
            missing=("cnv_reading_frame",),
        )

    def test_whole_gene_duplication_is_not_generic_intragenic_pvs1(self):
        self.assert_decision(
            self.cnv_evidence(
                cnv_event_type="duplication",
                affected_exons=["3"],
                duplication_contained_within_gene=True,
                whole_gene_duplication=True,
                tandem_duplication_status="proven_tandem",
                cnv_reading_frame="frameshift",
                nmd_prediction=True,
            ),
            "NOT_MET",
            "NONE",
        )

    def test_duplication_reports_the_first_missing_decision_critical_fact(self):
        common = {
            "cnv_event_type": "duplication",
            "affected_exons": ["3"],
            "cnv_reading_frame": "frameshift",
            "nmd_prediction": True,
        }
        with self.subTest("whole-gene duplication status unknown"):
            self.assert_decision(
                self.cnv_evidence(
                    **common,
                    duplication_contained_within_gene=True,
                    tandem_duplication_status="proven_tandem",
                ),
                "UNDETERMINED",
                "NONE",
                missing=("whole_gene_duplication",),
            )
        with self.subTest("intragenic status unknown"):
            self.assert_decision(
                self.cnv_evidence(
                    **common,
                    whole_gene_duplication=False,
                    tandem_duplication_status="proven_tandem",
                ),
                "UNDETERMINED",
                "NONE",
                missing=("duplication_contained_within_gene",),
            )
        with self.subTest("tandem status unknown"):
            self.assert_decision(
                self.cnv_evidence(
                    **common,
                    whole_gene_duplication=False,
                    duplication_contained_within_gene=True,
                ),
                "UNDETERMINED",
                "NONE",
                missing=("tandem_duplication_status",),
            )
        with self.subTest("NMD unknown"):
            incomplete = dict(common)
            incomplete.pop("nmd_prediction")
            self.assert_decision(
                self.cnv_evidence(
                    **incomplete,
                    whole_gene_duplication=False,
                    duplication_contained_within_gene=True,
                    tandem_duplication_status="proven_tandem",
                ),
                "UNDETERMINED",
                "NONE",
                missing=("nmd_prediction",),
            )

    def test_known_non_disruptive_duplication_is_not_met(self):
        self.assert_decision(
            self.cnv_evidence(
                cnv_event_type="duplication",
                affected_exons=["3"],
                whole_gene_duplication=False,
                duplication_contained_within_gene=True,
                tandem_duplication_status="proven_tandem",
                cnv_reading_frame="in_frame",
                nmd_prediction=False,
            ),
            "NOT_MET",
            "NONE",
        )


class InitiationCodonAndOverrideTests(PVS1TestCase):
    def test_initiation_codon_pathway(self):
        with self.subTest("alternative functional start"):
            self.assert_decision(
                self.start_loss_evidence(alternative_start_codon=True),
                "NOT_MET",
                "NONE",
            )
        with self.subTest("pathogenic evidence upstream of next start"):
            self.assert_decision(
                self.start_loss_evidence(upstream_pathogenic_variant_evidence=True),
                "ACTIVATED",
                "MODERATE",
            )
        with self.subTest("no upstream pathogenic evidence"):
            self.assert_decision(
                self.start_loss_evidence(),
                "ACTIVATED",
                "SUPPORTING",
            )

    def test_functional_alternative_start_excludes_without_relevance_input(self):
        evidence = self.start_loss_evidence(alternative_start_codon=True)
        evidence.pop("biologically_relevant_transcript")
        self.assert_decision(evidence, "NOT_MET", "NONE")

    def test_initiation_codon_missing_alternative_start_evidence_is_undetermined(self):
        evidence = self.start_loss_evidence()
        evidence.pop("alternative_start_codon")
        self.assert_decision(
            evidence,
            "UNDETERMINED",
            "NONE",
            missing=("alternative_start_codon",),
        )

    def test_initiation_codon_missing_upstream_pathogenic_evidence_is_undetermined(self):
        evidence = self.start_loss_evidence()
        evidence.pop("upstream_pathogenic_variant_evidence")
        self.assert_decision(
            evidence,
            "UNDETERMINED",
            "NONE",
            missing=("upstream_pathogenic_variant_evidence",),
        )

    @staticmethod
    def synthetic_start_loss_override():
        return {
            "name": "Synthetic start-loss VCEP specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/synthetic-start-loss-v1",
            "specification_id": "SYNTHETIC-START-LOSS-1",
            "gene_id": "GENE-VCEP",
            "disease_id": "DISEASE-VCEP",
            "transcript_id": "TX-VCEP.1",
            "rules": [
                {
                    "source_locator": "synthetic-start-loss-override",
                    "pathway": "INITIATION_CODON",
                    "when": {
                        "alternative_start_codon": False,
                        "upstream_pathogenic_variant_evidence": False,
                    },
                    "decision": {
                        "status": "ACTIVATED",
                        "strength": "VERY_STRONG",
                        "decision_path": "SYNTHETIC_VCEP_START_LOSS_VS",
                    },
                }
            ],
        }

    @staticmethod
    def synthetic_unclassified_splice_override():
        """A VCEP rule that can be selected before generic splice routing."""

        return {
            "name": "Synthetic unclassified-splice VCEP specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/synthetic-splice-v1",
            "specification_id": "SYNTHETIC-SPLICE-1",
            "gene_id": "GENE-SPLICE-VCEP",
            "disease_id": "DISEASE-SPLICE-VCEP",
            "transcript_id": "TX-SPLICE-VCEP.1",
            "rules": [
                {
                    # ANY permits a scoped rule to take precedence before the
                    # generic classifier can establish canonical +/-1/2
                    # versus a non-canonical splice-site position.
                    "source_locator": "synthetic-unclassified-splice-override",
                    "pathway": "ANY",
                    "when": {
                        "variant_consequence": "splice_site",
                        "splice_effect": "out_of_frame",
                    },
                    "decision": {
                        "status": "ACTIVATED",
                        "strength": "VERY_STRONG",
                        "decision_path": "SYNTHETIC_VCEP_UNCLASSIFIED_SPLICE_VS",
                    },
                }
            ],
        }

    @staticmethod
    def pten_gn003_v32_start_loss_override():
        """Pinned external VCEP fixture; it is supplied by the caller, not bundled."""

        return {
            "name": "ClinGen PTEN VCEP GN003",
            "version": "3.2",
            "source_url": "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN003",
            "specification_id": "GN003",
            "gene_id": "PTEN",
            "disease_id": "MONDO:0017623",
            "transcript_id": "NM_000314.8",
            "mode_of_inheritance": "AD",
            "rules": [
                {
                    "source_locator": "GN003 v3.2 initiation-codon PVS1 rule",
                    "pathway": "INITIATION_CODON",
                    "when": {
                        "alternative_start_codon": False,
                        "upstream_pathogenic_variant_evidence": False,
                    },
                    "decision": {
                        "status": "ACTIVATED",
                        "strength": "VERY_STRONG",
                        "decision_path": "GN003_V3_2.INITIATION_CODON.PVS1",
                    },
                }
            ],
        }

    def test_matching_versioned_vcep_rule_overrides_generic_start_loss_rule(self):
        evidence = self.start_loss_evidence(
            gene_id="GENE-VCEP",
            disease_id="DISEASE-VCEP",
            transcript_id="TX-VCEP.1",
            gene_specific_pvs1_rules=[self.synthetic_start_loss_override()],
        )
        decision = self.assert_decision(
            evidence,
            "ACTIVATED",
            "VERY_STRONG",
            require_audit=True,
        )
        payload = decision.to_dict()
        self.assertEqual("Synthetic start-loss VCEP specification", payload["guideline"]["name"])
        self.assertEqual("1.0.0", payload["guideline"]["version"])
        self.assertEqual("SYNTHETIC-START-LOSS-1", payload["guideline"]["rule_id"])
        self.assertEqual("SYNTHETIC_VCEP_START_LOSS_VS", payload["decision_path"])

    def test_matching_vcep_rule_precedes_generic_lof_noneligibility(self):
        """A scoped higher-authority VCEP rule is evaluated before generic LoF gates."""

        evidence = self.start_loss_evidence(
            gene_id="GENE-VCEP",
            disease_id="DISEASE-VCEP",
            transcript_id="TX-VCEP.1",
            lof_disease_mechanism="not_established",
            gene_specific_pvs1_rules=[self.synthetic_start_loss_override()],
        )
        self.assert_decision(
            evidence,
            "ACTIVATED",
            "VERY_STRONG",
            require_audit=True,
        )

    def test_path_specific_start_loss_vcep_rule_cannot_hijack_frameshift_pathway(self):
        """Only an applicable pathway (or ANY) may take early VCEP precedence."""

        evidence = self.eligible_null(
            gene_id="GENE-VCEP",
            disease_id="DISEASE-VCEP",
            transcript_id="TX-VCEP.1",
            alternative_start_codon=False,
            upstream_pathogenic_variant_evidence=False,
            gene_specific_pvs1_rules=[self.synthetic_start_loss_override()],
        )
        decision = self.assert_decision(
            evidence,
            "ACTIVATED",
            "VERY_STRONG",
            require_audit=True,
        )
        self.assertNotEqual(
            "SYNTHETIC_VCEP_START_LOSS_VS",
            decision.to_dict()["decision_path"],
        )

    def test_vcep_any_rule_precedes_generic_splice_position_classification(self):
        """A scoped override must not be blocked by absent generic +/-1/2 data."""

        evidence = self.eligible_null(
            gene_id="GENE-SPLICE-VCEP",
            disease_id="DISEASE-SPLICE-VCEP",
            transcript_id="TX-SPLICE-VCEP.1",
            variant_consequence="splice_site",
            splice_effect="frameshift_nmd",
            gene_specific_pvs1_rules=[self.synthetic_unclassified_splice_override()],
        )
        self.assertNotIn("splice_site_position", evidence)
        decision = self.assert_decision(
            evidence,
            "ACTIVATED",
            "VERY_STRONG",
            require_audit=True,
        )
        self.assertEqual(
            "SYNTHETIC_VCEP_UNCLASSIFIED_SPLICE_VS",
            decision.to_dict()["decision_path"],
        )

    def test_scoped_pten_gn003_v32_start_loss_fixture_overrides_generic_strength(self):
        """A realistic, version-pinned caller-supplied VCEP fixture is auditable."""

        evidence = self.start_loss_evidence(
            gene_id="PTEN",
            disease_id="MONDO:0017623",
            transcript_id="NM_000314.8",
            mode_of_inheritance="AD",
            gene_specific_pvs1_rules=[self.pten_gn003_v32_start_loss_override()],
        )
        decision = self.assert_decision(
            evidence,
            "ACTIVATED",
            "VERY_STRONG",
            require_audit=True,
        )
        guideline = decision.to_dict()["guideline"]
        self.assertEqual("GN003", guideline["rule_id"])
        self.assertEqual("3.2", guideline["version"])
        self.assertEqual(
            "https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN003",
            guideline["source_url"],
        )

    def test_vcep_manifest_range_operators_match_and_malformed_ones_are_rejected(self):
        valid = self.synthetic_start_loss_override()
        valid["rules"][0]["when"]["protein_truncation_percent"] = {
            "between": [5, 10]
        }
        matching_evidence = self.start_loss_evidence(
            gene_id="GENE-VCEP",
            disease_id="DISEASE-VCEP",
            transcript_id="TX-VCEP.1",
            protein_truncation_percent=5,
            gene_specific_pvs1_rules=[valid],
        )
        self.assert_decision(matching_evidence, "ACTIVATED", "VERY_STRONG")

        with self.subTest("reversed inclusive range"):
            malformed = self.synthetic_start_loss_override()
            malformed["rules"][0]["when"]["protein_truncation_percent"] = {
                "between": [10, 5]
            }
            evidence = self.start_loss_evidence(
                gene_id="GENE-VCEP",
                disease_id="DISEASE-VCEP",
                transcript_id="TX-VCEP.1",
                protein_truncation_percent=5,
                gene_specific_pvs1_rules=[malformed],
            )
            with self.assertRaises(ValueError):
                evaluate_pvs1(evidence)

        with self.subTest("unsupported operator"):
            malformed = self.synthetic_start_loss_override()
            malformed["rules"][0]["when"]["protein_truncation_percent"] = {
                "approximately": 5
            }
            evidence = self.start_loss_evidence(
                gene_id="GENE-VCEP",
                disease_id="DISEASE-VCEP",
                transcript_id="TX-VCEP.1",
                protein_truncation_percent=5,
                gene_specific_pvs1_rules=[malformed],
            )
            with self.assertRaises(ValueError):
                evaluate_pvs1(evidence)

    def test_wrong_scope_does_not_apply_a_gene_specific_override(self):
        with self.subTest("disease mismatch"):
            evidence = self.start_loss_evidence(
                gene_id="GENE-VCEP",
                disease_id="OTHER-DISEASE",
                transcript_id="TX-VCEP.1",
                gene_specific_pvs1_rules=[self.synthetic_start_loss_override()],
            )
            self.assert_decision(evidence, "ACTIVATED", "SUPPORTING")

        with self.subTest("transcript mismatch"):
            evidence = self.start_loss_evidence(
                gene_id="GENE-VCEP",
                disease_id="DISEASE-VCEP",
                transcript_id="TX-OTHER.1",
                gene_specific_pvs1_rules=[self.synthetic_start_loss_override()],
            )
            self.assert_decision(evidence, "ACTIVATED", "SUPPORTING")

    def test_unresolved_potential_vcep_scope_is_not_bypassed_by_generic_rule(self):
        evidence = self.start_loss_evidence(
            gene_id="GENE-VCEP",
            transcript_id="TX-VCEP.1",
            gene_specific_pvs1_rules=[self.synthetic_start_loss_override()],
        )
        evidence.pop("disease_id")
        decision = self.assert_decision(
            evidence,
            "UNDETERMINED",
            "NONE",
            missing=("disease_id",),
        )
        payload = decision.to_dict()
        self.assertEqual(
            "caller-supplied VCEP specification selection",
            payload["guideline_conditions"]["framework"],
        )
        self.assertEqual("UNRESOLVED", payload["guideline_conditions"]["selection"])
        self.assertEqual("UNKNOWN", payload["evidence_used"]["disease_id"]["state"])

    def test_vcep_specification_requires_versioned_metadata(self):
        unversioned = self.synthetic_start_loss_override()
        unversioned.pop("version")
        evidence = self.start_loss_evidence(
            gene_id="GENE-VCEP",
            disease_id="DISEASE-VCEP",
            transcript_id="TX-VCEP.1",
            gene_specific_pvs1_rules=[unversioned],
        )
        with self.assertRaises(ValueError):
            evaluate_pvs1(evidence)

    def test_conflicting_equally_scoped_vcep_rules_are_undetermined(self):
        first = self.synthetic_start_loss_override()
        second = deepcopy(first)
        second["name"] = "Conflicting synthetic specification"
        second["version"] = "2.0.0"
        second["rules"][0]["decision"] = {
            "status": "ACTIVATED",
            "strength": "MODERATE",
            "decision_path": "SYNTHETIC_CONFLICTING_START_LOSS_MODERATE",
        }
        evidence = self.start_loss_evidence(
            gene_id="GENE-VCEP",
            disease_id="DISEASE-VCEP",
            transcript_id="TX-VCEP.1",
            gene_specific_pvs1_rules=[first, second],
        )
        self.assert_decision(evidence, "UNDETERMINED", "NONE")


class OutputContractTests(PVS1TestCase):
    def test_unsupported_consequence_is_not_applicable_not_not_met(self):
        self.assert_decision(
            self.eligible_null(variant_consequence="missense_variant"),
            "NOT_APPLICABLE",
            "NONE",
        )

    def test_evidence_provenance_is_preserved_in_auditable_output(self):
        evidence = self.eligible_null(
            nmd_prediction={
                "value": True,
                "state": "known",
                "provenance": [
                    {
                        "source": "synthetic RNA assay",
                        "kind": "OBSERVED",
                        "identifier": "TEST-RNA-001",
                    }
                ],
            }
        )
        decision = self.assert_decision(
            evidence,
            "ACTIVATED",
            "VERY_STRONG",
            require_audit=True,
        )
        provenance = decision.to_dict()["evidence_provenance"]
        self.assertTrue(any(item.get("identifier") == "TEST-RNA-001" for item in provenance))

    def test_generic_and_vcep_results_include_condition_snapshots(self):
        generic = self.assert_decision(
            self.eligible_null(), "ACTIVATED", "VERY_STRONG", require_audit=True
        ).to_dict()
        self.assertEqual("FRAMESHIFT", generic["evidence_used"]["variant_consequence"]["value"])
        self.assertTrue(generic["guideline_conditions"])
        self.assertEqual(
            True,
            generic["guideline_conditions"]["evaluated_conditions"]["nmd_prediction"]["value"],
        )

        spec = {
            "name": "Synthetic auditable start-loss VCEP specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/audit-start-loss-v1",
            "specification_id": "SYNTHETIC-AUDIT-1",
            "gene_id": "GENE-AUDIT",
            "disease_id": "DISEASE-AUDIT",
            "rules": [{
                "source_locator": "synthetic audit start-loss rule",
                "pathway": "INITIATION_CODON",
                "when": {
                    "variant_consequence": "START_LOST",
                    "alternative_start_codon": False,
                    "upstream_pathogenic_variant_evidence": False,
                },
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "VERY_STRONG",
                    "decision_path": "SYNTHETIC_AUDIT_START_LOSS_VS",
                },
            }],
        }
        vcep = self.assert_decision(
            self.start_loss_evidence(
                gene_id="GENE-AUDIT",
                disease_id="DISEASE-AUDIT",
                gene_specific_pvs1_rules=[spec],
            ),
            "ACTIVATED",
            "VERY_STRONG",
        ).to_dict()
        self.assertEqual("INITIATION_CODON", vcep["guideline_conditions"]["pathway"])
        self.assertEqual("GENE-AUDIT", vcep["guideline_conditions"]["scope"]["gene_id"])
        self.assertIn("variant_consequence", vcep["evidence_used"])


class SafetyRegressionTests(PVS1TestCase):
    @staticmethod
    def deletion_override():
        return {
            "name": "Synthetic deletion VCEP specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/deletion-v1",
            "gene_id": "GENE-STRUCTURE",
            "disease_id": "DISEASE-STRUCTURE",
            "rules": [{
                "source_locator": "synthetic exon deletion rule",
                "pathway": "EXON_DELETION",
                "when": {"cnv_event_type": "DELETION"},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "VERY_STRONG",
                    "decision_path": "SYNTHETIC_DELETION_VS",
                },
            }],
        }

    def test_material_structural_conflict_cannot_be_bypassed_by_vcep(self):
        evidence = self.cnv_evidence(
            gene_id="GENE-STRUCTURE",
            disease_id="DISEASE-STRUCTURE",
            cnv_event_type="deletion",
            whole_gene_deletion=False,
            whole_gene_duplication=True,
            gene_specific_pvs1_rules=[self.deletion_override()],
        )
        self.assert_decision(evidence, "UNDETERMINED", "NONE")

    def test_malformed_vcep_manifest_is_rejected_before_structural_short_circuit(self):
        """Schema safety cannot depend on which clinical branch returns first."""

        malformed = {
            "name": "Synthetic malformed specification",
            # Deliberately omit the mandatory pinned version.
            "source_url": "https://example.invalid/pvs1/malformed-v1",
            "gene_id": "GENE-STRUCTURE",
            "disease_id": "DISEASE-STRUCTURE",
            "rules": [{
                "source_locator": "synthetic null rule",
                "pathway": "NONSENSE_OR_FRAMESHIFT",
                "when": {"variant_consequence": "FRAMESHIFT"},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "VERY_STRONG",
                    "decision_path": "SYNTHETIC_NULL_VS",
                },
            }],
        }
        evidence = self.eligible_null(
            cnv_event_type="deletion",
            gene_specific_pvs1_rules=[malformed],
        )
        with self.assertRaisesRegex(ValueError, r"version"):
            evaluate_pvs1(evidence)

    def test_conflicting_full_gene_evidence_cannot_be_ignored_on_a_small_variant_path(self):
        spec = {
            "name": "Synthetic null-variant VCEP specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/null-v1",
            "gene_id": "GENE-TEST",
            "disease_id": "DISEASE-TEST",
            "rules": [{
                "source_locator": "synthetic null rule",
                "pathway": "NONSENSE_OR_FRAMESHIFT",
                "when": {"variant_consequence": "FRAMESHIFT"},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "VERY_STRONG",
                    "decision_path": "SYNTHETIC_NULL_VS",
                },
            }],
        }
        for field in ("whole_gene_deletion", "whole_gene_duplication"):
            with self.subTest(field=field):
                decision = self.assert_decision(
                    self.eligible_null(
                        **{field: {"state": "conflicting", "value": True}},
                        gene_specific_pvs1_rules=[spec],
                    ),
                    "UNDETERMINED",
                    "NONE",
                )
                self.assertEqual("SVI2018.PATHWAY.STRUCTURAL_CONFLICT", decision.decision_path)

    def test_conflicting_full_gene_evidence_cannot_bypass_unannotated_start_loss(self):
        spec = {
            "name": "Synthetic initiation VCEP specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/initiation-v1",
            "gene_id": "GENE-START-CONFLICT",
            "disease_id": "DISEASE-START-CONFLICT",
            "rules": [{
                "source_locator": "synthetic initiation rule",
                "pathway": "INITIATION_CODON",
                "when": {"start_codon_affected": True},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "VERY_STRONG",
                    "decision_path": "SYNTHETIC_INITIATION_VS",
                },
            }],
        }
        for field in ("whole_gene_deletion", "whole_gene_duplication"):
            with self.subTest(field=field):
                decision = self.assert_decision(
                    {
                        "gene_id": "GENE-START-CONFLICT",
                        "disease_id": "DISEASE-START-CONFLICT",
                        "lof_disease_mechanism": "established",
                        "biologically_relevant_transcript": True,
                        "start_codon_affected": True,
                        "alternative_start_codon": False,
                        "upstream_pathogenic_variant_evidence": False,
                        field: {"state": "conflicting", "value": True},
                        "gene_specific_pvs1_rules": [spec],
                    },
                    "UNDETERMINED",
                    "NONE",
                )
                self.assertEqual("SVI2018.PATHWAY.STRUCTURAL_CONFLICT", decision.decision_path)

    def test_conflicting_whole_gene_duplication_cannot_be_bypassed_by_vcep(self):
        spec = {
            "name": "Synthetic duplication VCEP specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/duplication-v1",
            "gene_id": "GENE-DUP-CONFLICT",
            "disease_id": "DISEASE-DUP-CONFLICT",
            "rules": [{
                "source_locator": "synthetic duplication rule",
                "pathway": "DUPLICATION",
                "when": {"variant_consequence": "DUPLICATION"},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "VERY_STRONG",
                    "decision_path": "SYNTHETIC_DUPLICATION_VS",
                },
            }],
        }
        self.assert_decision(
            self.cnv_evidence(
                gene_id="GENE-DUP-CONFLICT",
                disease_id="DISEASE-DUP-CONFLICT",
                variant_consequence="duplication",
                whole_gene_duplication={"state": "conflicting", "value": True},
                gene_specific_pvs1_rules=[spec],
            ),
            "UNDETERMINED",
            "NONE",
            missing=("whole_gene_duplication",),
        )

    def test_full_gene_event_conflicts_cannot_be_bypassed_by_vcep(self):
        whole_deletion_spec = {
            "name": "Synthetic whole-gene deletion VCEP specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/full-deletion-v1",
            "gene_id": "GENE-FULL-CONFLICT",
            "disease_id": "DISEASE-FULL-CONFLICT",
            "rules": [{
                "source_locator": "synthetic full deletion rule",
                "pathway": "WHOLE_GENE_DELETION",
                "when": {"whole_gene_deletion": True},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "VERY_STRONG",
                    "decision_path": "SYNTHETIC_FULL_DELETION_VS",
                },
            }],
        }
        whole_duplication_spec = {
            "name": "Synthetic whole-gene duplication VCEP specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/full-duplication-v1",
            "gene_id": "GENE-FULL-CONFLICT",
            "disease_id": "DISEASE-FULL-CONFLICT",
            "rules": [{
                "source_locator": "synthetic full duplication rule",
                "pathway": "DUPLICATION",
                "when": {"whole_gene_duplication": True},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "VERY_STRONG",
                    "decision_path": "SYNTHETIC_FULL_DUPLICATION_VS",
                },
            }],
        }
        with self.subTest("deletion versus conflicting duplication evidence"):
            self.assert_decision(
                {
                    "gene_id": "GENE-FULL-CONFLICT",
                    "disease_id": "DISEASE-FULL-CONFLICT",
                    "whole_gene_deletion": True,
                    "whole_gene_duplication": {"state": "conflicting", "value": True},
                    "gene_specific_pvs1_rules": [whole_deletion_spec],
                },
                "UNDETERMINED",
                "NONE",
            )
        with self.subTest("duplication versus conflicting deletion evidence"):
            self.assert_decision(
                {
                    "gene_id": "GENE-FULL-CONFLICT",
                    "disease_id": "DISEASE-FULL-CONFLICT",
                    "whole_gene_duplication": True,
                    "whole_gene_deletion": {"state": "conflicting", "value": True},
                    "gene_specific_pvs1_rules": [whole_duplication_spec],
                },
                "UNDETERMINED",
                "NONE",
            )

    def test_exon_and_whole_gene_deletion_claims_are_not_silently_merged(self):
        for event in (None, "deletion"):
            with self.subTest(cnv_event_type=event):
                evidence = self.cnv_evidence(
                    variant_consequence="exon_deletion",
                    whole_gene_deletion=True,
                )
                if event is not None:
                    evidence["cnv_event_type"] = event
                self.assert_decision(evidence, "UNDETERMINED", "NONE")

    def test_start_codon_cross_fields_route_or_conflict_safely(self):
        self.assert_decision(
            self.start_loss_evidence(start_codon_affected=False),
            "UNDETERMINED",
            "NONE",
        )
        met1_stop = self.eligible_null(
            variant_consequence="stop_gain",
            start_codon_affected=True,
            alternative_start_codon=False,
            upstream_pathogenic_variant_evidence=False,
        )
        decision = self.assert_decision(met1_stop, "ACTIVATED", "SUPPORTING")
        self.assertIn("INITIATION_CODON", decision.decision_path)

    def test_false_rna_summary_and_observed_damaging_profile_are_conflicting(self):
        evidence = self.eligible_null(
            variant_consequence="splice_site",
            splice_site_position=5,
            rna_splicing_evidence={
                "state": "KNOWN",
                "value": False,
                "provenance": [{"source": "synthetic RNA", "kind": "OBSERVED"}],
            },
            rna_assay_interpretable=True,
            splice_effect={
                "state": "KNOWN",
                "value": "OUT_OF_FRAME",
                "provenance": [{"source": "synthetic RNA product", "kind": "OBSERVED"}],
            },
        )
        self.assert_decision(
            evidence,
            "UNDETERMINED",
            "NONE",
            missing=("rna_splicing_evidence", "splice_effect"),
        )

    def test_duplication_known_exclusions_do_not_require_activation_facts(self):
        with self.subTest("whole_gene"):
            self.assert_decision(
                self.cnv_evidence(
                    cnv_event_type="duplication",
                    whole_gene_duplication=True,
                ),
                "NOT_MET",
                "NONE",
            )
        with self.subTest("proven_not_tandem"):
            self.assert_decision(
                self.cnv_evidence(
                    cnv_event_type="duplication",
                    tandem_duplication_status="proven_not_tandem",
                ),
                "NOT_MET",
                "NONE",
            )

    def test_early_vcep_any_requires_a_variant_class_anchor(self):
        spec = {
            "name": "Synthetic unanchored ANY specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/unanchored-any-v1",
            "gene_id": "GENE-ANY",
            "disease_id": "DISEASE-ANY",
            "rules": [{
                "source_locator": "synthetic unanchored rule",
                "pathway": "ANY",
                "when": {"alternative_start_codon": False},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "VERY_STRONG",
                    "decision_path": "SYNTHETIC_UNANCHORED_ANY_VS",
                },
            }],
        }
        evidence = {
            "gene_id": "GENE-ANY",
            "disease_id": "DISEASE-ANY",
            "lof_disease_mechanism": "established",
            "alternative_start_codon": False,
            "gene_specific_pvs1_rules": [spec],
        }
        self.assert_decision(
            evidence,
            "UNDETERMINED",
            "NONE",
            missing=("variant_consequence",),
        )

    def test_early_vcep_anchors_must_match_the_actual_supported_pathway(self):
        mixed_initiation = {
            "name": "Synthetic mixed initiation specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/mixed-initiation-v1",
            "gene_id": "GENE-EARLY-ANCHOR",
            "disease_id": "DISEASE-EARLY-ANCHOR",
            "rules": [{
                "source_locator": "synthetic mixed initiation condition",
                "pathway": "INITIATION_CODON",
                "when": {"variant_consequence": {"in": ["START_LOST", "OTHER"]}},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "VERY_STRONG",
                    "decision_path": "SYNTHETIC_MIXED_INITIATION_VS",
                },
            }],
        }
        self.assert_decision(
            self.eligible_null(
                gene_id="GENE-EARLY-ANCHOR",
                disease_id="DISEASE-EARLY-ANCHOR",
                variant_consequence="missense_variant",
                gene_specific_pvs1_rules=[mixed_initiation],
            ),
            "NOT_APPLICABLE",
            "NONE",
        )

        unpositioned_noncanonical = {
            "name": "Synthetic noncanonical specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/noncanonical-v1",
            "gene_id": "GENE-EARLY-ANCHOR",
            "disease_id": "DISEASE-EARLY-ANCHOR",
            "rules": [{
                "source_locator": "synthetic noncanonical condition",
                "pathway": "NONCANONICAL_SPLICE",
                "when": {"variant_consequence": "SPLICE_SITE"},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "VERY_STRONG",
                    "decision_path": "SYNTHETIC_NONCANONICAL_VS",
                },
            }],
        }
        self.assert_decision(
            self.eligible_null(
                gene_id="GENE-EARLY-ANCHOR",
                disease_id="DISEASE-EARLY-ANCHOR",
                variant_consequence="splice_site",
                gene_specific_pvs1_rules=[unpositioned_noncanonical],
            ),
            "UNDETERMINED",
            "NONE",
            missing=("splice_site_position",),
        )

    def test_position_and_exon_conditioned_vcep_rules_can_anchor_missing_generic_class(self):
        cases = (
            (
                "canonical splice position",
                "CANONICAL_SPLICE",
                {"splice_site_position": 1},
                {"splice_site_position": 1},
                "SYNTHETIC_POSITIONAL_CANONICAL",
            ),
            (
                "affected exon interval",
                "EXON_DELETION",
                {"affected_exons": {"contains_any": ["3"]}},
                {"affected_exons": ["3"]},
                "SYNTHETIC_EXON_INTERVAL",
            ),
        )
        for label, pathway, when, evidence, decision_path in cases:
            with self.subTest(label=label):
                spec = {
                    "name": f"Synthetic {label} VCEP specification",
                    "version": "1.0.0",
                    "source_url": "https://example.invalid/pvs1/positional-v1",
                    "gene_id": "GENE-POSITIONAL",
                    "disease_id": "DISEASE-POSITIONAL",
                    "rules": [{
                        "source_locator": f"synthetic {label} rule",
                        "pathway": pathway,
                        "when": when,
                        "decision": {
                            "status": "ACTIVATED",
                            "strength": "STRONG",
                            "decision_path": decision_path,
                        },
                    }],
                }
                decision = self.assert_decision(
                    {
                        "gene_id": "GENE-POSITIONAL",
                        "disease_id": "DISEASE-POSITIONAL",
                        **evidence,
                        "gene_specific_pvs1_rules": [spec],
                    },
                    "ACTIVATED",
                    "STRONG",
                )
                self.assertEqual(decision_path, decision.decision_path)

    def test_full_gene_conflict_blocks_positional_and_exon_vcep_anchors(self):
        anchors = (
            (
                "splice position",
                "CANONICAL_SPLICE",
                {"splice_site_position": 1},
                {"splice_site_position": 1},
            ),
            (
                "affected exon",
                "EXON_DELETION",
                {"affected_exons": {"contains_any": ["3"]}},
                {"affected_exons": ["3"]},
            ),
        )
        for label, pathway, when, anchor_evidence in anchors:
            for conflict_field in ("whole_gene_deletion", "whole_gene_duplication"):
                with self.subTest(anchor=label, conflict_field=conflict_field):
                    spec = {
                        "name": "Synthetic structural-safety VCEP specification",
                        "version": "1.0.0",
                        "source_url": "https://example.invalid/pvs1/structural-safety-v1",
                        "gene_id": "GENE-POSITIONAL-CONFLICT",
                        "disease_id": "DISEASE-POSITIONAL-CONFLICT",
                        "rules": [{
                            "source_locator": "synthetic structural-safety rule",
                            "pathway": pathway,
                            "when": when,
                            "decision": {
                                "status": "ACTIVATED",
                                "strength": "VERY_STRONG",
                                "decision_path": "SYNTHETIC_STRUCTURAL_SAFETY_VS",
                            },
                        }],
                    }
                    decision = self.assert_decision(
                        {
                            "gene_id": "GENE-POSITIONAL-CONFLICT",
                            "disease_id": "DISEASE-POSITIONAL-CONFLICT",
                            **anchor_evidence,
                            conflict_field: {"state": "conflicting", "value": True},
                            "gene_specific_pvs1_rules": [spec],
                        },
                        "UNDETERMINED",
                        "NONE",
                    )
                    self.assertEqual("SVI2018.PATHWAY.STRUCTURAL_CONFLICT", decision.decision_path)

    def test_vcep_alias_conditions_are_normalized_before_matching(self):
        spec = {
            "name": "Synthetic alias condition specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/alias-condition-v1",
            "gene_id": "GENE-ALIASES",
            "disease_id": "DISEASE-ALIASES",
            "rules": [{
                "source_locator": "synthetic consequence aliases",
                "pathway": "NONSENSE_OR_FRAMESHIFT",
                "when": {
                    "variant_consequence": {
                        "in": ["nonsense", "frameshift_variant"]
                    },
                },
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "STRONG",
                    "decision_path": "SYNTHETIC_ALIAS_CONDITION_STRONG",
                },
            }],
        }
        decision = self.assert_decision(
            self.eligible_null(
                gene_id="GENE-ALIASES",
                disease_id="DISEASE-ALIASES",
                gene_specific_pvs1_rules=[spec],
            ),
            "ACTIVATED",
            "STRONG",
        )
        self.assertEqual("SYNTHETIC_ALIAS_CONDITION_STRONG", decision.decision_path)

    def test_unordered_vcep_condition_sets_have_canonical_audit_order(self):
        spec = {
            "name": "Synthetic unordered condition specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/unordered-condition-v1",
            "gene_id": "GENE-UNORDERED",
            "disease_id": "DISEASE-UNORDERED",
            "rules": [{
                "source_locator": "synthetic unordered condition",
                "pathway": "NONSENSE_OR_FRAMESHIFT",
                "when": {
                    "variant_consequence": {"frameshift_variant", "stop_gained"},
                },
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "STRONG",
                    "decision_path": "SYNTHETIC_UNORDERED_CONDITION_STRONG",
                },
            }],
        }
        decision = self.assert_decision(
            self.eligible_null(
                gene_id="GENE-UNORDERED",
                disease_id="DISEASE-UNORDERED",
                gene_specific_pvs1_rules=[spec],
            ),
            "ACTIVATED",
            "STRONG",
        )
        self.assertEqual(
            ["FRAMESHIFT", "STOP_GAIN"],
            decision.to_dict()["guideline_conditions"]["when"]["variant_consequence"],
        )

    def test_definite_vcep_mismatch_beats_earlier_unresolved_predicates(self):
        rule_mismatch = {
            "name": "Synthetic ordered-rule specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/ordered-rule-v1",
            "gene_id": "GENE-ORDERED",
            "disease_id": "DISEASE-ORDERED",
            "rules": [{
                "source_locator": "synthetic ordered rule",
                "pathway": "NONSENSE_OR_FRAMESHIFT",
                "when": {
                    "critical_protein_region": True,
                    "variant_consequence": "START_LOST",
                },
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "VERY_STRONG",
                    "decision_path": "SYNTHETIC_ORDERED_RULE_VS",
                },
            }],
        }
        self.assert_decision(
            self.eligible_null(
                gene_id="GENE-ORDERED",
                disease_id="DISEASE-ORDERED",
                critical_protein_region={"state": "conflicting", "value": True},
                gene_specific_pvs1_rules=[rule_mismatch],
            ),
            "ACTIVATED",
            "VERY_STRONG",
        )

        scope_mismatch = deepcopy(rule_mismatch)
        scope_mismatch["name"] = "Synthetic ordered scope specification"
        scope_mismatch["gene_id"] = "GENE-SCOPE-ONLY"
        scope_mismatch["disease_id"] = "DISEASE-SCOPE-ONLY"
        scope_mismatch["rules"][0]["when"] = {
            "variant_consequence": "FRAMESHIFT"
        }
        scope_evidence = self.eligible_null(
            disease_id="OTHER-DISEASE",
            gene_specific_pvs1_rules=[scope_mismatch],
        )
        scope_evidence.pop("gene_id")
        self.assert_decision(scope_evidence, "ACTIVATED", "VERY_STRONG")

    def test_known_branch_exclusions_are_not_hidden_by_missing_lof_eligibility(self):
        cases = [
            self.eligible_null(biologically_relevant_transcript=False),
            self.eligible_null(
                variant_consequence="canonical_splice_variant",
                splice_effect="frameshift_nmd",
                splice_rescue_evidence=True,
            ),
            self.cnv_evidence(
                variant_consequence="duplication",
                whole_gene_duplication=True,
            ),
            self.eligible_null(
                variant_consequence="splice_site",
                splice_site_position=5,
                rna_splicing_evidence={
                    "state": "known",
                    "value": False,
                    "provenance": [{"source": "synthetic RNA", "kind": "OBSERVED"}],
                },
                rna_assay_interpretable=True,
            ),
            self.start_loss_evidence(alternative_start_codon=True),
        ]
        for evidence in cases:
            with self.subTest(variant_consequence=evidence.get("variant_consequence")):
                evidence.pop("lof_disease_mechanism")
                self.assert_decision(evidence, "NOT_MET", "NONE")

    def test_vcep_specific_rule_precedes_any_and_unresolved_specific_blocks_it(self):
        spec = {
            "name": "Synthetic specificity precedence specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/specificity-v1",
            "gene_id": "GENE-SPECIFICITY",
            "disease_id": "DISEASE-SPECIFICITY",
            "rules": [
                {
                    "source_locator": "synthetic any fallback",
                    "pathway": "ANY",
                    "when": {"variant_consequence": "START_LOST"},
                    "decision": {
                        "status": "ACTIVATED",
                        "strength": "MODERATE",
                        "decision_path": "SYNTHETIC_ANY_MODERATE",
                    },
                },
                {
                    "source_locator": "synthetic start-loss specific rule",
                    "pathway": "INITIATION_CODON",
                    "when": {
                        "alternative_start_codon": False,
                        "upstream_pathogenic_variant_evidence": False,
                    },
                    "decision": {
                        "status": "ACTIVATED",
                        "strength": "VERY_STRONG",
                        "decision_path": "SYNTHETIC_SPECIFIC_VS",
                    },
                },
            ],
        }
        evidence = self.start_loss_evidence(
            gene_id="GENE-SPECIFICITY",
            disease_id="DISEASE-SPECIFICITY",
            gene_specific_pvs1_rules=[spec],
        )
        decision = self.assert_decision(evidence, "ACTIVATED", "VERY_STRONG")
        self.assertEqual("SYNTHETIC_SPECIFIC_VS", decision.decision_path)

        unresolved = deepcopy(evidence)
        unresolved.pop("alternative_start_codon")
        self.assert_decision(
            unresolved,
            "UNDETERMINED",
            "NONE",
            missing=("alternative_start_codon",),
        )

    def test_separate_vcep_specifications_do_not_gain_implicit_precedence(self):
        """A broad path-specific rule cannot silently outrank another VCEP's ANY rule."""

        broad_specific = {
            "name": "Synthetic broad initiation specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/broad-initiation-v1",
            "gene_id": "GENE-CROSS-SPEC",
            "disease_id": "DISEASE-CROSS-SPEC",
            "rules": [{
                "source_locator": "synthetic broad initiation rule",
                "pathway": "INITIATION_CODON",
                "when": {
                    "alternative_start_codon": False,
                    "upstream_pathogenic_variant_evidence": False,
                },
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "VERY_STRONG",
                    "decision_path": "SYNTHETIC_BROAD_INITIATION_VS",
                },
            }],
        }
        transcript_scoped_any = {
            "name": "Synthetic transcript ANY specification",
            "version": "1.0.0",
            "source_url": "https://example.invalid/pvs1/transcript-any-v1",
            "gene_id": "GENE-CROSS-SPEC",
            "disease_id": "DISEASE-CROSS-SPEC",
            "transcript_id": "TX-CROSS-SPEC.1",
            "rules": [{
                "source_locator": "synthetic transcript any rule",
                "pathway": "ANY",
                "when": {"variant_consequence": "START_LOST"},
                "decision": {
                    "status": "ACTIVATED",
                    "strength": "MODERATE",
                    "decision_path": "SYNTHETIC_TRANSCRIPT_ANY_MODERATE",
                },
            }],
        }
        self.assert_decision(
            self.start_loss_evidence(
                gene_id="GENE-CROSS-SPEC",
                disease_id="DISEASE-CROSS-SPEC",
                transcript_id="TX-CROSS-SPEC.1",
                gene_specific_pvs1_rules=[broad_specific, transcript_scoped_any],
            ),
            "UNDETERMINED",
            "NONE",
        )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
