"""Small pathway evaluator. Scores are never combined or imputed."""
from __future__ import annotations

import math
from .guidelines import BAYESDEL, PALB2, PALB2_DISEASE, PALB2_GENE, PALB2_TRANSCRIPT, PEJAVER, REVEL, RPE65, RPE65_DISEASE, RPE65_GENE, RPE65_TRANSCRIPT, SPLICEAI, WALKER, Rule
from .models import BP4Result, BiologicalFeatures, Prediction, State, Status, Strength


def _result(status: Status, reason: str, *, strength: Strength = Strength.NONE,
            path: str | None = None, prediction: Prediction | None = None,
            used: tuple[str, ...] = (), missing: tuple[str, ...] = (),
            exclusions: tuple[str, ...] = (), guideline=None) -> BP4Result:
    return BP4Result("BP4", status, strength, path, prediction.name if prediction else None,
                     used, missing, exclusions, guideline, reason)


def _valid(prediction: Prediction, rule: Rule) -> tuple[bool, str, bool]:
    """Return valid, reason, and whether missing information might resolve it."""
    if prediction.score is None:
        return False, "prediction score missing", True
    if not isinstance(prediction.score, (int, float)) or not math.isfinite(prediction.score):
        return False, "prediction score invalid", True
    if rule.name in {"REVEL", "SpliceAI"} and not 0 <= prediction.score <= 1:
        return False, "prediction score outside 0–1", True
    if rule.name == "BayesDel_noAF" and not -1 <= prediction.score <= 1:
        return False, "prediction score outside -1–1", True
    if rule.allowed_versions and prediction.version not in rule.allowed_versions:
        return False, "predictor version does not match calibration", True
    for field in ("variant_match", "transcript_match", "genome_build_match", "model_applicability", "quality"):
        state = getattr(prediction, field)
        if state != State.KNOWN:
            return False, f"{field} is {state}", state != State.NOT_APPLICABLE
    if prediction.out_of_distribution:
        return False, "prediction is out of distribution", True
    return True, "", False


def _find(features: BiologicalFeatures, name: str) -> Prediction | None:
    matches = [p for p in features.predictions if p.name.casefold() == name.casefold()]
    return matches[0] if len(matches) == 1 else None


def _evaluate_rule(features: BiologicalFeatures, rule: Rule, path: str) -> BP4Result:
    prediction = _find(features, rule.name)
    if prediction is None:
        return _result(Status.UNDETERMINED, f"{rule.name} prediction is missing or duplicated.", missing=(f"{rule.name} score",), guideline=rule.guideline)
    valid, why, unresolved = _valid(prediction, rule)
    if not valid:
        return _result(Status.UNDETERMINED if unresolved else Status.NOT_MET, why, prediction=prediction, missing=(why,) if unresolved else (), guideline=rule.guideline)
    strength = rule.strength(prediction.score)
    if strength == Strength.NONE:
        return _result(Status.NOT_MET, f"{rule.name} score {prediction.score} does not meet its benign threshold.", prediction=prediction, used=(f"{rule.name} score",), guideline=rule.guideline)
    return _result(Status.ACTIVATED, f"{rule.name} score {prediction.score} meets the published BP4 {strength.value.lower()} interval.", strength=strength, path=path, prediction=prediction, used=("consequence", "transcript_id_version", f"{rule.name} score", "prediction identity and mapping"), guideline=rule.guideline)


def _rpe65(features: BiologicalFeatures) -> BP4Result:
    if features.vcep_version != "1.0":
        return _result(Status.UNDETERMINED, "RPE65 VCEP version must be 1.0.", missing=("vcep_version",), guideline=RPE65)
    if features.disease_id not in RPE65_DISEASE:
        return _result(Status.UNDETERMINED, "RPE65 disease scope is unresolved.", missing=("disease_id",), guideline=RPE65)
    if features.transcript_id_version != RPE65_TRANSCRIPT:
        return _result(Status.UNDETERMINED, "RPE65 VCEP requires NM_000329.3.", missing=("transcript_id_version",), guideline=RPE65)
    if features.consequence == "missense":
        revel, splice = _find(features, "REVEL"), _find(features, "SpliceAI")
        if revel is None or splice is None:
            return _result(Status.UNDETERMINED, "RPE65 missense BP4 requires both REVEL and SpliceAI.", missing=tuple(x for x, p in (("REVEL score", revel), ("SpliceAI score", splice)) if p is None), guideline=RPE65)
        for p, rule in ((revel, REVEL), (splice, SPLICEAI)):
            valid, why, _ = _valid(p, rule)
            if not valid:
                return _result(Status.UNDETERMINED, why, prediction=p, missing=(why,), guideline=RPE65)
        if splice.score > 0.1 or revel.score > 0.290:
            return _result(Status.NOT_MET, "RPE65 requires SpliceAI ≤0.1 and REVEL ≤0.290.", used=("REVEL score", "SpliceAI score"), guideline=RPE65)
        strength = Strength.MODERATE if revel.score <= 0.183 else Strength.SUPPORTING
        return _result(Status.ACTIVATED, "RPE65 VCEP missense thresholds met by both predictors.", strength=strength, path="VCEP_RPE65_MISSENSE", prediction=revel, used=("gene_id", "disease_id", "transcript_id_version", "REVEL score", "SpliceAI score"), guideline=RPE65)
    if features.consequence in {"intronic", "synonymous", "noncanonical_splice"}:
        if features.splice_region_eligible is not True:
            return _result(Status.UNDETERMINED, "RPE65 splice-only site eligibility must be verified.", missing=("splice_region_eligible",), guideline=RPE65)
        if features.protein_effect_relevant is not False:
            return _result(Status.UNDETERMINED, "Protein effect must be excluded for splice-only BP4.", missing=("protein_effect_relevant",), guideline=RPE65)
        result = _evaluate_rule(features, SPLICEAI, "VCEP_RPE65_SPLICE")
        if result.status == Status.ACTIVATED:
            return _result(Status.ACTIVATED, "RPE65 splice-only BP4 threshold met.", strength=Strength.SUPPORTING, path="VCEP_RPE65_SPLICE", prediction=_find(features, "SpliceAI"), used=result.features_used, guideline=RPE65)
        return _result(result.status, result.reason, prediction=_find(features, "SpliceAI"), missing=result.missing_critical_features, guideline=RPE65)
    return _result(Status.NOT_APPLICABLE, "RPE65 BP4 has no implemented pathway for this consequence.", guideline=RPE65)


def _palb2(features: BiologicalFeatures) -> BP4Result:
    if features.vcep_version != "1.3" or features.disease_id not in PALB2_DISEASE or features.transcript_id_version != PALB2_TRANSCRIPT:
        return _result(Status.UNDETERMINED, "PALB2 VCEP 1.3 disease/transcript scope is unresolved.", missing=("vcep_version", "disease_id", "transcript_id_version"), guideline=PALB2)
    if features.consequence == "missense":
        return _result(Status.NOT_APPLICABLE, "PALB2 VCEP prohibits missense BP4.", guideline=PALB2)
    if features.consequence in {"noncanonical_splice", "intronic", "synonymous"}:
        if features.splice_region_eligible is not True or features.protein_effect_relevant is not False:
            return _result(Status.UNDETERMINED, "PALB2 splice eligibility or protein effect is unresolved.", missing=("splice_region_eligible", "protein_effect_relevant"), guideline=PALB2)
        result = _evaluate_rule(features, SPLICEAI, "VCEP_PALB2_SPLICE")
        if result.status == Status.ACTIVATED:
            return _result(Status.ACTIVATED, "PALB2 VCEP SpliceAI ≤0.1 supports BP4.", strength=Strength.SUPPORTING, path="VCEP_PALB2_SPLICE", prediction=_find(features, "SpliceAI"), used=result.features_used, guideline=PALB2)
        return _result(result.status, result.reason, prediction=_find(features, "SpliceAI"), missing=result.missing_critical_features, guideline=PALB2)
    return _result(Status.NOT_APPLICABLE, "PALB2 BP4 has no implemented pathway for this consequence.", guideline=PALB2)


def evaluate_bp4(biological_features: BiologicalFeatures) -> BP4Result:
    """Evaluate supplied evidence. Caller must establish variant identity and mapping."""
    f = biological_features
    if not isinstance(f, BiologicalFeatures):
        raise TypeError("biological_features must be BiologicalFeatures")
    if not f.variant_id or not f.gene_id or not f.consequence or not f.transcript_id_version:
        return _result(Status.UNDETERMINED, "Variant identity, gene, consequence and transcript are required.", missing=tuple(k for k in ("variant_id", "gene_id", "consequence", "transcript_id_version") if not getattr(f, k)))
    if f.biologically_relevant_transcript != State.KNOWN:
        return _result(Status.UNDETERMINED, "Relevant transcript is unresolved.", missing=("biologically_relevant_transcript",))
    exclusions = ("PP3 overlap", "BP7 overlap", "splice evidence overlap", "conflicting evidence")
    if f.pp3_evidence_conflict or f.bp7_same_prediction_prohibited or f.splice_evidence_double_counting or f.conflicting_evidence:
        return _result(Status.UNDETERMINED, "Computational evidence conflict or prohibited reuse must be resolved.", exclusions=exclusions)
    if f.gene_id in RPE65_GENE:
        return _rpe65(f)
    if f.gene_id in PALB2_GENE:
        return _palb2(f)
    if f.consequence == "missense":
        # A damaging calibrated SpliceAI result contradicts a benign missense call.
        splice = _find(f, "SpliceAI")
        if splice:
            valid, _, _ = _valid(splice, SPLICEAI)
            if valid and splice.score >= 0.2:
                return _result(Status.UNDETERMINED, "Credible SpliceAI spliceogenic prediction conflicts with missense benign evidence.", prediction=splice, exclusions=("conflicting splicing prediction",), guideline=WALKER)
        revel, bayes = _find(f, "REVEL"), _find(f, "BayesDel_noAF")
        valid = [(p, rule) for p, rule in ((revel, REVEL), (bayes, BAYESDEL)) if p is not None and _valid(p, rule)[0]]
        if valid:
            benign = [(p, rule) for p, rule in valid if rule.strength(p.score) != Strength.NONE]
            pathogenic = any((rule.name == "REVEL" and p.score >= 0.644) or (rule.name == "BayesDel_noAF" and p.score >= 0.13) for p, rule in valid)
            if benign and pathogenic:
                return _result(Status.UNDETERMINED, "Calibrated missense predictors give conflicting BP4/PP3 directions.", exclusions=("conflicting predictions",), guideline=PEJAVER)
            if benign:
                # Fixed preference avoids choosing whichever benign score is strongest.
                selected = next(((p, rule) for p, rule in benign if rule.name == "REVEL"), benign[0])
                return _evaluate_rule(f, selected[1], "CALIBRATED_MISSENSE_PREDICTOR")
            return _result(Status.NOT_MET, "Available calibrated missense scores do not meet a BP4 threshold.", used=tuple(f"{rule.name} score" for _, rule in valid), guideline=PEJAVER)
        return _result(Status.UNDETERMINED, "No valid calibrated missense predictor supplied.", missing=("valid calibrated missense score",), guideline=PEJAVER)
    if f.consequence in {"noncanonical_splice", "intronic", "synonymous"}:
        if f.splice_region_eligible is not True:
            return _result(Status.UNDETERMINED, "Noncanonical splice-site eligibility must be verified.", missing=("splice_region_eligible",), guideline=WALKER)
        if f.protein_effect_relevant is not False:
            return _result(Status.UNDETERMINED, "SpliceAI alone cannot exclude a relevant protein effect.", missing=("protein_effect_relevant",), guideline=WALKER)
        return _evaluate_rule(f, SPLICEAI, "CALIBRATED_SPLICING_PREDICTOR")
    return _result(Status.NOT_APPLICABLE, "No verified BP4 pathway for this consequence.")
