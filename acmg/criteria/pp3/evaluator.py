"""Deterministic, offline PP3 evaluator for supplied evidence."""

from dataclasses import fields, replace
from math import isfinite
from typing import Mapping

from .guidelines import PEJAVER, RUNX1, VCEPS, WALKER, revel_band, spliceai_band
from .models import Decision, EvidenceState, Features, Prediction, Result, Strength


def _out(status: Decision, reason: str, *, strength=Strength.NONE, path=None,
         predictor=None, used=(), missing=(), exclusions=(), guideline=None,
         bp4=False) -> Result:
    return Result(status, strength, path, predictor, tuple(used), tuple(missing),
                  tuple(exclusions), guideline, reason, bp4)


def _parse(value: Features | Mapping) -> Features:
    if isinstance(value, Features):
        raw = {field.name: getattr(value, field.name) for field in fields(Features)}
    elif isinstance(value, Mapping):
        raw = dict(value)
    else:
        raise TypeError("biological_features must be Features or a mapping")
    unknown = set(raw) - {field.name for field in fields(Features)}
    if unknown:
        raise ValueError(f"Unsupported input fields: {', '.join(sorted(unknown))}")
    raw["predictions"] = tuple(p if isinstance(p, Prediction) else Prediction(**p)
                               for p in raw.get("predictions", ()))
    raw["predictions"] = tuple(replace(p, state=EvidenceState(p.state)) for p in raw["predictions"])
    raw["missing_states"] = {k: EvidenceState(v) for k, v in raw.get("missing_states", {}).items()}
    return Features(**raw)


def _problem(p: Prediction, f: Features) -> str | None:
    if p.state != EvidenceState.KNOWN:
        return f"{p.name} prediction is {p.state.value}"
    if p.score is None:
        return f"{p.name} score is missing"
    if isinstance(p.score, bool) or not isinstance(p.score, (float, int)) or not isfinite(p.score) or not 0 <= p.score <= 1:
        return f"{p.name} score is invalid or outside [0, 1]"
    if p.quality_ok is False or p.applicable is False or p.out_of_distribution is True:
        return f"{p.name} prediction is technically unreliable or inapplicable"
    if f.variant_id and p.variant_id != f.variant_id:
        return f"{p.name} variant identity is missing or mismatched"
    if f.transcript_id_version and p.transcript_id_version != f.transcript_id_version:
        return f"{p.name} transcript is missing or mismatched"
    if f.genome_build and p.genome_build != f.genome_build:
        return f"{p.name} genome build is missing or mismatched"
    if not p.version:
        return f"{p.name} predictor version is missing"
    return None


def evaluate_pp3(biological_features: Features | Mapping) -> Result:
    """Apply pinned, variant-class-specific PP3 rules without stacking scores."""
    f = _parse(biological_features)
    consequence = (f.variant_consequence or "").upper()
    if not consequence:
        return _out(Decision.UNDETERMINED, "Variant consequence is required.", missing=("variant_consequence",))
    if f.missing_states.get("variant_consequence") == EvidenceState.CONFLICTING:
        return _out(Decision.UNDETERMINED, "Variant consequence conflicts.")
    if f.gene_disease_validity is False or f.mechanism_compatible is False or f.biologically_relevant_transcript is False:
        return _out(Decision.NOT_APPLICABLE, "Gene-disease relation, mechanism, or transcript is incompatible.")
    if f.conflicting_computational_evidence or f.bp4_evidence_conflict or any(
            state == EvidenceState.CONFLICTING for state in f.missing_states.values()):
        return _out(Decision.UNDETERMINED, "Applicable computational evidence conflicts without a resolution rule.")
    if f.biologically_relevant_transcript is not True or not f.transcript_id_version:
        return _out(Decision.UNDETERMINED, "A relevant versioned transcript is required.",
                    missing=("transcript_id_version", "biologically_relevant_transcript"))
    vcep = None
    if f.vcep_id or f.vcep_version:
        vcep = VCEPS.get((f.vcep_id, f.vcep_version))
        if vcep is None:
            return _out(Decision.UNDETERMINED, "VCEP version is not pinned.", missing=("supported_vcep_version",))
        if (f.gene_id or "").upper() != vcep.gene or f.disease_id != vcep.disease:
            return _out(Decision.NOT_APPLICABLE, "VCEP gene or disease scope does not match.")
        if f.transcript_id_version and f.transcript_id_version != vcep.transcript:
            return _out(Decision.UNDETERMINED, "VCEP transcript does not match.", missing=("RUNX1c_transcript",))
    elif (f.gene_id or "").upper() == "RUNX1":
        return _out(Decision.UNDETERMINED, "RUNX1 requires its pinned VCEP rule.", missing=("vcep_id", "vcep_version"))

    predictions: dict[str, Prediction] = {}
    for p in f.predictions:
        name = p.name.upper()
        if name in predictions:
            return _out(Decision.UNDETERMINED, f"Duplicate {name} predictions need reconciliation.")
        predictions[name] = p
    if not predictions:
        return _out(Decision.UNDETERMINED, "No calibrated prediction was supplied.", missing=("calibrated_prediction",))
    candidates: list[Result] = []
    problems: list[str] = []
    selected = (f.selected_missense_predictor or "").upper()
    if selected:
        if consequence != "MISSENSE":
            problems.append("Missense predictor is inapplicable to this variant class")
        elif selected != "REVEL":
            problems.append(f"No pinned calibration for {selected}")
        elif selected not in predictions:
            problems.append("Selected REVEL prediction is missing")
        else:
            p = predictions[selected]
            issue = _problem(p, f)
            if issue:
                problems.append(issue)
            elif p.calibration not in ({"RUNX1_GN008_3.1", "PEJAVER_2022"} if vcep else {"PEJAVER_2022"}):
                problems.append("REVEL calibration version mismatch")
            else:
                band, strength = revel_band(p.score)
                if vcep:
                    band, strength = ("PP3", "SUPPORTING") if p.score >= vcep.missense_min else ("NEUTRAL", "NONE")
                status = Decision.ACTIVATED if band == "PP3" else Decision.NOT_MET
                candidates.append(_out(status, f"REVEL {p.score} {'meets' if status == Decision.ACTIVATED else 'does not meet'} the PP3 threshold.",
                    strength=Strength(strength), path="CALIBRATED_MISSENSE_PREDICTOR", predictor="REVEL",
                    used=("variant_consequence", "transcript_id_version", "predictions.REVEL.score"),
                    exclusions=("variant_identity", "transcript_match", "technical_quality"),
                    guideline=RUNX1 if vcep else PEJAVER, bp4=band == "BP4_CANDIDATE"))

    elif "REVEL" in predictions:
        problems.append("Prospective missense predictor selection is missing")

    p = predictions.get("SPLICEAI")
    if p:
        if f.canonical_splice_site is not False:
            problems.append("Canonical ±1/2 splice sites are excluded")
        elif consequence not in {"MISSENSE", "SYNONYMOUS", "INTRONIC", "SPLICE_REGION", "NONCODING"}:
            problems.append("SpliceAI PP3 is inapplicable to this variant class")
        elif vcep and (consequence not in {"MISSENSE", "SYNONYMOUS", "INTRONIC"} or
                       (consequence == "INTRONIC" and f.intron_number not in range(4, 9))):
            problems.append("RUNX1 VCEP excludes this consequence or intron")
        elif f.biologically_relevant_transcript is not True or not f.transcript_id_version or f.splice_context_established is not True:
            problems.append("Splice context is unverified")
        elif f.pvs1_splice_overlap is not False or f.rna_splice_evidence_used is not False:
            problems.append("Splice prediction overlaps PVS1 or RNA evidence")
        else:
            issue = _problem(p, f)
            if issue:
                problems.append(issue)
            elif p.calibration not in ({"RUNX1_GN008_3.1", "WALKER_2023"} if vcep else {"WALKER_2023"}):
                problems.append("SpliceAI calibration version mismatch")
            else:
                band = ("PP3" if p.score >= vcep.splice_min else "BP4_CANDIDATE" if p.score <= 0.20 else "NEUTRAL") if vcep else spliceai_band(p.score)
                status = Decision.ACTIVATED if band == "PP3" else Decision.NOT_MET
                candidates.append(_out(status, f"SpliceAI {p.score} {'meets' if status == Decision.ACTIVATED else 'does not meet'} the PP3 threshold.",
                    strength=Strength.SUPPORTING if status == Decision.ACTIVATED else Strength.NONE,
                    path="CALIBRATED_SPLICING_PREDICTOR", predictor="SPLICEAI",
                    used=("variant_consequence", "transcript_id_version", "predictions.SPLICEAI.score"),
                    exclusions=("canonical_splice_site", "pvs1_splice_overlap", "rna_splice_evidence_used", "technical_quality"),
                    guideline=RUNX1 if vcep else WALKER,
                    bp4=band == "BP4_CANDIDATE" and (consequence in {"SYNONYMOUS", "INTRONIC"} or
                        bool(vcep and consequence == "MISSENSE" and predictions.get("REVEL") and
                             predictions["REVEL"].score is not None and predictions["REVEL"].score < 0.50))))

    positive = [candidate for candidate in candidates if candidate.status == Decision.ACTIVATED]
    if positive and any(candidate.bp4_candidate for candidate in candidates):
        return _out(Decision.UNDETERMINED, "Calibrated PP3 and BP4 evidence conflict.")
    if positive:
        return positive[0]  # Fixed pathway order: selected REVEL, then SpliceAI; one PP3.
    if problems:
        return _out(Decision.UNDETERMINED, "; ".join(problems), missing=problems)
    if candidates:
        return candidates[0]
    return _out(Decision.UNDETERMINED, "No applicable calibrated prediction was supplied.", missing=("calibrated_prediction",))
