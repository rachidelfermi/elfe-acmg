"""Rule-based PM5 assessment. Inputs are curated evidence, not annotations to infer."""

from dataclasses import fields
from typing import Mapping

from .guidelines import ACMG, Guideline, select_guideline
from .models import BiologicalFeatures, Candidate, PM5Result, Reference


def _record(cls, value):
    if isinstance(value, cls):
        return value
    if not isinstance(value, Mapping):
        raise TypeError(f"{cls.__name__} must be a mapping or {cls.__name__}")
    allowed = {f.name for f in fields(cls)}
    unknown = set(value) - allowed
    if unknown:
        raise ValueError(f"Unknown {cls.__name__} fields: {', '.join(sorted(unknown))}")
    return cls(**value)


def _features(value):
    if isinstance(value, BiologicalFeatures):
        record = value
    else:
        if not isinstance(value, Mapping):
            raise TypeError("biological_features must be a mapping or BiologicalFeatures")
        data = dict(value)
        data["candidate"] = _record(Candidate, data["candidate"])
        data["references"] = tuple(_record(Reference, r) for r in data.get("references", ()))
        record = _record(BiologicalFeatures, data)
    allowed_states = {"KNOWN", "UNKNOWN", "CONFLICTING", "NOT_APPLICABLE"}
    if any(state not in allowed_states for state in record.evidence_states.values()):
        raise ValueError("Evidence states must be KNOWN, UNKNOWN, CONFLICTING, or NOT_APPLICABLE")
    return record


def _result(status, strength, path, guideline, reason, used=(), features=(), missing=(), exclusions=()):
    return PM5Result("PM5", status, strength, path, tuple(used), tuple(features), tuple(sorted(set(missing))), tuple(exclusions), guideline.citation() if guideline else {}, reason)


def _base_check(c: Candidate, guideline: Guideline):
    missing = []
    if c.consequence is None:
        missing.append("candidate.consequence")
    elif c.consequence != "MISSENSE":
        return "NOT_APPLICABLE", "Candidate is not a missense substitution", missing
    for key in ("amino_acid_position", "reference_amino_acid", "alternate_amino_acid", "transcript_id_version", "transcript_relevant"):
        if getattr(c, key) is None:
            missing.append(f"candidate.{key}")
    if c.transcript_relevant is False:
        return "NOT_MET", "Candidate transcript is not disease relevant", missing
    if c.reference_amino_acid and c.alternate_amino_acid and c.reference_amino_acid == c.alternate_amino_acid:
        return "NOT_APPLICABLE", "Candidate has no amino acid substitution", missing
    if c.amino_acid_position is not None and c.amino_acid_position < 1:
        return "NOT_APPLICABLE", "Invalid protein position", missing
    if guideline.transcript and c.transcript_id_version and c.transcript_id_version != guideline.transcript:
        return "UNDETERMINED", "Candidate is not mapped to the pinned VCEP transcript", ["candidate.transcript_id_version"]
    return None, "", missing


def _reference_check(c: Candidate, r: Reference, guideline: Guideline):
    missing = []
    if r.variant_id and c.variant_id and r.variant_id == c.variant_id:
        return "NOT_MET", "Candidate cannot serve as its own reference", missing
    if r.consequence is None:
        missing.append("reference.consequence")
    elif r.consequence != "MISSENSE":
        return "NOT_MET", "Reference is not missense", missing
    for key in ("amino_acid_position", "reference_amino_acid", "alternate_amino_acid", "transcript_id_version"):
        if getattr(r, key) is None:
            missing.append(f"reference.{key}")
    if c.transcript_id_version and r.transcript_id_version and c.transcript_id_version != r.transcript_id_version:
        return "NOT_MET", "Transcript mismatch; protein coordinates cannot be presumed equivalent", missing
    if c.protein_accession_version and r.protein_accession_version and c.protein_accession_version != r.protein_accession_version:
        return "NOT_MET", "Protein isoform mismatch", missing
    if c.amino_acid_position and r.amino_acid_position and c.amino_acid_position != r.amino_acid_position:
        return "NOT_MET", "Different amino acid residue", missing
    if c.reference_amino_acid and r.reference_amino_acid and c.reference_amino_acid != r.reference_amino_acid:
        return "NOT_MET", "Reference amino acid differs", missing
    if c.alternate_amino_acid and r.alternate_amino_acid and c.alternate_amino_acid == r.alternate_amino_acid:
        return "NOT_MET", "Same amino acid substitution; assess PS1 separately", missing
    if r.classification_conflicting:
        return "UNDETERMINED", "Reference classification is conflicting", missing
    if r.classification is None:
        missing.append("reference.classification")
    elif r.classification not in ("PATHOGENIC", "LIKELY_PATHOGENIC"):
        return "NOT_MET", "Reference is not classified Pathogenic or Likely Pathogenic", missing
    for key in ("classification_source", "classification_evidence_reviewed", "classification_independent_of_pm5", "disease_match"):
        if getattr(r, key) is None:
            missing.append(f"reference.{key}")
    if r.classification_evidence_reviewed is False:
        return "UNDETERMINED", "Reference assertion has not been independently reviewed", missing
    if r.classification_independent_of_pm5 is False:
        return "NOT_MET", "Reference classification depends on PM5", missing
    if r.disease_match is False or r.mechanism_compatible is False:
        return "NOT_MET", "Reference disease or mechanism is incompatible", missing
    if r.splice_mechanism == "SPLICING":
        return "NOT_MET", "Reference pathogenicity is splice mediated, not transferable by residue", missing
    if r.splice_mechanism == "CONFLICTING":
        return "UNDETERMINED", "Reference splice mechanism is conflicting", missing
    if r.amino_acid_position is not None and r.amino_acid_position < 1:
        return "NOT_MET", "Invalid reference protein position", missing
    # Missing splice prediction is not treated as proof of normal splicing.
    return None, "", missing


def evaluate_pm5(biological_features: BiologicalFeatures | Mapping) -> PM5Result:
    """Evaluate one candidate against independently curated reference substitutions."""
    f = _features(biological_features)
    g = select_guideline(f.gene_id, f.disease_id, f.specification)
    if g is None:
        return _result("UNDETERMINED", "NONE", "GUIDELINE_SCOPE", ACMG, "No pinned guideline matches the requested gene, disease and version", missing=("applicable_guideline",))
    if f.evidence_states.get("candidate") == "NOT_APPLICABLE":
        return _result("NOT_APPLICABLE", "NONE", "CANDIDATE_ELIGIBILITY", g, "Candidate PM5 evidence is marked not applicable")
    if any(f.evidence_states.get(key) == "UNKNOWN" for key in ("candidate", "gene_disease")):
        return _result("UNDETERMINED", "NONE", "MISSING_CORE_EVIDENCE", g, "Core candidate or disease evidence is marked unknown", missing=("candidate_or_gene_disease",))
    if any(f.evidence_states.get(key) == "CONFLICTING" for key in ("candidate", "gene_disease")):
        return _result("UNDETERMINED", "NONE", "CONFLICTING_EVIDENCE", g, "Core candidate or disease evidence is conflicting")
    c = f.candidate
    status, reason, missing = _base_check(c, g)
    if status:
        return _result(status, "NONE", "CANDIDATE_ELIGIBILITY", g, reason, missing=missing)
    if missing:
        return _result("UNDETERMINED", "NONE", "CANDIDATE_ELIGIBILITY", g, "Candidate protein mapping or consequence is incomplete", missing=missing)
    context_missing = [name for name in ("gene_id", "disease_id", "gene_disease_validity") if getattr(f, name) is None]
    if context_missing or f.gene_disease_validity == "UNKNOWN":
        return _result("UNDETERMINED", "NONE", "DISEASE_CONTEXT", g, "Gene and disease context is incomplete", missing=context_missing or ("gene_disease_validity",))
    if f.gene_disease_validity in ("REFUTED", "DISPUTED"):
        return _result("NOT_MET", "NONE", "DISEASE_CONTEXT", g, "Gene-disease relationship is not established")
    if f.candidate_mechanism_compatible is False or c.splice_mechanism == "SPLICING":
        return _result("NOT_MET", "NONE", "MECHANISM_EXCLUSION", g, "Candidate mechanism is incompatible with protein-substitution PM5", exclusions=("candidate_mechanism",))
    if f.biological_evidence_conflicting or c.splice_mechanism == "CONFLICTING":
        return _result("UNDETERMINED", "NONE", "CONFLICTING_EVIDENCE", g, "Conflicting biological evidence needs resolution")
    if f.known_benign_substitution_at_residue:
        return _result("UNDETERMINED", "NONE", "BENIGN_RESIDUE_EVIDENCE", g, "Known benign substitution at residue requires variant-level review", exclusions=("benign_substitution",))
    if not f.references:
        return _result("UNDETERMINED", "NONE", "REFERENCE_MISSING", g, "No reference variant supplied", missing=("references",))
    eligible, unresolved, failures = [], list(missing), []
    for i, r in enumerate(f.references):
        state, why, absent = _reference_check(c, r, g)
        if state:
            failures.append(why)
            if state == "UNDETERMINED":
                unresolved.extend(absent or [f"references[{i}]"])
            continue
        if absent:
            unresolved.extend(absent)
            continue
        eligible.append(r)
    # A single fully qualified reference suffices; duplicate reports are collapsed by alternate amino acid.
    unique = {}
    for r in eligible:
        unique.setdefault(r.alternate_amino_acid, r)
    refs = list(unique.values())
    if not refs:
        if unresolved:
            return _result("UNDETERMINED", "NONE", "REFERENCE_REVIEW", g, "; ".join(sorted(set(failures))) or "Essential reference evidence is missing", missing=unresolved)
        return _result("NOT_MET", "NONE", "REFERENCE_EXCLUSION", g, "; ".join(sorted(set(failures))) or "No qualifying reference")
    strength = None
    selected = []
    if g.key == "HNF1A_1.1":
        pathogenic = [r for r in refs if r.classification == "PATHOGENIC"]
        if len(pathogenic) >= 2:
            strength, selected = "STRONG", pathogenic[:2]
        else:
            if pathogenic and (c.grantham_distance is None or any(r.grantham_distance is None for r in pathogenic)):
                return _result("UNDETERMINED", "NONE", "VCEP_SUBSTITUTION_COMPARISON", g, "HNF1A PM5 requires curated Grantham distances", missing=("candidate.grantham_distance", "reference.grantham_distance"))
            moderate = [r for r in pathogenic if c.grantham_distance >= r.grantham_distance]
            strength, selected = ("MODERATE", moderate[:1]) if moderate else ("SUPPORTING", refs[:1])
    elif g.key == "PLATELET_1":
        pathogenic = [r for r in refs if r.classification == "PATHOGENIC"]
        strength, selected = ("MODERATE", pathogenic[:1]) if pathogenic else ("SUPPORTING", refs[:1])
    else:
        pathogenic = [r for r in refs if r.classification == "PATHOGENIC"]
        if pathogenic:
            strength, selected = "MODERATE", pathogenic[:1]
    if strength is None:
        return _result("UNDETERMINED", "NONE", "REFERENCE_CLASSIFICATION", g, "Likely Pathogenic reference does not meet the original ACMG/AMP PM5 definition", missing=("pathogenic_reference",))
    ids = tuple(r.variant_id or f"reference:{r.alternate_amino_acid}" for r in selected)
    return _result("ACTIVATED", strength, "DIFFERENT_MISSENSE_SAME_RESIDUE", g, f"{len(selected)} distinct qualifying reference substitution(s); {g.rule_id} at {strength.lower()} strength", used=ids, features=("candidate missense", "matched transcript and residue", "different substitution", "reviewed independent reference classification", "disease match"), exclusions=("same substitution", "splicing mechanism", "circular classification", "known benign residue evidence"))
