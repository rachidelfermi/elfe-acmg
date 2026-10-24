"""Resolve documented reuse of the *same* observations, assays, or predictions."""


def _shared(a: dict, b: dict, field: str) -> bool:
    left, right = a.get(field), b.get(field)
    if not left or not right:
        return False
    if isinstance(left, str):
        left = [left]
    if isinstance(right, str):
        right = [right]
    return bool({str(item) for item in left} & {str(item) for item in right})


def resolve_dependencies(rows: list[dict], inputs: dict[str, object]) -> tuple[list[dict], list[dict], list[dict]]:
    """Return countable rows, recorded conflicts, and resolved overlaps.

    Separate opposing observations remain contradictory evidence for Table 5.
    Shared opposing observations are removed from both directions pending review.
    """
    active = {r["criterion"]: r for r in rows if r["status"] == "ACTIVATED"}
    suppressed: set[str] = set()
    conflicts: list[dict] = []
    resolved: list[dict] = []

    for first, second, field in (
        ("PS3", "BS3", "assays_used"),
        ("PP3", "BP4", "predictor_used"),
        ("PP1", "BS4", "families_used"),
    ):
        if first in active and second in active:
            if _shared(active[first], active[second], field):
                suppressed.update((first, second))
                resolved.append({"criteria": [first, second], "reason": f"Shared {field} cannot support opposing directions."})
            else:
                conflicts.append({"criteria": [first, second], "reason": "Opposing independent evidence requires review."})

    if "PM4" in active and "BP3" in active:
        # Both concern the same variant's protein-length event. No published
        # generic priority exists, so neither interpretation is counted.
        suppressed.update(("PM4", "BP3"))
        resolved.append({"criteria": ["PM4", "BP3"], "reason": "Protein-length interpretation is contradictory."})

    if "PVS1" in active and "PP3" in active:
        pp3_input = inputs.get("PP3")
        if isinstance(pp3_input, dict) and (pp3_input.get("pvs1_splice_overlap") is True or
                                             pp3_input.get("rna_splice_evidence_used") is True):
            suppressed.add("PP3")
            resolved.append({"criteria": ["PVS1", "PP3"], "reason": "Splice prediction/RNA evidence was already used for PVS1."})

    # The shared de novo evaluator yields one criterion, so PS2 and PM6
    # cannot count the same case twice by construction.
    if "PS2" in active or "PM6" in active:
        resolved.append({"criteria": ["PS2", "PM6"], "reason": "One shared de novo assessment assigned one criterion."})

    population = [c for c in ("BA1", "BS1", "PM2") if c in active]
    if len(population) > 1:
        # BA1 supersedes BS1; rarity against either benign frequency result
        # needs review because source and threshold reconciliation is absent.
        if "BA1" in population and "BS1" in population:
            suppressed.add("BS1")
            resolved.append({"criteria": ["BA1", "BS1"], "reason": "BA1 supersedes BS1 frequency evidence.", "unresolved": False})
        if "PM2" in population and ("BA1" in population or "BS1" in population):
            suppressed.update(population)
            conflicts.append({"criteria": population, "reason": "Contradictory validated population-frequency interpretations require source review."})

    countable = [r for r in rows if r["status"] == "ACTIVATED" and r["criterion"] not in suppressed]
    for r in rows:
        if r["criterion"] in suppressed:
            r["counted"] = False
            r["suppression_reason"] = "Evidence dependency or unresolved contradiction"
    return countable, conflicts, resolved
