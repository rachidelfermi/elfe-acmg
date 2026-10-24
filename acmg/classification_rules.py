"""Richards et al. 2015 Table 5 and ClinGen SVI 2020 PM2 update.

Strength modifiers are accepted as reported by existing ClinGen evaluators.
No point conversion, odds model, or probability is inferred.
"""

STRENGTHS = frozenset({"VERY_STRONG", "STRONG", "MODERATE", "SUPPORTING", "STAND_ALONE", "NONE"})


def classify_criteria(activated: list[dict]) -> tuple[str, str]:
    p = [row for row in activated if row["criterion"].startswith(("PVS", "PS", "PM", "PP")) and row["strength"] != "NONE"]
    b = [row for row in activated if row["criterion"].startswith(("BA", "BS", "BP")) and row["strength"] != "NONE"]
    # Table 5's Very Strong combination is explicitly PVS1. Other upgraded
    # criteria need their applicable specification's combination rule.
    vs = sum(row["criterion"] == "PVS1" and row["strength"] == "VERY_STRONG" for row in p)
    strong = sum(row["strength"] == "STRONG" for row in p)
    moderate = sum(row["strength"] == "MODERATE" for row in p)
    supporting = sum(row["strength"] == "SUPPORTING" for row in p)
    stand_alone = any(row["strength"] == "STAND_ALONE" for row in b)
    benign_strong = sum(row["strength"] == "STRONG" for row in b)
    benign_supporting = sum(row["strength"] == "SUPPORTING" for row in b)

    # Table 5 pathogenic (i), (ii), (iii); adjusted strengths retain their
    # categories but a non-PVS1 very-strong criterion needs a separate rule.
    pathogenic = (
        (vs >= 1 and (strong >= 1 or moderate >= 2 or
                      (moderate >= 1 and supporting >= 1) or supporting >= 2))
        or strong >= 2
        or (strong >= 1 and (moderate >= 3 or
                            (moderate >= 2 and supporting >= 2) or
                            (moderate >= 1 and supporting >= 4)))
    )
    likely_pathogenic = (
        (vs >= 1 and moderate >= 1)
        or (strong >= 1 and (moderate >= 1 or supporting >= 2))
        or moderate >= 3
        or (moderate >= 2 and supporting >= 2)
        or (moderate >= 1 and supporting >= 4)
    )
    benign = stand_alone or benign_strong >= 2
    likely_benign = (benign_strong >= 1 and benign_supporting >= 1) or benign_supporting >= 2

    # Table 5 explicitly sends contradictory pathogenic and benign evidence
    # to uncertain significance. BA1 is no exception when conflict is credible.
    if p and b:
        return "VUS", "Credible pathogenic and benign criteria conflict (ACMG/AMP 2015 Table 5)."
    if pathogenic:
        return "PATHOGENIC", "Pathogenic combination in ACMG/AMP 2015 Table 5."
    if likely_pathogenic:
        return "LIKELY_PATHOGENIC", "Likely pathogenic combination in ACMG/AMP 2015 Table 5."
    if vs >= 1 and any(row["criterion"] == "PM2" and row["strength"] == "SUPPORTING" for row in p):
        return "LIKELY_PATHOGENIC", "PVS1 plus PM2_Supporting: ClinGen SVI PM2 recommendation v1.0 (2020)."
    if benign:
        return "BENIGN", "Benign combination in ACMG/AMP 2015 Table 5."
    if likely_benign:
        return "LIKELY_BENIGN", "Likely benign combination in ACMG/AMP 2015 Table 5."
    return "VUS", "Available evidence does not meet an ACMG/AMP 2015 Table 5 combination."
