"""Synthetic PP3 test inputs; no example is a clinical assertion."""


def missense(score=0.8, **changes):
    data = {
        "variant_id": "v1", "variant_consequence": "MISSENSE",
        "transcript_id_version": "NM_TEST.1", "biologically_relevant_transcript": True,
        "selected_missense_predictor": "REVEL",
        "genome_build": "GRCh38", "predictions": [{
            "name": "REVEL", "score": score, "version": "1", "calibration": "PEJAVER_2022",
            "variant_id": "v1", "transcript_id_version": "NM_TEST.1",
            "genome_build": "GRCh38", "quality_ok": True,
        }],
    }
    data.update(changes)
    return data


def splicing(score=0.2, **changes):
    data = {
        "variant_id": "v1", "variant_consequence": "INTRONIC",
        "transcript_id_version": "NM_TEST.1", "biologically_relevant_transcript": True,
        "genome_build": "GRCh38", "splice_context_established": True,
        "canonical_splice_site": False, "pvs1_splice_overlap": False,
        "rna_splice_evidence_used": False,
        "predictions": [{"name": "SPLICEAI", "score": score, "version": "1.3", "calibration": "WALKER_2023",
                         "variant_id": "v1", "transcript_id_version": "NM_TEST.1",
                         "genome_build": "GRCh38", "quality_ok": True}],
    }
    data.update(changes)
    return data
