"""Typed input and output for supplied, already-curated population evidence."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import StrEnum
from typing import Any


class EvidenceState(StrEnum):
    KNOWN = "KNOWN"
    UNKNOWN = "UNKNOWN"
    CONFLICTING = "CONFLICTING"
    NOT_APPLICABLE = "NOT_APPLICABLE"


@dataclass(frozen=True)
class AncestryRecord:
    name: str
    continental: bool = True
    allele_frequency: float | None = None
    allele_count: int | None = None
    allele_number: int | None = None
    homozygote_count: int | None = None
    hemizygote_count: int | None = None
    founder_or_bottlenecked: bool = False
    technically_reliable: bool | None = None


@dataclass(frozen=True)
class PopulationRecord:
    name: str
    version: str
    genome_assembly: str
    normalized_variant_representation: str
    allele_count: int | None = None
    allele_number: int | None = None
    allele_frequency: float | None = None
    filtering_allele_frequency: float | None = None
    homozygote_count: int | None = None
    hemizygote_count: int | None = None
    male_allele_number: int | None = None
    ancestries: tuple[AncestryRecord, ...] = ()
    ancestries_complete: bool | None = None
    filter_pass: bool | None = None
    site_callable: bool | None = None
    identity_confirmed: bool | None = None
    technically_reliable: bool | None = None
    sample_size_adequate: bool | None = None
    detection_sensitivity_adequate: bool | None = None
    genotype_quality_adequate: bool | None = None
    low_complexity_region: bool | None = None
    segmental_duplication_region: bool | None = None
    mapping_ambiguity: bool | None = None
    sampling_bias_material: bool | None = None
    case_enrichment_material: bool | None = None
    population_substructure_material: bool | None = None
    frequency_conflict: bool = False
    frequency_evidence_state: EvidenceState | None = None
    frequency_missing_reason: str | None = None
    source_date: str | None = None
    technical_limitations: str | None = None


@dataclass(frozen=True)
class PM2Input:
    variant_id: str
    genome_assembly: str
    normalized_variant_representation: str
    gene_id: str | None = None
    disease_id: str | None = None
    inheritance_mode: str | None = None
    variant_type: str = "SNV"
    chromosome: str | None = None
    position: int | None = None
    reference_allele: str | None = None
    alternate_allele: str | None = None
    variant_mapping_confidence: str | None = None
    indel_normalization_adequate: bool | None = None
    multiallelic_allele_resolved: bool | None = None
    structural_variant_frequency_resolved: bool | None = None
    mitochondrial_heteroplasmy_resolved: bool | None = None
    ploidy_context_resolved: bool | None = None
    gene_disease_validity: bool | None = None
    disease_prevalence: float | None = None
    disease_penetrance: float | None = None
    allelic_heterogeneity: float | None = None
    genetic_heterogeneity: float | None = None
    population_records: tuple[PopulationRecord, ...] = ()
    ba1_established: bool = False
    bs1_established: bool = False
    bs2_established: bool = False
    frequency_conflict: bool = False
    vcep_specification_version: str | None = None


@dataclass(frozen=True)
class PM2Result:
    criterion: str
    status: str
    strength: str
    decision_path: str | None
    population_sources_used: tuple[str, ...]
    features_used: tuple[str, ...]
    missing_critical_features: tuple[str, ...]
    exclusions_evaluated: tuple[str, ...]
    guideline: dict[str, str]
    reason: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def parse_input(value: PM2Input | dict[str, Any]) -> PM2Input:
    if isinstance(value, PM2Input):
        return value
    if not isinstance(value, dict):
        raise TypeError("evaluate_pm2 expects PM2Input or a mapping")
    data = dict(value)
    records = []
    for record in data.get("population_records", ()):
        if isinstance(record, PopulationRecord):
            records.append(record)
        else:
            item = dict(record)
            if item.get("frequency_evidence_state") is not None:
                item["frequency_evidence_state"] = EvidenceState(item["frequency_evidence_state"])
            item["ancestries"] = tuple(
                ancestry if isinstance(ancestry, AncestryRecord) else AncestryRecord(**ancestry)
                for ancestry in item.get("ancestries", ())
            )
            records.append(PopulationRecord(**item))
    data["population_records"] = tuple(records)
    return PM2Input(**data)
