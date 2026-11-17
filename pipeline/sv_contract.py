"""Shared SV-table schema and VCF record parsing used by all pipeline stages."""

from __future__ import annotations

import hashlib


METADATA_COLUMNS = [
    "sv_record_id", "sv_id", "chrom", "start", "end", "sv_type", "length", "filter", "imprecise",
]


def canonical_sample_id(sample_id: str) -> str:
    """Map Coriell cell-line IDs (GM12345) to the 1000 Genomes form (NA12345)."""
    if sample_id.startswith("GM") and sample_id[2:].isdigit():
        return f"NA{sample_id[2:]}"
    return sample_id


def normalize_chrom(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def parse_info(info_text: str) -> dict[str, str | bool]:
    info: dict[str, str | bool] = {}
    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        elif item:
            info[item] = True
    return info


def parse_length(info: dict[str, str | bool], start: int, end: int) -> str:
    value = info.get("SVLEN")
    if isinstance(value, str) and value not in {"", "."}:
        try:
            return str(abs(int(value.split(",")[0])))
        except ValueError:
            pass
    interval_length = end - start
    return str(interval_length) if interval_length > 0 else ""


def simplify_sv_id(sv_id: str, chrom: str, start: int, end: int, sv_type: str, max_length: int = 80) -> str:
    if len(sv_id) <= max_length:
        return sv_id
    digest = hashlib.sha1(sv_id.encode()).hexdigest()[:10]
    return f"SV_{chrom}_{start}_{end}_{sv_type}_{digest}"


def parse_vcf_fields(
    fields: list[str], record_number: int, max_sv_id_length: int
) -> tuple[list, list[str], bool]:
    """Convert one split VCF data line into SV metadata and raw GT strings.

    ``record_number`` is the 1-based index of the record on its chromosome and
    becomes part of the unique ``sv_record_id``. Coordinates are converted to
    0-based, half-open. Returns (metadata row, genotypes, whether sv_id was shortened).
    """
    if len(fields) < 8:
        raise ValueError(f"Malformed VCF line with {len(fields)} columns: {' '.join(fields)[:120]}")
    chrom = normalize_chrom(fields[0])
    pos = int(fields[1])
    info = parse_info(fields[7])
    start = pos - 1
    end_value = info.get("END")
    end = int(end_value) if isinstance(end_value, str) and end_value != "." else pos
    sv_type = str(info.get("SVTYPE", "MISSING"))
    source_sv_id = fields[2] if fields[2] != "." else f"{chrom}:{pos}:{sv_type}:{record_number}"
    sv_id = simplify_sv_id(source_sv_id, chrom, start, end, sv_type, max_sv_id_length)
    metadata = [
        f"{chrom}_record_{record_number}", sv_id, chrom, start, end, sv_type,
        parse_length(info, start, end), fields[6], "IMPRECISE" in info,
    ]
    genotypes: list[str] = []
    if len(fields) > 9:
        format_keys = fields[8].split(":")
        if "GT" not in format_keys:
            raise ValueError(f"VCF record {source_sv_id} at {chrom}:{pos} has no GT field")
        gt_index = format_keys.index("GT")
        genotypes = [sample_field.split(":")[gt_index] for sample_field in fields[9:]]
    return metadata, genotypes, sv_id != source_sv_id
