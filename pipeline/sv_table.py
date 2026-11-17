"""The SV table: shared schema, VCF record parsing and VCF -> per-chromosome tables.

Every pipeline table starts with ``METADATA_COLUMNS``. Stage 1 creates the
``sv_genotypes.<chrom>.tsv`` tables here; later stages read them.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import logging
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pysam


log = logging.getLogger("sv_table")


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


class VcfRecordError(ValueError):
    """A VCF data line could not be parsed (as opposed to the file being unindexable)."""


def is_gzip_file(path: Path) -> bool:
    """True for gzip/BGZF content, whatever the file extension (.gz, .bgz, ...)."""
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def open_vcf_text(path: Path):
    return gzip.open(path, "rt") if is_gzip_file(path) else open(path, "rt")


def read_vcf_samples(path: Path) -> list[str]:
    """Sample names from the #CHROM header line, without reading any records."""
    with open_vcf_text(path) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    raise ValueError(f"{path} has no #CHROM header line; is it a VCF file?")


def _canonical_sample_ids(original_samples: list[str]) -> list[str]:
    canonical = [canonical_sample_id(sample) for sample in original_samples]
    duplicates = sorted(sample for sample, count in Counter(canonical).items() if count > 1)
    if duplicates:
        raise ValueError(f"Sample-ID canonicalization creates duplicate columns: {duplicates}")
    return canonical


def split_vcf_by_chromosome(
    vcf_path: Path,
    chroms: list[str],
    out_dir: Path,
    max_sv_id_length: int,
    chrom_workers: int = 1,
) -> dict[str, Path]:
    """Write one genotype table per chromosome, using a Tabix index when available."""
    qc_dir = out_dir / "debug_and_qc"
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    selected = set(chroms)
    paths = {chrom: out_dir / f"sv_genotypes.{chrom}.tsv" for chrom in chroms}
    handles: dict[str, object] = {}
    writers: dict[str, csv.writer] = {}
    record_counts: Counter[str] = Counter()
    shortened_counts: Counter[str] = Counter()
    filter_counts = {chrom: Counter() for chrom in chroms}
    genotype_counts = {chrom: Counter() for chrom in chroms}
    original_samples: list[str] | None = None
    canonical_samples: list[str] | None = None

    def open_writers(samples: list[str]) -> None:
        for chrom in chroms:
            handle = paths[chrom].open("w", newline="")
            handles[chrom] = handle
            writers[chrom] = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writers[chrom].writerow(METADATA_COLUMNS + samples)

    def write_record(fields: list[str]) -> None:
        try:
            metadata, genotypes, shortened = parse_vcf_fields(
                fields, record_counts[normalize_chrom(fields[0])] + 1, max_sv_id_length
            )
        except (ValueError, IndexError) as error:
            location = ":".join(fields[:2])
            raise VcfRecordError(f"Malformed VCF record at {location} in {vcf_path}: {error}") from error
        chrom = metadata[2]
        writers[chrom].writerow(metadata + genotypes)
        record_counts[chrom] += 1
        shortened_counts[chrom] += shortened
        filter_counts[chrom][fields[6]] += 1
        genotype_counts[chrom].update(genotypes)

    def reset_outputs() -> None:
        for handle in handles.values():
            handle.close()
        handles.clear()
        writers.clear()
        record_counts.clear()
        shortened_counts.clear()
        for chrom in chroms:
            filter_counts[chrom].clear()
            genotype_counts[chrom].clear()

    compressed = is_gzip_file(vcf_path)
    try:
        try:
            # Only compressed input is indexed: pysam would otherwise compress (and
            # replace) an uncompressed VCF in the user's input directory.
            if not compressed:
                raise ValueError("VCF is not compressed")
            if not Path(f"{vcf_path}.tbi").exists():
                log.info("Creating a Tabix index for %s", vcf_path)
                pysam.tabix_index(str(vcf_path), preset="vcf")
            with pysam.VariantFile(str(vcf_path)) as source:
                original_samples = list(source.header.samples)
            canonical_samples = _canonical_sample_ids(original_samples)
            open_writers(canonical_samples)

            def extract_chromosome(requested_chrom: str) -> None:
                with pysam.VariantFile(str(vcf_path)) as indexed_source:
                    for record in indexed_source.fetch(requested_chrom):
                        write_record(str(record).rstrip("\n").split("\t"))

            with ThreadPoolExecutor(max_workers=min(chrom_workers, len(chroms))) as executor:
                for future in [executor.submit(extract_chromosome, chrom) for chrom in chroms]:
                    future.result()
        except VcfRecordError:
            raise
        except (OSError, ValueError) as error:
            # Unindexable input (plain gzip, uncompressed, contig naming mismatch):
            # fall back to one sequential pass over the file.
            (log.warning if compressed else log.info)(
                "Indexed VCF access unavailable (%s); streaming the VCF instead", error
            )
            reset_outputs()
            with open_vcf_text(vcf_path) as source:
                for line in source:
                    if line.startswith("##"):
                        continue
                    if line.startswith("#CHROM"):
                        original_samples = line.rstrip("\n").split("\t")[9:]
                        canonical_samples = _canonical_sample_ids(original_samples)
                        open_writers(canonical_samples)
                        continue
                    if line.startswith("#"):
                        continue
                    if canonical_samples is None:
                        raise ValueError(f"{vcf_path} has no #CHROM header line")
                    fields = line.rstrip("\n").split("\t")
                    if normalize_chrom(fields[0]) in selected:
                        write_record(fields)
    finally:
        for handle in handles.values():
            handle.close()

    if original_samples is None or canonical_samples is None:
        raise ValueError(f"{vcf_path} has no #CHROM header line")
    with (out_dir / "samples.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "original_sample_id"])
        writer.writerows(zip(canonical_samples, original_samples))
    for chrom in chroms:
        qc = {
            "vcf": str(vcf_path.resolve()),
            "chrom": chrom,
            "records_written": record_counts[chrom],
            "sv_id_max_length": max_sv_id_length,
            "sv_ids_shortened": shortened_counts[chrom],
            "samples": len(canonical_samples),
            "sample_ids_changed": sum(a != b for a, b in zip(original_samples, canonical_samples)),
            "filter_counts": dict(sorted(filter_counts[chrom].items())),
            "genotype_counts": dict(sorted(genotype_counts[chrom].items())),
        }
        (qc_dir / f"vcf_qc.{chrom}.json").write_text(json.dumps(qc, indent=2) + "\n")
        log.info("Extracted %d records for %s", record_counts[chrom], chrom)
        if record_counts[chrom] == 0:
            log.warning("No VCF records found for %s; check chromosome naming in --chroms", chrom)
    return paths


def reusable_sv_tables(
    vcf_path: Path,
    chroms: list[str],
    out_dir: Path,
    max_sv_id_length: int,
) -> dict[str, Path] | None:
    """Return a complete interrupted-run VCF split when its QC matches this invocation."""
    qc_dir = out_dir / "debug_and_qc"
    tables = {chrom: out_dir / f"sv_genotypes.{chrom}.tsv" for chrom in chroms}
    if not (out_dir / "samples.tsv").exists():
        return None
    expected_vcf = str(vcf_path.resolve())
    for chrom, table in tables.items():
        qc_path = qc_dir / f"vcf_qc.{chrom}.json"
        if not table.exists() or not qc_path.exists():
            return None
        qc = json.loads(qc_path.read_text())
        if qc.get("vcf") != expected_vcf or qc.get("sv_id_max_length") != max_sv_id_length:
            return None
    log.info("Reusing the completed chromosome-specific VCF split from an interrupted run")
    return tables
