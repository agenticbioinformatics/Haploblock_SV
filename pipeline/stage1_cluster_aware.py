"""Stage 1: infer chromosome-specific SV-to-haploblock-cluster associations.

Reads a cohort SV VCF, retrieves (or reads locally) haploblock cluster
memberships, writes the standardized SV, haploblock and cluster tables, and
records every output path in ``config.yaml`` for the downstream stages.
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
from pathlib import Path

import yaml

from cluster_inference import ALL_CHROMS, InferenceSettings
from cluster_inference import run as run_cluster_aware
from haploblock_clusters import DEFAULT_BASE_URL, request_with_retries
from stage_io import DEFAULT_BOUNDARY_DISTANCE_BP, run_cli
from sv_table import canonical_sample_id, normalize_chrom, read_vcf_samples


log = logging.getLogger("stage1_cluster_aware")

DEFAULT_SAMPLE_METADATA_URL = (
    "https://s3.amazonaws.com/1000g-ont/PROCESSED_DATA/"
    "1kGP_LRSC_500_ONT_Metadata.tsv"
)
DEFAULT_GTF_URL = (
    "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/"
    "Homo_sapiens.GRCh38.115.gtf.gz"
)
# Column names of the 1000 Genomes ONT metadata table -> pipeline contract.
METADATA_COLUMN_NAMES = {
    "NHGRI_ID": "sample_id",
    "SubPopulation": "population",
    "SuperPopulation": "superpopulation",
}


CHROM_NAME_RE = re.compile(r"chr[A-Za-z0-9_]+")


def parse_chroms(value: str) -> list[str]:
    if value.strip().lower() == "all":
        return ALL_CHROMS
    chroms = list(dict.fromkeys(
        normalize_chrom(chrom.strip()) for chrom in value.split(",") if chrom.strip()
    ))
    if not chroms:
        raise ValueError("--chroms selected no chromosomes")
    invalid = [chrom for chrom in chroms if not CHROM_NAME_RE.fullmatch(chrom)]
    if invalid:
        # Chromosome names become file names, cache directories and URLs.
        raise ValueError(f"Invalid chromosome names in --chroms: {invalid}")
    return chroms


def load_sample_metadata(
    source: Path, original_to_canonical: dict[str, str]
) -> tuple[list[str], dict[str, dict[str, str]]]:
    """Read sample metadata with contract column names, keyed by canonical sample ID."""
    with source.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = [METADATA_COLUMN_NAMES.get(name, name) for name in reader.fieldnames or []]
        missing = {"sample_id", "population"} - set(fieldnames)
        if missing:
            raise ValueError(
                f"Sample metadata {source} lacks columns {sorted(missing)} "
                f"(or their 1000 Genomes ONT names {sorted(METADATA_COLUMN_NAMES)})"
            )
        metadata_by_sample = {}
        for source_row in reader:
            row = {
                METADATA_COLUMN_NAMES.get(name, name): value
                for name, value in source_row.items()
            }
            row["sample_id"] = original_to_canonical.get(
                row["sample_id"], canonical_sample_id(row["sample_id"])
            )
            metadata_by_sample[row["sample_id"]] = row
    return fieldnames, metadata_by_sample


def check_metadata_matches_vcf(metadata_source: Path, vcf: Path) -> None:
    """Fail before inference if no VCF sample has population metadata."""
    vcf_samples = {canonical_sample_id(sample) for sample in read_vcf_samples(vcf)}
    _, metadata_by_sample = load_sample_metadata(metadata_source, {})
    matched = len(vcf_samples & set(metadata_by_sample))
    if matched == 0:
        raise ValueError(
            f"No VCF sample in {vcf} has a row in the sample metadata {metadata_source}; "
            "check that the sample IDs match"
        )
    if matched < len(vcf_samples):
        log.warning(
            "%d of %d VCF samples have no population metadata and will be ignored by Stages 4, 6 and 7",
            len(vcf_samples) - matched, len(vcf_samples),
        )


def normalize_sample_metadata(source: Path, samples_path: Path, destination: Path) -> None:
    """Rename metadata columns to the contract and keep only VCF samples, in VCF order."""
    with samples_path.open() as handle:
        sample_rows = list(csv.DictReader(handle, delimiter="\t"))
    original_to_canonical = {row["original_sample_id"]: row["sample_id"] for row in sample_rows}
    fieldnames, metadata_by_sample = load_sample_metadata(source, original_to_canonical)

    matched = 0
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sample in sample_rows:
            if sample["sample_id"] in metadata_by_sample:
                writer.writerow(metadata_by_sample[sample["sample_id"]])
                matched += 1
    if matched == 0:
        raise ValueError(f"No VCF sample has a row in the sample metadata {source}")


def fetch_default(path: Path, url: str, retries: int, description: str) -> None:
    if path.exists():
        log.info("Reusing previously downloaded %s: %s", description, path)
        return
    log.info("Downloading %s from %s", description, url)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(request_with_retries(url, retries))
    temporary.replace(path)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, help="Cohort SV VCF (.vcf or .vcf.gz), GRCh38")
    parser.add_argument(
        "--sample-metadata", type=Path, default=None,
        help="TSV with sample_id and population columns; default downloads 1000 Genomes ONT metadata",
    )
    parser.add_argument(
        "--gtf", type=Path, default=None,
        help="Gene annotation GTF for Stage 8; default downloads Ensembl GRCh38 release 115",
    )
    parser.add_argument("--out-dir", type=Path, default=Path("stage1_output"))
    parser.add_argument("--chroms", default="all", help="Comma-separated chromosomes or 'all' for chr1-22,X")
    parser.add_argument("--cluster-base-url", default=DEFAULT_BASE_URL)
    parser.add_argument(
        "--cluster-root", type=Path, default=None,
        help="Use local cluster files from <root>/<chrom>/ instead of downloading",
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--retries", type=int, default=3, help="Download attempts per file")
    parser.add_argument("--max-sv-id-length", type=int, default=80)
    parser.add_argument("--association-threshold", type=float, default=0.75)
    parser.add_argument("--posterior-threshold", type=float, default=0.75)
    parser.add_argument("--max-iterations", type=int, default=25)
    parser.add_argument("--tolerance", type=float, default=1e-5)
    parser.add_argument(
        "--boundary-distance-bp", type=int, default=DEFAULT_BOUNDARY_DISTANCE_BP,
        help="Recorded in config.yaml as the Stage 2 near-boundary distance",
    )
    parser.add_argument(
        "--test", action="store_true",
        help="Run this stage end-to-end on a tiny synthetic dataset and exit",
    )
    args = parser.parse_args(argv)
    if not args.test and args.vcf is None:
        parser.error("--vcf is required (or use --test)")
    for name in ("association_threshold", "posterior_threshold"):
        if not 0 <= getattr(args, name) <= 1:
            parser.error(f"--{name.replace('_', '-')} must be between 0 and 1")
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    if args.test:
        from synthetic_data import self_test

        self_test(1)
        return

    for label, path in (("VCF", args.vcf), ("sample metadata", args.sample_metadata),
                        ("GTF", args.gtf), ("cluster root", args.cluster_root)):
        if path is not None and not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")
    chroms = parse_chroms(args.chroms)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    metadata_source = args.sample_metadata
    if metadata_source is None:
        metadata_source = args.out_dir / Path(DEFAULT_SAMPLE_METADATA_URL).name
        fetch_default(metadata_source, DEFAULT_SAMPLE_METADATA_URL, args.retries, "sample metadata")
    check_metadata_matches_vcf(metadata_source, args.vcf)

    gtf_source = args.gtf
    if gtf_source is None:
        gtf_source = args.out_dir / Path(DEFAULT_GTF_URL).name
        fetch_default(gtf_source, DEFAULT_GTF_URL, args.retries, "gene annotation")

    log.info("Running cluster-aware preprocessing for %d chromosome(s)", len(chroms))
    run_cluster_aware(
        InferenceSettings(
            vcf=args.vcf,
            chroms=chroms,
            out_dir=args.out_dir,
            cluster_base_url=args.cluster_base_url,
            cluster_root=args.cluster_root,
            threads=args.threads,
            retries=args.retries,
            max_sv_id_length=args.max_sv_id_length,
            association_threshold=args.association_threshold,
            posterior_threshold=args.posterior_threshold,
            max_iterations=args.max_iterations,
            tolerance=args.tolerance,
        )
    )

    normalized_metadata_path = args.out_dir / "sample_metadata.tsv"
    normalize_sample_metadata(metadata_source, args.out_dir / "samples.tsv", normalized_metadata_path)

    def per_chrom(name: str) -> dict[str, str]:
        return {chrom: str((args.out_dir / f"{name}.{chrom}.tsv").resolve()) for chrom in chroms}

    config = {
        "genome_build": "GRCh38",
        "data_sources": {
            "vcf": str(args.vcf.resolve()),
            "haploblock_clusters": (
                str(args.cluster_root.resolve()) if args.cluster_root is not None else args.cluster_base_url
            ),
            "gtf": str(gtf_source.resolve()) if args.gtf is not None else DEFAULT_GTF_URL,
            "sample_metadata": (
                str(args.sample_metadata.resolve())
                if args.sample_metadata is not None else DEFAULT_SAMPLE_METADATA_URL
            ),
        },
        "thresholds": {
            "association_probability": args.association_threshold,
            "heterozygote_assignment_probability": args.posterior_threshold,
            "boundary_distance_bp": args.boundary_distance_bp,
        },
        "settings": {
            "max_sv_id_length": args.max_sv_id_length,
            "max_em_iterations": args.max_iterations,
            "em_tolerance": args.tolerance,
        },
        "paths": {
            "vcf": str(args.vcf.resolve()),
            "samples": str((args.out_dir / "samples.tsv").resolve()),
            "sample_metadata": str(normalized_metadata_path.resolve()),
            "sv_genotypes": per_chrom("sv_genotypes"),
            "sv_to_clusters": per_chrom("sv_to_clusters"),
            "haploblocks": per_chrom("haploblocks"),
            "cluster_memberships": per_chrom("cluster_memberships"),
            "sv_block_summary": per_chrom("sv_block_summary"),
            "debug_and_qc": str((args.out_dir / "debug_and_qc").resolve()),
            "gtf": str(gtf_source.resolve()),
        },
    }
    config_path = args.out_dir / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    log.info("Stage 1 complete: %s", config_path.resolve())


if __name__ == "__main__":
    run_cli(main)
