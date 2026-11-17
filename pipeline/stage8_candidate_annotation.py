"""Stage 8: add locus context to SVs captured or missed by haploblock hashes.

The same genomic overlap has different interpretations for deletions,
duplications, inversions, and insertions. This stage keeps those consequence
labels separate from association, representation, and call-quality evidence.
"""

from __future__ import annotations

import argparse
import logging
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from stage_io import load_config, read_tsv, run_cli, write_stage_config
from sv_table import METADATA_COLUMNS, normalize_chrom, open_vcf_text


log = logging.getLogger("stage8_candidate_annotation")
GTF_COLUMNS = ["chrom", "feature", "start", "end", "gene_id", "gene_name", "gene_biotype"]


def gtf_attributes(text: str) -> dict[str, str]:
    return {
        key: value
        for key, value in re.findall(r'(\S+)\s+"([^"]+)"', text)
    }


def read_gtf(path: Path) -> pd.DataFrame:
    """Read gene and exon features from a GTF as 0-based, half-open intervals."""
    if not path.is_file():
        raise FileNotFoundError(f"GTF annotation not found: {path}")
    rows = []
    with open_vcf_text(path) as handle:  # detects gzip by content, not extension
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "exon"}:
                continue
            attributes = gtf_attributes(fields[8])
            rows.append({
                "chrom": normalize_chrom(fields[0]),
                "feature": fields[2],
                "start": int(fields[3]) - 1,
                "end": int(fields[4]),
                "gene_id": attributes.get("gene_id", ""),
                "gene_name": attributes.get("gene_name", attributes.get("gene_id", "")),
                "gene_biotype": attributes.get(
                    "gene_biotype", attributes.get("gene_type", "")
                ),
            })
    if not rows:
        raise ValueError(f"No gene or exon features found in {path}; is it a GTF file?")
    return pd.DataFrame(rows, columns=GTF_COLUMNS)


def overlaps(features: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    return features[(features["start"] < end) & (features["end"] > start)]


def genes_at_breakpoint(features: pd.DataFrame, position: int) -> set[str]:
    genes = features[
        (features["feature"] == "gene")
        & (features["start"] <= position)
        & (features["end"] > position)
    ]
    return set(genes["gene_name"])


def annotate_variant(variant: pd.Series, features: pd.DataFrame) -> tuple[str, str, str]:
    start, end = int(variant["start"]), int(variant["end"])
    sv_type = str(variant["sv_type"]).upper()
    span = overlaps(features, start, max(end, start + 1))
    genes = set(span[span["feature"] == "gene"]["gene_name"])
    exon_genes = set(span[span["feature"] == "exon"]["gene_name"])

    if sv_type == "DEL":
        complete_exon_genes = set(
            features[
                (features["feature"] == "exon")
                & (features["start"] >= start)
                & (features["end"] <= end)
            ]["gene_name"]
        )
        if complete_exon_genes:
            return "complete_exon_loss", ";".join(sorted(complete_exon_genes)), "exon"
        if exon_genes:
            return "exonic_deletion", ";".join(sorted(exon_genes)), "exon"
        if genes:
            return "intragenic_deletion", ";".join(sorted(genes)), "gene"
    elif sv_type == "DUP":
        full_genes = set(
            features[
                (features["feature"] == "gene")
                & (features["start"] >= start)
                & (features["end"] <= end)
            ]["gene_name"]
        )
        if full_genes:
            return "complete_gene_duplication", ";".join(sorted(full_genes)), "gene"
        if exon_genes or genes:
            affected = exon_genes or genes
            return "partial_gene_duplication", ";".join(sorted(affected)), "exon_or_gene"
    elif sv_type == "INV":
        breakpoint_genes = genes_at_breakpoint(features, start) | genes_at_breakpoint(
            features, max(start, end - 1)
        )
        if breakpoint_genes:
            return "inversion_breakpoint_disruption", ";".join(sorted(breakpoint_genes)), "breakpoint"
        if genes:
            return "gene_within_inversion", ";".join(sorted(genes)), "span_only"
    elif sv_type == "INS":
        point = overlaps(features, start, start + 1)
        point_exons = set(point[point["feature"] == "exon"]["gene_name"])
        point_genes = set(point[point["feature"] == "gene"]["gene_name"])
        if point_exons:
            return "exonic_insertion", ";".join(sorted(point_exons)), "breakpoint"
        if point_genes:
            return "intragenic_insertion", ";".join(sorted(point_genes)), "breakpoint"
    elif exon_genes or genes:
        affected = exon_genes or genes
        return "gene_overlap", ";".join(sorted(affected)), "exon_or_gene"
    return "intergenic_or_unannotated", "", "none"


def consequence_weight(consequence: str) -> int:
    return {
        "complete_exon_loss": 3,
        "exonic_deletion": 3,
        "exonic_insertion": 3,
        "inversion_breakpoint_disruption": 3,
        "complete_gene_duplication": 2,
        "partial_gene_duplication": 2,
        "intragenic_deletion": 2,
        "intragenic_insertion": 2,
        "gene_overlap": 1,
        "gene_within_inversion": 1,
    }.get(consequence, 0)


def annotate_candidates(
    candidates: pd.DataFrame,
    features: pd.DataFrame,
    classifications: pd.DataFrame | None = None,
    representation: pd.DataFrame | None = None,
    threads: int = 1,
) -> pd.DataFrame:
    result = candidates.copy()
    features_by_chrom = {chrom: group for chrom, group in features.groupby("chrom")}
    empty_features = features.iloc[0:0]
    chunk_size = max(1, (len(result) + threads * 4 - 1) // (threads * 4))
    chunks = [result.iloc[start:start + chunk_size] for start in range(0, len(result), chunk_size)]

    def annotate_chunk(chunk: pd.DataFrame) -> list[tuple[str, str, str]]:
        return [
            annotate_variant(row, features_by_chrom.get(row["chrom"], empty_features))
            for _, row in chunk.iterrows()
        ]

    if chunks:
        with ThreadPoolExecutor(max_workers=min(threads, len(chunks))) as executor:
            annotations = [
                annotation for chunk in executor.map(annotate_chunk, chunks) for annotation in chunk
            ]
    else:
        annotations = []
    result[["consequence", "genes", "overlap_basis"]] = pd.DataFrame(
        annotations, index=result.index, columns=["consequence", "genes", "overlap_basis"]
    )
    if classifications is not None:
        keys = ["sv_record_id"]
        classification_columns = keys + [
            column for column in ("sv_class", "specific_to_population")
            if column in classifications
        ]
        result = result.merge(
            classifications[classification_columns].drop_duplicates(keys), on=keys, how="left"
        )
    if representation is not None:
        keys = ["sv_record_id", "haploblock_id"]
        representation_columns = keys + [
            column for column in (
                "representation_pattern", "n_supported_carrier_clusters",
                "n_standard_evidence_carrier_clusters",
                "top_cluster_carrier_evidence_share",
                "effective_carrier_cluster_count",
                "top_standard_evidence_cluster_id",
                "top_standard_cluster_carrier_evidence_share",
                "effective_standard_carrier_cluster_count",
                "n_mixed_clusters_meeting_count_threshold",
                "n_mixed_diplotypes_meeting_count_threshold",
                "top_supported_cluster_carrier_rate",
                "top_standard_cluster_carrier_rate",
                "top_standard_cluster_population_count",
                "top_standard_cluster_populations",
                "population_context_pattern",
            ) if column in representation
        ]
        result = result.merge(
            representation[representation_columns].drop_duplicates(keys), on=keys, how="left"
        )
    if "representation_pattern" not in result:
        result["representation_pattern"] = "not_evaluated"
    else:
        result["representation_pattern"] = result["representation_pattern"].fillna(
            "not_evaluated"
        )

    imprecise = result["imprecise"].astype(str).str.lower().isin({"true", "1"})
    passed = result["filter"].astype(str).eq("PASS")
    result["call_quality"] = np.select(
        [passed & ~imprecise, passed & imprecise, ~passed & ~imprecise],
        ["pass_precise", "pass_imprecise", "nonpass_precise"],
        default="nonpass_imprecise",
    )
    result["consequence_priority"] = result["consequence"].map(consequence_weight)
    return result.reset_index(drop=True)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("stage7_output/config.yaml"),
                        help="Stage 7 config.yaml (carries Stage 1, 4 and 6 paths)")
    parser.add_argument("--gtf", type=Path, default=None,
                        help="Gene annotation GTF; defaults to the one registered by Stage 1")
    parser.add_argument("--out-dir", type=Path, default=Path("stage8_output"))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--test", action="store_true",
                        help="Run Stages 1, 4, 6, 7 and 8 on a tiny synthetic dataset and exit")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    if args.test:
        from synthetic_data import self_test

        self_test(8)
        return

    config = load_config(args.config, ["sv_cluster_summary"] + ([] if args.gtf else ["gtf"]))
    candidates = read_tsv(config["paths"]["sv_cluster_summary"], METADATA_COLUMNS + ["haploblock_id"])
    classifications = None
    if "sv_classification" in config["paths"]:
        classifications = read_tsv(config["paths"]["sv_classification"], ["sv_record_id"])
    representation = None
    if "sv_hash_representation" in config["paths"]:
        representation = read_tsv(config["paths"]["sv_hash_representation"], ["sv_record_id", "haploblock_id"])
    gtf_path = args.gtf or Path(config["paths"]["gtf"])
    result = annotate_candidates(
        candidates, read_gtf(gtf_path), classifications, representation, max(1, args.threads)
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.out_dir / "annotated_sv_candidates.tsv"
    result.to_csv(output_path, sep="\t", index=False)
    log.info("Annotated %d SV–haploblock candidates: %s", len(result),
             result["consequence"].value_counts().to_dict() if not result.empty else {})
    write_stage_config(config, args.out_dir, paths={"annotated_sv_candidates": str(output_path.resolve())})


if __name__ == "__main__":
    run_cli(main)
