"""Optional Stage 2: classify cohort SVs relative to haploblock boundaries.

Stage 1 already links SVs to associated haplotype clusters. This script is
only needed for the separate positional question. ``position_class`` records
exact overlap: zero blocks is outside, one is within, and two or more is
boundary_crossing. ``near_boundary`` separately records proximity to an edge.
The script reads the per-chromosome SV metadata and haploblock tables
registered in Stage 1's config.yaml.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from stage_io import DEFAULT_BOUNDARY_DISTANCE_BP, load_config, read_tsv, run_cli
from sv_table import METADATA_COLUMNS


log = logging.getLogger("stage2_boundary_classification")


def is_near_boundary(boundaries: np.ndarray, start: int, end: int, boundary_bp: int) -> bool:
    for position in (start, end):
        index = np.searchsorted(boundaries, position)
        if index < len(boundaries) and boundaries[index] - position < boundary_bp:
            return True
        if index and position - boundaries[index - 1] < boundary_bp:
            return True
    return False


def classify_sv_positions(sv: pd.DataFrame, hb: pd.DataFrame, boundary_bp: int) -> pd.DataFrame:
    position_class = np.full(len(sv), "outside_block", dtype=object)
    haploblock_ids = np.full(len(sv), "", dtype=object)
    near_boundary = np.full(len(sv), False, dtype=bool)

    for chrom, hb_group in hb.groupby("chrom", sort=False):
        hb_group = hb_group.sort_values("start")
        starts = hb_group["start"].to_numpy()
        ends = hb_group["end"].to_numpy()
        ids = hb_group["haploblock_id"].to_numpy()
        if len(starts) > 1 and np.any(ends[:-1] > starts[1:]):
            raise ValueError(f"Overlapping haploblocks on {chrom} are not supported")
        boundaries = np.sort(np.concatenate([starts, ends]))
        sv_mask = (sv["chrom"] == chrom).to_numpy()
        if not sv_mask.any():
            continue

        sv_starts = sv.loc[sv_mask, "start"].to_numpy()
        sv_ends = sv.loc[sv_mask, "end"].to_numpy()
        chrom_classes = np.full(sv_mask.sum(), "outside_block", dtype=object)
        chrom_ids = np.full(sv_mask.sum(), "", dtype=object)
        chrom_near_boundary = np.full(sv_mask.sum(), False, dtype=bool)

        for row_index, (sv_start, sv_end) in enumerate(zip(sv_starts, sv_ends)):
            first = np.searchsorted(ends, sv_start, side="right")
            last = np.searchsorted(starts, sv_end, side="left")
            if first < last:
                chrom_classes[row_index] = "within_block" if last - first == 1 else "boundary_crossing"
                chrom_ids[row_index] = ",".join(ids[first:last])
                chrom_near_boundary[row_index] = last - first > 1
            chrom_near_boundary[row_index] |= is_near_boundary(boundaries, sv_start, sv_end, boundary_bp)

        position_class[sv_mask] = chrom_classes
        haploblock_ids[sv_mask] = chrom_ids
        near_boundary[sv_mask] = chrom_near_boundary

    result = sv.copy()
    result["position_class"] = position_class
    result["haploblock_id"] = haploblock_ids
    result["near_boundary"] = near_boundary
    return result


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("stage1_output/config.yaml"),
                        help="Stage 1 config.yaml")
    parser.add_argument("--out-dir", type=Path, default=Path("stage2_output"))
    parser.add_argument(
        "--boundary-distance-bp", type=int, default=None,
        help="Near-boundary distance; defaults to the value recorded by Stage 1",
    )
    parser.add_argument("--test", action="store_true",
                        help="Run Stages 1-2 on a tiny synthetic dataset and exit")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    if args.test:
        from synthetic_data import self_test

        self_test(2)
        return

    config = load_config(args.config, ["haploblocks", "sv_genotypes"])
    chroms = list(config["paths"]["haploblocks"])
    boundary_bp = (
        args.boundary_distance_bp
        if args.boundary_distance_bp is not None
        else config.get("thresholds", {}).get("boundary_distance_bp", DEFAULT_BOUNDARY_DISTANCE_BP)
    )
    if boundary_bp < 0:
        raise ValueError("--boundary-distance-bp must be non-negative")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    qc = {"boundary_distance_bp": boundary_bp, "chromosomes": {}}
    for chrom in chroms:
        sv = read_tsv(config["paths"]["sv_genotypes"][chrom], METADATA_COLUMNS, usecols=METADATA_COLUMNS)
        haploblocks = read_tsv(
            config["paths"]["haploblocks"][chrom], ["haploblock_id", "chrom", "start", "end"]
        )
        annotated = classify_sv_positions(sv, haploblocks, boundary_bp)
        output_path = args.out_dir / f"boundary_svs.{chrom}.tsv"
        annotated.to_csv(output_path, sep="\t", index=False)
        counts = annotated["position_class"].value_counts().to_dict()
        qc["chromosomes"][chrom] = {
            "sv_records": len(annotated),
            "position_class_counts": counts,
            "near_boundary_records": int(annotated["near_boundary"].sum()),
            "output": str(output_path.resolve()),
        }
        log.info("%s: %s", chrom, counts)

    (args.out_dir / "boundary_qc.json").write_text(json.dumps(qc, indent=2) + "\n")


if __name__ == "__main__":
    run_cli(main)
