"""Stage 4: classify cohort SVs by population allele frequency.

Reads every chromosome-specific ``sv_genotypes`` table and the normalized
``sample_metadata`` table registered in Stage 1's config. Population labels
are independent of haploblock cluster labels so Stage 6 can compare them
without circularity.

Outputs:
  - sv_af_classification.tsv: one row per SV and population
  - sv_classification.tsv: one row per SV
  - sv_classification_haploblocks.tsv: one row per SV with haploblock assignments
  - haploblock_population_specific_summary.tsv: population-specific SVs per haploblock
  - population_specific_summary.tsv: population-specific SVs per population
  - stage4_summary.tsv: overall counts
  - stage4_plots/<chrom>/*.png (with --plots)
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from genotypes import genotype_counts, genotype_samples
from stage_io import load_config, read_sample_metadata, read_tsv, run_cli, write_stage_config
from sv_table import METADATA_COLUMNS


log = logging.getLogger("stage4_classify_af")
HAPLOBLOCK_COLUMNS = ["sv_record_id", "haploblock_id", "block_start", "block_end"]


def classify(
    sv: pd.DataFrame,
    metadata: pd.DataFrame,
    af_threshold: float,
    absent_af_threshold: float,
    min_samples_per_population: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return separate per-population AF and per-SV classification tables."""
    sv = sv.reset_index(drop=True)
    population_by_sample = dict(
        zip(metadata["sample_id"].astype(str), metadata["population"].astype(str))
    )
    samples = genotype_samples(sv, metadata)
    samples_by_population = {
        population: [sample for sample in samples if population_by_sample[sample] == population]
        for population in sorted({population_by_sample[sample] for sample in samples})
    }
    population_af = {}
    population_called_samples = {}
    population_called_alleles = {}
    for population, population_samples in samples_by_population.items():
        alt_sum = np.zeros(len(sv), dtype=int)
        allele_sum = np.zeros(len(sv), dtype=int)
        called_samples = np.zeros(len(sv), dtype=int)
        for sample in population_samples:
            alternate, called = genotype_counts(sv[sample])
            alt_sum += alternate
            allele_sum += called
            called_samples += called > 0
        with np.errstate(invalid="ignore", divide="ignore"):
            population_af[population] = np.where(allele_sum > 0, alt_sum / allele_sum, np.nan)
        population_called_samples[population] = called_samples
        population_called_alleles[population] = allele_sum

    classes = []
    specific_populations = []
    other_reasons = []
    populations = list(samples_by_population)
    for row_index in range(len(sv)):
        populations_with_data = [
            population
            for population in populations
            if population_called_samples[population][row_index] >= min_samples_per_population
        ]
        populations_above_threshold = [
            population
            for population in populations_with_data
            if population_af[population][row_index] >= af_threshold
        ]
        other_populations = [
            population
            for population in populations_with_data
            if population not in populations_above_threshold
        ]

        if len(populations_with_data) < 2:
            sv_class, specific_population, reason = "other", "", "insufficient_population_data"
        elif len(populations_above_threshold) >= 2:
            sv_class, specific_population, reason = "common", "", ""
        elif len(populations_above_threshold) == 1 and all(
            population_af[population][row_index] < absent_af_threshold
            for population in other_populations
        ):
            sv_class = "population_specific"
            specific_population = populations_above_threshold[0]
            reason = ""
        elif len(populations_above_threshold) == 1:
            sv_class, specific_population, reason = (
                "other", "", "one_population_high_plus_intermediate_elsewhere"
            )
        else:
            sv_class, specific_population, reason = "other", "", "absent_or_rare"

        classes.append(sv_class)
        specific_populations.append(specific_population)
        other_reasons.append(reason)

    classification_rows = []
    population_rows = []
    for row_index, variant in sv.iterrows():
        classification_rows.append(
            {
                **{column: variant[column] for column in METADATA_COLUMNS},
                "sv_class": classes[row_index],
                "specific_to_population": specific_populations[row_index],
                "other_reason": other_reasons[row_index],
            }
        )
        for population, population_samples in samples_by_population.items():
            af = population_af[population][row_index]
            population_rows.append(
                {
                    **{column: variant[column] for column in METADATA_COLUMNS},
                    "population": population,
                    "n_samples": len(population_samples),
                    "n_called": int(population_called_samples[population][row_index]),
                    "called_alleles": int(population_called_alleles[population][row_index]),
                    "pop_has_data": bool(
                        population_called_samples[population][row_index]
                        >= min_samples_per_population
                    ),
                    "af": float(af) if np.isfinite(af) else np.nan,
                }
            )
    return pd.DataFrame(population_rows), pd.DataFrame(classification_rows)

def load_haploblock_assignments(config: dict) -> pd.DataFrame:
    """Load one row per SV-haploblock overlap from Stage 1 summaries."""
    tables = [
        read_tsv(path, HAPLOBLOCK_COLUMNS)[HAPLOBLOCK_COLUMNS].drop_duplicates()
        for path in config["paths"]["sv_block_summary"].values()
    ]
    if not tables:
        return pd.DataFrame(columns=HAPLOBLOCK_COLUMNS)
    return pd.concat(tables, ignore_index=True).drop_duplicates()


def write_informative_summary(classification_haploblocks: pd.DataFrame, output_dir: Path) -> None:
    """Write overall, per-population and per-haploblock summaries of the SV classes."""
    data = classification_haploblocks
    sv_level = data.groupby("sv_record_id").agg(
        sv_class=("sv_class", "first"), has_haploblock=("haploblock_id", lambda ids: ids.notna().any())
    )
    overall_summary = pd.DataFrame([{
        "total_unique_svs": len(sv_level),
        "common_svs": int(sv_level["sv_class"].eq("common").sum()),
        "population_specific_svs": int(sv_level["sv_class"].eq("population_specific").sum()),
        "other_svs": int(sv_level["sv_class"].eq("other").sum()),
        "svs_overlapping_haploblocks": int(sv_level["has_haploblock"].sum()),
        "svs_outside_haploblocks": int((~sv_level["has_haploblock"]).sum()),
    }])
    overall_summary.to_csv(output_dir / "stage4_summary.tsv", sep="\t", index=False)

    # Each SV counted once per haploblock it overlaps.
    in_blocks = data.dropna(subset=["haploblock_id"]).drop_duplicates(["sv_record_id", "haploblock_id"])
    specific = in_blocks[in_blocks["sv_class"] == "population_specific"]

    population_summary = (
        specific.groupby("specific_to_population")
        .agg(population_specific_svs=("sv_record_id", "nunique"), haploblocks=("haploblock_id", "nunique"))
        .rename_axis("population").reset_index()
    )
    population_summary.to_csv(output_dir / "population_specific_summary.tsv", sep="\t", index=False)

    populations_per_block = specific.groupby("haploblock_id")["specific_to_population"].agg(
        lambda populations: ", ".join(sorted(set(populations)))
    )
    block_summary = (
        in_blocks.assign(
            common=in_blocks["sv_class"].eq("common"),
            population_specific=in_blocks["sv_class"].eq("population_specific"),
        )
        .groupby("haploblock_id")
        .agg(
            total_svs=("sv_record_id", "size"),
            common_svs=("common", "sum"),
            population_specific_svs=("population_specific", "sum"),
            block_start=("block_start", "first"),
            block_end=("block_end", "first"),
        )
        .join(populations_per_block.rename("populations_with_specific_svs"))
        .fillna({"populations_with_specific_svs": ""})
        .reset_index()
        [["haploblock_id", "total_svs", "common_svs", "population_specific_svs",
          "populations_with_specific_svs", "block_start", "block_end"]]
        .sort_values(["population_specific_svs", "total_svs"], ascending=False)
    )
    block_summary.to_csv(output_dir / "haploblock_population_specific_summary.tsv", sep="\t", index=False)

    log.info("Stage 4 summary: %s", overall_summary.iloc[0].to_dict())


def write_stage4_plots(
    classification_haploblocks: pd.DataFrame,
    output_dir: Path,
    chrom: str,
    top_haploblocks: int,
) -> None:
    """Write chromosome-specific plots for SV class and population specificity."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chrom_data = classification_haploblocks[
        classification_haploblocks["chrom"] == chrom
    ].copy()

    chrom_data = chrom_data.dropna(subset=["haploblock_id"])

    if chrom_data.empty:
        log.warning("No haploblock-overlapping SVs available for %s plots", chrom)
        return

    plot_dir = output_dir / chrom
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Count each SV once within each haploblock.
    unique_sv_block = chrom_data.drop_duplicates(
        ["sv_record_id", "haploblock_id"]
    )

    block_counts = (
        unique_sv_block.groupby(
            ["haploblock_id", "sv_class"],
            as_index=False,
        )
        .size()
        .pivot(
            index="haploblock_id",
            columns="sv_class",
            values="size",
        )
        .fillna(0)
    )

    for column in ["common", "population_specific", "other"]:
        if column not in block_counts:
            block_counts[column] = 0

    block_counts = block_counts[["common", "population_specific", "other"]]
    block_counts["total_svs"] = block_counts.sum(axis=1)
    block_counts = block_counts.sort_index()

    # Plot 1: class composition across all haploblocks.
    ax = block_counts[
        ["common", "population_specific", "other"]
    ].plot(
        kind="bar",
        stacked=True,
        figsize=(18, 7),
        color=["#4C78A8", "#E45756", "#B9B9B9"],
        width=0.9,
    )

    ax.set_title(f"{chrom}: SV classification by haploblock")
    ax.set_xlabel("Haploblock")
    ax.set_ylabel("Unique SV count")

    # Show only about 30 readable haploblock labels while retaining
    # every haploblock bar in the plot.
    label_step = max(1, len(block_counts) // 30)
    label_positions = list(range(0, len(block_counts), label_step))
    label_values = [
        str(block_counts.index[position])
        for position in label_positions
    ]

    ax.set_xticks(label_positions)
    ax.set_xticklabels(
        label_values,
        rotation=45,
        ha="right",
        fontsize=7,
    )

    ax.legend(title="SV class")
    ax.figure.tight_layout()
    ax.figure.savefig(
        plot_dir / "haploblock_sv_classification.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(ax.figure)

    # Select haploblocks with the largest population-specific SV counts.
    ranked_blocks = (
        block_counts.sort_values(
            ["population_specific", "total_svs"],
            ascending=[False, False],
        )
        .head(top_haploblocks)
        .index
    )

    specific_data = unique_sv_block[
        (unique_sv_block["sv_class"] == "population_specific")
        & unique_sv_block["haploblock_id"].isin(ranked_blocks)
    ]

    # Plot 2: population-specific SVs by haploblock and subpopulation.
    population_counts = (
        specific_data.dropna(subset=["specific_to_population"])
        .drop_duplicates(["sv_record_id", "haploblock_id"])
        .groupby(
            ["haploblock_id", "specific_to_population"],
            as_index=False,
        )
        .size()
        .pivot(
            index="haploblock_id",
            columns="specific_to_population",
            values="size",
        )
        .fillna(0)
    )

    if not population_counts.empty:
        population_counts = population_counts.loc[
            [block for block in ranked_blocks if block in population_counts.index]
        ]

        ax = population_counts.plot(
            kind="bar",
            stacked=True,
            figsize=(18, 8),
            colormap="tab20",
            width=0.9,
        )

        ax.set_title(
            f"{chrom}: population-specific SVs in top haploblocks"
        )
        ax.set_xlabel("Haploblock")
        ax.set_ylabel("Unique population-specific SV count")
        ax.tick_params(axis="x", labelrotation=90)
        ax.legend(
            title="Subpopulation",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
        )
        ax.figure.tight_layout()
        ax.figure.savefig(
            plot_dir / "population_specific_by_haploblock.png",
            dpi=200,
            bbox_inches="tight",
        )
        plt.close(ax.figure)

    # Plot 3: population-specific fraction by genomic position.
    fraction_data = block_counts.copy()
    block_coordinates = (
        unique_sv_block[
            ["haploblock_id", "block_start", "block_end"]
        ]
        .drop_duplicates("haploblock_id")
        .set_index("haploblock_id")
    )

    fraction_data = fraction_data.join(block_coordinates)
    fraction_data["population_specific_fraction"] = (
        fraction_data["population_specific"]
        / fraction_data["total_svs"].replace(0, np.nan)
    )
    fraction_data["midpoint_mb"] = (
        (fraction_data["block_start"] + fraction_data["block_end"]) / 2
    ) / 1_000_000

    ax = fraction_data.plot(
        x="midpoint_mb",
        y="population_specific_fraction",
        kind="scatter",
        figsize=(12, 6),
        s=np.maximum(fraction_data["total_svs"] * 2, 20),
        color="#E45756",
        alpha=0.75,
    )

    ax.set_title(
        f"{chrom}: population-specific SV fraction by haploblock position"
    )
    ax.set_xlabel("Haploblock midpoint (Mb)")
    ax.set_ylabel("Population-specific SVs / total SVs")
    ax.set_ylim(0, 1)
    ax.figure.tight_layout()
    ax.figure.savefig(
        plot_dir / "population_specific_fraction.png",
        dpi=200,
    )
    plt.close(ax.figure)

    block_counts.reset_index().to_csv(
        plot_dir / "haploblock_plot_data.tsv",
        sep="\t",
        index=False,
    )

    if not population_counts.empty:
        population_counts.reset_index().to_csv(
            plot_dir / "population_haploblock_plot_data.tsv",
            sep="\t",
            index=False,
        )

    log.info("Wrote Stage 4 plots for %s to %s", chrom, plot_dir)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, default=Path("stage1_output/config.yaml"),
                        help="Stage 1 config.yaml")
    parser.add_argument("--out-dir", type=Path, default=Path("stage4_output"))
    parser.add_argument("--af-threshold", type=float, default=0.05,
                        help="Population AF at or above which an SV counts as present")
    parser.add_argument("--absent-af-threshold", type=float, default=0.01,
                        help="AF below which an SV counts as absent from a population")
    parser.add_argument("--min-samples-per-population", type=int, default=2,
                        help="Called samples needed for a population to count as having data")
    parser.add_argument("--plots", action="store_true", help="Generate chromosome-specific plots")
    parser.add_argument("--plot-dir", type=Path, default=None,
                        help="Directory for plots; defaults to <out-dir>/stage4_plots")
    parser.add_argument("--top-haploblocks", type=int, default=20,
                        help="Number of haploblocks shown in ranked plots")
    parser.add_argument("--test", action="store_true",
                        help="Run Stages 1 and 4 on a tiny synthetic dataset and exit")
    args = parser.parse_args(argv)
    if not 0 <= args.absent_af_threshold <= args.af_threshold <= 1:
        parser.error("require 0 <= --absent-af-threshold <= --af-threshold <= 1")
    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    if args.test:
        from synthetic_data import self_test

        self_test(4)
        return

    config = load_config(args.config, ["sample_metadata", "sv_genotypes", "sv_block_summary"])
    metadata = read_sample_metadata(config)

    population_tables = []
    classification_tables = []
    for chrom, path in config["paths"]["sv_genotypes"].items():
        sv = read_tsv(path, METADATA_COLUMNS)
        by_population, classifications = classify(
            sv,
            metadata,
            args.af_threshold,
            args.absent_af_threshold,
            args.min_samples_per_population,
        )
        population_tables.append(by_population)
        classification_tables.append(classifications)
        log.info("%s: classified %d SVs", chrom, len(sv))

    by_population = pd.concat(population_tables, ignore_index=True)
    classifications = pd.concat(classification_tables, ignore_index=True)
    if classifications.empty:
        raise ValueError("No SV records found in the Stage 1 genotype tables")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    by_population_path = args.out_dir / "sv_af_classification.tsv"
    classification_path = args.out_dir / "sv_classification.tsv"
    by_population.to_csv(by_population_path, sep="\t", index=False)
    classifications.to_csv(classification_path, sep="\t", index=False)

    classification_haploblocks = classifications.merge(
        load_haploblock_assignments(config), on="sv_record_id", how="left"
    )
    classification_haploblocks_path = args.out_dir / "sv_classification_haploblocks.tsv"
    classification_haploblocks.to_csv(classification_haploblocks_path, sep="\t", index=False)
    write_informative_summary(classification_haploblocks, args.out_dir)

    if args.plots:
        plot_dir = args.plot_dir or args.out_dir / "stage4_plots"
        for chrom in config["paths"]["sv_genotypes"]:
            write_stage4_plots(classification_haploblocks, plot_dir, chrom, args.top_haploblocks)

    class_counts = classifications["sv_class"].value_counts().to_dict()
    log.info("Classified %d SVs: %s", len(classifications), class_counts)
    log.info("Other fraction: %.3f", class_counts.get("other", 0) / len(classifications))

    write_stage_config(
        config,
        args.out_dir,
        paths={
            "sv_af_classification": str(by_population_path.resolve()),
            "sv_classification": str(classification_path.resolve()),
            "sv_classification_haploblocks": str(classification_haploblocks_path.resolve()),
        },
        thresholds={
            "af_common_threshold": args.af_threshold,
            "af_absent_threshold": args.absent_af_threshold,
            "min_samples_per_population": args.min_samples_per_population,
        },
    )


if __name__ == "__main__":
    run_cli(main)
