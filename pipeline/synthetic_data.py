"""Tiny synthetic dataset and the ``--test`` self-test shared by all stage scripts.

The dataset mimics the real inputs at toy scale: one chromosome, two
haploblocks with three hash clusters each, 36 samples in three populations,
DEL/INS records whose carriers are concentrated in one cluster, local cluster
files in data.haploblocks.org layout, population metadata and a small GTF.
Everything is written to a temporary directory and removed afterwards.
"""

from __future__ import annotations

import importlib
import logging
import shutil
import tempfile
from pathlib import Path

import numpy as np


log = logging.getLogger("synthetic_data")

CHROM = "chr1"
CHROM_LENGTH = 100_000
BLOCKS = [(0, 20_000), (20_000, 40_000)]
CLUSTERS = ["c0", "c1", "c2"]
POPULATIONS = {"POP_A": "SUP_X", "POP_B": "SUP_X", "POP_C": "SUP_Y"}
SAMPLES_PER_POPULATION = 12
SVS_PER_BLOCK = 10

STAGE_MODULES = {
    0: "stage0_merge_vcfs",
    1: "stage1_cluster_aware",
    2: "stage2_boundary_classification",
    4: "stage4_classify_af",
    5: "stage5_type_enrichment",
    6: "stage6_cluster_association",
    7: "stage7_hash_representation",
    8: "stage8_candidate_annotation",
    9: "stage9_report",
}
# Which stage's config.yaml each stage consumes (the documented pipeline order).
CONFIG_FROM = {2: [1], 4: [1], 5: [1], 6: [4], 7: [6], 8: [7], 9: [5, 8]}
# Settings that keep each stage fast on the toy data.
TEST_ARGS = {
    4: ["--plots"],
    6: ["--permutations", "200", "--initial-permutations", "50",
        "--refinement-permutations", "0", "--threads", "1"],
    7: ["--threads", "1"],
    8: ["--threads", "1"],
    9: ["--agent", "off"],
}


def _sample_ids() -> list[tuple[str, str]]:
    return [
        (f"{population}_{index:02d}", population)
        for population in POPULATIONS
        for index in range(1, SAMPLES_PER_POPULATION + 1)
    ]


def write_synthetic_inputs(directory: Path, seed: int = 7) -> dict[str, Path]:
    """Write a cohort VCF, cluster files, sample metadata and GTF; return their paths."""
    import pysam

    rng = np.random.default_rng(seed)
    directory.mkdir(parents=True, exist_ok=True)
    samples = _sample_ids()

    metadata = directory / "sample_metadata.tsv"
    metadata.write_text(
        "sample_id\tpopulation\tsuperpopulation\n"
        + "".join(f"{sample}\t{pop}\t{POPULATIONS[pop]}\n" for sample, pop in samples)
    )

    # Each population has a different cluster mix so population structure is present.
    population_weights = {"POP_A": [0.6, 0.3, 0.1], "POP_B": [0.3, 0.4, 0.3], "POP_C": [0.1, 0.3, 0.6]}
    cluster_root = directory / "clusters"
    (cluster_root / CHROM).mkdir(parents=True, exist_ok=True)
    haplotype_clusters: dict[tuple[int, int], np.ndarray] = {}
    for block_index, (start, end) in enumerate(BLOCKS):
        lines = []
        for sample_index, (sample, population) in enumerate(samples):
            for haplotype in (0, 1):
                cluster = rng.choice(len(CLUSTERS), p=population_weights[population])
                haplotype_clusters[(block_index, sample_index * 2 + haplotype)] = cluster
                lines.append(f"{CLUSTERS[cluster]}\t{sample}_{CHROM}_region_{start}-{end}_hap{haplotype}\n")
        (cluster_root / CHROM / f"{CHROM}_{start}-{end}_cluster.tsv").write_text("".join(lines))

    records = []
    for block_index, (start, end) in enumerate(BLOCKS):
        positions = np.linspace(start + 500, end - 1500, SVS_PER_BLOCK).astype(int)
        for sv_index, position in enumerate(positions):
            sv_type = "DEL" if sv_index % 2 == 0 else "INS"
            length = int(rng.integers(60, 800))
            carrier_cluster = sv_index % len(CLUSTERS)
            genotypes = []
            for sample_index in range(len(samples)):
                if rng.random() < 0.03:
                    genotypes.append("./.")
                    continue
                alleles = [
                    int(rng.random() < (0.9 if haplotype_clusters[(block_index, sample_index * 2 + h)]
                                        == carrier_cluster else 0.02))
                    for h in (0, 1)
                ]
                genotypes.append(f"{min(alleles)}/{max(alleles)}")
            records.append((position, sv_type, length, genotypes))
    # One deletion spanning the block boundary, for Stage 2.
    records.append((BLOCKS[0][1] - 100, "DEL", 300, ["0/1" if i % 5 == 0 else "0/0" for i in range(len(samples))]))
    records.sort()

    vcf_text = directory / "cohort.vcf"
    header = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={CHROM},length={CHROM_LENGTH}>",
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(s for s, _ in samples),
    ]
    body = []
    for number, (position, sv_type, length, genotypes) in enumerate(records, start=1):
        end = position + length if sv_type == "DEL" else position + 1
        svlen = -length if sv_type == "DEL" else length
        body.append(
            f"{CHROM}\t{position}\tsyn_sv{number}\tN\t<{sv_type}>\t.\tPASS\t"
            f"SVTYPE={sv_type};SVLEN={svlen};END={end}\tGT\t" + "\t".join(genotypes)
        )
    vcf_text.write_text("\n".join(header + body) + "\n")
    vcf = Path(pysam.tabix_index(str(vcf_text), preset="vcf", force=True))

    gtf = directory / "genes.gtf"
    gtf_rows = []
    for gene_index, (gene_start, gene_end) in enumerate([(1_000, 9_000), (22_000, 30_000)], start=1):
        attributes = f'gene_id "GENE{gene_index}"; gene_name "GENE{gene_index}"; gene_biotype "protein_coding";'
        gtf_rows.append(f"{CHROM}\tsynthetic\tgene\t{gene_start}\t{gene_end}\t.\t+\t.\t{attributes}")
        for exon_start in range(gene_start, gene_end, 2_000):
            gtf_rows.append(f"{CHROM}\tsynthetic\texon\t{exon_start}\t{exon_start + 300}\t.\t+\t.\t{attributes}")
    gtf.write_text("\n".join(gtf_rows) + "\n")

    return {"vcf": vcf, "sample_metadata": metadata, "gtf": gtf, "cluster_root": cluster_root}


def write_single_sample_vcfs(directory: Path) -> Path:
    """Write two single-sample SV VCFs and a manifest listing them (Stage 0 input)."""
    directory.mkdir(parents=True, exist_ok=True)
    header = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={CHROM},length={CHROM_LENGTH}>",
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    manifest_lines = []
    # The second sample carries the same deletion with slightly shifted breakpoints,
    # which truvari collapse should merge.
    for sample, shift in (("SAMPLE_A", 0), ("SAMPLE_B", 5)):
        path = directory / f"{sample}.vcf"
        rows = [
            f"{CHROM}\t{1000 + shift}\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-500;END={1500 + shift}\tGT\t0/1",
            f"{CHROM}\t{5000 + shift}\t.\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=300;END={5001 + shift}\tGT\t1/1",
        ]
        path.write_text("\n".join(header + [
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"
        ] + rows) + "\n")
        manifest_lines.append(str(path))
    manifest = directory / "vcf_list.txt"
    manifest.write_text("\n".join(manifest_lines) + "\n")
    return manifest


def _run_stage(stage: int, workdir: Path, inputs: dict[str, Path], done: set[int]) -> None:
    if stage in done:
        return
    for upstream in CONFIG_FROM.get(stage, []):
        _run_stage(upstream, workdir, inputs, done)
    out_dir = workdir / f"stage{stage}_output"
    if stage == 1:
        argv = [
            "--vcf", str(inputs["vcf"]), "--sample-metadata", str(inputs["sample_metadata"]),
            "--gtf", str(inputs["gtf"]), "--cluster-root", str(inputs["cluster_root"]),
            "--chroms", CHROM, "--threads", "1",
        ]
    else:
        argv = []
        for upstream in CONFIG_FROM[stage]:
            argv += ["--config", str(workdir / f"stage{upstream}_output" / "config.yaml")]
        argv += TEST_ARGS.get(stage, [])
    log.info("Self-test: running Stage %d", stage)
    importlib.import_module(STAGE_MODULES[stage]).main(argv + ["--out-dir", str(out_dir)])
    done.add(stage)


def self_test(stage: int) -> None:
    """Run ``stage`` (and the upstream stages it needs) on synthetic data."""
    if stage not in STAGE_MODULES:
        raise ValueError(f"Unknown stage {stage}")
    workdir = Path(tempfile.mkdtemp(prefix=f"haploblock_sv_stage{stage}_test_"))
    try:
        if stage == 0:
            manifest = write_single_sample_vcfs(workdir / "inputs")
            importlib.import_module(STAGE_MODULES[0]).main(
                ["--vcf-list", str(manifest), "--out-dir", str(workdir / "stage0_output")]
            )
        else:
            inputs = write_synthetic_inputs(workdir / "inputs")
            _run_stage(stage, workdir, inputs, set())
        out_dir = workdir / f"stage{stage}_output"
        produced = sorted(str(path.relative_to(out_dir)) for path in out_dir.rglob("*") if path.is_file())
        if not produced:
            raise RuntimeError(f"Self-test for Stage {stage} produced no output files")
        empty = [name for name in produced if name.endswith(".tsv") and (out_dir / name).stat().st_size == 0]
        if empty:
            raise RuntimeError(f"Self-test for Stage {stage} wrote empty tables: {empty}")
        log.info("Self-test passed for Stage %d: %d output files (%s)",
                 stage, len(produced), ", ".join(produced[:8]) + (" ..." if len(produced) > 8 else ""))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
