# Stage 0: cohort SV merging

Run instructions are in the top-level [README](../README.md#usage); this page documents the method and output contract.

## Purpose

Stage 0 turns single-sample long-read SV calls into a cohort VCF in which records are intended to represent comparable biological alleles. This is the prerequisite for interpreting an SV on multiple haploblock backgrounds: without record reconciliation, that pattern could reflect different calls at a repeat rather than one allele.

The results in the README used a callset that was merged before the hackathon (`1kgp_ont_cohort.postfilter.full.vcf.gz`, produced with Sniffles, `bcftools merge`, `truvari collapse` and kanpig regenotyping). Kanpig requires resolved variant sequences, so that callset contains DEL and INS records but not INV or BND records. Stage 0 reproduces the merge-and-collapse part for a new cohort; it does not run kanpig.

## Inputs

- A manifest: a text file with one local single-sample VCF/VCF.GZ path per line. Blank lines and `#` comments are ignored.
- At least two VCFs, each with exactly one sample. Sample names must be unique unless `--allow-duplicate-sample-names` is given (passes `--force-samples` to `bcftools merge`).
- All VCFs must declare the same contig names and lengths (i.e. the same reference build) unless `--force-mismatched-references` is given.
- Optionally, an indexed (`samtools faidx`) reference FASTA via `--reference`, passed to `truvari collapse` so it can sequence-resolve symbolic alleles.
- `bcftools` and `truvari` on `PATH` (or given with `--bcftools` / `--truvari`).

## Method

1. Validate the manifest, sample names and reference contigs before doing any work.
2. `bcftools sort` every input.
3. Give every record a cohort-unique ID (`s<sample>_<CHROM>_<POS>_<REF>_<ALT>`) with `bcftools annotate`, then index.
4. `bcftools merge -m id`. Because no two records share an ID, bcftools only unions rows and pastes sample columns. With its default matching, bcftools would fuse different symbolic SVs that share `POS` (it ignores `END`/`SVLEN`) and keep only one extent; see the [truvari wiki](https://github.com/ACEnglish/truvari/wiki/collapse#symbolic-variants).
5. `truvari collapse` performs all SV matching. Tune its thresholds with repeated `--truvari-arg`, e.g. `--truvari-arg=--refdist --truvari-arg=500`.
6. Sort **both** truvari outputs before tabix indexing, because truvari does not guarantee coordinate order.

Stage 0 performs representation-level merging only. It does not apply analysis-specific size filters or silently remove `IMPRECISE` records; downstream stages decide which records to analyze.

## Outputs

```text
stage0_output/
├── staged/                              per-sample sorted and re-IDed inputs
├── merged.exact.vcf.gz(.tbi)            every post-merge record
├── merged.collapsed.vcf.gz(.tbi)        records kept by truvari -> Stage 1 --vcf
├── merged.collapsed_removed.vcf.gz(.tbi) records collapsed into another record
└── run.json                             inputs, options, staged files and outputs
```

A non-zero exit status means the run is incomplete and its outputs should not be used without inspection.

## Output contract expected by Stage 1

| VCF element | Requirement | Downstream use |
|---|---|---|
| Record identifier | Any value; need not be unique | Stage 1 creates its own unique `sv_record_id` and keeps this as `sv_id`. |
| Chromosome and coordinates | GRCh38, with an SV type and usable span | Defines haploblock overlap and gene context. |
| `SVTYPE` and `SVLEN`/`END` | Present where applicable | Used for type-specific analyses and annotation. |
| `FILTER` and `IMPRECISE` | Preserved rather than silently removed | Retained as call-quality context. |
| Per-sample `GT` | Present and interpretable | Used to calculate dosage, frequency and cluster association. |

## Major limitations

Calls in low-complexity or repetitive sequence may remain difficult to reconcile. An apparent multi-cluster SV should therefore be checked for local representation and read-level support before it is interpreted as recurrent mutation, recombination or an old allele.
