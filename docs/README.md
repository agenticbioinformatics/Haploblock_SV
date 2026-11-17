# Methods documentation

Installation and run commands are in the top-level [README](../README.md). The documents here are reference material: inputs, data contracts, calculations, parameters, output fields and limitations of each stage.

| Stage | Script | Methods document | Input config |
|---|---|---|---|
| 0 | `stage0_merge_vcfs.py` | [Cohort SV merging](stage0_cohort_sv_merging.md) | — (manifest of single-sample VCFs) |
| 1 | `stage1_cluster_aware.py` | [Cluster-aware preprocessing](stage1_cluster_aware_preprocessing.md) | — (cohort SV VCF) |
| 2 | `stage2_boundary_classification.py` | [Boundary classification](stage2_boundary_classification.md) | Stage 1 |
| 4 | `stage4_classify_af.py` | [Population frequency classification](stage4_population_frequency.md) | Stage 1 |
| 5 | `stage5_type_enrichment.py` | [SV-type enrichment](stage5_sv_type_enrichment.md) | Stage 1 |
| 6 | `stage6_cluster_association.py` | [Population-conditioned association](stage6_cluster_association.md) | Stage 4 |
| 7 | `stage7_hash_representation.py` | [Hash representation audit](stage7_hash_representation.md) | Stage 6 |
| 8 | `stage8_candidate_annotation.py` | [Candidate annotation](stage8_candidate_annotation.md) | Stage 7 |
| 9 | `stage9_report.py` | [Reporting](stage9_reporting.md) | Stage 5 and Stage 8 |

Stage numbers match the original project plan; Stage 3 (boundary enrichment) was dropped. Each stage copies its input `config.yaml`, adds its own output paths and thresholds, and writes the result to its output directory, so a later stage only needs the config of the stage before it.

The pipeline overview diagram is [flowchart_haploblock_SV.png](flowchart_haploblock_SV.png). Hackathon planning notes are kept in [PROMPTS.md](PROMPTS.md) and [CHANGELOG.md](CHANGELOG.md).

## The question the pipeline answers

When can a local, SNV-derived haploblock cluster stand in for direct SV genotyping, and when does the SV add information that the cluster misses? If an SV is consistently carried by one haploblock cluster, the cluster may be a useful proxy. If carriers and non-carriers are mixed within the same cluster, direct SV genotyping adds resolution. Population labels are used throughout so that broad ancestry differences are not mistaken for a local haploblock–SV relationship.

## What the results mean

Two complementary views describe whether a haploblock represents an SV well:

| View | Question | A result that supports the haploblock | A result that suggests the SV adds information |
|---|---|---|---|
| SV-centric | Where are carriers of this SV found? | Carrier evidence is concentrated in one cluster. | Carriers occur in several clusters. |
| Cluster-centric | Are members of this cluster alike for SV carriage? | Nearly all callable members are carriers, or nearly all are non-carriers. | The same sufficiently represented cluster contains both carriers and non-carriers. |

These are association results, not proof of causality or functional effect. A gene overlap identifies a locus worth examining; it does not establish that the SV changes that gene's function.

## Data caveats

- The VCF phase is not assumed to correspond to haploblock `hap0`/`hap1` labels. An unphased heterozygote is assigned to a cluster only when the surrounding evidence supports it; otherwise it stays ambiguous. Haploblock haplotype labels are not maternal and paternal labels.
- The 1000 Genomes ONT callset used for the Results contains only DEL and INS (kanpig regenotyping requires resolved sequences); this is a limitation of that callset, not of the method.
- Rare calls, missing genotypes, SV representation errors in repetitive sequence, sequencing batches and uneven population sampling can all create apparent associations. Treat candidate tags and within-cluster splits as prioritization results until they have adequate support and replication.
