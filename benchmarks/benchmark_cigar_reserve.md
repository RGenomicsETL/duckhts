CIGAR validation cost
================

This focused comparison measures the aligned-block reservation change against origin/develop in separate R processes on the same host and inputs. It measures
all seven metric helpers individually, `has_op('M')`, `has_op('P')`, and aligned
blocks, for both text and packed CIGAR. Candidate calls use the default overload
and explicit strict `FALSE` and `TRUE`. Every row is retained, including missing
CIGARs; any oracle or cross-extension output disagreement stops the run or render.

The nearest retained baseline is
[`benchmark_cigar_aligned_blocks.md`](benchmark_cigar_aligned_blocks.md), source
`1139f23ddde1`, with 27,377 ONT records / 2,817,619 operations / 1,408,917 blocks
and 9,569,553 HG00188 chr22 records / 10,964,914 operations / 9,971,254 blocks.
That report used minimap2 2.28-r1209 and samtools 1.21 for ONT. This focused report
does **not** establish an identical-input comparison with that older ONT
derivation. It does **not** rerun Riker; the approximately 17 GB acquisition and
short-read measurement remain outstanding. The full existing report is retained.

## Reproduction and scope

Run from the repository root with DBI, duckdb, digest, and duckhtsbench available
in the R library. ONT must already be staged through the duckhtsbench registry;
the driver only resolves cached paths and performs no download or staging.
Use one identical `taskset` CPU mask for both commands. Supply the full source
revision and `HEAD:src` IDs for each build, and `uncommitted` if the extension
was built from modified sources. Those IDs describe the base revision when the
source is uncommitted; the loaded binary’s SHA256 identifies the measured binary.
The driver records caller-supplied source metadata without certifying the build.
The CSV prefix for this rendering is benchmarks/data/cigar_reserve_257;
use it in place of `benchmarks/data/cigar_validation` below.

``` sh
# Set these to the local build paths and full Git object IDs for the measured builds.
taskset -c "$CPU" Rscript scripts/benchmark_cigar_validation.R baseline \
  "$BASELINE_EXTENSION" "$BASELINE_SOURCE" "$BASELINE_REVISION" "$BASELINE_SRC_TREE" \
  clean benchmarks/data/cigar_validation full 5
taskset -c "$CPU" Rscript scripts/benchmark_cigar_validation.R candidate \
  "$CANDIDATE_EXTENSION" "$CANDIDATE_SOURCE" "$CANDIDATE_REVISION" "$CANDIDATE_SRC_TREE" \
  "$CANDIDATE_SOURCE_STATE" benchmarks/data/cigar_validation full 5
Rscript -e 'rmarkdown::render("benchmarks/benchmark_cigar_validation.Rmd")'
```

Each invocation loads exactly one extension. For a smoke check use `smoke 2`
and prefix `benchmarks/data/cigar_validation_smoke`: this uses the first 256
sequential ONT records and 4,096 synthetic records with 50 ops each. Full runs
use every sequential ONT record and 8,192 synthetic records with 194 ops each.
Synthetic records have alternating H/S endpoints and repeated M/I/D/N/=/X ops;
lengths and positions vary with record identity. M matches near the start and
P is absent, making the cost of checking the full suffix in `has_op` visible.

Inputs are materialized before timing. Text is rendered directly from the same
packed words using htslib’s `MIDNSHP=X` mapping. SQL columns keep native calls
from being constant-folded. The seven metric oracles use independent packed-op
sums and literal first/last-op rules. Block expectations use the prefix-slice
SQL geometry in the existing aligned-block benchmark. Every result is compared
row by row with `IS DISTINCT FROM`, including NULLs and all three block lists,
before one excluded warm-up and the timed repetitions. Disagreements retain
counterexamples in a declared `*_mismatches.csv` output and stop execution.

Timed queries consume the full helper output in aggregates, returning one row
to R. Timing includes SQL execution, output hashing and the DBI aggregate
transfer; it excludes BAM I/O, input construction, oracle construction,
extension loading, warm-up and result comparisons. This measures resident
projection cost, not a pure C kernel, I/O throughput or memory use. Each pass
must match the independent oracle’s exact counts, totals, XOR and sum of hashes
linked to physical row ordinal, QNAME, FLAG, RNAME and POS. Rendering verifies
these results across processes and strict modes, together with the input
fingerprint and BAM SHA256. No disagreement is filtered or tolerated.

## Environment and inputs

| Property         | Baseline                                                         | Candidate                                                        |
|:-----------------|:-----------------------------------------------------------------|:-----------------------------------------------------------------|
| revision         | 81edf9bff85ab2daabeccd1daaf13e4bfa39ee47                         | 538c1a7b79a57c695c4ea3768c26cb2c672b0f41                         |
| src_tree         | eaefad637812712f4f261e09a6803c8b65b60e19                         | a0b0888be6ddd6e12b489b6450d114f735856176                         |
| source_state     | clean                                                            | clean                                                            |
| extension_sha256 | e122aee4aec5d55130bd40f74c379a4e83f1b0e72600a00a95a3d09ebd886689 | 8d0745ea78c07d3ae686e54ba90db35c776afbfb620850c18d5d62bc1b6af7d0 |
| source_dir       | /tmp/duckhts-codex257-develop-baseline                           | /root/duckhts-codex257                                           |
| run_utc          | 2026-09-25 03:09:42 UTC                                          | 2026-09-25 03:10:53 UTC                                          |
| host             | Ubuntu-2404-noble-amd64-base                                     | Ubuntu-2404-noble-amd64-base                                     |
| system           | Linux 6.8.0-78-generic x86_64                                    | Linux 6.8.0-78-generic x86_64                                    |
| cpu              | 13th Gen Intel(R) Core(TM) i5-13500                              | 13th Gen Intel(R) Core(TM) i5-13500                              |
| affinity         | 16                                                               | 16                                                               |
| threads          | 1                                                                | 1                                                                |
| r_version        | R version 4.6.0 (2026-04-24)                                     | R version 4.6.0 (2026-04-24)                                     |
| duckdb_version   | v1.5.5                                                           | v1.5.5                                                           |
| htslib_version   | 1.24                                                             | 1.24                                                             |
| driver_sha256    | 0209917eb42ba926986b51e6c118a4245f1e2d097fe968194d915285d2f174c9 | 0209917eb42ba926986b51e6c118a4245f1e2d097fe968194d915285d2f174c9 |
| size             | full                                                             | full                                                             |

| Property         | Value                                                                      |
|:-----------------|:---------------------------------------------------------------------------|
| ont_artifact     | ont_ecoli_k12_bam                                                          |
| ont_sha256       | 8bda421bd0a523b2b76ab3d9e984b9ba665c5ec799482378ef02965317c019ef           |
| ont_release      | ENA_PRJEB86481_ERR14686255                                                 |
| ont_sources      | artifact:ont_ecoli_k12_reference_fna;artifact:ont_ecoli_k12_reads_fastq_gz |
| ont_transform    | minimap2_map_ont;samtools_sort;samtools_index                              |
| ont_aligner      | minimap2 2.26-r1175                                                        |
| ont_sorter       | samtools 1.23                                                              |
| synthetic_rows   | 8192                                                                       |
| synthetic_cycles | 32                                                                         |
| riker            | Not staged; not measured (approximately 17 GB acquisition)                 |

| input      | size | records | operations | blocks  | null_cigars | empty_cigars | max_ops_per_record | input_fingerprint |
|:-----------|:-----|:--------|:-----------|:--------|:------------|:-------------|:-------------------|:------------------|
| ont        | full | 27377   | 2817619    | 1408917 | 0           | 1688         | 4602               | DA80E4E560F37079  |
| long_cigar | full | 8192    | 1589248    | 786432  | 0           | 0            | 194                | C8B966C931CFE7E1  |

Operation and block counts cover all input records. Missing CIGAR counts are
reported separately and remain in every helper’s output denominator. The
synthetic SQL is literal workload construction in the driver and requires no
artifact registration. ONT identifies ENA ERR14686255 on RefSeq GCF_000005845.2;
the receipt records its actual aligner and sorter above.

## Outputs and timing

| input      | workload                   | representation | output_records | nonnull_outputs | output_total | repetitions | baseline_s | default_s | false_s | true_s | default_over_baseline | false_over_baseline | true_over_baseline |
|:-----------|:---------------------------|:---------------|:---------------|:----------------|:-------------|:------------|-----------:|----------:|--------:|-------:|----------------------:|--------------------:|-------------------:|
| long_cigar | cigar_aligned_blocks       | packed         | 8192           | 8192            | 786432       | 5           |      0.017 |     0.018 |   0.018 |  0.018 |                1.0588 |              1.0588 |             1.0588 |
| long_cigar | cigar_aligned_blocks       | text           | 8192           | 8192            | 786432       | 5           |      0.020 |     0.028 |   0.028 |  0.028 |                1.4000 |              1.4000 |             1.4000 |
| long_cigar | cigar_aligned_query_length | packed         | 8192           | 8192            | 9565312      | 5           |      0.010 |     0.007 |   0.008 |  0.008 |                0.7000 |              0.8000 |             0.8000 |
| long_cigar | cigar_aligned_query_length | text           | 8192           | 8192            | 9565312      | 5           |      0.011 |     0.011 |   0.011 |  0.011 |                1.0000 |              1.0000 |             1.0000 |
| long_cigar | cigar_has_hard_clip        | packed         | 8192           | 8192            | 8192         | 5           |      0.008 |     0.007 |   0.008 |  0.008 |                0.8750 |              1.0000 |             1.0000 |
| long_cigar | cigar_has_hard_clip        | text           | 8192           | 8192            | 8192         | 5           |      0.011 |     0.011 |   0.012 |  0.011 |                1.0000 |              1.0909 |             1.0000 |
| long_cigar | cigar_has_soft_clip        | packed         | 8192           | 8192            | 8192         | 5           |      0.009 |     0.007 |   0.008 |  0.008 |                0.7778 |              0.8889 |             0.8889 |
| long_cigar | cigar_has_soft_clip        | text           | 8192           | 8192            | 8192         | 5           |      0.012 |     0.011 |   0.011 |  0.011 |                0.9167 |              0.9167 |             0.9167 |
| long_cigar | cigar_left_soft_clip       | packed         | 8192           | 8192            | 16381        | 5           |      0.008 |     0.007 |   0.009 |  0.008 |                0.8750 |              1.1250 |             1.0000 |
| long_cigar | cigar_left_soft_clip       | text           | 8192           | 8192            | 16381        | 5           |      0.011 |     0.011 |   0.012 |  0.011 |                1.0000 |              1.0909 |             1.0000 |
| long_cigar | cigar_query_length         | packed         | 8192           | 8192            | 9856124      | 5           |      0.009 |     0.007 |   0.008 |  0.008 |                0.7778 |              0.8889 |             0.8889 |
| long_cigar | cigar_query_length         | text           | 8192           | 8192            | 9856124      | 5           |      0.011 |     0.011 |   0.012 |  0.011 |                1.0000 |              1.0909 |             1.0000 |
| long_cigar | cigar_reference_length     | packed         | 8192           | 8192            | 10876032     | 5           |      0.010 |     0.007 |   0.007 |  0.007 |                0.7000 |              0.7000 |             0.7000 |
| long_cigar | cigar_reference_length     | text           | 8192           | 8192            | 10876032     | 5           |      0.011 |     0.011 |   0.011 |  0.012 |                1.0000 |              1.0000 |             1.0909 |
| long_cigar | cigar_right_soft_clip      | packed         | 8192           | 8192            | 12287        | 5           |      0.009 |     0.007 |   0.008 |  0.008 |                0.7778 |              0.8889 |             0.8889 |
| long_cigar | cigar_right_soft_clip      | text           | 8192           | 8192            | 12287        | 5           |      0.011 |     0.011 |   0.011 |  0.011 |                1.0000 |              1.0000 |             1.0000 |
| long_cigar | has_op_M                   | packed         | 8192           | 8192            | 8192         | 5           |      0.008 |     0.007 |   0.007 |  0.007 |                0.8750 |              0.8750 |             0.8750 |
| long_cigar | has_op_M                   | text           | 8192           | 8192            | 8192         | 5           |      0.012 |     0.011 |   0.012 |  0.012 |                0.9167 |              1.0000 |             1.0000 |
| long_cigar | has_op_P                   | packed         | 8192           | 8192            | 0            | 5           |      0.008 |     0.007 |   0.007 |  0.007 |                0.8750 |              0.8750 |             0.8750 |
| long_cigar | has_op_P                   | text           | 8192           | 8192            | 0            | 5           |      0.012 |     0.011 |   0.011 |  0.012 |                0.9167 |              0.9167 |             1.0000 |
| ont        | cigar_aligned_blocks       | packed         | 27377          | 25689           | 1408917      | 5           |      0.037 |     0.037 |   0.037 |  0.037 |                1.0000 |              1.0000 |             1.0000 |
| ont        | cigar_aligned_blocks       | text           | 27377          | 25689           | 1408917      | 5           |      0.052 |     0.074 |   0.074 |  0.076 |                1.4231 |              1.4231 |             1.4615 |
| ont        | cigar_aligned_query_length | packed         | 27377          | 25689           | 57342744     | 5           |      0.019 |     0.019 |   0.019 |  0.019 |                1.0000 |              1.0000 |             1.0000 |
| ont        | cigar_aligned_query_length | text           | 27377          | 25689           | 57342744     | 5           |      0.031 |     0.031 |   0.031 |  0.031 |                1.0000 |              1.0000 |             1.0000 |
| ont        | cigar_has_hard_clip        | packed         | 27377          | 25689           | 121          | 5           |      0.019 |     0.019 |   0.019 |  0.019 |                1.0000 |              1.0000 |             1.0000 |
| ont        | cigar_has_hard_clip        | text           | 27377          | 25689           | 121          | 5           |      0.032 |     0.031 |   0.031 |  0.031 |                0.9688 |              0.9688 |             0.9688 |
| ont        | cigar_has_soft_clip        | packed         | 27377          | 25689           | 20131        | 5           |      0.019 |     0.019 |   0.019 |  0.019 |                1.0000 |              1.0000 |             1.0000 |
| ont        | cigar_has_soft_clip        | text           | 27377          | 25689           | 20131        | 5           |      0.032 |     0.031 |   0.031 |  0.031 |                0.9688 |              0.9688 |             0.9688 |
| ont        | cigar_left_soft_clip       | packed         | 27377          | 25689           | 937452       | 5           |      0.019 |     0.019 |   0.019 |  0.019 |                1.0000 |              1.0000 |             1.0000 |
| ont        | cigar_left_soft_clip       | text           | 27377          | 25689           | 937452       | 5           |      0.031 |     0.031 |   0.031 |  0.031 |                1.0000 |              1.0000 |             1.0000 |
| ont        | cigar_query_length         | packed         | 27377          | 25689           | 60050517     | 5           |      0.019 |     0.019 |   0.019 |  0.019 |                1.0000 |              1.0000 |             1.0000 |
| ont        | cigar_query_length         | text           | 27377          | 25689           | 60050517     | 5           |      0.031 |     0.036 |   0.032 |  0.031 |                1.1613 |              1.0323 |             1.0000 |
| ont        | cigar_reference_length     | packed         | 27377          | 25689           | 58809044     | 5           |      0.019 |     0.019 |   0.019 |  0.019 |                1.0000 |              1.0000 |             1.0000 |
| ont        | cigar_reference_length     | text           | 27377          | 25689           | 58809044     | 5           |      0.031 |     0.031 |   0.031 |  0.031 |                1.0000 |              1.0000 |             1.0000 |
| ont        | cigar_right_soft_clip      | packed         | 27377          | 25689           | 958121       | 5           |      0.019 |     0.019 |   0.019 |  0.019 |                1.0000 |              1.0000 |             1.0000 |
| ont        | cigar_right_soft_clip      | text           | 27377          | 25689           | 958121       | 5           |      0.031 |     0.032 |   0.031 |  0.031 |                1.0323 |              1.0000 |             1.0000 |
| ont        | has_op_M                   | packed         | 27377          | 27377           | 25689        | 5           |      0.019 |     0.018 |   0.018 |  0.018 |                0.9474 |              0.9474 |             0.9474 |
| ont        | has_op_M                   | text           | 27377          | 27377           | 25689        | 5           |      0.032 |     0.032 |   0.032 |  0.032 |                1.0000 |              1.0000 |             1.0000 |
| ont        | has_op_P                   | packed         | 27377          | 27377           | 0            | 5           |      0.019 |     0.018 |   0.018 |  0.019 |                0.9474 |              0.9474 |             1.0000 |
| ont        | has_op_P                   | text           | 27377          | 27377           | 0            | 5           |      0.032 |     0.032 |   0.032 |  0.032 |                1.0000 |              1.0000 |             1.0000 |

All times are medians in seconds; ratios above 1 indicate extra cost relative
to the baseline. Output totals are summed lengths for numeric helpers, true
counts for Boolean helpers, and emitted block counts for aligned blocks.
The non-NULL count retains missing-output behavior. Candidate strict modes
have identical outputs to the default and baseline in this workload. The
default, explicit `FALSE`, and strict `TRUE` timings must be read separately:
full suffix validation can cost more than an early presence match. Small
timings include timer resolution and scheduling noise; a zero baseline median
has no meaningful ratio. These measurements make no scaling or malformed-input
performance claim.

The declared outputs use the prefix above, followed by `_baseline_` or
`_candidate_` and `metadata`, `inputs`, `results` or `timings` with a `.csv`
suffix; an oracle failure retains `*_mismatches.csv`. They retain per-repeat
elapsed times, minimum/median/maximum times,
complete output counts and both identity-linked fingerprints, input
denominators, and source/binary metadata. The report is rendered only after
both runs pass their independent oracle checks.
