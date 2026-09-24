CIGAR validation cost
================

This focused comparison measures the old CIGAR implementation and the shared
checked decoder in separate R processes on the same host and inputs. It measures
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
| revision         | 70d2d749c6e5b5f1df97117178d544efa9783f71                         | e36af9dd7a855f854b12979e3a2fe25812b06908                         |
| src_tree         | ac3e3678585711876f0fd2df2ddd0226722a9acb                         | 5981a12439da7e8b87cc3e9efac96fd3ecd29668                         |
| source_state     | clean                                                            | clean                                                            |
| extension_sha256 | 3a89ecce87ab9fa1b4604122fd3fdc288c99c441f1eaa23b34a2aecbc0265308 | ac66633907f18d5dcdb2a0f39a4794cbb50c2f69834ac5dcb4c8bc2358420e68 |
| source_dir       | /tmp/duckhts-pr243-codex.klkV1H/baseline                         | /tmp/duckhts-pr243-codex.klkV1H/repo                             |
| run_utc          | 2026-09-23 18:26:01 UTC                                          | 2026-09-23 18:26:40 UTC                                          |
| host             | Ubuntu-2404-noble-amd64-base                                     | Ubuntu-2404-noble-amd64-base                                     |
| system           | Linux 6.8.0-78-generic x86_64                                    | Linux 6.8.0-78-generic x86_64                                    |
| cpu              | 13th Gen Intel(R) Core(TM) i5-13500                              | 13th Gen Intel(R) Core(TM) i5-13500                              |
| affinity         | 0                                                                | 0                                                                |
| threads          | 1                                                                | 1                                                                |
| r_version        | R version 4.6.0 (2026-04-24)                                     | R version 4.6.0 (2026-04-24)                                     |
| duckdb_version   | v1.5.5                                                           | v1.5.5                                                           |
| htslib_version   | 1.24                                                             | 1.24                                                             |
| driver_sha256    | 0209917eb42ba926986b51e6c118a4245f1e2d097fe968194d915285d2f174c9 | 0209917eb42ba926986b51e6c118a4245f1e2d097fe968194d915285d2f174c9 |
| size             | full                                                             | full                                                             |

| Property         | Value                                                                      |
|:-----------------|:---------------------------------------------------------------------------|
| ont_artifact     | ont_ecoli_k12_bam                                                          |
| ont_sha256       | b9bf633c99fd27503f3413c9c3208dc8d39a05c675ff15395972b981fa2225b2           |
| ont_release      | ENA_PRJEB86481_ERR14686255                                                 |
| ont_sources      | artifact:ont_ecoli_k12_reference_fna;artifact:ont_ecoli_k12_reads_fastq_gz |
| ont_transform    | minimap2_map_ont;samtools_sort;samtools_index                              |
| ont_aligner      | minimap2 2.26-r1175                                                        |
| ont_sorter       | samtools 1.24                                                              |
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
| long_cigar | cigar_aligned_blocks       | packed         | 8192           | 8192            | 786432       | 5           |      0.009 |     0.009 |   0.009 |  0.009 |                1.0000 |              1.0000 |             1.0000 |
| long_cigar | cigar_aligned_blocks       | text           | 8192           | 8192            | 786432       | 5           |      0.011 |     0.010 |   0.010 |  0.011 |                0.9091 |              0.9091 |             1.0000 |
| long_cigar | cigar_aligned_query_length | packed         | 8192           | 8192            | 9565312      | 5           |      0.003 |     0.004 |   0.004 |  0.004 |                1.3333 |              1.3333 |             1.3333 |
| long_cigar | cigar_aligned_query_length | text           | 8192           | 8192            | 9565312      | 5           |      0.005 |     0.006 |   0.006 |  0.006 |                1.2000 |              1.2000 |             1.2000 |
| long_cigar | cigar_has_hard_clip        | packed         | 8192           | 8192            | 8192         | 5           |      0.003 |     0.004 |   0.004 |  0.004 |                1.3333 |              1.3333 |             1.3333 |
| long_cigar | cigar_has_hard_clip        | text           | 8192           | 8192            | 8192         | 5           |      0.006 |     0.006 |   0.006 |  0.006 |                1.0000 |              1.0000 |             1.0000 |
| long_cigar | cigar_has_soft_clip        | packed         | 8192           | 8192            | 8192         | 5           |      0.003 |     0.004 |   0.004 |  0.004 |                1.3333 |              1.3333 |             1.3333 |
| long_cigar | cigar_has_soft_clip        | text           | 8192           | 8192            | 8192         | 5           |      0.006 |     0.006 |   0.006 |  0.006 |                1.0000 |              1.0000 |             1.0000 |
| long_cigar | cigar_left_soft_clip       | packed         | 8192           | 8192            | 16381        | 5           |      0.003 |     0.004 |   0.004 |  0.004 |                1.3333 |              1.3333 |             1.3333 |
| long_cigar | cigar_left_soft_clip       | text           | 8192           | 8192            | 16381        | 5           |      0.005 |     0.006 |   0.006 |  0.006 |                1.2000 |              1.2000 |             1.2000 |
| long_cigar | cigar_query_length         | packed         | 8192           | 8192            | 9856124      | 5           |      0.003 |     0.004 |   0.004 |  0.004 |                1.3333 |              1.3333 |             1.3333 |
| long_cigar | cigar_query_length         | text           | 8192           | 8192            | 9856124      | 5           |      0.006 |     0.006 |   0.006 |  0.006 |                1.0000 |              1.0000 |             1.0000 |
| long_cigar | cigar_reference_length     | packed         | 8192           | 8192            | 10876032     | 5           |      0.003 |     0.004 |   0.004 |  0.005 |                1.3333 |              1.3333 |             1.6667 |
| long_cigar | cigar_reference_length     | text           | 8192           | 8192            | 10876032     | 5           |      0.006 |     0.006 |   0.006 |  0.006 |                1.0000 |              1.0000 |             1.0000 |
| long_cigar | cigar_right_soft_clip      | packed         | 8192           | 8192            | 12287        | 5           |      0.003 |     0.004 |   0.004 |  0.004 |                1.3333 |              1.3333 |             1.3333 |
| long_cigar | cigar_right_soft_clip      | text           | 8192           | 8192            | 12287        | 5           |      0.005 |     0.006 |   0.006 |  0.006 |                1.2000 |              1.2000 |             1.2000 |
| long_cigar | has_op_M                   | packed         | 8192           | 8192            | 8192         | 5           |      0.001 |     0.005 |   0.004 |  0.004 |                5.0000 |              4.0000 |             4.0000 |
| long_cigar | has_op_M                   | text           | 8192           | 8192            | 8192         | 5           |      0.001 |     0.006 |   0.006 |  0.006 |                6.0000 |              6.0000 |             6.0000 |
| long_cigar | has_op_P                   | packed         | 8192           | 8192            | 0            | 5           |      0.002 |     0.004 |   0.004 |  0.004 |                2.0000 |              2.0000 |             2.0000 |
| long_cigar | has_op_P                   | text           | 8192           | 8192            | 0            | 5           |      0.003 |     0.006 |   0.006 |  0.006 |                2.0000 |              2.0000 |             2.0000 |
| ont        | cigar_aligned_blocks       | packed         | 27377          | 25689           | 1408917      | 5           |      0.016 |     0.019 |   0.020 |  0.020 |                1.1875 |              1.2500 |             1.2500 |
| ont        | cigar_aligned_blocks       | text           | 27377          | 25689           | 1408917      | 5           |      0.026 |     0.030 |   0.029 |  0.030 |                1.1538 |              1.1154 |             1.1538 |
| ont        | cigar_aligned_query_length | packed         | 27377          | 25689           | 57342744     | 5           |      0.011 |     0.012 |   0.012 |  0.012 |                1.0909 |              1.0909 |             1.0909 |
| ont        | cigar_aligned_query_length | text           | 27377          | 25689           | 57342744     | 5           |      0.020 |     0.021 |   0.021 |  0.020 |                1.0500 |              1.0500 |             1.0000 |
| ont        | cigar_has_hard_clip        | packed         | 27377          | 25689           | 121          | 5           |      0.011 |     0.011 |   0.012 |  0.012 |                1.0000 |              1.0909 |             1.0909 |
| ont        | cigar_has_hard_clip        | text           | 27377          | 25689           | 121          | 5           |      0.020 |     0.020 |   0.021 |  0.020 |                1.0000 |              1.0500 |             1.0000 |
| ont        | cigar_has_soft_clip        | packed         | 27377          | 25689           | 20131        | 5           |      0.011 |     0.012 |   0.012 |  0.012 |                1.0909 |              1.0909 |             1.0909 |
| ont        | cigar_has_soft_clip        | text           | 27377          | 25689           | 20131        | 5           |      0.020 |     0.021 |   0.020 |  0.020 |                1.0500 |              1.0000 |             1.0000 |
| ont        | cigar_left_soft_clip       | packed         | 27377          | 25689           | 937452       | 5           |      0.011 |     0.012 |   0.012 |  0.013 |                1.0909 |              1.0909 |             1.1818 |
| ont        | cigar_left_soft_clip       | text           | 27377          | 25689           | 937452       | 5           |      0.020 |     0.021 |   0.022 |  0.020 |                1.0500 |              1.1000 |             1.0000 |
| ont        | cigar_query_length         | packed         | 27377          | 25689           | 60050517     | 5           |      0.011 |     0.012 |   0.012 |  0.012 |                1.0909 |              1.0909 |             1.0909 |
| ont        | cigar_query_length         | text           | 27377          | 25689           | 60050517     | 5           |      0.020 |     0.020 |   0.021 |  0.021 |                1.0000 |              1.0500 |             1.0500 |
| ont        | cigar_reference_length     | packed         | 27377          | 25689           | 58809044     | 5           |      0.011 |     0.012 |   0.012 |  0.012 |                1.0909 |              1.0909 |             1.0909 |
| ont        | cigar_reference_length     | text           | 27377          | 25689           | 58809044     | 5           |      0.021 |     0.020 |   0.021 |  0.020 |                0.9524 |              1.0000 |             0.9524 |
| ont        | cigar_right_soft_clip      | packed         | 27377          | 25689           | 958121       | 5           |      0.011 |     0.012 |   0.012 |  0.013 |                1.0909 |              1.0909 |             1.1818 |
| ont        | cigar_right_soft_clip      | text           | 27377          | 25689           | 958121       | 5           |      0.020 |     0.020 |   0.022 |  0.021 |                1.0000 |              1.1000 |             1.0500 |
| ont        | has_op_M                   | packed         | 27377          | 27377           | 25689        | 5           |      0.003 |     0.012 |   0.012 |  0.012 |                4.0000 |              4.0000 |             4.0000 |
| ont        | has_op_M                   | text           | 27377          | 27377           | 25689        | 5           |      0.003 |     0.021 |   0.022 |  0.021 |                7.0000 |              7.3333 |             7.0000 |
| ont        | has_op_P                   | packed         | 27377          | 27377           | 0            | 5           |      0.004 |     0.012 |   0.012 |  0.012 |                3.0000 |              3.0000 |             3.0000 |
| ont        | has_op_P                   | text           | 27377          | 27377           | 0            | 5           |      0.011 |     0.020 |   0.021 |  0.021 |                1.8182 |              1.9091 |             1.9091 |

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

The declared outputs are `data/cigar_validation_{baseline,candidate}_{metadata,inputs,results,timings}.csv`
(or the `cigar_validation_smoke` prefix), plus `*_mismatches.csv` on an oracle
failure. They retain per-repeat elapsed times, minimum/median/maximum times,
complete output counts and both identity-linked fingerprints, input
denominators, and source/binary metadata. The report is rendered only after
both runs pass their independent oracle checks.
