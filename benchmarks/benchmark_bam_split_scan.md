read_bam: byte-balanced scan pieces and named tag columns
================

This report measures two `read_bam` changes against the previous reader on real
alignments:

- **Scan pieces.** An indexed full-file BGZF BAM scan used to claim one reference
  per DuckDB thread, so a file whose data sits on one reference scanned on one
  thread. A reference with more than 16 MiB of indexed data is now cut into
  byte-balanced pieces; a record belongs to the piece that holds its start
  (`design/better_scans.md`, `src/bam_partition.c`).
- **Named tag columns.** `auxiliary_tags := ['XS', 'XA']` returns one text column
  per tag instead of building the `AUXILIARY_TAGS` map for every record.

`bam_bin_counts` does not use scan pieces; it still claims one reference per thread.
This change moved its scan from bind to execution on DuckDB worker threads, so its
rows measure that thread policy (see the note after the results).

## Inputs

All inputs resolve through the `duckhtsbench` registry; rendering does not download.

| Input            | Registry id                                   | Shape                                                                                              |
|------------------|-----------------------------------------------|----------------------------------------------------------------------------------------------------|
| exome            | `bam_split_hg00403_exome_bam`                 | 1000 Genomes phase 3 exome of HG00403 (GRCh37 with decoys), all chromosomes                        |
| exome unmapped   | `bam_split_hg00403_unmapped_bam`              | the same sample’s reads without coordinates only                                                   |
| chr20 1x, 2x, 4x | `insert_size_hg00403_chr20_bam`, `_2x`, `_4x` | 30x HG00403 chr20, a full GRCh38 header with one populated reference; 2x and 4x repeat the records |
| ONT              | `ont_ecoli_k12_bam`                           | 27,377 MinION reads on one 4.6 Mb reference                                                        |

Stage them once, with network access only in this step:

``` sh
Rscript -e 'pkgload::load_all("r/duckhtsbench"); for (id in duckhts_bench_stage_plan("bam-split-scan")$id) duckhts_bench_fetch(id)'
Rscript -e 'duckhtsbench::duckhts_bench_stage_ont_ecoli()'
```

The chr20 BAMs come from the `insert-size` workload’s staging functions.

## Method

- Engines: **baseline**, the extension built from the source before this change
  (`BAM_SPLIT_BASELINE_EXTENSION`), and **candidate**, the current build
  (`DUCKHTS_EXTENSION`, default `build/release/duckhts.duckdb_extension`).
- Every timed query runs in a fresh DuckDB CLI process (`DUCKDB`, default `duckdb`)
  with `SET threads`, `SET preserve_insertion_order = false` and a 16 GiB
  memory limit. The time is the CLI `.timer` wall time of the measured statement
  only; peak RSS and CPU seconds (user plus system, all threads, including htslib
  decompression threads) are GNU `time` for the whole process.
- Three fresh processes per engine, input, workload and thread count; the table
  shows the median and range.
- Declared before measuring: a 900 s wall limit per process and a 20 GiB peak-RSS
  ceiling. A breach is a failed cell, not a reason to raise the limit.
- Workloads:
  - **scan**: `count(*)`, `sum(POS)` and `sum(hash(QNAME))` over the file.
  - **ctas**: `CREATE TEMP TABLE` of six core columns (not on the exome, whose
    rows would not fit the memory limit without spilling).
  - **tags map** and **tags named**: a fingerprint of two non-standard tags read
    from the map, then as named columns (candidate only).
  - **bin counts**: `bam_bin_counts(path, 100000)`, 100 kb bins.
- Correctness before timing: for each input and engine, a four-thread scan’s
  record fingerprint (count, distinct and summed `FILE_OFFSET`, summed hash of
  name, flag, reference, position and CIGAR) must equal the candidate’s
  sequential scan (`scan_mode := 'sequential'`, one thread), an independent
  single-stream path. The named-tag fingerprint must equal the map fingerprint.

## Correctness

| input          | engine         | equal |
|:---------------|:---------------|:------|
| exome          | baseline       | TRUE  |
| exome          | candidate      | TRUE  |
| exome          | candidate tags | TRUE  |
| exome unmapped | baseline       | TRUE  |
| exome unmapped | candidate      | TRUE  |
| chr20 1x       | baseline       | TRUE  |
| chr20 1x       | candidate      | TRUE  |
| chr20 1x       | candidate tags | TRUE  |
| chr20 2x       | baseline       | TRUE  |
| chr20 2x       | candidate      | TRUE  |
| chr20 2x       | candidate tags | TRUE  |
| chr20 4x       | baseline       | TRUE  |
| chr20 4x       | candidate      | TRUE  |
| chr20 4x       | candidate tags | TRUE  |
| ONT            | baseline       | TRUE  |
| ONT            | candidate      | TRUE  |
| ONT            | candidate tags | TRUE  |

## Results

Median of three fresh processes, with the range; seconds are the measured
statement only, and peak RSS is the largest of the three. A cell with any failed,
timed-out or over-ceiling process is marked failed and has no median or speedup.

| input          | workload   | threads | status.baseline | seconds.baseline | range.baseline | peak_rss_mib.baseline | cpu_seconds.baseline | status.candidate | seconds.candidate | range.candidate | peak_rss_mib.candidate | cpu_seconds.candidate | speedup |
|:---------------|:-----------|--------:|:----------------|-----------------:|:---------------|----------------------:|---------------------:|:-----------------|------------------:|:----------------|-----------------------:|----------------------:|--------:|
| exome          | bin_counts |       1 | PASS            |             2.77 | 2.66-3.03      |                   302 |                 43.9 | PASS             |             10.31 | 10.26-10.31     |                     61 |                  25.3 |    0.27 |
| exome          | bin_counts |       4 | PASS            |             2.94 | 2.75-2.95      |                   302 |                 48.6 | PASS             |              3.65 | 3.65-3.66       |                    105 |                  37.2 |    0.81 |
| exome          | scan       |       1 | PASS            |            10.30 | 10.28-10.30    |                    57 |                 27.2 | PASS             |             10.39 | 10.33-10.40     |                     56 |                  27.4 |    0.99 |
| exome          | scan       |       4 | PASS            |             3.84 | 3.22-3.84      |                   103 |                 41.6 | PASS             |              3.82 | 3.81-3.83       |                    102 |                  42.2 |    1.01 |
| exome          | tags_map   |       1 | PASS            |            49.95 | 49.52-50.07    |                    58 |                 72.9 | PASS             |             51.11 | 50.93-51.28     |                     57 |                  74.0 |    0.98 |
| exome          | tags_map   |       4 | PASS            |            15.78 | 14.34-15.83    |                   110 |                 95.3 | PASS             |             15.94 | 15.92-15.99     |                    108 |                  96.1 |    0.99 |
| exome          | tags_named |       1 | NA              |               NA | NA             |                    NA |                   NA | PASS             |             12.34 | 12.25-12.38     |                     56 |                  35.0 |      NA |
| exome          | tags_named |       4 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              4.71 | 4.31-4.79       |                    101 |                  52.9 |      NA |
| exome unmapped | bin_counts |       1 | PASS            |             0.00 | 0.00-0.00      |                    46 |                  0.1 | PASS             |              0.00 | 0.00-0.00       |                     44 |                   0.0 |     NaN |
| exome unmapped | bin_counts |       4 | PASS            |             0.00 | 0.00-0.00      |                    47 |                  0.1 | PASS             |              0.00 | 0.00-0.00       |                     45 |                   0.1 |     NaN |
| exome unmapped | ctas       |       1 | PASS            |             0.41 | 0.40-0.41      |                   194 |                  0.8 | PASS             |              0.42 | 0.41-0.42       |                    194 |                   0.9 |    0.98 |
| exome unmapped | ctas       |       4 | PASS            |             0.41 | 0.40-0.42      |                   194 |                  0.9 | PASS             |              0.40 | 0.40-0.41       |                    194 |                   0.8 |    1.02 |
| exome unmapped | scan       |       1 | PASS            |             0.17 | 0.17-0.17      |                    41 |                  0.6 | PASS             |              0.17 | 0.17-0.17       |                     41 |                   0.6 |    1.00 |
| exome unmapped | scan       |       4 | PASS            |             0.17 | 0.17-0.17      |                    42 |                  0.6 | PASS             |              0.17 | 0.17-0.17       |                     41 |                   0.6 |    1.00 |
| chr20 1x       | bin_counts |       1 | PASS            |             1.76 | 1.75-1.76      |                    79 |                  4.8 | PASS             |              1.76 | 1.75-1.76       |                     47 |                   4.7 |    1.00 |
| chr20 1x       | bin_counts |       4 | PASS            |             1.76 | 1.73-1.76      |                    80 |                  4.8 | PASS             |              1.75 | 1.75-1.75       |                     50 |                   4.7 |    1.01 |
| chr20 1x       | ctas       |       1 | PASS            |             3.36 | 3.33-3.39      |                  1428 |                  7.4 | PASS             |              3.40 | 3.36-3.41       |                   1428 |                   7.5 |    0.99 |
| chr20 1x       | ctas       |       4 | PASS            |             3.30 | 3.28-3.31      |                  1429 |                  7.3 | PASS             |              1.03 | 1.02-1.04       |                   1436 |                   9.3 |    3.20 |
| chr20 1x       | scan       |       1 | PASS            |             1.75 | 1.75-1.78      |                    43 |                  5.1 | PASS             |              1.77 | 1.76-1.77       |                     42 |                   5.1 |    0.99 |
| chr20 1x       | scan       |       4 | PASS            |             1.76 | 1.74-1.76      |                    49 |                  5.1 | PASS             |              0.57 | 0.56-0.59       |                     50 |                   6.7 |    3.09 |
| chr20 1x       | tags_map   |       1 | PASS            |             8.89 | 8.84-8.91      |                    48 |                 13.0 | PASS             |              9.14 | 9.08-9.18       |                     48 |                  13.3 |    0.97 |
| chr20 1x       | tags_map   |       4 | PASS            |             8.83 | 8.80-8.87      |                    50 |                 12.9 | PASS             |              2.44 | 2.43-2.46       |                     62 |                  14.8 |    3.62 |
| chr20 1x       | tags_named |       1 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              2.71 | 2.71-2.73       |                     45 |                   6.8 |      NA |
| chr20 1x       | tags_named |       4 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              0.83 | 0.81-0.85       |                     57 |                   8.4 |      NA |
| chr20 2x       | bin_counts |       1 | PASS            |             2.54 | 2.52-2.55      |                    79 |                  7.2 | PASS             |              2.54 | 2.53-2.55       |                     48 |                   7.2 |    1.00 |
| chr20 2x       | bin_counts |       4 | PASS            |             2.54 | 2.53-2.55      |                    80 |                  7.3 | PASS             |              2.53 | 2.52-2.54       |                     50 |                   7.1 |    1.00 |
| chr20 2x       | ctas       |       1 | PASS            |             6.48 | 6.44-6.51      |                  2812 |                 12.7 | PASS             |              6.54 | 6.47-6.60       |                   2812 |                  12.8 |    0.99 |
| chr20 2x       | ctas       |       4 | PASS            |             6.33 | 6.32-6.44      |                  2813 |                 12.6 | PASS             |              1.85 | 1.82-1.90       |                   2820 |                  15.3 |    3.42 |
| chr20 2x       | scan       |       1 | PASS            |             2.53 | 2.52-2.55      |                    42 |                  7.9 | PASS             |              2.58 | 2.58-2.59       |                     42 |                   8.0 |    0.98 |
| chr20 2x       | scan       |       4 | PASS            |             2.54 | 2.52-2.60      |                    49 |                  8.0 | PASS             |              0.96 | 0.93-1.17       |                     51 |                  10.9 |    2.65 |
| chr20 2x       | tags_map   |       1 | PASS            |            17.65 | 17.64-17.70    |                    49 |                 23.8 | PASS             |             18.22 | 18.05-18.29     |                     49 |                  24.5 |    0.97 |
| chr20 2x       | tags_map   |       4 | PASS            |            17.77 | 17.74-17.81    |                    52 |                 41.8 | PASS             |              4.98 | 4.86-5.04       |                     64 |                  27.6 |    3.57 |
| chr20 2x       | tags_named |       1 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              5.24 | 5.24-5.26       |                     47 |                  11.5 |      NA |
| chr20 2x       | tags_named |       4 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              1.54 | 1.52-1.87       |                     59 |                  14.2 |      NA |
| chr20 4x       | bin_counts |       1 | PASS            |             3.96 | 3.95-3.97      |                    79 |                 11.7 | PASS             |              3.95 | 3.93-3.96       |                     47 |                  11.6 |    1.00 |
| chr20 4x       | bin_counts |       4 | PASS            |             3.96 | 3.94-3.96      |                    80 |                 11.6 | PASS             |              3.95 | 3.94-3.99       |                     50 |                  11.5 |    1.00 |
| chr20 4x       | ctas       |       1 | PASS            |            12.77 | 12.76-12.81    |                  5581 |                 23.1 | PASS             |             12.86 | 12.83-13.01     |                   5581 |                  23.2 |    0.99 |
| chr20 4x       | ctas       |       4 | PASS            |            12.54 | 12.51-12.58    |                  5582 |                 22.9 | PASS             |              3.51 | 3.45-3.60       |                   5588 |                  27.0 |    3.57 |
| chr20 4x       | scan       |       1 | PASS            |             4.51 | 4.45-4.74      |                    42 |                 13.6 | PASS             |              4.62 | 4.56-4.63       |                     42 |                  13.7 |    0.98 |
| chr20 4x       | scan       |       4 | PASS            |             4.50 | 4.46-4.51      |                    49 |                 13.6 | PASS             |              1.95 | 1.91-2.01       |                     51 |                  21.5 |    2.31 |
| chr20 4x       | tags_map   |       1 | PASS            |            34.49 | 34.44-34.60    |                    51 |                 44.4 | PASS             |             35.62 | 35.54-35.64     |                     50 |                  45.8 |    0.97 |
| chr20 4x       | tags_map   |       4 | PASS            |            35.03 | 34.46-35.43    |                    53 |                 80.2 | PASS             |             10.03 | 9.56-10.32      |                     68 |                  53.3 |    3.49 |
| chr20 4x       | tags_named |       1 | NA              |               NA | NA             |                    NA |                   NA | PASS             |             10.18 | 10.17-10.19     |                     48 |                  20.1 |      NA |
| chr20 4x       | tags_named |       4 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              3.27 | 2.89-3.50       |                     64 |                  26.5 |      NA |
| ONT            | bin_counts |       1 | PASS            |             0.10 | 0.10-0.10      |                    41 |                  0.3 | PASS             |              0.10 | 0.10-0.10       |                     41 |                   0.3 |    1.00 |
| ONT            | bin_counts |       4 | PASS            |             0.10 | 0.10-0.10      |                    41 |                  0.3 | PASS             |              0.10 | 0.10-0.10       |                     42 |                   0.3 |    1.00 |
| ONT            | ctas       |       1 | PASS            |             0.10 | 0.10-0.10      |                    52 |                  0.3 | PASS             |              0.10 | 0.10-0.11       |                     52 |                   0.3 |    1.00 |
| ONT            | ctas       |       4 | PASS            |             0.11 | 0.11-0.11      |                    61 |                  0.3 | PASS             |              0.04 | 0.04-0.04       |                     63 |                   0.4 |    2.75 |
| ONT            | scan       |       1 | PASS            |             0.10 | 0.10-0.10      |                    41 |                  0.3 | PASS             |              0.10 | 0.10-0.10       |                     41 |                   0.3 |    1.00 |
| ONT            | scan       |       4 | PASS            |             0.10 | 0.10-0.10      |                    41 |                  0.3 | PASS             |              0.03 | 0.03-0.04       |                     45 |                   0.3 |    3.33 |
| ONT            | tags_map   |       1 | PASS            |             0.10 | 0.10-0.10      |                    43 |                  0.3 | PASS             |              0.10 | 0.10-0.10       |                     42 |                   0.3 |    1.00 |
| ONT            | tags_map   |       4 | PASS            |             0.10 | 0.10-0.10      |                    43 |                  0.3 | PASS             |              0.03 | 0.03-0.04       |                     50 |                   0.4 |    3.33 |
| ONT            | tags_named |       1 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              0.10 | 0.10-0.10       |                     41 |                   0.3 |      NA |
| ONT            | tags_named |       4 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              0.03 | 0.03-0.03       |                     45 |                   0.3 |      NA |

| key                      | value                                                            |
|:-------------------------|:-----------------------------------------------------------------|
| source_revision          | b8c36db3ff1b76c093ff601e29d1e44e70e178fd                         |
| candidate_sha256         | 772d019b6816d8651ccd08420191a390d882619d304b345d9bbf7e25e1fb1930 |
| baseline_sha256          | 2cbf19455343e1993d1d2f7df8fc47538f1c21b3bef9330a86611253f3c48678 |
| baseline_source_revision | 7474748d68a87ae102d72c58ec631eae7ba4b005                         |
| duckdb                   | v1.5.1 (Variegata) 7dbb2e646f                                    |
| nproc                    | 20                                                               |
| cpu                      | 13th Gen Intel(R) Core(TM) i5-13500                              |
| rendered_utc             | 2026-11-03 13:49:56 UTC                                          |

The chr20 rows are the 1x/2x/4x scaling ladder on one populated reference.
`tags_named` has no baseline cell because the baseline cannot name tags; compare
it with the candidate’s `tags_map` row.

The `bam_bin_counts` rows compare different thread policies, not the same work at
the same concurrency. The baseline ran its scan at bind with up to 16 threads of its
own whatever `SET threads` said, so its CPU seconds far exceed its wall time even
in the one-thread cells. The candidate runs on DuckDB worker threads and follows the
setting. On the exome, with 86 references, the baseline’s one- and four-thread cells
are therefore faster in wall time than the candidate’s; on chr20, with one populated
reference, both run one reader and match. The wall-time comparison at equal
concurrency is a thread count of 16 or more, which these cells do not include.
