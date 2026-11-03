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
| exome          | bin_counts |       1 | PASS            |             2.75 | 2.57-3.04      |                   301 |                 45.6 | PASS             |             10.29 | 10.22-10.31     |                     61 |                  25.4 |    0.27 |
| exome          | bin_counts |       4 | PASS            |             2.96 | 2.73-2.99      |                   302 |                 48.5 | PASS             |              3.66 | 3.63-3.67       |                    105 |                  37.3 |    0.81 |
| exome          | scan       |       1 | PASS            |            10.28 | 10.26-10.33    |                    56 |                 27.1 | PASS             |             10.40 | 10.36-10.41     |                     56 |                  27.5 |    0.99 |
| exome          | scan       |       4 | PASS            |             3.83 | 3.19-3.84      |                   103 |                 41.8 | PASS             |              3.83 | 3.81-3.84       |                    102 |                  42.3 |    1.00 |
| exome          | tags_map   |       1 | PASS            |            49.74 | 49.57-49.93    |                    58 |                 72.6 | PASS             |             50.28 | 50.14-50.31     |                     58 |                  73.3 |    0.99 |
| exome          | tags_map   |       4 | PASS            |            15.80 | 14.52-15.99    |                   110 |                 95.3 | PASS             |             15.75 | 15.73-15.76     |                    108 |                  95.4 |    1.00 |
| exome          | tags_named |       1 | NA              |               NA | NA             |                    NA |                   NA | PASS             |             12.26 | 12.22-12.29     |                     56 |                  34.8 |      NA |
| exome          | tags_named |       4 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              4.78 | 4.43-4.83       |                    101 |                  53.0 |      NA |
| exome unmapped | bin_counts |       1 | PASS            |             0.00 | 0.00-0.00      |                    47 |                  0.1 | PASS             |              0.00 | 0.00-0.00       |                     44 |                   0.0 |     NaN |
| exome unmapped | bin_counts |       4 | PASS            |             0.00 | 0.00-0.00      |                    47 |                  0.1 | PASS             |              0.00 | 0.00-0.00       |                     45 |                   0.1 |     NaN |
| exome unmapped | ctas       |       1 | PASS            |             0.41 | 0.40-0.41      |                   194 |                  0.8 | PASS             |              0.42 | 0.41-0.42       |                    194 |                   0.9 |    0.98 |
| exome unmapped | ctas       |       4 | PASS            |             0.41 | 0.40-0.41      |                   194 |                  0.8 | PASS             |              0.41 | 0.40-0.50       |                    194 |                   0.9 |    1.00 |
| exome unmapped | scan       |       1 | PASS            |             0.17 | 0.17-0.17      |                    41 |                  0.6 | PASS             |              0.17 | 0.17-0.17       |                     41 |                   0.6 |    1.00 |
| exome unmapped | scan       |       4 | PASS            |             0.17 | 0.17-0.17      |                    42 |                  0.6 | PASS             |              0.17 | 0.17-0.17       |                     42 |                   0.6 |    1.00 |
| chr20 1x       | bin_counts |       1 | PASS            |             1.76 | 1.76-1.79      |                    81 |                  4.8 | PASS             |              1.75 | 1.75-1.76       |                     48 |                   4.7 |    1.01 |
| chr20 1x       | bin_counts |       4 | PASS            |             1.74 | 1.74-1.75      |                    80 |                  4.8 | PASS             |              1.75 | 1.75-1.76       |                     50 |                   4.8 |    0.99 |
| chr20 1x       | ctas       |       1 | PASS            |             3.36 | 3.33-3.38      |                  1428 |                  7.4 | PASS             |              3.40 | 3.37-3.42       |                   1428 |                   7.5 |    0.99 |
| chr20 1x       | ctas       |       4 | PASS            |             3.31 | 3.25-3.31      |                  1429 |                  7.3 | PASS             |              1.03 | 1.01-1.04       |                   1438 |                   9.1 |    3.21 |
| chr20 1x       | scan       |       1 | PASS            |             1.75 | 1.73-1.75      |                    42 |                  5.0 | PASS             |              1.76 | 1.75-1.77       |                     42 |                   5.1 |    0.99 |
| chr20 1x       | scan       |       4 | PASS            |             1.76 | 1.76-1.80      |                    49 |                  6.9 | PASS             |              0.58 | 0.58-0.59       |                     51 |                   6.7 |    3.03 |
| chr20 1x       | tags_map   |       1 | PASS            |             8.89 | 8.86-8.93      |                    48 |                 13.0 | PASS             |              9.03 | 8.98-9.05       |                     48 |                  13.2 |    0.98 |
| chr20 1x       | tags_map   |       4 | PASS            |             8.88 | 8.86-8.97      |                    50 |                 13.0 | PASS             |              2.43 | 2.40-2.44       |                     62 |                  14.8 |    3.65 |
| chr20 1x       | tags_named |       1 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              2.68 | 2.67-2.70       |                     45 |                   6.8 |      NA |
| chr20 1x       | tags_named |       4 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              0.83 | 0.82-0.84       |                     57 |                   8.4 |      NA |
| chr20 2x       | bin_counts |       1 | PASS            |             2.56 | 2.54-2.56      |                    79 |                  7.3 | PASS             |              2.54 | 2.54-2.55       |                     48 |                   7.2 |    1.01 |
| chr20 2x       | bin_counts |       4 | PASS            |             2.54 | 2.53-2.54      |                    79 |                  7.2 | PASS             |              2.54 | 2.54-2.55       |                     50 |                   7.2 |    1.00 |
| chr20 2x       | ctas       |       1 | PASS            |             6.54 | 6.45-6.54      |                  2812 |                 12.7 | PASS             |              6.54 | 6.54-6.57       |                   2812 |                  12.9 |    1.00 |
| chr20 2x       | ctas       |       4 | PASS            |             6.41 | 6.33-6.46      |                  2814 |                 12.7 | PASS             |              1.83 | 1.81-1.86       |                   2821 |                  15.2 |    3.50 |
| chr20 2x       | scan       |       1 | PASS            |             2.52 | 2.50-2.53      |                    42 |                  7.9 | PASS             |              2.57 | 2.57-2.59       |                     43 |                   8.0 |    0.98 |
| chr20 2x       | scan       |       4 | PASS            |             2.54 | 2.53-2.62      |                    49 |                 10.6 | PASS             |              1.13 | 0.96-1.14       |                     50 |                  13.0 |    2.25 |
| chr20 2x       | tags_map   |       1 | PASS            |            17.59 | 17.58-17.61    |                    49 |                 23.8 | PASS             |             17.86 | 17.82-17.97     |                     50 |                  24.1 |    0.98 |
| chr20 2x       | tags_map   |       4 | PASS            |            17.70 | 17.59-17.72    |                    52 |                 23.9 | PASS             |              4.76 | 4.75-4.78       |                     64 |                  26.7 |    3.72 |
| chr20 2x       | tags_named |       1 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              5.26 | 5.26-5.28       |                     47 |                  11.6 |      NA |
| chr20 2x       | tags_named |       4 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              1.54 | 1.54-1.76       |                     59 |                  14.2 |      NA |
| chr20 4x       | bin_counts |       1 | PASS            |             3.97 | 3.94-3.99      |                    80 |                 11.7 | PASS             |              3.97 | 3.95-3.98       |                     48 |                  11.7 |    1.00 |
| chr20 4x       | bin_counts |       4 | PASS            |             3.95 | 3.92-3.97      |                    80 |                 11.6 | PASS             |              3.95 | 3.94-3.96       |                     50 |                  11.6 |    1.00 |
| chr20 4x       | ctas       |       1 | PASS            |            12.72 | 12.68-12.82    |                  5581 |                 23.1 | PASS             |             12.86 | 12.75-12.86     |                   5580 |                  23.2 |    0.99 |
| chr20 4x       | ctas       |       4 | PASS            |            12.47 | 12.42-12.51    |                  5582 |                 22.6 | PASS             |              3.57 | 3.55-3.60       |                   5588 |                  27.2 |    3.49 |
| chr20 4x       | scan       |       1 | PASS            |             4.46 | 4.45-4.48      |                    43 |                 13.4 | PASS             |              4.60 | 4.59-4.62       |                     43 |                  13.8 |    0.97 |
| chr20 4x       | scan       |       4 | PASS            |             4.44 | 4.40-4.48      |                    49 |                 13.3 | PASS             |              1.98 | 1.80-2.00       |                     51 |                  21.5 |    2.24 |
| chr20 4x       | tags_map   |       1 | PASS            |            34.58 | 34.52-34.72    |                    50 |                 44.6 | PASS             |             34.97 | 34.87-35.10     |                     50 |                  45.1 |    0.99 |
| chr20 4x       | tags_map   |       4 | PASS            |            34.53 | 34.32-34.57    |                    53 |                 44.5 | PASS             |              9.36 | 9.34-9.45       |                     68 |                  50.0 |    3.69 |
| chr20 4x       | tags_named |       1 | NA              |               NA | NA             |                    NA |                   NA | PASS             |             10.19 | 10.13-10.20     |                     47 |                  20.1 |      NA |
| chr20 4x       | tags_named |       4 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              3.23 | 2.88-3.50       |                     63 |                  26.4 |      NA |
| ONT            | bin_counts |       1 | PASS            |             0.10 | 0.10-0.10      |                    41 |                  0.3 | PASS             |              0.10 | 0.10-0.10       |                     41 |                   0.3 |    1.00 |
| ONT            | bin_counts |       4 | PASS            |             0.10 | 0.10-0.10      |                    41 |                  0.3 | PASS             |              0.10 | 0.10-0.10       |                     42 |                   0.3 |    1.00 |
| ONT            | ctas       |       1 | PASS            |             0.10 | 0.10-0.10      |                    52 |                  0.3 | PASS             |              0.11 | 0.11-0.12       |                     52 |                   0.3 |    0.91 |
| ONT            | ctas       |       4 | PASS            |             0.11 | 0.11-0.12      |                    62 |                  0.3 | PASS             |              0.04 | 0.04-0.04       |                     63 |                   0.4 |    2.75 |
| ONT            | scan       |       1 | PASS            |             0.10 | 0.10-0.10      |                    41 |                  0.3 | PASS             |              0.10 | 0.10-0.11       |                     41 |                   0.3 |    1.00 |
| ONT            | scan       |       4 | PASS            |             0.10 | 0.10-0.10      |                    41 |                  0.3 | PASS             |              0.03 | 0.03-0.04       |                     45 |                   0.4 |    3.33 |
| ONT            | tags_map   |       1 | PASS            |             0.10 | 0.10-0.10      |                    42 |                  0.3 | PASS             |              0.10 | 0.10-0.11       |                     43 |                   0.3 |    1.00 |
| ONT            | tags_map   |       4 | PASS            |             0.10 | 0.10-0.10      |                    43 |                  0.3 | PASS             |              0.03 | 0.03-0.04       |                     49 |                   0.4 |    3.33 |
| ONT            | tags_named |       1 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              0.10 | 0.10-0.10       |                     41 |                   0.3 |      NA |
| ONT            | tags_named |       4 | NA              |               NA | NA             |                    NA |                   NA | PASS             |              0.03 | 0.03-0.03       |                     45 |                   0.3 |      NA |

| key                      | value                                                            |
|:-------------------------|:-----------------------------------------------------------------|
| source_revision          | 28cef804dc337cd803ae6ffd76d1a5e6b3ac9f43                         |
| candidate_sha256         | 172641b77e2893829bf5bcc849b0d6385156079e9bc4c03174b2263a83deed56 |
| baseline_sha256          | 2cbf19455343e1993d1d2f7df8fc47538f1c21b3bef9330a86611253f3c48678 |
| baseline_source_revision | 7474748d68a87ae102d72c58ec631eae7ba4b005                         |
| duckdb                   | v1.5.1 (Variegata) 7dbb2e646f                                    |
| nproc                    | 20                                                               |
| cpu                      | 13th Gen Intel(R) Core(TM) i5-13500                              |
| rendered_utc             | 2026-11-03 09:30:13 UTC                                          |

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
