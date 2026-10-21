Low-pass query compatibility: DuckDB 1.5 and 2.0
================

## Workload and checks

This repeats the synthetic scaling workload from [the low-pass
calibration report](benchmark_count_mixture_lowpass.md). The existing
`count_mixture_lowpass/memory_worker.R` generates the input view. No
external corpus is needed.

The site axis has 100,000, 200,000, and 400,000 sites, with one sample
and candidate. The group axis has 100,000 sites and (samples,
candidates) of (1,1), (2,1), and (2,2). Each query returns one
unrelated-model row per sample/candidate pair. All positions fit within
one 1 Mb block per group. This tests input scaling, not scaling in the
number of jackknife blocks.

Each configuration runs three times, sequentially, in fresh R processes.
Each query uses one DuckDB thread, a 3 GB memory limit, and CPU affinity
0–5. Query time excludes startup, extension loading, and view creation.
Peak process RSS includes those costs. DuckDB profiling supplies peak
buffers and spill.

The memory ceiling is unchanged: 512 MiB plus three times modeled live
state (768 bytes per SQL input row and 1,024 bytes per native
block/group). Every run must satisfy that ceiling for RSS and buffers,
without spill. Results are aligned by
`(sample, candidate, relationship_model)`; SQL row order is not
specified. Duplicate keys fail the check. Every column must then be
identical to baseline, including counters, fit values, statuses, and
missing values. No rows or results are filtered.

The baseline is built from source revision
`0d765acb8f8eaf16931852ac37079c5668cd6c51`. The candidate runs on both R
DuckDB engines. Comparing baseline and candidate on DuckDB 1.5 isolates
the SQL change; comparing engines also includes engine changes. The
baseline query is rejected by DuckDB 2.0’s window binder, so it has no
2.0 timing.

## Reproduction

Build the baseline from an archive of its revision and build the
candidate with `make release -j2`. Run from the repository root after
other build and test jobs finish:

``` r
rmarkdown::render("benchmarks/benchmark_duckdb2_lowpass.Rmd", params = list(
  run = TRUE,
  baseline_extension = Sys.getenv("BASELINE_EXTENSION"),
  candidate_extension = Sys.getenv("CANDIDATE_EXTENSION"),
  duckdb2_library = Sys.getenv("DUCKDB2_R_LIB")
))
```

Both extensions use the same release build settings. The measurement
chunk starts each worker with `taskset -c 0-5`. Set `run = FALSE` to
render stored results.

## Results

|     | variant       | engine                                   | source_revision                          | extension_sha256                                                 |
|:----|:--------------|:-----------------------------------------|:-----------------------------------------|:-----------------------------------------------------------------|
| 1   | baseline_1.5  | v1.5.5 d8cdaa33fda R-package-1.5.5       | 0d765acb8f8eaf16931852ac37079c5668cd6c51 | 7a549d1ea0a9b6b07ec499cdce651f2c4d89f107dea6fefca9565b55d0d18394 |
| 19  | candidate_1.5 | v1.5.5 d8cdaa33fda R-package-1.5.5       | 7111af111e98d36eb387bde0f037461693a55648 | bc128e51f31a3d89ebf58402f53c0409f660467d29a55f36d934f1814d10320d |
| 37  | candidate_2.0 | v2.0.0 d54cbbe3df R-package-1.99.99.9001 | 7111af111e98d36eb387bde0f037461693a55648 | bc128e51f31a3d89ebf58402f53c0409f660467d29a55f36d934f1814d10320d |

| R_version                    | cpu_model                           |
|:-----------------------------|:------------------------------------|
| R version 4.6.0 (2026-04-24) | 13th Gen Intel(R) Core(TM) i5-13500 |

| variant       | axis   | factor | input_rows | output_rows | threads | median_seconds | min_seconds | max_seconds | peak_rss_MiB | peak_buffers_MiB | ceiling_MiB |
|:--------------|:-------|-------:|-----------:|------------:|--------:|---------------:|------------:|------------:|-------------:|-----------------:|------------:|
| baseline_1.5  | groups |      1 |      1e+05 |           1 |       1 |          0.168 |       0.165 |       0.169 |      178.785 |           90.547 |     731.729 |
| candidate_1.5 | groups |      1 |      1e+05 |           1 |       1 |          0.163 |       0.162 |       0.164 |      178.801 |           90.531 |     731.729 |
| candidate_2.0 | groups |      1 |      1e+05 |           1 |       1 |          0.171 |       0.168 |       0.172 |      193.629 |           44.600 |     731.729 |
| baseline_1.5  | sites  |      1 |      1e+05 |           1 |       1 |          0.168 |       0.168 |       0.168 |      178.785 |           90.547 |     731.729 |
| candidate_1.5 | sites  |      1 |      1e+05 |           1 |       1 |          0.163 |       0.162 |       0.163 |      178.637 |           90.531 |     731.729 |
| candidate_2.0 | sites  |      1 |      1e+05 |           1 |       1 |          0.170 |       0.170 |       0.173 |      193.484 |           44.600 |     731.729 |
| baseline_1.5  | groups |      2 |      2e+05 |           2 |       1 |          0.322 |       0.321 |       0.326 |      251.309 |          205.031 |     951.459 |
| candidate_1.5 | groups |      2 |      2e+05 |           2 |       1 |          0.326 |       0.318 |       0.326 |      251.285 |          205.031 |     951.459 |
| candidate_2.0 | groups |      2 |      2e+05 |           2 |       1 |          0.322 |       0.321 |       0.322 |      223.180 |           85.162 |     951.459 |
| baseline_1.5  | sites  |      2 |      2e+05 |           1 |       1 |          0.340 |       0.339 |       0.342 |      235.137 |          174.623 |     951.456 |
| candidate_1.5 | sites  |      2 |      2e+05 |           1 |       1 |          0.328 |       0.328 |       0.332 |      234.949 |          174.584 |     951.456 |
| candidate_2.0 | sites  |      2 |      2e+05 |           1 |       1 |          0.333 |       0.333 |       0.336 |      233.332 |           85.162 |     951.456 |
| baseline_1.5  | groups |      4 |      4e+05 |           4 |       1 |          0.545 |       0.544 |       0.546 |      271.180 |          237.344 |    1390.918 |
| candidate_1.5 | groups |      4 |      4e+05 |           4 |       1 |          0.535 |       0.532 |       0.538 |      271.305 |          237.402 |    1390.918 |
| candidate_2.0 | groups |      4 |      4e+05 |           4 |       1 |          0.555 |       0.553 |       0.555 |      281.477 |          118.149 |    1390.918 |
| baseline_1.5  | sites  |      4 |      4e+05 |           1 |       1 |          0.687 |       0.685 |       0.689 |      347.160 |          338.234 |    1390.909 |
| candidate_1.5 | sites  |      4 |      4e+05 |           1 |       1 |          0.667 |       0.659 |       0.668 |      347.301 |          338.219 |    1390.909 |
| candidate_2.0 | sites  |      4 |      4e+05 |           1 |       1 |          0.665 |       0.665 |       0.674 |      311.520 |          166.848 |    1390.909 |

All 54 results match their baseline exactly. All runs satisfy the memory
ceilings without spill. The one-times workloads are shorter than five
seconds; under the existing declaration, this does not support a
timing-exponent verdict. The table reports all timing ranges rather than
selecting favorable runs.

The historical identical workload is retained in
[count_mixture_lowpass/results/scaling.csv](count_mixture_lowpass/results/scaling.csv).
Its environment differs from this run, so it is context rather than a
regression control. This report measures only the low-pass query. It
does not measure cgranges throughput or BGZF compression/decompression
throughput.
