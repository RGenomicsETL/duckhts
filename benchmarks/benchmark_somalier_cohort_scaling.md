Somalier cohort relatedness scaling
================

This report gives the 1×/2×/4× scaling evidence `STYLE.md` asks for, for the
multi-sample Somalier path from count evidence to pairwise relatedness. It
measures four stages, each as its own fresh process:

- `sketches`: `duckhts_somalier_prepare_sketches` over a count-evidence
  relation, written with `COPY ... TO` Parquet. One sketch row per sample.
- `pairs_reduce`: `duckhts_somalier_relatedness_all_pairs` over the stored
  sketches, reduced to a count and sums of its fields.
- `pairs_copy`: the same all-pairs result, every field, written with
  `COPY ... TO` Parquet.
- `pairs_scalar`: the per-pair function `duckhts_somalier_relatedness` over a
  self-join of the stored sketches, reduced like `pairs_reduce`. This was the
  only all-pairs form before `duckhts_somalier_relatedness_all_pairs`.

The all-pairs macro checks the contents of each sketch once. The per-pair
function checks both sketches of every pair. [Before and after](#before-and-after)
compares both with the observations of revision `a5f3b47a`.

## What was already measured

| path                                                    | existing evidence                                                                                                         | gap                                                                |
|:--------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------|
| `duckhts_somalier_vcf_counts`                           | `benchmark_sql_scaling_audit.md`: file, panel and joint doublings, three fresh processes at 1 and 4 threads, 256 MiB caps | none for memory; all cells are below the 5-second floor            |
| `duckhts_somalier_import_sites`                         | `benchmark_sql_scaling_audit.md`: three real sizes, repeated, with an RSS budget                                          | no timing verdict (below the floor)                                |
| `duckhts_somalier_bam_counts`                           | `benchmark_somalier_site_extraction.md`: one 17,000-site synthetic workload                                               | no 1×/2×/4× series, no budget                                      |
| `duckhts_somalier_sex`                                  | `benchmark_somalier_sex.md`: 500/1,000/2,000 samples, budget gated                                                        | sites are fixed; synthetic only                                    |
| panel selection (`rduckhts_somalier_find_sites`)        | `benchmark_somalier_find_sites_chr22.md`: whole chr22 and chrX, memory recorded                                           | one size per chromosome                                            |
| sketches and all-pair relatedness                       | `benchmark_somalier_reviewed.md`: one size (250 samples, 17,000 sites), 4 threads                                         | no series, no budget, no buffer or spill record: **measured here** |
| `COPY ... TO` of the relatedness result                 | `benchmark_somalier_reviewed.md`: inside one end-to-end case at one size                                                  | **measured here**                                                  |
| sketch verification, CHARR, matched contamination       | `benchmark_somalier_reviewed.md`: one size                                                                                | not measured here                                                  |
| `*_convert_parquet_sql` COPY statements (BAM, VCF, GFF) | none; listed as open finding 6 of `benchmark_sql_scaling_audit.md`                                                        | not measured here; they are not Somalier paths                     |

Open findings 5 and 6 of `benchmark_sql_scaling_audit.md` are the source of the
“multi-sample Somalier cohort and COPY workloads” item. Finding 5 asks for a
capped quadratic contract for all-pairs relatedness, which this report gives.
Finding 6 concerns the reader `*_convert_parquet_sql` COPY statements, which
this report does not measure.

## Growth dimensions and size contract

The dimensions are samples `N`, panel sites `S` and pairs `P = N (N - 1) / 2`.

- Count evidence has `N × S` rows.
- A sketch holds three bit masks of `ceiling(S / 64)` 64-bit words: `3 S / 8`
  bytes per sample. Sketch output has `N` rows.
- The all-pairs output has `P` rows of 30 fixed-width fields. It is quadratic
  in `N` and does not depend on `S`. No intermediate has `N × S × P` rows: a
  pair is computed from two sketches, not from site rows.
- Pair work is `P × S / 64` word operations.
- Caps: `max_sites` bounds `S` per sketch (default 1,000,000; hard limit
  100,000,000). The SQL has no cap on `P`. The caller chooses the pair
  relation; a `pairs_table` join gives a linear selection instead of all pairs.

The series are:

- samples: 125, 250, 500 at 17,000 sites (7,750, 31,125 and 124,750 pairs);
- sites: 17,000, 34,000, 68,000 at 250 samples;
- joint: (125 samples, 4,250 sites), (250, 8,500), (500, 17,000).

A configuration shared by two series is measured once and listed in each.

## Inputs

The synthetic count evidence is deterministic arithmetic data. It is generated
at render time into a temporary directory, which is untimed staging. It uses
the formulas of the registered `somalier-synthetic` generator
(`r/duckhtsbench/R/somalier.R`) with the number of sites as a parameter. Rows
are ordered by sample and then site, as a concatenation of per-sample count
extractions is. About 1% of the rows are unavailable (NULL counts).

The real public workload uses the registry artifacts
`roh_counts_chr20_na18507`, `roh_counts_chr20_hg00403` and
`roh_counts_chr20_hg00188`. They hold read counts that
`duckhts_somalier_bam_counts` made from three 1000 Genomes 30× CRAMs at 108,757
chr20 sites. Each file is checked against its registered row count and counts
digest before anything is measured. The render derives a panel from the sites
(dense `site_index` in position order, alleles in lexical A/B order) and
orients the reference and alternate counts into A and B counts. This
derivation is untimed and is not written to the cache.

No registered input holds count evidence for more samples. The registered
1000 Genomes phased VCFs have GT but no `FORMAT/AD`, and counts are never
inferred from GT. The real workload therefore covers the site dimension
(108,757 sites) with three samples and three pairs. It does not cover a real
cohort. Staging one would need range reads of many 30× CRAMs, which this
report does not do.

## Method

Every observation is one fresh DuckDB CLI process under `/usr/bin/time`. The
JSON profiler gives the latency of the measured statement and DuckDB’s peak
buffer and temporary bytes. `/usr/bin/time` gives peak RSS for the whole
process. `max_temp_directory_size` is zero, so a spill fails the render.
`memory_limit` is the DuckDB default. Each cell has three repetitions at one
thread and at four threads, in an order that alternates sizes. The process is
pinned to CPUs 0-3. Input Parquet files are read from a warm page cache.

Each statement consumes its output. `sketches` and `pairs_copy` write every
output column to Parquet; a second, untimed process reads the file back.
`pairs_reduce` and `pairs_scalar` sum fields of every pair.

The pair stages are checked against a reference: every field of
`duckhts_somalier_relatedness` for every pair, written once per configuration
by an untimed process. The render stops unless:

- the rows of `duckhts_somalier_relatedness_all_pairs` equal the reference
  rows, every field, for every configuration;
- the Parquet file of every `pairs_copy` observation equals the reference
  rows, every field;
- every `pairs_reduce` and `pairs_scalar` observation reports the reference’s
  row count and its sums of `jointly_called`, `ibs0` and `ibs2` exactly, and
  its sum of `relatedness` within a relative 1e-9. That sum adds doubles in an
  order the plan chooses, so it is not compared bit for bit; the row
  comparisons above compare each `relatedness` value exactly.

The pair stages read sketches that an untimed process staged with the same
`sketches` statement.

## Budgets

The budgets were declared before the first measurement, at revision
`a5f3b47a`. That run stopped at its own gate: the real cell of three pairs
used 42.5 to 44.1 MiB against a ceiling of 41.5 MiB. The first model scaled
the join output by the number of pairs, so three pairs were budgeted three
rows. DuckDB allocates whole 2,048-row vectors for the 30-field result in each
operator, whatever the row count. The model below adds that fixed term, and
the sketch relation the macro materializes. No other term changed. The
observations of the first run are kept in
[`benchmark_somalier_cohort_scaling_runs_a5f3b47a.csv`](benchmark_somalier_cohort_scaling_runs_a5f3b47a.csv).

The overhead is the median peak RSS of three empty runs per thread count. An
empty run starts the CLI, loads the extension, and copies a one-row Parquet
file to a one-row Parquet file.

Each ceiling is overhead + 3 × the live state below, and each is a gate: the
render stops if any observation exceeds its ceiling or spills. `W` is
`ceiling(S / 64)`, and a sketch is `3 × W × 8 + 256` bytes.

- **`sketches`.** The live state is:
  - the panel, which the macro and the panel digest materialize several times,
    with its row hashes and the hash table that joins it to the evidence:
    600 bytes per site;
  - the aggregate states: four masks of `W` words per sample, in one hash table
    per thread and in the combined table: `(threads + 1) × N × 4 × W × 8` bytes;
  - the sketches, in the aggregate result and in the Parquet writer:
    `2 × N` sketches;
  - the evidence scan: one Parquet row group (122,880 rows) of 120 decoded
    bytes per row and thread. The evidence itself is streamed.
- **`pairs_reduce`** and **`pairs_scalar`.** The live state is:
  - the sketch relation, compressed and decoded, the checked copy the macro
    materializes, and both sides of the join: `6 × N` sketches;
  - one join output vector per thread. The join compares `sample_id`, and the
    pair function receives both sketches of each row:
    `min(P, threads × 2,048) × 2` sketches;
  - the result vectors, fixed per thread: four operators hold one 2,048-row
    vector of the 30-field result at 16 bytes per field,
    `threads × 4 × 2,048 × 30 × 16` bytes.
- **`pairs_copy`.** The live state of `pairs_reduce`, and the Parquet writer’s
  row group: `min(P, threads × 122,880)` rows of 400 bytes.

A time exponent is the log of the median time ratio divided by the log of the
work ratio, per doubling. The work is `N × S` for `sketches` and `P × S` for
the pair stages, so a pair stage that is quadratic in samples and linear in
pairs has exponent 1. An exponent above 1.25 needs an investigation. A series
whose one-thread 1× run is under 5 seconds gets a memory verdict only.

``` sh
taskset -c 0-3 Rscript -e 'rmarkdown::render("benchmarks/benchmark_somalier_cohort_scaling.Rmd")'
```

    #> NULL

Source revision: 68d865e2c9044aa6362bcf13e7832ac957323697; `src` tree 16e9a40aad728d9cc6396f71a80600581523bc3e. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: cd1a0c4320967084a42c6f00d11f9a83cdf0dc2476068c7a5bf557d32ad5aefe. CPU affinity: 0-3. Host load before the run: 3.70, 2.05, 2.20; after: 4.97, 4.54, 3.24. Empty-run overhead: 39.5 MiB at one thread and 39.5 MiB at four. Every observation is in [`benchmark_somalier_cohort_scaling_runs.csv`](benchmark_somalier_cohort_scaling_runs.csv).

## Inputs as measured

`evidence_kib` and `sketch_kib` are the Parquet file sizes. The decoded evidence is about 120 bytes per row. `macro_differing_rows` counts the rows, over every field, that differ between `duckhts_somalier_relatedness_all_pairs` and the per-pair function.

| source    | samples |   sites | evidence_rows |   pairs | evidence_kib | sketch_kib | macro_differing_rows |
|:----------|--------:|--------:|--------------:|--------:|-------------:|-----------:|---------------------:|
| synthetic |     125 |  17,000 |     2,125,000 |   7,750 |      6,326.6 |       83.9 |                    0 |
| synthetic |     250 |  17,000 |     4,250,000 |  31,125 |     12,637.5 |      163.0 |                    0 |
| synthetic |     500 |  17,000 |     8,500,000 | 124,750 |     25,116.3 |      318.8 |                    0 |
| synthetic |     250 |  34,000 |     8,500,000 |  31,125 |     36,766.0 |      292.9 |                    0 |
| synthetic |     250 |  68,000 |    17,000,000 |  31,125 |     74,662.7 |      552.2 |                    0 |
| synthetic |     125 |   4,250 |       531,250 |   7,750 |      1,130.5 |       32.8 |                    0 |
| synthetic |     250 |   8,500 |     2,125,000 |  31,125 |      4,199.5 |       95.1 |                    0 |
| real      |       3 | 108,757 |       326,271 |       3 |      3,527.4 |      113.1 |                    0 |

## `sketches`: count evidence to sketches, written to Parquet

| dimension | scale | threads | samples |  sites | input_rows | output_rows | output_kib | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | peak_buffer_mib | spill_bytes | budget_mib | within_budget |
|:----------|------:|--------:|--------:|-------:|-----------:|------------:|-----------:|---------------:|------------:|------------:|---------------:|------------:|----------------:|------------:|-----------:|:--------------|
| samples   |     1 |       1 |     125 |  17000 |    2125000 |         125 |       86.7 |          0.306 |       0.299 |       0.306 |           62.6 |        62.9 |            41.8 |           0 |      121.7 | TRUE          |
| samples   |     2 |       1 |     250 |  17000 |    4250000 |         250 |      169.4 |          0.579 |       0.576 |       0.581 |           68.4 |        68.6 |            42.2 |           0 |      132.6 | TRUE          |
| samples   |     4 |       1 |     500 |  17000 |    8500000 |         500 |      333.9 |          1.139 |       1.130 |       2.022 |           81.4 |        81.8 |            52.8 |           0 |      154.3 | TRUE          |
| samples   |     1 |       4 |     125 |  17000 |    2125000 |         125 |       83.8 |          0.150 |       0.144 |       0.181 |           75.7 |        79.8 |            85.5 |           0 |      257.4 | TRUE          |
| samples   |     2 |       4 |     250 |  17000 |    4250000 |         250 |      157.8 |          0.258 |       0.256 |       0.260 |           82.6 |        82.7 |            86.8 |           0 |      277.4 | TRUE          |
| samples   |     4 |       4 |     500 |  17000 |    8500000 |         500 |      301.0 |          0.488 |       0.485 |       0.497 |           96.5 |        96.5 |            87.4 |           0 |      317.3 | TRUE          |
| sites     |     1 |       1 |     250 |  17000 |    4250000 |         250 |      169.4 |          0.579 |       0.576 |       0.581 |           68.4 |        68.6 |            42.2 |           0 |      132.6 | TRUE          |
| sites     |     2 |       1 |     250 |  34000 |    8500000 |         250 |      283.5 |          1.185 |       1.183 |       1.435 |           88.2 |        89.8 |            68.1 |           0 |      183.1 | TRUE          |
| sites     |     4 |       1 |     250 |  68000 |   17000000 |         250 |      542.2 |          2.399 |       2.350 |       2.417 |          128.7 |       128.9 |           124.3 |           0 |      284.0 | TRUE          |
| sites     |     1 |       4 |     250 |  17000 |    4250000 |         250 |      157.8 |          0.258 |       0.256 |       0.260 |           82.6 |        82.7 |            86.8 |           0 |      277.4 | TRUE          |
| sites     |     2 |       4 |     250 |  34000 |    8500000 |         250 |      294.4 |          0.501 |       0.501 |       0.501 |          107.8 |       108.8 |           109.3 |           0 |      346.2 | TRUE          |
| sites     |     4 |       4 |     250 |  68000 |   17000000 |         250 |      555.3 |          0.994 |       0.975 |       1.019 |          152.3 |       153.3 |           143.0 |           0 |      483.5 | TRUE          |
| joint     |     1 |       1 |     125 |   4250 |     531250 |         125 |       34.5 |          0.083 |       0.082 |       0.083 |           52.5 |        52.6 |            29.8 |           0 |       91.9 | TRUE          |
| joint     |     2 |       1 |     250 |   8500 |    2125000 |         250 |      100.9 |          0.290 |       0.285 |       0.293 |           59.0 |        59.1 |            35.3 |           0 |      107.3 | TRUE          |
| joint     |     4 |       1 |     500 |  17000 |    8500000 |         500 |      333.9 |          1.139 |       1.130 |       2.022 |           81.4 |        81.8 |            52.8 |           0 |      154.3 | TRUE          |
| joint     |     1 |       4 |     125 |   4250 |     531250 |         125 |       32.8 |          0.047 |       0.047 |       0.048 |           61.2 |        62.3 |            60.3 |           0 |      220.7 | TRUE          |
| joint     |     2 |       4 |     250 |   8500 |    2125000 |         250 |       96.0 |          0.141 |       0.139 |       0.143 |           69.7 |        69.8 |            80.0 |           0 |      243.0 | TRUE          |
| joint     |     4 |       4 |     500 |  17000 |    8500000 |         500 |      301.0 |          0.488 |       0.485 |       0.497 |           96.5 |        96.5 |            87.4 |           0 |      317.3 | TRUE          |
| real      |     1 |       1 |       3 | 108757 |     326271 |           3 |      113.1 |          0.172 |       0.170 |       0.174 |          128.9 |       129.1 |           146.7 |           0 |      270.1 | TRUE          |
| real      |     1 |       4 |       3 | 108757 |     326271 |           3 |      113.1 |          0.106 |       0.100 |       0.117 |          136.5 |       136.9 |           174.3 |           0 |      398.0 | TRUE          |

## `pairs_reduce`: all pairs, reduced

|     | dimension | scale | threads | samples |  sites | input_rows | output_rows | output_kib | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | peak_buffer_mib | spill_bytes | budget_mib | within_budget |
|:----|:----------|------:|--------:|--------:|-------:|-----------:|------------:|-----------:|---------------:|------------:|------------:|---------------:|------------:|----------------:|------------:|-----------:|:--------------|
| 21  | samples   |     1 |       1 |     125 |  17000 |        125 |        7750 |        0.1 |          0.072 |       0.071 |       0.073 |           86.9 |        87.0 |            18.2 |           0 |      142.8 | TRUE          |
| 22  | samples   |     2 |       1 |     250 |  17000 |        250 |       31125 |        0.1 |          0.208 |       0.206 |       0.211 |          100.2 |       100.4 |            24.7 |           0 |      157.1 | TRUE          |
| 23  | samples   |     4 |       1 |     500 |  17000 |        500 |      124750 |        0.1 |          0.733 |       0.725 |       0.736 |          114.0 |       114.2 |            36.8 |           0 |      185.6 | TRUE          |
| 24  | samples   |     1 |       4 |     125 |  17000 |        125 |        7750 |        0.1 |          0.073 |       0.069 |       0.075 |           95.8 |        95.8 |            18.3 |           0 |      393.2 | TRUE          |
| 25  | samples   |     2 |       4 |     250 |  17000 |        250 |       31125 |        0.1 |          0.213 |       0.206 |       0.225 |          113.2 |       115.4 |            24.9 |           0 |      424.3 | TRUE          |
| 26  | samples   |     4 |       4 |     500 |  17000 |        500 |      124750 |        0.1 |          0.785 |       0.739 |       0.842 |          135.1 |       136.4 |            36.9 |           0 |      452.8 | TRUE          |
| 27  | sites     |     1 |       1 |     250 |  17000 |        250 |       31125 |        0.1 |          0.208 |       0.206 |       0.211 |          100.2 |       100.4 |            24.7 |           0 |      157.1 | TRUE          |
| 28  | sites     |     2 |       1 |     250 |  34000 |        250 |       31125 |        0.1 |          0.385 |       0.382 |       0.389 |          160.8 |       161.0 |            36.8 |           0 |      259.3 | TRUE          |
| 29  | sites     |     4 |       1 |     250 |  68000 |        250 |       31125 |        0.1 |          1.034 |       1.022 |       1.036 |          205.8 |       205.9 |            61.4 |           0 |      463.3 | TRUE          |
| 30  | sites     |     1 |       4 |     250 |  17000 |        250 |       31125 |        0.1 |          0.213 |       0.206 |       0.225 |          113.2 |       115.4 |            24.9 |           0 |      424.3 | TRUE          |
| 31  | sites     |     2 |       4 |     250 |  34000 |        250 |       31125 |        0.1 |          0.396 |       0.386 |       0.396 |          159.7 |       165.6 |            36.9 |           0 |      750.9 | TRUE          |
| 32  | sites     |     4 |       4 |     250 |  68000 |        250 |       31125 |        0.1 |          1.045 |       1.044 |       1.134 |          211.4 |       217.5 |            61.5 |           0 |     1403.0 | TRUE          |
| 33  | joint     |     1 |       1 |     125 |   4250 |        125 |        7750 |        0.1 |          0.023 |       0.023 |       0.023 |           59.8 |        59.8 |            13.1 |           0 |       76.6 | TRUE          |
| 34  | joint     |     2 |       1 |     250 |   8500 |        250 |       31125 |        0.1 |          0.116 |       0.115 |       0.116 |           77.9 |        78.3 |            18.2 |           0 |      106.0 | TRUE          |
| 35  | joint     |     4 |       1 |     500 |  17000 |        500 |      124750 |        0.1 |          0.733 |       0.725 |       0.736 |          114.0 |       114.2 |            36.8 |           0 |      185.6 | TRUE          |
| 36  | joint     |     1 |       4 |     125 |   4250 |        125 |        7750 |        0.1 |          0.023 |       0.023 |       0.024 |           58.5 |        58.6 |            13.2 |           0 |      171.2 | TRUE          |
| 37  | joint     |     2 |       4 |     250 |   8500 |        250 |       31125 |        0.1 |          0.113 |       0.112 |       0.113 |           77.4 |        79.8 |            18.3 |           0 |      260.9 | TRUE          |
| 38  | joint     |     4 |       4 |     500 |  17000 |        500 |      124750 |        0.1 |          0.785 |       0.739 |       0.842 |          135.1 |       136.4 |            36.9 |           0 |      452.8 | TRUE          |
| 39  | real      |     1 |       1 |       3 | 108757 |          3 |           3 |        0.1 |          0.005 |       0.005 |       0.005 |           44.9 |        45.0 |            11.6 |           0 |       53.6 | TRUE          |
| 40  | real      |     1 |       4 |       3 | 108757 |          3 |           3 |        0.1 |          0.004 |       0.004 |       0.005 |           45.4 |        45.5 |            11.5 |           0 |       87.3 | TRUE          |

## `pairs_copy`: all pairs, written to Parquet

|     | dimension | scale | threads | samples |  sites | input_rows | output_rows | output_kib | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | peak_buffer_mib | spill_bytes | budget_mib | within_budget |
|:----|:----------|------:|--------:|--------:|-------:|-----------:|------------:|-----------:|---------------:|------------:|------------:|---------------:|------------:|----------------:|------------:|-----------:|:--------------|
| 41  | samples   |     1 |       1 |     125 |  17000 |        125 |        7750 |      111.9 |          0.076 |       0.075 |       0.076 |           90.9 |        90.9 |            20.2 |           0 |      151.7 | TRUE          |
| 42  | samples   |     2 |       1 |     250 |  17000 |        250 |       31125 |      366.1 |          0.226 |       0.224 |       0.230 |          103.2 |       108.2 |            34.6 |           0 |      192.7 | TRUE          |
| 43  | samples   |     4 |       1 |     500 |  17000 |        500 |      124750 |     1425.3 |          0.829 |       0.803 |       0.880 |          202.0 |       202.9 |            78.2 |           0 |      326.2 | TRUE          |
| 44  | samples   |     1 |       4 |     125 |  17000 |        125 |        7750 |      111.9 |          0.079 |       0.077 |       0.109 |           95.4 |        97.4 |            20.4 |           0 |      402.1 | TRUE          |
| 45  | samples   |     2 |       4 |     250 |  17000 |        250 |       31125 |      366.1 |          0.232 |       0.214 |       0.252 |          116.5 |       125.1 |            34.8 |           0 |      459.9 | TRUE          |
| 46  | samples   |     4 |       4 |     500 |  17000 |        500 |      124750 |     1425.3 |          0.801 |       0.796 |       0.826 |          196.1 |       203.3 |            78.3 |           0 |      595.5 | TRUE          |
| 47  | sites     |     1 |       1 |     250 |  17000 |        250 |       31125 |      366.1 |          0.226 |       0.224 |       0.230 |          103.2 |       108.2 |            34.6 |           0 |      192.7 | TRUE          |
| 48  | sites     |     2 |       1 |     250 |  34000 |        250 |       31125 |      355.1 |          0.434 |       0.432 |       0.541 |          148.6 |       148.9 |            45.2 |           0 |      294.9 | TRUE          |
| 49  | sites     |     4 |       1 |     250 |  68000 |        250 |       31125 |      351.3 |          1.053 |       1.043 |       1.093 |          214.1 |       214.2 |            66.8 |           0 |      499.0 | TRUE          |
| 50  | sites     |     1 |       4 |     250 |  17000 |        250 |       31125 |      366.1 |          0.232 |       0.214 |       0.252 |          116.5 |       125.1 |            34.8 |           0 |      459.9 | TRUE          |
| 51  | sites     |     2 |       4 |     250 |  34000 |        250 |       31125 |      355.1 |          0.422 |       0.412 |       0.437 |          194.7 |       203.3 |            45.3 |           0 |      786.5 | TRUE          |
| 52  | sites     |     4 |       4 |     250 |  68000 |        250 |       31125 |      351.3 |          1.062 |       1.056 |       1.129 |          214.4 |       219.4 |            66.9 |           0 |     1438.6 | TRUE          |
| 53  | joint     |     1 |       1 |     125 |   4250 |        125 |        7750 |      116.3 |          0.025 |       0.025 |       0.026 |           60.2 |        61.0 |            15.1 |           0 |       85.5 | TRUE          |
| 54  | joint     |     2 |       1 |     250 |   8500 |        250 |       31125 |      356.2 |          0.122 |       0.121 |       0.127 |           87.0 |        87.3 |            41.3 |           0 |      141.6 | TRUE          |
| 55  | joint     |     4 |       1 |     500 |  17000 |        500 |      124750 |     1425.3 |          0.829 |       0.803 |       0.880 |          202.0 |       202.9 |            78.2 |           0 |      326.2 | TRUE          |
| 56  | joint     |     1 |       4 |     125 |   4250 |        125 |        7750 |      116.3 |          0.028 |       0.028 |       0.039 |           62.0 |        62.2 |            15.2 |           0 |      180.1 | TRUE          |
| 57  | joint     |     2 |       4 |     250 |   8500 |        250 |       31125 |      356.2 |          0.127 |       0.126 |       0.134 |          100.4 |       103.4 |            41.4 |           0 |      296.6 | TRUE          |
| 58  | joint     |     4 |       4 |     500 |  17000 |        500 |      124750 |     1425.3 |          0.801 |       0.796 |       0.826 |          196.1 |       203.3 |            78.3 |           0 |      595.5 | TRUE          |
| 59  | real      |     1 |       1 |       3 | 108757 |          3 |           3 |        4.6 |          0.005 |       0.005 |       0.005 |           44.9 |        45.1 |             9.6 |           0 |       53.6 | TRUE          |
| 60  | real      |     1 |       4 |       3 | 108757 |          3 |           3 |        4.6 |          0.004 |       0.004 |       0.004 |           46.2 |        46.3 |            10.0 |           0 |       87.3 | TRUE          |

## `pairs_scalar`: all pairs through the per-pair function, reduced

|     | dimension | scale | threads | samples |  sites | input_rows | output_rows | output_kib | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | peak_buffer_mib | spill_bytes | budget_mib | within_budget |
|:----|:----------|------:|--------:|--------:|-------:|-----------:|------------:|-----------:|---------------:|------------:|------------:|---------------:|------------:|----------------:|------------:|-----------:|:--------------|
| 61  | samples   |     1 |       1 |     125 |  17000 |        125 |        7750 |        0.1 |          0.164 |       0.163 |       0.164 |           85.6 |        85.8 |            18.4 |           0 |      142.8 | TRUE          |
| 62  | samples   |     2 |       1 |     250 |  17000 |        250 |       31125 |        0.1 |          0.583 |       0.575 |       0.587 |           99.7 |       104.5 |            27.4 |           0 |      157.1 | TRUE          |
| 63  | samples   |     4 |       1 |     500 |  17000 |        500 |      124750 |        0.1 |          2.253 |       2.250 |       2.274 |          114.7 |       114.7 |            45.5 |           0 |      185.6 | TRUE          |
| 64  | samples   |     1 |       4 |     125 |  17000 |        125 |        7750 |        0.1 |          0.164 |       0.162 |       0.170 |           94.9 |        95.0 |            18.4 |           0 |      393.2 | TRUE          |
| 65  | samples   |     2 |       4 |     250 |  17000 |        250 |       31125 |        0.1 |          0.580 |       0.579 |       0.584 |           96.9 |        97.5 |            27.8 |           0 |      424.3 | TRUE          |
| 66  | samples   |     4 |       4 |     500 |  17000 |        500 |      124750 |        0.1 |          2.256 |       2.253 |       2.259 |          116.9 |       121.5 |            45.7 |           0 |      452.8 | TRUE          |
| 67  | sites     |     1 |       1 |     250 |  17000 |        250 |       31125 |        0.1 |          0.583 |       0.575 |       0.587 |           99.7 |       104.5 |            27.4 |           0 |      157.1 | TRUE          |
| 68  | sites     |     2 |       1 |     250 |  34000 |        250 |       31125 |        0.1 |          1.137 |       1.135 |       1.138 |          161.3 |       161.4 |            45.7 |           0 |      259.3 | TRUE          |
| 69  | sites     |     4 |       1 |     250 |  68000 |        250 |       31125 |        0.1 |          2.530 |       2.513 |       2.534 |          202.7 |       202.7 |            81.7 |           0 |      463.3 | TRUE          |
| 70  | sites     |     1 |       4 |     250 |  17000 |        250 |       31125 |        0.1 |          0.580 |       0.579 |       0.584 |           96.9 |        97.5 |            27.8 |           0 |      424.3 | TRUE          |
| 71  | sites     |     2 |       4 |     250 |  34000 |        250 |       31125 |        0.1 |          1.140 |       1.136 |       1.143 |          146.9 |       154.1 |            45.4 |           0 |      750.9 | TRUE          |
| 72  | sites     |     4 |       4 |     250 |  68000 |        250 |       31125 |        0.1 |          2.709 |       2.699 |       2.715 |          209.1 |       209.5 |            81.7 |           0 |     1403.0 | TRUE          |
| 73  | joint     |     1 |       1 |     125 |   4250 |        125 |        7750 |        0.1 |          0.047 |       0.047 |       0.048 |           54.4 |        54.5 |            10.8 |           0 |       76.6 | TRUE          |
| 74  | joint     |     2 |       1 |     250 |   8500 |        250 |       31125 |        0.1 |          0.302 |       0.302 |       0.303 |           73.2 |        73.3 |            18.4 |           0 |      106.0 | TRUE          |
| 75  | joint     |     4 |       1 |     500 |  17000 |        500 |      124750 |        0.1 |          2.253 |       2.250 |       2.274 |          114.7 |       114.7 |            45.5 |           0 |      185.6 | TRUE          |
| 76  | joint     |     1 |       4 |     125 |   4250 |        125 |        7750 |        0.1 |          0.047 |       0.047 |       0.048 |           55.6 |        56.0 |            10.9 |           0 |      171.2 | TRUE          |
| 77  | joint     |     2 |       4 |     250 |   8500 |        250 |       31125 |        0.1 |          0.302 |       0.302 |       0.307 |           75.2 |        83.6 |            18.4 |           0 |      260.9 | TRUE          |
| 78  | joint     |     4 |       4 |     500 |  17000 |        500 |      124750 |        0.1 |          2.256 |       2.253 |       2.259 |          116.9 |       121.5 |            45.7 |           0 |      452.8 | TRUE          |
| 79  | real      |     1 |       1 |       3 | 108757 |          3 |           3 |        0.1 |          0.003 |       0.003 |       0.003 |           42.2 |        42.6 |             7.4 |           0 |       53.6 | TRUE          |
| 80  | real      |     1 |       4 |       3 | 108757 |          3 |           3 |        0.1 |          0.003 |       0.003 |       0.003 |           43.1 |        43.2 |             7.4 |           0 |       87.3 | TRUE          |

## Before and after

`before` is the `pairs_reduce` stage of revision `a5f3b47a`, which ran the
per-pair function over a self-join; its observations are in
[`benchmark_somalier_cohort_scaling_runs_a5f3b47a.csv`](benchmark_somalier_cohort_scaling_runs_a5f3b47a.csv).
`per_pair_function` and `all_pairs_macro` are the `pairs_scalar` and
`pairs_reduce` stages of this revision. All three produce the same rows from
the same sketches, with the same statement shape, host, CPU set and DuckDB
CLI. The two runs were made at different times on a shared host, so the
ratios carry that noise. Times are median seconds of three fresh processes.

| samples |  sites |   pairs | threads | before | per_pair_function | all_pairs_macro | macro_speedup |
|--------:|-------:|--------:|--------:|-------:|------------------:|----------------:|--------------:|
|     125 |  4,250 |   7,750 |       1 |  0.081 |             0.047 |           0.023 |         3.550 |
|     250 |  8,500 |  31,125 |       1 |  0.574 |             0.302 |           0.116 |         4.969 |
|     125 | 17,000 |   7,750 |       1 |  0.298 |             0.164 |           0.072 |         4.150 |
|     250 | 17,000 |  31,125 |       1 |  1.121 |             0.583 |           0.208 |         5.399 |
|     500 | 17,000 | 124,750 |       1 |  4.399 |             2.253 |           0.733 |         6.002 |
|     250 | 34,000 |  31,125 |       1 |  2.201 |             1.137 |           0.385 |         5.710 |
|     250 | 68,000 |  31,125 |       1 |  4.633 |             2.530 |           1.034 |         4.481 |
|     125 |  4,250 |   7,750 |       4 |  0.081 |             0.047 |           0.023 |         3.469 |
|     250 |  8,500 |  31,125 |       4 |  0.579 |             0.302 |           0.113 |         5.123 |
|     125 | 17,000 |   7,750 |       4 |  0.295 |             0.164 |           0.073 |         4.045 |
|     250 | 17,000 |  31,125 |       4 |  1.125 |             0.580 |           0.213 |         5.285 |
|     500 | 17,000 | 124,750 |       4 |  4.413 |             2.256 |           0.785 |         5.619 |
|     250 | 34,000 |  31,125 |       4 |  2.204 |             1.140 |           0.396 |         5.568 |
|     250 | 68,000 |  31,125 |       4 |  4.834 |             2.709 |           1.045 |         4.627 |

## Time exponents and memory ratios

Time exponents per doubling, against the work of each stage, and peak-RSS ratios to 1×:

| stage        | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | timing_verdict     | exponent_above_1.25 |
|:-------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|:-------------------|:--------------------|
| sketches     | samples   |       1 |            0.923 |            0.975 |        1.093 |        1.300 | none: 1x under 5 s | FALSE               |
| sketches     | samples   |       4 |            0.779 |            0.919 |        1.091 |        1.275 | none: 1x under 5 s | FALSE               |
| sketches     | sites     |       1 |            1.032 |            1.018 |        1.289 |        1.882 | none: 1x under 5 s | FALSE               |
| sketches     | sites     |       4 |            0.956 |            0.990 |        1.305 |        1.844 | none: 1x under 5 s | FALSE               |
| sketches     | joint     |       1 |            0.905 |            0.986 |        1.124 |        1.550 | none: 1x under 5 s | FALSE               |
| sketches     | joint     |       4 |            0.789 |            0.897 |        1.139 |        1.577 | none: 1x under 5 s | FALSE               |
| pairs_reduce | samples   |       1 |            0.764 |            0.909 |        1.153 |        1.312 | none: 1x under 5 s | FALSE               |
| pairs_reduce | samples   |       4 |            0.771 |            0.940 |        1.182 |        1.410 | none: 1x under 5 s | FALSE               |
| pairs_reduce | sites     |       1 |            0.892 |            1.424 |        1.605 |        2.054 | none: 1x under 5 s | TRUE                |
| pairs_reduce | sites     |       4 |            0.894 |            1.400 |        1.411 |        1.867 | none: 1x under 5 s | TRUE                |
| pairs_reduce | joint     |       1 |            0.777 |            0.888 |        1.303 |        1.906 | none: 1x under 5 s | FALSE               |
| pairs_reduce | joint     |       4 |            0.756 |            0.931 |        1.323 |        2.309 | none: 1x under 5 s | FALSE               |
| pairs_copy   | samples   |       1 |            0.787 |            0.936 |        1.135 |        2.222 | none: 1x under 5 s | FALSE               |
| pairs_copy   | samples   |       4 |            0.776 |            0.893 |        1.221 |        2.056 | none: 1x under 5 s | FALSE               |
| pairs_copy   | sites     |       1 |            0.940 |            1.279 |        1.440 |        2.075 | none: 1x under 5 s | TRUE                |
| pairs_copy   | sites     |       4 |            0.865 |            1.332 |        1.671 |        1.840 | none: 1x under 5 s | TRUE                |
| pairs_copy   | joint     |       1 |            0.751 |            0.922 |        1.445 |        3.355 | none: 1x under 5 s | FALSE               |
| pairs_copy   | joint     |       4 |            0.719 |            0.886 |        1.619 |        3.163 | none: 1x under 5 s | FALSE               |
| pairs_scalar | samples   |       1 |            0.913 |            0.974 |        1.165 |        1.340 | none: 1x under 5 s | FALSE               |
| pairs_scalar | samples   |       4 |            0.908 |            0.978 |        1.021 |        1.232 | none: 1x under 5 s | FALSE               |
| pairs_scalar | sites     |       1 |            0.964 |            1.154 |        1.618 |        2.033 | none: 1x under 5 s | FALSE               |
| pairs_scalar | sites     |       4 |            0.974 |            1.249 |        1.516 |        2.158 | none: 1x under 5 s | FALSE               |
| pairs_scalar | joint     |       1 |            0.893 |            0.966 |        1.346 |        2.108 | none: 1x under 5 s | FALSE               |
| pairs_scalar | joint     |       4 |            0.889 |            0.966 |        1.353 |        2.103 | none: 1x under 5 s | FALSE               |

## Findings

Every observation is within its declared budget and none spills; the render
enforces both. The largest observation uses
84% of its ceiling.

The largest all-pairs cell (500 samples, 17,000 sites, 124,750 pairs) takes
0.83 s on one thread and
0.8 s on four when written to Parquet, at
202 and
196.1 MiB peak RSS. The Parquet output is
1425.3 KiB.

On the real counts, one thread builds three sketches of 108,757 sites in
0.172 s at
128.9 MiB peak RSS, and writes the three pairs in
0.005 s at
44.9 MiB.

For 124,750 pairs of 17,000 sites on one thread, the per-pair function took
4.4 s at revision `a5f3b47a` and takes
2.25 s now; the all-pairs macro takes
0.73 s, 6 times faster than before.
The pair stages still grow with pairs times sites, which is the size of the
comparison itself. They do not use more than one thread:
0.73 s on one thread and
0.79 s on four. DuckDB splits a
scan by row groups, and a relation of a few hundred sketches is one row group.

24 of the 24 series have a one-thread 1× run under 5 seconds and get a memory verdict only.
4 series have a time exponent above 1.25 at either doubling.
The pair stages along sites exceed 1.25 at the second doubling, from 34,000 to
68,000 sites. The cause is not established. One candidate, not tested, is that
250 sketches of 68,000 sites no longer fit the CPU cache. Flags on the
`sketches` stage follow the host load: the minimum times of that stage scale
linearly in the tables above.

## Limitations

- The host was shared with other work during the run. The load averages are
  recorded above. Times can be inflated and their spread widened; peak RSS,
  buffer and spill values do not depend on the load.
- The series are synthetic. The real workload has three samples and three
  pairs, so no real cohort is measured.
- The operating scale this report supports is what it measures: up to 500
  samples, 68,000 synthetic sites, 108,757 real sites and 124,750 pairs.
- Evidence is ordered by sample. Evidence that interleaves samples is not
  measured.
- The pair comparison runs in one thread, and the number of output pairs is
  not capped.
- Count extraction, sketch verification, CHARR, matched contamination, the R
  wrappers and the reader `*_convert_parquet_sql` COPY statements are not
  measured here.
- The DuckDB CLI is v1.5.1 (Variegata) 7dbb2e646f; the R package and other hosts can
  differ.
