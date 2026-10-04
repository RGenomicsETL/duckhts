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
output column to Parquet; a second, untimed process reads the file back and
checks the row count. `pairs_reduce` and `pairs_scalar` sum fields of every
pair. The three pair stages must agree on the row count and on the sum of
`jointly_called`.

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

Source revision: 8a17605021950429e3d32bfc572390122880bd8b; `src` tree f0076f48eb7c93058db828c1afba6ab14d367086. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: f6e3ba6cb2d19fe6f303ac950978b77672a86420bfcc9482047738a0065cb948. CPU affinity: 0-3. Host load before the run: 7.29, 6.68, 5.44; after: 2.16, 4.85, 4.92. Empty-run overhead: 39.4 MiB at one thread and 39.6 MiB at four. Every observation is in [`benchmark_somalier_cohort_scaling_runs.csv`](benchmark_somalier_cohort_scaling_runs.csv).

## Inputs as measured

`evidence_kib` and `sketch_kib` are the Parquet file sizes. The decoded evidence is about 120 bytes per row.

| source    | samples |   sites | evidence_rows |   pairs | evidence_kib | sketch_kib |
|:----------|--------:|--------:|--------------:|--------:|-------------:|-----------:|
| synthetic |     125 |  17,000 |     2,125,000 |   7,750 |      6,326.6 |       84.5 |
| synthetic |     250 |  17,000 |     4,250,000 |  31,125 |     12,637.5 |      163.3 |
| synthetic |     500 |  17,000 |     8,500,000 | 124,750 |     25,116.3 |      321.5 |
| synthetic |     250 |  34,000 |     8,500,000 |  31,125 |     36,766.0 |      290.5 |
| synthetic |     250 |  68,000 |    17,000,000 |  31,125 |     74,662.7 |      544.1 |
| synthetic |     125 |   4,250 |       531,250 |   7,750 |      1,130.5 |       32.8 |
| synthetic |     250 |   8,500 |     2,125,000 |  31,125 |      4,199.5 |       96.0 |
| real      |       3 | 108,757 |       326,271 |       3 |      3,527.4 |      113.1 |

## `sketches`: count evidence to sketches, written to Parquet

| dimension | scale | threads | samples |  sites | input_rows | output_rows | output_kib | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | peak_buffer_mib | spill_bytes | budget_mib | within_budget |
|:----------|------:|--------:|--------:|-------:|-----------:|------------:|-----------:|---------------:|------------:|------------:|---------------:|------------:|----------------:|------------:|-----------:|:--------------|
| samples   |     1 |       1 |     125 |  17000 |    2125000 |         125 |       86.7 |          0.292 |       0.290 |       0.319 |           62.5 |        62.5 |            41.8 |           0 |      121.6 | TRUE          |
| samples   |     2 |       1 |     250 |  17000 |    4250000 |         250 |      169.4 |          0.552 |       0.552 |       0.565 |           68.1 |        68.8 |            42.2 |           0 |      132.5 | TRUE          |
| samples   |     4 |       1 |     500 |  17000 |    8500000 |         500 |      333.9 |          1.092 |       1.081 |       1.095 |           82.6 |        82.8 |            52.8 |           0 |      154.1 | TRUE          |
| samples   |     1 |       4 |     125 |  17000 |    2125000 |         125 |       84.2 |          0.148 |       0.142 |       0.150 |           75.7 |        75.8 |            86.5 |           0 |      257.5 | TRUE          |
| samples   |     2 |       4 |     250 |  17000 |    4250000 |         250 |      159.1 |          0.251 |       0.250 |       0.254 |           81.9 |        81.9 |            86.8 |           0 |      277.4 | TRUE          |
| samples   |     4 |       4 |     500 |  17000 |    8500000 |         500 |      295.3 |          0.480 |       0.470 |       0.486 |           95.5 |        96.2 |            87.4 |           0 |      317.4 | TRUE          |
| sites     |     1 |       1 |     250 |  17000 |    4250000 |         250 |      169.4 |          0.552 |       0.552 |       0.565 |           68.1 |        68.8 |            42.2 |           0 |      132.5 | TRUE          |
| sites     |     2 |       1 |     250 |  34000 |    8500000 |         250 |      283.5 |          1.127 |       1.125 |       1.159 |           89.1 |        89.3 |            68.2 |           0 |      182.9 | TRUE          |
| sites     |     4 |       1 |     250 |  68000 |   17000000 |         250 |      542.2 |          2.366 |       2.278 |       2.443 |          129.0 |       129.7 |           124.3 |           0 |      283.8 | TRUE          |
| sites     |     1 |       4 |     250 |  17000 |    4250000 |         250 |      159.1 |          0.251 |       0.250 |       0.254 |           81.9 |        81.9 |            86.8 |           0 |      277.4 | TRUE          |
| sites     |     2 |       4 |     250 |  34000 |    8500000 |         250 |      292.0 |          0.482 |       0.477 |       0.495 |          107.5 |       108.7 |           109.3 |           0 |      346.2 | TRUE          |
| sites     |     4 |       4 |     250 |  68000 |   17000000 |         250 |      554.8 |          0.965 |       0.953 |       0.978 |          153.2 |       153.5 |           143.0 |           0 |      483.5 | TRUE          |
| joint     |     1 |       1 |     125 |   4250 |     531250 |         125 |       34.5 |          0.080 |       0.079 |       0.083 |           52.7 |        52.8 |            29.8 |           0 |       91.8 | TRUE          |
| joint     |     2 |       1 |     250 |   8500 |    2125000 |         250 |      100.9 |          0.284 |       0.282 |       0.325 |           58.7 |        58.8 |            35.3 |           0 |      107.2 | TRUE          |
| joint     |     4 |       1 |     500 |  17000 |    8500000 |         500 |      333.9 |          1.092 |       1.081 |       1.095 |           82.6 |        82.8 |            52.8 |           0 |      154.1 | TRUE          |
| joint     |     1 |       4 |     125 |   4250 |     531250 |         125 |       33.0 |          0.046 |       0.043 |       0.046 |           61.1 |        61.5 |            68.4 |           0 |      220.8 | TRUE          |
| joint     |     2 |       4 |     250 |   8500 |    2125000 |         250 |       95.2 |          0.140 |       0.137 |       0.146 |           69.9 |        70.6 |            80.0 |           0 |      243.0 | TRUE          |
| joint     |     4 |       4 |     500 |  17000 |    8500000 |         500 |      295.3 |          0.480 |       0.470 |       0.486 |           95.5 |        96.2 |            87.4 |           0 |      317.4 | TRUE          |
| real      |     1 |       1 |       3 | 108757 |     326271 |           3 |      113.1 |          0.166 |       0.163 |       0.172 |          129.0 |       129.3 |           146.7 |           0 |      269.9 | TRUE          |
| real      |     1 |       4 |       3 | 108757 |     326271 |           3 |      113.0 |          0.094 |       0.092 |       0.098 |          136.2 |       136.4 |           175.7 |           0 |      398.0 | TRUE          |

## `pairs_reduce`: all pairs, reduced

|     | dimension | scale | threads | samples |  sites | input_rows | output_rows | output_kib | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | peak_buffer_mib | spill_bytes | budget_mib | within_budget |
|:----|:----------|------:|--------:|--------:|-------:|-----------:|------------:|-----------:|---------------:|------------:|------------:|---------------:|------------:|----------------:|------------:|-----------:|:--------------|
| 21  | samples   |     1 |       1 |     125 |  17000 |        125 |        7750 |        0.1 |          0.068 |       0.067 |       0.070 |           86.5 |        86.7 |            18.3 |           0 |      142.7 | TRUE          |
| 22  | samples   |     2 |       1 |     250 |  17000 |        250 |       31125 |        0.1 |          0.191 |       0.187 |       0.199 |          102.3 |       102.5 |            24.7 |           0 |      157.0 | TRUE          |
| 23  | samples   |     4 |       1 |     500 |  17000 |        500 |      124750 |        0.1 |          0.720 |       0.690 |       0.750 |          127.4 |       127.7 |            36.8 |           0 |      185.5 | TRUE          |
| 24  | samples   |     1 |       4 |     125 |  17000 |        125 |        7750 |        0.1 |          0.068 |       0.066 |       0.069 |           94.3 |        94.7 |            18.3 |           0 |      393.3 | TRUE          |
| 25  | samples   |     2 |       4 |     250 |  17000 |        250 |       31125 |        0.1 |          0.202 |       0.197 |       0.208 |          114.8 |       119.2 |            24.8 |           0 |      424.3 | TRUE          |
| 26  | samples   |     4 |       4 |     500 |  17000 |        500 |      124750 |        0.1 |          0.772 |       0.710 |       0.776 |          128.9 |       135.5 |            36.9 |           0 |      452.8 | TRUE          |
| 27  | sites     |     1 |       1 |     250 |  17000 |        250 |       31125 |        0.1 |          0.191 |       0.187 |       0.199 |          102.3 |       102.5 |            24.7 |           0 |      157.0 | TRUE          |
| 28  | sites     |     2 |       1 |     250 |  34000 |        250 |       31125 |        0.1 |          0.373 |       0.365 |       0.391 |          160.7 |       161.0 |            36.8 |           0 |      259.2 | TRUE          |
| 29  | sites     |     4 |       1 |     250 |  68000 |        250 |       31125 |        0.1 |          0.996 |       0.995 |       1.030 |          205.4 |       205.6 |            61.4 |           0 |      463.2 | TRUE          |
| 30  | sites     |     1 |       4 |     250 |  17000 |        250 |       31125 |        0.1 |          0.202 |       0.197 |       0.208 |          114.8 |       119.2 |            24.8 |           0 |      424.3 | TRUE          |
| 31  | sites     |     2 |       4 |     250 |  34000 |        250 |       31125 |        0.1 |          0.393 |       0.389 |       0.393 |          184.2 |       199.6 |            36.8 |           0 |      750.9 | TRUE          |
| 32  | sites     |     4 |       4 |     250 |  68000 |        250 |       31125 |        0.1 |          0.990 |       0.990 |       0.992 |          216.7 |       221.5 |            61.4 |           0 |     1403.0 | TRUE          |
| 33  | joint     |     1 |       1 |     125 |   4250 |        125 |        7750 |        0.1 |          0.022 |       0.021 |       0.022 |           59.4 |        59.5 |            13.1 |           0 |       76.5 | TRUE          |
| 34  | joint     |     2 |       1 |     250 |   8500 |        250 |       31125 |        0.1 |          0.112 |       0.109 |       0.115 |           77.8 |        77.8 |            18.3 |           0 |      105.9 | TRUE          |
| 35  | joint     |     4 |       1 |     500 |  17000 |        500 |      124750 |        0.1 |          0.720 |       0.690 |       0.750 |          127.4 |       127.7 |            36.8 |           0 |      185.5 | TRUE          |
| 36  | joint     |     1 |       4 |     125 |   4250 |        125 |        7750 |        0.1 |          0.023 |       0.022 |       0.023 |           57.6 |        57.7 |            13.2 |           0 |      171.2 | TRUE          |
| 37  | joint     |     2 |       4 |     250 |   8500 |        250 |       31125 |        0.1 |          0.106 |       0.106 |       0.112 |           77.4 |        77.6 |            18.3 |           0 |      261.0 | TRUE          |
| 38  | joint     |     4 |       4 |     500 |  17000 |        500 |      124750 |        0.1 |          0.772 |       0.710 |       0.776 |          128.9 |       135.5 |            36.9 |           0 |      452.8 | TRUE          |
| 39  | real      |     1 |       1 |       3 | 108757 |          3 |           3 |        0.1 |          0.004 |       0.004 |       0.004 |           43.7 |        44.0 |             8.8 |           0 |       53.5 | TRUE          |
| 40  | real      |     1 |       4 |       3 | 108757 |          3 |           3 |        0.1 |          0.004 |       0.004 |       0.004 |           44.5 |        44.6 |             9.3 |           0 |       87.4 | TRUE          |

## `pairs_copy`: all pairs, written to Parquet

|     | dimension | scale | threads | samples |  sites | input_rows | output_rows | output_kib | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | peak_buffer_mib | spill_bytes | budget_mib | within_budget |
|:----|:----------|------:|--------:|--------:|-------:|-----------:|------------:|-----------:|---------------:|------------:|------------:|---------------:|------------:|----------------:|------------:|-----------:|:--------------|
| 41  | samples   |     1 |       1 |     125 |  17000 |        125 |        7750 |      111.9 |          0.073 |       0.072 |       0.075 |           90.4 |        90.5 |            20.2 |           0 |      151.6 | TRUE          |
| 42  | samples   |     2 |       1 |     250 |  17000 |        250 |       31125 |      366.1 |          0.216 |       0.216 |       0.229 |          108.1 |       108.3 |            34.7 |           0 |      192.6 | TRUE          |
| 43  | samples   |     4 |       1 |     500 |  17000 |        500 |      124750 |     1425.3 |          0.793 |       0.786 |       0.823 |          196.0 |       201.1 |            78.1 |           0 |      326.1 | TRUE          |
| 44  | samples   |     1 |       4 |     125 |  17000 |        125 |        7750 |      111.9 |          0.073 |       0.070 |       0.075 |           95.6 |        96.4 |            20.3 |           0 |      402.1 | TRUE          |
| 45  | samples   |     2 |       4 |     250 |  17000 |        250 |       31125 |      366.1 |          0.210 |       0.203 |       0.219 |          109.1 |       117.1 |            34.8 |           0 |      459.9 | TRUE          |
| 46  | samples   |     4 |       4 |     500 |  17000 |        500 |      124750 |     1425.3 |          0.770 |       0.759 |       0.779 |          199.9 |       201.1 |            78.2 |           0 |      595.6 | TRUE          |
| 47  | sites     |     1 |       1 |     250 |  17000 |        250 |       31125 |      366.1 |          0.216 |       0.216 |       0.229 |          108.1 |       108.3 |            34.7 |           0 |      192.6 | TRUE          |
| 48  | sites     |     2 |       1 |     250 |  34000 |        250 |       31125 |      355.1 |          0.413 |       0.402 |       0.419 |          150.1 |       150.3 |            45.2 |           0 |      294.8 | TRUE          |
| 49  | sites     |     4 |       1 |     250 |  68000 |        250 |       31125 |      351.3 |          1.008 |       1.007 |       1.037 |          213.7 |       213.9 |            66.9 |           0 |      498.8 | TRUE          |
| 50  | sites     |     1 |       4 |     250 |  17000 |        250 |       31125 |      366.1 |          0.210 |       0.203 |       0.219 |          109.1 |       117.1 |            34.8 |           0 |      459.9 | TRUE          |
| 51  | sites     |     2 |       4 |     250 |  34000 |        250 |       31125 |      355.1 |          0.392 |       0.391 |       0.421 |          197.3 |       200.5 |            45.3 |           0 |      786.6 | TRUE          |
| 52  | sites     |     4 |       4 |     250 |  68000 |        250 |       31125 |      351.3 |          1.008 |       1.001 |       1.017 |          216.0 |       218.3 |            66.8 |           0 |     1438.6 | TRUE          |
| 53  | joint     |     1 |       1 |     125 |   4250 |        125 |        7750 |      116.3 |          0.024 |       0.024 |       0.024 |           60.7 |        60.7 |            15.1 |           0 |       85.4 | TRUE          |
| 54  | joint     |     2 |       1 |     250 |   8500 |        250 |       31125 |      356.2 |          0.116 |       0.116 |       0.121 |           87.7 |        87.7 |            40.5 |           0 |      141.5 | TRUE          |
| 55  | joint     |     4 |       1 |     500 |  17000 |        500 |      124750 |     1425.3 |          0.793 |       0.786 |       0.823 |          196.0 |       201.1 |            78.1 |           0 |      326.1 | TRUE          |
| 56  | joint     |     1 |       4 |     125 |   4250 |        125 |        7750 |      116.3 |          0.026 |       0.026 |       0.027 |           60.9 |        61.3 |            15.2 |           0 |      180.1 | TRUE          |
| 57  | joint     |     2 |       4 |     250 |   8500 |        250 |       31125 |      356.2 |          0.118 |       0.117 |       0.120 |           97.0 |        98.1 |            41.3 |           0 |      296.6 | TRUE          |
| 58  | joint     |     4 |       4 |     500 |  17000 |        500 |      124750 |     1425.3 |          0.770 |       0.759 |       0.779 |          199.9 |       201.1 |            78.2 |           0 |      595.6 | TRUE          |
| 59  | real      |     1 |       1 |       3 | 108757 |          3 |           3 |        4.6 |          0.004 |       0.004 |       0.004 |           43.6 |        43.7 |             8.7 |           0 |       53.5 | TRUE          |
| 60  | real      |     1 |       4 |       3 | 108757 |          3 |           3 |        4.6 |          0.004 |       0.004 |       0.004 |           45.6 |        45.8 |             9.2 |           0 |       87.4 | TRUE          |

## `pairs_scalar`: all pairs through the per-pair function, reduced

|     | dimension | scale | threads | samples |  sites | input_rows | output_rows | output_kib | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | peak_buffer_mib | spill_bytes | budget_mib | within_budget |
|:----|:----------|------:|--------:|--------:|-------:|-----------:|------------:|-----------:|---------------:|------------:|------------:|---------------:|------------:|----------------:|------------:|-----------:|:--------------|
| 61  | samples   |     1 |       1 |     125 |  17000 |        125 |        7750 |        0.1 |          0.155 |       0.155 |       0.158 |           86.7 |        86.7 |            18.3 |           0 |      142.7 | TRUE          |
| 62  | samples   |     2 |       1 |     250 |  17000 |        250 |       31125 |        0.1 |          0.569 |       0.562 |       0.569 |           99.8 |       100.0 |            27.4 |           0 |      157.0 | TRUE          |
| 63  | samples   |     4 |       1 |     500 |  17000 |        500 |      124750 |        0.1 |          2.158 |       2.152 |       2.254 |          112.8 |       112.8 |            45.4 |           0 |      185.5 | TRUE          |
| 64  | samples   |     1 |       4 |     125 |  17000 |        125 |        7750 |        0.1 |          0.157 |       0.156 |       0.159 |           93.9 |        94.8 |            18.4 |           0 |      393.3 | TRUE          |
| 65  | samples   |     2 |       4 |     250 |  17000 |        250 |       31125 |        0.1 |          0.565 |       0.564 |       0.574 |           97.4 |       102.8 |            27.7 |           0 |      424.3 | TRUE          |
| 66  | samples   |     4 |       4 |     500 |  17000 |        500 |      124750 |        0.1 |          2.182 |       2.180 |       2.186 |          121.5 |       121.7 |            45.7 |           0 |      452.8 | TRUE          |
| 67  | sites     |     1 |       1 |     250 |  17000 |        250 |       31125 |        0.1 |          0.569 |       0.562 |       0.569 |           99.8 |       100.0 |            27.4 |           0 |      157.0 | TRUE          |
| 68  | sites     |     2 |       1 |     250 |  34000 |        250 |       31125 |        0.1 |          1.090 |       1.082 |       1.138 |          159.7 |       163.6 |            45.4 |           0 |      259.2 | TRUE          |
| 69  | sites     |     4 |       1 |     250 |  68000 |        250 |       31125 |        0.1 |          2.435 |       2.433 |       2.498 |          202.7 |       202.8 |            81.7 |           0 |      463.2 | TRUE          |
| 70  | sites     |     1 |       4 |     250 |  17000 |        250 |       31125 |        0.1 |          0.565 |       0.564 |       0.574 |           97.4 |       102.8 |            27.7 |           0 |      424.3 | TRUE          |
| 71  | sites     |     2 |       4 |     250 |  34000 |        250 |       31125 |        0.1 |          1.112 |       1.111 |       1.116 |          147.9 |       154.8 |            45.7 |           0 |      750.9 | TRUE          |
| 72  | sites     |     4 |       4 |     250 |  68000 |        250 |       31125 |        0.1 |          2.662 |       2.639 |       2.777 |          209.4 |       210.0 |            81.5 |           0 |     1403.0 | TRUE          |
| 73  | joint     |     1 |       1 |     125 |   4250 |        125 |        7750 |        0.1 |          0.047 |       0.045 |       0.047 |           54.5 |        54.5 |            10.9 |           0 |       76.5 | TRUE          |
| 74  | joint     |     2 |       1 |     250 |   8500 |        250 |       31125 |        0.1 |          0.291 |       0.290 |       0.299 |           73.2 |        73.3 |            18.4 |           0 |      105.9 | TRUE          |
| 75  | joint     |     4 |       1 |     500 |  17000 |        500 |      124750 |        0.1 |          2.158 |       2.152 |       2.254 |          112.8 |       112.8 |            45.4 |           0 |      185.5 | TRUE          |
| 76  | joint     |     1 |       4 |     125 |   4250 |        125 |        7750 |        0.1 |          0.047 |       0.044 |       0.051 |           55.5 |        55.8 |            11.0 |           0 |      171.2 | TRUE          |
| 77  | joint     |     2 |       4 |     250 |   8500 |        250 |       31125 |        0.1 |          0.295 |       0.295 |       0.304 |           75.4 |        82.6 |            18.4 |           0 |      261.0 | TRUE          |
| 78  | joint     |     4 |       4 |     500 |  17000 |        500 |      124750 |        0.1 |          2.182 |       2.180 |       2.186 |          121.5 |       121.7 |            45.7 |           0 |      452.8 | TRUE          |
| 79  | real      |     1 |       1 |       3 | 108757 |          3 |           3 |        0.1 |          0.003 |       0.003 |       0.003 |           42.6 |        42.8 |             7.4 |           0 |       53.5 | TRUE          |
| 80  | real      |     1 |       4 |       3 | 108757 |          3 |           3 |        0.1 |          0.003 |       0.003 |       0.004 |           42.6 |        43.1 |             7.4 |           0 |       87.4 | TRUE          |

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
|     125 |  4,250 |   7,750 |       1 |  0.081 |             0.047 |           0.022 |         3.764 |
|     250 |  8,500 |  31,125 |       1 |  0.574 |             0.291 |           0.112 |         5.135 |
|     125 | 17,000 |   7,750 |       1 |  0.298 |             0.155 |           0.068 |         4.393 |
|     250 | 17,000 |  31,125 |       1 |  1.121 |             0.569 |           0.191 |         5.859 |
|     500 | 17,000 | 124,750 |       1 |  4.399 |             2.158 |           0.720 |         6.111 |
|     250 | 34,000 |  31,125 |       1 |  2.201 |             1.090 |           0.373 |         5.907 |
|     250 | 68,000 |  31,125 |       1 |  4.633 |             2.435 |           0.996 |         4.652 |
|     125 |  4,250 |   7,750 |       4 |  0.081 |             0.047 |           0.023 |         3.560 |
|     250 |  8,500 |  31,125 |       4 |  0.579 |             0.295 |           0.106 |         5.449 |
|     125 | 17,000 |   7,750 |       4 |  0.295 |             0.157 |           0.068 |         4.356 |
|     250 | 17,000 |  31,125 |       4 |  1.125 |             0.565 |           0.202 |         5.578 |
|     500 | 17,000 | 124,750 |       4 |  4.413 |             2.182 |           0.772 |         5.714 |
|     250 | 34,000 |  31,125 |       4 |  2.204 |             1.112 |           0.393 |         5.603 |
|     250 | 68,000 |  31,125 |       4 |  4.834 |             2.662 |           0.990 |         4.882 |

## Time exponents and memory ratios

Time exponents per doubling, against the work of each stage, and peak-RSS ratios to 1×:

| stage        | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | timing_verdict     | exponent_above_1.25 |
|:-------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|:-------------------|:--------------------|
| sketches     | samples   |       1 |            0.919 |            0.984 |        1.090 |        1.322 | none: 1x under 5 s | FALSE               |
| sketches     | samples   |       4 |            0.765 |            0.937 |        1.082 |        1.262 | none: 1x under 5 s | FALSE               |
| sketches     | sites     |       1 |            1.029 |            1.070 |        1.308 |        1.894 | none: 1x under 5 s | FALSE               |
| sketches     | sites     |       4 |            0.943 |            1.001 |        1.313 |        1.871 | none: 1x under 5 s | FALSE               |
| sketches     | joint     |       1 |            0.911 |            0.970 |        1.114 |        1.567 | none: 1x under 5 s | FALSE               |
| sketches     | joint     |       4 |            0.805 |            0.890 |        1.144 |        1.563 | none: 1x under 5 s | FALSE               |
| pairs_reduce | samples   |       1 |            0.746 |            0.954 |        1.183 |        1.473 | none: 1x under 5 s | FALSE               |
| pairs_reduce | samples   |       4 |            0.786 |            0.967 |        1.217 |        1.367 | none: 1x under 5 s | FALSE               |
| pairs_reduce | sites     |       1 |            0.961 |            1.418 |        1.571 |        2.008 | none: 1x under 5 s | TRUE                |
| pairs_reduce | sites     |       4 |            0.963 |            1.332 |        1.605 |        1.888 | none: 1x under 5 s | TRUE                |
| pairs_reduce | joint     |       1 |            0.790 |            0.895 |        1.310 |        2.145 | none: 1x under 5 s | FALSE               |
| pairs_reduce | joint     |       4 |            0.739 |            0.953 |        1.344 |        2.238 | none: 1x under 5 s | FALSE               |
| pairs_copy   | samples   |       1 |            0.779 |            0.936 |        1.196 |        2.168 | none: 1x under 5 s | FALSE               |
| pairs_copy   | samples   |       4 |            0.756 |            0.937 |        1.141 |        2.091 | none: 1x under 5 s | FALSE               |
| pairs_copy   | sites     |       1 |            0.931 |            1.287 |        1.389 |        1.977 | none: 1x under 5 s | TRUE                |
| pairs_copy   | sites     |       4 |            0.903 |            1.362 |        1.808 |        1.980 | none: 1x under 5 s | TRUE                |
| pairs_copy   | joint     |       1 |            0.759 |            0.923 |        1.445 |        3.229 | none: 1x under 5 s | FALSE               |
| pairs_copy   | joint     |       4 |            0.720 |            0.900 |        1.593 |        3.282 | none: 1x under 5 s | FALSE               |
| pairs_scalar | samples   |       1 |            0.933 |            0.961 |        1.151 |        1.301 | none: 1x under 5 s | FALSE               |
| pairs_scalar | samples   |       4 |            0.922 |            0.973 |        1.037 |        1.294 | none: 1x under 5 s | FALSE               |
| pairs_scalar | sites     |       1 |            0.938 |            1.160 |        1.600 |        2.031 | none: 1x under 5 s | FALSE               |
| pairs_scalar | sites     |       4 |            0.975 |            1.260 |        1.518 |        2.150 | none: 1x under 5 s | TRUE                |
| pairs_scalar | joint     |       1 |            0.880 |            0.962 |        1.343 |        2.070 | none: 1x under 5 s | FALSE               |
| pairs_scalar | joint     |       4 |            0.884 |            0.962 |        1.359 |        2.189 | none: 1x under 5 s | FALSE               |

## Findings

Every observation is within its declared budget and none spills; the render
enforces both. The largest observation uses
82% of its ceiling.

The largest all-pairs cell (500 samples, 17,000 sites, 124,750 pairs) takes
0.79 s on one thread and
0.77 s on four when written to Parquet, at
196 and
199.9 MiB peak RSS. The Parquet output is
1425.3 KiB.

On the real counts, one thread builds three sketches of 108,757 sites in
0.166 s at
129 MiB peak RSS, and writes the three pairs in
0.004 s at
43.6 MiB.

For 124,750 pairs of 17,000 sites on one thread, the per-pair function took
4.4 s at revision `a5f3b47a` and takes
2.16 s now; the all-pairs macro takes
0.72 s, 6.1 times faster than before.
The pair stages still grow with pairs times sites, which is the size of the
comparison itself. They do not use more than one thread:
0.72 s on one thread and
0.77 s on four. DuckDB splits a
scan by row groups, and a relation of a few hundred sketches is one row group.

24 of the 24 series have a one-thread 1× run under 5 seconds and get a memory verdict only.
5 series have a time exponent above 1.25 at either doubling.
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
