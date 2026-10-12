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

The all-pairs macro checks the contents of each sketch once and compares the
pairs on DuckDB’s threads. The per-pair function checks both sketches of every
pair. [Before and after](#before-and-after) compares both with the observations
of revision `a5f3b47a`. [Change from the previous render](#change-from-the-previous-render)
compares this render with the previous one, which ran the macro’s comparisons
in one thread.

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
- The macro holds the checked sketches as one list. The pairs are numbered and
  split into at most 65,536 blocks of consecutive pairs. Each thread that
  compares pairs copies the list once for each chunk of up to 2,048 blocks
  that it scans, and reads both sketches of a pair in place from that copy
  (see [Budgets](#budgets)). Across all threads, a query makes at most 32
  such copies, whatever `N` is.
- Caps: `max_sites` bounds `S` per sketch (default 1,000,000; hard limit
  100,000,000). `max_pairs` bounds `P` (default 50,000,000): a relation with
  more pairs is an error before any pair is compared. The caller chooses the
  pair relation; a `pairs_table` join gives a linear selection instead of all
  pairs.

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
thread and at four threads, in an order that alternates sizes. The render runs
at `nice -n -15` under `taskset`, and every CLI process inherits both. This
host mixes performance and efficiency cores, so the header below records the
CPU affinity, the kind of core and the nice value; every row of the runs file
records the affinity. Input Parquet files are read from a warm page cache.

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

The budgets were first declared before the first measurement, at revision
`a5f3b47a`. That run stopped at its own gate: the real cell of three pairs
used 42.5 to 44.1 MiB against a ceiling of 41.5 MiB. The first model scaled
the join output by the number of pairs, so three pairs were budgeted three
rows. DuckDB allocates whole 2,048-row vectors for the 30-field result in each
operator, whatever the row count. The model gained that fixed term, and the
sketch relation that the macro then materialized. The observations of the
first run are kept in
[`benchmark_somalier_cohort_scaling_runs_a5f3b47a.csv`](benchmark_somalier_cohort_scaling_runs_a5f3b47a.csv).

The all-pairs macro now compares pairs on DuckDB’s threads, with a new plan.
The terms of `pairs_reduce` and `pairs_copy` below follow that plan. They come
from the macro and the DuckDB 1.5 source, and they were declared before this
render. The terms of `sketches` and `pairs_scalar` did not change.

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
- **`pairs_reduce`.** The macro collects the checked sketches into one list in
  `sample_id` order, in a materialized one-row relation. A cross product joins
  the block rows with that list, and the pair indexes of each block are
  unnested above the join; the unnest passes the list on by reference. The
  pair function reads both sketches of a pair in place from the list. DuckDB
  1.5 frees an operator’s state only when the query ends, so three copies of
  the list stay live while the pairs are compared, `3 × N` sketches:
  - the materialized relation that holds the list;
  - the copy of the list that the cross product keeps;
  - the segments of the ordered `list` aggregate: after its sort, the
    aggregate copies each sketch into segments of its operator’s arena.

  Each thread that compares pairs adds:
  - two copies of the list’s mask words, `2 × N` sketches. The cross product
    reads its copy of the list again for each chunk of up to 2,048 block rows.
    The mask words fill many 2,048-value vectors, so the reader copies them
    into one new vector. It allocates the new copy before it frees the
    previous one;
  - the result vectors, fixed per thread: four operators hold one 2,048-row
    vector of the 30-field result at 16 bytes per field,
    `4 × 2,048 × 30 × 16` bytes.

  The model counts every thread, also when fewer threads compare pairs. The
  cohort stage, which builds the list, ends before the first pair is compared.
  It holds at most four lists at once, at the end of the ordered aggregate:
  the sorted sketches, the sorted rows that it reads back, the list segments
  and the result list. That is less than the pair stage holds with one thread.
- **`pairs_copy`.** The live state of `pairs_reduce`, and the Parquet writer’s
  row group: `min(P, threads × 122,880)` rows of 400 bytes.
- **`pairs_scalar`.** The per-pair query did not change, and its terms are the
  ones declared at `a5f3b47a`: `6 × N` sketches for the sketch relation and
  both sides of the join; one join output vector per thread, which carries
  both sketches of each row, `min(P, threads × 2,048) × 2` sketches; and the
  result vectors of `pairs_reduce`.

A time exponent is the log of the median time ratio divided by the log of the
work ratio, per doubling. The work is `N × S` for `sketches` and `P × S` for
the pair stages, so a pair stage that is quadratic in samples and linear in
pairs has exponent 1. An exponent above 1.25 needs an investigation. A series
whose one-thread 1× run is under 5 seconds gets a memory verdict only.

``` sh
nice -n -15 taskset -c 0-3 Rscript -e 'rmarkdown::render("benchmarks/benchmark_somalier_cohort_scaling.Rmd")'
```

Source revision: 57c27de3446a8482d03eae179ce490294f06f4d4; `src` tree fda943b6dcbb090d387fc4a934091e6d17589e62. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: ea551be30ca37c6e5a672b7cf7c2c48396cc762c3611c2c4a005b6a5ebc5dbdc. CPU affinity: 0-3 (2 physical cores, 2 CPUs per core, 4.8 GHz maximum). Nice value: -15. Host load before the run: 1.39, 1.63, 1.96; after: 1.83, 1.66, 1.92. Empty-run overhead: 39.8 MiB at one thread and 39.8 MiB at four. Every observation is in [`benchmark_somalier_cohort_scaling_runs.csv`](benchmark_somalier_cohort_scaling_runs.csv).

## Inputs as measured

`evidence_kib` and `sketch_kib` are the Parquet file sizes. The decoded evidence is about 120 bytes per row. `macro_differing_rows` counts the rows, over every field, that differ between `duckhts_somalier_relatedness_all_pairs` and the per-pair function.

| source    | samples |   sites | evidence_rows |   pairs | evidence_kib | sketch_kib | macro_differing_rows |
|:----------|--------:|--------:|--------------:|--------:|-------------:|-----------:|---------------------:|
| synthetic |     125 |  17,000 |     2,125,000 |   7,750 |      6,326.6 |       83.9 |                    0 |
| synthetic |     250 |  17,000 |     4,250,000 |  31,125 |     12,637.5 |      161.4 |                    0 |
| synthetic |     500 |  17,000 |     8,500,000 | 124,750 |     25,116.3 |      318.0 |                    0 |
| synthetic |     250 |  34,000 |     8,500,000 |  31,125 |     36,766.0 |      290.4 |                    0 |
| synthetic |     250 |  68,000 |    17,000,000 |  31,125 |     74,662.7 |      542.8 |                    0 |
| synthetic |     125 |   4,250 |       531,250 |   7,750 |      1,130.5 |       32.4 |                    0 |
| synthetic |     250 |   8,500 |     2,125,000 |  31,125 |      4,199.5 |       95.2 |                    0 |
| real      |       3 | 108,757 |       326,271 |       3 |      3,527.4 |      113.1 |                    0 |

## `sketches`: count evidence to sketches, written to Parquet

| dimension | scale | threads | samples |  sites | input_rows | output_rows | output_kib | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | peak_buffer_mib | spill_bytes | budget_mib | within_budget |
|:----------|------:|--------:|--------:|-------:|-----------:|------------:|-----------:|---------------:|------------:|------------:|---------------:|------------:|----------------:|------------:|-----------:|:--------------|
| samples   |     1 |       1 |     125 |  17000 |    2125000 |         125 |       86.7 |          0.289 |       0.289 |       0.290 |           62.8 |        63.1 |            41.9 |           0 |      122.0 | TRUE          |
| samples   |     2 |       1 |     250 |  17000 |    4250000 |         250 |      169.4 |          0.555 |       0.553 |       0.555 |           68.6 |        68.6 |            42.2 |           0 |      132.8 | TRUE          |
| samples   |     4 |       1 |     500 |  17000 |    8500000 |         500 |      333.9 |          1.079 |       1.077 |       1.084 |           82.5 |        82.6 |            52.8 |           0 |      154.5 | TRUE          |
| samples   |     1 |       4 |     125 |  17000 |    2125000 |         125 |       84.3 |          0.141 |       0.139 |       0.147 |           76.5 |        76.6 |            86.5 |           0 |      257.7 | TRUE          |
| samples   |     2 |       4 |     250 |  17000 |    4250000 |         250 |      158.9 |          0.251 |       0.250 |       0.262 |           82.8 |        83.0 |            86.8 |           0 |      277.7 | TRUE          |
| samples   |     4 |       4 |     500 |  17000 |    8500000 |         500 |      308.8 |          0.471 |       0.462 |       0.481 |           95.5 |        96.2 |            87.4 |           0 |      317.6 | TRUE          |
| sites     |     1 |       1 |     250 |  17000 |    4250000 |         250 |      169.4 |          0.555 |       0.553 |       0.555 |           68.6 |        68.6 |            42.2 |           0 |      132.8 | TRUE          |
| sites     |     2 |       1 |     250 |  34000 |    8500000 |         250 |      283.5 |          1.124 |       1.117 |       1.127 |           88.3 |        88.4 |            68.2 |           0 |      183.3 | TRUE          |
| sites     |     4 |       1 |     250 |  68000 |   17000000 |         250 |      542.2 |          2.246 |       2.245 |       2.255 |          129.4 |       129.6 |           124.3 |           0 |      284.2 | TRUE          |
| sites     |     1 |       4 |     250 |  17000 |    4250000 |         250 |      158.9 |          0.251 |       0.250 |       0.262 |           82.8 |        83.0 |            86.8 |           0 |      277.7 | TRUE          |
| sites     |     2 |       4 |     250 |  34000 |    8500000 |         250 |      293.4 |          0.481 |       0.480 |       0.483 |          108.2 |       109.8 |           109.3 |           0 |      346.4 | TRUE          |
| sites     |     4 |       4 |     250 |  68000 |   17000000 |         250 |      554.1 |          0.944 |       0.938 |       0.948 |          151.9 |       154.1 |           143.0 |           0 |      483.8 | TRUE          |
| joint     |     1 |       1 |     125 |   4250 |     531250 |         125 |       34.5 |          0.079 |       0.078 |       0.079 |           53.0 |        53.1 |            29.8 |           0 |       92.2 | TRUE          |
| joint     |     2 |       1 |     250 |   8500 |    2125000 |         250 |      100.9 |          0.279 |       0.277 |       0.280 |           58.5 |        58.6 |            35.3 |           0 |      107.6 | TRUE          |
| joint     |     4 |       1 |     500 |  17000 |    8500000 |         500 |      333.9 |          1.079 |       1.077 |       1.084 |           82.5 |        82.6 |            52.8 |           0 |      154.5 | TRUE          |
| joint     |     1 |       4 |     125 |   4250 |     531250 |         125 |       33.0 |          0.043 |       0.043 |       0.045 |           61.1 |        61.5 |            68.3 |           0 |      221.0 | TRUE          |
| joint     |     2 |       4 |     250 |   8500 |    2125000 |         250 |       95.6 |          0.135 |       0.134 |       0.136 |           70.0 |        70.4 |            80.0 |           0 |      243.3 | TRUE          |
| joint     |     4 |       4 |     500 |  17000 |    8500000 |         500 |      308.8 |          0.471 |       0.462 |       0.481 |           95.5 |        96.2 |            87.4 |           0 |      317.6 | TRUE          |
| real      |     1 |       1 |       3 | 108757 |     326271 |           3 |      113.1 |          0.162 |       0.162 |       0.162 |          129.7 |       130.1 |           146.7 |           0 |      270.3 | TRUE          |
| real      |     1 |       4 |       3 | 108757 |     326271 |           3 |      113.1 |          0.097 |       0.089 |       0.098 |          136.8 |       136.9 |           175.7 |           0 |      398.3 | TRUE          |

## `pairs_reduce`: all pairs, reduced

|     | dimension | scale | threads | samples |  sites | input_rows | output_rows | output_kib | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | peak_buffer_mib | spill_bytes | budget_mib | within_budget |
|:----|:----------|------:|--------:|--------:|-------:|-----------:|------------:|-----------:|---------------:|------------:|------------:|---------------:|------------:|----------------:|------------:|-----------:|:--------------|
| 21  | samples   |     1 |       1 |     125 |  17000 |        125 |        7750 |        0.1 |          0.038 |       0.038 |       0.039 |           50.5 |        50.9 |            22.1 |           0 |       62.9 | TRUE          |
| 22  | samples   |     2 |       1 |     250 |  17000 |        250 |       31125 |        0.1 |          0.131 |       0.131 |       0.131 |           54.7 |        54.8 |            30.1 |           0 |       74.8 | TRUE          |
| 23  | samples   |     4 |       1 |     500 |  17000 |        500 |      124750 |        0.1 |          0.486 |       0.486 |       0.489 |           61.9 |        62.1 |            48.2 |           0 |       98.5 | TRUE          |
| 24  | samples   |     1 |       4 |     125 |  17000 |        125 |        7750 |        0.1 |          0.024 |       0.024 |       0.025 |           58.4 |        59.6 |            44.2 |           0 |      110.9 | TRUE          |
| 25  | samples   |     2 |       4 |     250 |  17000 |        250 |       31125 |        0.1 |          0.075 |       0.073 |       0.076 |           68.8 |        70.5 |            61.2 |           0 |      137.1 | TRUE          |
| 26  | samples   |     4 |       4 |     500 |  17000 |        500 |      124750 |        0.1 |          0.256 |       0.255 |       0.264 |           85.7 |        86.1 |            87.8 |           0 |      189.3 | TRUE          |
| 27  | sites     |     1 |       1 |     250 |  17000 |        250 |       31125 |        0.1 |          0.131 |       0.131 |       0.131 |           54.7 |        54.8 |            30.1 |           0 |       74.8 | TRUE          |
| 28  | sites     |     2 |       1 |     250 |  34000 |        250 |       31125 |        0.1 |          0.238 |       0.238 |       0.240 |           60.6 |        60.7 |            47.1 |           0 |       97.6 | TRUE          |
| 29  | sites     |     4 |       1 |     250 |  68000 |        250 |       31125 |        0.1 |          0.456 |       0.454 |       0.457 |           77.8 |        78.1 |            82.2 |           0 |      143.2 | TRUE          |
| 30  | sites     |     1 |       4 |     250 |  17000 |        250 |       31125 |        0.1 |          0.075 |       0.073 |       0.076 |           68.8 |        70.5 |            61.2 |           0 |      137.1 | TRUE          |
| 31  | sites     |     2 |       4 |     250 |  34000 |        250 |       31125 |        0.1 |          0.137 |       0.135 |       0.139 |           86.5 |        87.5 |            87.1 |           0 |      187.3 | TRUE          |
| 32  | sites     |     4 |       4 |     250 |  68000 |        250 |       31125 |        0.1 |          0.262 |       0.250 |       0.266 |          118.0 |       118.1 |           140.4 |           0 |      287.5 | TRUE          |
| 33  | joint     |     1 |       1 |     125 |   4250 |        125 |        7750 |        0.1 |          0.017 |       0.017 |       0.017 |           47.5 |        47.7 |            16.5 |           0 |       54.4 | TRUE          |
| 34  | joint     |     2 |       1 |     250 |   8500 |        250 |       31125 |        0.1 |          0.076 |       0.075 |       0.077 |           50.6 |        50.9 |            22.1 |           0 |       63.4 | TRUE          |
| 35  | joint     |     4 |       1 |     500 |  17000 |        500 |      124750 |        0.1 |          0.486 |       0.486 |       0.489 |           61.9 |        62.1 |            48.2 |           0 |       98.5 | TRUE          |
| 36  | joint     |     1 |       4 |     125 |   4250 |        125 |        7750 |        0.1 |          0.012 |       0.012 |       0.012 |           54.1 |        54.2 |            33.9 |           0 |       92.1 | TRUE          |
| 37  | joint     |     2 |       4 |     250 |   8500 |        250 |       31125 |        0.1 |          0.043 |       0.043 |       0.045 |           60.9 |        60.9 |            47.8 |           0 |      111.9 | TRUE          |
| 38  | joint     |     4 |       4 |     500 |  17000 |        500 |      124750 |        0.1 |          0.256 |       0.255 |       0.264 |           85.7 |        86.1 |            87.8 |           0 |      189.3 | TRUE          |
| 39  | real      |     1 |       1 |       3 | 108757 |          3 |           3 |        0.1 |          0.006 |       0.006 |       0.006 |           46.1 |        46.2 |            15.1 |           0 |       52.8 | TRUE          |
| 40  | real      |     1 |       4 |       3 | 108757 |          3 |           3 |        0.1 |          0.006 |       0.005 |       0.006 |           47.9 |        48.0 |            12.3 |           0 |       88.7 | TRUE          |

## `pairs_copy`: all pairs, written to Parquet

|     | dimension | scale | threads | samples |  sites | input_rows | output_rows | output_kib | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | peak_buffer_mib | spill_bytes | budget_mib | within_budget |
|:----|:----------|------:|--------:|--------:|-------:|-----------:|------------:|-----------:|---------------:|------------:|------------:|---------------:|------------:|----------------:|------------:|-----------:|:--------------|
| 41  | samples   |     1 |       1 |     125 |  17000 |        125 |        7750 |      103.4 |          0.042 |       0.042 |       0.043 |           53.1 |        53.3 |            22.6 |           0 |       71.8 | TRUE          |
| 42  | samples   |     2 |       1 |     250 |  17000 |        250 |       31125 |      330.8 |          0.145 |       0.145 |       0.145 |           74.0 |        74.0 |            32.2 |           0 |      110.4 | TRUE          |
| 43  | samples   |     4 |       1 |     500 |  17000 |        500 |      124750 |     1266.1 |          0.548 |       0.545 |       0.548 |          150.0 |       150.1 |           100.8 |           0 |      239.2 | TRUE          |
| 44  | samples   |     1 |       4 |     125 |  17000 |        125 |        7750 |      103.4 |          0.045 |       0.043 |       0.053 |           58.4 |        58.5 |            22.7 |           0 |      119.8 | TRUE          |
| 45  | samples   |     2 |       4 |     250 |  17000 |        250 |       31125 |      330.8 |          0.152 |       0.150 |       0.152 |           80.1 |        83.9 |            32.2 |           0 |      172.7 | TRUE          |
| 46  | samples   |     4 |       4 |     500 |  17000 |        500 |      124750 |     1266.1 |          0.547 |       0.543 |       0.570 |          153.1 |       157.2 |           101.0 |           0 |      332.1 | TRUE          |
| 47  | sites     |     1 |       1 |     250 |  17000 |        250 |       31125 |      330.8 |          0.145 |       0.145 |       0.145 |           74.0 |        74.0 |            32.2 |           0 |      110.4 | TRUE          |
| 48  | sites     |     2 |       1 |     250 |  34000 |        250 |       31125 |      319.8 |          0.253 |       0.252 |       0.254 |           81.7 |        81.9 |            49.2 |           0 |      133.2 | TRUE          |
| 49  | sites     |     4 |       1 |     250 |  68000 |        250 |       31125 |      311.3 |          0.470 |       0.469 |       0.474 |           94.7 |        94.8 |            84.2 |           0 |      178.8 | TRUE          |
| 50  | sites     |     1 |       4 |     250 |  17000 |        250 |       31125 |      330.8 |          0.152 |       0.150 |       0.152 |           80.1 |        83.9 |            32.2 |           0 |      172.7 | TRUE          |
| 51  | sites     |     2 |       4 |     250 |  34000 |        250 |       31125 |      319.8 |          0.255 |       0.253 |       0.271 |           94.1 |        95.4 |            49.3 |           0 |      222.9 | TRUE          |
| 52  | sites     |     4 |       4 |     250 |  68000 |        250 |       31125 |      311.3 |          0.473 |       0.468 |       0.485 |          110.2 |       118.7 |            84.4 |           0 |      323.2 | TRUE          |
| 53  | joint     |     1 |       1 |     125 |   4250 |        125 |        7750 |      107.1 |          0.020 |       0.020 |       0.021 |           50.3 |        50.3 |            16.6 |           0 |       63.3 | TRUE          |
| 54  | joint     |     2 |       1 |     250 |   8500 |        250 |       31125 |      317.8 |          0.092 |       0.090 |       0.092 |           70.2 |        70.4 |            45.0 |           0 |       99.0 | TRUE          |
| 55  | joint     |     4 |       1 |     500 |  17000 |        500 |      124750 |     1266.1 |          0.548 |       0.545 |       0.548 |          150.0 |       150.1 |           100.8 |           0 |      239.2 | TRUE          |
| 56  | joint     |     1 |       4 |     125 |   4250 |        125 |        7750 |      107.1 |          0.021 |       0.021 |       0.022 |           54.6 |        55.2 |            15.6 |           0 |      101.0 | TRUE          |
| 57  | joint     |     2 |       4 |     250 |   8500 |        250 |       31125 |      317.8 |          0.095 |       0.092 |       0.105 |           74.5 |        75.9 |            45.1 |           0 |      147.6 | TRUE          |
| 58  | joint     |     4 |       4 |     500 |  17000 |        500 |      124750 |     1266.1 |          0.547 |       0.543 |       0.570 |          153.1 |       157.2 |           101.0 |           0 |      332.1 | TRUE          |
| 59  | real      |     1 |       1 |       3 | 108757 |          3 |           3 |        4.6 |          0.006 |       0.006 |       0.006 |           45.9 |        46.1 |            15.0 |           0 |       52.8 | TRUE          |
| 60  | real      |     1 |       4 |       3 | 108757 |          3 |           3 |        4.6 |          0.006 |       0.006 |       0.006 |           47.8 |        48.3 |            13.4 |           0 |       88.7 | TRUE          |

## `pairs_scalar`: all pairs through the per-pair function, reduced

|     | dimension | scale | threads | samples |  sites | input_rows | output_rows | output_kib | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | peak_buffer_mib | spill_bytes | budget_mib | within_budget |
|:----|:----------|------:|--------:|--------:|-------:|-----------:|------------:|-----------:|---------------:|------------:|------------:|---------------:|------------:|----------------:|------------:|-----------:|:--------------|
| 61  | samples   |     1 |       1 |     125 |  17000 |        125 |        7750 |        0.1 |          0.156 |       0.155 |       0.159 |           86.2 |        86.2 |            18.4 |           0 |      143.1 | TRUE          |
| 62  | samples   |     2 |       1 |     250 |  17000 |        250 |       31125 |        0.1 |          0.558 |       0.554 |       0.561 |          100.3 |       100.4 |            27.4 |           0 |      157.4 | TRUE          |
| 63  | samples   |     4 |       1 |     500 |  17000 |        500 |      124750 |        0.1 |          2.147 |       2.134 |       2.156 |          113.0 |       113.0 |            45.7 |           0 |      185.9 | TRUE          |
| 64  | samples   |     1 |       4 |     125 |  17000 |        125 |        7750 |        0.1 |          0.158 |       0.153 |       0.158 |           94.9 |        95.5 |            18.4 |           0 |      393.5 | TRUE          |
| 65  | samples   |     2 |       4 |     250 |  17000 |        250 |       31125 |        0.1 |          0.558 |       0.556 |       0.559 |          100.9 |       100.9 |            27.7 |           0 |      424.6 | TRUE          |
| 66  | samples   |     4 |       4 |     500 |  17000 |        500 |      124750 |        0.1 |          2.147 |       2.146 |       2.153 |          113.3 |       121.8 |            45.5 |           0 |      453.1 | TRUE          |
| 67  | sites     |     1 |       1 |     250 |  17000 |        250 |       31125 |        0.1 |          0.558 |       0.554 |       0.561 |          100.3 |       100.4 |            27.4 |           0 |      157.4 | TRUE          |
| 68  | sites     |     2 |       1 |     250 |  34000 |        250 |       31125 |        0.1 |          1.080 |       1.078 |       1.082 |          164.4 |       164.7 |            45.4 |           0 |      259.6 | TRUE          |
| 69  | sites     |     4 |       1 |     250 |  68000 |        250 |       31125 |        0.1 |          2.409 |       2.401 |       2.411 |          203.0 |       203.1 |            81.7 |           0 |      463.6 | TRUE          |
| 70  | sites     |     1 |       4 |     250 |  17000 |        250 |       31125 |        0.1 |          0.558 |       0.556 |       0.559 |          100.9 |       100.9 |            27.7 |           0 |      424.6 | TRUE          |
| 71  | sites     |     2 |       4 |     250 |  34000 |        250 |       31125 |        0.1 |          1.098 |       1.090 |       1.099 |          148.1 |       153.6 |            45.7 |           0 |      751.2 | TRUE          |
| 72  | sites     |     4 |       4 |     250 |  68000 |        250 |       31125 |        0.1 |          2.607 |       2.410 |       2.624 |          205.9 |       210.0 |            81.5 |           0 |     1403.3 | TRUE          |
| 73  | joint     |     1 |       1 |     125 |   4250 |        125 |        7750 |        0.1 |          0.045 |       0.044 |       0.045 |           55.1 |        55.1 |            10.9 |           0 |       76.9 | TRUE          |
| 74  | joint     |     2 |       1 |     250 |   8500 |        250 |       31125 |        0.1 |          0.286 |       0.283 |       0.286 |           73.4 |        73.9 |            18.4 |           0 |      106.3 | TRUE          |
| 75  | joint     |     4 |       1 |     500 |  17000 |        500 |      124750 |        0.1 |          2.147 |       2.134 |       2.156 |          113.0 |       113.0 |            45.7 |           0 |      185.9 | TRUE          |
| 76  | joint     |     1 |       4 |     125 |   4250 |        125 |        7750 |        0.1 |          0.046 |       0.045 |       0.046 |           55.4 |        56.0 |            10.9 |           0 |      171.5 | TRUE          |
| 77  | joint     |     2 |       4 |     250 |   8500 |        250 |       31125 |        0.1 |          0.289 |       0.285 |       0.290 |           75.4 |        76.1 |            18.5 |           0 |      261.2 | TRUE          |
| 78  | joint     |     4 |       4 |     500 |  17000 |        500 |      124750 |        0.1 |          2.147 |       2.146 |       2.153 |          113.3 |       121.8 |            45.5 |           0 |      453.1 | TRUE          |
| 79  | real      |     1 |       1 |       3 | 108757 |          3 |           3 |        0.1 |          0.003 |       0.003 |       0.003 |           42.9 |        43.2 |             7.4 |           0 |       53.9 | TRUE          |
| 80  | real      |     1 |       4 |       3 | 108757 |          3 |           3 |        0.1 |          0.003 |       0.003 |       0.003 |           43.3 |        43.4 |             7.4 |           0 |       87.6 | TRUE          |

## Before and after

`before` is the `pairs_reduce` stage of revision `a5f3b47a`, which ran the
per-pair function over a self-join; its observations are in
[`benchmark_somalier_cohort_scaling_runs_a5f3b47a.csv`](benchmark_somalier_cohort_scaling_runs_a5f3b47a.csv).
`per_pair_function` and `all_pairs_macro` are the `pairs_scalar` and
`pairs_reduce` stages of this revision. All three produce the same rows from
the same sketches, with the same statement shape, host and DuckDB CLI.
`macro_speedup` is `before` divided by `all_pairs_macro`. Times are median
seconds of three fresh processes.

| samples |  sites |   pairs | threads | before | per_pair_function | all_pairs_macro | macro_speedup |
|--------:|-------:|--------:|--------:|-------:|------------------:|----------------:|--------------:|
|     125 |  4,250 |   7,750 |       1 |  0.081 |             0.045 |           0.017 |         4.783 |
|     250 |  8,500 |  31,125 |       1 |  0.574 |             0.286 |           0.076 |         7.524 |
|     125 | 17,000 |   7,750 |       1 |  0.298 |             0.156 |           0.038 |         7.748 |
|     250 | 17,000 |  31,125 |       1 |  1.121 |             0.558 |           0.131 |         8.530 |
|     500 | 17,000 | 124,750 |       1 |  4.399 |             2.147 |           0.486 |         9.045 |
|     250 | 34,000 |  31,125 |       1 |  2.201 |             1.080 |           0.238 |         9.229 |
|     250 | 68,000 |  31,125 |       1 |  4.633 |             2.409 |           0.456 |        10.149 |
|     125 |  4,250 |   7,750 |       4 |  0.081 |             0.046 |           0.012 |         6.905 |
|     250 |  8,500 |  31,125 |       4 |  0.579 |             0.289 |           0.043 |        13.365 |
|     125 | 17,000 |   7,750 |       4 |  0.295 |             0.158 |           0.024 |        12.050 |
|     250 | 17,000 |  31,125 |       4 |  1.125 |             0.558 |           0.075 |        14.931 |
|     500 | 17,000 | 124,750 |       4 |  4.413 |             2.147 |           0.256 |        17.219 |
|     250 | 34,000 |  31,125 |       4 |  2.204 |             1.098 |           0.137 |        16.143 |
|     250 | 68,000 |  31,125 |       4 |  4.834 |             2.607 |           0.262 |        18.416 |

The `a5f3b47a` run and this render both ran on CPUs 0-3 (2 physical cores, 2 CPUs per core, 4.8 GHz maximum). The two runs were made at different times on a shared
host, so the ratios also carry that noise.

## Change from the previous render

The previous render, of revision 68d865e2 (commit
bba5b238), measured the same cells with the same driver, inputs and
DuckDB CLI on the same host. The measured statements and the process driver
did not change. This render changes the budget terms of the macro stages, the
recorded host fields and the report text. Host load differs: the previous run
recorded
3.70, 2.05, 2.20 before and 4.97, 4.54, 3.24 after, against the
loads recorded above. Both renders ran on CPUs 0-3 (2 physical cores, 2 CPUs per core, 4.8 GHz maximum). The previous render command sets no nice value; this render ran at nice value -15.

The `sketches` and `pairs_scalar` stages do the same work in both revisions.
Their median time ratio, this render over the previous one, is
0.95 at one thread and
0.96 at four: that is the change of
load and nice value. `vs_same_work` divides a stage’s median time ratio by that
value. Peak RSS does not depend on host load or on the kind of core. Each ratio compares the
medians of three fresh processes in one cell; the summary gives the median,
minimum and maximum over the cells of a stage:

| stage        | threads | code      | median_time_ratio | min_time_ratio | max_time_ratio | vs_same_work | median_rss_ratio |
|:-------------|--------:|:----------|------------------:|---------------:|---------------:|-------------:|-----------------:|
| sketches     |       1 | same work |             0.948 |          0.936 |          0.961 |        0.997 |            1.004 |
| sketches     |       4 | same work |             0.955 |          0.916 |          0.972 |        0.995 |            1.002 |
| pairs_reduce |       1 | changed   |             0.647 |          0.441 |          1.218 |        0.680 |            0.564 |
| pairs_reduce |       4 | changed   |             0.349 |          0.251 |          1.323 |        0.364 |            0.622 |
| pairs_copy   |       1 | changed   |             0.651 |          0.446 |          1.168 |        0.685 |            0.730 |
| pairs_copy   |       4 | changed   |             0.669 |          0.445 |          1.306 |        0.696 |            0.715 |
| pairs_scalar |       1 | same work |             0.952 |          0.946 |          0.965 |        1.002 |            1.007 |
| pairs_scalar |       4 | same work |             0.961 |          0.927 |          0.965 |        1.000 |            1.002 |

`pairs_reduce`, cell by cell:

| samples |   sites | threads | seconds_before | rss_mib_before | seconds_after | rss_mib_after | time_ratio | rss_ratio |
|--------:|--------:|--------:|---------------:|---------------:|--------------:|--------------:|-----------:|----------:|
|     125 |   4,250 |       1 |          0.023 |         59.750 |         0.017 |        47.547 |      0.742 |     0.796 |
|     250 |   8,500 |       1 |          0.116 |         77.863 |         0.076 |        50.570 |      0.660 |     0.649 |
|     125 |  17,000 |       1 |          0.072 |         86.930 |         0.038 |        50.535 |      0.536 |     0.581 |
|     250 |  17,000 |       1 |          0.208 |        100.215 |         0.131 |        54.715 |      0.633 |     0.546 |
|     500 |  17,000 |       1 |          0.733 |        113.953 |         0.486 |        61.945 |      0.664 |     0.544 |
|     250 |  34,000 |       1 |          0.385 |        160.840 |         0.238 |        60.633 |      0.619 |     0.377 |
|     250 |  68,000 |       1 |          1.034 |        205.785 |         0.456 |        77.812 |      0.441 |     0.378 |
|       3 | 108,757 |       1 |          0.005 |         44.879 |         0.006 |        46.062 |      1.218 |     1.026 |
|     125 |   4,250 |       4 |          0.023 |         58.453 |         0.012 |        54.059 |      0.502 |     0.925 |
|     250 |   8,500 |       4 |          0.113 |         77.410 |         0.043 |        60.867 |      0.383 |     0.786 |
|     125 |  17,000 |       4 |          0.073 |         95.812 |         0.024 |        58.371 |      0.336 |     0.609 |
|     250 |  17,000 |       4 |          0.213 |        113.250 |         0.075 |        68.789 |      0.354 |     0.607 |
|     500 |  17,000 |       4 |          0.785 |        135.070 |         0.256 |        85.691 |      0.326 |     0.634 |
|     250 |  34,000 |       4 |          0.396 |        159.672 |         0.137 |        86.453 |      0.345 |     0.541 |
|     250 |  68,000 |       4 |          1.045 |        211.387 |         0.262 |       117.996 |      0.251 |     0.558 |
|       3 | 108,757 |       4 |          0.004 |         45.367 |         0.006 |        47.863 |      1.323 |     1.055 |

On the largest shared cell, 500 samples of
17,000 sites, `pairs_reduce` went from
114 to 61.9 MiB
peak RSS on one thread and from 135.1 to
85.7 MiB on four. Its time went from
0.73 to 0.49 s on one
thread and from 0.79 to
0.26 s on four.

## Time exponents and memory ratios

Time exponents per doubling, against the work of each stage, and peak-RSS ratios to 1×:

| stage        | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | timing_verdict     | exponent_above_1.25 |
|:-------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|:-------------------|:--------------------|
| sketches     | samples   |       1 |            0.941 |            0.960 |        1.092 |        1.314 | none: 1x under 5 s | FALSE               |
| sketches     | samples   |       4 |            0.834 |            0.909 |        1.082 |        1.248 | none: 1x under 5 s | FALSE               |
| sketches     | sites     |       1 |            1.019 |            0.998 |        1.287 |        1.886 | none: 1x under 5 s | FALSE               |
| sketches     | sites     |       4 |            0.941 |            0.972 |        1.307 |        1.835 | none: 1x under 5 s | FALSE               |
| sketches     | joint     |       1 |            0.914 |            0.976 |        1.104 |        1.557 | none: 1x under 5 s | FALSE               |
| sketches     | joint     |       4 |            0.824 |            0.900 |        1.146 |        1.563 | none: 1x under 5 s | FALSE               |
| pairs_reduce | samples   |       1 |            0.884 |            0.943 |        1.083 |        1.226 | none: 1x under 5 s | FALSE               |
| pairs_reduce | samples   |       4 |            0.810 |            0.882 |        1.178 |        1.467 | none: 1x under 5 s | FALSE               |
| pairs_reduce | sites     |       1 |            0.859 |            0.937 |        1.108 |        1.422 | none: 1x under 5 s | FALSE               |
| pairs_reduce | sites     |       4 |            0.857 |            0.943 |        1.257 |        1.715 | none: 1x under 5 s | FALSE               |
| pairs_reduce | joint     |       1 |            0.721 |            0.890 |        1.065 |        1.303 | none: 1x under 5 s | FALSE               |
| pairs_reduce | joint     |       4 |            0.627 |            0.854 |        1.126 |        1.584 | none: 1x under 5 s | FALSE               |
| pairs_copy   | samples   |       1 |            0.897 |            0.957 |        1.394 |        2.825 | none: 1x under 5 s | FALSE               |
| pairs_copy   | samples   |       4 |            0.876 |            0.924 |        1.372 |        2.622 | none: 1x under 5 s | FALSE               |
| pairs_copy   | sites     |       1 |            0.800 |            0.894 |        1.104 |        1.280 | none: 1x under 5 s | FALSE               |
| pairs_copy   | sites     |       4 |            0.754 |            0.888 |        1.175 |        1.376 | none: 1x under 5 s | FALSE               |
| pairs_copy   | joint     |       1 |            0.724 |            0.858 |        1.396 |        2.982 | none: 1x under 5 s | FALSE               |
| pairs_copy   | joint     |       4 |            0.718 |            0.842 |        1.364 |        2.804 | none: 1x under 5 s | FALSE               |
| pairs_scalar | samples   |       1 |            0.915 |            0.971 |        1.164 |        1.311 | none: 1x under 5 s | FALSE               |
| pairs_scalar | samples   |       4 |            0.910 |            0.970 |        1.063 |        1.194 | none: 1x under 5 s | FALSE               |
| pairs_scalar | sites     |       1 |            0.953 |            1.157 |        1.639 |        2.024 | none: 1x under 5 s | FALSE               |
| pairs_scalar | sites     |       4 |            0.976 |            1.247 |        1.468 |        2.041 | none: 1x under 5 s | FALSE               |
| pairs_scalar | joint     |       1 |            0.891 |            0.969 |        1.332 |        2.051 | none: 1x under 5 s | FALSE               |
| pairs_scalar | joint     |       4 |            0.884 |            0.964 |        1.361 |        2.045 | none: 1x under 5 s | FALSE               |

## Findings

Every observation is within its declared budget and none spills; the render
enforces both. The largest observation uses
88% of its ceiling.

The largest all-pairs cell (500 samples,
17,000 sites, 124,750 pairs) takes
0.55 s on one thread and
0.55 s on four when written to Parquet, at
150 and
153.1 MiB peak RSS. The Parquet output is
1266.1 KiB.

On the real counts, one thread builds three sketches of 108,757 sites in
0.162 s at
129.7 MiB peak RSS, and writes the three pairs in
0.006 s at
45.9 MiB.

For 124,750 pairs of 17,000
sites on one thread, the per-pair function takes
2.15 s and the all-pairs macro
0.49 s, 4.4
times faster on the same cores. At revision `a5f3b47a`, on CPUs 0-3, the
per-pair function took 4.4 s. The pair stages still
grow with pairs times sites, which is the size of the comparison itself.

In `pairs_reduce`, the macro compares pairs on DuckDB’s threads: the largest cell takes
0.49 s on one thread and
0.26 s on four,
1.9
times faster. In `pairs_copy`, written to Parquet, it does not: `COPY` keeps
insertion order by default, and the relations that feed the pairs give DuckDB
no batch index, so the copy runs as one task. The macro guarantees no row order. DuckDB’s Parquet
writer runs in parallel when insertion order is not kept
(`SET preserve_insertion_order = false`); this report does not measure that
setting.

24 of the 24 series have a one-thread 1× run under 5 seconds and get a memory verdict only.
0 series have a time exponent above 1.25 at either doubling.

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
- The macro compares pairs on up to 32 threads. This report measures one and
  four threads, on the cores named above. `COPY` to Parquet with the default
  insertion order compares the pairs in one task.
- `max_pairs` caps the output at 50,000,000 pairs by default, 400 times the
  largest cell here. Cells between 124,750 pairs and the cap are not
  measured.
- Count extraction, sketch verification, CHARR, matched contamination, the R
  wrappers and the reader `*_convert_parquet_sql` COPY statements are not
  measured here.
- The DuckDB CLI is v1.5.1 (Variegata) 7dbb2e646f; the R package and other hosts can
  differ.
