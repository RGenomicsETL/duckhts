ROH macro scaling
================

This report gives the 1×/2×/4× scaling evidence `STYLE.md` asks for, for the two
runs-of-homozygosity entry points added for \#318:

- `duckhts_roh_counts` over a relation of read counts. Its growth dimensions are
  the sites per sample and chromosome (the length of one kernel list) and the
  number of samples.
- `duckhts_roh_ancestry` over the registered chr20 children BCF of the ROH
  ancestry evaluation (849,143 records). It is measured along three series:
  - the number of samples decoded together (8, 16, 32) with the full reference;
  - the reference rows, with 8 samples and the chr20 reference cut at position
    quantiles to a quarter, a half and all of its sites (the BCF is unchanged,
    so only the reference stage grows);
  - both together: (8, quarter), (16, half) and (32, all).

  A configuration shared by several series is measured once and listed in each.

Every observation is one fresh DuckDB CLI process under `/usr/bin/time`. The
JSON profiler gives query latency and DuckDB’s peak buffer and temporary
bytes, and `/usr/bin/time` gives peak RSS. `max_temp_directory_size` is zero,
so a spill fails the render. Each cell has three repetitions at one thread and
at four threads, run in an order that alternates sizes. Each query reduces the
runs it returns to a count and a total length, so the macro’s output is
consumed.

The read counts are synthetic and deterministic. They are generated at render
time into a temporary directory, which is untimed staging. Each sample has a
2,000-site homozygous stretch every 10,000 sites over Hardy-Weinberg
background, with genotypes drawn at each site’s frequency and 30 reads per site,
so every input yields runs. The ancestry
inputs are registry artifacts `roh_ancestry_chr20_children_bcf`,
`roh_ancestry_chr20_reference_long` and `roh_ancestry_chr20_q`. Their samples
are the first 8, 16 or 32 sample IDs of the registered proportions, sorted.
The reference cuts are the 25% and 50% position quantiles of its distinct sites.

## Budgets, declared before measuring

The overhead is the median peak RSS of three empty runs that start the CLI,
load the extension and run no scan.

- **`duckhts_roh_counts`.** The live state is every sample’s site list, held by
  the per-(sample, chromosome) list aggregate until the scan ends. That is 24
  bytes of struct data per site (position, frequency and two counts), taken as
  48 bytes with list-segment and validity overhead. On top of that come the
  kernel arrays for the longest list (emission, path and position storage, 40
  bytes per site). The ceiling is overhead + 3 × rows × 88 bytes. This budget
  is a gate: the render stops if any counts observation exceeds it or spills.
- **`duckhts_roh_ancestry`.** The same list aggregate holds samples × records
  elements. Each element is 28 bytes of struct data (position, frequency and
  genotype evidence), taken as 56 bytes. The reference join adds one hash
  entry per reference row (64 bytes). The ceiling is overhead + 3 × (samples ×
  records × 56 bytes + reference rows × 64 bytes). The report states whether
  each cell meets it. See *Findings* for why this budget is reported rather
  than gated.

``` sh
taskset -c 8-15 Rscript -e 'rmarkdown::render("benchmarks/benchmark_roh_scaling.Rmd")'
```

Source revision: 3818c6a9dd1260632da1b8ef79af7443089606c3; `src` tree 47d02e582ca0974166a5dba434945157899d8383. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: 45aea60192672fcdc0831e81aa4beeb7da6da9eec5b7967f33b6b48afadf75d9. Host load before the run: 0.10, 0.25, 0.42. Empty-run overhead: 38.3 MiB. Every observation is in [`benchmark_roh_scaling_runs.csv`](benchmark_roh_scaling_runs.csv).

| macro                | dimension | scale | threads | samples | input_rows | runs_out | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib | within_budget |
|:---------------------|:----------|------:|--------:|--------:|-----------:|:---------|---------------:|------------:|------------:|---------------:|------------:|-----------:|:--------------|
| duckhts_roh_ancestry | joint     |     1 |       1 |       8 |    6793144 | 356      |          6.222 |       6.221 |       6.370 |         1780.3 |      1786.5 |     1233.4 | FALSE         |
| duckhts_roh_ancestry | joint     |     2 |       1 |      16 |   13586288 | 1256     |         10.084 |      10.035 |      10.392 |         3437.9 |      3439.5 |     2428.5 | FALSE         |
| duckhts_roh_ancestry | joint     |     4 |       1 |      32 |   27172576 | 5120     |         18.249 |      18.212 |      18.771 |         6599.0 |      6603.6 |     4818.8 | FALSE         |
| duckhts_roh_ancestry | joint     |     1 |       4 |       8 |    6793144 | 356      |          5.023 |       5.021 |       5.203 |         1808.6 |      1836.2 |     1233.4 | FALSE         |
| duckhts_roh_ancestry | joint     |     2 |       4 |      16 |   13586288 | 1256     |          7.591 |       7.582 |       7.647 |         3474.4 |      3522.0 |     2428.5 | FALSE         |
| duckhts_roh_ancestry | joint     |     4 |       4 |      32 |   27172576 | 5120     |         12.507 |      12.332 |      12.951 |         6388.3 |      7007.9 |     4818.8 | FALSE         |
| duckhts_roh_ancestry | reference |     1 |       1 |       8 |    6793144 | 356      |          6.222 |       6.221 |       6.370 |         1780.3 |      1786.5 |     1233.4 | FALSE         |
| duckhts_roh_ancestry | reference |     2 |       1 |       8 |    6793144 | 641      |          6.564 |       6.555 |       6.790 |         2158.2 |      2166.7 |     1340.2 | FALSE         |
| duckhts_roh_ancestry | reference |     4 |       1 |       8 |    6793144 | 1337     |          7.316 |       7.267 |       7.571 |         2719.3 |      2723.4 |     1553.7 | FALSE         |
| duckhts_roh_ancestry | reference |     1 |       4 |       8 |    6793144 | 356      |          5.023 |       5.021 |       5.203 |         1808.6 |      1836.2 |     1233.4 | FALSE         |
| duckhts_roh_ancestry | reference |     2 |       4 |       8 |    6793144 | 641      |          5.367 |       5.363 |       5.400 |         1950.8 |      1996.7 |     1340.2 | FALSE         |
| duckhts_roh_ancestry | reference |     4 |       4 |       8 |    6793144 | 1337     |          5.610 |       5.574 |       5.678 |         2111.9 |      2124.4 |     1553.7 | FALSE         |
| duckhts_roh_ancestry | samples   |     1 |       1 |       8 |    6793144 | 1337     |          7.316 |       7.267 |       7.571 |         2719.3 |      2723.4 |     1553.7 | FALSE         |
| duckhts_roh_ancestry | samples   |     2 |       1 |      16 |   13586288 | 2637     |         10.837 |      10.815 |      11.007 |         4040.3 |      4062.1 |     2642.1 | FALSE         |
| duckhts_roh_ancestry | samples   |     4 |       1 |      32 |   27172576 | 5120     |         18.249 |      18.212 |      18.771 |         6599.0 |      6603.6 |     4818.8 | FALSE         |
| duckhts_roh_ancestry | samples   |     1 |       4 |       8 |    6793144 | 1337     |          5.610 |       5.574 |       5.678 |         2111.9 |      2124.4 |     1553.7 | FALSE         |
| duckhts_roh_ancestry | samples   |     2 |       4 |      16 |   13586288 | 2637     |          7.868 |       7.749 |       7.948 |         3709.4 |      3767.8 |     2642.1 | FALSE         |
| duckhts_roh_ancestry | samples   |     4 |       4 |      32 |   27172576 | 5120     |         12.507 |      12.332 |      12.951 |         6388.3 |      7007.9 |     4818.8 | FALSE         |
| duckhts_roh_counts   | samples   |     1 |       1 |       8 |    2000000 | 200      |          0.333 |       0.326 |       0.337 |          254.8 |       254.8 |      541.8 | TRUE          |
| duckhts_roh_counts   | samples   |     2 |       1 |      16 |    4000000 | 400      |          0.664 |       0.662 |       0.679 |          466.4 |       466.6 |     1045.3 | TRUE          |
| duckhts_roh_counts   | samples   |     4 |       1 |      32 |    8000000 | 800      |          1.350 |       1.345 |       1.401 |          912.4 |       912.8 |     2052.4 | TRUE          |
| duckhts_roh_counts   | samples   |     1 |       4 |       8 |    2000000 | 200      |          0.178 |       0.167 |       0.209 |          268.8 |       277.3 |      541.8 | TRUE          |
| duckhts_roh_counts   | samples   |     2 |       4 |      16 |    4000000 | 400      |          0.331 |       0.328 |       0.352 |          456.8 |       466.4 |     1045.3 | TRUE          |
| duckhts_roh_counts   | samples   |     4 |       4 |      32 |    8000000 | 800      |          0.630 |       0.583 |       0.634 |          912.8 |       919.2 |     2052.4 | TRUE          |
| duckhts_roh_counts   | sites     |     1 |       1 |       1 |    1000000 | 100      |          0.166 |       0.166 |       0.188 |          167.9 |       168.0 |      290.0 | TRUE          |
| duckhts_roh_counts   | sites     |     2 |       1 |       1 |    2000000 | 200      |          0.337 |       0.333 |       0.352 |          310.5 |       310.5 |      541.8 | TRUE          |
| duckhts_roh_counts   | sites     |     4 |       1 |       1 |    4000000 | 400      |          0.701 |       0.700 |       0.709 |          702.3 |       702.3 |     1045.3 | TRUE          |
| duckhts_roh_counts   | sites     |     1 |       4 |       1 |    1000000 | 100      |          0.180 |       0.175 |       0.183 |          181.6 |       181.6 |      290.0 | TRUE          |
| duckhts_roh_counts   | sites     |     2 |       4 |       1 |    2000000 | 200      |          0.344 |       0.339 |       0.352 |          355.8 |       356.5 |      541.8 | TRUE          |
| duckhts_roh_counts   | sites     |     4 |       4 |       1 |    4000000 | 400      |          0.691 |       0.691 |       0.697 |          697.4 |       701.4 |     1045.3 | TRUE          |

Time exponents per doubling (log2 of the median ratio) and peak-RSS ratios to 1×:

| macro                | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | one_x_seconds |
|:---------------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|--------------:|
| duckhts_roh_ancestry | joint     |       1 |            0.697 |            0.856 |        1.931 |        3.707 |         6.222 |
| duckhts_roh_ancestry | reference |       1 |            0.077 |            0.157 |        1.212 |        1.527 |         6.222 |
| duckhts_roh_ancestry | samples   |       1 |            0.567 |            0.752 |        1.486 |        2.427 |         7.316 |
| duckhts_roh_counts   | samples   |       1 |            0.997 |            1.025 |        1.830 |        3.581 |         0.333 |
| duckhts_roh_counts   | sites     |       1 |            1.022 |            1.056 |        1.849 |        4.183 |         0.166 |
| duckhts_roh_ancestry | joint     |       4 |            0.596 |            0.720 |        1.921 |        3.532 |         5.023 |
| duckhts_roh_ancestry | reference |       4 |            0.095 |            0.064 |        1.079 |        1.168 |         5.023 |
| duckhts_roh_ancestry | samples   |       4 |            0.488 |            0.669 |        1.756 |        3.025 |         5.610 |
| duckhts_roh_counts   | samples   |       4 |            0.894 |            0.931 |        1.699 |        3.396 |         0.178 |
| duckhts_roh_counts   | sites     |       4 |            0.932 |            1.007 |        1.959 |        3.840 |         0.180 |

## Findings

`duckhts_roh_counts` stays within its declared budget at every size and
thread count, and the render enforces that. Its one-thread 1× runs are under
the 5-second floor, so it gets a memory verdict, not a timing verdict.

`duckhts_roh_ancestry` exceeds its declared model. Memory grows with the number of samples decoded together, because the per-(sample, chromosome) list aggregate keeps every sample’s site list until the BCF scan ends, and DuckDB’s list-aggregate state costs several times the 28 bytes of struct data per element. The reference stage adds a fixed cost that grows with reference rows (see the reference series). `duckhts_roh` and `duckhts_roh_af_table` share the list stage and are already on `main` (#323): the ancestry evaluation’s eight-child, 21-autosome queries peaked at 16.5 GB.

The maintainer approved this as a reviewed exception on 2026-10-03: the budget is reported, not gated, until \#329 makes the decode memory-bounded. Until then, `samples :=` is the memory control. Decode in sample batches, as the evaluation does.
