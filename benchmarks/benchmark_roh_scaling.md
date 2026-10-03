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

Source revision: cb3b36cd8e95322c695b347901adfedfd76731d6; `src` tree fbed0118c09af45ac5de804e5da4b9b99cde5c4d. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: beadbb8aefd3adcb7da3ea15eee3b6836830e06d3b70b146145ee651fbccca17. Host load before the run: 1.07, 0.76, 0.68. Empty-run overhead: 38.2 MiB. Every observation is in [`benchmark_roh_scaling_runs.csv`](benchmark_roh_scaling_runs.csv).

| macro                | dimension | scale | threads | samples | input_rows | runs_out | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib | within_budget |
|:---------------------|:----------|------:|--------:|--------:|-----------:|:---------|---------------:|------------:|------------:|---------------:|------------:|-----------:|:--------------|
| duckhts_roh_ancestry | joint     |     1 |       1 |       8 |    6793144 | 356      |          6.194 |       6.182 |       6.227 |         1781.1 |      1782.3 |     1233.3 | FALSE         |
| duckhts_roh_ancestry | joint     |     2 |       1 |      16 |   13586288 | 1256     |         10.051 |      10.035 |      10.084 |         3448.1 |      3449.0 |     2428.5 | FALSE         |
| duckhts_roh_ancestry | joint     |     4 |       1 |      32 |   27172576 | 5120     |         18.120 |      18.063 |      18.255 |         6595.2 |      6595.4 |     4818.7 | FALSE         |
| duckhts_roh_ancestry | joint     |     1 |       4 |       8 |    6793144 | 356      |          4.960 |       4.953 |       5.295 |         1822.7 |      1826.3 |     1233.3 | FALSE         |
| duckhts_roh_ancestry | joint     |     2 |       4 |      16 |   13586288 | 1256     |          7.534 |       7.532 |       8.260 |         3553.8 |      3674.6 |     2428.5 | FALSE         |
| duckhts_roh_ancestry | joint     |     4 |       4 |      32 |   27172576 | 5120     |         12.527 |      12.383 |      13.244 |         6659.4 |      6809.7 |     4818.7 | FALSE         |
| duckhts_roh_ancestry | reference |     1 |       1 |       8 |    6793144 | 356      |          6.194 |       6.182 |       6.227 |         1781.1 |      1782.3 |     1233.3 | FALSE         |
| duckhts_roh_ancestry | reference |     2 |       1 |       8 |    6793144 | 641      |          6.522 |       6.502 |       6.523 |         2158.4 |      2165.4 |     1340.1 | FALSE         |
| duckhts_roh_ancestry | reference |     4 |       1 |       8 |    6793144 | 1337     |          7.219 |       7.216 |       7.241 |         2716.7 |      2740.7 |     1553.6 | FALSE         |
| duckhts_roh_ancestry | reference |     1 |       4 |       8 |    6793144 | 356      |          4.960 |       4.953 |       5.295 |         1822.7 |      1826.3 |     1233.3 | FALSE         |
| duckhts_roh_ancestry | reference |     2 |       4 |       8 |    6793144 | 641      |          5.390 |       5.192 |       5.675 |         2012.9 |      2023.0 |     1340.1 | FALSE         |
| duckhts_roh_ancestry | reference |     4 |       4 |       8 |    6793144 | 1337     |          5.581 |       5.507 |       5.589 |         2078.7 |      2104.6 |     1553.6 | FALSE         |
| duckhts_roh_ancestry | samples   |     1 |       1 |       8 |    6793144 | 1337     |          7.219 |       7.216 |       7.241 |         2716.7 |      2740.7 |     1553.6 | FALSE         |
| duckhts_roh_ancestry | samples   |     2 |       1 |      16 |   13586288 | 2637     |         10.812 |      10.798 |      10.822 |         4011.8 |      4021.2 |     2642.0 | FALSE         |
| duckhts_roh_ancestry | samples   |     4 |       1 |      32 |   27172576 | 5120     |         18.120 |      18.063 |      18.255 |         6595.2 |      6595.4 |     4818.7 | FALSE         |
| duckhts_roh_ancestry | samples   |     1 |       4 |       8 |    6793144 | 1337     |          5.581 |       5.507 |       5.589 |         2078.7 |      2104.6 |     1553.6 | FALSE         |
| duckhts_roh_ancestry | samples   |     2 |       4 |      16 |   13586288 | 2637     |          7.908 |       7.872 |       8.022 |         3620.2 |      3673.0 |     2642.0 | FALSE         |
| duckhts_roh_ancestry | samples   |     4 |       4 |      32 |   27172576 | 5120     |         12.527 |      12.383 |      13.244 |         6659.4 |      6809.7 |     4818.7 | FALSE         |
| duckhts_roh_counts   | samples   |     1 |       1 |       8 |    2000000 | 200      |          0.335 |       0.331 |       0.338 |          255.0 |       255.1 |      541.7 | TRUE          |
| duckhts_roh_counts   | samples   |     2 |       1 |      16 |    4000000 | 400      |          0.662 |       0.657 |       0.663 |          466.6 |       466.8 |     1045.3 | TRUE          |
| duckhts_roh_counts   | samples   |     4 |       1 |      32 |    8000000 | 800      |          1.352 |       1.331 |       1.374 |          912.9 |       913.0 |     2052.3 | TRUE          |
| duckhts_roh_counts   | samples   |     1 |       4 |       8 |    2000000 | 200      |          0.193 |       0.169 |       0.198 |          256.3 |       273.8 |      541.7 | TRUE          |
| duckhts_roh_counts   | samples   |     2 |       4 |      16 |    4000000 | 400      |          0.344 |       0.343 |       0.346 |          440.4 |       459.0 |     1045.3 | TRUE          |
| duckhts_roh_counts   | samples   |     4 |       4 |      32 |    8000000 | 800      |          0.674 |       0.655 |       0.682 |          907.7 |       939.7 |     2052.3 | TRUE          |
| duckhts_roh_counts   | sites     |     1 |       1 |       1 |    1000000 | 100      |          0.178 |       0.174 |       0.184 |          162.5 |       182.9 |      289.9 | TRUE          |
| duckhts_roh_counts   | sites     |     2 |       1 |       1 |    2000000 | 200      |          0.339 |       0.337 |       0.363 |          310.6 |       310.9 |      541.7 | TRUE          |
| duckhts_roh_counts   | sites     |     4 |       1 |       1 |    4000000 | 400      |          0.704 |       0.704 |       0.753 |          702.2 |       732.7 |     1045.3 | TRUE          |
| duckhts_roh_counts   | sites     |     1 |       4 |       1 |    1000000 | 100      |          0.179 |       0.179 |       0.184 |          195.5 |       197.8 |      289.9 | TRUE          |
| duckhts_roh_counts   | sites     |     2 |       4 |       1 |    2000000 | 200      |          0.344 |       0.342 |       0.346 |          358.8 |       368.1 |      541.7 | TRUE          |
| duckhts_roh_counts   | sites     |     4 |       4 |       1 |    4000000 | 400      |          0.687 |       0.685 |       0.687 |          692.3 |       702.7 |     1045.3 | TRUE          |

Time exponents per doubling (log2 of the median ratio) and peak-RSS ratios to 1×:

| macro                | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | one_x_seconds |
|:---------------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|--------------:|
| duckhts_roh_ancestry | joint     |       1 |            0.698 |            0.850 |        1.936 |        3.703 |         6.194 |
| duckhts_roh_ancestry | reference |       1 |            0.075 |            0.146 |        1.212 |        1.525 |         6.194 |
| duckhts_roh_ancestry | samples   |       1 |            0.583 |            0.745 |        1.477 |        2.428 |         7.219 |
| duckhts_roh_counts   | samples   |       1 |            0.983 |            1.030 |        1.830 |        3.580 |         0.335 |
| duckhts_roh_counts   | sites     |       1 |            0.928 |            1.054 |        1.911 |        4.321 |         0.178 |
| duckhts_roh_ancestry | joint     |       4 |            0.603 |            0.733 |        1.950 |        3.654 |         4.960 |
| duckhts_roh_ancestry | reference |       4 |            0.120 |            0.050 |        1.104 |        1.140 |         4.960 |
| duckhts_roh_ancestry | samples   |       4 |            0.503 |            0.664 |        1.742 |        3.204 |         5.581 |
| duckhts_roh_counts   | samples   |       4 |            0.835 |            0.971 |        1.718 |        3.542 |         0.193 |
| duckhts_roh_counts   | sites     |       4 |            0.941 |            0.996 |        1.835 |        3.541 |         0.179 |

## Findings

`duckhts_roh_counts` stays within its declared budget at every size and
thread count, and the render enforces that. Its one-thread 1× runs are under
the 5-second floor, so it gets a memory verdict, not a timing verdict.

`duckhts_roh_ancestry` exceeds its declared model. Memory grows with the number of samples decoded together, because the per-(sample, chromosome) list aggregate keeps every sample’s site list until the BCF scan ends, and DuckDB’s list-aggregate state costs several times the 28 bytes of struct data per element. The reference stage adds a fixed cost that grows with reference rows (see the reference series). `duckhts_roh` and `duckhts_roh_af_table` share the list stage and are already on `main` (#323): the ancestry evaluation’s eight-child, 21-autosome queries peaked at 16.5 GB.

The maintainer approved this as a reviewed exception on 2026-10-03: the budget is reported, not gated, until \#329 makes the decode memory-bounded. Until then, `samples :=` is the memory control. Decode in sample batches, as the evaluation does.
