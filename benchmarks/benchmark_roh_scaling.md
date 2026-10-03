ROH macro scaling
================

This report gives the 1×/2×/4× scaling evidence `STYLE.md` asks for, for the two
runs-of-homozygosity entry points added for \#318:

- `duckhts_roh_counts` over a relation of read counts. Its growth dimensions are
  the sites per sample and chromosome (the length of one kernel list) and the
  number of samples. Synthetic counts give the 1×/2×/4× series. A real public
  workload adds read counts from three 1000 Genomes 30× CRAMs at the 108,757
  chr20 sites of the ROH ancestry evaluation, decoded for 1, 2 and 3 samples.
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

The synthetic read counts are deterministic. They are generated at render
time into a temporary directory, which is untimed staging. Each sample has a
2,000-site homozygous stretch every 10,000 sites over Hardy-Weinberg
background, with genotypes drawn at each site’s frequency and 30 reads per site,
so every input yields runs. The real counts are registry artifacts
`roh_counts_chr20_na18507`, `roh_counts_chr20_hg00403` and
`roh_counts_chr20_hg00188`: `duckhts_somalier_bam_counts` at the chr20
evaluation sites of the registered public CRAMs, read by range requests, with
MAPQ ≥ 20, base quality ≥ 13 and mate-overlap suppression
(`benchmarks/roh_counts_stage.R`). Each is checked against its registered row
count and counts digest before anything is measured. The ancestry
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

Source revision: ce79e1f9651156b6405342fe866ea97980d299ad; `src` tree 47d02e582ca0974166a5dba434945157899d8383. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: 45aea60192672fcdc0831e81aa4beeb7da6da9eec5b7967f33b6b48afadf75d9. Host load before the run: 1.15, 1.20, 0.87. Empty-run overhead: 38.2 MiB. Every observation is in [`benchmark_roh_scaling_runs.csv`](benchmark_roh_scaling_runs.csv).

| macro                | dimension    | scale | threads | samples | input_rows | runs_out | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib | within_budget |
|:---------------------|:-------------|------:|--------:|--------:|-----------:|:---------|---------------:|------------:|------------:|---------------:|------------:|-----------:|:--------------|
| duckhts_roh_ancestry | joint        |     1 |       1 |       8 |    6793144 | 356      |          6.247 |       6.200 |       6.249 |         1780.5 |      1821.3 |     1233.4 | FALSE         |
| duckhts_roh_ancestry | joint        |     2 |       1 |      16 |   13586288 | 1256     |         10.128 |      10.086 |      10.152 |         3440.4 |      3449.8 |     2428.5 | FALSE         |
| duckhts_roh_ancestry | joint        |     4 |       1 |      32 |   27172576 | 5120     |         18.219 |      18.214 |      18.279 |         6601.3 |      6601.8 |     4818.8 | FALSE         |
| duckhts_roh_ancestry | joint        |     1 |       4 |       8 |    6793144 | 356      |          5.072 |       5.067 |       5.122 |         1940.7 |      1943.0 |     1233.4 | FALSE         |
| duckhts_roh_ancestry | joint        |     2 |       4 |      16 |   13586288 | 1256     |          7.630 |       7.597 |       7.653 |         3490.6 |      3553.7 |     2428.5 | FALSE         |
| duckhts_roh_ancestry | joint        |     4 |       4 |      32 |   27172576 | 5120     |         12.621 |      12.346 |      12.797 |         6760.4 |      6887.1 |     4818.8 | FALSE         |
| duckhts_roh_ancestry | reference    |     1 |       1 |       8 |    6793144 | 356      |          6.247 |       6.200 |       6.249 |         1780.5 |      1821.3 |     1233.4 | FALSE         |
| duckhts_roh_ancestry | reference    |     2 |       1 |       8 |    6793144 | 641      |          6.551 |       6.529 |       6.558 |         2156.9 |      2158.8 |     1340.1 | FALSE         |
| duckhts_roh_ancestry | reference    |     4 |       1 |       8 |    6793144 | 1337     |          7.254 |       7.166 |       7.320 |         2717.0 |      2729.2 |     1553.6 | FALSE         |
| duckhts_roh_ancestry | reference    |     1 |       4 |       8 |    6793144 | 356      |          5.072 |       5.067 |       5.122 |         1940.7 |      1943.0 |     1233.4 | FALSE         |
| duckhts_roh_ancestry | reference    |     2 |       4 |       8 |    6793144 | 641      |          5.373 |       5.318 |       5.420 |         1929.1 |      2028.7 |     1340.1 | FALSE         |
| duckhts_roh_ancestry | reference    |     4 |       4 |       8 |    6793144 | 1337     |          5.611 |       5.550 |       5.691 |         2063.0 |      2108.3 |     1553.6 | FALSE         |
| duckhts_roh_ancestry | samples      |     1 |       1 |       8 |    6793144 | 1337     |          7.254 |       7.166 |       7.320 |         2717.0 |      2729.2 |     1553.6 | FALSE         |
| duckhts_roh_ancestry | samples      |     2 |       1 |      16 |   13586288 | 2637     |         10.858 |      10.803 |      10.993 |         4040.8 |      4043.2 |     2642.0 | FALSE         |
| duckhts_roh_ancestry | samples      |     4 |       1 |      32 |   27172576 | 5120     |         18.219 |      18.214 |      18.279 |         6601.3 |      6601.8 |     4818.8 | FALSE         |
| duckhts_roh_ancestry | samples      |     1 |       4 |       8 |    6793144 | 1337     |          5.611 |       5.550 |       5.691 |         2063.0 |      2108.3 |     1553.6 | FALSE         |
| duckhts_roh_ancestry | samples      |     2 |       4 |      16 |   13586288 | 2637     |          7.909 |       7.721 |       7.931 |         3715.2 |      3935.0 |     2642.0 | FALSE         |
| duckhts_roh_ancestry | samples      |     4 |       4 |      32 |   27172576 | 5120     |         12.621 |      12.346 |      12.797 |         6760.4 |      6887.1 |     4818.8 | FALSE         |
| duckhts_roh_counts   | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.022 |       0.021 |       0.023 |           54.4 |        54.5 |       65.6 | TRUE          |
| duckhts_roh_counts   | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.040 |       0.039 |       0.041 |           65.3 |        65.7 |       93.0 | TRUE          |
| duckhts_roh_counts   | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.060 |       0.059 |       0.061 |           79.3 |        79.6 |      120.4 | TRUE          |
| duckhts_roh_counts   | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.024 |       0.023 |       0.028 |           60.7 |        61.6 |       65.6 | TRUE          |
| duckhts_roh_counts   | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.030 |       0.029 |       0.033 |           77.3 |        78.1 |       93.0 | TRUE          |
| duckhts_roh_counts   | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.044 |       0.044 |       0.045 |           87.7 |        92.8 |      120.4 | TRUE          |
| duckhts_roh_counts   | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.329 |       0.327 |       0.334 |          254.9 |       254.9 |      541.8 | TRUE          |
| duckhts_roh_counts   | samples      |     2 |       1 |      16 |    4000000 | 400      |          0.661 |       0.659 |       0.669 |          466.6 |       466.7 |     1045.3 | TRUE          |
| duckhts_roh_counts   | samples      |     4 |       1 |      32 |    8000000 | 800      |          1.337 |       1.331 |       1.339 |          903.5 |       912.8 |     2052.4 | TRUE          |
| duckhts_roh_counts   | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.193 |       0.171 |       0.196 |          260.4 |       272.2 |      541.8 | TRUE          |
| duckhts_roh_counts   | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.329 |       0.328 |       0.346 |          463.1 |       489.4 |     1045.3 | TRUE          |
| duckhts_roh_counts   | samples      |     4 |       4 |      32 |    8000000 | 800      |          0.592 |       0.580 |       0.669 |          915.4 |       917.7 |     2052.4 | TRUE          |
| duckhts_roh_counts   | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.177 |       0.167 |       0.182 |          162.4 |       167.8 |      290.0 | TRUE          |
| duckhts_roh_counts   | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.338 |       0.331 |       0.338 |          310.4 |       310.7 |      541.8 | TRUE          |
| duckhts_roh_counts   | sites        |     4 |       1 |       1 |    4000000 | 400      |          0.702 |       0.699 |       0.704 |          702.2 |       702.4 |     1045.3 | TRUE          |
| duckhts_roh_counts   | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.176 |       0.174 |       0.177 |          175.8 |       183.7 |      290.0 | TRUE          |
| duckhts_roh_counts   | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.358 |       0.342 |       0.363 |          346.9 |       351.8 |      541.8 | TRUE          |
| duckhts_roh_counts   | sites        |     4 |       4 |       1 |    4000000 | 400      |          0.685 |       0.681 |       0.689 |          702.8 |       716.7 |     1045.3 | TRUE          |

Time exponents per doubling (log2 of the median ratio) and peak-RSS ratios to 1×:

| macro                | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | one_x_seconds |
|:---------------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|--------------:|
| duckhts_roh_ancestry | joint     |       1 |            0.697 |            0.847 |        1.932 |        3.708 |         6.247 |
| duckhts_roh_ancestry | reference |       1 |            0.069 |            0.147 |        1.211 |        1.526 |         6.247 |
| duckhts_roh_ancestry | samples   |       1 |            0.582 |            0.747 |        1.487 |        2.430 |         7.254 |
| duckhts_roh_counts   | samples   |       1 |            1.007 |            1.017 |        1.831 |        3.545 |         0.329 |
| duckhts_roh_counts   | sites     |       1 |            0.929 |            1.056 |        1.911 |        4.324 |         0.177 |
| duckhts_roh_ancestry | joint     |       4 |            0.589 |            0.726 |        1.799 |        3.483 |         5.072 |
| duckhts_roh_ancestry | reference |       4 |            0.083 |            0.063 |        0.994 |        1.063 |         5.072 |
| duckhts_roh_ancestry | samples   |       4 |            0.495 |            0.674 |        1.801 |        3.277 |         5.611 |
| duckhts_roh_counts   | samples   |       4 |            0.772 |            0.846 |        1.778 |        3.515 |         0.193 |
| duckhts_roh_counts   | sites     |       4 |            1.020 |            0.936 |        1.973 |        3.998 |         0.176 |

## Findings

`duckhts_roh_counts` stays within its declared budget at every size and
thread count, synthetic and real, and the render enforces that. Its one-thread
1× runs are under the 5-second floor, so it gets a memory verdict, not a timing
verdict. On the real workload, one thread decodes one sample’s 108,757 chr20
sites in 0.022 s at 54.4 MiB
peak RSS and three samples in 0.06 s at
79.3 MiB. The operating scale this report supports
is what it measures: per-chromosome lists of up to 4 million synthetic sites,
up to 32 synthetic samples, and three real 30× samples at 108,757 sites;
larger panels and cohorts are not measured here.

`duckhts_roh_ancestry` exceeds its declared model. Memory grows with the number of samples decoded together, because the per-(sample, chromosome) list aggregate keeps every sample’s site list until the BCF scan ends, and DuckDB’s list-aggregate state costs several times the 28 bytes of struct data per element. The reference stage adds a fixed cost that grows with reference rows (see the reference series). `duckhts_roh` and `duckhts_roh_af_table` share the list stage and are already on `main` (#323): the ancestry evaluation’s eight-child, 21-autosome queries peaked at 16.5 GB.

The maintainer approved this as a reviewed exception on 2026-10-03: the budget is reported, not gated, until \#329 makes the decode memory-bounded. Until then, `samples :=` is the memory control. Decode in sample batches, as the evaluation does.
