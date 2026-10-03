ROH macro scaling
================

This report gives the 1×/2×/4× scaling evidence `STYLE.md` asks for, for the two
runs-of-homozygosity entry points added for \#318. Since \#329 the sites of each
sample and chromosome are held in native buffers that DuckHTS bounds, instead
of DuckDB lists; the last section compares this revision with the one before.

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

Each ceiling is overhead + 3 × the live state below, and each is a gate: the
render stops if any observation exceeds its ceiling or spills.

- **`duckhts_roh_counts`.** The live state is:
  - every sample’s packed site list, 24 bytes per row (position, frequency and
    two counts), held natively until the scan ends and then by DuckDB as one
    BLOB per sample and chromosome;
  - the kernel workspace of each thread for the longest list, 38 bytes per site
    (position, two emissions, the Viterbi path and the forward-backward values).

  That is rows × 24 + threads × sites per list × 38 bytes.
- **`duckhts_roh_ancestry`.** A record without a reference site is dropped
  before the sites are collected, so the lists hold reference sites only: the
  reference rows divided by its number of groups. The live state is:
  - the packed site lists, 16 bytes per sample and reference site (position,
    frequency and genotype evidence);
  - the distinct called sites and the hash table that joins the reference to
    them, 64 bytes per BCF record each;
  - the reference join and the per-site frequency lists, 64 bytes per reference
    row;
  - the kernel workspace, 38 bytes per reference site and thread.

  That is samples × reference sites × 16 + records × 128 + reference rows × 64 +
  threads × reference sites × 38 bytes.

``` sh
taskset -c 8-15 Rscript -e 'rmarkdown::render("benchmarks/benchmark_roh_scaling.Rmd")'
```

Source revision: 4cb72dc8967d0eeb4e93a659552ea6c56f22eb55; `src` tree 5908c1addd715de1121c9da6bf1480372b441e94. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: e49af8c9c196f8655f6f74e5b8604915d964a9b5f165ee30a39c2af7f284b01a. Host load before the run: 2.85, 3.11, 2.41. Empty-run overhead: 38.5 MiB. Every observation is in [`benchmark_roh_scaling_runs.csv`](benchmark_roh_scaling_runs.csv).

| macro                | dimension    | scale | threads | samples | input_rows | runs_out | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib | within_budget |
|:---------------------|:-------------|------:|--------:|--------:|-----------:|:---------|---------------:|------------:|------------:|---------------:|------------:|-----------:|:--------------|
| duckhts_roh_ancestry | joint        |     1 |       1 |       8 |    6793144 | 356      |          3.746 |       3.729 |       3.765 |          245.2 |       245.9 |      469.4 | TRUE          |
| duckhts_roh_ancestry | joint        |     2 |       1 |      16 |   13586288 | 1256     |          5.121 |       5.103 |       5.123 |          387.3 |       388.9 |      609.6 | TRUE          |
| duckhts_roh_ancestry | joint        |     4 |       1 |      32 |   27172576 | 5120     |          8.249 |       8.247 |       8.316 |          693.9 |       705.9 |      951.2 | TRUE          |
| duckhts_roh_ancestry | joint        |     1 |       4 |       8 |    6793144 | 356      |          3.679 |       3.674 |       3.730 |          275.7 |       276.6 |      478.4 | TRUE          |
| duckhts_roh_ancestry | joint        |     2 |       4 |      16 |   13586288 | 1256     |          5.124 |       5.067 |       5.173 |          473.9 |       485.3 |      627.8 | TRUE          |
| duckhts_roh_ancestry | joint        |     4 |       4 |      32 |   27172576 | 5120     |          8.237 |       8.183 |       8.299 |          791.2 |       799.7 |      987.4 | TRUE          |
| duckhts_roh_ancestry | reference    |     1 |       1 |       8 |    6793144 | 356      |          3.746 |       3.729 |       3.765 |          245.2 |       245.9 |      469.4 | TRUE          |
| duckhts_roh_ancestry | reference    |     2 |       1 |       8 |    6793144 | 641      |          3.996 |       3.972 |       4.005 |          376.6 |       377.3 |      589.3 | TRUE          |
| duckhts_roh_ancestry | reference    |     4 |       1 |       8 |    6793144 | 1337     |          4.530 |       4.509 |       4.540 |          632.9 |       634.5 |      829.2 | TRUE          |
| duckhts_roh_ancestry | reference    |     1 |       4 |       8 |    6793144 | 356      |          3.679 |       3.674 |       3.730 |          275.7 |       276.6 |      478.4 | TRUE          |
| duckhts_roh_ancestry | reference    |     2 |       4 |       8 |    6793144 | 641      |          3.964 |       3.962 |       3.968 |          479.0 |       479.9 |      607.4 | TRUE          |
| duckhts_roh_ancestry | reference    |     4 |       4 |       8 |    6793144 | 1337     |          4.528 |       4.475 |       4.528 |          767.0 |       774.6 |      865.4 | TRUE          |
| duckhts_roh_ancestry | samples      |     1 |       1 |       8 |    6793144 | 1337     |          4.530 |       4.509 |       4.540 |          632.9 |       634.5 |      829.2 | TRUE          |
| duckhts_roh_ancestry | samples      |     2 |       1 |      16 |   13586288 | 2637     |          5.799 |       5.776 |       5.826 |          662.2 |       662.9 |      869.9 | TRUE          |
| duckhts_roh_ancestry | samples      |     4 |       1 |      32 |   27172576 | 5120     |          8.249 |       8.247 |       8.316 |          693.9 |       705.9 |      951.2 | TRUE          |
| duckhts_roh_ancestry | samples      |     1 |       4 |       8 |    6793144 | 1337     |          4.528 |       4.475 |       4.528 |          767.0 |       774.6 |      865.4 | TRUE          |
| duckhts_roh_ancestry | samples      |     2 |       4 |      16 |   13586288 | 2637     |          5.751 |       5.700 |       5.791 |          751.9 |       766.1 |      906.1 | TRUE          |
| duckhts_roh_ancestry | samples      |     4 |       4 |      32 |   27172576 | 5120     |          8.237 |       8.183 |       8.299 |          791.2 |       799.7 |      987.4 | TRUE          |
| duckhts_roh_counts   | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.019 |       0.018 |       0.019 |           49.0 |        49.0 |       57.7 | TRUE          |
| duckhts_roh_counts   | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.032 |       0.032 |       0.033 |           49.7 |        49.8 |       65.2 | TRUE          |
| duckhts_roh_counts   | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.047 |       0.047 |       0.048 |           52.0 |        52.2 |       72.7 | TRUE          |
| duckhts_roh_counts   | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.025 |       0.020 |       0.034 |           51.4 |        51.8 |       93.2 | TRUE          |
| duckhts_roh_counts   | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.032 |       0.028 |       0.035 |           57.7 |        59.3 |      100.7 | TRUE          |
| duckhts_roh_counts   | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.035 |       0.034 |       0.037 |           62.5 |        62.9 |      108.2 | TRUE          |
| duckhts_roh_counts   | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.252 |       0.249 |       0.253 |           96.2 |        96.2 |      203.0 | TRUE          |
| duckhts_roh_counts   | samples      |     2 |       1 |      16 |    4000000 | 400      |          0.494 |       0.490 |       0.495 |          143.2 |       143.4 |      340.3 | TRUE          |
| duckhts_roh_counts   | samples      |     4 |       1 |      32 |    8000000 | 800      |          0.976 |       0.972 |       0.984 |          237.3 |       238.1 |      614.9 | TRUE          |
| duckhts_roh_counts   | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.137 |       0.136 |       0.141 |          126.7 |       127.7 |      284.5 | TRUE          |
| duckhts_roh_counts   | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.261 |       0.247 |       0.292 |          184.4 |       199.5 |      421.8 | TRUE          |
| duckhts_roh_counts   | samples      |     4 |       4 |      32 |    8000000 | 800      |          0.478 |       0.465 |       0.503 |          298.7 |       299.9 |      696.5 | TRUE          |
| duckhts_roh_counts   | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.146 |       0.146 |       0.149 |          123.5 |       123.6 |      215.8 | TRUE          |
| duckhts_roh_counts   | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.303 |       0.303 |       0.306 |          179.7 |       180.0 |      393.2 | TRUE          |
| duckhts_roh_counts   | sites        |     4 |       1 |       1 |    4000000 | 400      |          0.605 |       0.604 |       0.614 |          318.0 |       318.1 |      748.0 | TRUE          |
| duckhts_roh_counts   | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.146 |       0.143 |       0.153 |          130.7 |       131.0 |      542.0 | TRUE          |
| duckhts_roh_counts   | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.291 |       0.284 |       0.292 |          198.7 |       198.8 |     1045.5 | TRUE          |
| duckhts_roh_counts   | sites        |     4 |       4 |       1 |    4000000 | 400      |          0.568 |       0.565 |       0.572 |          348.5 |       348.7 |     2052.6 | TRUE          |

Time exponents per doubling (log2 of the median ratio) and peak-RSS ratios to 1×:

| macro                | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | one_x_seconds |
|:---------------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|--------------:|
| duckhts_roh_ancestry | joint     |       1 |            0.451 |            0.688 |        1.580 |        2.830 |         3.746 |
| duckhts_roh_ancestry | reference |       1 |            0.093 |            0.181 |        1.536 |        2.581 |         3.746 |
| duckhts_roh_ancestry | samples   |       1 |            0.356 |            0.508 |        1.046 |        1.096 |         4.530 |
| duckhts_roh_counts   | samples   |       1 |            0.972 |            0.983 |        1.489 |        2.467 |         0.252 |
| duckhts_roh_counts   | sites     |       1 |            1.050 |            0.998 |        1.455 |        2.575 |         0.146 |
| duckhts_roh_ancestry | joint     |       4 |            0.478 |            0.685 |        1.719 |        2.870 |         3.679 |
| duckhts_roh_ancestry | reference |       4 |            0.108 |            0.192 |        1.737 |        2.782 |         3.679 |
| duckhts_roh_ancestry | samples   |       4 |            0.345 |            0.518 |        0.980 |        1.032 |         4.528 |
| duckhts_roh_counts   | samples   |       4 |            0.929 |            0.873 |        1.455 |        2.358 |         0.137 |
| duckhts_roh_counts   | sites     |       4 |            0.999 |            0.964 |        1.520 |        2.666 |         0.146 |

## Findings

Both macros stay within their declared budgets at every size and thread count,
synthetic and real, and the render enforces that.

`duckhts_roh_counts`: its one-thread 1× runs are under the 5-second floor, so
it gets a memory verdict, not a timing verdict. On the real workload, one
thread decodes one sample’s 108,757 chr20 sites in
0.019 s at 49 MiB
peak RSS and three samples in 0.047 s at
52 MiB. The operating scale this report supports
is what it measures: per-chromosome lists of up to 4 million synthetic sites,
up to 32 synthetic samples, and three real 30× samples at 108,757 sites;
larger panels and cohorts are not measured here.

`duckhts_roh_ancestry`: memory still grows with the number of samples decoded
together, because one pass over the BCF needs every sample’s sites. Each
sample now adds 16 bytes per reference site. `samples :=` chooses the batch, and
`max_sites` and `max_site_bytes` turn an oversized one into an explicit error.

## Change from the previous revision

The previous revision, a5f3b47a, held each sample’s sites in a
DuckDB list aggregate, kept the records without a reference site in those
lists, built the ancestry reference with ordered list aggregates and gave the
planner no row estimate for `read_bcf`. Its
observations are read from git; they were made with the same driver and inputs
on the same host. One thread, median of three:

| macro                | dimension    | scale | seconds_before | rss_mib_before | seconds_after | rss_mib_after | time_ratio | rss_ratio |
|:---------------------|:-------------|------:|---------------:|---------------:|--------------:|--------------:|-----------:|----------:|
| duckhts_roh_ancestry | joint        |     1 |          6.247 |       1780.504 |         3.746 |       245.227 |      0.600 |     0.138 |
| duckhts_roh_ancestry | joint        |     2 |         10.128 |       3440.414 |         5.121 |       387.270 |      0.506 |     0.113 |
| duckhts_roh_ancestry | joint        |     4 |         18.219 |       6601.320 |         8.249 |       693.922 |      0.453 |     0.105 |
| duckhts_roh_ancestry | reference    |     1 |          6.247 |       1780.504 |         3.746 |       245.227 |      0.600 |     0.138 |
| duckhts_roh_ancestry | reference    |     2 |          6.551 |       2156.871 |         3.996 |       376.586 |      0.610 |     0.175 |
| duckhts_roh_ancestry | reference    |     4 |          7.254 |       2717.016 |         4.530 |       632.918 |      0.625 |     0.233 |
| duckhts_roh_ancestry | samples      |     1 |          7.254 |       2717.016 |         4.530 |       632.918 |      0.625 |     0.233 |
| duckhts_roh_ancestry | samples      |     2 |         10.858 |       4040.809 |         5.799 |       662.188 |      0.534 |     0.164 |
| duckhts_roh_ancestry | samples      |     4 |         18.219 |       6601.320 |         8.249 |       693.922 |      0.453 |     0.105 |
| duckhts_roh_counts   | real_samples |     1 |          0.022 |         54.445 |         0.019 |        48.996 |      0.861 |     0.900 |
| duckhts_roh_counts   | real_samples |     2 |          0.040 |         65.332 |         0.032 |        49.668 |      0.807 |     0.760 |
| duckhts_roh_counts   | real_samples |     3 |          0.060 |         79.270 |         0.047 |        52.047 |      0.784 |     0.657 |
| duckhts_roh_counts   | samples      |     1 |          0.329 |        254.898 |         0.252 |        96.160 |      0.766 |     0.377 |
| duckhts_roh_counts   | samples      |     2 |          0.661 |        466.645 |         0.494 |       143.176 |      0.747 |     0.307 |
| duckhts_roh_counts   | samples      |     4 |          1.337 |        903.488 |         0.976 |       237.281 |      0.730 |     0.263 |
| duckhts_roh_counts   | sites        |     1 |          0.177 |        162.422 |         0.146 |       123.512 |      0.825 |     0.760 |
| duckhts_roh_counts   | sites        |     2 |          0.338 |        310.410 |         0.303 |       179.703 |      0.897 |     0.579 |
| duckhts_roh_counts   | sites        |     4 |          0.702 |        702.176 |         0.605 |       317.961 |      0.861 |     0.453 |

With 32 samples and the full reference, `duckhts_roh_ancestry` went from
6601 MiB to 694 MiB and
from 18.2 s to 8.2 s.
The segments are unchanged: the chr20 evaluation arms reproduce their 51,083
(ancestry) and 60,204 (pooled) segments exactly.
