ROH macro scaling
================

This report gives the 1×/2×/4× scaling evidence `STYLE.md` asks for, for the
runs-of-homozygosity entry points added for \#318 and the per-sample fitted
count decode. Since \#329 the sites of each sample and chromosome are held in
native buffers that DuckHTS bounds, instead of DuckDB lists; the last sections
compare this revision with the list-aggregate revision and with the previous
render of this report.

- `duckhts_roh_counts` over a relation of read counts. Its growth dimensions are
  the sites per sample and chromosome (the length of one kernel list) and the
  number of samples. Synthetic counts give the 1×/2×/4× series. A real public
  workload adds read counts from three 1000 Genomes 30× CRAMs at the 108,757
  chr20 sites of the ROH ancestry evaluation, decoded for 1, 2 and 3 samples.
- `duckhts_roh_counts_fitted` over the same inputs, with a fit relation that
  gives every sample `seq_error` 0.001 and `contamination` 0, the scalar
  defaults. It is measured on the same series as `duckhts_roh_counts`.
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
  - the kernel workspace of each thread for the longest list, 36.25 bytes per
    site (position, two emissions, the forward-backward values and two
    traceback bits).

  That is rows × 24 + threads × sites per list × 36.25 bytes.
- **`duckhts_roh_counts_fitted`.** The same live state as `duckhts_roh_counts`.
  The fit lookup holds one row per sample, and the sample set is a second scan
  of the input that keeps one row per sample, so neither grows with the sites.
- **`duckhts_roh_ancestry`.** A record without a reference site is dropped
  before the sites are collected, so the lists hold reference sites only: the
  reference rows divided by its number of groups. The live state is:
  - the packed site lists, 16 bytes per sample and reference site (position,
    frequency and genotype evidence);
  - the distinct called sites and the hash table that joins the reference to
    them, 64 bytes per BCF record each;
  - the reference join and the per-site frequency lists, 64 bytes per reference
    row;
  - the kernel workspace, 36.25 bytes per reference site and thread.

  That is samples × reference sites × 16 + records × 128 + reference rows × 64 +
  threads × reference sites × 36.25 bytes.

``` sh
taskset -c 8-15 Rscript -e 'rmarkdown::render("benchmarks/benchmark_roh_scaling.Rmd")'
```

Source revision: b6a14e3be4f31e18356b9adde1f8e1d085363ae2; `src` tree b0474026510f692da693fe57fd8d1c56d79cd70b. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: 5fd33730e8af91c365ac3dcf58091e3faf3b2d49c7849ef283a8db8450bf4bfe. Host load before the run: 1.42, 1.29, 1.10. Empty-run overhead: 38.6 MiB. Every observation is in [`benchmark_roh_scaling_runs.csv`](benchmark_roh_scaling_runs.csv).

| macro                     | dimension    | scale | threads | samples | input_rows | runs_out | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib | within_budget |
|:--------------------------|:-------------|------:|--------:|--------:|-----------:|:---------|---------------:|------------:|------------:|---------------:|------------:|-----------:|:--------------|
| duckhts_roh_ancestry      | joint        |     1 |       1 |       8 |    6793144 | 356      |          3.698 |       3.669 |       3.718 |          244.5 |       245.3 |      469.3 | TRUE          |
| duckhts_roh_ancestry      | joint        |     2 |       1 |      16 |   13586288 | 1256     |          5.060 |       5.016 |       5.097 |          385.2 |       387.2 |      609.5 | TRUE          |
| duckhts_roh_ancestry      | joint        |     4 |       1 |      32 |   27172576 | 5120     |          8.205 |       8.166 |       8.248 |          698.5 |       700.0 |      950.8 | TRUE          |
| duckhts_roh_ancestry      | joint        |     1 |       4 |       8 |    6793144 | 356      |          3.661 |       3.656 |       3.673 |          275.2 |       276.6 |      478.0 | TRUE          |
| duckhts_roh_ancestry      | joint        |     2 |       4 |      16 |   13586288 | 1256     |          5.024 |       5.021 |       5.086 |          478.1 |       496.2 |      626.7 | TRUE          |
| duckhts_roh_ancestry      | joint        |     4 |       4 |      32 |   27172576 | 5120     |          7.992 |       7.964 |       8.003 |          786.9 |       800.3 |      985.3 | TRUE          |
| duckhts_roh_ancestry      | reference    |     1 |       1 |       8 |    6793144 | 356      |          3.698 |       3.669 |       3.718 |          244.5 |       245.3 |      469.3 | TRUE          |
| duckhts_roh_ancestry      | reference    |     2 |       1 |       8 |    6793144 | 641      |          3.966 |       3.960 |       3.978 |          372.7 |       372.9 |      589.1 | TRUE          |
| duckhts_roh_ancestry      | reference    |     4 |       1 |       8 |    6793144 | 1337     |          4.487 |       4.468 |       4.491 |          629.4 |       632.5 |      828.8 | TRUE          |
| duckhts_roh_ancestry      | reference    |     1 |       4 |       8 |    6793144 | 356      |          3.661 |       3.656 |       3.673 |          275.2 |       276.6 |      478.0 | TRUE          |
| duckhts_roh_ancestry      | reference    |     2 |       4 |       8 |    6793144 | 641      |          3.924 |       3.910 |       3.926 |          475.7 |       478.8 |      606.4 | TRUE          |
| duckhts_roh_ancestry      | reference    |     4 |       4 |       8 |    6793144 | 1337     |          4.448 |       4.421 |       4.470 |          769.4 |       771.3 |      863.3 | TRUE          |
| duckhts_roh_ancestry      | samples      |     1 |       1 |       8 |    6793144 | 1337     |          4.487 |       4.468 |       4.491 |          629.4 |       632.5 |      828.8 | TRUE          |
| duckhts_roh_ancestry      | samples      |     2 |       1 |      16 |   13586288 | 2637     |          5.729 |       5.702 |       5.737 |          647.4 |       658.1 |      869.4 | TRUE          |
| duckhts_roh_ancestry      | samples      |     4 |       1 |      32 |   27172576 | 5120     |          8.205 |       8.166 |       8.248 |          698.5 |       700.0 |      950.8 | TRUE          |
| duckhts_roh_ancestry      | samples      |     1 |       4 |       8 |    6793144 | 1337     |          4.448 |       4.421 |       4.470 |          769.4 |       771.3 |      863.3 | TRUE          |
| duckhts_roh_ancestry      | samples      |     2 |       4 |      16 |   13586288 | 2637     |          5.665 |       5.655 |       5.677 |          748.7 |       762.9 |      904.0 | TRUE          |
| duckhts_roh_ancestry      | samples      |     4 |       4 |      32 |   27172576 | 5120     |          7.992 |       7.964 |       8.003 |          786.9 |       800.3 |      985.3 | TRUE          |
| duckhts_roh_counts        | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.019 |       0.019 |       0.019 |           49.3 |        49.3 |       57.3 | TRUE          |
| duckhts_roh_counts        | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.033 |       0.032 |       0.033 |           51.4 |        51.6 |       64.8 | TRUE          |
| duckhts_roh_counts        | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.048 |       0.047 |       0.049 |           52.9 |        53.1 |       72.2 | TRUE          |
| duckhts_roh_counts        | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.034 |       0.020 |       0.035 |           50.4 |        51.7 |       91.1 | TRUE          |
| duckhts_roh_counts        | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.030 |       0.020 |       0.033 |           57.7 |        60.2 |       98.6 | TRUE          |
| duckhts_roh_counts        | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.034 |       0.033 |       0.035 |           63.1 |        65.1 |      106.1 | TRUE          |
| duckhts_roh_counts        | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.253 |       0.253 |       0.253 |           96.0 |        96.1 |      201.8 | TRUE          |
| duckhts_roh_counts        | samples      |     2 |       1 |      16 |    4000000 | 400      |          0.499 |       0.496 |       0.503 |          142.9 |       143.1 |      339.1 | TRUE          |
| duckhts_roh_counts        | samples      |     4 |       1 |      32 |    8000000 | 800      |          0.991 |       0.986 |       0.992 |          236.6 |       236.8 |      613.8 | TRUE          |
| duckhts_roh_counts        | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.140 |       0.139 |       0.145 |          136.1 |       142.3 |      279.6 | TRUE          |
| duckhts_roh_counts        | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.270 |       0.269 |       0.271 |          210.2 |       210.5 |      416.9 | TRUE          |
| duckhts_roh_counts        | samples      |     4 |       4 |      32 |    8000000 | 800      |          0.473 |       0.464 |       0.496 |          335.5 |       341.0 |      691.6 | TRUE          |
| duckhts_roh_counts        | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.150 |       0.144 |       0.151 |          121.8 |       121.9 |      210.9 | TRUE          |
| duckhts_roh_counts        | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.306 |       0.305 |       0.307 |          186.8 |       186.8 |      383.3 | TRUE          |
| duckhts_roh_counts        | sites        |     4 |       1 |       1 |    4000000 | 400      |          0.612 |       0.611 |       0.613 |          333.7 |       333.8 |      728.1 | TRUE          |
| duckhts_roh_counts        | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.144 |       0.144 |       0.148 |          129.2 |       129.3 |      522.1 | TRUE          |
| duckhts_roh_counts        | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.304 |       0.303 |       0.309 |          195.6 |       198.5 |     1005.6 | TRUE          |
| duckhts_roh_counts        | sites        |     4 |       4 |       1 |    4000000 | 400      |          0.592 |       0.587 |       0.599 |          344.8 |       345.4 |     1972.6 | TRUE          |
| duckhts_roh_counts_fitted | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.022 |       0.022 |       0.023 |           55.7 |        56.0 |       57.3 | TRUE          |
| duckhts_roh_counts_fitted | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.038 |       0.038 |       0.038 |           58.6 |        58.7 |       64.8 | TRUE          |
| duckhts_roh_counts_fitted | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.054 |       0.053 |       0.055 |           61.1 |        61.3 |       72.2 | TRUE          |
| duckhts_roh_counts_fitted | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.025 |       0.024 |       0.028 |           55.4 |        55.8 |       91.1 | TRUE          |
| duckhts_roh_counts_fitted | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.026 |       0.025 |       0.029 |           65.9 |        66.1 |       98.6 | TRUE          |
| duckhts_roh_counts_fitted | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.038 |       0.033 |       0.041 |           70.0 |        73.3 |      106.1 | TRUE          |
| duckhts_roh_counts_fitted | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.291 |       0.289 |       0.292 |          104.0 |       105.3 |      201.8 | TRUE          |
| duckhts_roh_counts_fitted | samples      |     2 |       1 |      16 |    4000000 | 400      |          0.570 |       0.563 |       0.571 |          159.0 |       161.4 |      339.1 | TRUE          |
| duckhts_roh_counts_fitted | samples      |     4 |       1 |      32 |    8000000 | 800      |          1.126 |       1.125 |       1.128 |          266.1 |       266.1 |      613.8 | TRUE          |
| duckhts_roh_counts_fitted | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.194 |       0.191 |       0.194 |          150.3 |       156.1 |      279.6 | TRUE          |
| duckhts_roh_counts_fitted | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.332 |       0.301 |       0.348 |          215.5 |       219.1 |      416.9 | TRUE          |
| duckhts_roh_counts_fitted | samples      |     4 |       4 |      32 |    8000000 | 800      |          0.595 |       0.529 |       0.604 |          391.5 |       393.2 |      691.6 | TRUE          |
| duckhts_roh_counts_fitted | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.170 |       0.169 |       0.175 |          127.1 |       127.4 |      210.9 | TRUE          |
| duckhts_roh_counts_fitted | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.348 |       0.343 |       0.349 |          192.7 |       192.8 |      383.3 | TRUE          |
| duckhts_roh_counts_fitted | sites        |     4 |       1 |       1 |    4000000 | 400      |          0.688 |       0.687 |       0.698 |          338.6 |       338.7 |      728.1 | TRUE          |
| duckhts_roh_counts_fitted | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.156 |       0.154 |       0.160 |          137.1 |       137.8 |      522.1 | TRUE          |
| duckhts_roh_counts_fitted | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.322 |       0.309 |       0.331 |          203.0 |       204.4 |     1005.6 | TRUE          |
| duckhts_roh_counts_fitted | sites        |     4 |       4 |       1 |    4000000 | 400      |          0.618 |       0.612 |       0.633 |          345.9 |       348.2 |     1972.6 | TRUE          |

Time exponents per doubling (log2 of the median ratio) and peak-RSS ratios to 1×:

| macro                     | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | one_x_seconds |
|:--------------------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|--------------:|
| duckhts_roh_ancestry      | joint     |       1 |            0.452 |            0.698 |        1.575 |        2.857 |         3.698 |
| duckhts_roh_ancestry      | reference |       1 |            0.101 |            0.178 |        1.524 |        2.574 |         3.698 |
| duckhts_roh_ancestry      | samples   |       1 |            0.353 |            0.518 |        1.029 |        1.110 |         4.487 |
| duckhts_roh_counts        | samples   |       1 |            0.981 |            0.990 |        1.489 |        2.465 |         0.253 |
| duckhts_roh_counts_fitted | samples   |       1 |            0.971 |            0.984 |        1.529 |        2.559 |         0.291 |
| duckhts_roh_counts        | sites     |       1 |            1.033 |            0.999 |        1.534 |        2.740 |         0.150 |
| duckhts_roh_counts_fitted | sites     |       1 |            1.037 |            0.982 |        1.516 |        2.664 |         0.170 |
| duckhts_roh_ancestry      | joint     |       4 |            0.457 |            0.670 |        1.737 |        2.859 |         3.661 |
| duckhts_roh_ancestry      | reference |       4 |            0.100 |            0.181 |        1.729 |        2.796 |         3.661 |
| duckhts_roh_ancestry      | samples   |       4 |            0.349 |            0.496 |        0.973 |        1.023 |         4.448 |
| duckhts_roh_counts        | samples   |       4 |            0.952 |            0.808 |        1.544 |        2.465 |         0.140 |
| duckhts_roh_counts_fitted | samples   |       4 |            0.779 |            0.842 |        1.434 |        2.605 |         0.194 |
| duckhts_roh_counts        | sites     |       4 |            1.078 |            0.960 |        1.514 |        2.669 |         0.144 |
| duckhts_roh_counts_fitted | sites     |       4 |            1.050 |            0.940 |        1.481 |        2.523 |         0.156 |

## Findings

Both macros stay within their declared budgets at every size and thread count,
synthetic and real, and the render enforces that.

`duckhts_roh_counts`: its one-thread 1× runs are under the 5-second floor, so
it gets a memory verdict, not a timing verdict. On the real workload, one
thread decodes one sample’s 108,757 chr20 sites in
0.019 s at 49.3 MiB
peak RSS and three samples in 0.048 s at
52.9 MiB. The operating scale this report supports
is what it measures: per-chromosome lists of up to 4 million synthetic sites,
up to 32 synthetic samples, and three real 30× samples at 108,757 sites;
larger panels and cohorts are not measured here.

`duckhts_roh_counts_fitted` returns the same runs as `duckhts_roh_counts` in
every cell, which the render checks, and stays within the same budget. Its
peak RSS over the scalar decode at the largest inputs, one thread, median of
three: 4.9 MiB with 4 million sites of one sample and
29.5 MiB with 32 samples of 250,000 sites.

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
| duckhts_roh_ancestry | joint        |     1 |          6.247 |       1780.504 |         3.698 |       244.484 |      0.592 |     0.137 |
| duckhts_roh_ancestry | joint        |     2 |         10.128 |       3440.414 |         5.060 |       385.195 |      0.500 |     0.112 |
| duckhts_roh_ancestry | joint        |     4 |         18.219 |       6601.320 |         8.205 |       698.496 |      0.450 |     0.106 |
| duckhts_roh_ancestry | reference    |     1 |          6.247 |       1780.504 |         3.698 |       244.484 |      0.592 |     0.137 |
| duckhts_roh_ancestry | reference    |     2 |          6.551 |       2156.871 |         3.966 |       372.727 |      0.605 |     0.173 |
| duckhts_roh_ancestry | reference    |     4 |          7.254 |       2717.016 |         4.487 |       629.414 |      0.619 |     0.232 |
| duckhts_roh_ancestry | samples      |     1 |          7.254 |       2717.016 |         4.487 |       629.414 |      0.619 |     0.232 |
| duckhts_roh_ancestry | samples      |     2 |         10.858 |       4040.809 |         5.729 |       647.402 |      0.528 |     0.160 |
| duckhts_roh_ancestry | samples      |     4 |         18.219 |       6601.320 |         8.205 |       698.496 |      0.450 |     0.106 |
| duckhts_roh_counts   | real_samples |     1 |          0.022 |         54.445 |         0.019 |        49.262 |      0.881 |     0.905 |
| duckhts_roh_counts   | real_samples |     2 |          0.040 |         65.332 |         0.033 |        51.418 |      0.820 |     0.787 |
| duckhts_roh_counts   | real_samples |     3 |          0.060 |         79.270 |         0.048 |        52.871 |      0.805 |     0.667 |
| duckhts_roh_counts   | samples      |     1 |          0.329 |        254.898 |         0.253 |        96.043 |      0.769 |     0.377 |
| duckhts_roh_counts   | samples      |     2 |          0.661 |        466.645 |         0.499 |       142.938 |      0.755 |     0.306 |
| duckhts_roh_counts   | samples      |     4 |          1.337 |        903.488 |         0.991 |       236.582 |      0.741 |     0.262 |
| duckhts_roh_counts   | sites        |     1 |          0.177 |        162.422 |         0.150 |       121.785 |      0.844 |     0.750 |
| duckhts_roh_counts   | sites        |     2 |          0.338 |        310.410 |         0.306 |       186.762 |      0.907 |     0.602 |
| duckhts_roh_counts   | sites        |     4 |          0.702 |        702.176 |         0.612 |       333.688 |      0.872 |     0.475 |

With 32 samples and the full reference, `duckhts_roh_ancestry` went from
6601 MiB to 698 MiB and
from 18.2 s to 8.2 s.
The segments are unchanged: the chr20 evaluation arms reproduce their 51,083
(ancestry) and 60,204 (pooled) segments exactly.

## Change from the previous render

The previous render, d0af7e98, measured the same cells with the
same driver on the same host, with the traceback held in two bytes per site.
The fitted decode is new and has no earlier cells. Host load differs between
renders. One thread, median of three:

| macro                | dimension    | scale | seconds_before | rss_mib_before | seconds_after | rss_mib_after | time_ratio | rss_ratio |
|:---------------------|:-------------|------:|---------------:|---------------:|--------------:|--------------:|-----------:|----------:|
| duckhts_roh_ancestry | joint        |     1 |          3.916 |        244.262 |         3.698 |       244.484 |      0.944 |     1.001 |
| duckhts_roh_ancestry | joint        |     2 |          5.378 |        381.645 |         5.060 |       385.195 |      0.941 |     1.009 |
| duckhts_roh_ancestry | joint        |     4 |          8.689 |        697.676 |         8.205 |       698.496 |      0.944 |     1.001 |
| duckhts_roh_ancestry | reference    |     1 |          3.916 |        244.262 |         3.698 |       244.484 |      0.944 |     1.001 |
| duckhts_roh_ancestry | reference    |     2 |          4.242 |        372.129 |         3.966 |       372.727 |      0.935 |     1.002 |
| duckhts_roh_ancestry | reference    |     4 |          4.777 |        627.699 |         4.487 |       629.414 |      0.939 |     1.003 |
| duckhts_roh_ancestry | samples      |     1 |          4.777 |        627.699 |         4.487 |       629.414 |      0.939 |     1.003 |
| duckhts_roh_ancestry | samples      |     2 |          6.081 |        656.496 |         5.729 |       647.402 |      0.942 |     0.986 |
| duckhts_roh_ancestry | samples      |     4 |          8.689 |        697.676 |         8.205 |       698.496 |      0.944 |     1.001 |
| duckhts_roh_counts   | real_samples |     1 |          0.020 |         49.211 |         0.019 |        49.262 |      0.964 |     1.001 |
| duckhts_roh_counts   | real_samples |     2 |          0.035 |         51.223 |         0.033 |        51.418 |      0.927 |     1.004 |
| duckhts_roh_counts   | real_samples |     3 |          0.051 |         52.949 |         0.048 |        52.871 |      0.934 |     0.999 |
| duckhts_roh_counts   | samples      |     1 |          0.277 |         96.129 |         0.253 |        96.043 |      0.913 |     0.999 |
| duckhts_roh_counts   | samples      |     2 |          0.541 |        142.773 |         0.499 |       142.938 |      0.923 |     1.001 |
| duckhts_roh_counts   | samples      |     4 |          1.072 |        236.762 |         0.991 |       236.582 |      0.924 |     0.999 |
| duckhts_roh_counts   | sites        |     1 |          0.160 |        123.488 |         0.150 |       121.785 |      0.936 |     0.986 |
| duckhts_roh_counts   | sites        |     2 |          0.330 |        190.465 |         0.306 |       186.762 |      0.927 |     0.981 |
| duckhts_roh_counts   | sites        |     4 |          0.667 |        340.273 |         0.612 |       333.688 |      0.918 |     0.981 |
