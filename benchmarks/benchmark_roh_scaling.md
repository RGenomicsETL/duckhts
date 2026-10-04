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

Source revision: 57b54cec465a1ff67d64ae9cbaf204c1cef1abaa; `src` tree d11a52b6ef4052aa5cb12df9d11978d4fc6998dc. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: 18aeaf0535464a7794a085c3464102615098cc9d261a36bf2cec64243dd22e0f. Host load before the run: 5.14, 2.97, 1.34. Empty-run overhead: 37.9 MiB. Every observation is in [`benchmark_roh_scaling_runs.csv`](benchmark_roh_scaling_runs.csv).

| macro                | dimension    | scale | threads | samples | input_rows | runs_out | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib | within_budget |
|:---------------------|:-------------|------:|--------:|--------:|-----------:|:---------|---------------:|------------:|------------:|---------------:|------------:|-----------:|:--------------|
| duckhts_roh_ancestry | joint        |     1 |       1 |       8 |    6793144 | 356      |          3.733 |       3.715 |       3.741 |          244.1 |       244.5 |      468.8 | TRUE          |
| duckhts_roh_ancestry | joint        |     2 |       1 |      16 |   13586288 | 1256     |          5.170 |       5.155 |       5.323 |          381.5 |       382.2 |      609.1 | TRUE          |
| duckhts_roh_ancestry | joint        |     4 |       1 |      32 |   27172576 | 5120     |          8.239 |       8.190 |       8.269 |          692.4 |       694.6 |      950.7 | TRUE          |
| duckhts_roh_ancestry | joint        |     1 |       4 |       8 |    6793144 | 356      |          3.692 |       3.692 |       3.700 |          277.1 |       279.5 |      477.9 | TRUE          |
| duckhts_roh_ancestry | joint        |     2 |       4 |      16 |   13586288 | 1256     |          5.113 |       5.067 |       5.125 |          479.6 |       488.7 |      627.2 | TRUE          |
| duckhts_roh_ancestry | joint        |     4 |       4 |      32 |   27172576 | 5120     |          8.188 |       8.125 |       8.282 |          791.8 |       809.8 |      986.9 | TRUE          |
| duckhts_roh_ancestry | reference    |     1 |       1 |       8 |    6793144 | 356      |          3.733 |       3.715 |       3.741 |          244.1 |       244.5 |      468.8 | TRUE          |
| duckhts_roh_ancestry | reference    |     2 |       1 |       8 |    6793144 | 641      |          3.989 |       3.969 |       4.065 |          370.9 |       374.3 |      588.8 | TRUE          |
| duckhts_roh_ancestry | reference    |     4 |       1 |       8 |    6793144 | 1337     |          4.505 |       4.500 |       4.621 |          630.4 |       631.2 |      828.7 | TRUE          |
| duckhts_roh_ancestry | reference    |     1 |       4 |       8 |    6793144 | 356      |          3.692 |       3.692 |       3.700 |          277.1 |       279.5 |      477.9 | TRUE          |
| duckhts_roh_ancestry | reference    |     2 |       4 |       8 |    6793144 | 641      |          3.958 |       3.919 |       4.027 |          477.9 |       488.8 |      606.9 | TRUE          |
| duckhts_roh_ancestry | reference    |     4 |       4 |       8 |    6793144 | 1337     |          4.505 |       4.479 |       4.551 |          813.9 |       816.4 |      864.9 | TRUE          |
| duckhts_roh_ancestry | samples      |     1 |       1 |       8 |    6793144 | 1337     |          4.505 |       4.500 |       4.621 |          630.4 |       631.2 |      828.7 | TRUE          |
| duckhts_roh_ancestry | samples      |     2 |       1 |      16 |   13586288 | 2637     |          5.758 |       5.750 |       5.951 |          656.9 |       657.7 |      869.3 | TRUE          |
| duckhts_roh_ancestry | samples      |     4 |       1 |      32 |   27172576 | 5120     |          8.239 |       8.190 |       8.269 |          692.4 |       694.6 |      950.7 | TRUE          |
| duckhts_roh_ancestry | samples      |     1 |       4 |       8 |    6793144 | 1337     |          4.505 |       4.479 |       4.551 |          813.9 |       816.4 |      864.9 | TRUE          |
| duckhts_roh_ancestry | samples      |     2 |       4 |      16 |   13586288 | 2637     |          5.829 |       5.717 |       5.863 |          781.3 |       789.2 |      905.5 | TRUE          |
| duckhts_roh_ancestry | samples      |     4 |       4 |      32 |   27172576 | 5120     |          8.188 |       8.125 |       8.282 |          791.8 |       809.8 |      986.9 | TRUE          |
| duckhts_roh_counts   | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.019 |       0.019 |       0.020 |           48.9 |        49.0 |       57.2 | TRUE          |
| duckhts_roh_counts   | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.033 |       0.031 |       0.033 |           49.6 |        49.7 |       64.7 | TRUE          |
| duckhts_roh_counts   | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.047 |       0.047 |       0.048 |           52.2 |        52.2 |       72.1 | TRUE          |
| duckhts_roh_counts   | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.026 |       0.021 |       0.030 |           49.9 |        50.7 |       92.7 | TRUE          |
| duckhts_roh_counts   | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.027 |       0.027 |       0.034 |           56.9 |        61.6 |      100.1 | TRUE          |
| duckhts_roh_counts   | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.034 |       0.032 |       0.035 |           62.1 |        63.1 |      107.6 | TRUE          |
| duckhts_roh_counts   | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.251 |       0.249 |       0.254 |           96.2 |        96.4 |      202.4 | TRUE          |
| duckhts_roh_counts   | samples      |     2 |       1 |      16 |    4000000 | 400      |          0.494 |       0.489 |       0.503 |          142.9 |       143.1 |      339.7 | TRUE          |
| duckhts_roh_counts   | samples      |     4 |       1 |      32 |    8000000 | 800      |          0.983 |       0.979 |       1.015 |          237.6 |       237.7 |      614.4 | TRUE          |
| duckhts_roh_counts   | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.147 |       0.136 |       0.153 |          127.1 |       127.4 |      284.0 | TRUE          |
| duckhts_roh_counts   | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.262 |       0.262 |       0.265 |          184.7 |       184.8 |      421.3 | TRUE          |
| duckhts_roh_counts   | samples      |     4 |       4 |      32 |    8000000 | 800      |          0.519 |       0.455 |       0.533 |          304.7 |       316.2 |      695.9 | TRUE          |
| duckhts_roh_counts   | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.148 |       0.146 |       0.150 |          123.3 |       123.3 |      215.3 | TRUE          |
| duckhts_roh_counts   | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.304 |       0.301 |       0.312 |          179.2 |       179.2 |      392.7 | TRUE          |
| duckhts_roh_counts   | sites        |     4 |       1 |       1 |    4000000 | 400      |          0.606 |       0.604 |       0.615 |          317.8 |       318.0 |      747.4 | TRUE          |
| duckhts_roh_counts   | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.147 |       0.146 |       0.151 |          130.8 |       130.9 |      541.4 | TRUE          |
| duckhts_roh_counts   | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.292 |       0.281 |       0.301 |          198.7 |       198.7 |     1045.0 | TRUE          |
| duckhts_roh_counts   | sites        |     4 |       4 |       1 |    4000000 | 400      |          0.578 |       0.565 |       0.613 |          348.3 |       348.4 |     2052.1 | TRUE          |

Time exponents per doubling (log2 of the median ratio) and peak-RSS ratios to 1×:

| macro                | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | one_x_seconds |
|:---------------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|--------------:|
| duckhts_roh_ancestry | joint     |       1 |            0.470 |            0.672 |        1.563 |        2.837 |         3.733 |
| duckhts_roh_ancestry | reference |       1 |            0.096 |            0.176 |        1.519 |        2.583 |         3.733 |
| duckhts_roh_ancestry | samples   |       1 |            0.354 |            0.517 |        1.042 |        1.098 |         4.505 |
| duckhts_roh_counts   | samples   |       1 |            0.976 |            0.992 |        1.485 |        2.470 |         0.251 |
| duckhts_roh_counts   | sites     |       1 |            1.038 |            0.993 |        1.453 |        2.577 |         0.148 |
| duckhts_roh_ancestry | joint     |       4 |            0.470 |            0.679 |        1.731 |        2.857 |         3.692 |
| duckhts_roh_ancestry | reference |       4 |            0.101 |            0.187 |        1.725 |        2.937 |         3.692 |
| duckhts_roh_ancestry | samples   |       4 |            0.372 |            0.490 |        0.960 |        0.973 |         4.505 |
| duckhts_roh_counts   | samples   |       4 |            0.831 |            0.987 |        1.453 |        2.397 |         0.147 |
| duckhts_roh_counts   | sites     |       4 |            0.993 |            0.985 |        1.519 |        2.663 |         0.147 |

## Findings

Both macros stay within their declared budgets at every size and thread count,
synthetic and real, and the render enforces that.

`duckhts_roh_counts`: its one-thread 1× runs are under the 5-second floor, so
it gets a memory verdict, not a timing verdict. On the real workload, one
thread decodes one sample’s 108,757 chr20 sites in
0.019 s at 48.9 MiB
peak RSS and three samples in 0.047 s at
52.2 MiB. The operating scale this report supports
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
| duckhts_roh_ancestry | joint        |     1 |          6.247 |       1780.504 |         3.733 |       244.129 |      0.597 |     0.137 |
| duckhts_roh_ancestry | joint        |     2 |         10.128 |       3440.414 |         5.170 |       381.492 |      0.510 |     0.111 |
| duckhts_roh_ancestry | joint        |     4 |         18.219 |       6601.320 |         8.239 |       692.359 |      0.452 |     0.105 |
| duckhts_roh_ancestry | reference    |     1 |          6.247 |       1780.504 |         3.733 |       244.129 |      0.597 |     0.137 |
| duckhts_roh_ancestry | reference    |     2 |          6.551 |       2156.871 |         3.989 |       370.883 |      0.609 |     0.172 |
| duckhts_roh_ancestry | reference    |     4 |          7.254 |       2717.016 |         4.505 |       630.449 |      0.621 |     0.232 |
| duckhts_roh_ancestry | samples      |     1 |          7.254 |       2717.016 |         4.505 |       630.449 |      0.621 |     0.232 |
| duckhts_roh_ancestry | samples      |     2 |         10.858 |       4040.809 |         5.758 |       656.945 |      0.530 |     0.163 |
| duckhts_roh_ancestry | samples      |     4 |         18.219 |       6601.320 |         8.239 |       692.359 |      0.452 |     0.105 |
| duckhts_roh_counts   | real_samples |     1 |          0.022 |         54.445 |         0.019 |        48.852 |      0.892 |     0.897 |
| duckhts_roh_counts   | real_samples |     2 |          0.040 |         65.332 |         0.033 |        49.641 |      0.820 |     0.760 |
| duckhts_roh_counts   | real_samples |     3 |          0.060 |         79.270 |         0.047 |        52.160 |      0.794 |     0.658 |
| duckhts_roh_counts   | samples      |     1 |          0.329 |        254.898 |         0.251 |        96.219 |      0.765 |     0.377 |
| duckhts_roh_counts   | samples      |     2 |          0.661 |        466.645 |         0.494 |       142.945 |      0.748 |     0.306 |
| duckhts_roh_counts   | samples      |     4 |          1.337 |        903.488 |         0.983 |       237.617 |      0.736 |     0.263 |
| duckhts_roh_counts   | sites        |     1 |          0.177 |        162.422 |         0.148 |       123.285 |      0.835 |     0.759 |
| duckhts_roh_counts   | sites        |     2 |          0.338 |        310.410 |         0.304 |       179.238 |      0.901 |     0.577 |
| duckhts_roh_counts   | sites        |     4 |          0.702 |        702.176 |         0.606 |       317.848 |      0.862 |     0.453 |

With 32 samples and the full reference, `duckhts_roh_ancestry` went from
6601 MiB to 692 MiB and
from 18.2 s to 8.2 s.
The segments are unchanged: the chr20 evaluation arms reproduce their 51,083
(ancestry) and 60,204 (pooled) segments exactly.
