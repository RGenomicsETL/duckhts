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

Source revision: 3efd5720fbea15f0e73eb2695cd6f361cc18b8a9; `src` tree dbc1c3106c860ff6e7e7134a1ffb526c10dd8854. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: d6ac4b7fd9a45b856e84e24062ac2758c5aa4f633438bdececcd8695fdcbabe4. Host load before the run: 6.46, 4.42, 2.84. Empty-run overhead: 38.3 MiB. Every observation is in [`benchmark_roh_scaling_runs.csv`](benchmark_roh_scaling_runs.csv).

| macro                | dimension    | scale | threads | samples | input_rows | runs_out | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib | within_budget |
|:---------------------|:-------------|------:|--------:|--------:|-----------:|:---------|---------------:|------------:|------------:|---------------:|------------:|-----------:|:--------------|
| duckhts_roh_ancestry | joint        |     1 |       1 |       8 |    6793144 | 356      |          3.916 |       3.782 |       3.976 |          244.3 |       244.7 |      469.2 | TRUE          |
| duckhts_roh_ancestry | joint        |     2 |       1 |      16 |   13586288 | 1256     |          5.378 |       5.216 |       5.408 |          381.6 |       383.4 |      609.5 | TRUE          |
| duckhts_roh_ancestry | joint        |     4 |       1 |      32 |   27172576 | 5120     |          8.689 |       8.374 |       8.725 |          697.7 |       699.3 |      951.1 | TRUE          |
| duckhts_roh_ancestry | joint        |     1 |       4 |       8 |    6793144 | 356      |          4.000 |       3.781 |       4.088 |          277.8 |       284.9 |      478.3 | TRUE          |
| duckhts_roh_ancestry | joint        |     2 |       4 |      16 |   13586288 | 1256     |          5.408 |       5.144 |       5.624 |          477.9 |       490.7 |      627.6 | TRUE          |
| duckhts_roh_ancestry | joint        |     4 |       4 |      32 |   27172576 | 5120     |          8.657 |       8.186 |       8.852 |          810.8 |       840.1 |      987.3 | TRUE          |
| duckhts_roh_ancestry | reference    |     1 |       1 |       8 |    6793144 | 356      |          3.916 |       3.782 |       3.976 |          244.3 |       244.7 |      469.2 | TRUE          |
| duckhts_roh_ancestry | reference    |     2 |       1 |       8 |    6793144 | 641      |          4.242 |       4.053 |       4.253 |          372.1 |       372.4 |      589.2 | TRUE          |
| duckhts_roh_ancestry | reference    |     4 |       1 |       8 |    6793144 | 1337     |          4.777 |       4.770 |       4.814 |          627.7 |       627.7 |      829.0 | TRUE          |
| duckhts_roh_ancestry | reference    |     1 |       4 |       8 |    6793144 | 356      |          4.000 |       3.781 |       4.088 |          277.8 |       284.9 |      478.3 | TRUE          |
| duckhts_roh_ancestry | reference    |     2 |       4 |       8 |    6793144 | 641      |          4.186 |       4.026 |       4.203 |          472.6 |       476.6 |      607.3 | TRUE          |
| duckhts_roh_ancestry | reference    |     4 |       4 |       8 |    6793144 | 1337     |          4.848 |       4.577 |       4.863 |          784.9 |       820.0 |      865.3 | TRUE          |
| duckhts_roh_ancestry | samples      |     1 |       1 |       8 |    6793144 | 1337     |          4.777 |       4.770 |       4.814 |          627.7 |       627.7 |      829.0 | TRUE          |
| duckhts_roh_ancestry | samples      |     2 |       1 |      16 |   13586288 | 2637     |          6.081 |       6.069 |       6.125 |          656.5 |       657.8 |      869.7 | TRUE          |
| duckhts_roh_ancestry | samples      |     4 |       1 |      32 |   27172576 | 5120     |          8.689 |       8.374 |       8.725 |          697.7 |       699.3 |      951.1 | TRUE          |
| duckhts_roh_ancestry | samples      |     1 |       4 |       8 |    6793144 | 1337     |          4.848 |       4.577 |       4.863 |          784.9 |       820.0 |      865.3 | TRUE          |
| duckhts_roh_ancestry | samples      |     2 |       4 |      16 |   13586288 | 2637     |          6.110 |       5.752 |       6.159 |          783.0 |       805.0 |      905.9 | TRUE          |
| duckhts_roh_ancestry | samples      |     4 |       4 |      32 |   27172576 | 5120     |          8.657 |       8.186 |       8.852 |          810.8 |       840.1 |      987.3 | TRUE          |
| duckhts_roh_counts   | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.020 |       0.019 |       0.020 |           49.2 |        49.3 |       57.6 | TRUE          |
| duckhts_roh_counts   | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.035 |       0.035 |       0.036 |           51.2 |        51.4 |       65.1 | TRUE          |
| duckhts_roh_counts   | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.051 |       0.051 |       0.057 |           52.9 |        53.1 |       72.5 | TRUE          |
| duckhts_roh_counts   | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.028 |       0.027 |       0.036 |           51.7 |        51.7 |       93.1 | TRUE          |
| duckhts_roh_counts   | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.035 |       0.034 |       0.036 |           58.4 |        58.8 |      100.5 | TRUE          |
| duckhts_roh_counts   | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.036 |       0.036 |       0.036 |           64.9 |        65.7 |      108.0 | TRUE          |
| duckhts_roh_counts   | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.277 |       0.273 |       0.278 |           96.1 |        96.5 |      202.8 | TRUE          |
| duckhts_roh_counts   | samples      |     2 |       1 |      16 |    4000000 | 400      |          0.541 |       0.540 |       0.549 |          142.8 |       143.1 |      340.1 | TRUE          |
| duckhts_roh_counts   | samples      |     4 |       1 |      32 |    8000000 | 800      |          1.072 |       1.067 |       1.077 |          236.8 |       237.2 |      614.8 | TRUE          |
| duckhts_roh_counts   | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.159 |       0.147 |       0.172 |          130.4 |       134.1 |      284.3 | TRUE          |
| duckhts_roh_counts   | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.297 |       0.268 |       0.299 |          228.2 |       232.8 |      421.7 | TRUE          |
| duckhts_roh_counts   | samples      |     4 |       4 |      32 |    8000000 | 800      |          0.542 |       0.506 |       0.564 |          333.1 |       335.3 |      696.3 | TRUE          |
| duckhts_roh_counts   | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.160 |       0.160 |       0.162 |          123.5 |       123.8 |      215.7 | TRUE          |
| duckhts_roh_counts   | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.330 |       0.330 |       0.335 |          190.5 |       190.5 |      393.1 | TRUE          |
| duckhts_roh_counts   | sites        |     4 |       1 |       1 |    4000000 | 400      |          0.667 |       0.665 |       0.667 |          340.3 |       340.5 |      747.8 | TRUE          |
| duckhts_roh_counts   | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.153 |       0.148 |       0.154 |          131.0 |       131.0 |      541.8 | TRUE          |
| duckhts_roh_counts   | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.325 |       0.306 |       0.345 |          199.2 |       201.3 |     1045.4 | TRUE          |
| duckhts_roh_counts   | sites        |     4 |       4 |       1 |    4000000 | 400      |          0.626 |       0.594 |       0.636 |          341.7 |       348.7 |     2052.5 | TRUE          |

Time exponents per doubling (log2 of the median ratio) and peak-RSS ratios to 1×:

| macro                | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | one_x_seconds |
|:---------------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|--------------:|
| duckhts_roh_ancestry | joint     |       1 |            0.458 |            0.692 |        1.562 |        2.856 |         3.916 |
| duckhts_roh_ancestry | reference |       1 |            0.116 |            0.171 |        1.523 |        2.569 |         3.916 |
| duckhts_roh_ancestry | samples   |       1 |            0.348 |            0.515 |        1.046 |        1.112 |         4.777 |
| duckhts_roh_counts   | samples   |       1 |            0.964 |            0.988 |        1.486 |        2.464 |         0.277 |
| duckhts_roh_counts   | sites     |       1 |            1.045 |            1.013 |        1.543 |        2.755 |         0.160 |
| duckhts_roh_ancestry | joint     |       4 |            0.435 |            0.679 |        1.720 |        2.919 |         4.000 |
| duckhts_roh_ancestry | reference |       4 |            0.065 |            0.212 |        1.701 |        2.825 |         4.000 |
| duckhts_roh_ancestry | samples   |       4 |            0.334 |            0.503 |        0.998 |        1.033 |         4.848 |
| duckhts_roh_counts   | samples   |       4 |            0.900 |            0.870 |        1.750 |        2.554 |         0.159 |
| duckhts_roh_counts   | sites     |       4 |            1.082 |            0.945 |        1.521 |        2.608 |         0.153 |

## Findings

Both macros stay within their declared budgets at every size and thread count,
synthetic and real, and the render enforces that.

`duckhts_roh_counts`: its one-thread 1× runs are under the 5-second floor, so
it gets a memory verdict, not a timing verdict. On the real workload, one
thread decodes one sample’s 108,757 chr20 sites in
0.02 s at 49.2 MiB
peak RSS and three samples in 0.051 s at
52.9 MiB. The operating scale this report supports
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
| duckhts_roh_ancestry | joint        |     1 |          6.247 |       1780.504 |         3.916 |       244.262 |      0.627 |     0.137 |
| duckhts_roh_ancestry | joint        |     2 |         10.128 |       3440.414 |         5.378 |       381.645 |      0.531 |     0.111 |
| duckhts_roh_ancestry | joint        |     4 |         18.219 |       6601.320 |         8.689 |       697.676 |      0.477 |     0.106 |
| duckhts_roh_ancestry | reference    |     1 |          6.247 |       1780.504 |         3.916 |       244.262 |      0.627 |     0.137 |
| duckhts_roh_ancestry | reference    |     2 |          6.551 |       2156.871 |         4.242 |       372.129 |      0.648 |     0.173 |
| duckhts_roh_ancestry | reference    |     4 |          7.254 |       2717.016 |         4.777 |       627.699 |      0.659 |     0.231 |
| duckhts_roh_ancestry | samples      |     1 |          7.254 |       2717.016 |         4.777 |       627.699 |      0.659 |     0.231 |
| duckhts_roh_ancestry | samples      |     2 |         10.858 |       4040.809 |         6.081 |       656.496 |      0.560 |     0.162 |
| duckhts_roh_ancestry | samples      |     4 |         18.219 |       6601.320 |         8.689 |       697.676 |      0.477 |     0.106 |
| duckhts_roh_counts   | real_samples |     1 |          0.022 |         54.445 |         0.020 |        49.211 |      0.914 |     0.904 |
| duckhts_roh_counts   | real_samples |     2 |          0.040 |         65.332 |         0.035 |        51.223 |      0.885 |     0.784 |
| duckhts_roh_counts   | real_samples |     3 |          0.060 |         79.270 |         0.051 |        52.949 |      0.862 |     0.668 |
| duckhts_roh_counts   | samples      |     1 |          0.329 |        254.898 |         0.277 |        96.129 |      0.843 |     0.377 |
| duckhts_roh_counts   | samples      |     2 |          0.661 |        466.645 |         0.541 |       142.773 |      0.818 |     0.306 |
| duckhts_roh_counts   | samples      |     4 |          1.337 |        903.488 |         1.072 |       236.762 |      0.802 |     0.262 |
| duckhts_roh_counts   | sites        |     1 |          0.177 |        162.422 |         0.160 |       123.488 |      0.902 |     0.760 |
| duckhts_roh_counts   | sites        |     2 |          0.338 |        310.410 |         0.330 |       190.465 |      0.978 |     0.614 |
| duckhts_roh_counts   | sites        |     4 |          0.702 |        702.176 |         0.667 |       340.273 |      0.949 |     0.485 |

With 32 samples and the full reference, `duckhts_roh_ancestry` went from
6601 MiB to 698 MiB and
from 18.2 s to 8.7 s.
The segments are unchanged: the chr20 evaluation arms reproduce their 51,083
(ancestry) and 60,204 (pooled) segments exactly.
