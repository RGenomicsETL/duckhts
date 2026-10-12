Count error fit scaling
================

This report gives the 1×/2×/4× scaling evidence `STYLE.md` asks for, for
`duckhts_count_error_fit()`, and its first measurement on real counts. The
function fits read error, contamination and five related numbers of each
sample from allele read counts at sites of known frequency. DuckDB builds a
histogram of (block, frequency bin, depth, alt count) cells; the native code
packs a sample’s cells into one bounded buffer and fits the model on it. So
the work has two parts: the histogram, which grows with the rows, and the
fit, which grows with the cells and not with the rows. No earlier report
measures this function; this is its first baseline.

## Workloads

- **Sites.** One synthetic sample of 2, 4 and 8 million sites, depth 20 to
  40, Hardy-Weinberg genotypes at frequencies 0.15 to 0.85 with a second
  genome at 2%, read error 1e-3. The positions are 100 bases apart, so the
  sample spans 200, 400 and 800 Mb and has 20, 40 and 80 blocks of 10 Mb.
- **Samples.** 8, 16 and 32 synthetic samples of 250,000 sites each, the same
  model with the second genome at 2%, each sample spanning 25 Mb and three
  blocks. The buffers of the samples share one byte budget.
- **Real.** The three staged chr20 count tables of 1000 Genomes samples
  (`roh_counts_chr20_na18507`, `roh_counts_chr20_hg00403`,
  `roh_counts_chr20_hg00188`: `duckhts_somalier_bam_counts` at the chr20
  sites of the ancestry panel, about 32× each), checked against their
  registered identities before anything is measured. One chromosome of one
  sample fits in under the 5 s floor that `STYLE.md` sets for a timing
  verdict, so this workload is measured at its one size, three samples, and
  carries no scaling verdict: it records what the function costs on real
  counts and what it finds in them. A genome is 45 times one chromosome and
  is not measured here.

The synthetic inputs are generated untimed at render time, with every draw a
hash of its coordinates, and the number of histogram cells of each input is
counted untimed: it is the input denominator of the fit.

Every observation is one fresh DuckDB CLI process under `/usr/bin/time`,
pinned to eight cores by the render command and run at scheduling priority
−15, with three repetitions at one thread and at four threads in an order
that alternates sizes. The JSON profiler gives the query latency and DuckDB’s
peak buffer and temporary bytes; `/usr/bin/time` gives the peak RSS.
`max_temp_directory_size` is zero, so a spill fails the render; the memory
limit is DuckDB’s default, which the report records. Each query consumes
every row of the fit. One fit runs on one thread; the fits of a query’s
samples run in the threads that hold their rows, so the four-thread cells
measure both the histogram and that spread of the samples over threads.

## Budget, declared before measuring

The overhead is the median peak RSS of three empty runs that start the CLI,
load the extension and scan nothing.

The ceiling is overhead + 3 × the live state of the largest input of each
workload, and it is a gate on the peak RSS: the render stops if any
observation’s peak RSS exceeds it or if any observation spills. DuckDB’s
peak buffer count is its own accounting of what its buffer manager held at
most; it is reported beside the RSS and not gated. The live state of one
query is, summed over its samples:

- the packed cells, 16 bytes each: held natively by the aggregate, up to
  three times that while its buffer doubles, then by DuckDB as one BLOB per
  sample, then as one sorted copy in the fit with its pair index, 24 bytes
  per cell, plus 88 bytes per distinct (depth, alt) pair and 24 bytes per
  depth up to the deepest, 1,000 at most; the native part of this is what
  `max_cell_bytes` charges;
- DuckDB’s hash aggregate that builds the cells from the rows, 160 bytes
  per cell, and the window that numbers the blocks, 64 bytes per cell, both
  of which DuckDB manages and could spill. The rows themselves stream from
  Parquet into the aggregate and are not held.

The cells are counted from the input before measuring, so the budget is a
number for each workload, written below with the measurements.

The render also stops if a time exponent of a synthetic workload is above
1.25 at either doubling, or if its 1× run takes less than 5 s on one thread,
which is the least that `STYLE.md` accepts for a timing verdict.

``` sh
taskset -c 8-15 Rscript -e 'rmarkdown::render("benchmarks/benchmark_count_error_fit.Rmd")'
```

Source revision: e7667ff353725f27c1dc1c64732239fdd5a9ec21; `src` tree 0340e562109701f322a733450b16be60a5daba10. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: cfcb1394601761f2a5ebd1a0d6745fea6099fc6b13232b388accb88938decf97. Host load at start: 1.27, 1.24, 1.27. DuckDB memory limit: 50.0 GiB (the default); temporary directory limit: 0.

Overhead: 39.5 MiB. Ceilings: real 153.6 MiB, samples 1343.3 MiB, sites 1190.9 MiB.

## Inputs

`rows_used` are the rows the fit keeps, `cells` the histogram cells they make and `pairs` the distinct (depth, alt count) pairs, all with the default options.

| workload | scale | samples |     sites | rows_used |     cells |  pairs | budget_mib |
|:---------|------:|--------:|----------:|----------:|----------:|-------:|-----------:|
| sites    |     1 |       1 | 2,000,000 | 2,000,000 |   322,483 |    642 |    1,190.9 |
| sites    |     2 |       1 | 4,000,000 | 4,000,000 |   644,876 |    649 |    1,190.9 |
| sites    |     4 |       1 | 8,000,000 | 8,000,000 | 1,289,659 |    650 |    1,190.9 |
| samples  |     1 |       8 |   250,000 | 2,000,000 |   363,259 |  4,878 |    1,343.3 |
| samples  |     2 |      16 |   250,000 | 4,000,000 |   726,273 |  9,770 |    1,343.3 |
| samples  |     4 |      32 |   250,000 | 8,000,000 | 1,452,666 | 19,553 |    1,343.3 |
| real     |     1 |       3 |   108,738 |   326,214 |   126,676 |  3,124 |      153.6 |

## Measurements

`mean_contamination` is the mean fitted share over the samples of the query, the same in every repetition. Times are the query latency of three repetitions.

| workload | scale | threads | samples |      rows |     cells | mean_contamination | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | max_buffer_mib | budget_mib |
|:---------|------:|--------:|--------:|----------:|----------:|-------------------:|---------------:|------------:|------------:|---------------:|------------:|---------------:|-----------:|
| real     |     1 |       1 |       3 |   326,214 |   126,676 |               0.00 |          3.291 |       3.278 |       3.302 |           61.3 |        61.4 |           33.3 |      153.6 |
| real     |     1 |       4 |       3 |   326,214 |   126,676 |               0.00 |          2.265 |       2.263 |       2.296 |           79.9 |        82.0 |           60.1 |      153.6 |
| samples  |     1 |       1 |       8 | 2,000,000 |   363,259 |               0.02 |          9.287 |       9.285 |       9.326 |           92.7 |        92.8 |           78.2 |    1,343.3 |
| samples  |     2 |       1 |      16 | 4,000,000 |   726,273 |               0.02 |         18.703 |      18.565 |      18.749 |          134.2 |       134.2 |          150.3 |    1,343.3 |
| samples  |     4 |       1 |      32 | 8,000,000 | 1,452,666 |               0.02 |         36.882 |      36.837 |      37.006 |          225.0 |       225.4 |          297.9 |    1,343.3 |
| samples  |     1 |       4 |       8 | 2,000,000 |   363,259 |               0.02 |          7.086 |       7.004 |       7.127 |          218.7 |       220.7 |          369.9 |    1,343.3 |
| samples  |     2 |       4 |      16 | 4,000,000 |   726,273 |               0.02 |         10.291 |       9.736 |      12.106 |          403.9 |       404.4 |          742.4 |    1,343.3 |
| samples  |     4 |       4 |      32 | 8,000,000 | 1,452,666 |               0.02 |         17.265 |      17.258 |      20.157 |          633.7 |       705.4 |        1,302.0 |    1,343.3 |
| sites    |     1 |       1 |       1 | 2,000,000 |   322,483 |               0.02 |          9.139 |       9.118 |       9.178 |          101.3 |       101.4 |           77.7 |    1,190.9 |
| sites    |     2 |       1 |       1 | 4,000,000 |   644,876 |               0.02 |         16.029 |      16.020 |      16.031 |          154.3 |       154.4 |          140.9 |    1,190.9 |
| sites    |     4 |       1 |       1 | 8,000,000 | 1,289,659 |               0.02 |         29.903 |      29.796 |      29.995 |          260.7 |       260.8 |          280.4 |    1,190.9 |
| sites    |     1 |       4 |       1 | 2,000,000 |   322,483 |               0.02 |          9.076 |       8.992 |       9.140 |          136.7 |       144.5 |          189.7 |    1,190.9 |
| sites    |     2 |       4 |       1 | 4,000,000 |   644,876 |               0.02 |         15.830 |      15.698 |      16.048 |          229.4 |       229.5 |          352.0 |    1,190.9 |
| sites    |     4 |       4 |       1 | 8,000,000 | 1,289,659 |               0.02 |         29.435 |      29.241 |      29.489 |          393.1 |       409.1 |          559.6 |    1,190.9 |

Growth of the median time of the synthetic workloads against the rows (sites) or the samples (an exponent of 1 is linear), and peak-RSS ratios to 1×:

| workload | threads | one_x_seconds | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x |
|:---------|--------:|--------------:|-----------------:|-----------------:|-------------:|-------------:|
| samples  |       1 |         9.287 |            1.010 |            0.980 |        1.448 |        2.427 |
| samples  |       4 |         7.086 |            0.538 |            0.747 |        1.847 |        2.898 |
| sites    |       1 |         9.139 |            0.811 |            0.900 |        1.523 |        2.574 |
| sites    |       4 |         9.076 |            0.803 |            0.895 |        1.678 |        2.876 |

## Findings

Every observation is within its ceiling and none spills; every synthetic 1× run on one thread takes at least 5 s and no time exponent is above 1.25. The render enforces all of this, and that every fit of every run reports `ok` and the same estimates in every repetition.

One sample of 2,000,000 sites takes 9.1 s on one thread (29.9 s at 8,000,000 sites) at 260.8 MiB peak RSS; 32 samples of 250,000 sites take 36.9 s on one thread and 17.3 s on four, at 225.4 MiB and 705.4 MiB. Four threads do not speed up one sample’s fit, and spread the fits of 32 samples over the threads. The three real chr20 samples take 3.3 s on one thread at 61.4 MiB.

## The real samples

The fit of the three staged chr20 samples with the default options, untimed. These are single 30× samples with no second genome known, so `contamination` is the share the model attributes to one; `seq_error` is the chance that a read shows the other panel allele at a homozygous site, with site artefacts that `artefact_weight` does not absorb.

| sample  |  sites | depth | blocks | seq_error | contamination | contamination_sd | allele_balance | spread_hom | spread_het | artefact_weight | contamination_relative | relative_minus_unrelated | status |
|:--------|-------:|------:|-------:|----------:|--------------:|-----------------:|---------------:|-----------:|-----------:|----------------:|-----------------------:|-------------------------:|:-------|
| HG00188 | 108755 |  31.6 |      7 |  0.000119 |      0.000326 |          6.6e-05 |       0.497437 |   0.009021 |   0.000835 |        0.001825 |               0.000654 |                      0.1 | ok     |
| HG00403 | 108736 |  31.9 |      7 |  0.000129 |      0.000319 |          6.6e-05 |       0.497970 |   0.012590 |   0.000186 |        0.002955 |               0.000639 |                      0.1 | ok     |
| NA18507 | 108723 |  33.0 |      7 |  0.000141 |      0.000263 |          7.4e-05 |       0.498413 |   0.013956 |   0.000914 |        0.002657 |               0.000527 |                      0.1 | ok     |

## A second genome of another population

The same function on a count-level mixture of two of those samples, untimed: HG00403 (CHS) receives reads of NA18507 (YRI) at a share of 0, 0.02, 0.05 and 0.20, each read of the receiver kept with chance 1 − share and each read of the other sample added with chance share, with the receiver’s panel frequencies. The fitted `contamination` is the number to compare with `share`. The other sample’s genotypes do not follow the receiver’s frequencies, which the unrelated model takes them to, so the fit attributes more reads to it than the share: this is the bias the catalog entry of `duckhts_count_error_fit()` names. The sites common to both tables are used, so the share-0 row is the receiver alone on those sites and differs slightly from its row above.

| share |  sites | depth | seq_error | contamination | contamination_sd | homozygosity_excess | allele_balance | artefact_weight | contamination_relative | relative_minus_unrelated | status |
|------:|-------:|------:|----------:|--------------:|-----------------:|--------------------:|---------------:|----------------:|-----------------------:|-------------------------:|:-------|
|  0.00 | 108736 |  31.9 |  0.000129 |      0.000319 |         0.000066 |            0.121353 |       0.497970 |        0.002955 |               0.000639 |                      0.1 | ok     |
|  0.02 | 108742 |  31.9 |  0.000147 |      0.023461 |         0.000614 |            0.121221 |       0.498188 |        0.003131 |               0.040744 |                   -198.0 | ok     |
|  0.05 | 108751 |  32.0 |  0.000295 |      0.057636 |         0.001187 |            0.120789 |       0.498108 |        0.003545 |               0.092059 |                   -853.7 | ok     |
|  0.20 | 108756 |  32.1 |  0.000343 |      0.220912 |         0.001562 |            0.109416 |       0.496043 |        0.004612 |               0.243859 |                  -4610.0 | ok     |

## Change from the previous render

The previous render, 14c810a2, measured the same cells with the
same driver on the same host and cores. Host load differs between renders.
Median of the repetitions:

| workload | scale | threads | seconds_before | rss_mib_before | seconds_after | rss_mib_after | time_ratio | rss_ratio |
|:---------|------:|--------:|---------------:|---------------:|--------------:|--------------:|-----------:|----------:|
| real     |     1 |       1 |          3.267 |         60.387 |         3.291 |        61.273 |      1.008 |     1.015 |
| real     |     1 |       4 |          2.264 |         77.156 |         2.265 |        79.863 |      1.000 |     1.035 |
| samples  |     1 |       1 |          8.596 |         91.445 |         9.287 |        92.699 |      1.080 |     1.014 |
| samples  |     1 |       4 |          6.410 |        223.066 |         7.086 |       218.695 |      1.105 |     0.980 |
| samples  |     2 |       1 |         17.301 |        133.109 |        18.703 |       134.227 |      1.081 |     1.008 |
| samples  |     2 |       4 |         11.272 |        401.102 |        10.291 |       403.914 |      0.913 |     1.007 |
| samples  |     4 |       1 |         34.681 |        224.156 |        36.882 |       225.027 |      1.063 |     1.004 |
| samples  |     4 |       4 |         19.251 |        617.016 |        17.265 |       633.691 |      0.897 |     1.027 |
| sites    |     1 |       1 |          8.185 |         99.961 |         9.139 |       101.332 |      1.116 |     1.014 |
| sites    |     1 |       4 |          8.076 |        134.055 |         9.076 |       136.668 |      1.124 |     1.019 |
| sites    |     2 |       1 |         13.992 |        153.297 |        16.029 |       154.344 |      1.146 |     1.007 |
| sites    |     2 |       4 |         13.867 |        233.477 |        15.830 |       229.359 |      1.142 |     0.982 |
| sites    |     4 |       1 |         26.689 |        259.844 |        29.903 |       260.664 |      1.120 |     1.003 |
| sites    |     4 |       4 |         26.156 |        393.004 |        29.435 |       393.098 |      1.125 |     1.000 |
