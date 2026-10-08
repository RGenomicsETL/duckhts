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

Source revision: 3fa45fdfb29a23449fd90b488d82313abea68da3; `src` tree 3014c1ba19dcedff6a22f9a6b19233a25362696f. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: 419c4c2d72f1c786f5f2b6dc221d00f8b98af67236bc24a6b751a9ce3f70f71a. Host load at start: 0.88, 1.04, 0.79. DuckDB memory limit: 50.0 GiB (the default); temporary directory limit: 0.

Overhead: 38.2 MiB. Ceilings: real 152.2 MiB, samples 1342 MiB, sites 1189.6 MiB.

## Inputs

`rows_used` are the rows the fit keeps, `cells` the histogram cells they make and `pairs` the distinct (depth, alt count) pairs, all with the default options.

| workload | scale | samples |     sites | rows_used |     cells |  pairs | budget_mib |
|:---------|------:|--------:|----------:|----------:|----------:|-------:|-----------:|
| sites    |     1 |       1 | 2,000,000 | 2,000,000 |   322,483 |    642 |    1,189.6 |
| sites    |     2 |       1 | 4,000,000 | 4,000,000 |   644,876 |    649 |    1,189.6 |
| sites    |     4 |       1 | 8,000,000 | 8,000,000 | 1,289,659 |    650 |    1,189.6 |
| samples  |     1 |       8 |   250,000 | 2,000,000 |   363,259 |  4,878 |    1,342.0 |
| samples  |     2 |      16 |   250,000 | 4,000,000 |   726,273 |  9,770 |    1,342.0 |
| samples  |     4 |      32 |   250,000 | 8,000,000 | 1,452,666 | 19,553 |    1,342.0 |
| real     |     1 |       3 |   108,738 |   326,214 |   126,676 |  3,124 |      152.2 |

## Measurements

`mean_contamination` is the mean fitted share over the samples of the query, the same in every repetition. Times are the query latency of three repetitions.

| workload | scale | threads | samples |      rows |     cells | mean_contamination | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | max_buffer_mib | budget_mib |
|:---------|------:|--------:|--------:|----------:|----------:|-------------------:|---------------:|------------:|------------:|---------------:|------------:|---------------:|-----------:|
| real     |     1 |       1 |       3 |   326,214 |   126,676 |               0.00 |          3.245 |       3.241 |       3.253 |           57.6 |        57.7 |           32.7 |      152.2 |
| real     |     1 |       4 |       3 |   326,214 |   126,676 |               0.00 |          2.238 |       2.222 |       2.254 |           71.3 |        71.9 |           54.9 |      152.2 |
| samples  |     1 |       1 |       8 | 2,000,000 |   363,259 |               0.02 |          8.519 |       8.514 |       8.533 |           86.3 |        86.5 |           70.9 |    1,342.0 |
| samples  |     2 |       1 |      16 | 4,000,000 |   726,273 |               0.02 |         17.206 |      17.191 |      17.220 |          127.1 |       127.2 |          139.2 |    1,342.0 |
| samples  |     4 |       1 |      32 | 8,000,000 | 1,452,666 |               0.02 |         34.445 |      34.353 |      34.502 |          212.1 |       212.4 |          275.8 |    1,342.0 |
| samples  |     1 |       4 |       8 | 2,000,000 |   363,259 |               0.02 |          6.398 |       6.390 |       6.413 |          202.3 |       214.6 |          360.7 |    1,342.0 |
| samples  |     2 |       4 |      16 | 4,000,000 |   726,273 |               0.02 |         10.261 |       8.915 |      10.501 |          370.7 |       371.3 |          666.9 |    1,342.0 |
| samples  |     4 |       4 |      32 | 8,000,000 | 1,452,666 |               0.02 |         17.359 |      16.308 |      17.491 |          630.9 |       633.8 |        1,175.8 |    1,342.0 |
| sites    |     1 |       1 |       1 | 2,000,000 |   322,483 |               0.02 |          8.149 |       8.139 |       8.162 |           94.0 |        94.1 |           70.8 |    1,189.6 |
| sites    |     2 |       1 |       1 | 4,000,000 |   644,876 |               0.02 |         13.849 |      13.829 |      13.888 |          142.4 |       142.7 |          123.2 |    1,189.6 |
| sites    |     4 |       1 |       1 | 8,000,000 | 1,289,659 |               0.02 |         26.430 |      26.366 |      26.469 |          239.1 |       239.2 |          244.2 |    1,189.6 |
| sites    |     1 |       4 |       1 | 2,000,000 |   322,483 |               0.02 |          8.018 |       8.006 |       8.037 |          130.1 |       136.6 |          185.6 |    1,189.6 |
| sites    |     2 |       4 |       1 | 4,000,000 |   644,876 |               0.02 |         13.654 |      13.627 |      13.684 |          214.5 |       224.5 |          333.6 |    1,189.6 |
| sites    |     4 |       4 |       1 | 8,000,000 | 1,289,659 |               0.02 |         26.001 |      25.950 |      26.317 |          369.1 |       378.3 |          559.2 |    1,189.6 |

Growth of the median time of the synthetic workloads against the rows (sites) or the samples (an exponent of 1 is linear), and peak-RSS ratios to 1×:

| workload | threads | one_x_seconds | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x |
|:---------|--------:|--------------:|-----------------:|-----------------:|-------------:|-------------:|
| samples  |       1 |         8.519 |            1.014 |            1.001 |        1.473 |        2.458 |
| samples  |       4 |         6.398 |            0.682 |            0.758 |        1.832 |        3.119 |
| sites    |       1 |         8.149 |            0.765 |            0.932 |        1.515 |        2.544 |
| sites    |       4 |         8.018 |            0.768 |            0.929 |        1.649 |        2.837 |

## Findings

Every observation is within its ceiling and none spills; every synthetic 1× run on one thread takes at least 5 s and no time exponent is above 1.25. The render enforces all of this, and that every fit of every run reports `ok` and the same estimates in every repetition.

One sample of 2,000,000 sites takes 8.1 s on one thread (26.4 s at 8,000,000 sites) at 239.2 MiB peak RSS; 32 samples of 250,000 sites take 34.4 s on one thread and 17.4 s on four, at 212.4 MiB and 633.8 MiB. Four threads do not speed up one sample’s fit, and spread the fits of 32 samples over the threads. The three real chr20 samples take 3.2 s on one thread at 57.7 MiB.

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
|  0.02 | 108742 |  31.9 |  0.000147 |      0.023466 |         0.000612 |            0.121175 |       0.498203 |        0.003116 |               0.040766 |                   -198.0 | ok     |
|  0.05 | 108751 |  32.0 |  0.000293 |      0.057623 |         0.001196 |            0.120443 |       0.498053 |        0.003568 |               0.092063 |                   -853.8 | ok     |
|  0.20 | 108756 |  32.1 |  0.000343 |      0.220912 |         0.001562 |            0.109416 |       0.496043 |        0.004612 |               0.243869 |                  -4610.0 | ok     |
