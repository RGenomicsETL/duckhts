Depth fraction lane: scaling and a real-data check
================

This report gives the 1×/2×/4× scaling evidence `STYLE.md` asks for, for the
three macros of the depth fraction lane (`design/depth_fraction_lane.md`,
<https://github.com/RGenomicsETL/duckhts/issues/343>): `duckhts_depth_bins()`,
`duckhts_depth_region_fraction()` and `duckhts_depth_sex_fraction()`. It also
checks them on real low-pass data made from three registered 30× genomes by
read-name hash subsampling. No earlier report measures these macros; this is
their first baseline.

## Workloads

- **Synthetic scaling.** A genome of the GRCh38 shape (chr1 to chr22, X and
  Y lengths) in fixed bins, with a GC share per bin and reads drawn from a
  GC-biased mean, every draw a hash of its coordinates. Samples:
  64, 128 and 256 at 25 kb bins (123,440 bins). Bins: 25, 12.5 and 6.25 kb
  (123,440, 246,880 and 493,760 bins) at 64 samples. Joint: 128 samples at
  12.5 kb. The estimators read the same bins with `expected` set to the mean
  and every bin kept, so their inputs do not depend on the correction. The
  region estimator is also measured at 22, 44 and 88 regions (whole
  autosomes, halves and quarters) on the 1× input.
- **Real.** `depth_fraction_hg00403_buckets`,
  `depth_fraction_na18507_buckets` and `depth_fraction_hg00188_buckets`:
  reads (MAPQ at least 20, no unmapped, secondary, failed, duplicate or
  supplementary flag) in 100 kb bins of chr1 to chr22, X and Y of three
  1000 Genomes 30× CRAMs, split into 300 of 3,000 read-name buckets. Ten
  buckets are a 0.1× read-name subsample, 40 are 0.4× and 100 are 1×. The
  three inputs are checked against their registered row and read counts, the
  reference FASTA against its registered size and the sample panel against
  its registered SHA-256 before anything is computed. GC comes from
  `fasta_nuc()` on the registered GRCh38 reference.

Every timed observation is one fresh DuckDB CLI process under
`/usr/bin/time`, pinned to cores 12 to 19 by the render command and run at
scheduling priority −15, with three repetitions at one thread and at four
threads in an order that alternates sizes. The JSON profiler gives the query
latency and DuckDB’s peak buffer and temporary bytes; `/usr/bin/time` gives
the peak RSS. `max_temp_directory_size` is zero, so a spill fails the
render. Each query consumes every output row in an aggregate of all output
columns, so no column can be pruned.

## Budget, declared before measuring

The overhead of a macro is the median peak RSS of three runs of its query on
an empty input of the same schema at four threads: the CLI, the extension and
the plan’s fixed allocations. The ceiling of a workload is overhead + 3 ×
the live state of its largest input, a gate on the peak RSS: the render stops
if any observation exceeds it or spills. The live state, from the plans:

- `duckhts_depth_bins`: the input materialized once, 104 bytes a row (two
  strings of 16, eight numbers of 8 or less and the reason), and the holistic
  aggregates, 8 bytes a value for each of the two quantiles, the stratum
  medians and the centre over the autosomal rows: 136 bytes a row.
- `duckhts_depth_region_fraction`: the narrow columns materialized once, 96
  bytes a row; the joins build on regions and samples × regions.
- `duckhts_depth_sex_fraction`: it streams; its state is a few sums per
  sample. Its live state is taken as the decoded row groups being read, one
  of 122,880 rows of 64 bytes for each thread.

The render also stops if a time exponent is above 1.25 at either doubling of
a series whose 1× run takes at least 5 s on one thread. A series whose 1× run
is shorter gets its memory verdict and its times without a timing verdict, as
`STYLE.md` allows.

## Real-data checks, declared before running

The three genomes are male in the 1000 Genomes sample panel. Composites are
sums of disjoint bucket sets at 0.1×, 0.4× and 1× (10, 40 and 100 buckets):

- **Pure:** each genome alone. Expected: major genome XY, a second-genome
  fraction of 0, and no change of any autosome.
- **Planted gain:** HG00403 with extra buckets on chr18 only, so that chr18
  has the depth of a one-copy gain in a share m of cells (one extra bucket in
  ten is m = 0.2). Expected: shift m on chr18.
- **XY major with an XX second genome:** HG00403 at share 1 − f, and NA18507
  at share f on the autosomes with twice as many buckets on X and none on Y,
  which emulates an XX genome from a male one. Expected: major XY, fraction f.
- **XX major with an XY second genome:** NA18507 emulated as XX at share
  1 − f (autosomes, X at twice the buckets, no Y) with HG00403 at share f.
  Expected: major XX, fraction f. The emulated XX genome has no Y reads at
  all, while a real XX genome has some mismapped Y reads; this check is
  optimistic for the Y estimate.

Bins: 1 Mb, with fewer than 5% N bases, outside the pseudo-autosomal regions
of X (GRCh38 chrX:10,001-2,781,479 and 155,701,383-156,030,895). Bins are
masked and estimates calibrated with a panel of the genomes not used in a
composite (leave-one-out): each genome’s 300 buckets (3×) give corrected
bins; a bin is kept for a composite when the median over the panel genomes of
its ratio divided by its expected level (1 for an autosome, 0.5 for X or Y in
an XY genome) is within \[0.75, 1.25\]. This removes centromere models,
unmappable and multi-copy bins without a separate mappability track. The
same masked panel bins are the `reference_table` of the estimators, all with
sex XY. Acceptance, with the reference:

- every major genome as expected;
- pure genomes: the absolute fraction at most 0.03, and the absolute shift of
  every autosome at most 0.08;
- planted gains: the shift within 0.06 of m at 0.1× and within 0.04 at 0.4×
  and 1×;
- mixtures: the fraction within 0.04 of f at 0.1× and within 0.03 at 0.4× and
  1×.

The bin-based estimates without a reference are reported beside them and are
not gated: residual bias between chromosomes and between X and the autosomes
is what the reference set exists to remove.

A first run of these checks staged chr18 to chr22, X and Y only. It recovered
0.6 to 0.7 m for the planted gains: chr18 then held a third of the autosomal
bins, and most of the bins of the low-GC strata, so the GC correction of the
sample absorbed much of the gain. This render stages the whole genome, and
also reports the planted gains corrected on chr18 to chr22 alone (mode
`five_chromosomes`, not gated) to show the effect.

``` sh
taskset -c 12-19 Rscript -e 'rmarkdown::render("benchmarks/benchmark_depth_fraction.Rmd")'
```

Source revision: 949f5ab0e9e19531ac7c2a83e20107088d1e769c; `src` tree d1a7d32f271763dbbc7a9c4fe693bdcff7c041c1. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: 2494437603f3ecca63dbe8b5b64ef7fddc6bab433d9d47aad06bf5bbb78e4e50. Host load at start: 2.78, 3.53, 4.89. DuckDB memory limit: 50.0 GiB (the default); temporary directory limit: 0. Cores 12 to 19, priority −15.

## Real data

The staged inputs match their registered identities:

| sample_id |      rows |      reads | max_bucket |
|:----------|----------:|-----------:|-----------:|
| HG00188   | 8,294,029 | 58,745,362 |        299 |
| HG00403   | 8,306,508 | 59,501,002 |        299 |
| NA18507   | 8,323,627 | 61,353,704 |        299 |

1 Mb bins after the N and pseudo-autosomal filters, and the bins each leave-one-out mask keeps (group = the genomes a composite uses):

| group           | autosomes | chrX | chrY |
|:----------------|----------:|-----:|-----:|
| before the mask |     2,660 |  143 |   19 |
| HG00403         |     2,564 |  138 |   11 |
| NA18507         |     2,574 |  138 |   12 |
| HG00188         |     2,571 |  138 |   12 |
| HG00403+NA18507 |     2,565 |  138 |   11 |

Acceptance checks, declared above:

| check                                         | pass |
|:----------------------------------------------|:-----|
| major genome as expected (with reference)     | TRUE |
| pure: \|fraction\| \<= 0.03                   | TRUE |
| pure: \|shift\| \<= 0.08 on every autosome    | TRUE |
| gain: \|shift - m\| within 0.06 / 0.04        | TRUE |
| mixtures: \|fraction - f\| within 0.04 / 0.03 | TRUE |

### Sex-chromosome fraction

`truth` is the share of the second genome; `major` the expected major genome. `mode` is `bins` (the ideal copy levels and the bin interval) or `reference` (levels from the masked panel of the other genomes, all XY).

| kind     | depth | truth | major | mode      | x_bins | y_bins | x_copies | y_copies | major_sex | fraction_x | fraction_y | fraction | fraction_low | fraction_high | agreement_z | status       |
|:---------|:------|------:|:------|:----------|-------:|-------:|---------:|---------:|:----------|-----------:|-----------:|---------:|-------------:|--------------:|------------:|:-------------|
| gain     | 0.1x  |  0.20 | XY    | bins      |    138 |     11 |    0.965 |    1.011 | XY        |     -0.035 |     -0.011 |   -0.033 |       -0.047 |        -0.019 |      -0.818 | ok           |
| gain     | 0.1x  |  0.20 | XY    | reference |    138 |     11 |    0.992 |    1.069 | XY        |     -0.008 |     -0.069 |   -0.012 |       -0.026 |         0.003 |       2.050 | ok           |
| gain     | 0.1x  |  0.40 | XY    | bins      |    138 |     11 |    0.963 |    1.009 | XY        |     -0.037 |     -0.009 |   -0.035 |       -0.051 |        -0.020 |      -0.928 | ok           |
| gain     | 0.1x  |  0.40 | XY    | reference |    138 |     11 |    0.990 |    1.066 | XY        |     -0.010 |     -0.066 |   -0.014 |       -0.030 |         0.002 |       1.714 | ok           |
| gain     | 0.4x  |  0.10 | XY    | bins      |    138 |     11 |    0.973 |    1.010 | XY        |     -0.027 |     -0.010 |   -0.026 |       -0.035 |        -0.017 |      -0.969 | ok           |
| gain     | 0.4x  |  0.10 | XY    | reference |    138 |     11 |    1.000 |    1.067 | XY        |      0.000 |     -0.067 |   -0.004 |       -0.014 |         0.005 |       3.528 | x_y_disagree |
| gain     | 0.4x  |  0.20 | XY    | bins      |    138 |     11 |    0.971 |    1.008 | XY        |     -0.029 |     -0.008 |   -0.027 |       -0.037 |        -0.017 |      -1.066 | ok           |
| gain     | 0.4x  |  0.20 | XY    | reference |    138 |     11 |    0.999 |    1.065 | XY        |     -0.001 |     -0.065 |   -0.006 |       -0.016 |         0.004 |       3.207 | x_y_disagree |
| gain     | 1x    |  0.04 | XY    | bins      |    138 |     11 |    0.975 |    1.001 | XY        |     -0.025 |     -0.001 |   -0.023 |       -0.031 |        -0.016 |      -1.618 | ok           |
| gain     | 1x    |  0.04 | XY    | reference |    138 |     11 |    1.002 |    1.057 | XY        |      0.002 |     -0.057 |   -0.002 |       -0.010 |         0.006 |       3.675 | x_y_disagree |
| gain     | 1x    |  0.10 | XY    | bins      |    138 |     11 |    0.974 |    1.000 | XY        |     -0.026 |      0.000 |   -0.024 |       -0.032 |        -0.016 |      -1.691 | ok           |
| gain     | 1x    |  0.10 | XY    | reference |    138 |     11 |    1.001 |    1.056 | XY        |      0.001 |     -0.056 |   -0.003 |       -0.011 |         0.005 |       3.495 | x_y_disagree |
| gain     | 1x    |  0.20 | XY    | bins      |    138 |     11 |    0.972 |    0.998 | XY        |     -0.028 |      0.002 |   -0.026 |       -0.034 |        -0.017 |      -1.764 | ok           |
| gain     | 1x    |  0.20 | XY    | reference |    138 |     11 |    0.999 |    1.055 | XY        |     -0.001 |     -0.055 |   -0.004 |       -0.013 |         0.005 |       3.071 | x_y_disagree |
| pure     | 0.1x  |  0.00 | XY    | bins      |    138 |     12 |    0.971 |    0.945 | XY        |     -0.029 |      0.055 |   -0.022 |       -0.035 |        -0.008 |      -3.282 | x_y_disagree |
| pure     | 0.1x  |  0.00 | XY    | reference |    138 |     12 |    0.996 |    0.989 | XY        |     -0.004 |      0.011 |   -0.003 |       -0.017 |         0.011 |      -0.563 | ok           |
| pure     | 0.1x  |  0.00 | XY    | bins      |    138 |     11 |    0.969 |    1.014 | XY        |     -0.031 |     -0.014 |   -0.030 |       -0.044 |        -0.016 |      -0.630 | ok           |
| pure     | 0.1x  |  0.00 | XY    | reference |    138 |     11 |    0.996 |    1.071 | XY        |     -0.004 |     -0.071 |   -0.009 |       -0.023 |         0.006 |       2.304 | ok           |
| pure     | 0.1x  |  0.00 | XY    | bins      |    138 |     12 |    0.967 |    0.938 | XY        |     -0.033 |      0.062 |   -0.025 |       -0.040 |        -0.009 |      -3.345 | x_y_disagree |
| pure     | 0.1x  |  0.00 | XY    | reference |    138 |     12 |    0.991 |    0.980 | XY        |     -0.009 |      0.020 |   -0.007 |       -0.023 |         0.009 |      -0.994 | ok           |
| pure     | 0.4x  |  0.00 | XY    | bins      |    138 |     12 |    0.972 |    0.924 | XY        |     -0.028 |      0.076 |   -0.019 |       -0.028 |        -0.010 |      -6.291 | x_y_disagree |
| pure     | 0.4x  |  0.00 | XY    | reference |    138 |     12 |    0.997 |    0.968 | XY        |     -0.003 |      0.032 |    0.000 |       -0.009 |         0.009 |      -2.053 | ok           |
| pure     | 0.4x  |  0.00 | XY    | bins      |    138 |     11 |    0.974 |    1.011 | XY        |     -0.026 |     -0.011 |   -0.025 |       -0.034 |        -0.016 |      -0.845 | ok           |
| pure     | 0.4x  |  0.00 | XY    | reference |    138 |     11 |    1.001 |    1.068 | XY        |      0.001 |     -0.068 |   -0.003 |       -0.013 |         0.006 |       3.712 | x_y_disagree |
| pure     | 0.4x  |  0.00 | XY    | bins      |    138 |     12 |    0.973 |    0.929 | XY        |     -0.027 |      0.071 |   -0.019 |       -0.030 |        -0.008 |      -4.763 | x_y_disagree |
| pure     | 0.4x  |  0.00 | XY    | reference |    138 |     12 |    0.996 |    0.971 | XY        |     -0.004 |      0.029 |   -0.001 |       -0.013 |         0.010 |      -1.553 | ok           |
| pure     | 1x    |  0.00 | XY    | bins      |    138 |     12 |    0.973 |    0.930 | XY        |     -0.027 |      0.070 |   -0.019 |       -0.027 |        -0.011 |      -6.920 | x_y_disagree |
| pure     | 1x    |  0.00 | XY    | reference |    138 |     12 |    0.997 |    0.974 | XY        |     -0.003 |      0.026 |    0.000 |       -0.008 |         0.008 |      -1.966 | ok           |
| pure     | 1x    |  0.00 | XY    | bins      |    138 |     11 |    0.975 |    1.001 | XY        |     -0.025 |     -0.001 |   -0.023 |       -0.031 |        -0.015 |      -1.542 | ok           |
| pure     | 1x    |  0.00 | XY    | reference |    138 |     11 |    1.003 |    1.058 | XY        |      0.003 |     -0.058 |   -0.002 |       -0.010 |         0.006 |       3.763 | x_y_disagree |
| pure     | 1x    |  0.00 | XY    | bins      |    138 |     12 |    0.971 |    0.930 | XY        |     -0.029 |      0.070 |   -0.021 |       -0.031 |        -0.010 |      -5.236 | x_y_disagree |
| pure     | 1x    |  0.00 | XY    | reference |    138 |     12 |    0.994 |    0.971 | XY        |     -0.006 |      0.029 |   -0.003 |       -0.014 |         0.007 |      -1.746 | ok           |
| xx_major | 0.1x  |  0.10 | XX    | bins      |    138 |     11 |    1.854 |    0.104 | XX        |      0.146 |      0.104 |    0.121 |        0.107 |         0.135 |       2.884 | ok           |
| xx_major | 0.1x  |  0.10 | XX    | reference |    138 |     11 |    1.904 |    0.110 | XX        |      0.096 |      0.110 |    0.104 |        0.090 |         0.119 |      -0.912 | ok           |
| xx_major | 0.1x  |  0.20 | XX    | bins      |    138 |     11 |    1.762 |    0.201 | XX        |      0.238 |      0.201 |    0.223 |        0.207 |         0.239 |       2.207 | ok           |
| xx_major | 0.1x  |  0.20 | XX    | reference |    138 |     11 |    1.809 |    0.212 | XX        |      0.191 |      0.212 |    0.199 |        0.183 |         0.216 |      -1.263 | ok           |
| xx_major | 0.1x  |  0.30 | XX    | bins      |    138 |     11 |    1.660 |    0.313 | XX        |      0.340 |      0.313 |    0.332 |        0.315 |         0.349 |       1.456 | ok           |
| xx_major | 0.1x  |  0.30 | XX    | reference |    138 |     11 |    1.705 |    0.330 | XX        |      0.295 |      0.330 |    0.306 |        0.288 |         0.323 |      -1.791 | ok           |
| xx_major | 0.4x  |  0.05 | XX    | bins      |    138 |     11 |    1.891 |    0.050 | XX        |      0.109 |      0.050 |    0.065 |        0.057 |         0.073 |       6.172 | x_y_disagree |
| xx_major | 0.4x  |  0.05 | XX    | reference |    138 |     11 |    1.941 |    0.053 | XX        |      0.059 |      0.053 |    0.054 |        0.046 |         0.063 |       0.563 | ok           |
| xx_major | 0.4x  |  0.10 | XX    | bins      |    138 |     11 |    1.845 |    0.104 | XX        |      0.155 |      0.104 |    0.125 |        0.115 |         0.135 |       4.802 | x_y_disagree |
| xx_major | 0.4x  |  0.10 | XX    | reference |    138 |     11 |    1.894 |    0.110 | XX        |      0.106 |      0.110 |    0.108 |        0.098 |         0.119 |      -0.423 | ok           |
| xx_major | 0.4x  |  0.20 | XX    | bins      |    138 |     11 |    1.750 |    0.202 | XX        |      0.250 |      0.202 |    0.230 |        0.218 |         0.241 |       4.017 | x_y_disagree |
| xx_major | 0.4x  |  0.20 | XX    | reference |    138 |     11 |    1.797 |    0.213 | XX        |      0.203 |      0.213 |    0.207 |        0.195 |         0.219 |      -0.832 | ok           |
| xx_major | 1x    |  0.02 | XX    | bins      |    138 |     11 |    1.922 |    0.020 | XX        |      0.078 |      0.020 |    0.027 |        0.021 |         0.032 |       6.964 | x_y_disagree |
| xx_major | 1x    |  0.02 | XX    | reference |    138 |     11 |    1.974 |    0.021 | XX        |      0.026 |      0.021 |    0.022 |        0.016 |         0.027 |       0.562 | ok           |
| xx_major | 1x    |  0.05 | XX    | bins      |    138 |     11 |    1.894 |    0.051 | XX        |      0.106 |      0.051 |    0.065 |        0.057 |         0.072 |       6.183 | x_y_disagree |
| xx_major | 1x    |  0.05 | XX    | reference |    138 |     11 |    1.945 |    0.054 | XX        |      0.055 |      0.054 |    0.054 |        0.046 |         0.062 |       0.069 | ok           |
| xx_major | 1x    |  0.10 | XX    | bins      |    138 |     11 |    1.847 |    0.099 | XX        |      0.153 |      0.099 |    0.120 |        0.111 |         0.130 |       5.590 | x_y_disagree |
| xx_major | 1x    |  0.10 | XX    | reference |    138 |     11 |    1.897 |    0.105 | XX        |      0.103 |      0.105 |    0.104 |        0.095 |         0.114 |      -0.179 | ok           |
| xx_major | 1x    |  0.20 | XX    | bins      |    138 |     11 |    1.751 |    0.199 | XX        |      0.249 |      0.199 |    0.228 |        0.218 |         0.239 |       4.615 | x_y_disagree |
| xx_major | 1x    |  0.20 | XX    | reference |    138 |     11 |    1.798 |    0.210 | XX        |      0.202 |      0.210 |    0.205 |        0.195 |         0.216 |      -0.754 | ok           |
| xy_major | 0.1x  |  0.10 | XY    | bins      |    138 |     11 |    1.068 |    0.910 | XY        |      0.068 |      0.090 |    0.070 |        0.055 |         0.084 |      -0.840 | ok           |
| xy_major | 0.1x  |  0.10 | XY    | reference |    138 |     11 |    1.097 |    0.960 | XY        |      0.097 |      0.040 |    0.092 |        0.077 |         0.107 |       2.052 | ok           |
| xy_major | 0.1x  |  0.20 | XY    | bins      |    138 |     11 |    1.167 |    0.820 | XY        |      0.167 |      0.180 |    0.168 |        0.153 |         0.183 |      -0.531 | ok           |
| xy_major | 0.1x  |  0.20 | XY    | reference |    138 |     11 |    1.198 |    0.865 | XY        |      0.198 |      0.135 |    0.192 |        0.177 |         0.208 |       2.390 | ok           |
| xy_major | 0.1x  |  0.30 | XY    | bins      |    138 |     11 |    1.263 |    0.716 | XY        |      0.263 |      0.284 |    0.265 |        0.250 |         0.281 |      -0.896 | ok           |
| xy_major | 0.1x  |  0.30 | XY    | reference |    138 |     11 |    1.297 |    0.755 | XY        |      0.297 |      0.245 |    0.290 |        0.274 |         0.306 |       2.071 | ok           |
| xy_major | 0.4x  |  0.05 | XY    | bins      |    138 |     11 |    1.024 |    0.958 | XY        |      0.024 |      0.042 |    0.026 |        0.016 |         0.035 |      -0.992 | ok           |
| xy_major | 0.4x  |  0.05 | XY    | reference |    138 |     11 |    1.052 |    1.011 | XY        |      0.052 |     -0.011 |    0.047 |        0.038 |         0.057 |       3.438 | x_y_disagree |
| xy_major | 0.4x  |  0.10 | XY    | bins      |    138 |     11 |    1.073 |    0.903 | XY        |      0.073 |      0.097 |    0.075 |        0.065 |         0.084 |      -1.418 | ok           |
| xy_major | 0.4x  |  0.10 | XY    | reference |    138 |     11 |    1.102 |    0.953 | XY        |      0.102 |      0.047 |    0.097 |        0.087 |         0.107 |       3.015 | x_y_disagree |
| xy_major | 0.4x  |  0.20 | XY    | bins      |    138 |     11 |    1.175 |    0.804 | XY        |      0.175 |      0.196 |    0.177 |        0.167 |         0.187 |      -1.244 | ok           |
| xy_major | 0.4x  |  0.20 | XY    | reference |    138 |     11 |    1.207 |    0.849 | XY        |      0.207 |      0.151 |    0.201 |        0.191 |         0.212 |       3.172 | x_y_disagree |
| xy_major | 1x    |  0.02 | XY    | bins      |    138 |     11 |    0.995 |    0.982 | XY        |     -0.005 |      0.018 |   -0.003 |       -0.011 |         0.005 |      -1.528 | ok           |
| xy_major | 1x    |  0.02 | XY    | reference |    138 |     11 |    1.022 |    1.036 | XY        |      0.022 |     -0.036 |    0.018 |        0.010 |         0.026 |       3.613 | x_y_disagree |
| xy_major | 1x    |  0.05 | XY    | bins      |    138 |     11 |    1.025 |    0.953 | XY        |      0.025 |      0.047 |    0.027 |        0.019 |         0.035 |      -1.471 | ok           |
| xy_major | 1x    |  0.05 | XY    | reference |    138 |     11 |    1.053 |    1.005 | XY        |      0.053 |     -0.005 |    0.048 |        0.040 |         0.056 |       3.683 | x_y_disagree |
| xy_major | 1x    |  0.10 | XY    | bins      |    138 |     11 |    1.076 |    0.901 | XY        |      0.076 |      0.099 |    0.078 |        0.070 |         0.086 |      -1.583 | ok           |
| xy_major | 1x    |  0.10 | XY    | reference |    138 |     11 |    1.105 |    0.950 | XY        |      0.105 |      0.050 |    0.100 |        0.092 |         0.109 |       3.554 | x_y_disagree |
| xy_major | 1x    |  0.20 | XY    | bins      |    138 |     11 |    1.174 |    0.804 | XY        |      0.174 |      0.196 |    0.176 |        0.168 |         0.185 |      -1.549 | ok           |
| xy_major | 1x    |  0.20 | XY    | reference |    138 |     11 |    1.206 |    0.848 | XY        |      0.206 |      0.152 |    0.200 |        0.191 |         0.209 |       3.532 | x_y_disagree |

### Region fraction

Every autosome of the pure genomes, with the reference, summarised by depth and mode: the mean and the largest absolute shift and the mean standard error.

| depth | mode      | regions | mean_shift | max_abs_shift | mean_shift_se |
|:------|:----------|--------:|-----------:|--------------:|--------------:|
| 0.1x  | bins      |      66 |     -0.002 |         0.063 |         0.013 |
| 0.1x  | reference |      66 |      0.001 |         0.055 |         0.014 |
| 0.4x  | bins      |      66 |     -0.003 |         0.046 |         0.009 |
| 0.4x  | reference |      66 |      0.000 |         0.032 |         0.010 |
| 1x    | bins      |      66 |     -0.003 |         0.042 |         0.008 |
| 1x    | reference |      66 |      0.000 |         0.025 |         0.009 |

chr18 of the planted gains; `truth` is m. `five_chromosomes` is the same composite with the correction and the reference trained on chr18 to chr22, X and Y alone.

| depth | truth | mode             | bins | shift | shift_se | shift_low | shift_high | dispersion | reference_shift | reference_samples | status |
|:------|------:|:-----------------|-----:|------:|---------:|----------:|-----------:|-----------:|----------------:|------------------:|:-------|
| 0.1x  |  0.20 | bins             |   73 | 0.171 |    0.015 |     0.141 |      0.201 |      2.672 |              NA |                NA | ok     |
| 0.1x  |  0.20 | five_chromosomes |   70 | 0.122 |    0.021 |    -0.143 |      0.387 |      3.306 |           0.006 |                 2 | ok     |
| 0.1x  |  0.20 | reference        |   73 | 0.178 |    0.016 |    -0.029 |      0.385 |      2.672 |          -0.007 |                 2 | ok     |
| 0.1x  |  0.40 | bins             |   73 | 0.371 |    0.016 |     0.340 |      0.402 |      2.665 |              NA |                NA | ok     |
| 0.1x  |  0.40 | five_chromosomes |   70 | 0.225 |    0.026 |    -0.105 |      0.555 |      5.007 |           0.006 |                 2 | ok     |
| 0.1x  |  0.40 | reference        |   73 | 0.378 |    0.017 |     0.163 |      0.593 |      2.665 |          -0.007 |                 2 | ok     |
| 0.4x  |  0.10 | bins             |   73 | 0.082 |    0.010 |     0.063 |      0.101 |      4.441 |              NA |                NA | ok     |
| 0.4x  |  0.10 | five_chromosomes |   70 | 0.059 |    0.013 |    -0.104 |      0.222 |      4.774 |           0.006 |                 2 | ok     |
| 0.4x  |  0.10 | reference        |   73 | 0.089 |    0.011 |    -0.053 |      0.231 |      4.441 |          -0.007 |                 2 | ok     |
| 0.4x  |  0.20 | bins             |   73 | 0.183 |    0.010 |     0.164 |      0.202 |      4.430 |              NA |                NA | ok     |
| 0.4x  |  0.20 | five_chromosomes |   70 | 0.121 |    0.014 |    -0.062 |      0.304 |      6.038 |           0.006 |                 2 | ok     |
| 0.4x  |  0.20 | reference        |   73 | 0.190 |    0.011 |     0.046 |      0.335 |      4.430 |          -0.007 |                 2 | ok     |
| 1x    |  0.04 | bins             |   73 | 0.034 |    0.008 |     0.018 |      0.050 |      8.301 |              NA |                NA | ok     |
| 1x    |  0.04 | five_chromosomes |   70 | 0.027 |    0.011 |    -0.115 |      0.170 |      8.948 |           0.006 |                 2 | ok     |
| 1x    |  0.04 | reference        |   73 | 0.041 |    0.010 |    -0.086 |      0.168 |      8.301 |          -0.007 |                 2 | ok     |
| 1x    |  0.10 | bins             |   73 | 0.094 |    0.008 |     0.077 |      0.110 |      8.294 |              NA |                NA | ok     |
| 1x    |  0.10 | five_chromosomes |   70 | 0.066 |    0.012 |    -0.083 |      0.214 |      9.627 |           0.006 |                 2 | ok     |
| 1x    |  0.10 | reference        |   73 | 0.101 |    0.010 |    -0.027 |      0.229 |      8.294 |          -0.007 |                 2 | ok     |
| 1x    |  0.20 | bins             |   73 | 0.193 |    0.009 |     0.176 |      0.210 |      8.244 |              NA |                NA | ok     |
| 1x    |  0.20 | five_chromosomes |   70 | 0.120 |    0.014 |    -0.055 |      0.295 |     13.685 |           0.006 |                 2 | ok     |
| 1x    |  0.20 | reference        |   73 | 0.200 |    0.010 |     0.070 |      0.330 |      8.244 |          -0.007 |                 2 | ok     |

## Scaling

Overheads, live state of the largest input and ceilings (MiB):

| macro           | overhead_mib | live_state_mib | budget_mib |
|:----------------|-------------:|---------------:|-----------:|
| depth_bins      |         45.7 |         4098.6 |    12341.5 |
| region_fraction |         52.6 |         2893.1 |     8732.0 |
| sex_fraction    |         44.2 |           30.0 |      134.2 |

Inputs (rows = samples × bins; the joint input doubles both):

| series  | scale | samples | bin_bp |    bins |       rows |
|:--------|------:|--------:|-------:|--------:|-----------:|
| samples |     1 |      64 | 25,000 | 123,440 |  7,900,160 |
| samples |     2 |     128 | 25,000 | 123,440 | 15,800,320 |
| samples |     4 |     256 | 25,000 | 123,440 | 31,600,640 |
| bins    |     2 |      64 | 12,500 | 246,880 | 15,800,320 |
| bins    |     4 |      64 |  6,250 | 493,760 | 31,600,640 |
| joint   |     4 |     128 | 12,500 | 246,880 | 31,600,640 |

Median, minimum and maximum query latency of three repetitions, peak RSS and DuckDB’s peak buffer memory:

| macro           | series  | scale | threads | samples |    bins | regions |       rows | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | max_buffer_mib | budget_mib |
|:----------------|:--------|------:|--------:|--------:|--------:|--------:|-----------:|---------------:|------------:|------------:|---------------:|------------:|---------------:|-----------:|
| depth_bins      | bins    |     1 |       1 |      64 | 123,440 |      NA |  7,900,160 |          6.108 |       6.102 |       6.162 |        2,230.0 |     2,230.4 |        4,709.1 |   12,341.5 |
| depth_bins      | bins    |     2 |       1 |      64 | 246,880 |      NA | 15,800,320 |         12.240 |      12.240 |      12.316 |        4,406.6 |     4,406.6 |        9,392.7 |   12,341.5 |
| depth_bins      | bins    |     4 |       1 |      64 | 493,760 |      NA | 31,600,640 |         24.736 |      24.587 |      24.736 |        8,758.4 |     8,758.5 |       18,763.5 |   12,341.5 |
| depth_bins      | bins    |     1 |       4 |      64 | 123,440 |      NA |  7,900,160 |          1.831 |       1.816 |       2.002 |        2,242.6 |     2,242.8 |        4,754.6 |   12,341.5 |
| depth_bins      | bins    |     2 |       4 |      64 | 246,880 |      NA | 15,800,320 |          3.526 |       3.510 |       3.565 |        4,419.9 |     4,420.2 |        9,440.7 |   12,341.5 |
| depth_bins      | bins    |     4 |       4 |      64 | 493,760 |      NA | 31,600,640 |          6.907 |       6.860 |       6.909 |        8,769.2 |     8,784.0 |       18,810.4 |   12,341.5 |
| depth_bins      | joint   |     1 |       1 |      64 | 123,440 |      NA |  7,900,160 |          6.108 |       6.102 |       6.162 |        2,230.0 |     2,230.4 |        4,709.1 |   12,341.5 |
| depth_bins      | joint   |     4 |       1 |     128 | 246,880 |      NA | 31,600,640 |         24.787 |      24.723 |      24.902 |        8,758.6 |     8,758.8 |       18,762.0 |   12,341.5 |
| depth_bins      | joint   |     1 |       4 |      64 | 123,440 |      NA |  7,900,160 |          1.831 |       1.816 |       2.002 |        2,242.6 |     2,242.8 |        4,754.6 |   12,341.5 |
| depth_bins      | joint   |     4 |       4 |     128 | 246,880 |      NA | 31,600,640 |          6.949 |       6.938 |       7.002 |        8,773.1 |     8,777.4 |       18,810.1 |   12,341.5 |
| depth_bins      | samples |     1 |       1 |      64 | 123,440 |      NA |  7,900,160 |          6.108 |       6.102 |       6.162 |        2,230.0 |     2,230.4 |        4,709.1 |   12,341.5 |
| depth_bins      | samples |     2 |       1 |     128 | 123,440 |      NA | 15,800,320 |         12.418 |      12.383 |      12.421 |        4,409.2 |     4,409.8 |        9,393.2 |   12,341.5 |
| depth_bins      | samples |     4 |       1 |     256 | 123,440 |      NA | 31,600,640 |         24.995 |      24.827 |      25.142 |        8,758.9 |     8,759.0 |       18,762.6 |   12,341.5 |
| depth_bins      | samples |     1 |       4 |      64 | 123,440 |      NA |  7,900,160 |          1.831 |       1.816 |       2.002 |        2,242.6 |     2,242.8 |        4,754.6 |   12,341.5 |
| depth_bins      | samples |     2 |       4 |     128 | 123,440 |      NA | 15,800,320 |          3.532 |       3.521 |       3.589 |        4,421.0 |     4,423.2 |        9,440.1 |   12,341.5 |
| depth_bins      | samples |     4 |       4 |     256 | 123,440 |      NA | 31,600,640 |          7.391 |       7.283 |       7.516 |        8,775.9 |     8,784.0 |       18,817.0 |   12,341.5 |
| region_fraction | bins    |     1 |       1 |      64 | 123,440 |      22 |  7,900,160 |          1.712 |       1.706 |       1.714 |          635.9 |       636.4 |        1,231.1 |    8,732.0 |
| region_fraction | bins    |     2 |       1 |      64 | 246,880 |      22 | 15,800,320 |          3.394 |       3.394 |       3.404 |        1,219.4 |     1,219.6 |        2,442.7 |    8,732.0 |
| region_fraction | bins    |     4 |       1 |      64 | 493,760 |      22 | 31,600,640 |          6.754 |       6.750 |       6.760 |        2,385.3 |     2,385.3 |        4,867.0 |    8,732.0 |
| region_fraction | bins    |     1 |       4 |      64 | 123,440 |      22 |  7,900,160 |          0.499 |       0.498 |       0.500 |          643.1 |       643.2 |        1,251.6 |    8,732.0 |
| region_fraction | bins    |     2 |       4 |      64 | 246,880 |      22 | 15,800,320 |          0.937 |       0.937 |       1.006 |        1,225.0 |     1,225.8 |        2,463.9 |    8,732.0 |
| region_fraction | bins    |     4 |       4 |      64 | 493,760 |      22 | 31,600,640 |          1.808 |       1.808 |       1.813 |        2,390.4 |     2,390.9 |        4,887.4 |    8,732.0 |
| region_fraction | joint   |     1 |       1 |      64 | 123,440 |      22 |  7,900,160 |          1.712 |       1.706 |       1.714 |          635.9 |       636.4 |        1,231.1 |    8,732.0 |
| region_fraction | joint   |     4 |       1 |     128 | 246,880 |      22 | 31,600,640 |          6.758 |       6.739 |       6.764 |        2,385.2 |     2,385.2 |        4,867.8 |    8,732.0 |
| region_fraction | joint   |     1 |       4 |      64 | 123,440 |      22 |  7,900,160 |          0.499 |       0.498 |       0.500 |          643.1 |       643.2 |        1,251.6 |    8,732.0 |
| region_fraction | joint   |     4 |       4 |     128 | 246,880 |      22 | 31,600,640 |          1.807 |       1.802 |       1.813 |        2,391.0 |     2,391.1 |        4,887.7 |    8,732.0 |
| region_fraction | regions |     1 |       1 |      64 | 123,440 |      22 |  7,900,160 |          1.712 |       1.706 |       1.714 |          635.9 |       636.4 |        1,231.1 |    8,732.0 |
| region_fraction | regions |     2 |       1 |      64 | 123,440 |      44 |  7,900,160 |          1.793 |       1.792 |       1.794 |          637.8 |       638.1 |        1,232.0 |    8,732.0 |
| region_fraction | regions |     4 |       1 |      64 | 123,440 |      88 |  7,900,160 |          1.943 |       1.937 |       1.947 |          638.0 |       638.3 |        1,232.4 |    8,732.0 |
| region_fraction | regions |     1 |       4 |      64 | 123,440 |      22 |  7,900,160 |          0.499 |       0.498 |       0.500 |          643.1 |       643.2 |        1,251.6 |    8,732.0 |
| region_fraction | regions |     2 |       4 |      64 | 123,440 |      44 |  7,900,160 |          0.524 |       0.523 |       0.524 |          643.9 |       644.1 |        1,249.5 |    8,732.0 |
| region_fraction | regions |     4 |       4 |      64 | 123,440 |      88 |  7,900,160 |          0.563 |       0.561 |       0.568 |          645.1 |       645.9 |        1,257.2 |    8,732.0 |
| region_fraction | samples |     1 |       1 |      64 | 123,440 |      22 |  7,900,160 |          1.712 |       1.706 |       1.714 |          635.9 |       636.4 |        1,231.1 |    8,732.0 |
| region_fraction | samples |     2 |       1 |     128 | 123,440 |      22 | 15,800,320 |          3.389 |       3.385 |       3.403 |        1,219.8 |     1,220.0 |        2,443.5 |    8,732.0 |
| region_fraction | samples |     4 |       1 |     256 | 123,440 |      22 | 31,600,640 |          6.740 |       6.730 |       6.747 |        2,385.4 |     2,385.5 |        4,867.9 |    8,732.0 |
| region_fraction | samples |     1 |       4 |      64 | 123,440 |      22 |  7,900,160 |          0.499 |       0.498 |       0.500 |          643.1 |       643.2 |        1,251.6 |    8,732.0 |
| region_fraction | samples |     2 |       4 |     128 | 123,440 |      22 | 15,800,320 |          0.935 |       0.932 |       0.936 |        1,226.1 |     1,226.4 |        2,463.5 |    8,732.0 |
| region_fraction | samples |     4 |       4 |     256 | 123,440 |      22 | 31,600,640 |          1.809 |       1.807 |       1.809 |        2,393.6 |     2,393.8 |        4,888.9 |    8,732.0 |
| sex_fraction    | bins    |     1 |       1 |      64 | 123,440 |      NA |  7,900,160 |          0.640 |       0.639 |       0.641 |           49.2 |        49.3 |           10.7 |      134.2 |
| sex_fraction    | bins    |     2 |       1 |      64 | 246,880 |      NA | 15,800,320 |          1.268 |       1.267 |       1.270 |           52.6 |        52.8 |           15.1 |      134.2 |
| sex_fraction    | bins    |     4 |       1 |      64 | 493,760 |      NA | 31,600,640 |          2.539 |       2.535 |       2.553 |           59.8 |        59.9 |           21.9 |      134.2 |
| sex_fraction    | bins    |     1 |       4 |      64 | 123,440 |      NA |  7,900,160 |          0.203 |       0.202 |       0.210 |           63.2 |        64.1 |           25.6 |      134.2 |
| sex_fraction    | bins    |     2 |       4 |      64 | 246,880 |      NA | 15,800,320 |          0.376 |       0.373 |       0.381 |           65.7 |        65.8 |           29.6 |      134.2 |
| sex_fraction    | bins    |     4 |       4 |      64 | 493,760 |      NA | 31,600,640 |          0.708 |       0.705 |       0.716 |           71.7 |        72.6 |           56.1 |      134.2 |
| sex_fraction    | joint   |     1 |       1 |      64 | 123,440 |      NA |  7,900,160 |          0.640 |       0.639 |       0.641 |           49.2 |        49.3 |           10.7 |      134.2 |
| sex_fraction    | joint   |     4 |       1 |     128 | 246,880 |      NA | 31,600,640 |          2.537 |       2.519 |       2.538 |           59.9 |        60.0 |           21.9 |      134.2 |
| sex_fraction    | joint   |     1 |       4 |      64 | 123,440 |      NA |  7,900,160 |          0.203 |       0.202 |       0.210 |           63.2 |        64.1 |           25.6 |      134.2 |
| sex_fraction    | joint   |     4 |       4 |     128 | 246,880 |      NA | 31,600,640 |          0.706 |       0.701 |       0.716 |           71.8 |        72.0 |           55.6 |      134.2 |
| sex_fraction    | samples |     1 |       1 |      64 | 123,440 |      NA |  7,900,160 |          0.640 |       0.639 |       0.641 |           49.2 |        49.3 |           10.7 |      134.2 |
| sex_fraction    | samples |     2 |       1 |     128 | 123,440 |      NA | 15,800,320 |          1.277 |       1.277 |       1.277 |           52.6 |        52.6 |           15.4 |      134.2 |
| sex_fraction    | samples |     4 |       1 |     256 | 123,440 |      NA | 31,600,640 |          2.538 |       2.530 |       2.544 |           59.8 |        59.9 |           22.4 |      134.2 |
| sex_fraction    | samples |     1 |       4 |      64 | 123,440 |      NA |  7,900,160 |          0.203 |       0.202 |       0.210 |           63.2 |        64.1 |           25.6 |      134.2 |
| sex_fraction    | samples |     2 |       4 |     128 | 123,440 |      NA | 15,800,320 |          0.371 |       0.371 |       0.377 |           66.3 |        67.1 |           30.0 |      134.2 |
| sex_fraction    | samples |     4 |       4 |     256 | 123,440 |      NA | 31,600,640 |          0.710 |       0.706 |       0.710 |           72.1 |        72.1 |           55.8 |      134.2 |

The real workload, every composite of the real-data checks in 1 Mb bins, corrected by `duckhts_depth_bins` (no timing verdict: it is far below 5 s):

| threads |    rows | median_seconds | max_rss_mib |
|--------:|--------:|---------------:|------------:|
|       1 | 101,592 |          0.092 |        86.4 |
|       4 | 101,592 |          0.049 |        92.0 |

Growth of the median time against the series’ dimension (an exponent of 1 is linear; for the joint series, against the rows) and the peak-RSS ratio of the largest to the 1× input:

| macro           | series  | threads | one_x_seconds | time_exponent_2x | time_exponent_4x | rss_ratio_4x | timing_verdict    |
|:----------------|:--------|--------:|--------------:|-----------------:|-----------------:|-------------:|:------------------|
| depth_bins      | bins    |       1 |         6.108 |            1.003 |            1.015 |        3.928 | yes               |
| region_fraction | bins    |       1 |         1.712 |            0.987 |            0.993 |        3.751 | no (1x under 5 s) |
| sex_fraction    | bins    |       1 |         0.640 |            0.986 |            1.001 |        1.215 | no (1x under 5 s) |
| depth_bins      | joint   |       1 |         6.108 |               NA |            1.010 |        3.928 | yes               |
| region_fraction | joint   |       1 |         1.712 |               NA |            0.990 |        3.751 | no (1x under 5 s) |
| sex_fraction    | joint   |       1 |         0.640 |               NA |            0.993 |        1.217 | no (1x under 5 s) |
| region_fraction | regions |       1 |         1.712 |            0.066 |            0.116 |        1.003 | no (1x under 5 s) |
| depth_bins      | samples |       1 |         6.108 |            1.024 |            1.009 |        3.928 | yes               |
| region_fraction | samples |       1 |         1.712 |            0.985 |            0.992 |        3.751 | no (1x under 5 s) |
| sex_fraction    | samples |       1 |         0.640 |            0.996 |            0.991 |        1.215 | no (1x under 5 s) |
| depth_bins      | bins    |       4 |         1.831 |            0.945 |            0.970 |        3.910 | yes               |
| region_fraction | bins    |       4 |         0.499 |            0.910 |            0.948 |        3.717 | no (1x under 5 s) |
| sex_fraction    | bins    |       4 |         0.203 |            0.886 |            0.914 |        1.134 | no (1x under 5 s) |
| depth_bins      | joint   |       4 |         1.831 |               NA |            0.962 |        3.912 | yes               |
| region_fraction | joint   |       4 |         0.499 |               NA |            0.928 |        3.718 | no (1x under 5 s) |
| sex_fraction    | joint   |       4 |         0.203 |               NA |            0.898 |        1.136 | no (1x under 5 s) |
| region_fraction | regions |       4 |         0.499 |            0.070 |            0.104 |        1.003 | no (1x under 5 s) |
| depth_bins      | samples |       4 |         1.831 |            0.948 |            1.065 |        3.913 | yes               |
| region_fraction | samples |       4 |         0.499 |            0.906 |            0.953 |        3.722 | no (1x under 5 s) |
| sex_fraction    | samples |       4 |         0.203 |            0.869 |            0.936 |        1.141 | no (1x under 5 s) |

## Findings

Every observation is within its ceiling and none spills. Timing verdicts apply to the series whose 1x run takes at least 5 s on one thread (depth_bins bins, depth_bins joint, depth_bins samples), and none of their exponents is above 1.07; the other series get a memory verdict and their times without a timing verdict. `duckhts_depth_bins` corrects 7,900,160 rows in 6.1 s on one thread and 31,600,640 rows in 25.0 s (7.4 s on four threads), at 8759 MiB peak RSS, 289 bytes a row against the 136 declared. The region estimator takes 6.7 s and 2385.5 MiB at 31,600,640 rows, and the sex estimator, which streams, 2.5 s and 59.9 MiB.

On the real genomes, without a reference the X of the pure XY genomes reads as 0.967 to 0.975 copies, so the X estimate of the fraction is -0.033 to -0.025, while the reference brings every pure fraction within 0.009 of 0. The residual is a difference in depth between X and the autosomes that GC does not explain (mappability, read filters); a reference set of the same processing removes it.

HG00403 reads as 1.06 to 1.07 Y copies against the other two genomes, so its Y estimate of the fraction is negative and `agreement_z` flags it (`x_y_disagree` at 0.4x and 1x): the Y chromosomes of different men differ in copy number, which a reference of other men does not remove. The combined fraction stays within the acceptance because the X estimate carries most of the weight. The planted chr18 gains are recovered with the reference within 0.022 of m on the whole genome, and at 0.56 to 0.68 of m when the correction and the reference see chr18 to chr22 alone.
