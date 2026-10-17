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
  three inputs are checked against their registered row and read counts and
  receipt SHA-256 values. Staging validates each full CRAM against its supplier
  MD5; receipts retain the producer and actual input hashes. The reference
  FASTA and sample panel are checked against their registered checksums before
  anything is computed. GC comes from
  `fasta_nuc()` on the registered GRCh38 reference.

Every timed observation is one fresh DuckDB CLI process under
`/usr/bin/time`, with its launcher’s actual CPU affinity and scheduling
priority recorded below, with three repetitions at one thread and at four
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

| artifact                       | producer_sha256                                                  | alignments_sha256                                                | output_sha256                                                    |
|:-------------------------------|:-----------------------------------------------------------------|:-----------------------------------------------------------------|:-----------------------------------------------------------------|
| depth_fraction_hg00403_buckets | ef335aaf51f413b09365cecec55417c7b1efe860cd2db58ff5243c997976daa8 | fd85e89918bf47a7378761eb925fd654557b782e9eab896f290af43b3c2c7b15 | 6ee4a25bb1ce08469aae1003761ad7123860a924dc10f44b21d84fb7bb63c55b |
| depth_fraction_na18507_buckets | d8d601db566efafd7ae5e09650f32bd7c930a4d060e715405eb36209e2aaab88 | d4ae010de22e991805b91e95ad4352fe1658970fa577b11e696c94b737313cb8 | 7cf2438544588eaf8a9707b434f2b613de4ae692650c401ca67f79fca6513147 |
| depth_fraction_hg00188_buckets | 17b1a43a2216db7505d1323467900b5638608086664ab44d9aaef5213097c22e | 1ac3e97343325b497929e2e65a21ff96c7168f559cc28c915274e41aae0116ca | 1dbd8f70ad3d7b7e5ce8e72e1d608ac445de51a3509b472e9fbcff3c54819103 |

Source revision: 53ced276f90af2ec4352e6b528aa5cc51a6bc4f5; `src` tree d79a118e8c83e3977b7503f1f90d722019028bb6. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: 7683031d3ad8bc08b9e489011e43efb60bca55b2cd08308b382684ddc56507f6. Host load at start: 0.36, 0.81, 0.69. DuckDB memory limit: 50.0 GiB (the default); temporary directory limit: 0. Launcher affinity: pid 3670909’s current affinity list: 12-19; niceness: 0.

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
| depth_bins      |         46.2 |         4098.6 |    12342.0 |
| region_fraction |         53.8 |         2893.1 |     8733.1 |
| sex_fraction    |         44.1 |           30.0 |      134.1 |

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
| depth_bins      | bins    |     1 |       1 |      64 | 123,440 |      NA |  7,900,160 |          6.204 |       6.198 |       6.279 |        2,230.8 |     2,230.9 |        4,709.1 |   12,342.0 |
| depth_bins      | bins    |     2 |       1 |      64 | 246,880 |      NA | 15,800,320 |         12.460 |      12.459 |      12.503 |        4,407.3 |     4,407.4 |        9,392.7 |   12,342.0 |
| depth_bins      | bins    |     4 |       1 |      64 | 493,760 |      NA | 31,600,640 |         25.002 |      24.999 |      25.359 |        8,758.8 |     8,758.9 |       18,763.5 |   12,342.0 |
| depth_bins      | bins    |     1 |       4 |      64 | 123,440 |      NA |  7,900,160 |          1.856 |       1.853 |       1.856 |        2,242.7 |     2,243.1 |        4,754.7 |   12,342.0 |
| depth_bins      | bins    |     2 |       4 |      64 | 246,880 |      NA | 15,800,320 |          3.588 |       3.565 |       3.594 |        4,427.2 |     4,428.8 |        9,438.2 |   12,342.0 |
| depth_bins      | bins    |     4 |       4 |      64 | 493,760 |      NA | 31,600,640 |          6.975 |       6.945 |       7.041 |        8,773.2 |     8,784.4 |       18,806.6 |   12,342.0 |
| depth_bins      | joint   |     1 |       1 |      64 | 123,440 |      NA |  7,900,160 |          6.204 |       6.198 |       6.279 |        2,230.8 |     2,230.9 |        4,709.1 |   12,342.0 |
| depth_bins      | joint   |     4 |       1 |     128 | 246,880 |      NA | 31,600,640 |         25.289 |      24.981 |      25.378 |        8,759.2 |     8,759.3 |       18,762.0 |   12,342.0 |
| depth_bins      | joint   |     1 |       4 |      64 | 123,440 |      NA |  7,900,160 |          1.856 |       1.853 |       1.856 |        2,242.7 |     2,243.1 |        4,754.7 |   12,342.0 |
| depth_bins      | joint   |     4 |       4 |     128 | 246,880 |      NA | 31,600,640 |          7.089 |       7.030 |       7.106 |        8,787.8 |     8,796.1 |       18,806.9 |   12,342.0 |
| depth_bins      | samples |     1 |       1 |      64 | 123,440 |      NA |  7,900,160 |          6.204 |       6.198 |       6.279 |        2,230.8 |     2,230.9 |        4,709.1 |   12,342.0 |
| depth_bins      | samples |     2 |       1 |     128 | 123,440 |      NA | 15,800,320 |         12.637 |      12.564 |      12.643 |        4,409.3 |     4,409.5 |        9,393.2 |   12,342.0 |
| depth_bins      | samples |     4 |       1 |     256 | 123,440 |      NA | 31,600,640 |         25.415 |      25.394 |      25.591 |        8,759.4 |     8,759.5 |       18,762.6 |   12,342.0 |
| depth_bins      | samples |     1 |       4 |      64 | 123,440 |      NA |  7,900,160 |          1.856 |       1.853 |       1.856 |        2,242.7 |     2,243.1 |        4,754.7 |   12,342.0 |
| depth_bins      | samples |     2 |       4 |     128 | 123,440 |      NA | 15,800,320 |          3.596 |       3.585 |       3.605 |        4,423.6 |     4,425.9 |        9,439.5 |   12,342.0 |
| depth_bins      | samples |     4 |       4 |     256 | 123,440 |      NA | 31,600,640 |          7.437 |       7.427 |       7.544 |        8,774.1 |     8,776.4 |       18,810.3 |   12,342.0 |
| region_fraction | bins    |     1 |       1 |      64 | 123,440 |      22 |  7,900,160 |          1.689 |       1.685 |       1.693 |          637.0 |       637.0 |        1,231.1 |    8,733.1 |
| region_fraction | bins    |     2 |       1 |      64 | 246,880 |      22 | 15,800,320 |          3.366 |       3.357 |       3.366 |        1,220.1 |     1,220.1 |        2,442.7 |    8,733.1 |
| region_fraction | bins    |     4 |       1 |      64 | 493,760 |      22 | 31,600,640 |          6.683 |       6.683 |       6.685 |        2,385.8 |     2,385.9 |        4,867.0 |    8,733.1 |
| region_fraction | bins    |     1 |       4 |      64 | 123,440 |      22 |  7,900,160 |          0.493 |       0.493 |       0.493 |          642.9 |       644.2 |        1,251.8 |    8,733.1 |
| region_fraction | bins    |     2 |       4 |      64 | 246,880 |      22 | 15,800,320 |          0.923 |       0.922 |       0.928 |        1,225.4 |     1,225.6 |        2,463.9 |    8,733.1 |
| region_fraction | bins    |     4 |       4 |      64 | 493,760 |      22 | 31,600,640 |          1.789 |       1.785 |       1.790 |        2,390.8 |     2,391.1 |        4,886.6 |    8,733.1 |
| region_fraction | joint   |     1 |       1 |      64 | 123,440 |      22 |  7,900,160 |          1.689 |       1.685 |       1.693 |          637.0 |       637.0 |        1,231.1 |    8,733.1 |
| region_fraction | joint   |     4 |       1 |     128 | 246,880 |      22 | 31,600,640 |          6.685 |       6.671 |       6.842 |        2,384.9 |     2,385.4 |        4,867.8 |    8,733.1 |
| region_fraction | joint   |     1 |       4 |      64 | 123,440 |      22 |  7,900,160 |          0.493 |       0.493 |       0.493 |          642.9 |       644.2 |        1,251.8 |    8,733.1 |
| region_fraction | joint   |     4 |       4 |     128 | 246,880 |      22 | 31,600,640 |          1.786 |       1.782 |       1.787 |        2,391.5 |     2,392.4 |        4,887.3 |    8,733.1 |
| region_fraction | regions |     1 |       1 |      64 | 123,440 |      22 |  7,900,160 |          1.689 |       1.685 |       1.693 |          637.0 |       637.0 |        1,231.1 |    8,733.1 |
| region_fraction | regions |     2 |       1 |      64 | 123,440 |      44 |  7,900,160 |          1.725 |       1.724 |       1.740 |          638.1 |       638.2 |        1,232.0 |    8,733.1 |
| region_fraction | regions |     4 |       1 |      64 | 123,440 |      88 |  7,900,160 |          1.789 |       1.789 |       1.792 |          637.7 |       637.8 |        1,232.4 |    8,733.1 |
| region_fraction | regions |     1 |       4 |      64 | 123,440 |      22 |  7,900,160 |          0.493 |       0.493 |       0.493 |          642.9 |       644.2 |        1,251.8 |    8,733.1 |
| region_fraction | regions |     2 |       4 |      64 | 123,440 |      44 |  7,900,160 |          0.507 |       0.504 |       0.507 |          644.7 |       644.9 |        1,250.8 |    8,733.1 |
| region_fraction | regions |     4 |       4 |      64 | 123,440 |      88 |  7,900,160 |          0.524 |       0.523 |       0.524 |          645.4 |       645.9 |        1,252.7 |    8,733.1 |
| region_fraction | samples |     1 |       1 |      64 | 123,440 |      22 |  7,900,160 |          1.689 |       1.685 |       1.693 |          637.0 |       637.0 |        1,231.1 |    8,733.1 |
| region_fraction | samples |     2 |       1 |     128 | 123,440 |      22 | 15,800,320 |          3.359 |       3.356 |       3.364 |        1,220.2 |     1,220.3 |        2,443.5 |    8,733.1 |
| region_fraction | samples |     4 |       1 |     256 | 123,440 |      22 | 31,600,640 |          6.685 |       6.682 |       6.707 |        2,385.0 |     2,385.9 |        4,867.9 |    8,733.1 |
| region_fraction | samples |     1 |       4 |      64 | 123,440 |      22 |  7,900,160 |          0.493 |       0.493 |       0.493 |          642.9 |       644.2 |        1,251.8 |    8,733.1 |
| region_fraction | samples |     2 |       4 |     128 | 123,440 |      22 | 15,800,320 |          0.926 |       0.925 |       0.926 |        1,226.5 |     1,227.5 |        2,458.3 |    8,733.1 |
| region_fraction | samples |     4 |       4 |     256 | 123,440 |      22 | 31,600,640 |          1.787 |       1.784 |       1.789 |        2,393.7 |     2,393.8 |        4,888.5 |    8,733.1 |
| sex_fraction    | bins    |     1 |       1 |      64 | 123,440 |      NA |  7,900,160 |          0.642 |       0.641 |       0.645 |           49.6 |        49.6 |           10.7 |      134.1 |
| sex_fraction    | bins    |     2 |       1 |      64 | 246,880 |      NA | 15,800,320 |          1.276 |       1.272 |       1.276 |           52.9 |        52.9 |           15.1 |      134.1 |
| sex_fraction    | bins    |     4 |       1 |      64 | 493,760 |      NA | 31,600,640 |          2.536 |       2.530 |       2.540 |           60.2 |        60.4 |           21.9 |      134.1 |
| sex_fraction    | bins    |     1 |       4 |      64 | 123,440 |      NA |  7,900,160 |          0.207 |       0.206 |       0.213 |           64.2 |        64.4 |           25.4 |      134.1 |
| sex_fraction    | bins    |     2 |       4 |      64 | 246,880 |      NA | 15,800,320 |          0.374 |       0.372 |       0.377 |           66.6 |        66.8 |           29.6 |      134.1 |
| sex_fraction    | bins    |     4 |       4 |      64 | 493,760 |      NA | 31,600,640 |          0.716 |       0.713 |       0.736 |           72.7 |        72.8 |           55.9 |      134.1 |
| sex_fraction    | joint   |     1 |       1 |      64 | 123,440 |      NA |  7,900,160 |          0.642 |       0.641 |       0.645 |           49.6 |        49.6 |           10.7 |      134.1 |
| sex_fraction    | joint   |     4 |       1 |     128 | 246,880 |      NA | 31,600,640 |          2.532 |       2.527 |       2.539 |           60.1 |        60.1 |           21.9 |      134.1 |
| sex_fraction    | joint   |     1 |       4 |      64 | 123,440 |      NA |  7,900,160 |          0.207 |       0.206 |       0.213 |           64.2 |        64.4 |           25.4 |      134.1 |
| sex_fraction    | joint   |     4 |       4 |     128 | 246,880 |      NA | 31,600,640 |          0.710 |       0.708 |       0.715 |           72.4 |        72.6 |           56.2 |      134.1 |
| sex_fraction    | samples |     1 |       1 |      64 | 123,440 |      NA |  7,900,160 |          0.642 |       0.641 |       0.645 |           49.6 |        49.6 |           10.7 |      134.1 |
| sex_fraction    | samples |     2 |       1 |     128 | 123,440 |      NA | 15,800,320 |          1.274 |       1.271 |       1.277 |           53.0 |        53.0 |           15.4 |      134.1 |
| sex_fraction    | samples |     4 |       1 |     256 | 123,440 |      NA | 31,600,640 |          2.534 |       2.527 |       2.551 |           60.1 |        60.2 |           22.4 |      134.1 |
| sex_fraction    | samples |     1 |       4 |      64 | 123,440 |      NA |  7,900,160 |          0.207 |       0.206 |       0.213 |           64.2 |        64.4 |           25.4 |      134.1 |
| sex_fraction    | samples |     2 |       4 |     128 | 123,440 |      NA | 15,800,320 |          0.377 |       0.369 |       0.379 |           66.3 |        67.1 |           29.8 |      134.1 |
| sex_fraction    | samples |     4 |       4 |     256 | 123,440 |      NA | 31,600,640 |          0.716 |       0.710 |       0.716 |           72.6 |        73.4 |           56.0 |      134.1 |

The real workload, every composite of the real-data checks in 1 Mb bins, corrected by `duckhts_depth_bins` (no timing verdict: it is far below 5 s):

| threads |    rows | median_seconds | max_rss_mib |
|--------:|--------:|---------------:|------------:|
|       1 | 101,592 |          0.093 |        86.6 |
|       4 | 101,592 |          0.050 |        93.0 |

Growth of the median time against the series’ dimension (an exponent of 1 is linear; for the joint series, against the rows) and the peak-RSS ratio of the largest to the 1× input:

| macro           | series  | threads | one_x_seconds | time_exponent_2x | time_exponent_4x | rss_ratio_4x | timing_verdict    |
|:----------------|:--------|--------:|--------------:|-----------------:|-----------------:|-------------:|:------------------|
| depth_bins      | bins    |       1 |         6.204 |            1.006 |            1.005 |        3.926 | yes               |
| region_fraction | bins    |       1 |         1.689 |            0.995 |            0.990 |        3.745 | no (1x under 5 s) |
| sex_fraction    | bins    |       1 |         0.642 |            0.991 |            0.992 |        1.214 | no (1x under 5 s) |
| depth_bins      | joint   |       1 |         6.204 |               NA |            1.014 |        3.926 | yes               |
| region_fraction | joint   |       1 |         1.689 |               NA |            0.992 |        3.744 | no (1x under 5 s) |
| sex_fraction    | joint   |       1 |         0.642 |               NA |            0.990 |        1.212 | no (1x under 5 s) |
| region_fraction | regions |       1 |         1.689 |            0.030 |            0.053 |        1.001 | no (1x under 5 s) |
| depth_bins      | samples |       1 |         6.204 |            1.026 |            1.008 |        3.927 | yes               |
| region_fraction | samples |       1 |         1.689 |            0.991 |            0.993 |        3.744 | no (1x under 5 s) |
| sex_fraction    | samples |       1 |         0.642 |            0.989 |            0.992 |        1.212 | no (1x under 5 s) |
| depth_bins      | bins    |       4 |         1.856 |            0.951 |            0.959 |        3.912 | yes               |
| region_fraction | bins    |       4 |         0.493 |            0.904 |            0.954 |        3.719 | no (1x under 5 s) |
| sex_fraction    | bins    |       4 |         0.207 |            0.855 |            0.938 |        1.132 | no (1x under 5 s) |
| depth_bins      | joint   |       4 |         1.856 |               NA |            0.967 |        3.918 | yes               |
| region_fraction | joint   |       4 |         0.493 |               NA |            0.928 |        3.720 | no (1x under 5 s) |
| sex_fraction    | joint   |       4 |         0.207 |               NA |            0.890 |        1.128 | no (1x under 5 s) |
| region_fraction | regions |       4 |         0.493 |            0.038 |            0.049 |        1.004 | no (1x under 5 s) |
| depth_bins      | samples |       4 |         1.856 |            0.954 |            1.048 |        3.912 | yes               |
| region_fraction | samples |       4 |         0.493 |            0.908 |            0.949 |        3.723 | no (1x under 5 s) |
| sex_fraction    | samples |       4 |         0.207 |            0.868 |            0.923 |        1.131 | no (1x under 5 s) |

## Findings

Every observation is within its ceiling and none spills. Timing verdicts apply to the series whose 1x run takes at least 5 s on one thread (depth_bins bins, depth_bins joint, depth_bins samples), and none of their exponents is above 1.05; the other series get a memory verdict and their times without a timing verdict. `duckhts_depth_bins` corrects 7,900,160 rows in 6.2 s on one thread and 31,600,640 rows in 25.4 s (7.4 s on four threads), at 8759.5 MiB peak RSS, 289 bytes a row against the 136 declared. The region estimator takes 6.7 s and 2385.9 MiB at 31,600,640 rows, and the sex estimator, which streams, 2.5 s and 60.2 MiB.

On the real genomes, without a reference the X of the pure XY genomes reads as 0.967 to 0.975 copies, so the X estimate of the fraction is -0.033 to -0.025, while the reference brings every pure fraction within 0.009 of 0. The residual is a difference in depth between X and the autosomes that GC does not explain (mappability, read filters); a reference set of the same processing removes it.

HG00403 reads as 1.06 to 1.07 Y copies against the other two genomes, so its Y estimate of the fraction is negative and `agreement_z` flags it (`x_y_disagree` at 0.4x and 1x): the Y chromosomes of different men differ in copy number, which a reference of other men does not remove. The combined fraction stays within the acceptance because the X estimate carries most of the weight. The planted chr18 gains are recovered with the reference within 0.022 of m on the whole genome, and at 0.56 to 0.68 of m when the correction and the reference see chr18 to chr22 alone.
