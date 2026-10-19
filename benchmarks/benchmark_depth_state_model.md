Depth state model: agreement with ichorCNA v0.6.1, planted gains and scaling
================

This report validates phase 2 of the depth fraction lane
(`design/depth_fraction_lane.md`, <https://github.com/RGenomicsETL/duckhts/issues/343>):
`duckhts_depth_state_solutions()` and `duckhts_depth_state_segments()`, a port of the
ichorCNA copy-number state model. It runs the port and the unmodified upstream R code
(GavinHaLab/ichorCNA tag v0.6.1 with HMMcopy 1.54.0, `benchmarks/depth_state/upstream_ichorcna.R`)
on the same corrected bins and compares fraction, ploidy, guard and segments; checks planted
chromosome-18 gains and samples without change at 0.1×, 0.4× and 1×; and gives the 1×/2×/4×
scaling evidence `STYLE.md` asks for. No earlier report measures these macros; this is their
first baseline.

The real bins are the composites of `benchmarks/benchmark_depth_fraction.md`: sums of
read-name buckets of three registered 30× genomes, corrected with `duckhts_depth_bins()` and
masked by the leave-one-out panel rule of that report. Upstream cannot take corrected bins
through its command line, which applies its own loess correction, so its model is entered at
the function level that `run_ichorCNA` calls after correction: `getDefaultParameters`,
`HMMsegment` and `correctIntegerCN` with `run_ichorCNA`’s default arguments, plus its guard.

``` sh
taskset -c 12-19 Rscript -e 'rmarkdown::render("benchmarks/benchmark_depth_state_model.Rmd")'
```

Extension at revision `ab28ee89025896d4a4e176d6acac0448f10d48dd` (source tree `ca89acc14a96abbe16c8b5971950485b069950e6`), DuckDB CLI `v1.5.1 (Variegata) 7dbb2e646f`,
host load `0.49, 0.98, 1.27` at start, launcher `pid 1393077's current affinity list: 12-19`.

## Declaration, before any run

Composites, from the design of the depth fraction report: for each of the genomes
HG00403, NA18507 and HG00188, a pure sample at 0.1×, 0.4× and 1× (10, 40 and 100 of the
300 read-name buckets); and HG00403 with extra chromosome-18 buckets so that chromosome 18
has the depth of a one-copy gain in a share `m` of cells, `m` = 0.2 and 0.4 at 0.1×, 0.1 and
0.2 at 0.4×, 0.04, 0.1 and 0.2 at 1×. The gain series is the dilution series: the same
change at decreasing share. Each composite is corrected alone with `duckhts_depth_bins()`
(1 Mb bins, GC from the GRCh38 reference, bins with 5% or more N and the X
pseudo-autosomal regions removed) and keeps the bins that the other two genomes, at
30×, correct to a ratio within 0.75 and 1.25 (the panel mask of the depth fraction report).
Y is dropped; X is decoded and not trained on. Both models run the default grid, normal
fractions 0.5 to 0.9 at ploidy starts 2 and 3, with `max_cn` 7.

Criteria:

- **A1, best solution.** For every composite, the port and upstream rank the same start first,
  or their best solutions agree within the A2 tolerances.
- **A2, parameters.** For the best solution of every composite, the normal fractions differ by at
  most 0.02 and the ploidies by at most 0.05; the guard decision is the same.
- **A3, segments.** For the best solution of every composite, the Viterbi states agree on at
  least 98% of bins.
- **A4, every start.** Over all composites and starts, at least 90% of the starts have normal
  fractions within 0.02, ploidies within 0.05 and objectives within 5 of upstream. The rest are
  listed; a difference on a flat objective is the expected effect of the two optimizers.
- **G1, gains.** For each gain composite the port’s best solution has a ploidy within 0.15 of 2,
  chromosome 18 at three copies over at least 90% of its kept bins, no other autosome with more
  than 10% of its bins altered, and a tumour fraction within 0.05 of `m` at 0.4× and 1× and within
  0.08 at 0.1×.
- **N1, no change.** For each pure composite, the best solution of the port and of upstream has a
  tumour fraction of at most 0.02.
- **S1, scaling.** On synthetic bins in fresh processes with four threads: with one sample, the
  latency at 2× and 4× bins is at most 2.6× and 6.8× that at 1× (exponent at most 1.4); with
  four samples of the 1× size, at most 2× the one-sample latency; peak RSS stays under 1 GiB
  throughout. The model’s declared memory is about `max_cn * 28 + 40` bytes a bin per sample
  in flight.

| test_id             | kind | depth | truth | genomes |
|:--------------------|:-----|:------|------:|:--------|
| pure HG00403 0.1x   | pure | 0.1x  |  0.00 | HG00403 |
| pure NA18507 0.1x   | pure | 0.1x  |  0.00 | NA18507 |
| pure HG00188 0.1x   | pure | 0.1x  |  0.00 | HG00188 |
| gain chr18 0.2 0.1x | gain | 0.1x  |  0.20 | HG00403 |
| gain chr18 0.4 0.1x | gain | 0.1x  |  0.40 | HG00403 |
| pure HG00403 0.4x   | pure | 0.4x  |  0.00 | HG00403 |
| pure NA18507 0.4x   | pure | 0.4x  |  0.00 | NA18507 |
| pure HG00188 0.4x   | pure | 0.4x  |  0.00 | HG00188 |
| gain chr18 0.1 0.4x | gain | 0.4x  |  0.10 | HG00403 |
| gain chr18 0.2 0.4x | gain | 0.4x  |  0.20 | HG00403 |
| pure HG00403 1x     | pure | 1x    |  0.00 | HG00403 |
| pure NA18507 1x     | pure | 1x    |  0.00 | NA18507 |
| pure HG00188 1x     | pure | 1x    |  0.00 | HG00188 |
| gain chr18 0.04 1x  | gain | 1x    |  0.04 | HG00403 |
| gain chr18 0.1 1x   | gain | 1x    |  0.10 | HG00403 |
| gain chr18 0.2 1x   | gain | 1x    |  0.20 | HG00403 |

The composites

## Execution and provenance

| artifact                       | producer_sha256                                                  | output_sha256                                                    |
|:-------------------------------|:-----------------------------------------------------------------|:-----------------------------------------------------------------|
| depth_fraction_hg00403_buckets | ef335aaf51f413b09365cecec55417c7b1efe860cd2db58ff5243c997976daa8 | 6ee4a25bb1ce08469aae1003761ad7123860a924dc10f44b21d84fb7bb63c55b |
| depth_fraction_na18507_buckets | d8d601db566efafd7ae5e09650f32bd7c930a4d060e715405eb36209e2aaab88 | 7cf2438544588eaf8a9707b434f2b613de4ae692650c401ca67f79fca6513147 |
| depth_fraction_hg00188_buckets | 17b1a43a2216db7505d1323467900b5638608086664ab44d9aaef5213097c22e | 1dbd8f70ad3d7b7e5ce8e72e1d608ac445de51a3509b472e9fbcff3c54819103 |

Kept bins per composite after the panel mask

| sample_id           | kept | bins |
|:--------------------|-----:|-----:|
| gain chr18 0.04 1x  | 2686 | 2702 |
| gain chr18 0.1 0.4x | 2686 | 2702 |
| gain chr18 0.1 1x   | 2686 | 2702 |
| gain chr18 0.2 0.1x | 2702 | 2702 |
| gain chr18 0.2 0.4x | 2702 | 2702 |
| gain chr18 0.2 1x   | 2702 | 2702 |
| gain chr18 0.4 0.1x | 2702 | 2702 |
| pure HG00188 0.1x   | 2699 | 2709 |
| pure HG00188 0.4x   | 2709 | 2709 |
| pure HG00188 1x     | 2683 | 2709 |
| pure HG00403 0.1x   | 2686 | 2702 |
| pure HG00403 0.4x   | 2686 | 2702 |
| pure HG00403 1x     | 2686 | 2702 |
| pure NA18507 0.1x   | 2702 | 2712 |
| pure NA18507 0.4x   | 2686 | 2712 |
| pure NA18507 1x     | 2686 | 2712 |

## Agreement with upstream

| test_id             | same_start | normal_duckhts | normal_upstream | ploidy_duckhts | ploidy_upstream | guard_duckhts | guard_upstream | state_agreement | A1    | A2    | A3    |
|:--------------------|:-----------|---------------:|----------------:|---------------:|----------------:|:--------------|:---------------|----------------:|:------|:------|:------|
| pure HG00403 0.1x   | FALSE      |         0.9898 |          0.9898 |         3.2226 |          3.1951 | FALSE         | FALSE          |          0.8850 | TRUE  | TRUE  | FALSE |
| pure NA18507 0.1x   | TRUE       |         0.9926 |          0.9926 |         2.7804 |          2.7804 | FALSE         | FALSE          |          1.0000 | TRUE  | TRUE  | TRUE  |
| pure HG00188 0.1x   | TRUE       |         0.9949 |          0.9949 |         3.3457 |          3.3456 | FALSE         | FALSE          |          1.0000 | TRUE  | TRUE  | TRUE  |
| gain chr18 0.2 0.1x | TRUE       |         0.9283 |          0.9279 |         2.0058 |          2.0048 | FALSE         | FALSE          |          1.0000 | TRUE  | TRUE  | TRUE  |
| gain chr18 0.4 0.1x | FALSE      |         0.8748 |          0.6244 |         2.0311 |          2.0103 | FALSE         | FALSE          |          0.9219 | FALSE | FALSE | FALSE |
| pure HG00403 0.4x   | FALSE      |         0.9997 |          1.0000 |         3.9961 |          4.0009 | FALSE         | FALSE          |          0.8805 | TRUE  | TRUE  | FALSE |
| pure NA18507 0.4x   | FALSE      |         0.9929 |          0.9919 |         3.2499 |          2.9381 | FALSE         | FALSE          |          0.7617 | FALSE | FALSE | FALSE |
| pure HG00188 0.4x   | TRUE       |         0.9814 |          0.9809 |         1.9671 |          1.9612 | FALSE         | FALSE          |          1.0000 | TRUE  | TRUE  | TRUE  |
| gain chr18 0.1 0.4x | FALSE      |         0.9821 |          0.9662 |         2.2880 |          1.9930 | FALSE         | FALSE          |          0.8637 | FALSE | FALSE | FALSE |
| gain chr18 0.2 0.4x | TRUE       |         0.9574 |          0.9574 |         2.0018 |          2.0017 | FALSE         | FALSE          |          1.0000 | TRUE  | TRUE  | TRUE  |
| pure HG00403 1x     | TRUE       |         0.9825 |          0.9817 |         1.9270 |          1.9355 | FALSE         | FALSE          |          0.9792 | TRUE  | TRUE  | FALSE |
| pure NA18507 1x     | FALSE      |         0.9922 |          0.9751 |         3.2029 |          2.1840 | FALSE         | FALSE          |          0.4423 | FALSE | FALSE | FALSE |
| pure HG00188 1x     | TRUE       |         1.0000 |          1.0000 |         1.9576 |          1.9575 | TRUE          | TRUE           |          1.0000 | TRUE  | TRUE  | TRUE  |
| gain chr18 0.04 1x  | FALSE      |         0.9821 |          0.9779 |         1.9886 |          1.9899 | FALSE         | FALSE          |          0.9181 | TRUE  | TRUE  | FALSE |
| gain chr18 0.1 1x   | FALSE      |         0.9141 |          0.9141 |         2.0160 |          2.0160 | FALSE         | FALSE          |          1.0000 | TRUE  | TRUE  | TRUE  |
| gain chr18 0.2 1x   | TRUE       |         0.9071 |          0.9071 |         2.0120 |          2.0120 | FALSE         | FALSE          |          1.0000 | TRUE  | TRUE  | TRUE  |

Best solutions: the port against upstream

| criterion | composites |  of | starts_within | pass  |
|:----------|-----------:|----:|:--------------|:------|
| A1        |         12 |  16 | NA            | FALSE |
| A2        |         12 |  16 | NA            | FALSE |
| A3        |          8 |  16 | NA            | FALSE |
| A4        |         NA |  NA | 88.1% of 160  | FALSE |

Starts outside the A4 tolerances

|     | sample_id           | start_normal | start_ploidy | normal_fraction | normal_fraction_upstream | ploidy | ploidy_upstream |    dlik |
|:----|:--------------------|-------------:|-------------:|----------------:|-------------------------:|-------:|----------------:|--------:|
| 2   | gain chr18 0.04 1x  |          0.5 |            3 |          0.3509 |                   0.5521 | 2.9989 |          2.9985 |  0.0001 |
| 7   | gain chr18 0.04 1x  |          0.8 |            2 |          0.9767 |                   0.9779 | 1.9752 |          1.9899 | 15.2263 |
| 11  | gain chr18 0.1 0.4x |          0.5 |            2 |          0.9821 |                   1.0000 | 2.2880 |          4.0005 | 36.8186 |
| 13  | gain chr18 0.1 0.4x |          0.6 |            2 |          0.9821 |                   1.0000 | 2.2890 |          4.0005 | 36.8104 |
| 14  | gain chr18 0.1 0.4x |          0.6 |            3 |          0.7301 |                   0.6074 | 3.0005 |          3.0003 |  0.0058 |
| 15  | gain chr18 0.1 0.4x |          0.7 |            2 |          0.9821 |                   1.0000 | 2.2940 |          4.0014 | 33.4478 |
| 22  | gain chr18 0.1 1x   |          0.5 |            3 |          0.3903 |                   0.5725 | 2.9992 |          2.9990 |  0.0001 |
| 76  | pure HG00188 0.1x   |          0.7 |            3 |          0.2624 |                   0.0621 | 3.0023 |          3.0020 |  0.0001 |
| 78  | pure HG00188 0.1x   |          0.8 |            3 |          0.4584 |                   0.4312 | 3.0030 |          3.0028 |  0.0001 |
| 84  | pure HG00188 0.4x   |          0.6 |            3 |          0.7315 |                   0.6607 | 2.9937 |          2.9948 |  0.0003 |
| 94  | pure HG00188 1x     |          0.6 |            3 |          0.7378 |                   0.6635 | 2.9945 |          2.9956 |  0.0003 |
| 100 | pure HG00188 1x     |          0.9 |            3 |          0.9223 |                   0.9449 | 2.9829 |          2.9741 |  0.2182 |
| 104 | pure HG00403 0.1x   |          0.6 |            3 |          0.8013 |                   0.8231 | 2.9946 |          2.9940 |  0.0004 |
| 122 | pure HG00403 1x     |          0.5 |            3 |          0.3253 |                   0.5383 | 2.9985 |          2.9980 |  0.0002 |
| 142 | pure NA18507 0.4x   |          0.5 |            3 |          0.5807 |                   0.5021 | 3.0027 |          3.0023 |  0.0008 |
| 147 | pure NA18507 0.4x   |          0.8 |            2 |          0.9912 |                   0.9903 | 2.8243 |          2.7081 |  0.9255 |
| 149 | pure NA18507 0.4x   |          0.9 |            2 |          0.9929 |                   0.9906 | 3.2499 |          2.7354 |  2.9208 |
| 153 | pure NA18507 1x     |          0.6 |            2 |          0.9922 |                   0.9848 | 3.2029 |          2.1754 |  9.3371 |
| 154 | pure NA18507 1x     |          0.6 |            3 |          0.8515 |                   0.8766 | 2.9976 |          2.9972 |  0.0430 |

Observation: agreement among the composites with a reported fraction, and which optimizer ended higher where the starts differ

| composites_with_a_fraction |  A1 |  A2 |  A3 | starts_outside_A4 | port_objective_higher_by_1 | upstream_objective_higher_by_1 |
|---------------------------:|----:|----:|----:|------------------:|---------------------------:|-------------------------------:|
|                          8 |   5 |   5 |   4 |                19 |                          5 |                              1 |

## Planted gains and samples without change

| test_id             | depth | truth | designed_log2 | chr18_median_log2 | tumour_fraction | tumour_fraction_upstream | ploidy | status | chr18_copies | chr18_copies_upstream | extra_copies_times_fraction | chr18_three_copies | worst_other_altered | G1    |
|:--------------------|:------|------:|--------------:|------------------:|----------------:|-------------------------:|-------:|:-------|-------------:|----------------------:|----------------------------:|-------------------:|--------------------:|:------|
| gain chr18 0.2 0.1x | 0.1x  |  0.20 |        0.1375 |            0.1076 |          0.0717 |                   0.0721 | 2.0058 | ok     |            4 |                     4 |                      0.1434 |                  0 |                   1 | FALSE |
| gain chr18 0.4 0.1x | 0.1x  |  0.40 |        0.2630 |            0.2514 |          0.1252 |                   0.3756 | 2.0311 | ok     |            5 |                     3 |                      0.3756 |                  0 |                   0 | FALSE |
| gain chr18 0.1 0.4x | 0.4x  |  0.10 |        0.0704 |            0.0538 |          0.0179 |                   0.0338 | 2.2880 | ok     |            6 |                     4 |                      0.0717 |                  0 |                   1 | FALSE |
| gain chr18 0.2 0.4x | 0.4x  |  0.20 |        0.1375 |            0.1197 |          0.0426 |                   0.0426 | 2.0018 | ok     |            6 |                     6 |                      0.1703 |                  0 |                   1 | FALSE |
| gain chr18 0.04 1x  | 1x    |  0.04 |        0.0286 |            0.0162 |          0.0179 |                   0.0221 | 1.9886 | ok     |            3 |                     3 |                      0.0179 |                  1 |                   1 | FALSE |
| gain chr18 0.1 1x   | 1x    |  0.10 |        0.0704 |            0.0580 |          0.0859 |                   0.0859 | 2.0160 | ok     |            3 |                     3 |                      0.0859 |                  1 |                   0 | TRUE  |
| gain chr18 0.2 1x   | 1x    |  0.20 |        0.1375 |            0.1266 |          0.0929 |                   0.0929 | 2.0120 | ok     |            4 |                     4 |                      0.1858 |                  0 |                   0 | FALSE |

Planted chromosome-18 gains: the best solution of the port, with the copy number both models give chromosome 18

| test_id           | tumour_fraction | status     | tumour_fraction_upstream | guard_upstream | N1    |
|:------------------|----------------:|:-----------|-------------------------:|:---------------|:------|
| pure HG00403 0.1x |          0.0102 | ok         |                   0.0102 | FALSE          | TRUE  |
| pure NA18507 0.1x |          0.0074 | ok         |                   0.0074 | FALSE          | TRUE  |
| pure HG00188 0.1x |          0.0051 | ok         |                   0.0051 | FALSE          | TRUE  |
| pure HG00403 0.4x |          0.0003 | ok         |                   0.0000 | FALSE          | TRUE  |
| pure NA18507 0.4x |          0.0071 | ok         |                   0.0081 | FALSE          | TRUE  |
| pure HG00188 0.4x |          0.0186 | ok         |                   0.0191 | FALSE          | TRUE  |
| pure HG00403 1x   |          0.0175 | ok         |                   0.0183 | FALSE          | TRUE  |
| pure NA18507 1x   |          0.0078 | ok         |                   0.0249 | FALSE          | FALSE |
| pure HG00188 1x   |          0.0000 | guard_zero |                   0.0000 | TRUE           | TRUE  |

Samples without change

## Scaling

| label     | bins_per_contig | samples |   bins | solutions | seconds | max_rss_mib |
|:----------|----------------:|--------:|-------:|----------:|--------:|------------:|
| 1x bins   |             136 |       1 |  29920 |        10 |   6.360 |      49.867 |
| 2x bins   |             272 |       1 |  59840 |        10 |  11.575 |      51.398 |
| 4x bins   |             544 |       1 | 119680 |        10 |  23.383 |      55.090 |
| 2 samples |             136 |       2 |  59840 |        20 |   6.365 |      52.938 |
| 4 samples |             136 |       4 | 119680 |        40 |  18.505 |      56.230 |

Fresh processes, four threads, spill forbidden

The latency exponent over 1×, 2× and 4× bins is 0.94; four samples of the
1× size take 2.91× the one-sample latency; S1
fails. One decode runs per sample on the thread that holds that
sample’s row, so parallelism across samples depends on how DuckDB spreads the per-sample rows
over threads; it is not a property of the kernel.

## Findings

- **Declared criteria.** A1 holds for 12 of 16 composites, A2
  for 12, A3 for 8; A4 has
  88.1% of starts within tolerance
  (fails). G1 holds for 1 of
  7 gain composites; N1 for 8 of 9 pure
  composites. S1 fails (exponent 0.94, peak
  RSS 56 MiB, four samples at
  2.91×).
- **Where the two models land on the same optimum they agree exactly.** In
  5 of 7 gain
  composites the port and upstream give chromosome 18 the same copy number and the same fraction to
  four decimals, with every bin in the same state. The disagreements sit on flat objectives: pure
  samples and the smallest gains, where the normal fraction is near 1 and ploidy and states are not
  identified, or starts where one optimizer climbed to a different local maximum. Among the starts
  outside A4, the port’s objective ended higher than upstream’s by more than 1 in 5 and
  lower in 1; both optimizers find local maxima, as the design note says.
- **A single altered chromosome does not identify the fraction.** Both models see the planted
  shift (chromosome 18 median log2 ratio against the designed `log2(1 + m/2)` in the table) and
  explain it as a higher copy number in a smaller share of cells: extra copies times fraction is
  close to `m`, the called copy number is 3 only when `m` is small enough for the prior to prefer it.
  This is the model’s known degeneracy on genomes with one event, shared with upstream, so G1 as
  declared (three copies at fraction `m`) fails for both; a fraction is identified when a genome
  carries changes of more than one magnitude, which this design does not plant.
- **Samples without change report no fraction** (at most
  0.019 from the port, 0.025 from upstream); the guard fired
  on 1 of them, the rest converged to a normal fraction near
  1, whose ploidy and copy numbers are not evidence.

Every table above is written under `benchmarks/results/depth-state/`; the per-segment tables are
gzip-compressed.
