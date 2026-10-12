ROH from read counts: the fitted emission on chr20
================

`duckhts_roh_counts` takes `allele_balance` and `artefact_weight`, and
`duckhts_roh_counts_fitted` reads them, with `seq_error` and
`contamination`, from the output of `duckhts_count_error_fit`
(<https://github.com/RGenomicsETL/duckhts/issues/329>,
<https://github.com/RGenomicsETL/duckhts/issues/332>). The SQL tests check the
emission against an independent oracle and on a synthetic fixture. This
report checks it on real data, with two questions:

- **Agreement.** Does the full fit bring the read-count decode closer to the
  genotype decode than the scalar defaults or the fitted error and
  contamination alone?
- **Sensitivity.** Does a change of one fitted number, within a range tied to
  its uncertainty, move the runs?

The scope is small on purpose: chr20, three samples. The report makes no
timing or scaling claim; it measures memory only to gate it.

## Inputs

The inputs are those of `benchmark_roh_counts_validation.md`, from the
`r/duckhtsbench` registry:

- the read counts of NA18507, HG00403 and HG00188 at the 108,757 chr20
  evaluation sites (`roh_counts_chr20_na18507`, `roh_counts_chr20_hg00403`,
  `roh_counts_chr20_hg00188`), from the 1000 Genomes 30× CRAMs, with the
  site’s `af`;
- the genotypes of the same samples at the same sites
  (`roh_counts_validation_chr20_genotypes`), FORMAT/GT only.

Before anything is decoded, the render compares each file with the identity
its registry row records (`benchmarks/roh_counts_validation/inputs.R`). A
file that differs stops the render. The genotype decode takes the `af` of the
first sample’s counts as its frequency relation, as in the validation report.

## Design, declared before any fitted decode

The values below are in `benchmarks/roh_counts_validation/fitted_declaration.R`,
committed on its own before any decode with a fitted `allele_balance` or
`artefact_weight` was run, and before this report. This section states them
in words.

**Fit.** `duckhts_count_error_fit` with its default options, on the counts of
all decoded samples. A sample whose `status` is not `ok` is not decoded with
its fit.

**Decodes.** Every count decode is `duckhts_roh_counts_fitted` with a fit
relation of four numbers per sample:

| Decode                               | `seq_error` | `contamination` | `allele_balance` | `artefact_weight` |
|:-------------------------------------|:------------|:----------------|:-----------------|:------------------|
| \(a\) scalar defaults                | 0.001       | 0               | 0.5              | 0                 |
| \(b\) fitted error and contamination | fit         | fit             | 0.5              | 0                 |
| \(c\) full fit                       | fit         | fit             | fit              | fit               |

Decode (a) must equal `duckhts_roh_counts` with its scalar defaults, run for
run; the render stops otherwise. The genotype decodes are
`duckhts_roh_af_table(..., gt_error := 30)` and `gt_error := 100`, as in the
validation report.

**Metrics**, per sample, all runs, as in the validation report: FROH (bases in
runs over the span from the first to the last site), the FROH difference test
minus reference, and the Jaccard index of run bases (bases in the runs of both
decodes over bases in the runs of either). Runs of at least 1 Mb are in the
CSV files.

**Criterion.** On every sample, the full fit (c) is at least as close to the
genotype decode at `gt_error := 30` as (a) and as (b), by the absolute
difference of bases in runs (the absolute FROH difference times the span). It
is a verdict of the report, not a stop of the render. The comparison at
`gt_error := 100` is reported without a criterion. Each row also shows the
validation report’s tolerance (absolute FROH difference at most 0.01, Jaccard
index at least 0.7); it is not the criterion here.

**Sensitivity.** Decode (c) is repeated with one fitted number changed and the
three others kept:

- `allele_balance` −0.01, −0.005, +0.005 and +0.01;
- `artefact_weight` × 0.5 and × 2;
- `seq_error` × 0.5 and × 2;
- `contamination` − and + `contamination_sd` (0.001 when it is NULL).

A changed value outside the range that `duckhts_roh_counts` accepts is set to
the end of that range and the row is marked. For each sample and change, the
reference is decode (c):

- the FROH change and the change of bases in runs;
- the changed share: bases in the runs of exactly one of the two decodes,
  over the bases in the runs of (c).

**Bound**, declared: on every sample and change, the absolute FROH change is
below 0.005 and the changed share is below 2%. A row that exceeds it is
reported with the number that drives it. The bound is a verdict, not a stop of
the render, and it is not changed after the run.

**Stability**, declared after the first render (it found the runs most
sensitive to `artefact_weight`) in its own commit, and before any decode with
an `artefact_weight_sd`. `duckhts_count_error_fit` now reports
`artefact_weight_sd`, and `duckhts_roh_counts_fitted` gives each run
`stable_fraction`: the share of its bases that are also in a run when the
sample is decoded at `artefact_weight` − and + `artefact_weight_sd` (kept in
\[0, 0.05\]). Decode (c) passes each sample’s fitted `artefact_weight_sd`; every
other count decode, the perturbed ones included, passes 0. The report gives,
per sample and without a verdict:

- `artefact_weight_sd` and the two weights it gives, beside the declared × 0.5
  and × 2;
- the `stable_fraction` of the runs of (c): the runs, the runs at 1, its
  quantiles, and the stable bases (run length times `stable_fraction`,
  summed) over the run bases;
- for the `artefact_weight` × 0.5 and × 2 rows: the lost bases (in runs of (c)
  that the changed decode does not call), split by whether the run that holds
  them has `stable_fraction` 1 or below 1, and the gained bases (called by the
  changed decode outside any run of (c)), which no `stable_fraction` covers.
  The FROH change is gained minus lost over the span; the changed bases of
  the sensitivity table are gained plus lost.

**Gained bases**, declared after the stability render and before any decode
that reports `gained_bases`, in its own commit. `duckhts_roh_counts_fitted`
now also gives each run `gained_bases`: the bases that its decode at
w + `artefact_weight_sd` calls in runs that meet the run (share a base or
abut) but outside every run at w. Runs at w + sd that meet no run at w have
no row. The driver decodes once more at w + `artefact_weight_sd` with
`artefact_weight_sd` 0 to count them. The report gives, per sample and
without a verdict:

- the `gained_bases` of the runs of (c), summed, and the bases of the w + sd
  decode in runs that meet no run of (c). The summed column must equal the
  bases of that extra decode outside the runs of (c) in runs that meet one,
  or the render stops;
- for the `artefact_weight` × 2 rows: the gained bases of the decomposition,
  split into those in runs that meet a run of (c) and those in new runs that
  meet none.

**Gates.** Each query (the fit, every decode, three empty runs) runs in a
fresh DuckDB CLI process under `/usr/bin/time`, with the JSON profiler and
`max_temp_directory_size` set to zero, at two threads. The render stops if an
input differs from its registry identity, if any query spills, or if its peak
RSS exceeds its ceiling. Each ceiling is the overhead (the median peak RSS of
the three empty runs) plus 3 × the live state of the query, as the scaling
reports declare it:

- **count decode**, `benchmark_roh_scaling.md`: decoded samples × sites × 24
  bytes for the packed site lists, plus threads × sites × 36.25 bytes for the
  kernel workspace;
- **genotype decode** (`duckhts_roh_af_table`), the same report: samples ×
  frequency rows × 16 + frequency rows × 128 + threads × frequency rows ×
  36.25 bytes;
- **fit**, `benchmark_count_error_fit.md`: histogram cells × 312 + distinct
  (depth, alt count) pairs × 88 + samples × 1,001 × 24 bytes, with the cells
  and pairs counted untimed from the input.

`ROH_FITTED_SIZE=reduced` checks the driver: it decodes one sample
(NA18507) instead of three and writes its observations to a temporary
directory instead of `benchmarks/results/roh-fitted-emission/`. Every input
check and gate still runs. A reduced run is not evidence. The default is
`full`.

``` sh
nice -n -15 taskset -c 8-15 Rscript -e 'rmarkdown::render("benchmarks/benchmark_roh_fitted_emission.Rmd")'
# Driver check only; writes the page outside benchmarks/:
ROH_FITTED_SIZE=reduced taskset -c 0-3 Rscript -e \
  'rmarkdown::render("benchmarks/benchmark_roh_fitted_emission.Rmd", output_dir = tempdir())'
```

## Run

| Item                         | Value                                                            |
|:-----------------------------|:-----------------------------------------------------------------|
| revision                     | bc2bc025074df976ebcdacafd458e3c54197f358                         |
| src_tree                     | d8d9a75ea2378d58e4524438db2e204272c6b76f                         |
| input_identities             | match the registry                                               |
| extension_sha256             | e983edbbe3c312f6e8eac9fef495a0bb060cfa711bf2f9e1c2b2d6f885de48d3 |
| duckdb_cli                   | v1.5.1 (Variegata) 7dbb2e646f                                    |
| r_version                    | R version 4.6.0 (2026-04-24)                                     |
| host_load                    | 1.11, 1.18, 0.88                                                 |
| cpu_affinity                 | 8-15                                                             |
| threads                      | 2                                                                |
| samples                      | NA18507,HG00403,HG00188                                          |
| decoded_samples              | HG00188,HG00403,NA18507                                          |
| sites                        | 108757                                                           |
| span_bases                   | 64251060                                                         |
| fit_cells                    | 126676                                                           |
| fit_pairs                    | 3124                                                             |
| overhead_mib                 | 39.8                                                             |
| defaults_identical_to_scalar | TRUE                                                             |

## Fit

| Sample  |   Sites | Depth |   seq_error | contamination | contamination_sd | allele_balance | artefact_weight | artefact_weight_sd | status |
|:--------|--------:|------:|------------:|--------------:|-----------------:|---------------:|----------------:|-------------------:|:-------|
| HG00188 | 108,755 |  31.6 | 0.000119226 |   0.000325977 |      6.59711e-05 |       0.497437 |      0.00182514 |        0.000939117 | ok     |
| HG00403 | 108,736 |  31.9 | 0.000129028 |   0.000318597 |      6.58372e-05 |        0.49797 |      0.00295505 |        0.000563165 | ok     |
| NA18507 | 108,723 |  33.0 | 0.000140831 |   0.000262749 |      7.39054e-05 |       0.498413 |       0.0026572 |         0.00134403 | ok     |

## Agreement with the genotype decode

`reference` is the genotype decode at the `gt_error` of the row and `test` is
the count decode; all runs.

| gt_error | Sample  | Decode                               | Genotype FROH | Count FROH | FROH difference | Jaccard | Validation tolerance |
|---------:|:--------|:-------------------------------------|--------------:|-----------:|----------------:|--------:|:---------------------|
|       30 | NA18507 | \(a\) scalar defaults                |        0.1772 |     0.1469 |         -0.0304 |   0.824 | fail                 |
|       30 | NA18507 | \(b\) fitted error and contamination |        0.1772 |     0.1460 |         -0.0313 |   0.819 | fail                 |
|       30 | NA18507 | \(c\) full fit                       |        0.1772 |     0.1771 |         -0.0002 |   0.975 | pass                 |
|       30 | HG00403 | \(a\) scalar defaults                |        0.3411 |     0.3030 |         -0.0381 |   0.883 | fail                 |
|       30 | HG00403 | \(b\) fitted error and contamination |        0.3411 |     0.3020 |         -0.0391 |   0.881 | fail                 |
|       30 | HG00403 | \(c\) full fit                       |        0.3411 |     0.3446 |         +0.0035 |   0.987 | pass                 |
|       30 | HG00188 | \(a\) scalar defaults                |        0.3382 |     0.2822 |         -0.0559 |   0.834 | fail                 |
|       30 | HG00188 | \(b\) fitted error and contamination |        0.3382 |     0.2813 |         -0.0569 |   0.831 | fail                 |
|       30 | HG00188 | \(c\) full fit                       |        0.3382 |     0.3353 |         -0.0029 |   0.989 | pass                 |
|      100 | NA18507 | \(a\) scalar defaults                |        0.1477 |     0.1469 |         -0.0009 |   0.986 | pass                 |
|      100 | NA18507 | \(b\) fitted error and contamination |        0.1477 |     0.1460 |         -0.0018 |   0.980 | pass                 |
|      100 | NA18507 | \(c\) full fit                       |        0.1477 |     0.1771 |         +0.0293 |   0.834 | fail                 |
|      100 | HG00403 | \(a\) scalar defaults                |        0.3054 |     0.3030 |         -0.0024 |   0.975 | pass                 |
|      100 | HG00403 | \(b\) fitted error and contamination |        0.3054 |     0.3020 |         -0.0034 |   0.973 | pass                 |
|      100 | HG00403 | \(c\) full fit                       |        0.3054 |     0.3446 |         +0.0392 |   0.885 | fail                 |
|      100 | HG00188 | \(a\) scalar defaults                |        0.2856 |     0.2822 |         -0.0034 |   0.981 | pass                 |
|      100 | HG00188 | \(b\) fitted error and contamination |        0.2856 |     0.2813 |         -0.0043 |   0.978 | pass                 |
|      100 | HG00188 | \(c\) full fit                       |        0.2856 |     0.3353 |         +0.0497 |   0.849 | fail                 |

### Criterion

The gap is the absolute difference of bases in runs between the count decode
and the genotype decode at `gt_error := 30`.

| Sample  |   Gap (a) |   Gap (b) | Gap (c) | Verdict |
|:--------|----------:|----------:|--------:|:--------|
| HG00188 | 3,593,146 | 3,654,279 | 184,073 | pass    |
| HG00403 | 2,448,479 | 2,509,848 | 224,108 | pass    |
| NA18507 | 1,951,952 | 2,009,195 |  10,731 | pass    |

The criterion holds on 3 of 3
samples: it holds.

## Sensitivity of the full fit

Each row changes one number of decode (c) for one sample. The reference is
decode (c) itself.

| Number          | Change | Sample  |      Fitted |     Changed | At guard | FROH change | Bases change | Changed bases | Changed share | Verdict |
|:----------------|:-------|:--------|------------:|------------:|:---------|------------:|-------------:|--------------:|--------------:|:--------|
| allele_balance  | -0.01  | NA18507 |    0.498413 |    0.488413 |          |    +0.00049 |       31,595 |        66,031 |         0.58% | pass    |
| allele_balance  | -0.01  | HG00403 |     0.49797 |     0.48797 |          |    -0.00104 |      -66,906 |        88,180 |         0.40% | pass    |
| allele_balance  | -0.01  | HG00188 |    0.497437 |    0.487437 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |
| allele_balance  | -0.005 | NA18507 |    0.498413 |    0.493413 |          |    +0.00000 |          251 |         7,983 |         0.07% | pass    |
| allele_balance  | -0.005 | HG00403 |     0.49797 |     0.49297 |          |    +0.00010 |        6,202 |         6,202 |         0.03% | pass    |
| allele_balance  | -0.005 | HG00188 |    0.497437 |    0.492437 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |
| allele_balance  | +0.005 | NA18507 |    0.498413 |    0.503413 |          |    -0.00135 |      -86,990 |        86,990 |         0.76% | pass    |
| allele_balance  | +0.005 | HG00403 |     0.49797 |     0.50297 |          |    -0.00004 |       -2,870 |         2,870 |         0.01% | pass    |
| allele_balance  | +0.005 | HG00188 |    0.497437 |    0.502437 |          |    +0.00012 |        7,717 |         7,717 |         0.04% | pass    |
| allele_balance  | +0.01  | NA18507 |    0.498413 |    0.508413 |          |    -0.00135 |      -86,990 |        86,990 |         0.76% | pass    |
| allele_balance  | +0.01  | HG00403 |     0.49797 |     0.50797 |          |    -0.00004 |       -2,870 |         2,870 |         0.01% | pass    |
| allele_balance  | +0.01  | HG00188 |    0.497437 |    0.507437 |          |    +0.00032 |       20,269 |        20,269 |         0.09% | pass    |
| artefact_weight | × 0.5  | NA18507 |   0.0026572 |   0.0013286 |          |    -0.00652 |     -419,209 |       419,209 |         3.68% | fail    |
| artefact_weight | × 0.5  | HG00403 |  0.00295505 |  0.00147753 |          |    -0.00310 |     -199,147 |       199,147 |         0.90% | pass    |
| artefact_weight | × 0.5  | HG00188 |  0.00182514 | 0.000912572 |          |    -0.00375 |     -240,877 |       240,877 |         1.12% | pass    |
| artefact_weight | × 2    | NA18507 |   0.0026572 |  0.00531439 |          |    +0.00495 |      318,363 |       318,363 |         2.80% | fail    |
| artefact_weight | × 2    | HG00403 |  0.00295505 |  0.00591011 |          |    +0.00280 |      179,992 |       179,992 |         0.81% | pass    |
| artefact_weight | × 2    | HG00188 |  0.00182514 |  0.00365029 |          |    +0.00294 |      188,954 |       188,954 |         0.88% | pass    |
| seq_error       | × 0.5  | NA18507 | 0.000140831 | 7.04154e-05 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |
| seq_error       | × 0.5  | HG00403 | 0.000129028 | 6.45139e-05 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |
| seq_error       | × 0.5  | HG00188 | 0.000119226 |  5.9613e-05 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |
| seq_error       | × 2    | NA18507 | 0.000140831 | 0.000281662 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |
| seq_error       | × 2    | HG00403 | 0.000129028 | 0.000258056 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |
| seq_error       | × 2    | HG00188 | 0.000119226 | 0.000238452 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |
| contamination   | -1 sd  | NA18507 | 0.000262749 | 0.000188844 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |
| contamination   | -1 sd  | HG00403 | 0.000318597 |  0.00025276 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |
| contamination   | -1 sd  | HG00188 | 0.000325977 | 0.000260006 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |
| contamination   | +1 sd  | NA18507 | 0.000262749 | 0.000336655 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |
| contamination   | +1 sd  | HG00403 | 0.000318597 | 0.000384434 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |
| contamination   | +1 sd  | HG00188 | 0.000325977 | 0.000391948 |          |    +0.00000 |            0 |             0 |         0.00% | pass    |

| Number          | Largest FROH change | Largest changed share | Rows within the bound |
|:----------------|--------------------:|----------------------:|----------------------:|
| allele_balance  |             0.00135 |                 0.76% |              12 of 12 |
| artefact_weight |             0.00652 |                 3.68% |                4 of 6 |
| seq_error       |             0.00000 |                 0.00% |                6 of 6 |
| contamination   |             0.00000 |                 0.00% |                6 of 6 |

28 of 30 rows are within the
declared bound (absolute FROH change below 0.005,
changed share below 2%).
The rows outside it change: `artefact_weight`.

## Stability of the full fit

`stable_fraction` of the runs of decode (c), from its fitted
`artefact_weight_sd`. The two weights are those of the stability decodes; the
declared perturbations are × 0.5 and × 2.

| Sample  | artefact_weight | artefact_weight_sd | w − sd      | w + sd     | × 0.5       | × 2        | Runs | Runs at 1 | q0    | q0.25 | q0.5  | q0.75 | q1    | Stable bases |
|:--------|:----------------|:-------------------|:------------|:-----------|:------------|:-----------|-----:|----------:|:------|:------|:------|:------|:------|:-------------|
| HG00188 | 0.00182514      | 0.000939117        | 0.000886026 | 0.00276426 | 0.000912572 | 0.00365029 |  178 |       170 | 0.154 | 1.000 | 1.000 | 1.000 | 1.000 | 98.8%        |
| HG00403 | 0.00295505      | 0.000563165        | 0.00239189  | 0.00351822 | 0.00147753  | 0.00591011 |  189 |       185 | 0.000 | 1.000 | 1.000 | 1.000 | 1.000 | 99.3%        |
| NA18507 | 0.0026572       | 0.00134403         | 0.00131317  | 0.00400122 | 0.0013286   | 0.00531439 |  149 |       131 | 0.000 | 1.000 | 1.000 | 1.000 | 1.000 | 96.3%        |

How the changed bases of the `artefact_weight` rows of the sensitivity table
split. Lost bases are in runs of (c) that the changed decode does not call;
gained bases are called outside any run of (c).

| Change | Sample  | FROH change |    Lost | Lost, stable_fraction 1 | Lost, stable_fraction \< 1 |  Gained | Lost share below 1 | Changed share below 1 |
|:-------|:--------|------------:|--------:|------------------------:|---------------------------:|--------:|-------------------:|----------------------:|
| × 0.5  | NA18507 |    -0.00652 | 419,209 |                       0 |                    419,209 |       0 |             100.0% |                100.0% |
| × 0.5  | HG00403 |    -0.00310 | 199,147 |                  53,309 |                    145,838 |       0 |              73.2% |                 73.2% |
| × 0.5  | HG00188 |    -0.00375 | 240,877 |                       0 |                    240,877 |       0 |             100.0% |                100.0% |
| × 2    | NA18507 |    +0.00495 |       0 |                       0 |                          0 | 318,363 |               none |                  0.0% |
| × 2    | HG00403 |    +0.00280 |       0 |                       0 |                          0 | 179,992 |               none |                  0.0% |
| × 2    | HG00188 |    +0.00294 |       0 |                       0 |                          0 | 188,954 |               none |                  0.0% |

### Gained bases

`gained_bases` of the runs of (c), from its decode at w + `artefact_weight_sd`,
and the bases of that decode in runs that meet no run of (c), which no row
reports.

| Sample  |     w + sd | Runs with a gain | gained_bases, summed | In new runs, not reported |
|:--------|-----------:|-----------------:|---------------------:|--------------------------:|
| HG00188 | 0.00276426 |                6 |               66,465 |                         0 |
| HG00403 | 0.00351822 |                3 |               19,185 |                         0 |
| NA18507 | 0.00400122 |                8 |              111,976 |                    65,969 |

The gained bases of the `artefact_weight` × 2 rows, by whether their run at
× 2 meets a run of (c).

| Change | Sample  |     Weight |  Gained | Next to a run of (c) | In new runs | Share next to a run |
|:-------|:--------|-----------:|--------:|---------------------:|------------:|--------------------:|
| × 2    | NA18507 | 0.00531439 | 318,363 |              189,073 |     129,290 |               59.4% |
| × 2    | HG00403 | 0.00591011 | 179,992 |              179,992 |           0 |              100.0% |
| × 2    | HG00188 | 0.00365029 | 188,954 |              124,361 |      64,593 |               65.8% |

## Gates

Every query, with its latency, peak RSS, DuckDB’s peak buffer memory, its
ceiling and its temporary bytes. The render stops on any spill or any peak
RSS above its ceiling, and it reached this table.

| Query           | Kind            | Seconds | Peak RSS MiB | Peak buffer MiB | Ceiling MiB | Temporary bytes |
|:----------------|:----------------|--------:|-------------:|----------------:|------------:|----------------:|
| empty_1         | empty           |   0.000 |         39.7 |             0.2 |             |               0 |
| empty_2         | empty           |   0.007 |         39.8 |             0.2 |             |               0 |
| empty_3         | empty           |   0.000 |         39.9 |             0.2 |             |               0 |
| fit             | fit             |   3.328 |         69.3 |            39.8 |       153.9 |               0 |
| genotype_gt30   | genotype_decode |   0.178 |         93.2 |            97.1 |       117.1 |               0 |
| genotype_gt100  | genotype_decode |   0.165 |         93.5 |            97.1 |       117.1 |               0 |
| decode_a        | count_decode    |   0.077 |         66.7 |            40.5 |        84.7 |               0 |
| decode_b        | count_decode    |   0.078 |         70.7 |            40.5 |        84.7 |               0 |
| decode_c        | count_decode    |   0.109 |         70.9 |            39.9 |        84.7 |               0 |
| scalar_defaults | count_decode    |   0.037 |         61.0 |            36.9 |        84.7 |               0 |
| perturbed_1     | count_decode    |   0.113 |         67.5 |            40.5 |        84.7 |               0 |
| perturbed_2     | count_decode    |   0.113 |         65.8 |            41.0 |        84.7 |               0 |
| perturbed_3     | count_decode    |   0.109 |         67.2 |            40.5 |        84.7 |               0 |
| perturbed_4     | count_decode    |   0.148 |         69.6 |            29.4 |        84.7 |               0 |
| perturbed_5     | count_decode    |   0.155 |         65.7 |            29.4 |        84.7 |               0 |
| perturbed_6     | count_decode    |   0.155 |         65.2 |            29.4 |        84.7 |               0 |
| perturbed_7     | count_decode    |   0.109 |         70.0 |            40.5 |        84.7 |               0 |
| perturbed_8     | count_decode    |   0.109 |         71.5 |            40.5 |        84.7 |               0 |
| perturbed_9     | count_decode    |   0.112 |         70.8 |            40.5 |        84.7 |               0 |
| perturbed_10    | count_decode    |   0.149 |         67.7 |            41.1 |        84.7 |               0 |
| stability_high  | count_decode    |   0.151 |         68.4 |            29.4 |        84.7 |               0 |

Every observation is in `benchmarks/results/roh-fitted-emission/`: the fit, the four numbers of every decode, every run, the metrics and the queries.

## Limitations

- One chromosome and three samples. The rows are observations, not rates.
- The genotype decode is a reference, not the truth: its genotypes come from
  the same sequencing, with joint calling and phasing behind them. The
  criterion asks whether the fit moves the read-count decode toward it, not
  whether either is right.
- The fit and the decode use the same counts. The fitted numbers are not
  tested on data they were not fitted on.
- The perturbations change one number at a time. They do not test changes of
  several numbers together.
- `contamination` and `artefact_weight` have intervals from the fit
  (`contamination_sd`, `artefact_weight_sd`). The ranges of the other two
  numbers (± 0.005 and ± 0.01, × 0.5 and × 2) are declared, not derived from
  an interval of the fit, and the × 0.5 and × 2 of `artefact_weight` were
  declared before its interval existed.
- `stable_fraction` measures only the bases of a run that a change of the
  artefact weight would lose. `gained_bases` counts the bases that the higher
  weight of the interval adds next to a run; bases it calls in new runs that
  meet no run have no row, and this report counts them separately.
- The spreads of the fit (`spread_hom`, `spread_het`) are not used by the
  decode.
