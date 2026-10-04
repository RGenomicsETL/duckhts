ROH from read counts: real-data validation on chr20
================

`duckhts_roh_counts` finds runs of homozygosity from allele read counts, with a
per-read error and a contamination fraction. Its SQL tests compare it with the
same model computed in SQL and with a synthetic titration. This report adds two
checks on real data, as section 2 of
<https://github.com/RGenomicsETL/duckhts/issues/329> asks:

- **(a) Agreement.** The read-count path against the VCF genotype path, on the
  same samples, sites and allele frequencies.
- **(b) Titration.** Real read counts of two samples mixed at a known fraction,
  decoded without and with the contamination term.

The scope is small on purpose: chr20, three samples, no timing. It makes no
performance claim.

## Inputs

All inputs are `r/duckhtsbench` registry artifacts.

- Read counts of NA18507, HG00403 and HG00188 at the 108,757 chr20 evaluation
  sites of the ROH ancestry evaluation: `roh_counts_chr20_na18507`,
  `roh_counts_chr20_hg00403` and `roh_counts_chr20_hg00188`. They come from the
  1000 Genomes 30× CRAMs (`benchmarks/roh_counts_stage.R`). Each row has the
  site’s `af`, the source VCF’s INFO/AF of the ALT allele, clamped to
  \[0.001, 0.999\].
- Genotypes of the same three samples at the same sites:
  `roh_counts_validation_chr20_genotypes`. It is derived from the phased chr20
  VCF of the 3,202 samples (`roh_ancestry_chr20_source`) and the site list
  (`roh_ancestry_chr20_af_sites`) by
  `benchmarks/roh_counts_validation/stage_genotypes.R`: three samples, FORMAT/GT
  only, exactly one record per site. Staging checks the SHA-256 of the cached
  source and the record digest of the output.

The read counts and the genotypes come from the same 30× sequencing of the
1000 Genomes samples. The VCF path is an independent decode path, with joint calling and phasing behind its
genotypes. It is not an independent measurement of the samples.

## Design, declared before any run was decoded

The values below are in `benchmarks/roh_counts_validation/declaration.R`. That
file and this section were committed before the first decode.

**Denominator.** FROH is the bases in runs divided by the span from the first
to the last evaluation site on chr20, inclusive. A run’s bases are
`end - start + 1`.

**Run classes.** Every metric is given for all runs and for long runs, those of
at least 1,000,000 bases. Each decode is filtered by its own run lengths.

**Metrics**, per sample and run class, for a test decode against a reference:

- FROH of each decode and the difference, test minus reference;
- the base-level overlap as a Jaccard index: bases in the runs of both decodes
  divided by bases in the runs of either. It is undefined when neither decode
  has a run in the class.

**Tolerance.**

- The absolute FROH difference is at most 0.01.
- The Jaccard index is at least 0.9 for long runs and at least 0.7 for all runs.
- A class in which neither decode has a run agrees. A class in which only one
  decode has a run has a Jaccard index of 0 and fails.

Justification: both paths use the same sites, frequencies and hidden Markov
model, and differ only in the genotype evidence (a called genotype at phred 30,
or about 30 reads at the default `seq_error`). A long run rests on hundreds of
sites and should not depend on that difference. A short run rests on a few
sites, so its ends can move.

**(a) Agreement.** The reference is `duckhts_roh_af_table(..., gt_error := 30)`
on the genotype BCF, with the `af` of the counts as the frequency relation. The
test is `duckhts_roh_counts` on the counts with the default `seq_error` (0.001)
and `contamination := 0`. A sample passes a run class when it meets the
tolerance.

**(b) Titration.** For each ordered pair (receiver, contaminant) of the three
samples and each α in 1%, 2%, 5% and 10%, the mixed count of each allele at
each site is

    Binomial(receiver count, 1 − α) + Binomial(contaminant count, α).

This thins real counts and adds them. It is not a mixture of reads in a CRAM:
no read is realigned, filtered or counted again, and the two samples’ counts
are drawn independently at each site. The mean depth stays close to that of
one sample.

Cell *k* of the grid, in the order receiver, contaminant, α, draws with
`set.seed(20261004 + k)` and R’s default generators. The draws are, in order:
receiver REF, contaminant REF, receiver ALT, contaminant ALT.

Each mixture is decoded twice: with `contamination := 0` and with
`contamination := α`. The reference is the decode of the receiver’s unmixed
counts from (a). The same metrics and tolerance apply. The declared expectation
is for the decode with `contamination := α`; the decode without the term is the
control that shows what the term changes.

The contaminant is one individual, and `duckhts_roh_counts` models a
contaminant read by the population frequency `af`. That approximation is part
of what (b) measures.

All observations are kept: every run of every decode is in
`benchmarks/results/roh-counts-validation/`.

``` sh
Rscript benchmarks/roh_counts_validation/stage.R
taskset -c 4-7 Rscript -e 'rmarkdown::render("benchmarks/benchmark_roh_counts_validation.Rmd")'
```

## Run

| Item                | Value                                    |
|:--------------------|:-----------------------------------------|
| revision            | bdddcaeeeddb5dfc4b582961dfcd13dc3ac510d8 |
| tracked_changes     | no                                       |
| extension_version   | 1.5.2.9011                               |
| duckdb_version      | v1.5.5                                   |
| r_version           | R version 4.6.0 (2026-04-24)             |
| rng_kind            | Mersenne-Twister/Inversion/Rejection     |
| seed                | 20261004                                 |
| threads             | 2                                        |
| sites               | 108757                                   |
| first_pos           | 80457                                    |
| last_pos            | 64331516                                 |
| span_bases          | 64251060                                 |
| gt_error            | 30                                       |
| long_run_bases      | 1000000                                  |
| max_froh_difference | 0.01                                     |
| min_jaccard_long    | 0.9                                      |
| min_jaccard_all     | 0.7                                      |

## (a) Read counts against VCF genotypes

`reference` is the VCF genotype path and `test` is the read-count path.

| Sample  | Run class   | VCF runs | Count runs | VCF FROH | Count FROH | FROH difference | Shared bases | Bases in either | Jaccard | FROH | Overlap | Verdict |
|:--------|:------------|---------:|-----------:|---------:|-----------:|----------------:|-------------:|----------------:|--------:|:-----|:--------|:--------|
| NA18507 | all runs    |      150 |        154 |   0.1772 |     0.1469 |         -0.0304 |    9,406,217 |      11,416,643 |   0.824 | fail | pass    | fail    |
| NA18507 | runs ≥ 1 Mb |        0 |          0 |   0.0000 |     0.0000 |         +0.0000 |            0 |               0 |    none | pass | pass    | pass    |
| HG00403 | all runs    |      188 |        233 |   0.3411 |     0.3030 |         -0.0381 |   19,408,034 |      21,971,637 |   0.883 | fail | pass    | fail    |
| HG00403 | runs ≥ 1 Mb |        0 |          0 |   0.0000 |     0.0000 |         +0.0000 |            0 |               0 |    none | pass | pass    | pass    |
| HG00188 | all runs    |      178 |        215 |   0.3382 |     0.2822 |         -0.0559 |   18,128,775 |      21,733,049 |   0.834 | fail | pass    | fail    |
| HG00188 | runs ≥ 1 Mb |        1 |          0 |   0.0195 |     0.0000 |         -0.0195 |            0 |       1,252,351 |   0.000 | fail | fail    | fail    |

2 of 6 sample and run-class rows meet
the declared tolerance. “none” means that neither decode has a run in the
class.

## (b) Count-level contamination titration

Each row summarises the six ordered pairs at one α and one decode. The FROH
difference and the Jaccard index are against the receiver’s unmixed decode;
they are shown as the median and the range over the pairs.

| Run class   | α   | Contamination term | FROH difference              | Jaccard                | Pairs without runs | Pairs passing |
|:------------|:----|:-------------------|:-----------------------------|:-----------------------|-------------------:|:--------------|
| all runs    | 1%  | 0                  | -0.0003 (-0.0008 to +0.0000) | 0.998 (0.996 to 1.000) |                  0 | 6 of 6        |
| all runs    | 1%  | α                  | +0.0018 (+0.0008 to +0.0027) | 0.994 (0.990 to 0.995) |                  0 | 6 of 6        |
| all runs    | 2%  | 0                  | -0.0022 (-0.0055 to -0.0013) | 0.991 (0.975 to 0.995) |                  0 | 6 of 6        |
| all runs    | 2%  | α                  | +0.0016 (+0.0002 to +0.0032) | 0.990 (0.987 to 0.995) |                  0 | 6 of 6        |
| all runs    | 5%  | 0                  | -0.0445 (-0.0826 to -0.0219) | 0.832 (0.727 to 0.852) |                  0 | 0 of 6        |
| all runs    | 5%  | α                  | +0.0018 (+0.0005 to +0.0050) | 0.984 (0.972 to 0.986) |                  0 | 6 of 6        |
| all runs    | 10% | 0                  | -0.1804 (-0.2692 to -0.1170) | 0.197 (0.111 to 0.402) |                  0 | 0 of 6        |
| all runs    | 10% | α                  | -0.0002 (-0.0087 to +0.0020) | 0.939 (0.919 to 0.953) |                  0 | 6 of 6        |
| runs ≥ 1 Mb | 1%  | 0                  | +0.0000 (+0.0000 to +0.0000) | none                   |                  6 | 6 of 6        |
| runs ≥ 1 Mb | 1%  | α                  | +0.0000 (+0.0000 to +0.0000) | none                   |                  6 | 6 of 6        |
| runs ≥ 1 Mb | 2%  | 0                  | +0.0000 (+0.0000 to +0.0000) | none                   |                  6 | 6 of 6        |
| runs ≥ 1 Mb | 2%  | α                  | +0.0000 (+0.0000 to +0.0000) | none                   |                  6 | 6 of 6        |
| runs ≥ 1 Mb | 5%  | 0                  | +0.0000 (+0.0000 to +0.0000) | none                   |                  6 | 6 of 6        |
| runs ≥ 1 Mb | 5%  | α                  | +0.0000 (+0.0000 to +0.0000) | none                   |                  6 | 6 of 6        |
| runs ≥ 1 Mb | 10% | 0                  | +0.0000 (+0.0000 to +0.0000) | none                   |                  6 | 6 of 6        |
| runs ≥ 1 Mb | 10% | α                  | +0.0000 (+0.0000 to +0.0000) | none                   |                  6 | 6 of 6        |

With `contamination := α`, 48 of 48 pair,
α and run-class rows meet the declared tolerance. With `contamination := 0`,
36 of 48 do.

The realised contaminant read fraction of the mixtures, over all sites, is
within 0.0037
of α in every cell, and the mean depth of a mixture is between
31.6 and
33.0 reads per site.

### Every pair

| Run class   | Receiver | Contaminant | alpha | Contamination term | Reference runs | Runs | Reference FROH |   FROH | FROH difference | Jaccard | Verdict |
|:------------|:---------|:------------|:------|:-------------------|---------------:|-----:|---------------:|-------:|----------------:|--------:|:--------|
| all runs    | NA18507  | HG00403     | 1%    | 0                  |            154 |  153 |         0.1469 | 0.1463 |         -0.0005 |   0.996 | pass    |
| all runs    | NA18507  | HG00403     | 1%    | α                  |            154 |  150 |         0.1469 | 0.1476 |         +0.0008 |   0.995 | pass    |
| all runs    | NA18507  | HG00403     | 2%    | 0                  |            154 |  147 |         0.1469 | 0.1431 |         -0.0037 |   0.975 | pass    |
| all runs    | NA18507  | HG00403     | 2%    | α                  |            154 |  148 |         0.1469 | 0.1471 |         +0.0002 |   0.991 | pass    |
| all runs    | NA18507  | HG00403     | 5%    | 0                  |            154 |  113 |         0.1469 | 0.1250 |         -0.0219 |   0.850 | fail    |
| all runs    | NA18507  | HG00403     | 5%    | α                  |            154 |  144 |         0.1469 | 0.1475 |         +0.0006 |   0.986 | pass    |
| all runs    | NA18507  | HG00403     | 10%   | 0                  |            154 |   34 |         0.1469 | 0.0299 |         -0.1170 |   0.203 | fail    |
| all runs    | NA18507  | HG00403     | 10%   | α                  |            154 |  138 |         0.1469 | 0.1487 |         +0.0018 |   0.953 | pass    |
| all runs    | NA18507  | HG00188     | 1%    | 0                  |            154 |  154 |         0.1469 | 0.1469 |         +0.0000 |   1.000 | pass    |
| all runs    | NA18507  | HG00188     | 1%    | α                  |            154 |  150 |         0.1469 | 0.1476 |         +0.0008 |   0.995 | pass    |
| all runs    | NA18507  | HG00188     | 2%    | 0                  |            154 |  152 |         0.1469 | 0.1454 |         -0.0015 |   0.990 | pass    |
| all runs    | NA18507  | HG00188     | 2%    | α                  |            154 |  150 |         0.1469 | 0.1476 |         +0.0007 |   0.995 | pass    |
| all runs    | NA18507  | HG00188     | 5%    | 0                  |            154 |  121 |         0.1469 | 0.1205 |         -0.0264 |   0.820 | fail    |
| all runs    | NA18507  | HG00188     | 5%    | α                  |            154 |  145 |         0.1469 | 0.1481 |         +0.0012 |   0.984 | pass    |
| all runs    | NA18507  | HG00188     | 10%   | 0                  |            154 |   36 |         0.1469 | 0.0285 |         -0.1183 |   0.192 | fail    |
| all runs    | NA18507  | HG00188     | 10%   | α                  |            154 |  132 |         0.1469 | 0.1474 |         +0.0005 |   0.939 | pass    |
| all runs    | HG00403  | NA18507     | 1%    | 0                  |            233 |  233 |         0.3030 | 0.3030 |         -0.0000 |   0.999 | pass    |
| all runs    | HG00403  | NA18507     | 1%    | α                  |            233 |  228 |         0.3030 | 0.3047 |         +0.0017 |   0.994 | pass    |
| all runs    | HG00403  | NA18507     | 2%    | 0                  |            233 |  226 |         0.3030 | 0.2975 |         -0.0055 |   0.980 | pass    |
| all runs    | HG00403  | NA18507     | 2%    | α                  |            233 |  227 |         0.3030 | 0.3058 |         +0.0028 |   0.990 | pass    |
| all runs    | HG00403  | NA18507     | 5%    | 0                  |            233 |  169 |         0.3030 | 0.2203 |         -0.0826 |   0.727 | fail    |
| all runs    | HG00403  | NA18507     | 5%    | α                  |            233 |  219 |         0.3030 | 0.3066 |         +0.0036 |   0.985 | pass    |
| all runs    | HG00403  | NA18507     | 10%   | 0                  |            233 |   42 |         0.3030 | 0.0337 |         -0.2692 |   0.111 | fail    |
| all runs    | HG00403  | NA18507     | 10%   | α                  |            233 |  200 |         0.3030 | 0.2946 |         -0.0084 |   0.929 | pass    |
| all runs    | HG00403  | HG00188     | 1%    | 0                  |            233 |  232 |         0.3030 | 0.3021 |         -0.0008 |   0.997 | pass    |
| all runs    | HG00403  | HG00188     | 1%    | α                  |            233 |  228 |         0.3030 | 0.3048 |         +0.0018 |   0.994 | pass    |
| all runs    | HG00403  | HG00188     | 2%    | 0                  |            233 |  230 |         0.3030 | 0.3008 |         -0.0022 |   0.991 | pass    |
| all runs    | HG00403  | HG00188     | 2%    | α                  |            233 |  225 |         0.3030 | 0.3042 |         +0.0013 |   0.991 | pass    |
| all runs    | HG00403  | HG00188     | 5%    | 0                  |            233 |  201 |         0.3030 | 0.2582 |         -0.0447 |   0.852 | fail    |
| all runs    | HG00403  | HG00188     | 5%    | α                  |            233 |  222 |         0.3030 | 0.3080 |         +0.0050 |   0.983 | pass    |
| all runs    | HG00403  | HG00188     | 10%   | 0                  |            233 |  101 |         0.3030 | 0.1219 |         -0.1811 |   0.402 | fail    |
| all runs    | HG00403  | HG00188     | 10%   | α                  |            233 |  201 |         0.3030 | 0.3050 |         +0.0020 |   0.951 | pass    |
| all runs    | HG00188  | NA18507     | 1%    | 0                  |            215 |  216 |         0.2822 | 0.2817 |         -0.0006 |   0.998 | pass    |
| all runs    | HG00188  | NA18507     | 1%    | α                  |            215 |  212 |         0.2822 | 0.2840 |         +0.0018 |   0.994 | pass    |
| all runs    | HG00188  | NA18507     | 2%    | 0                  |            215 |  211 |         0.2822 | 0.2802 |         -0.0021 |   0.993 | pass    |
| all runs    | HG00188  | NA18507     | 2%    | α                  |            215 |  207 |         0.2822 | 0.2842 |         +0.0020 |   0.987 | pass    |
| all runs    | HG00188  | NA18507     | 5%    | 0                  |            215 |  162 |         0.2822 | 0.2181 |         -0.0641 |   0.773 | fail    |
| all runs    | HG00188  | NA18507     | 5%    | α                  |            215 |  200 |         0.2822 | 0.2828 |         +0.0005 |   0.972 | pass    |
| all runs    | HG00188  | NA18507     | 10%   | 0                  |            215 |   44 |         0.2822 | 0.0521 |         -0.2302 |   0.185 | fail    |
| all runs    | HG00188  | NA18507     | 10%   | α                  |            215 |  181 |         0.2822 | 0.2736 |         -0.0087 |   0.919 | pass    |
| all runs    | HG00188  | HG00403     | 1%    | 0                  |            215 |  214 |         0.2822 | 0.2822 |         -0.0000 |   0.998 | pass    |
| all runs    | HG00188  | HG00403     | 1%    | α                  |            215 |  211 |         0.2822 | 0.2849 |         +0.0027 |   0.990 | pass    |
| all runs    | HG00188  | HG00403     | 2%    | 0                  |            215 |  215 |         0.2822 | 0.2810 |         -0.0013 |   0.995 | pass    |
| all runs    | HG00188  | HG00403     | 2%    | α                  |            215 |  208 |         0.2822 | 0.2854 |         +0.0032 |   0.989 | pass    |
| all runs    | HG00188  | HG00403     | 5%    | 0                  |            215 |  181 |         0.2822 | 0.2380 |         -0.0442 |   0.843 | fail    |
| all runs    | HG00188  | HG00403     | 5%    | α                  |            215 |  204 |         0.2822 | 0.2847 |         +0.0025 |   0.984 | pass    |
| all runs    | HG00188  | HG00403     | 10%   | 0                  |            215 |   87 |         0.2822 | 0.1026 |         -0.1797 |   0.363 | fail    |
| all runs    | HG00188  | HG00403     | 10%   | α                  |            215 |  189 |         0.2822 | 0.2813 |         -0.0010 |   0.939 | pass    |
| runs ≥ 1 Mb | NA18507  | HG00403     | 1%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00403     | 1%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00403     | 2%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00403     | 2%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00403     | 5%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00403     | 5%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00403     | 10%   | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00403     | 10%   | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00188     | 1%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00188     | 1%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00188     | 2%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00188     | 2%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00188     | 5%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00188     | 5%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00188     | 10%   | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | NA18507  | HG00188     | 10%   | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | NA18507     | 1%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | NA18507     | 1%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | NA18507     | 2%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | NA18507     | 2%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | NA18507     | 5%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | NA18507     | 5%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | NA18507     | 10%   | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | NA18507     | 10%   | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | HG00188     | 1%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | HG00188     | 1%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | HG00188     | 2%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | HG00188     | 2%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | HG00188     | 5%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | HG00188     | 5%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | HG00188     | 10%   | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00403  | HG00188     | 10%   | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | NA18507     | 1%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | NA18507     | 1%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | NA18507     | 2%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | NA18507     | 2%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | NA18507     | 5%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | NA18507     | 5%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | NA18507     | 10%   | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | NA18507     | 10%   | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | HG00403     | 1%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | HG00403     | 1%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | HG00403     | 2%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | HG00403     | 2%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | HG00403     | 5%    | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | HG00403     | 5%    | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | HG00403     | 10%   | 0                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |
| runs ≥ 1 Mb | HG00188  | HG00403     | 10%   | α                  |              0 |    0 |         0.0000 | 0.0000 |         +0.0000 |    none | pass    |

## Limitations

- One chromosome and three samples. The rows above are observations, not an
  estimate of a rate.

- The two paths of (a) share their sequencing data, sites and frequencies.
  Agreement shows that the read model and the genotype model decode the same
  evidence alike. It does not show that either set of runs is true.

- 2)  thins and adds counts. It does not reproduce what a mixture of reads
      would change upstream: alignment, duplicate marking, base and mapping quality
      filters, or the correlation of reads within a fragment.

- The contamination fraction is supplied to the decode. Nothing here estimates
  it.

- The reference of (b) is the receiver’s own decode, so (b) measures how much
  a mixture moves the runs, not whether the runs are true.
