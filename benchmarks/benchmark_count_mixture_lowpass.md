Low-pass conditional mixture moments: chr20 calibration and scaling
================

- [Declaration](#declaration)
  - [Frozen criteria](#frozen-criteria)
- [Execution and provenance](#execution-and-provenance)
- [Criteria and numerical
  invariants](#criteria-and-numerical-invariants)
- [Linked chr20 recovery](#linked-chr20-recovery)
- [Registered count thinning](#registered-count-thinning)
- [Scaling and memory](#scaling-and-memory)
  - [Exploratory aggregate
    comparison](#exploratory-aggregate-comparison)
- [Retained evidence and scope](#retained-evidence-and-scope)

## Declaration

The executable declaration is `count_mixture_lowpass/declaration.R`. It
is committed before the campaign it judges. This is a **chr20
point-estimation and interval-calibration check**, not whole-genome
release evidence.

The registered phased 1000G chr20 source, founder pedigree and matched
sites supply linked haplotypes. Sixty CEU and sixty YRI founders are
drawn without replacement with seed 332202610. Each replicate draws
distinct founders: two CEU haplotypes form the sample, two other CEU
haplotypes form a same-population second genome, or two YRI haplotypes
form an other-population second genome. An IBD1 second genome transmits
one complete phased haplotype of the sample and takes its other
haplotype from a distinct CEU founder. Frequencies are the called phased
allele frequencies among all founders of the corresponding population,
clamped to \[0.001, 0.999\]. Founder labels exclude pedigree
descendants; residual relatedness is not independently screened. Linkage
is retained across sites; no genotype-independent SNP simulation
substitutes for phased genomes.

There are 50 replicates in each of 108 primary strata: three
relationships, shares 0/0.01/0.02/0.05/0.10/0.20, depths 0.1/0.2/0.4,
and fixed errors 0.0001/0.0002. Depth is Poisson. A fragment
independently chooses a source, then one of its alleles, then has
symmetric sequencing error. Equivalently, ALT is binomial conditional on
the two phased dosages. With probability 0.002 a site’s ALT count is
uniform on 0..depth. Excess homozygosity is zero. The fitting
calibration equals these generating parameters. Matched IBD1 and
unrelated hypotheses remain separate conditional fits, not relationship
classification. Blocks are **1 Mb on chr20**, to reach the unchanged
30-informative-block interval floor; no 10 Mb whole-genome precision
claim is made. All statuses, missing intervals and failed fits remain
denominators.

The three registered 30x chr20 count tables are thinned at 1/300, 1/150
and 1/75, both as single sources and as ordered two-source combinations
at each nonzero share. Ref and ALT counts are independently
binomial-thinned, not reads or templates. The nominal share is the
thinning probability share, not necessarily the realized fragment share
when source depths differ. Each stratum has 50 draws, fixed error
0.00015, artefact weight 0.002 and zero excess homozygosity. These are
calibration-transfer diagnostics, not known-truth accuracy gates.
Single-source baseline offsets are retained.

### Frozen criteria

Release targets from the design require at least 200 independent
replicates per primary stratum on a genome-wide linked reference. With
50 chr20 draws these targets are **not evaluated as release gates**,
even where a numerical threshold is met. Observed thresholds are shown
separately:

- For correctly specified interior shares at least 0.05: absolute mean
  bias \<= 0.005/0.0035/0.0025 at 0.1/0.2/0.4x.
- At least 90% have absolute error \<= 0.020/0.014/0.010 for anchored
  other-population singles, \<= 0.025/0.015/0.008 for matched unrelated,
  \<= 0.050/0.025/0.015 for matched IBD1. Missing fits count as
  failures.
- Interior 95% interval coverage is 90–98%, counting missing intervals
  as uncovered. Zero-share positive claims use lower limit \> 0; their
  one-sided exact 95% binomial upper bound must be \<= 0.075. Missing
  intervals also have a separately reported denominator, not an inferred
  null claim.
- At least one replicate in every synthetic and thinning stratum has all
  native-input block moment sums independently reconstructed in R, with
  exact integer counters and relative error \<= 1e-10 (absolute 1e-10
  near zero). Every primary stratum’s first replicate is also fitted at
  one and four threads: statuses must agree and shares and moments must
  agree to the same tolerance.
- Fresh processes, three repetitions each, scale sites 100k/200k/400k
  and sample-by-candidate groups 1/2/4. Peak RSS and DuckDB peak buffers
  must be \<= 512 MiB overhead + 3 times modeled live state: 768 bytes
  per SQL input row plus 1024 bytes per physical native block and group.
  This explicitly includes SQL validation/join/sort state, not just
  bounded native moments. Memory limit is 3 GB, no spill. Timing
  exponent \<= 1.25 at each doubling is evaluated only if the 1x query
  takes at least five seconds; smaller workloads report timing without a
  verdict. Input rows, output rows, elapsed time, spreads, buffers, RSS
  and spill remain visible.

The exploratory ordered/unordered/materialization comparison is reported
separately from the declared campaign. Same-population relationship
ambiguity, weak pairs, frequency mismatch and calibration transfer are
limitations, not reasons to drop observations.

## Execution and provenance

The campaign uses registered inputs only. From the repository root, with
a private library containing the current `duckhtsbench` registry:

``` sh
taskset -c 0-5 Rscript benchmarks/count_mixture_lowpass/stage.R build/release/duckhts.duckdb_extension
taskset -c 0-5 Rscript benchmarks/count_mixture_lowpass/run.R build/release/duckhts.duckdb_extension /tmp/lowpass-campaign
```

`collect.R` reads the campaign, retained compensated-sum campaign and
exploratory memory probes, then writes the checked-in results. Recovery
workers run in three fresh R/DuckDB processes at a time, with one DuckDB
thread except the explicit four-thread comparison. Recovery timings are
not performance evidence. Scaling runs alone in fresh processes, one
DuckDB thread, three repetitions, six assigned CPU cores (0–5). Inputs
are warm cached; query elapsed time excludes startup, extension loading
and view creation. RSS includes those costs. DuckDB 3 GB memory limit
and no spill apply to the measured scaling queries. Source acquisition
and haplotype staging are outside measured queries.

Campaign source revision: `5eb4dc079d5756caab436d902d5779fd079826a5`.

R version 4.6.0 (2026-04-24); DuckDB R 1.5.5.

CPU(s): 20 Model name: 13th Gen Intel(R) Core(TM) i5-13500 BIOS Model
name: 13th Gen Intel(R) Core(TM) i5-13500 To Be Filled By O.E.M. CPU @
2.4GHz Socket(s): 1 NUMA node0 CPU(s): 0-19

Native source and driver SHA-256 hashes and the loaded extension hash
are retained in
[metadata.txt](count_mixture_lowpass/results/metadata.txt). The
declaration was committed as `60ce83730fa16dcb5bcc89f1ecfc6174ca95480c`
before the campaign.

Synthetic: 5400 replicates / 108 strata; count thinning: 4950 replicates
/ 99 strata.

The phased source is the 1000G 3202-sample GRCh38 release of 2020-10-28,
chr20. Raw source and founder pedigree identities are pinned in the
registry. There are 108,757 matched sites and 120 selected founders in
the derived artifact `count_mixture_lowpass_chr20_haplotypes`. The three
count tables each contain 87,867 common chr20 sites; their canonical
count hashes are pinned in the registry. The network-free staging test
checks temporary cache isolation, phased decoding, population
frequencies, repeat staging, and failure on a missing raw input.
Pedigree and phase are checked; cryptic relatedness and
individual-frequency agreement are not established.

## Criteria and numerical invariants

| criterion                                     | denominator | verdict                                    |
|:----------------------------------------------|------------:|:-------------------------------------------|
| R moment oracle                               |         207 | met                                        |
| Exact integer counters                        |         207 | met                                        |
| 1/4-thread statuses                           |         207 | met                                        |
| 1/4-thread shares                             |         207 | met                                        |
| 1/4-thread moments                            |         207 | met                                        |
| Scaling RSS ceiling                           |          18 | met                                        |
| Scaling DuckDB buffer ceiling                 |          18 | met                                        |
| No scaling spill                              |          18 | met                                        |
| Timing exponent \<= 1.25                      |           2 | not evaluated: 1x \< 5 seconds             |
| Genome-wide release recovery/interval targets |         108 | not evaluated: chr20, 50 \< 200 replicates |

The R oracle forms weighted cross-products and ALT/discordance moments
from the generated counts and frequencies, independently of the macro’s
SQL expressions. It compares every block and every moment supplied to
the native fitter. Integer QC counters compare exactly. Native
integer-limb sum tests independently cover dyadic integer oracles,
subnormals, merges, round-to-even ties, sticky rounding, invalid inputs
and capacity errors.

| datasets | max_R_relative_moment_error | integer_counter_failures | status_failures | max_thread_share_error | max_thread_moment_error |
|---------:|----------------------------:|-------------------------:|----------------:|-----------------------:|------------------------:|
|      207 |                           0 |                        0 |               0 |                      0 |                       0 |

## Linked chr20 recovery

The primary point below is the conditional whole-fit `share`, not a
relationship selected by its score. Every stratum has 50 independent
random draws conditional on this fixed panel. Opposite conditional
models are also retained in the compressed replicate table.
`conditional` counts `ok_conditional` statuses; other statuses remain
denominators. Mean bias uses all finite reported shares (its count is
retained in `strata.csv`). A missing or non-conditional fit is a
point-accuracy failure, and a missing interval is uncovered. Wide or
disconnected confidence sets are not removed. `singles_mean`,
`singles_bias` and `singles_within_limit_fraction` retain the separate
anchored read-only diagnostic for other-population strata; they do not
replace the declared whole-fit comparator.

| criterion | strata | met | failed | incomplete | release_verdict                 |
|:----------|-------:|----:|-------:|-----------:|:--------------------------------|
| bias      |     54 |  13 |     41 |          0 | not evaluated: chr20, 50 \< 200 |
| accuracy  |     90 |   1 |     89 |          0 | not evaluated: chr20, 50 \< 200 |
| coverage  |     90 |  73 |     17 |          0 | not evaluated: chr20, 50 \< 200 |
| null      |     18 |   5 |     11 |          2 | not evaluated: chr20, 50 \< 200 |

These numerical thresholds mostly **fail for point accuracy** on chr20.
Adequate-looking coverage can reflect wide intervals, not accurate
shares. At 0.1x this evidence does not establish dependable 1–2% share
detection. No universal detection limit is inferred. In particular, the
matched IBD1 branch can have uninformative confidence sets and positive
null estimates. The 200-replicate genome-wide criteria are not passed by
pooling depths, errors, relationships or small-share strata. The
whole-genome reference has more sites and linkage blocks; chr20 results
neither certify nor refute its numerical targets.

|  id | relationship | depth | seq_error | truth | conditional | mean_share |    bias | accurate_fraction | interval_coverage | positive_claims |
|----:|:-------------|------:|----------:|------:|------------:|-----------:|--------:|------------------:|------------------:|----------------:|
|   1 | same         |   0.1 |     1e-04 |  0.00 |          50 |     0.0268 |  0.0268 |              0.62 |              0.98 |               1 |
|   2 | other        |   0.1 |     1e-04 |  0.00 |          50 |     0.0080 |  0.0080 |              0.82 |              1.00 |               0 |
|   3 | ibd1         |   0.1 |     1e-04 |  0.00 |          37 |     0.0430 |  0.0430 |              0.70 |              0.98 |               1 |
|   4 | same         |   0.1 |     1e-04 |  0.01 |          49 |     0.0388 |  0.0288 |              0.64 |              0.96 |               2 |
|   5 | other        |   0.1 |     1e-04 |  0.01 |          50 |     0.0110 |  0.0010 |              0.92 |              1.00 |               1 |
|   6 | ibd1         |   0.1 |     1e-04 |  0.01 |          25 |     0.0881 |  0.0781 |              0.46 |              0.96 |               2 |
|   7 | same         |   0.1 |     1e-04 |  0.02 |          49 |     0.0301 |  0.0101 |              0.72 |              1.00 |               1 |
|   8 | other        |   0.1 |     1e-04 |  0.02 |          50 |     0.0252 |  0.0052 |              0.72 |              0.94 |               5 |
|   9 | ibd1         |   0.1 |     1e-04 |  0.02 |          35 |     0.0598 |  0.0398 |              0.66 |              0.96 |               2 |
|  10 | same         |   0.1 |     1e-04 |  0.05 |          50 |     0.0447 | -0.0053 |              0.26 |              0.94 |               3 |
|  11 | other        |   0.1 |     1e-04 |  0.05 |          50 |     0.0412 | -0.0088 |              0.46 |              0.96 |              13 |
|  12 | ibd1         |   0.1 |     1e-04 |  0.05 |          27 |     0.0913 |  0.0413 |              0.50 |              1.00 |               2 |
|  13 | same         |   0.1 |     1e-04 |  0.10 |          49 |     0.1024 |  0.0024 |              0.26 |              0.96 |              14 |
|  14 | other        |   0.1 |     1e-04 |  0.10 |          50 |     0.0961 | -0.0039 |              0.42 |              0.98 |              40 |
|  15 | ibd1         |   0.1 |     1e-04 |  0.10 |          16 |     0.1395 |  0.0395 |              0.00 |              0.92 |               4 |
|  16 | same         |   0.1 |     1e-04 |  0.20 |          49 |     0.2288 |  0.0288 |              0.18 |              0.90 |              38 |
|  17 | other        |   0.1 |     1e-04 |  0.20 |          50 |     0.1946 | -0.0054 |              0.36 |              0.98 |              50 |
|  18 | ibd1         |   0.1 |     1e-04 |  0.20 |          24 |     0.2505 |  0.0505 |              0.00 |              0.78 |              16 |
|  19 | same         |   0.2 |     1e-04 |  0.00 |          50 |     0.0102 |  0.0102 |              0.72 |              1.00 |               0 |
|  20 | other        |   0.2 |     1e-04 |  0.00 |          50 |     0.0064 |  0.0064 |              0.82 |              0.98 |               1 |
|  21 | ibd1         |   0.2 |     1e-04 |  0.00 |          47 |     0.0267 |  0.0267 |              0.64 |              0.98 |               0 |
|  22 | same         |   0.2 |     1e-04 |  0.01 |          50 |     0.0170 |  0.0070 |              0.70 |              1.00 |               1 |
|  23 | other        |   0.2 |     1e-04 |  0.01 |          50 |     0.0122 |  0.0022 |              0.80 |              0.94 |               5 |
|  24 | ibd1         |   0.2 |     1e-04 |  0.01 |          49 |     0.0443 |  0.0343 |              0.50 |              0.96 |               2 |
|  25 | same         |   0.2 |     1e-04 |  0.02 |          50 |     0.0223 |  0.0023 |              0.24 |              0.94 |               4 |
|  26 | other        |   0.2 |     1e-04 |  0.02 |          49 |     0.0190 | -0.0010 |              0.58 |              0.94 |               5 |
|  27 | ibd1         |   0.2 |     1e-04 |  0.02 |          50 |     0.0380 |  0.0180 |              0.64 |              1.00 |               2 |
|  28 | same         |   0.2 |     1e-04 |  0.05 |          50 |     0.0449 | -0.0051 |              0.44 |              0.98 |               8 |
|  29 | other        |   0.2 |     1e-04 |  0.05 |          50 |     0.0444 | -0.0056 |              0.46 |              0.98 |              22 |
|  30 | ibd1         |   0.2 |     1e-04 |  0.05 |          48 |     0.0718 |  0.0218 |              0.20 |              0.92 |              10 |
|  31 | same         |   0.2 |     1e-04 |  0.10 |          50 |     0.1019 |  0.0019 |              0.34 |              0.94 |              39 |
|  32 | other        |   0.2 |     1e-04 |  0.10 |          50 |     0.1016 |  0.0016 |              0.40 |              0.92 |              50 |
|  33 | ibd1         |   0.2 |     1e-04 |  0.10 |          44 |     0.1228 |  0.0228 |              0.18 |              0.96 |              17 |
|  34 | same         |   0.2 |     1e-04 |  0.20 |          50 |     0.2126 |  0.0126 |              0.26 |              0.98 |              50 |
|  35 | other        |   0.2 |     1e-04 |  0.20 |          50 |     0.1961 | -0.0039 |              0.40 |              0.98 |              50 |
|  36 | ibd1         |   0.2 |     1e-04 |  0.20 |          45 |     0.2441 |  0.0441 |              0.12 |              0.90 |              35 |
|  37 | same         |   0.4 |     1e-04 |  0.00 |          50 |     0.0089 |  0.0089 |              0.62 |              0.98 |               1 |
|  38 | other        |   0.4 |     1e-04 |  0.00 |          49 |     0.0055 |  0.0055 |              0.74 |              0.98 |               0 |
|  39 | ibd1         |   0.4 |     1e-04 |  0.00 |          50 |     0.0221 |  0.0221 |              0.60 |              0.94 |               3 |
|  40 | same         |   0.4 |     1e-04 |  0.01 |          50 |     0.0163 |  0.0063 |              0.32 |              0.94 |               6 |
|  41 | other        |   0.4 |     1e-04 |  0.01 |          50 |     0.0123 |  0.0023 |              0.72 |              0.92 |               6 |
|  42 | ibd1         |   0.4 |     1e-04 |  0.01 |          50 |     0.0249 |  0.0149 |              0.56 |              1.00 |               0 |
|  43 | same         |   0.4 |     1e-04 |  0.02 |          50 |     0.0219 |  0.0019 |              0.28 |              0.92 |               7 |
|  44 | other        |   0.4 |     1e-04 |  0.02 |          50 |     0.0234 |  0.0034 |              0.50 |              0.94 |              16 |
|  45 | ibd1         |   0.4 |     1e-04 |  0.02 |          49 |     0.0370 |  0.0170 |              0.26 |              0.98 |               4 |
|  46 | same         |   0.4 |     1e-04 |  0.05 |          50 |     0.0503 |  0.0003 |              0.18 |              0.88 |              31 |
|  47 | other        |   0.4 |     1e-04 |  0.05 |          50 |     0.0516 |  0.0016 |              0.54 |              0.94 |              48 |
|  48 | ibd1         |   0.4 |     1e-04 |  0.05 |          50 |     0.0499 | -0.0001 |              0.26 |              1.00 |               7 |
|  49 | same         |   0.4 |     1e-04 |  0.10 |          50 |     0.0951 | -0.0049 |              0.22 |              0.96 |              49 |
|  50 | other        |   0.4 |     1e-04 |  0.10 |          50 |     0.1015 |  0.0015 |              0.44 |              0.92 |              50 |
|  51 | ibd1         |   0.4 |     1e-04 |  0.10 |          50 |     0.0994 | -0.0006 |              0.22 |              0.98 |              28 |
|  52 | same         |   0.4 |     1e-04 |  0.20 |          50 |     0.2052 |  0.0052 |              0.22 |              0.96 |              50 |
|  53 | other        |   0.4 |     1e-04 |  0.20 |          50 |     0.2046 |  0.0046 |              0.34 |              0.98 |              50 |
|  54 | ibd1         |   0.4 |     1e-04 |  0.20 |          50 |     0.2251 |  0.0251 |              0.20 |              0.96 |              50 |
|  55 | same         |   0.1 |     2e-04 |  0.00 |          49 |     0.0258 |  0.0258 |              0.66 |              0.98 |               1 |
|  56 | other        |   0.1 |     2e-04 |  0.00 |          50 |     0.0072 |  0.0072 |              0.86 |              1.00 |               0 |
|  57 | ibd1         |   0.1 |     2e-04 |  0.00 |          34 |     0.0494 |  0.0494 |              0.66 |              0.98 |               1 |
|  58 | same         |   0.1 |     2e-04 |  0.01 |          48 |     0.0297 |  0.0197 |              0.70 |              0.94 |               2 |
|  59 | other        |   0.1 |     2e-04 |  0.01 |          50 |     0.0138 |  0.0038 |              0.80 |              1.00 |               1 |
|  60 | ibd1         |   0.1 |     2e-04 |  0.01 |          39 |     0.0389 |  0.0289 |              0.74 |              0.98 |               1 |
|  61 | same         |   0.1 |     2e-04 |  0.02 |          47 |     0.0357 |  0.0157 |              0.68 |              1.00 |               2 |
|  62 | other        |   0.1 |     2e-04 |  0.02 |          50 |     0.0221 |  0.0021 |              0.76 |              0.94 |               5 |
|  63 | ibd1         |   0.1 |     2e-04 |  0.02 |          29 |     0.0792 |  0.0592 |              0.56 |              0.98 |               1 |
|  64 | same         |   0.1 |     2e-04 |  0.05 |          49 |     0.0597 |  0.0097 |              0.30 |              0.96 |               6 |
|  65 | other        |   0.1 |     2e-04 |  0.05 |          50 |     0.0423 | -0.0077 |              0.48 |              0.90 |              12 |
|  66 | ibd1         |   0.1 |     2e-04 |  0.05 |          32 |     0.0701 |  0.0201 |              0.60 |              0.96 |               2 |
|  67 | same         |   0.1 |     2e-04 |  0.10 |          47 |     0.1124 |  0.0124 |              0.12 |              0.88 |              21 |
|  68 | other        |   0.1 |     2e-04 |  0.10 |          50 |     0.1004 |  0.0004 |              0.58 |              0.94 |              43 |
|  69 | ibd1         |   0.1 |     2e-04 |  0.10 |          24 |     0.1357 |  0.0357 |              0.00 |              0.94 |               4 |
|  70 | same         |   0.1 |     2e-04 |  0.20 |          47 |     0.1900 | -0.0100 |              0.40 |              0.96 |              37 |
|  71 | other        |   0.1 |     2e-04 |  0.20 |          50 |     0.1975 | -0.0025 |              0.38 |              0.96 |              50 |
|  72 | ibd1         |   0.1 |     2e-04 |  0.20 |          15 |     0.1974 | -0.0026 |              0.00 |              0.92 |               8 |
|  73 | same         |   0.2 |     2e-04 |  0.00 |          50 |     0.0137 |  0.0137 |              0.72 |              0.94 |               3 |
|  74 | other        |   0.2 |     2e-04 |  0.00 |          50 |     0.0056 |  0.0056 |              0.86 |              1.00 |               0 |
|  75 | ibd1         |   0.2 |     2e-04 |  0.00 |          50 |     0.0328 |  0.0328 |              0.64 |              0.92 |               4 |
|  76 | same         |   0.2 |     2e-04 |  0.01 |          50 |     0.0223 |  0.0123 |              0.64 |              0.98 |               4 |
|  77 | other        |   0.2 |     2e-04 |  0.01 |          50 |     0.0087 | -0.0013 |              0.88 |              0.94 |               2 |
|  78 | ibd1         |   0.2 |     2e-04 |  0.01 |          46 |     0.0438 |  0.0338 |              0.62 |              0.94 |               5 |
|  79 | same         |   0.2 |     2e-04 |  0.02 |          50 |     0.0251 |  0.0051 |              0.54 |              0.98 |               3 |
|  80 | other        |   0.2 |     2e-04 |  0.02 |          50 |     0.0213 |  0.0013 |              0.50 |              0.94 |               7 |
|  81 | ibd1         |   0.2 |     2e-04 |  0.02 |          49 |     0.0431 |  0.0231 |              0.62 |              1.00 |               2 |
|  82 | same         |   0.2 |     2e-04 |  0.05 |          50 |     0.0594 |  0.0094 |              0.36 |              0.94 |              21 |
|  83 | other        |   0.2 |     2e-04 |  0.05 |          50 |     0.0397 | -0.0103 |              0.40 |              0.90 |              24 |
|  84 | ibd1         |   0.2 |     2e-04 |  0.05 |          48 |     0.0657 |  0.0157 |              0.34 |              0.94 |               4 |
|  85 | same         |   0.2 |     2e-04 |  0.10 |          50 |     0.1041 |  0.0041 |              0.30 |              1.00 |              41 |
|  86 | other        |   0.2 |     2e-04 |  0.10 |          50 |     0.0991 | -0.0009 |              0.56 |              0.98 |              49 |
|  87 | ibd1         |   0.2 |     2e-04 |  0.10 |          42 |     0.1052 |  0.0052 |              0.22 |              0.90 |               9 |
|  88 | same         |   0.2 |     2e-04 |  0.20 |          50 |     0.2086 |  0.0086 |              0.10 |              0.94 |              50 |
|  89 | other        |   0.2 |     2e-04 |  0.20 |          50 |     0.1967 | -0.0033 |              0.52 |              0.92 |              50 |
|  90 | ibd1         |   0.2 |     2e-04 |  0.20 |          45 |     0.1972 | -0.0028 |              0.16 |              0.92 |              31 |
|  91 | same         |   0.4 |     2e-04 |  0.00 |          50 |     0.0087 |  0.0087 |              0.52 |              1.00 |               0 |
|  92 | other        |   0.4 |     2e-04 |  0.00 |          50 |     0.0074 |  0.0074 |              0.66 |              0.98 |               1 |
|  93 | ibd1         |   0.4 |     2e-04 |  0.00 |          49 |     0.0228 |  0.0228 |              0.56 |              0.92 |               3 |
|  94 | same         |   0.4 |     2e-04 |  0.01 |          50 |     0.0164 |  0.0064 |              0.32 |              0.98 |               2 |
|  95 | other        |   0.4 |     2e-04 |  0.01 |          49 |     0.0121 |  0.0021 |              0.74 |              0.88 |               6 |
|  96 | ibd1         |   0.4 |     2e-04 |  0.01 |          50 |     0.0210 |  0.0110 |              0.62 |              1.00 |               0 |
|  97 | same         |   0.4 |     2e-04 |  0.02 |          50 |     0.0235 |  0.0035 |              0.24 |              0.96 |              12 |
|  98 | other        |   0.4 |     2e-04 |  0.02 |          50 |     0.0213 |  0.0013 |              0.54 |              0.96 |              11 |
|  99 | ibd1         |   0.4 |     2e-04 |  0.02 |          50 |     0.0357 |  0.0157 |              0.18 |              0.96 |               5 |
| 100 | same         |   0.4 |     2e-04 |  0.05 |          50 |     0.0585 |  0.0085 |              0.20 |              0.94 |              35 |
| 101 | other        |   0.4 |     2e-04 |  0.05 |          50 |     0.0539 |  0.0039 |              0.56 |              0.96 |              49 |
| 102 | ibd1         |   0.4 |     2e-04 |  0.05 |          50 |     0.0650 |  0.0150 |              0.22 |              0.96 |              11 |
| 103 | same         |   0.4 |     2e-04 |  0.10 |          50 |     0.1077 |  0.0077 |              0.14 |              0.90 |              50 |
| 104 | other        |   0.4 |     2e-04 |  0.10 |          50 |     0.0994 | -0.0006 |              0.38 |              0.94 |              50 |
| 105 | ibd1         |   0.4 |     2e-04 |  0.10 |          50 |     0.0926 | -0.0074 |              0.18 |              1.00 |              23 |
| 106 | same         |   0.4 |     2e-04 |  0.20 |          50 |     0.2093 |  0.0093 |              0.30 |              0.96 |              50 |
| 107 | other        |   0.4 |     2e-04 |  0.20 |          50 |     0.2037 |  0.0037 |              0.46 |              0.94 |              50 |
| 108 | ibd1         |   0.4 |     2e-04 |  0.20 |          50 |     0.2142 |  0.0142 |              0.16 |              0.92 |              48 |

All per-stratum verdicts, missing-interval counts and exact one-sided
null bounds are in
[strata.csv](count_mixture_lowpass/results/strata.csv). The 50-draw
binomial uncertainty is material: even zero observed positive null
claims has a one-sided 95% upper bound of about 0.058. Missing intervals
make the null verdict incomplete, rather than certifying the null
criterion.

## Registered count thinning

|     | receiver | depth | replicates | mean_share | intervals | positive_claims |
|:----|:---------|------:|-----------:|-----------:|----------:|----------------:|
| 109 | NA18507  |   0.1 |         50 |    0.00030 |        49 |               0 |
| 110 | HG00403  |   0.1 |         50 |    0.00310 |        49 |               0 |
| 111 | HG00188  |   0.1 |         50 |    0.00141 |        46 |               0 |
| 112 | NA18507  |   0.2 |         50 |    0.00000 |        47 |               0 |
| 113 | HG00403  |   0.2 |         50 |    0.00000 |        49 |               0 |
| 114 | HG00188  |   0.2 |         50 |    0.00034 |        42 |               0 |
| 115 | NA18507  |   0.4 |         50 |    0.00000 |        27 |               0 |
| 116 | HG00403  |   0.4 |         50 |    0.00000 |        44 |               0 |
| 117 | HG00188  |   0.4 |         50 |    0.00000 |        40 |               0 |

Single-source offsets and missing intervals are retained. The largest
single-source mean share is about 0.0031 at nominal 0.1x. At higher
depth, several baseline fits concentrate at zero while intervals are
unavailable. These counts have empirical error and allele-specific
depth, not the synthetic generating calibration. The two-source share is
nominal and cannot be treated as exact truth without correcting for
source coverage. Thus two-source bias and coverage are diagnostics, not
release gate passes.

|     |  id | receiver | donor   | depth | truth | mean_share | intervals | positive_claims |
|:----|----:|:---------|:--------|------:|------:|-----------:|----------:|----------------:|
| 118 |  10 | HG00403  | NA18507 |   0.1 |  0.01 |     0.0051 |        49 |               0 |
| 119 |  11 | HG00188  | NA18507 |   0.1 |  0.01 |     0.0062 |        47 |               0 |
| 120 |  12 | NA18507  | HG00403 |   0.1 |  0.01 |     0.0003 |        49 |               0 |
| 121 |  13 | HG00188  | HG00403 |   0.1 |  0.01 |     0.0020 |        50 |               0 |
| 122 |  14 | NA18507  | HG00188 |   0.1 |  0.01 |     0.0011 |        49 |               0 |
| 123 |  15 | HG00403  | HG00188 |   0.1 |  0.01 |     0.0014 |        50 |               0 |
| 124 |  16 | HG00403  | NA18507 |   0.2 |  0.01 |     0.0002 |        49 |               0 |
| 125 |  17 | HG00188  | NA18507 |   0.2 |  0.01 |     0.0007 |        46 |               0 |
| 126 |  18 | NA18507  | HG00403 |   0.2 |  0.01 |     0.0000 |        48 |               0 |
| 127 |  19 | HG00188  | HG00403 |   0.2 |  0.01 |     0.0007 |        47 |               0 |
| 128 |  20 | NA18507  | HG00188 |   0.2 |  0.01 |     0.0001 |        48 |               0 |
| 129 |  21 | HG00403  | HG00188 |   0.2 |  0.01 |     0.0005 |        49 |               0 |
| 130 |  22 | HG00403  | NA18507 |   0.4 |  0.01 |     0.0003 |        47 |               0 |
| 131 |  23 | HG00188  | NA18507 |   0.4 |  0.01 |     0.0003 |        42 |               0 |
| 132 |  24 | NA18507  | HG00403 |   0.4 |  0.01 |     0.0000 |        40 |               0 |
| 133 |  25 | HG00188  | HG00403 |   0.4 |  0.01 |     0.0000 |        47 |               0 |
| 134 |  26 | NA18507  | HG00188 |   0.4 |  0.01 |     0.0000 |        44 |               0 |
| 135 |  27 | HG00403  | HG00188 |   0.4 |  0.01 |     0.0000 |        48 |               0 |
| 136 |  28 | HG00403  | NA18507 |   0.1 |  0.02 |     0.0100 |        50 |               0 |
| 137 |  29 | HG00188  | NA18507 |   0.1 |  0.02 |     0.0054 |        44 |               0 |
| 138 |  30 | NA18507  | HG00403 |   0.1 |  0.02 |     0.0007 |        50 |               0 |
| 139 |  31 | HG00188  | HG00403 |   0.1 |  0.02 |     0.0070 |        49 |               0 |
| 140 |  32 | NA18507  | HG00188 |   0.1 |  0.02 |     0.0024 |        50 |               0 |
| 141 |  33 | HG00403  | HG00188 |   0.1 |  0.02 |     0.0082 |        50 |               0 |
| 142 |  34 | HG00403  | NA18507 |   0.2 |  0.02 |     0.0013 |        50 |               0 |
| 143 |  35 | HG00188  | NA18507 |   0.2 |  0.02 |     0.0012 |        45 |               0 |
| 144 |  36 | NA18507  | HG00403 |   0.2 |  0.02 |     0.0001 |        50 |               0 |
| 145 |  37 | HG00188  | HG00403 |   0.2 |  0.02 |     0.0023 |        49 |               0 |
| 146 |  38 | NA18507  | HG00188 |   0.2 |  0.02 |     0.0003 |        49 |               0 |
| 147 |  39 | HG00403  | HG00188 |   0.2 |  0.02 |     0.0008 |        50 |               0 |
| 148 |  40 | HG00403  | NA18507 |   0.4 |  0.02 |     0.0002 |        48 |               0 |
| 149 |  41 | HG00188  | NA18507 |   0.4 |  0.02 |     0.0006 |        46 |               0 |
| 150 |  42 | NA18507  | HG00403 |   0.4 |  0.02 |     0.0000 |        45 |               0 |
| 151 |  43 | HG00188  | HG00403 |   0.4 |  0.02 |     0.0001 |        49 |               0 |
| 152 |  44 | NA18507  | HG00188 |   0.4 |  0.02 |     0.0000 |        48 |               0 |
| 153 |  45 | HG00403  | HG00188 |   0.4 |  0.02 |     0.0000 |        48 |               0 |
| 154 |  46 | HG00403  | NA18507 |   0.1 |  0.05 |     0.0163 |        49 |               0 |
| 155 |  47 | HG00188  | NA18507 |   0.1 |  0.05 |     0.0296 |        49 |               1 |
| 156 |  48 | NA18507  | HG00403 |   0.1 |  0.05 |     0.0082 |        50 |               0 |
| 157 |  49 | HG00188  | HG00403 |   0.1 |  0.05 |     0.0239 |        49 |               0 |
| 158 |  50 | NA18507  | HG00188 |   0.1 |  0.05 |     0.0069 |        50 |               0 |
| 159 |  51 | HG00403  | HG00188 |   0.1 |  0.05 |     0.0124 |        50 |               0 |
| 160 |  52 | HG00403  | NA18507 |   0.2 |  0.05 |     0.0170 |        50 |               0 |
| 161 |  53 | HG00188  | NA18507 |   0.2 |  0.05 |     0.0176 |        48 |               0 |
| 162 |  54 | NA18507  | HG00403 |   0.2 |  0.05 |     0.0041 |        50 |               0 |
| 163 |  55 | HG00188  | HG00403 |   0.2 |  0.05 |     0.0108 |        50 |               0 |
| 164 |  56 | NA18507  | HG00188 |   0.2 |  0.05 |     0.0021 |        50 |               0 |
| 165 |  57 | HG00403  | HG00188 |   0.2 |  0.05 |     0.0108 |        50 |               1 |
| 166 |  58 | HG00403  | NA18507 |   0.4 |  0.05 |     0.0146 |        50 |               2 |
| 167 |  59 | HG00188  | NA18507 |   0.4 |  0.05 |     0.0133 |        48 |               0 |
| 168 |  60 | NA18507  | HG00403 |   0.4 |  0.05 |     0.0006 |        50 |               0 |
| 169 |  61 | HG00188  | HG00403 |   0.4 |  0.05 |     0.0093 |        49 |               0 |
| 170 |  62 | NA18507  | HG00188 |   0.4 |  0.05 |     0.0013 |        50 |               0 |
| 171 |  63 | HG00403  | HG00188 |   0.4 |  0.05 |     0.0061 |        50 |               0 |
| 172 |  64 | HG00403  | NA18507 |   0.1 |  0.10 |     0.0686 |        50 |               9 |
| 173 |  65 | HG00188  | NA18507 |   0.1 |  0.10 |     0.0678 |        50 |               6 |
| 174 |  66 | NA18507  | HG00403 |   0.1 |  0.10 |     0.0572 |        50 |              10 |
| 175 |  67 | HG00188  | HG00403 |   0.1 |  0.10 |     0.0519 |        49 |               3 |
| 176 |  68 | NA18507  | HG00188 |   0.1 |  0.10 |     0.0381 |        50 |               1 |
| 177 |  69 | HG00403  | HG00188 |   0.1 |  0.10 |     0.0453 |        50 |               3 |
| 178 |  70 | HG00403  | NA18507 |   0.2 |  0.10 |     0.0674 |        50 |              32 |
| 179 |  71 | HG00188  | NA18507 |   0.2 |  0.10 |     0.0695 |        47 |              21 |
| 180 |  72 | NA18507  | HG00403 |   0.2 |  0.10 |     0.0369 |        50 |              13 |
| 181 |  73 | HG00188  | HG00403 |   0.2 |  0.10 |     0.0511 |        48 |              10 |
| 182 |  74 | NA18507  | HG00188 |   0.2 |  0.10 |     0.0390 |        50 |              12 |
| 183 |  75 | HG00403  | HG00188 |   0.2 |  0.10 |     0.0470 |        50 |              16 |
| 184 |  76 | HG00403  | NA18507 |   0.4 |  0.10 |     0.0667 |        50 |              47 |
| 185 |  77 | HG00188  | NA18507 |   0.4 |  0.10 |     0.0665 |        48 |              26 |
| 186 |  78 | NA18507  | HG00403 |   0.4 |  0.10 |     0.0406 |        50 |              35 |
| 187 |  79 | HG00188  | HG00403 |   0.4 |  0.10 |     0.0634 |        49 |              26 |
| 188 |  80 | NA18507  | HG00188 |   0.4 |  0.10 |     0.0301 |        50 |              21 |
| 189 |  81 | HG00403  | HG00188 |   0.4 |  0.10 |     0.0482 |        50 |              35 |
| 190 |  82 | HG00403  | NA18507 |   0.1 |  0.20 |     0.1978 |        50 |              44 |
| 191 |  83 | HG00188  | NA18507 |   0.1 |  0.20 |     0.1630 |        49 |              31 |
| 192 |  84 | NA18507  | HG00403 |   0.1 |  0.20 |     0.1465 |        50 |              31 |
| 193 |  85 | HG00188  | HG00403 |   0.1 |  0.20 |     0.1554 |        50 |              27 |
| 194 |  86 | NA18507  | HG00188 |   0.1 |  0.20 |     0.1086 |        50 |              20 |
| 195 |  87 | HG00403  | HG00188 |   0.1 |  0.20 |     0.1418 |        50 |              33 |
| 196 |  88 | HG00403  | NA18507 |   0.2 |  0.20 |     0.1856 |        50 |              50 |
| 197 |  89 | HG00188  | NA18507 |   0.2 |  0.20 |     0.1962 |        48 |              48 |
| 198 |  90 | NA18507  | HG00403 |   0.2 |  0.20 |     0.1463 |        50 |              50 |
| 199 |  91 | HG00188  | HG00403 |   0.2 |  0.20 |     0.1544 |        47 |              47 |
| 200 |  92 | NA18507  | HG00188 |   0.2 |  0.20 |     0.1501 |        50 |              50 |
| 201 |  93 | HG00403  | HG00188 |   0.2 |  0.20 |     0.1305 |        50 |              48 |
| 202 |  94 | HG00403  | NA18507 |   0.4 |  0.20 |     0.1845 |        50 |              50 |
| 203 |  95 | HG00188  | NA18507 |   0.4 |  0.20 |     0.1880 |        49 |              49 |
| 204 |  96 | NA18507  | HG00403 |   0.4 |  0.20 |     0.1474 |        50 |              50 |
| 205 |  97 | HG00188  | HG00403 |   0.4 |  0.20 |     0.1670 |        50 |              50 |
| 206 |  98 | NA18507  | HG00188 |   0.4 |  0.20 |     0.1437 |        50 |              50 |
| 207 |  99 | HG00403  | HG00188 |   0.4 |  0.20 |     0.1467 |        50 |              50 |

## Scaling and memory

Each input row has two frequencies, three counts, position, sample,
candidate and chromosome. Sites use a generated numeric panel. Group
scaling is (samples,candidates) = (1,1), (2,1), (2,2), with shared
counts validated across candidates. The ranges fit within one 1 Mb block
per group: this measures SQL input scaling and bounded aggregation,
**not** chr20 interval fitting or scaling in the number of jackknife
blocks. Output denominators are conditional rows, not raw sites. The
live-state model and ceiling were frozen before the run. The 512 MiB
overhead allowance is conservative relative to the measured
empty-process RSS:

| rss_MiB | buffer_bytes | seconds |
|--------:|-------------:|--------:|
| 121.883 |            0 |   0.000 |
| 122.031 |            0 |   0.001 |
| 121.871 |            0 |   0.000 |

| axis   | factor | input_rows | output_rows | median_seconds | min_seconds | max_seconds | max_rss_MiB | max_buffer_MiB | ceiling_MiB | spill_MiB | exponent | exponent_verdict               |
|:-------|-------:|-----------:|------------:|---------------:|------------:|------------:|------------:|---------------:|------------:|----------:|---------:|:-------------------------------|
| groups |      1 |      1e+05 |           1 |          0.171 |       0.170 |       0.175 |     179.715 |         90.547 |     731.729 |         0 |       NA | not evaluated: 1x \< 5 seconds |
| groups |      2 |      2e+05 |           2 |          0.329 |       0.327 |       0.333 |     252.559 |        205.031 |     951.459 |         0 |    0.944 | not evaluated: 1x \< 5 seconds |
| groups |      4 |      4e+05 |           4 |          0.546 |       0.543 |       0.549 |     272.301 |        237.344 |    1390.918 |         0 |    0.731 | not evaluated: 1x \< 5 seconds |
| sites  |      1 |      1e+05 |           1 |          0.172 |       0.171 |       0.174 |     180.020 |         90.605 |     731.729 |         0 |       NA | not evaluated: 1x \< 5 seconds |
| sites  |      2 |      2e+05 |           1 |          0.346 |       0.343 |       0.351 |     236.270 |        174.623 |     951.456 |         0 |    1.008 | not evaluated: 1x \< 5 seconds |
| sites  |      4 |      4e+05 |           1 |          0.686 |       0.682 |       0.688 |     348.301 |        338.234 |    1390.909 |         0 |    0.987 | not evaluated: 1x \< 5 seconds |

All 18 scaling runs satisfy RSS and DuckDB buffer ceilings without
spill. The 1x runs are shorter than five seconds, so the exponent gate
is **not evaluated**. The apparent exponents are descriptive, not a
passed timing criterion. Longer identical workloads and a multi-block
timing workload are still missing. Individual repetitions and modeled
live-state bytes are in
[scaling.csv](count_mixture_lowpass/results/scaling.csv).

### Exploratory aggregate comparison

| implementation | sites | threads | seconds | rss_MiB | buffer_MiB | spill_bytes |
|:---------------|------:|--------:|--------:|--------:|-----------:|------------:|
| materialized   | 1e+05 |       1 |   0.173 | 196.426 |    123.859 |           0 |
| materialized   | 2e+05 |       1 |   0.331 | 269.625 |    242.734 |           0 |
| materialized   | 4e+05 |       1 |   0.666 | 416.742 |    481.422 |           0 |
| ordered        | 1e+05 |       1 |   0.355 | 259.855 |    157.594 |           0 |
| ordered        | 2e+05 |       1 |   0.713 | 365.625 |    428.402 |           0 |
| ordered        | 4e+05 |       1 |   1.577 | 577.973 |    855.062 |           0 |
| unordered      | 1e+05 |       1 |   0.159 | 188.148 |    103.453 |           0 |
| unordered      | 2e+05 |       1 |   0.316 | 252.184 |    202.266 |           0 |
| unordered      | 4e+05 |       1 |   0.644 | 379.914 |    399.453 |           0 |

The ordered implementation was measured at `1c30c875`; its native
sources are unchanged in the current rebased counterpart `6e017d1c`.
These one-thread, single-run probes use the same numeric range workload
as scaling, but short identifier labels and one repetition rather than
three; they are exploratory comparisons, not a no-regression timing
claim. `unordered` is compensated `fsum` without raw materialization;
`materialized` is `fsum` with it. Explicit materialization increased
peak buffers and RSS in this comparison, so the macro uses
`NOT MATERIALIZED`. It still scans the relation repeatedly and SQL
memory remains input-sized.

Ordered aggregates retained rows: at 400k sites their peak buffers were
855 MiB and RSS 578 MiB. The final exact-sum range workload uses fixed
152-byte integer-limb accumulators, rounded once to double; final peak
buffers and RSS are in the three-repeat table above. Native moments do
not retain sites. This is **not** a constant-memory claim for
validation, joins, window ordering or the entire query.

Compensated sums pass the moment tolerance but did not preserve the
share criterion near a floating objective minimum. The retained
[compensated-sum
counterexamples](count_mixture_lowpass/results/compensated_sum_counterexamples.csv)
contain all 207 dataset checks, including 69 failed share comparisons
(maximum about 1.6e-7). They are not discarded to certify the final
sums. Exact sums preserve the fitter and the unchanged bitwise SQL
regression.

## Retained evidence and scope

- [replicates.csv.gz](count_mixture_lowpass/results/replicates.csv.gz):
  all conditional fits, statuses, primary-model flags and principal
  diagnostics.
- [replicate_qc.csv.gz](count_mixture_lowpass/results/replicate_qc.csv.gz):
  all replicate input/output counters, pair denominators and shared
  singles, pair-residual and leverage diagnostics. Join on
  `(kind,id,replicate)`.
- [oracle_and_threads.csv](count_mixture_lowpass/results/oracle_and_threads.csv):
  every independently checked stratum, not only passing examples.

The approximately 1.3 MB result bundle retains failures as well as
successes. This campaign checks calibrated chr20 simulations and real
count-thinning transfer. It does not establish a genome-wide precision
or detection limit, arbitrary calibration robustness, sample-specific
excess homozygosity, relationship classification, or error calibration
for arbitrary low-pass inputs. No criterion is widened to obtain a pass.
