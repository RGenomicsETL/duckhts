Genome-wide low-pass second-genome campaign
================

- [Declaration and scope](#declaration-and-scope)
  - [Frozen release targets](#frozen-release-targets)
- [Inputs, simulation and provenance](#inputs-simulation-and-provenance)
  - [Execution](#execution)
  - [Covered-site and microbatch
    checks](#covered-site-and-microbatch-checks)
- [Release outcomes](#release-outcomes)
  - [Detection and recovery](#detection-and-recovery)
  - [Intervals, statuses and ambiguous
    relationships](#intervals-statuses-and-ambiguous-relationships)
  - [Secondary depth](#secondary-depth)
- [Independent numerical checks](#independent-numerical-checks)
- [Performance and scaling](#performance-and-scaling)
- [Interpretation](#interpretation)

## Declaration and scope

The frozen declaration is
[`count_mixture_lowpass_genome/declaration.R`](count_mixture_lowpass_genome/declaration.R),
with its pre-run [execution
declaration](count_mixture_lowpass_genome/execution_declaration.md). It
commits 200 replicates in each of 90 strata: shares
0/0.02/0.05/0.10/0.20, same-superpopulation unrelated,
other-superpopulation unrelated and IBD1, at 0.1x, 0.2x and 0.4x, with
sequencing error 0.0001 and 0.0002. Seventy-five strata are primary:
0.1x/0.4x at error 0.0001, and all depths at error 0.0002. The 15
error-0.0001 strata at 0.2x are secondary. All 18,000 replicates use
site artefact probability 0.002 and fitted excess homozygosity zero.
Statistical criteria and the calibration grid are unchanged; an
incomplete stratum is not evaluated against release targets.

This report is separate from the [chr20
report](benchmark_count_mixture_lowpass.md); its sections and verdicts
are unchanged. These are phased-genome simulations, not real count-table
calibration-transfer measurements. Correct generating relationships
select which conditional fit is judged; the algorithm does not learn
that relationship.

### Frozen release targets

- At shares \>= 0.05, absolute mean bias \<= 0.005/0.0035/0.0025 at
  0.1/0.2/0.4x. At least 90% of replicates have absolute error \<=
  0.025/0.015/0.008 for matched unrelated, 0.020/0.014/0.010 for
  anchored other-population unrelated, or 0.050/0.025/0.015 for matched
  IBD1.
- Coverage \>= 92%, with the exact two-sided 95% binomial interval
  containing the nominal 95% coverage. Null positive claims use
  `share_lower > 0`; their exact one-sided 95% binomial upper bound must
  be \<= 0.075.
- Failed statuses, missing shares and missing/invalid intervals remain
  in denominators. A recovery pass also requires every fit to be
  `ok_conditional` with a finite share and interval. Strata below 200
  completions have no verdict.
- Interval width mean/median/90th percentile/maximum and full-range
  intervals are reported. Matched-frequency diagnostics must remain
  `unresolved`, including IBD1; relative scores do not choose a
  relationship.
- Replicate 1 of each stratum has an independent R moment oracle, exact
  integer counters, and relative moment error \<= 1e-10 (absolute near
  zero). One/four DuckDB threads must agree on statuses, missingness and
  diagnostics; share and interval differences must be \<= 1e-5, moment
  differences \<= 1e-10.
- No spill. RSS and reported DuckDB peak buffer each have the frozen
  ceiling 768 MiB + 3 times (512 bytes per input counts row + live R
  panels/counts
  - native 256 MiB cap). For the fresh-process scaling series, live R
    state is the largest chromosome panel/count pair. This is a modeled
    SQL preparation allowance, not a measurement of native live state.
    Native moments are bounded independently; joins, validation and
    ordered preparation dominate these SQL workloads. The DuckDB memory
    setting is 6 GB; measured peak buffer is reported separately, not
    assumed to equal that setting.
- Fixed-block site scaling uses all autosomes at strides 4/2/1
  (1x/2x/4x), in fresh processes. The end-to-end endpoint timing
  exponent must be \<= 1.25 if the smallest run takes \>= 5 seconds.
  Query-only timings are also retained.

## Inputs, simulation and provenance

Each autosome’s children BCF and matched AF Parquet is resolved and
SHA-256 verified through `r/duckhtsbench`. The chr20 AF identity is the
ordered rebuild pinned to `b862fc7d`. Each derived chromosome artifact
is registered as `count_mixture_lowpass_chr<N>_children_haplotypes`; its
source identities and settings are embedded, and the run manifest pins
the derived bytes again. No download is needed. The network-free staging
fixture covers allele order, all four diploid encodings, repeat staging,
family-disjoint donors, rejection of unphased GTs and checksum failure.
A full BCF audit checks all children for unphased GTs, independently of
the panel join.

Thirty-six CEU and 56 YRI children remain after selecting one child per
family. Every replicate draws one CEU child as the sample, a distinct
CEU child for same-superpopulation or IBD1, or a YRI child for
other-superpopulation. Different families supply the unrelated genomes;
cryptic relatedness is not independently screened. These are actual
phased diploid children, rather than the chr20 driver’s diploids
assembled from separate founders. Sample and donor identity is fixed
across all 22 autosomes. IBD1 takes one phased sample haplotype and one
phased donor haplotype, with fixed haplotype labels within each
chromosome and no crossover simulation, matching the chr20
simplification.

The provided `AF_EUR`/`AF_AFR` anchors are used directly, without
re-estimation or clamping. These 1000G AFs include the simulated genomes
in their full 2504/3202-genome frequency source: there is no
leave-two-out exclusion. This induces finite-panel dependence between
genotypes and frequencies and can make alignment optimistic relative to
held-out transfer. It can alter baseline heterozygosity, bias and
interval calibration; its direction and size are not isolated here.
CEU/YRI are narrower populations than EUR/AFR, so this campaign does not
remove within-superpopulation frequency mismatch.

Depth is Poisson; ALT conditional on depth is binomial with probability
`e + (1 - 2e) * ((1-share)*sample_dosage/2 + share*donor_dosage/2)`.
With probability 0.002 a site’s ALT count is uniform on 0..depth. Every
site’s depth and the same random draws are retained. Only covered sites
enter the bulk-loaded counts table; full-panel sites, zero-depth sites,
excluded sites and unavailable sites are reported separately per
replicate. Macro-native QC therefore describes covered input rows, not
all simulated sites. The macro uses `block_bases := 10000000`,
`min_blocks := 30`, default maximum depth 8, and the same error/artefact
calibration used to generate counts. The native extension and every
`src/` file are unchanged.

Declaration commit: `1735221a64bea0f49813090625e75cd7ea0d6acb`; batching
amendment: `957d03c91bfb3a49a0ae48286ae3275af6f446f2`; execution
amendment: `7fa78dc6ebdc371166e1551e15d3c1a46ac23112`. Campaign driver
revision: `fd2573193d2d1b1f82eb2098a0214352445fe0e4`.

R version 4.6.0 (2026-04-24); DuckDB R 1.5.5.

Input/extension/driver hashes:
[metadata.json](count_mixture_lowpass_genome/results/metadata.json).

Completed 18000 replicates in 90 strata; 15000 primary replicates.
Native output rows retained: 30000.

All 22 autosome BCFs have zero unphased records. Derived artifacts total
130.88 MiB; the limit is 4 GiB. All full-depth inputs contain 4883749
panel sites.

### Execution

Run from the worktree root with RBCFTools’ bundled bcftools in the
private R library, not a PATH executable. Staging and identity
verification are outside timed simulations. Recovery uses fourteen
long-lived R/DuckDB workers, `nice -n 10`, on CPU IDs **0–5,12–19**,
with one DuckDB thread per worker; explicit identity probes also use
four threads. Each worker loads chromosome panels once. Microbatches
contain 4/2/1 replicates at depths 0.1/0.2/0.4, with one macro call
grouping by sample identity. BLAS/OpenMP are limited to one thread.
Input identities and all-child phase are verified once at startup.
Per-seed outputs are atomically written under the manifest SHA-256, and
resumed only when source/input identities and settings match.

``` sh
taskset -c 0-5,12-19 env R_LIBS=/tmp/RBCFTools-release-lib \
  Rscript benchmarks/count_mixture_lowpass_genome/stage.R 1 build/release/duckhts.duckdb_extension
# Repeat staging for chromosomes 2–22.
taskset -c 0-5,12-19 nice -n 10 env R_LIBS=/tmp/RBCFTools-release-lib \
  Rscript benchmarks/count_mixture_lowpass_genome/run.R build/release/duckhts.duckdb_extension /tmp/lowpass-genome
```

Full per-replicate native outputs and process logs are retained by the
runner; checked-in gzip RDS tables preserve every seed, conditional
share/interval/status, family-disjoint identity and performance
denominator, with lossless doubles and categorical text labels. Read
them with `readRDS()`. Missing results are not replaced by successful
draws. Native outputs can be reconstructed from the registered
artifacts, manifest, declaration and seeds.

### Covered-site and microbatch checks

| depth | replicate |      seed | exact |   sites | zero_depth_sites | excluded_sites | unavailable_sites |
|------:|----------:|----------:|:------|--------:|-----------------:|---------------:|------------------:|
|   0.1 |         1 | 332206612 | TRUE  | 4883749 |          4419026 |              0 |                 0 |
|   0.1 |         2 | 332206613 | TRUE  | 4883749 |          4420419 |              0 |                 0 |
|   0.1 |         3 | 332206614 | TRUE  | 4883749 |          4419524 |              0 |                 0 |
|   0.4 |         1 | 332221612 | TRUE  | 4883749 |          3272685 |              0 |                 0 |
|   0.4 |         2 | 332221613 | TRUE  | 4883749 |          3274442 |              0 |                 0 |
|   0.4 |         3 | 332221614 | TRUE  | 4883749 |          3273258 |              0 |                 0 |
|   0.2 |         1 | 332236612 | TRUE  | 4883749 |          3996931 |              0 |                 0 |
|   0.2 |         2 | 332236613 | TRUE  | 4883749 |          3997608 |              0 |                 0 |
|   0.2 |         3 | 332236614 | TRUE  | 4883749 |          3998892 |              0 |                 0 |

| depth | batch_size | input_rows | simulation_seconds | loading_seconds | fit_seconds | buffer_bytes | spill_bytes | live_R_bytes |
|------:|-----------:|-----------:|-------------------:|----------------:|------------:|-------------:|------------:|-------------:|
|   0.1 |          4 |    1857208 |              1.468 |           0.323 |       4.227 |   1660424192 |           0 |    859326800 |
|   0.2 |          2 |    1772959 |              0.701 |           0.311 |       4.189 |   1660424192 |           0 |    847858456 |
|   0.4 |          1 |    1611064 |              0.558 |           0.294 |       3.304 |   1660424192 |           0 |    825835488 |

All-site and covered-site shares, interval endpoints, standard errors,
singles shares, pair residuals, statuses, diagnostics, blocks and pair
blocks must be bit-identical in all nine controls before using the
shortcut. Covered count rows and QC also match the original simulator
exactly. The declaration retains the full-site controls’ peak-buffer
caveat; campaign batches are chosen from the covered-input measurements,
not from the configured memory setting. The native control tables are
retained for direct IEEE-bit comparison.

## Release outcomes

| criterion                           | evaluated_strata | passing_strata |
|:------------------------------------|-----------------:|---------------:|
| Interior bias                       |               45 |             40 |
| 90% absolute error                  |               45 |             37 |
| Coverage and binomial compatibility |               75 |             53 |
| Null false-positive upper bound     |               15 |             10 |
| Matched diagnostic                  |               50 |             31 |
| All release criteria                |               75 |             30 |

### Detection and recovery

Positive-call rates below mean `share_lower > 0` under the supplied
generating relationship. They are observations, not a separate release
gate or relationship classification. Null positive calls are subject to
the declared binomial gate.

|     | relationship | depth | error | truth | replicates | mean_share |    bias | absolute_error_90 | accuracy | coverage | detection | release_pass |
|:----|:-------------|------:|------:|------:|-----------:|-----------:|--------:|------------------:|---------:|---------:|----------:|:-------------|
| 1   | same         |   0.1 | 1e-04 |  0.00 |        200 |     0.0046 |  0.0046 |            0.0124 |    1.000 |    0.945 |     0.055 | FALSE        |
| 2   | same         |   0.1 | 1e-04 |  0.02 |        200 |     0.0210 |  0.0010 |            0.0129 |    1.000 |    0.970 |     0.680 | FALSE        |
| 3   | same         |   0.1 | 1e-04 |  0.05 |        200 |     0.0501 |  0.0001 |            0.0145 |    1.000 |    0.960 |     1.000 | TRUE         |
| 4   | same         |   0.1 | 1e-04 |  0.10 |        200 |     0.1009 |  0.0009 |            0.0188 |    0.990 |    0.955 |     1.000 | TRUE         |
| 5   | same         |   0.1 | 1e-04 |  0.20 |        200 |     0.2018 |  0.0018 |            0.0252 |    0.895 |    0.950 |     1.000 | FALSE        |
| 6   | other        |   0.1 | 1e-04 |  0.00 |        200 |     0.0001 |  0.0001 |            0.0000 |    1.000 |    1.000 |     0.000 | FALSE        |
| 7   | other        |   0.1 | 1e-04 |  0.02 |        200 |     0.0144 | -0.0056 |            0.0114 |    1.000 |    0.775 |     0.885 | FALSE        |
| 8   | other        |   0.1 | 1e-04 |  0.05 |        200 |     0.0453 | -0.0047 |            0.0103 |    0.990 |    0.835 |     1.000 | FALSE        |
| 9   | other        |   0.1 | 1e-04 |  0.10 |        200 |     0.0971 | -0.0029 |            0.0092 |    1.000 |    0.900 |     1.000 | FALSE        |
| 10  | other        |   0.1 | 1e-04 |  0.20 |        200 |     0.2014 |  0.0014 |            0.0089 |    1.000 |    0.920 |     1.000 | TRUE         |
| 11  | ibd1         |   0.1 | 1e-04 |  0.00 |        200 |     0.0073 |  0.0073 |            0.0219 |    1.000 |    0.975 |     0.025 | FALSE        |
| 12  | ibd1         |   0.1 | 1e-04 |  0.02 |        200 |     0.0219 |  0.0019 |            0.0231 |    1.000 |    0.970 |     0.250 | FALSE        |
| 13  | ibd1         |   0.1 | 1e-04 |  0.05 |        200 |     0.0516 |  0.0016 |            0.0306 |    1.000 |    0.975 |     0.825 | FALSE        |
| 14  | ibd1         |   0.1 | 1e-04 |  0.10 |        200 |     0.1023 |  0.0023 |            0.0312 |    0.995 |    0.960 |     1.000 | TRUE         |
| 15  | ibd1         |   0.1 | 1e-04 |  0.20 |        200 |     0.2026 |  0.0026 |            0.0470 |    0.925 |    0.955 |     1.000 | TRUE         |
| 16  | same         |   0.4 | 1e-04 |  0.00 |        200 |     0.0016 |  0.0016 |            0.0054 |    0.980 |    0.950 |     0.050 | FALSE        |
| 17  | same         |   0.4 | 1e-04 |  0.02 |        200 |     0.0200 |  0.0000 |            0.0064 |    0.970 |    0.945 |     0.995 | FALSE        |
| 18  | same         |   0.4 | 1e-04 |  0.05 |        200 |     0.0503 |  0.0003 |            0.0061 |    0.965 |    0.960 |     1.000 | TRUE         |
| 19  | same         |   0.4 | 1e-04 |  0.10 |        200 |     0.1010 |  0.0010 |            0.0065 |    0.945 |    0.940 |     1.000 | TRUE         |
| 20  | same         |   0.4 | 1e-04 |  0.20 |        200 |     0.2008 |  0.0008 |            0.0085 |    0.870 |    0.950 |     1.000 | FALSE        |
| 21  | other        |   0.4 | 1e-04 |  0.00 |        200 |     0.0003 |  0.0003 |            0.0013 |    1.000 |    1.000 |     0.000 | FALSE        |
| 22  | other        |   0.4 | 1e-04 |  0.02 |        200 |     0.0184 | -0.0016 |            0.0047 |    1.000 |    0.890 |     1.000 | FALSE        |
| 23  | other        |   0.4 | 1e-04 |  0.05 |        200 |     0.0488 | -0.0012 |            0.0048 |    1.000 |    0.925 |     1.000 | TRUE         |
| 24  | other        |   0.4 | 1e-04 |  0.10 |        200 |     0.1003 |  0.0003 |            0.0049 |    1.000 |    0.920 |     1.000 | TRUE         |
| 25  | other        |   0.4 | 1e-04 |  0.20 |        200 |     0.2027 |  0.0027 |            0.0067 |    0.985 |    0.835 |     1.000 | FALSE        |
| 26  | ibd1         |   0.4 | 1e-04 |  0.00 |        200 |     0.0032 |  0.0032 |            0.0098 |    0.975 |    0.965 |     0.035 | FALSE        |
| 27  | ibd1         |   0.4 | 1e-04 |  0.02 |        200 |     0.0202 |  0.0002 |            0.0123 |    0.955 |    0.945 |     0.810 | FALSE        |
| 28  | ibd1         |   0.4 | 1e-04 |  0.05 |        200 |     0.0501 |  0.0001 |            0.0126 |    0.970 |    0.965 |     1.000 | TRUE         |
| 29  | ibd1         |   0.4 | 1e-04 |  0.10 |        200 |     0.1024 |  0.0024 |            0.0143 |    0.910 |    0.935 |     1.000 | TRUE         |
| 30  | ibd1         |   0.4 | 1e-04 |  0.20 |        200 |     0.2018 |  0.0018 |            0.0192 |    0.820 |    0.925 |     1.000 | FALSE        |
| 46  | same         |   0.1 | 2e-04 |  0.00 |        200 |     0.0040 |  0.0040 |            0.0131 |    1.000 |    0.965 |     0.035 | FALSE        |
| 47  | same         |   0.1 | 2e-04 |  0.02 |        200 |     0.0214 |  0.0014 |            0.0142 |    1.000 |    0.950 |     0.705 | FALSE        |
| 48  | same         |   0.1 | 2e-04 |  0.05 |        200 |     0.0505 |  0.0005 |            0.0148 |    0.990 |    0.955 |     1.000 | TRUE         |
| 49  | same         |   0.1 | 2e-04 |  0.10 |        200 |     0.1015 |  0.0015 |            0.0176 |    0.985 |    0.960 |     1.000 | TRUE         |
| 50  | same         |   0.1 | 2e-04 |  0.20 |        200 |     0.2006 |  0.0006 |            0.0247 |    0.900 |    0.935 |     1.000 | TRUE         |
| 51  | other        |   0.1 | 2e-04 |  0.00 |        200 |     0.0002 |  0.0002 |            0.0000 |    1.000 |    1.000 |     0.000 | FALSE        |
| 52  | other        |   0.1 | 2e-04 |  0.02 |        200 |     0.0137 | -0.0063 |            0.0117 |    1.000 |    0.745 |     0.845 | FALSE        |
| 53  | other        |   0.1 | 2e-04 |  0.05 |        200 |     0.0447 | -0.0053 |            0.0116 |    0.995 |    0.815 |     0.995 | FALSE        |
| 54  | other        |   0.1 | 2e-04 |  0.10 |        200 |     0.0974 | -0.0026 |            0.0098 |    1.000 |    0.905 |     1.000 | FALSE        |
| 55  | other        |   0.1 | 2e-04 |  0.20 |        200 |     0.2009 |  0.0009 |            0.0085 |    1.000 |    0.955 |     1.000 | TRUE         |
| 56  | ibd1         |   0.1 | 2e-04 |  0.00 |        200 |     0.0077 |  0.0077 |            0.0240 |    0.995 |    0.950 |     0.050 | FALSE        |
| 57  | ibd1         |   0.1 | 2e-04 |  0.02 |        200 |     0.0211 |  0.0011 |            0.0204 |    1.000 |    0.970 |     0.180 | FALSE        |
| 58  | ibd1         |   0.1 | 2e-04 |  0.05 |        200 |     0.0506 |  0.0006 |            0.0330 |    1.000 |    0.915 |     0.800 | FALSE        |
| 59  | ibd1         |   0.1 | 2e-04 |  0.10 |        200 |     0.1005 |  0.0005 |            0.0354 |    0.980 |    0.940 |     1.000 | TRUE         |
| 60  | ibd1         |   0.1 | 2e-04 |  0.20 |        200 |     0.2040 |  0.0040 |            0.0446 |    0.930 |    0.950 |     1.000 | TRUE         |
| 61  | same         |   0.4 | 2e-04 |  0.00 |        200 |     0.0015 |  0.0015 |            0.0047 |    0.995 |    0.975 |     0.025 | FALSE        |
| 62  | same         |   0.4 | 2e-04 |  0.02 |        200 |     0.0203 |  0.0003 |            0.0053 |    0.990 |    0.965 |     1.000 | TRUE         |
| 63  | same         |   0.4 | 2e-04 |  0.05 |        200 |     0.0502 |  0.0002 |            0.0071 |    0.935 |    0.900 |     1.000 | FALSE        |
| 64  | same         |   0.4 | 2e-04 |  0.10 |        200 |     0.1003 |  0.0003 |            0.0060 |    0.970 |    0.975 |     1.000 | TRUE         |
| 65  | same         |   0.4 | 2e-04 |  0.20 |        200 |     0.2009 |  0.0009 |            0.0083 |    0.885 |    0.945 |     0.995 | FALSE        |
| 66  | other        |   0.4 | 2e-04 |  0.00 |        200 |     0.0003 |  0.0003 |            0.0013 |    1.000 |    0.990 |     0.005 | FALSE        |
| 67  | other        |   0.4 | 2e-04 |  0.02 |        200 |     0.0181 | -0.0019 |            0.0052 |    1.000 |    0.870 |     1.000 | FALSE        |
| 68  | other        |   0.4 | 2e-04 |  0.05 |        200 |     0.0491 | -0.0009 |            0.0047 |    1.000 |    0.925 |     1.000 | TRUE         |
| 69  | other        |   0.4 | 2e-04 |  0.10 |        200 |     0.0998 | -0.0002 |            0.0046 |    1.000 |    0.935 |     1.000 | TRUE         |
| 70  | other        |   0.4 | 2e-04 |  0.20 |        200 |     0.2029 |  0.0029 |            0.0067 |    0.980 |    0.845 |     1.000 | FALSE        |
| 71  | ibd1         |   0.4 | 2e-04 |  0.00 |        200 |     0.0027 |  0.0027 |            0.0076 |    0.995 |    0.985 |     0.015 | FALSE        |
| 72  | ibd1         |   0.4 | 2e-04 |  0.02 |        200 |     0.0204 |  0.0004 |            0.0133 |    0.945 |    0.925 |     0.790 | FALSE        |
| 73  | ibd1         |   0.4 | 2e-04 |  0.05 |        200 |     0.0507 |  0.0007 |            0.0135 |    0.940 |    0.940 |     1.000 | TRUE         |
| 74  | ibd1         |   0.4 | 2e-04 |  0.10 |        200 |     0.1011 |  0.0011 |            0.0140 |    0.920 |    0.955 |     1.000 | TRUE         |
| 75  | ibd1         |   0.4 | 2e-04 |  0.20 |        200 |     0.2008 |  0.0008 |            0.0165 |    0.870 |    0.970 |     1.000 | FALSE        |
| 76  | same         |   0.2 | 2e-04 |  0.00 |        200 |     0.0023 |  0.0023 |            0.0073 |    0.990 |    0.970 |     0.030 | FALSE        |
| 77  | same         |   0.2 | 2e-04 |  0.02 |        200 |     0.0203 |  0.0003 |            0.0076 |    0.995 |    0.975 |     0.970 | TRUE         |
| 78  | same         |   0.2 | 2e-04 |  0.05 |        200 |     0.0507 |  0.0007 |            0.0075 |    0.990 |    0.940 |     1.000 | TRUE         |
| 79  | same         |   0.2 | 2e-04 |  0.10 |        200 |     0.1010 |  0.0010 |            0.0100 |    0.985 |    0.960 |     1.000 | TRUE         |
| 80  | same         |   0.2 | 2e-04 |  0.20 |        200 |     0.2020 |  0.0020 |            0.0144 |    0.905 |    0.915 |     0.995 | FALSE        |
| 81  | other        |   0.2 | 2e-04 |  0.00 |        200 |     0.0004 |  0.0004 |            0.0013 |    1.000 |    0.995 |     0.000 | FALSE        |
| 82  | other        |   0.2 | 2e-04 |  0.02 |        200 |     0.0169 | -0.0031 |            0.0068 |    1.000 |    0.860 |     1.000 | FALSE        |
| 83  | other        |   0.2 | 2e-04 |  0.05 |        200 |     0.0478 | -0.0022 |            0.0064 |    1.000 |    0.910 |     1.000 | FALSE        |
| 84  | other        |   0.2 | 2e-04 |  0.10 |        200 |     0.0991 | -0.0009 |            0.0067 |    1.000 |    0.940 |     1.000 | TRUE         |
| 85  | other        |   0.2 | 2e-04 |  0.20 |        200 |     0.2021 |  0.0021 |            0.0078 |    1.000 |    0.920 |     1.000 | TRUE         |
| 86  | ibd1         |   0.2 | 2e-04 |  0.00 |        200 |     0.0042 |  0.0042 |            0.0124 |    1.000 |    0.985 |     0.015 | FALSE        |
| 87  | ibd1         |   0.2 | 2e-04 |  0.02 |        200 |     0.0215 |  0.0015 |            0.0179 |    0.990 |    0.935 |     0.605 | FALSE        |
| 88  | ibd1         |   0.2 | 2e-04 |  0.05 |        200 |     0.0509 |  0.0009 |            0.0162 |    0.975 |    0.965 |     0.990 | TRUE         |
| 89  | ibd1         |   0.2 | 2e-04 |  0.10 |        200 |     0.1018 |  0.0018 |            0.0204 |    0.955 |    0.930 |     1.000 | TRUE         |
| 90  | ibd1         |   0.2 | 2e-04 |  0.20 |        200 |     0.2005 |  0.0005 |            0.0256 |    0.890 |    0.960 |     1.000 | FALSE        |

|     | relationship | depth | error | replicates | mean_share | detection | null_upper | null_pass |
|:----|:-------------|------:|------:|-----------:|-----------:|----------:|-----------:|:----------|
| 1   | same         |   0.1 | 1e-04 |        200 |     0.0046 |     0.055 |     0.0894 | FALSE     |
| 6   | other        |   0.1 | 1e-04 |        200 |     0.0001 |     0.000 |     0.0149 | TRUE      |
| 11  | ibd1         |   0.1 | 1e-04 |        200 |     0.0073 |     0.025 |     0.0518 | TRUE      |
| 16  | same         |   0.4 | 1e-04 |        200 |     0.0016 |     0.050 |     0.0833 | FALSE     |
| 21  | other        |   0.4 | 1e-04 |        200 |     0.0003 |     0.000 |     0.0149 | TRUE      |
| 26  | ibd1         |   0.4 | 1e-04 |        200 |     0.0032 |     0.035 |     0.0647 | TRUE      |
| 46  | same         |   0.1 | 2e-04 |        200 |     0.0040 |     0.035 |     0.0647 | TRUE      |
| 51  | other        |   0.1 | 2e-04 |        200 |     0.0002 |     0.000 |     0.0149 | TRUE      |
| 56  | ibd1         |   0.1 | 2e-04 |        200 |     0.0077 |     0.050 |     0.0833 | FALSE     |
| 61  | same         |   0.4 | 2e-04 |        200 |     0.0015 |     0.025 |     0.0518 | TRUE      |
| 66  | other        |   0.4 | 2e-04 |        200 |     0.0003 |     0.005 |     0.0235 | FALSE     |
| 71  | ibd1         |   0.4 | 2e-04 |        200 |     0.0027 |     0.015 |     0.0383 | TRUE      |
| 76  | same         |   0.2 | 2e-04 |        200 |     0.0023 |     0.030 |     0.0583 | TRUE      |
| 81  | other        |   0.2 | 2e-04 |        200 |     0.0004 |     0.000 |     0.0149 | FALSE     |
| 86  | ibd1         |   0.2 | 2e-04 |        200 |     0.0042 |     0.015 |     0.0383 | TRUE      |

### Intervals, statuses and ambiguous relationships

|     | relationship | depth | error | truth | interval_mean | interval_median | interval_90 | interval_max | full_range_intervals | missing_interval |
|:----|:-------------|------:|------:|------:|--------------:|----------------:|------------:|-------------:|---------------------:|-----------------:|
| 1   | same         |   0.1 | 1e-04 |  0.00 |        0.0191 |          0.0187 |      0.0302 |       0.0359 |                    0 |                0 |
| 2   | same         |   0.1 | 1e-04 |  0.02 |        0.0327 |          0.0338 |      0.0365 |       0.0386 |                    0 |                0 |
| 3   | same         |   0.1 | 1e-04 |  0.05 |        0.0367 |          0.0368 |      0.0389 |       0.0428 |                    0 |                0 |
| 4   | same         |   0.1 | 1e-04 |  0.10 |        0.0419 |          0.0420 |      0.0453 |       0.0498 |                    0 |                0 |
| 5   | same         |   0.1 | 1e-04 |  0.20 |        0.0579 |          0.0578 |      0.0635 |       0.0676 |                    0 |                0 |
| 6   | other        |   0.1 | 1e-04 |  0.00 |        0.0048 |          0.0042 |      0.0074 |       0.0139 |                    0 |                0 |
| 7   | other        |   0.1 | 1e-04 |  0.02 |        0.0178 |          0.0180 |      0.0193 |       0.0223 |                    0 |                0 |
| 8   | other        |   0.1 | 1e-04 |  0.05 |        0.0188 |          0.0187 |      0.0205 |       0.0235 |                    0 |                0 |
| 9   | other        |   0.1 | 1e-04 |  0.10 |        0.0195 |          0.0194 |      0.0213 |       0.0244 |                    0 |                0 |
| 10  | other        |   0.1 | 1e-04 |  0.20 |        0.0200 |          0.0199 |      0.0217 |       0.0239 |                    0 |                0 |
| 11  | ibd1         |   0.1 | 1e-04 |  0.00 |        0.0359 |          0.0332 |      0.0572 |       0.0756 |                    0 |                0 |
| 12  | ibd1         |   0.1 | 1e-04 |  0.02 |        0.0545 |          0.0584 |      0.0713 |       0.0816 |                    0 |                0 |
| 13  | ibd1         |   0.1 | 1e-04 |  0.05 |        0.0721 |          0.0731 |      0.0807 |       0.0892 |                    0 |                0 |
| 14  | ibd1         |   0.1 | 1e-04 |  0.10 |        0.0839 |          0.0837 |      0.0928 |       0.1019 |                    0 |                0 |
| 15  | ibd1         |   0.1 | 1e-04 |  0.20 |        0.1166 |          0.1141 |      0.1383 |       0.1650 |                    0 |                0 |
| 16  | same         |   0.4 | 1e-04 |  0.00 |        0.0072 |          0.0061 |      0.0123 |       0.0143 |                    0 |                0 |
| 17  | same         |   0.4 | 1e-04 |  0.02 |        0.0144 |          0.0141 |      0.0157 |       0.0189 |                    0 |                1 |
| 18  | same         |   0.4 | 1e-04 |  0.05 |        0.0148 |          0.0146 |      0.0159 |       0.0202 |                    0 |                0 |
| 19  | same         |   0.4 | 1e-04 |  0.10 |        0.0159 |          0.0158 |      0.0174 |       0.0200 |                    0 |                0 |
| 20  | same         |   0.4 | 1e-04 |  0.20 |        0.0203 |          0.0203 |      0.0222 |       0.0238 |                    0 |                0 |
| 21  | other        |   0.4 | 1e-04 |  0.00 |        0.0032 |          0.0026 |      0.0057 |       0.0093 |                    0 |                0 |
| 22  | other        |   0.4 | 1e-04 |  0.02 |        0.0095 |          0.0095 |      0.0102 |       0.0118 |                    0 |                0 |
| 23  | other        |   0.4 | 1e-04 |  0.05 |        0.0100 |          0.0099 |      0.0107 |       0.0119 |                    0 |                0 |
| 24  | other        |   0.4 | 1e-04 |  0.10 |        0.0105 |          0.0104 |      0.0113 |       0.0122 |                    0 |                0 |
| 25  | other        |   0.4 | 1e-04 |  0.20 |        0.0121 |          0.0121 |      0.0132 |       0.0143 |                    0 |                0 |
| 26  | ibd1         |   0.4 | 1e-04 |  0.00 |        0.0151 |          0.0149 |      0.0230 |       0.0298 |                    0 |                0 |
| 27  | ibd1         |   0.4 | 1e-04 |  0.02 |        0.0278 |          0.0281 |      0.0307 |       0.0338 |                    0 |                0 |
| 28  | ibd1         |   0.4 | 1e-04 |  0.05 |        0.0302 |          0.0299 |      0.0326 |       0.0383 |                    0 |                0 |
| 29  | ibd1         |   0.4 | 1e-04 |  0.10 |        0.0336 |          0.0334 |      0.0363 |       0.0435 |                    0 |                0 |
| 30  | ibd1         |   0.4 | 1e-04 |  0.20 |        0.0432 |          0.0429 |      0.0467 |       0.0526 |                    0 |                0 |
| 46  | same         |   0.1 | 2e-04 |  0.00 |        0.0179 |          0.0154 |      0.0302 |       0.0383 |                    0 |                0 |
| 47  | same         |   0.1 | 2e-04 |  0.02 |        0.0328 |          0.0340 |      0.0364 |       0.0391 |                    0 |                0 |
| 48  | same         |   0.1 | 2e-04 |  0.05 |        0.0369 |          0.0369 |      0.0396 |       0.0426 |                    0 |                0 |
| 49  | same         |   0.1 | 2e-04 |  0.10 |        0.0422 |          0.0423 |      0.0453 |       0.0518 |                    0 |                0 |
| 50  | same         |   0.1 | 2e-04 |  0.20 |        0.0578 |          0.0577 |      0.0638 |       0.0691 |                    0 |                0 |
| 51  | other        |   0.1 | 2e-04 |  0.00 |        0.0051 |          0.0043 |      0.0085 |       0.0157 |                    0 |                0 |
| 52  | other        |   0.1 | 2e-04 |  0.02 |        0.0177 |          0.0179 |      0.0194 |       0.0212 |                    0 |                0 |
| 53  | other        |   0.1 | 2e-04 |  0.05 |        0.0187 |          0.0188 |      0.0199 |       0.0220 |                    0 |                1 |
| 54  | other        |   0.1 | 2e-04 |  0.10 |        0.0195 |          0.0194 |      0.0209 |       0.0232 |                    0 |                0 |
| 55  | other        |   0.1 | 2e-04 |  0.20 |        0.0199 |          0.0199 |      0.0215 |       0.0235 |                    0 |                0 |
| 56  | ibd1         |   0.1 | 2e-04 |  0.00 |        0.0361 |          0.0328 |      0.0610 |       0.0746 |                    0 |                0 |
| 57  | ibd1         |   0.1 | 2e-04 |  0.02 |        0.0540 |          0.0540 |      0.0725 |       0.0828 |                    0 |                0 |
| 58  | ibd1         |   0.1 | 2e-04 |  0.05 |        0.0716 |          0.0728 |      0.0805 |       0.0916 |                    0 |                0 |
| 59  | ibd1         |   0.1 | 2e-04 |  0.10 |        0.0839 |          0.0835 |      0.0919 |       0.1036 |                    0 |                0 |
| 60  | ibd1         |   0.1 | 2e-04 |  0.20 |        0.1168 |          0.1136 |      0.1363 |       0.1778 |                    0 |                0 |
| 61  | same         |   0.4 | 2e-04 |  0.00 |        0.0073 |          0.0067 |      0.0113 |       0.0148 |                    0 |                0 |
| 62  | same         |   0.4 | 2e-04 |  0.02 |        0.0142 |          0.0140 |      0.0155 |       0.0184 |                    0 |                0 |
| 63  | same         |   0.4 | 2e-04 |  0.05 |        0.0148 |          0.0146 |      0.0161 |       0.0185 |                    0 |                0 |
| 64  | same         |   0.4 | 2e-04 |  0.10 |        0.0161 |          0.0160 |      0.0173 |       0.0205 |                    0 |                0 |
| 65  | same         |   0.4 | 2e-04 |  0.20 |        0.0203 |          0.0202 |      0.0219 |       0.0241 |                    0 |                1 |
| 66  | other        |   0.4 | 2e-04 |  0.00 |        0.0032 |          0.0027 |      0.0058 |       0.0088 |                    0 |                1 |
| 67  | other        |   0.4 | 2e-04 |  0.02 |        0.0096 |          0.0095 |      0.0104 |       0.0124 |                    0 |                0 |
| 68  | other        |   0.4 | 2e-04 |  0.05 |        0.0100 |          0.0099 |      0.0107 |       0.0123 |                    0 |                0 |
| 69  | other        |   0.4 | 2e-04 |  0.10 |        0.0106 |          0.0105 |      0.0114 |       0.0125 |                    0 |                0 |
| 70  | other        |   0.4 | 2e-04 |  0.20 |        0.0122 |          0.0120 |      0.0131 |       0.0145 |                    0 |                0 |
| 71  | ibd1         |   0.4 | 2e-04 |  0.00 |        0.0144 |          0.0138 |      0.0212 |       0.0304 |                    0 |                0 |
| 72  | ibd1         |   0.4 | 2e-04 |  0.02 |        0.0278 |          0.0279 |      0.0306 |       0.0358 |                    0 |                0 |
| 73  | ibd1         |   0.4 | 2e-04 |  0.05 |        0.0303 |          0.0300 |      0.0331 |       0.0396 |                    0 |                0 |
| 74  | ibd1         |   0.4 | 2e-04 |  0.10 |        0.0334 |          0.0332 |      0.0364 |       0.0412 |                    0 |                0 |
| 75  | ibd1         |   0.4 | 2e-04 |  0.20 |        0.0434 |          0.0432 |      0.0465 |       0.0535 |                    0 |                0 |
| 76  | same         |   0.2 | 2e-04 |  0.00 |        0.0105 |          0.0099 |      0.0172 |       0.0211 |                    0 |                0 |
| 77  | same         |   0.2 | 2e-04 |  0.02 |        0.0202 |          0.0201 |      0.0218 |       0.0242 |                    0 |                0 |
| 78  | same         |   0.2 | 2e-04 |  0.05 |        0.0214 |          0.0214 |      0.0230 |       0.0260 |                    0 |                0 |
| 79  | same         |   0.2 | 2e-04 |  0.10 |        0.0239 |          0.0238 |      0.0255 |       0.0278 |                    0 |                0 |
| 80  | same         |   0.2 | 2e-04 |  0.20 |        0.0319 |          0.0320 |      0.0342 |       0.0382 |                    0 |                1 |
| 81  | other        |   0.2 | 2e-04 |  0.00 |        0.0043 |          0.0038 |      0.0077 |       0.0133 |                    0 |                1 |
| 82  | other        |   0.2 | 2e-04 |  0.02 |        0.0129 |          0.0128 |      0.0139 |       0.0150 |                    0 |                0 |
| 83  | other        |   0.2 | 2e-04 |  0.05 |        0.0135 |          0.0134 |      0.0143 |       0.0159 |                    0 |                0 |
| 84  | other        |   0.2 | 2e-04 |  0.10 |        0.0145 |          0.0144 |      0.0156 |       0.0169 |                    0 |                0 |
| 85  | other        |   0.2 | 2e-04 |  0.20 |        0.0159 |          0.0158 |      0.0171 |       0.0187 |                    0 |                0 |
| 86  | ibd1         |   0.2 | 2e-04 |  0.00 |        0.0215 |          0.0214 |      0.0322 |       0.0404 |                    0 |                0 |
| 87  | ibd1         |   0.2 | 2e-04 |  0.02 |        0.0371 |          0.0393 |      0.0432 |       0.0465 |                    0 |                0 |
| 88  | ibd1         |   0.2 | 2e-04 |  0.05 |        0.0430 |          0.0429 |      0.0461 |       0.0535 |                    0 |                0 |
| 89  | ibd1         |   0.2 | 2e-04 |  0.10 |        0.0485 |          0.0484 |      0.0522 |       0.0605 |                    0 |                0 |
| 90  | ibd1         |   0.2 | 2e-04 |  0.20 |        0.0643 |          0.0638 |      0.0710 |       0.0857 |                    0 |                0 |

| relationship_model | status         | relation_diagnostic  | count |
|:-------------------|:---------------|:---------------------|------:|
| unrelated          | ok_conditional | compatible_unrelated |  4169 |
| unrelated          | model_mismatch | model_mismatch       |     4 |
| unrelated          | ok_conditional | model_mismatch       |   791 |
| unrelated          | ok_conditional | shared_allele_signal |    31 |
| ibd1               | ok_conditional | uninformative        |  1158 |
| unrelated          | ok_conditional | uninformative        |  2163 |
| ibd1               | model_mismatch | unresolved           |   949 |
| unrelated          | model_mismatch | unresolved           |     6 |
| ibd1               | ok_conditional | unresolved           |  9893 |
| unrelated          | ok_conditional | unresolved           | 10836 |

| failed_statuses | missing_shares | missing_intervals | failed_jackknife_refits |
|----------------:|---------------:|------------------:|------------------------:|
|              10 |              0 |                10 |                       0 |

Both matched conditional models remain in the gzip table. Their
existence and numerically successful statuses do not resolve the
same-superpopulation relationship or establish calibrated confidence
intervals. Across all matched native rows, 21,684 diagnostics are
`unresolved` and 2,316 are `uninformative`; the strict
all-replicates-`unresolved` criterion passes in 31/50 primary matched
strata. The criterion is not broadened to include `uninformative`.

### Secondary depth

|     | relationship | depth | error | truth | replicates | mean_share |    bias | accuracy | coverage | detection | release_pass |
|:----|:-------------|------:|------:|------:|-----------:|-----------:|--------:|---------:|---------:|----------:|:-------------|
| 31  | same         |   0.2 | 1e-04 |  0.00 |        200 |     0.0023 |  0.0023 |    1.000 |    0.970 |     0.030 | FALSE        |
| 32  | same         |   0.2 | 1e-04 |  0.02 |        200 |     0.0206 |  0.0006 |    0.990 |    0.955 |     0.980 | TRUE         |
| 33  | same         |   0.2 | 1e-04 |  0.05 |        200 |     0.0507 |  0.0007 |    0.995 |    0.940 |     1.000 | TRUE         |
| 34  | same         |   0.2 | 1e-04 |  0.10 |        200 |     0.1010 |  0.0010 |    0.980 |    0.945 |     0.995 | FALSE        |
| 35  | same         |   0.2 | 1e-04 |  0.20 |        200 |     0.2018 |  0.0018 |    0.930 |    0.940 |     0.995 | FALSE        |
| 36  | other        |   0.2 | 1e-04 |  0.00 |        200 |     0.0002 |  0.0002 |    1.000 |    0.995 |     0.000 | FALSE        |
| 37  | other        |   0.2 | 1e-04 |  0.02 |        200 |     0.0172 | -0.0028 |    1.000 |    0.860 |     0.995 | FALSE        |
| 38  | other        |   0.2 | 1e-04 |  0.05 |        200 |     0.0477 | -0.0023 |    1.000 |    0.895 |     1.000 | FALSE        |
| 39  | other        |   0.2 | 1e-04 |  0.10 |        200 |     0.0990 | -0.0010 |    1.000 |    0.925 |     1.000 | TRUE         |
| 40  | other        |   0.2 | 1e-04 |  0.20 |        200 |     0.2033 |  0.0033 |    0.995 |    0.910 |     1.000 | FALSE        |
| 41  | ibd1         |   0.2 | 1e-04 |  0.00 |        200 |     0.0047 |  0.0047 |    0.985 |    0.955 |     0.040 | FALSE        |
| 42  | ibd1         |   0.2 | 1e-04 |  0.02 |        200 |     0.0226 |  0.0026 |    0.990 |    0.945 |     0.630 | FALSE        |
| 43  | ibd1         |   0.2 | 1e-04 |  0.05 |        200 |     0.0520 |  0.0020 |    0.985 |    0.965 |     1.000 | TRUE         |
| 44  | ibd1         |   0.2 | 1e-04 |  0.10 |        200 |     0.1022 |  0.0022 |    0.950 |    0.950 |     1.000 | TRUE         |
| 45  | ibd1         |   0.2 | 1e-04 |  0.20 |        200 |     0.2027 |  0.0027 |    0.865 |    0.940 |     1.000 | FALSE        |

## Independent numerical checks

| checked_strata | max_R_moment_error | counter_failures | status_failures | share_missingness_failures | interval_missingness_failures | diagnostic_failures | max_thread_share_error | max_thread_interval_error | max_thread_moment_error | solo_batch_failures | bit_identical_thread_checks |
|---------------:|-------------------:|-----------------:|----------------:|---------------------------:|------------------------------:|--------------------:|-----------------------:|--------------------------:|------------------------:|--------------------:|----------------------------:|
|             90 |                  0 |                0 |               0 |                          0 |                             0 |                   0 |                      0 |                         0 |                       0 |                   0 |                          90 |

The R oracle reconstructs all moments and covered-input QC counters for
every physical chromosome/10 Mb block, using direct weighted
cross-products, not the SQL expressions. Full-panel and zero-depth QC
are retained separately for every seed; zero-depth rows contribute no
moments. The SQL probe exposes the actual moment rows passed to the
native aggregate. No implementation output is used as its own oracle.

## Performance and scaling

| component                           |  seconds | replicates | workers | threads_per_worker |
|:------------------------------------|---------:|-----------:|--------:|-------------------:|
| Startup and orchestration remainder |   700.22 |          0 |       1 |                  1 |
| 200-replicate pilot                 |   243.28 |        200 |       1 |                  1 |
| Concurrent campaign                 | 10373.50 |      17800 |      14 |                  1 |
| Total command                       | 11317.00 |      18000 |      14 |                  1 |

    ## Simulated site exposures: 87907482000; covered counts rows: 17760643873; output rows: 30000.

The 200-replicate pilot took 243.3 seconds (1.216 seconds per
replicate). The concurrent campaign took 10,373.5 seconds (2.88 hours);
the whole command, including startup identity/phase audits and the
pilot, took 3:08:37. The original eight-hour budget is met, but the
roughly 1.5-hour target is missed. The assigned fourteen logical CPUs
comprise three SMT core pairs and eight other cores on an Intel
i5-13500; a single-worker probe did not predict concurrent throughput.
No strata or replicates are dropped to meet the target.

Campaign timings are amortized within each microbatch and include
neither per-replicate process startup nor repeated panel reading. Worker
RSS is a long-lived-process high-water mark, not an isolated
per-replicate measurement.

| stratum | simulation_seconds | loading_seconds | fit_seconds | wall_seconds |
|--------:|-------------------:|----------------:|------------:|-------------:|
|       1 |             0.7491 |          0.2003 |      3.0498 |       4.1241 |
|       2 |             0.7583 |          0.1882 |      3.0972 |       4.1894 |
|       3 |             0.7646 |          0.2037 |      3.0890 |       4.2082 |
|       4 |             0.2192 |          0.0672 |      0.8581 |       1.1693 |
|       5 |             0.7377 |          0.1947 |      3.1079 |       4.1929 |
|       6 |             0.7415 |          0.2084 |      3.0460 |       4.1382 |
|       7 |             0.7523 |          0.1876 |      3.0863 |       4.1692 |
|       8 |             0.7742 |          0.1998 |      3.1512 |       4.2408 |
|       9 |             0.7372 |          0.1934 |      3.1034 |       4.1578 |
|      10 |             0.7409 |          0.1997 |      3.1028 |       4.1815 |
|      11 |             0.7808 |          0.1944 |      3.1290 |       4.2219 |
|      12 |             0.7910 |          0.1908 |      3.0788 |       4.2111 |
|      13 |             0.7441 |          0.2065 |      3.1272 |       4.2117 |
|      14 |             0.7814 |          0.2087 |      3.1688 |       4.2874 |
|      15 |             0.7852 |          0.1957 |      3.0728 |       4.1845 |
|      16 |             1.4360 |          0.5690 |      9.3339 |      11.6601 |
|      17 |             1.4012 |          0.5735 |      9.1666 |      11.4808 |
|      18 |             1.4240 |          0.5734 |      9.3888 |      11.7638 |
|      19 |             1.4444 |          0.6015 |      9.6178 |      12.0490 |
|      20 |             1.4177 |          0.5953 |      9.3947 |      11.7821 |
|      21 |             1.4488 |          0.5891 |      9.3375 |      11.7096 |
|      22 |             1.4236 |          0.5514 |      9.1794 |      11.5071 |
|      23 |             1.4263 |          0.5868 |      9.4423 |      11.8065 |
|      24 |             1.4433 |          0.5737 |      9.3110 |      11.7026 |
|      25 |             1.4323 |          0.5865 |      9.3381 |      11.7166 |
|      26 |             1.3236 |          0.5515 |      9.1147 |      11.3808 |
|      27 |             1.4179 |          0.5716 |      9.3991 |      11.7484 |
|      28 |             1.4072 |          0.5972 |      9.4443 |      11.8097 |
|      29 |             1.4183 |          0.5829 |      9.3343 |      11.6890 |
|      30 |             1.2926 |          0.5271 |      8.8047 |      10.9591 |
|      31 |             0.9805 |          0.3155 |      5.0080 |       6.4793 |
|      32 |             0.9515 |          0.3124 |      4.9517 |       6.4023 |
|      33 |             0.9211 |          0.3050 |      4.9145 |       6.3237 |
|      34 |             0.9014 |          0.3192 |      4.9816 |       6.4009 |
|      35 |             0.9220 |          0.2980 |      4.8809 |       6.2876 |
|      36 |             0.8975 |          0.2929 |      4.9048 |       6.2737 |
|      37 |             0.9022 |          0.2940 |      4.8871 |       6.2803 |
|      38 |             0.9768 |          0.3191 |      4.9954 |       6.4633 |
|      39 |             0.9373 |          0.3116 |      5.0161 |       6.4479 |
|      40 |             0.9291 |          0.3067 |      4.9526 |       6.3752 |
|      41 |             0.9134 |          0.3028 |      4.9990 |       6.3967 |
|      42 |             0.9359 |          0.3288 |      5.0871 |       6.5410 |
|      43 |             0.8916 |          0.2888 |      4.8837 |       6.2514 |
|      44 |             0.9806 |          0.3291 |      5.3689 |       6.8481 |
|      45 |             0.9161 |          0.3071 |      4.9996 |       6.3934 |
|      46 |             0.6295 |          0.1545 |      2.6143 |       3.4884 |
|      47 |             0.6149 |          0.1578 |      2.6314 |       3.4830 |
|      48 |             0.6472 |          0.1709 |      2.6136 |       3.5226 |
|      49 |             0.6371 |          0.1566 |      2.6044 |       3.4874 |
|      50 |             0.6396 |          0.1582 |      2.6575 |       3.5393 |
|      51 |             0.6271 |          0.1484 |      2.4872 |       3.3447 |
|      52 |             0.6359 |          0.1576 |      2.5863 |       3.4770 |
|      53 |             0.6416 |          0.1571 |      2.5983 |       3.4923 |
|      54 |             0.6547 |          0.1613 |      2.5606 |       3.4628 |
|      55 |             0.6393 |          0.1624 |      2.6071 |       3.5021 |
|      56 |             0.6230 |          0.1558 |      2.5976 |       3.4750 |
|      57 |             0.6366 |          0.1576 |      2.6195 |       3.5052 |
|      58 |             0.7426 |          0.1805 |      3.0477 |       4.0728 |
|      59 |             0.7171 |          0.1869 |      3.0012 |       3.9898 |
|      60 |             0.6700 |          0.1632 |      2.7000 |       3.6157 |
|      61 |             1.4848 |          0.6363 |      9.9435 |      12.3728 |
|      62 |             1.5109 |          0.6695 |     10.1266 |      12.6309 |
|      63 |             1.5008 |          0.6300 |      9.8351 |      12.2663 |
|      64 |             1.5287 |          0.6321 |      9.8680 |      12.3394 |
|      65 |             1.5021 |          0.6343 |      9.8463 |      12.2915 |
|      66 |             1.6106 |          0.6253 |      9.9856 |      12.5372 |
|      67 |             1.4521 |          0.6130 |      9.9229 |      12.3279 |
|      68 |             1.4364 |          0.5969 |      9.8602 |      12.1965 |
|      69 |             1.5014 |          0.6169 |      9.7771 |      12.2191 |
|      70 |             1.5403 |          0.6648 |      9.9478 |      12.4896 |
|      71 |             1.5229 |          0.6248 |     10.0007 |      12.4704 |
|      72 |             1.2088 |          0.5460 |      7.7639 |       9.8519 |
|      73 |             1.5148 |          0.6755 |     10.1656 |      12.6967 |
|      74 |             1.5405 |          0.6806 |     10.4645 |      13.0190 |
|      75 |             1.0852 |          0.5149 |      7.2859 |       9.2484 |
|      76 |             0.9701 |          0.3806 |      5.4571 |       7.0100 |
|      77 |             0.9899 |          0.3761 |      5.6237 |       7.1788 |
|      78 |             1.0133 |          0.3751 |      5.6164 |       7.2202 |
|      79 |             1.0219 |          0.3834 |      5.5821 |       7.1998 |
|      80 |             0.9623 |          0.3591 |      5.4892 |       7.0060 |
|      81 |             0.9959 |          0.3497 |      5.4954 |       7.0367 |
|      82 |             1.0581 |          0.3835 |      5.5687 |       7.1996 |
|      83 |             0.9881 |          0.3641 |      5.5942 |       7.1386 |
|      84 |             0.9655 |          0.3754 |      5.3556 |       6.8843 |
|      85 |             1.0223 |          0.3699 |      5.5223 |       7.1153 |
|      86 |             0.3639 |          0.1509 |      1.9215 |       2.4982 |
|      87 |             0.4530 |          0.1892 |      2.3894 |       3.1238 |
|      88 |             0.6909 |          0.2634 |      3.6199 |       4.7595 |
|      89 |             0.3724 |          0.1565 |      1.9311 |       2.5250 |
|      90 |             0.4402 |          0.1891 |      2.3205 |       3.0281 |

| stratum |   sites | covered_sites | zero_depth_sites | excluded_sites | unavailable_sites |
|--------:|--------:|--------------:|-----------------:|---------------:|------------------:|
|       1 | 4883749 |      464803.3 |          4418946 |          0.000 |                 0 |
|       2 | 4883749 |      464677.9 |          4419071 |          0.000 |                 0 |
|       3 | 4883749 |      464796.5 |          4418952 |          0.000 |                 0 |
|       4 | 4883749 |      464771.6 |          4418977 |          0.000 |                 0 |
|       5 | 4883749 |      464687.1 |          4419062 |          0.000 |                 0 |
|       6 | 4883749 |      464697.8 |          4419051 |          0.000 |                 0 |
|       7 | 4883749 |      464723.5 |          4419026 |          0.000 |                 0 |
|       8 | 4883749 |      464755.7 |          4418993 |          0.000 |                 0 |
|       9 | 4883749 |      464795.2 |          4418954 |          0.000 |                 0 |
|      10 | 4883749 |      464814.2 |          4418935 |          0.000 |                 0 |
|      11 | 4883749 |      464809.5 |          4418939 |          0.000 |                 0 |
|      12 | 4883749 |      464787.6 |          4418961 |          0.000 |                 0 |
|      13 | 4883749 |      464735.8 |          4419013 |          0.000 |                 0 |
|      14 | 4883749 |      464816.6 |          4418932 |          0.000 |                 0 |
|      15 | 4883749 |      464797.0 |          4418952 |          0.000 |                 0 |
|      16 | 4883749 |     1610012.3 |          3273737 |          0.000 |                 0 |
|      17 | 4883749 |     1610084.2 |          3273665 |          0.000 |                 0 |
|      18 | 4883749 |     1610055.0 |          3273694 |          0.000 |                 0 |
|      19 | 4883749 |     1610138.0 |          3273611 |          0.000 |                 0 |
|      20 | 4883749 |     1610060.1 |          3273689 |          0.000 |                 0 |
|      21 | 4883749 |     1609982.7 |          3273766 |          0.010 |                 0 |
|      22 | 4883749 |     1609994.0 |          3273755 |          0.000 |                 0 |
|      23 | 4883749 |     1610204.8 |          3273544 |          0.000 |                 0 |
|      24 | 4883749 |     1610106.0 |          3273643 |          0.000 |                 0 |
|      25 | 4883749 |     1610077.9 |          3273671 |          0.000 |                 0 |
|      26 | 4883749 |     1610073.3 |          3273676 |          0.000 |                 0 |
|      27 | 4883749 |     1610133.4 |          3273616 |          0.000 |                 0 |
|      28 | 4883749 |     1610084.7 |          3273664 |          0.005 |                 0 |
|      29 | 4883749 |     1610097.0 |          3273652 |          0.005 |                 0 |
|      30 | 4883749 |     1610155.4 |          3273594 |          0.000 |                 0 |
|      31 | 4883749 |      885299.5 |          3998449 |          0.000 |                 0 |
|      32 | 4883749 |      885184.8 |          3998564 |          0.000 |                 0 |
|      33 | 4883749 |      885269.2 |          3998480 |          0.000 |                 0 |
|      34 | 4883749 |      885160.7 |          3998588 |          0.000 |                 0 |
|      35 | 4883749 |      885380.6 |          3998368 |          0.000 |                 0 |
|      36 | 4883749 |      885155.8 |          3998593 |          0.000 |                 0 |
|      37 | 4883749 |      885276.2 |          3998473 |          0.000 |                 0 |
|      38 | 4883749 |      885313.0 |          3998436 |          0.000 |                 0 |
|      39 | 4883749 |      885158.2 |          3998591 |          0.000 |                 0 |
|      40 | 4883749 |      885132.1 |          3998617 |          0.000 |                 0 |
|      41 | 4883749 |      885301.3 |          3998448 |          0.000 |                 0 |
|      42 | 4883749 |      885217.4 |          3998532 |          0.000 |                 0 |
|      43 | 4883749 |      885267.0 |          3998482 |          0.000 |                 0 |
|      44 | 4883749 |      885307.0 |          3998442 |          0.000 |                 0 |
|      45 | 4883749 |      885203.4 |          3998546 |          0.000 |                 0 |
|      46 | 4883749 |      464707.7 |          4419041 |          0.000 |                 0 |
|      47 | 4883749 |      464769.2 |          4418980 |          0.000 |                 0 |
|      48 | 4883749 |      464731.1 |          4419018 |          0.000 |                 0 |
|      49 | 4883749 |      464775.6 |          4418973 |          0.000 |                 0 |
|      50 | 4883749 |      464809.6 |          4418939 |          0.000 |                 0 |
|      51 | 4883749 |      464650.1 |          4419099 |          0.000 |                 0 |
|      52 | 4883749 |      464758.6 |          4418990 |          0.000 |                 0 |
|      53 | 4883749 |      464781.5 |          4418968 |          0.000 |                 0 |
|      54 | 4883749 |      464679.2 |          4419070 |          0.000 |                 0 |
|      55 | 4883749 |      464797.7 |          4418951 |          0.000 |                 0 |
|      56 | 4883749 |      464715.5 |          4419034 |          0.000 |                 0 |
|      57 | 4883749 |      464807.9 |          4418941 |          0.000 |                 0 |
|      58 | 4883749 |      464699.5 |          4419050 |          0.000 |                 0 |
|      59 | 4883749 |      464779.6 |          4418969 |          0.000 |                 0 |
|      60 | 4883749 |      464814.2 |          4418935 |          0.000 |                 0 |
|      61 | 4883749 |     1610159.4 |          3273590 |          0.005 |                 0 |
|      62 | 4883749 |     1609974.3 |          3273775 |          0.010 |                 0 |
|      63 | 4883749 |     1610075.8 |          3273673 |          0.005 |                 0 |
|      64 | 4883749 |     1610181.0 |          3273568 |          0.000 |                 0 |
|      65 | 4883749 |     1610133.5 |          3273616 |          0.000 |                 0 |
|      66 | 4883749 |     1610041.5 |          3273708 |          0.005 |                 0 |
|      67 | 4883749 |     1610010.5 |          3273739 |          0.005 |                 0 |
|      68 | 4883749 |     1610005.1 |          3273744 |          0.000 |                 0 |
|      69 | 4883749 |     1610136.4 |          3273613 |          0.010 |                 0 |
|      70 | 4883749 |     1610088.2 |          3273661 |          0.000 |                 0 |
|      71 | 4883749 |     1609958.4 |          3273791 |          0.000 |                 0 |
|      72 | 4883749 |     1610160.9 |          3273588 |          0.005 |                 0 |
|      73 | 4883749 |     1610048.4 |          3273701 |          0.005 |                 0 |
|      74 | 4883749 |     1610095.7 |          3273653 |          0.000 |                 0 |
|      75 | 4883749 |     1610044.5 |          3273705 |          0.000 |                 0 |
|      76 | 4883749 |      885279.1 |          3998470 |          0.000 |                 0 |
|      77 | 4883749 |      885393.7 |          3998355 |          0.000 |                 0 |
|      78 | 4883749 |      885266.7 |          3998482 |          0.000 |                 0 |
|      79 | 4883749 |      885339.8 |          3998409 |          0.000 |                 0 |
|      80 | 4883749 |      885265.4 |          3998484 |          0.000 |                 0 |
|      81 | 4883749 |      885280.6 |          3998468 |          0.000 |                 0 |
|      82 | 4883749 |      885368.6 |          3998380 |          0.000 |                 0 |
|      83 | 4883749 |      885232.4 |          3998517 |          0.000 |                 0 |
|      84 | 4883749 |      885206.0 |          3998543 |          0.000 |                 0 |
|      85 | 4883749 |      885234.0 |          3998515 |          0.000 |                 0 |
|      86 | 4883749 |      885360.8 |          3998388 |          0.000 |                 0 |
|      87 | 4883749 |      885315.0 |          3998434 |          0.000 |                 0 |
|      88 | 4883749 |      885194.6 |          3998554 |          0.000 |                 0 |
|      89 | 4883749 |      885313.0 |          3998436 |          0.000 |                 0 |
|      90 | 4883749 |      885424.8 |          3998324 |          0.000 |                 0 |

Fresh isolated process timings below include process startup, extension
load, artifact reading, count simulation, ingestion and fitting in the
wall total. The simulation column covers chromosome
reading/simulation/ingestion; fit time covers only the macro query. Peak
RSS covers the entire process. These runs exclude the R oracle and
thread probes. Scaling retains all autosomes and physical blocks while
using strides 4/2/1, repeated three times at each depth. Input artifacts
are warm cached. No disk spill is allowed.

| depth | stride | repetition |   sites | output_rows | wall_seconds | simulation_seconds | fit_seconds |  rss_bytes | peak_buffer_bytes | spill_bytes |
|------:|-------:|-----------:|--------:|------------:|-------------:|-------------------:|------------:|-----------:|------------------:|------------:|
|   0.1 |      4 |          1 | 1220948 |           2 |         6.12 |              3.622 |       1.947 | 1082208256 |        1017073664 |           0 |
|   0.1 |      4 |          2 | 1220948 |           2 |         6.06 |              3.570 |       1.946 | 1082294272 |        1017135104 |           0 |
|   0.1 |      4 |          3 | 1220948 |           2 |         6.09 |              3.558 |       1.989 | 1082204160 |        1017057280 |           0 |
|   0.1 |      2 |          1 | 2441881 |           2 |         8.60 |              4.079 |       3.939 | 1965490176 |        2064629760 |           0 |
|   0.1 |      2 |          2 | 2441881 |           2 |         8.62 |              4.097 |       3.946 | 1965408256 |        2064711680 |           0 |
|   0.1 |      2 |          3 | 2441881 |           2 |         8.61 |              4.093 |       3.956 | 1965649920 |        2064728064 |           0 |
|   0.1 |      1 |          1 | 4883749 |           2 |        14.21 |              4.408 |       9.108 | 4305403904 |        8377683968 |           0 |
|   0.1 |      1 |          2 | 4883749 |           2 |        13.97 |              4.410 |       8.913 | 4307398656 |        8377638912 |           0 |
|   0.1 |      1 |          3 | 4883749 |           2 |        13.85 |              4.307 |       8.898 | 4309794816 |        8377675776 |           0 |
|   0.4 |      4 |          1 | 1220948 |           2 |         5.92 |              3.461 |       1.945 | 1082728448 |        1017597952 |           0 |
|   0.4 |      4 |          2 | 1220948 |           2 |         5.96 |              3.486 |       1.946 | 1082658816 |        1017659390 |           0 |
|   0.4 |      4 |          3 | 1220948 |           2 |         6.00 |              3.499 |       1.979 | 1082699776 |        1017606144 |           0 |
|   0.4 |      2 |          1 | 2441881 |           2 |         8.57 |              4.022 |       3.999 | 1965842432 |        2065186816 |           0 |
|   0.4 |      2 |          2 | 2441881 |           2 |         8.49 |              3.981 |       3.970 | 1965301760 |        2065256448 |           0 |
|   0.4 |      2 |          3 | 2441881 |           2 |         8.52 |              4.014 |       3.952 | 1965346816 |        2065235968 |           0 |
|   0.4 |      1 |          1 | 4883749 |           2 |        14.17 |              4.347 |       9.182 | 4311035904 |        8378200064 |           0 |
|   0.4 |      1 |          2 | 4883749 |           2 |        14.13 |              4.346 |       9.137 | 4304101376 |        8378167296 |           0 |
|   0.4 |      1 |          3 | 4883749 |           2 |        14.17 |              4.372 |       9.146 | 4310396928 |        8378167296 |           0 |

    ## 
    ## Depth 0.1: endpoint exponent 0.5989; gate met.
    ## 
    ## Depth 0.4: endpoint exponent 0.6247; gate met.

| workload         | observations | max_RSS_GiB | max_buffer_GiB | max_spill_bytes | RSS_ceiling_failures | buffer_ceiling_failures |
|:-----------------|-------------:|------------:|---------------:|----------------:|---------------------:|------------------------:|
| Campaign/probes  |        18000 |      3.0767 |         2.6444 |               0 |                    0 |                       0 |
| Isolated scaling |           18 |      4.0150 |         7.8028 |               0 |                    0 |                       0 |

The nearest checked-in baseline is the chr20 report’s
fixed-site/fixed-group scaling experiment. It has different site counts,
blocks, population-frequency construction and workload; these results do
not establish a performance regression or absence of one. The production
SQL/native path is unchanged.

## Interpretation

Only 30/75 primary strata meet every declared release criterion. At a 5%
share, 0.1x gives positive interval lower limits in 199–200/200
unrelated replicates and 160–165/200 IBD1 replicates; at 0.4x all tested
relationships give 200/200. Those call rates are not a calibrated
detection limit. Five of fifteen primary null strata fail the
false-positive upper-bound criterion, and no primary null stratum meets
every criterion. Bias passes in 40/45 interior strata, point accuracy in
37/45, and coverage in 53/75 primary strata. Ten generating-model fits
lack intervals; no missing intervals or failed statuses are filtered.

Release criteria are assessed only for completed 200-replicate strata,
under the provided frequency anchors and generating relationship.
Detection, point recovery and calibrated coverage are distinct outcomes.
An `ok_conditional` status is a numerical fit status, not a
release-validation verdict. In-panel frequencies, linked finite genomes,
within-superpopulation mismatch and fixed IBD1 haplotypes limit transfer
to independent low-pass datasets. Real-count genome-wide calibration
transfer is not evaluated.
