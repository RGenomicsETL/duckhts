# Copy-number fractions from read depth

Status: current implementation guidance for phase 1 (corrected bins and the two
estimators without a state model); open design for phase 2 (the state model). Public
issue: https://github.com/RGenomicsETL/duckhts/issues/343.

## Purpose

Read depth in fixed bins is the signal that survives low-pass sequencing (about 0.1x).
A share `m` of cells with one more or one fewer copy of a region moves the depth of that
region by `m`/2 relative to two copies. A second genome of the other sex at share `f`
moves the depth of X and Y by a known amount. This lane turns `bam_bin_counts` output
into those fractions as ordinary relations, with the method and its assumptions in each
row. Uses: cell-free DNA, mosaic fractions in low-pass genomes, cohort sex and
whole-chromosome checks, and the share of a second genome of the other sex.

## Decisions

The issue lists seven open choices. They are settled as follows.

1. **Contract.** The phase 2 state model implements the published ichorCNA model
   (Adalsteinsson, Ha, Freeman et al. 2017, https://doi.org/10.1038/s41467-017-00965-y)
   and is checked against `GavinHaLab/ichorCNA` tag v0.6.1 run as a black box. The
   correction step is DuckHTS's own declared method, not a loess match. ichorCNA and
   HMMcopy are GPL-3.0 and the extension is MIT: no upstream code or data is copied.
2. **First deliverable.** Corrected bins and the two estimators that need no state model.
3. **Smoother.** A median for each GC stratum, computed in SQL. No loess kernel. A native
   loess kernel is added only if phase 2 validation shows that the stratum median limits
   agreement with upstream.
4. **Home of the decoder.** Inside DuckHTS first, as a native kernel behind a table macro
   like the ROH decode, with a kernel interface that the segmentation engine of
   https://github.com/RGenomicsETL/duckhts/issues/255 can take over.
5. **Reference set.** Optional. Each estimator has a `_ref` macro taking `reference_table`.
   A per-bin panel median is one SQL aggregate the caller can write; it is not a separate
   function.
6. **Fragment length.** Later. Phase 1 counts read starts only.
7. **Y evidence.** Aligned depth only. Y-specific k-mer counts from unaligned reads are
   later.

## Phase 1 relations

Three table macros. Each reads caller-named relations through `query_table`, so any bin
source composes: `bam_bin_counts`, `fasta_nuc` for GC, `read_bigwig` for mappability,
`read_bed` and interval joins for exclusion lists.

### `duckhts_depth_bins(bins_table, ...)`: corrected bins

Input: one row for each sample and bin with `sample_id`, `chrom`, `start`, `end`
(zero-based, half-open), `reads` (reads in the bin) and `gc` (G+C share of the A/C/G/T
reference bases of the bin, from 0 to 1). An optional `mappability` column (0 to 1) is
used when present.

Model. The expected count of a bin is its length times the median read rate (count per
base) of the training bins of the same sample and GC stratum, scaled so that the median
ratio of the kept positive-read autosomal bins is 1. Empty bins do not set this scale
and retain zero ratios when their bin and stratum rules pass. The corrected ratio is
`count / expected`; it is 1 at the sample's typical autosomal depth, about 0.5 for one
copy and 1.5 for three.

- GC stratum: `floor(gc / gc_step)`; `gc_step` is in `[1/2147483647, 1]`
  so the stratum fits the public integer column at `gc = 1`.
- Training bins: autosomal, mappable, with a GC value and a read rate inside the
  `outlier_quantile` and `1 - outlier_quantile` quantiles of the sample's autosomal rates,
  and at least one read. Only training bins set the stratum medians.
- Kept bins: a count, a GC value, mappability at least `min_mappability` when given, a
  stratum with at least `min_stratum_bins` training bins, and a ratio of at most
  `max_ratio`. Every other rule depends on the bin and not on its count, so a bin with no
  reads (Y in an XX genome) stays kept. `reason` names why a bin is not kept.
- Autosomal means a contig whose `duckhts_contig_key` is not X, Y or MT. The caller
  restricts contigs (for example to the primary assembly) and removes known problem bins
  (centromeres, gaps, pseudo-autosomal regions on X where that matters) with SQL first.

Mappability only filters bins; it is not a second correction. Bins are not normalised
against a reference set here. Because the strata are fitted on the sample itself, a change
that holds a large share of some strata's training bins is partly absorbed: correct a whole
genome, or give the estimators `expected` counts from a reference panel, which they accept
from any source.

### Regional estimators

`duckhts_depth_region_fraction(bins_table, regions_table)` uses the bin model;
`duckhts_depth_region_fraction_ref(bins_table, regions_table, reference_table)`
uses reference samples.

`bins_table` is corrected bins: `sample_id`, `chrom`, `start`, `end`, `reads`,
`expected`, `kept`. `regions_table` has `region`, `chrom`, `start`, `end`; several rows
with one name form one region, and a bin belongs to a region when its midpoint lies in
one of the region's intervals.

Estimator. With `C` and `E` the sums of count and expected over kept bins, the region's
depth ratio is `(C_r / E_r) / (C_b / E_b)`, where `b` is the baseline: kept autosomal
bins outside the region. `shift = 2 * (ratio - 1)` is the share of cells with a one-copy
gain (positive) or loss (negative), relative to two copies, if the baseline is unchanged.
A constitutional one-copy gain gives 1.

Interval without a reference set: quasi-Poisson bins. The variance of a bin count is
`phi` times its mean, with the dispersion `phi` estimated from the baseline bins, so the
variance of a sum of counts is `phi` times that sum (at least one read). Bins are treated
as independent. The interval is the estimate plus or minus 1.959964 standard errors.

With a reference set (corrected bins of reference samples): each reference sample's
shift for the region is computed the same way, leaving out the tested sample when it is
also a reference. The reported shift is the sample's shift minus the reference mean; the
standard error is the reference standard deviation times `sqrt(1 + 1/n)`, and the
interval uses the 97.5% t quantile with `n - 1` degrees of freedom. The reference set
absorbs region-specific bias and correlated noise that the bin model ignores.

### Sex-copy estimators

`duckhts_depth_sex_fraction(bins_table)` uses ideal copy levels;
`duckhts_depth_sex_fraction_ref(bins_table, reference_table)` uses reference samples.

Model. A major genome and a second genome of the other sex with DNA share `f` < 0.5.
With `c_X` and `c_Y` the X and Y copy numbers implied by depth (twice the depth ratio to
the autosomes):

| Major genome | `c_X` | `c_Y` | `f` from X | `f` from Y |
| --- | --- | --- | --- | --- |
| XX | `2 - f` | `f` | `2 - c_X` | `c_Y` |
| XY | `1 + f` | `1 - f` | `c_X - 1` | `1 - c_Y` |

The major genome is XX when `c_X >= 1.5`. The two estimates are combined by inverse
variance; `agreement_z` compares them. An X or Y copy change of the major genome (loss of
Y in part of the cells, X aneuploidy) moves the two estimates differently, and a
disagreement says so.

Reference set: corrected bins with a `sex` column (`XX` or `XY`). XY references set the
depth of one X copy and one Y copy; XX references, when present, set the depth of two X
copies and the Y background of a genome without Y. Without XX references two X copies
are twice one, and the Y background is 0. The tested sample is left out of its own
reference.

## Bounds

- Corrected bins: the input is materialized once, and the quantile and median aggregates
  hold one value per autosomal bin and sample; the output has one row per input row. Joins
  build on samples x GC strata and on samples.
- Region fraction: the narrow input columns are materialized once. A region's intervals
  are merged, and the join matches each kept bin to the merged intervals of its contig;
  the state per sample and region is a few sums.
- Sex fraction: one pass over the input into a few sums per sample.

## Phase 2 contract: the state model

Not implemented. The decoder takes the corrected log2 ratios of one sample's kept bins
in contig and position order and returns:

- **Solutions:** one row per sample and start: normal fraction `n`, ploidy `phi`,
  log-likelihood, share of the genome altered, rank, and status (`ok`, `guard_zero`,
  `close_to_best`). Starts form a grid of normal fractions and ploidies, as upstream.
- **Segments:** one row per sample, solution and segment: interval, copy number, median
  log ratio, number of bins.

Emission: Student's t on the log ratio with mean
`log(((1 - n) c + 2 n) / ((1 - n) phi + 2 n))` for tumour copy number `c` up to
`max_cn` (default 7), in the log base the pinned source uses for its ratios; EM with priors; transitions with self-probability `txn_e`
(default 0.9999999). Guard: when the largest segment has at most 50 bins and at most 5%
of bins are altered, the fraction is set to zero and the status says so. Method version
and the pinned upstream tag go into every row. Validation is declared before any run:
fraction, ploidy and segments against upstream on the same bins with stated tolerances,
a dilution series by pooling reads, samples without changes, and depth at 0.1x, 0.4x and
1x.

## Limits

- A sample without copy-number change gives no fraction; depth measures altered genome.
- Fraction and ploidy trade against each other in the state model; all solutions are
  returned.
- Depth does not see copy-neutral changes. A constitutional copy-number variant looks
  like a change in every cell.
- The bin model treats bins as independent; real libraries have correlated residual
  bias, so intervals without a reference set are optimistic for large regions.
