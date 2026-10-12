# The declared design of the planted-truth and read-mixture evaluation of the
# ROH decodes: inputs, simulation, arms, metrics, comparisons, criteria and
# gates. benchmark_roh_truth_and_mixtures.Rmd states the same values in prose
# and reads them from here. They are fixed before any arm is decoded.
#
# Coordinates are one-based and inclusive.
roh_truth_declaration <- list(
  # Registered inputs.
  source_artifact = "roh_ancestry_chr20_source",
  pedigree_artifact = "roh_ancestry_pedigree",
  sites_artifact = "roh_ancestry_chr20_af_sites",
  # Children: each population gives `children` children. Founders are the
  # pedigree rows with FatherID and MotherID "0" of that population, in
  # SampleID order. For each population in this order, sample.int draws
  # 3 x children founders without replacement: the first third are the
  # first parents (A), the next third the second parents (B), the last third
  # the second genomes of the read counts (C). Child i takes one haplotype of
  # A_i and one of B_i, each drawn with chance 1/2, so A_i and B_i are two
  # individuals and no two children share a founder.
  populations = c("YRI", "CEU", "CHS", "GIH", "PEL"),
  children = 20L,
  reduced_populations = "CEU",
  reduced_children = 4L,
  seed = 20261012L,
  # Planted autozygous segments, drawn per child in this order of lengths.
  # A start is uniform over the positions that keep the segment between the
  # first and the last panel site. A draw is rejected when the segment
  # overlaps a window of window_bases (floor((pos - 1) / window_bases)) that
  # holds fewer than min_window_sites panel sites, or when it lies within
  # min_gap_bases of a segment already planted in the child. Inside a
  # segment the second haplotype is replaced by the first, so the child is
  # homozygous for A_i's haplotype there.
  planted_lengths = c(5e6, 3e6, 2e6, 1e6),
  min_gap_bases = 1e6,
  window_bases = 1e5,
  min_window_sites = 20L,
  max_draws = 10000L,
  # Frequencies: the ALT allele frequency among all founders of the child's
  # population at each panel site, clamped to [af_clamp, 1 - af_clamp] as the
  # staged read counts clamp theirs.
  af_clamp = 1e-3,
  # Genotype evidence: each site of each child independently, with chance
  # gt_error_rate, reports one of the two other genotypes, each with chance
  # 1/2. Unphased GT.
  gt_error_rate = 0.01,
  # Read evidence: depth is Poisson with mean `depth` at every site. With
  # chance artefact_share a site's alt count is uniform on 0..depth. Otherwise
  # it is binomial with chance (1 - contamination) q(g) + contamination q(g_C),
  # g the child's and g_C founder C's ALT dosage, q = seq_error,
  # allele_balance, 1 - seq_error for 0, 1, 2 copies. ref_count is the rest.
  depth = 30,
  truth = c(seq_error = 5e-4, contamination = 0.02, allele_balance = 0.48,
            artefact_weight = 0.003),
  # Arms, all on the same children and panel sites.
  #   gt30, gt100   duckhts_roh_af_table, the child's GT with planted errors,
  #                 the population frequencies, gt_error := 30 and 100; one
  #                 query per population and gt_error.
  #   counts_default duckhts_roh_counts with its scalar defaults.
  #   counts_truth  duckhts_roh_counts_fitted with the four simulated numbers
  #                 and artefact_weight_sd 0.
  #   counts_fit    duckhts_roh_counts_fitted with each child's row of
  #                 duckhts_count_error_fit (default options, every child in
  #                 one query), its artefact_weight_sd included, so its runs
  #                 carry stable_fraction and gained_bases. A child whose fit
  #                 has a NULL among the four numbers is not decoded and
  #                 counts as a child with no runs in this arm.
  #   fit_low, fit_high counts_fit at artefact_weight w - sd and w + sd (kept
  #                 in [0, 0.05]) with artefact_weight_sd 0: the decodes whose
  #                 runs locate the bases stable_fraction and gained_bases
  #                 count.
  arms = c("gt30", "gt100", "counts_default", "counts_truth", "counts_fit"),
  gt_errors = c(gt30 = 30, gt100 = 100),
  # Metrics per child and arm. Planted bases are the bases of the planted
  # segments. A called base is a base of a run of the arm; a true called base
  # is also planted, a false one is not. sensitivity = true called bases /
  # planted bases; false discovery = false called bases / called bases (NA
  # with no called base). Pooled over children: the sums before the ratio.
  # The criteria use the runs of at least min_call_bases; every run is
  # reported too.
  min_call_bases = 5e5,
  # Paired bootstrap over children, resampled within each population, with
  # set.seed(seed + 1): bootstrap draws of the pooled difference of two arms,
  # 95% percentile interval.
  bootstrap_draws = 2000L,
  bootstrap_level = 0.95,
  # Comparisons, first arm minus second, with intervals.
  comparisons = data.frame(
    test = c("counts_fit", "counts_fit", "gt30", "counts_fit"),
    reference = c("counts_default", "counts_truth", "gt100", "gt30"),
    criterion = c("T1", "T2", "", ""),
    stringsAsFactors = FALSE),
  # T1, the artefact and fit terms help: on runs of at least min_call_bases,
  # the interval of counts_fit - counts_default lies above 0 for sensitivity
  # and below t1_max_fdr_increase for false discovery.
  t1_max_fdr_increase = 0.01,
  # T2, the fit decodes as well as the truth: both intervals of
  # counts_fit - counts_truth lie inside [-t2_margin, t2_margin].
  t2_margin = 0.01,
  # T3, the fit recovers the simulated numbers: in every population, the
  # median over its children of each fitted number lies within its tolerance
  # of the simulated value, and at least t3_min_ok of all children have
  # status 'ok'.
  t3_tolerance = c(seq_error = 2e-4, contamination = 3e-3, allele_balance = 3e-3,
                   artefact_weight = 1e-3),
  t3_min_ok = 0.9,
  # DuckDB threads of every decode and of the fit.
  threads = 2L)

roh_mixture_declaration <- list(
  # Two registered 30x samples of different populations; no related pair is
  # staged. The slices are chr20 CRAMs of the registered public CRAMs.
  receiver = "HG00403",
  donor = "NA18507",
  slice_artifacts = c(HG00403 = "roh_mixture_hg00403_chr20_cram",
                      NA18507 = "roh_mixture_na18507_chr20_cram"),
  # The read counts of each slice alone must reproduce these registered
  # counts, which also give the counting settings.
  count_artifacts = c(HG00403 = "roh_counts_chr20_hg00403",
                      NA18507 = "roh_counts_chr20_na18507"),
  sites_artifact = "roh_ancestry_chr20_af_sites",
  reference_artifact = "liftover_grch38_fasta",
  # Read mixtures. A template (read pair) of the receiver is kept with chance
  # 1 - alpha and one of the donor with chance alpha, decided by samtools
  # view --subsample with --subsample-seed 0, which hashes the read name
  # (samtools 1.23 sam_view.c lines 177-179), so mates stay together. The two
  # kept sets are merged into one coordinate-sorted, indexed CRAM. alpha 0 is
  # the receiver slice alone.
  alphas = c(0.01, 0.02, 0.05, 0.10),
  reduced_alphas = 0.05,
  subsample_seed = 0L,
  worker_count = 4L,
  count_threads = 4L,
  # Decodes of each mixture and of alpha 0: duckhts_roh_counts_fitted with the
  # mixture's fit (duckhts_count_error_fit, default options, all mixtures in
  # one query; artefact_weight_sd 0), with the fitted contamination
  # ("with_term") and with contamination 0 and the other three numbers fitted
  # ("without_term"). FROH is run bases over the span of the panel sites, all
  # runs. Erosion is FROH of a decode of the mixture minus FROH of the
  # with_term decode of alpha 0.
  #
  # M1: the fitted contamination increases strictly with alpha over 0 and the
  # declared alphas.
  # M2: at every declared alpha, |erosion with_term| <= |erosion without_term|
  # + m2_tolerance.
  m2_tolerance = 0.002,
  # Before each mixture is written, the free space of its filesystem minus
  # twice the bytes of the two slices must stay at least min_free_bytes.
  min_free_bytes = 8 * 2^30,
  threads = 2L)

# Gates. The render stops when one fails.
#   - Every registered input has its registered identity: the source VCF, the
#     pedigree and the panel sites their SHA-256, the reference its bytes, each
#     slice its alignments, position sum and base sum.
#   - The children are what the plan says: inside every planted segment no
#     true genotype is heterozygous, and outside them every true genotype is
#     the sum of the two drawn founder haplotypes.
#   - The read counts of each slice alone equal its registered counts.
#   - Each DuckDB query runs in a fresh CLI process under /usr/bin/time with
#     the JSON profiler and max_temp_directory_size 0: none spills, and none
#     has a peak RSS above overhead + 3 x its live state, the overhead being
#     the median peak RSS of three empty runs. Live state, in bytes:
#       count decode     rows x 24 + threads x sites x 36.25
#                        (benchmark_roh_scaling.md, also for the fitted macro)
#       genotype decode  samples x sites x 16 + sites x 128
#                        + threads x sites x 36.25
#                        (duckhts_roh_af_table, the same report)
#       fit              cells x 312 + pairs x 88 + samples x 1001 x 24
#                        (benchmark_count_error_fit.md), cells and pairs
#                        counted untimed from the input
#       read counting    worker_count x (max_region_bytes
#                        + reference_cache_bytes + remote_cache_bytes
#                        + remote_block_bytes + the longest reference contig
#                        of the panel) + panel sites x 512, from the per-job
#                        bounds in the catalog entry of
#                        duckhts_somalier_bam_counts at its defaults, the
#                        contig being the CRAM decoder's reference copy
#   - The bases that stable_fraction and gained_bases count equal the bases
#     that fit_low and fit_high locate: for each run of counts_fit,
#     length x (1 - stable_fraction) differs from its located bases by less
#     than 0.5, and the summed gained_bases equal the located gained bases.
roh_truth_ceiling_bytes <- list(
  count_decode = function(rows, sites, threads) rows * 24 + threads * sites * 36.25,
  genotype_decode = function(samples, sites, threads) {
    samples * sites * 16 + sites * 128 + threads * sites * 36.25
  },
  fit = function(cells, pairs, samples) cells * 312 + pairs * 88 + samples * 1001 * 24,
  read_counting = function(worker_count, contig_bytes, sites,
                           max_region_bytes = 67108864, reference_cache_bytes = 67108864,
                           remote_cache_bytes = 67108864, remote_block_bytes = 1048576) {
    worker_count * (max_region_bytes + reference_cache_bytes + remote_cache_bytes +
                      remote_block_bytes + contig_bytes) + sites * 512
  })
