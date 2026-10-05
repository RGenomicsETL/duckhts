# The declared design of the read-count validation: inputs, metrics, tolerance
# and seed. benchmark_roh_counts_validation.Rmd states the same values in prose
# and reads them from here. They are fixed before any run is decoded.
roh_validation_declaration <- list(
  samples = c("NA18507", "HG00403", "HG00188"),
  count_artifacts = c(NA18507 = "roh_counts_chr20_na18507",
                      HG00403 = "roh_counts_chr20_hg00403",
                      HG00188 = "roh_counts_chr20_hg00188"),
  genotype_artifact = "roh_counts_validation_chr20_genotypes",
  # Phred genotype error of the VCF path, as bcftools roh -G 30.
  gt_error = 30,
  # A long run has at least this many bases (end - start + 1).
  long_run_bases = 1e6,
  # Tolerance: absolute FROH difference, and minimum Jaccard index of run bases.
  max_froh_difference = 0.01,
  min_jaccard_long = 0.9,
  min_jaccard_all = 0.7,
  # Titration: contaminant read fractions and the base seed. Cell k of the
  # (receiver, contaminant, alpha) grid, in the order run.R builds it, draws
  # with set.seed(seed + k).
  alphas = c(0.01, 0.02, 0.05, 0.10),
  seed = 20261004L)
