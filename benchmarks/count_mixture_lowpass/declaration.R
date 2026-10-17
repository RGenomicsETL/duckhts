lowpass_declaration <- list(
  seed = 332202610L,
  populations = c(sample = "CEU", other = "YRI"),
  founders_per_population = 60L,
  replicates = 50L,
  release_replicates = 200L,
  relationships = c("same", "other", "ibd1"),
  shares = c(0, .01, .02, .05, .10, .20),
  depths = c(.1, .2, .4),
  errors = c(.0001, .0002),
  artefact_weight = .002,
  homozygosity_excess = 0,
  block_bases = 1e6,
  min_blocks = 30L,
  max_depth = 8L,
  threads = 1L,
  moment_tolerance = 1e-10,
  bias_limits = c(`0.1` = .005, `0.2` = .0035, `0.4` = .0025),
  error_limits = list(other = c(.020, .014, .010), same = c(.025, .015, .008),
                      ibd1 = c(.050, .025, .015)),
  error_fraction = .90,
  coverage_limits = c(.90, .98),
  null_upper_limit = .075,
  thinning_probabilities = c(1/300, 1/150, 1/75),
  count_samples = c("NA18507", "HG00403", "HG00188"),
  thinning_replicates = 50L,
  thinning_error = .00015,
  scaling_sites = c(100000L, 200000L, 400000L),
  scaling_groups = c(1L, 2L, 4L),
  scaling_repeats = 3L,
  scaling_seconds = 5,
  exponent_limit = 1.25,
  memory_overhead = 512 * 2^20,
  sql_live_bytes_per_row = 768,
  native_bytes_per_block = 1024,
  memory_limit = "3GB",
  spill_allowed = FALSE,
  source_artifacts = c("roh_ancestry_chr20_source", "roh_ancestry_pedigree",
                       "roh_ancestry_chr20_af_sites"),
  count_artifacts = c("roh_counts_chr20_na18507", "roh_counts_chr20_hg00403",
                      "roh_counts_chr20_hg00188")
)

lowpass_strata <- function() {
  x <- expand.grid(relationship = lowpass_declaration$relationships,
                   share = lowpass_declaration$shares,
                   depth = lowpass_declaration$depths,
                   seq_error = lowpass_declaration$errors,
                   stringsAsFactors = FALSE)
  x$id <- seq_len(nrow(x))
  x$seed <- lowpass_declaration$seed + x$id
  x
}

lowpass_thinning_strata <- function() {
  samples <- lowpass_declaration$count_samples
  single <- expand.grid(receiver = samples, probability = lowpass_declaration$thinning_probabilities,
                        stringsAsFactors = FALSE)
  single$donor <- single$receiver
  single$share <- 0
  mixture <- expand.grid(receiver = samples, donor = samples,
                         probability = lowpass_declaration$thinning_probabilities,
                         share = lowpass_declaration$shares[-1L], stringsAsFactors = FALSE)
  mixture <- mixture[mixture$receiver != mixture$donor, ]
  x <- rbind(single[names(mixture)], mixture)
  x$id <- seq_len(nrow(x))
  x
}
