# Genome-wide declaration for the low-pass second-genome campaign.
# Commit this file before running any replicate that it judges.
genome_declaration <- list(
  seed = 332202611L,
  shares = c(0, .02, .05, .10, .20),
  depths = c(.1, .4, .2),
  relationships = c("same", "other", "ibd1"),
  errors = c(.0001, .0002),
  artefact_weight = .002,
  replicates = 200L,
  workers = 2L,
  threads = 7L,
  cores = "0-5,12-19",
  wall_budget_hours = 1.5,
  memory_limit = "6GB",
  block_bases = 10000000L,
  min_blocks = 30L,
  max_depth = 8L,
  max_state_bytes = 268435456L,
  spill_bytes = 0,
  memory_overhead_bytes = 768 * 1024^2,
  sql_live_bytes_per_site = 512,
  live_memory_multiplier = 3,
  scaling_strides = c(4L, 2L, 1L),
  scaling_max_exponent = 1.25,
  scaling_min_seconds = 5,
  oracle_relative_tolerance = 1e-10,
  thread_share_tolerance = 1e-5,
  primary_depths = c(.1, .4),
  secondary_depths = .2,
  additional_primary_error = .0002,
  bias_limits = c(`0.1` = .005, `0.2` = .0035, `0.4` = .0025),
  accuracy_limits = list(
    same = c(`0.1` = .025, `0.2` = .015, `0.4` = .008),
    other = c(`0.1` = .020, `0.2` = .014, `0.4` = .010),
    ibd1 = c(`0.1` = .050, `0.2` = .025, `0.4` = .015)),
  minimum_accuracy = .90,
  minimum_coverage = .92,
  nominal_coverage = .95,
  null_false_positive_upper = .075,
  matched_diagnostic = "unresolved"
)

genome_strata <- function() {
  grid <- expand.grid(share = genome_declaration$shares,
    relationship = genome_declaration$relationships,
    depth = genome_declaration$depths, error = genome_declaration$errors,
    stringsAsFactors = FALSE)
  grid$stratum <- seq_len(nrow(grid))
  grid$primary <- grid$depth %in% genome_declaration$primary_depths |
    grid$error == genome_declaration$additional_primary_error
  grid$replicates <- genome_declaration$replicates
  grid
}

# Bias and 90% absolute-error gates apply to correctly specified interior
# shares >= 5%; coverage applies to every stratum. The coverage binomial 95%
# interval must include 95%. Null claims require share_lower > 0; their exact
# one-sided 95% binomial upper bound must be <= .075. Missing/failed fits count
# against accuracy and coverage, and prevent a recovery pass. No stratum with
# fewer than 200 completed replicates is evaluated against release targets.
# Report every status, diagnostic, interval mean/median/90th percentile/max,
# full-range interval frequency, detection rate and binomial uncertainty.
# Matched-frequency rows, including IBD1, must remain unresolved; scores do not
# select relationships. Inspect both output models, retain all replicates.
# Check independent R moments for replicate 1 in each stratum, plus 1/4-thread
# statuses, missingness, diagnostic, shares, moments and interval endpoints.
# Memory ceiling: overhead + 3 * (512 bytes * sites + largest live R panel/count
# pair + native 256 MiB cap), applied separately to RSS and reported DuckDB
# peak buffer. Fixed-block 1x/2x/4x timings use strides 4/2/1 across all autosomes.
# Test the endpoint exponent only if the smallest end-to-end timing is >= 5 s.
# Campaign execution uses long-lived workers, one input-identity verification
# at startup, and one panel load per chromosome per worker. All 200 replicates
# of a stratum share one bulk counts relation and one macro invocation, with
# sample_id '<stratum>-<replicate>'. Simulation remains independent of the
# fitting implementation; each replicate's seed is derived from campaign seed,
# stratum and replicate. Sampled solo refits must be bit-identical to batch rows.
# Independent moments and one/four-thread checks retain the tolerances above.
# Fresh-process timing and memory measurements belong to the separate scaling
# series. First measure one batched stratum; the campaign target is 1.5 hours
# on the declared cores, not an assumed measured throughput.
# The existing 45 strata, their 200 replicates, simulation model, calibration
# and release criteria are unchanged. The same 45-stratum grid at e=2e-4 is
# additionally primary, with 200 replicates each and identical release gates.
# Total: 18,000 replicates, 90 strata, of which 75 are primary.
