source("benchmarks/count_mixture_lowpass/common.R")
args <- commandArgs(TRUE)
scratch <- args[[1L]]
baseline <- args[[2L]]
probe <- args[[3L]]
output <- "benchmarks/count_mixture_lowpass/results"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

write_result <- function(x, name, compressed = FALSE) {
  path <- file.path(output, name)
  if (compressed) {
    connection <- gzfile(path, "wt", compression = 9L)
    on.exit(close(connection))
    write.csv(x, connection, row.names = FALSE)
  } else write.csv(x, path, row.names = FALSE)
}

load_campaign <- function(path) {
  files <- list.files(path, pattern = "^(synthetic|thinning)-[0-9]+\\.rds$", full.names = TRUE)
  expected <- c(sprintf("synthetic-%d.rds", lowpass_strata()$id),
                 sprintf("thinning-%d.rds", lowpass_thinning_strata()$id))
  stopifnot(setequal(basename(files), expected))
  fits <- checks <- strata <- vector("list", length(files))
  for (i in seq_along(files)) {
    x <- readRDS(files[[i]])
    kind <- x$checks$kind
    fits[[i]] <- cbind(kind = kind, x$fits)
    checks[[i]] <- x$checks
    strata[[i]] <- cbind(kind = kind, x$stratum)
  }
  list(fits = do.call(rbind, fits), checks = do.call(rbind, checks), strata = strata)
}

campaign <- load_campaign(scratch)
old <- load_campaign(baseline)
write_result(campaign$checks, "oracle_and_threads.csv")
write_result(old$checks, "compensated_sum_counterexamples.csv")
# All conditional models and statuses are retained; fixed calibration and
# method strings are in the declaration rather than repeated on every row.
fit_columns <- c("kind", "id", "replicate", "primary", "relationship_model", "share", "share_lower",
  "share_upper", "share_se", "moment_score", "score_rank", "jackknife_failed_refits",
  "at_bound", "interval_hull", "status", "relation_diagnostic")
write_result(campaign$fits[fit_columns], "replicates.csv.gz", TRUE)
qc_columns <- c("kind", "id", "replicate", "sites", "usable_sites", "zero_depth_sites", "singleton_sites",
  "reads", "other_reads", "positive_pair_sites", "excluded_sites", "unavailable_sites", "blocks", "pair_blocks",
  "concordant_ref_pairs", "discordant_pairs", "concordant_alt_pairs", "pair_denominator",
  "singles_share", "singles_se", "pair_residual", "pair_residual_se", "discordance_ratio", "max_block_leverage")
qc <- unique(campaign$fits[qc_columns])
stopifnot(!anyDuplicated(qc[c("kind", "id", "replicate")]))
write_result(qc, "replicate_qc.csv.gz", TRUE)

summaries <- lapply(campaign$strata, function(stratum) {
  kind <- stratum$kind
  x <- campaign$fits[campaign$fits$kind == kind & campaign$fits$id == stratum$id & campaign$fits$primary, ]
  n <- nrow(x)
  expected_replicates <- if (kind == "synthetic") lowpass_declaration$replicates else lowpass_declaration$thinning_replicates
  stopifnot(n == expected_replicates, identical(sort(x$replicate), seq_len(expected_replicates)))
  truth <- stratum$share
  finite <- is.finite(x$share)
  conditional <- x$status == "ok_conditional"
  intervals <- is.finite(x$share_lower) & is.finite(x$share_upper)
  coverage <- sum(intervals & x$share_lower <= truth & x$share_upper >= truth, na.rm = TRUE) / n
  claimed <- sum(intervals & x$share_lower > 0, na.rm = TRUE)
  null_upper <- if (claimed == n) 1 else qbeta(.95, claimed + 1, n - claimed)
  bias <- if (any(finite)) mean(x$share[finite]) - truth else NA_real_
  relationship <- if (kind == "synthetic") stratum$relationship else "counts"
  depth <- if (kind == "synthetic") stratum$depth else stratum$probability*30
  error <- if (kind == "synthetic") stratum$seq_error else lowpass_declaration$thinning_error
  depth_index <- match(depth, lowpass_declaration$depths)
  limit <- if (kind == "synthetic") lowpass_declaration$error_limits[[relationship]][[depth_index]] else NA_real_
  bias_limit <- if (kind == "synthetic") lowpass_declaration$bias_limits[[depth_index]] else NA_real_
  accurate <- if (kind == "synthetic") sum(conditional & finite & abs(x$share-truth) <= limit, na.rm = TRUE)/n else NA_real_
  bias_verdict <- if (kind != "synthetic" || truth < .05) "not applicable" else
    if (all(conditional & finite) && abs(bias) <= bias_limit) "met" else "fail"
  error_verdict <- if (kind != "synthetic" || truth == 0) "not applicable" else
    if (accurate >= lowpass_declaration$error_fraction) "met" else "fail"
  coverage_verdict <- if (kind != "synthetic" || truth == 0) "not applicable" else
    if (coverage >= lowpass_declaration$coverage_limits[[1L]] &&
        coverage <= lowpass_declaration$coverage_limits[[2L]]) "met" else "fail"
  null_verdict <- if (kind != "synthetic" || truth != 0) "not applicable" else
    if (null_upper > lowpass_declaration$null_upper_limit) "fail" else
      if (!all(intervals)) "incomplete" else "met"
  data.frame(kind = kind, id = stratum$id, relationship = relationship, depth = depth, seq_error = error,
    truth = truth, receiver = if (kind == "thinning") stratum$receiver else "CEU",
    donor = if (kind == "thinning") stratum$donor else relationship,
    replicates = n, conditional = sum(conditional), finite_shares = sum(finite), intervals = sum(intervals),
    mean_share = if (any(finite)) mean(x$share[finite]) else NA_real_, bias = bias,
    mean_absolute_error = if (any(finite)) mean(abs(x$share[finite]-truth)) else NA_real_,
    singles_mean = if (any(is.finite(x$singles_share))) mean(x$singles_share[is.finite(x$singles_share)]) else NA_real_,
    singles_bias = if (any(is.finite(x$singles_share))) mean(x$singles_share[is.finite(x$singles_share)])-truth else NA_real_,
    singles_within_limit_fraction = if (kind == "synthetic" && relationship == "other")
      sum(is.finite(x$singles_share) & abs(x$singles_share-truth) <= limit, na.rm = TRUE)/n else NA_real_,
    accurate_fraction = accurate, interval_coverage = coverage, positive_claims = claimed,
    positive_claim_fraction = claimed/n, null_upper_95 = if (truth == 0) null_upper else NA_real_,
    bias_observed = bias_verdict, accuracy_observed = error_verdict, coverage_observed = coverage_verdict,
    null_observed = null_verdict, release_targets = "not evaluated: chr20, 50 < 200 replicates")
})
summary <- do.call(rbind, summaries)
summary <- summary[order(summary$kind, summary$id), ]
write_result(summary, "strata.csv")

scaling_files <- list.files(scratch, pattern = "^scale-.*\\.rds$", full.names = TRUE)
scaling <- do.call(rbind, lapply(scaling_files, function(file) {
  x <- readRDS(file)
  fields <- strsplit(sub("\\.rds$", "", basename(file)), "-", fixed = TRUE)[[1L]]
  live <- lowpass_declaration$sql_live_bytes_per_row*x$input_rows +
    lowpass_declaration$native_bytes_per_block*ceiling(x$sites/1e6)*x$groups
  ceiling <- lowpass_declaration$memory_overhead + 3*live
  rss <- as.numeric(readLines(paste0(file, ".rss")))*1024
  data.frame(axis = fields[[2L]], factor = as.integer(fields[[3L]]), repeat_id = as.integer(fields[[4L]]),
    sites = x$sites, samples = x$samples, candidates = x$candidates, input_rows = x$input_rows,
    output_rows = x$output_rows, threads = x$threads, seconds = x$seconds, rss_bytes = rss,
    buffer_bytes = x$buffer_bytes, spill_bytes = x$temp_bytes, live_bytes = live, ceiling_bytes = ceiling,
    rss_ok = rss <= ceiling, buffer_ok = x$buffer_bytes <= ceiling, spill_ok = x$temp_bytes == 0)
}))
expected_scaling <- expand.grid(axis = c("sites", "groups"), factor = c(1L, 2L, 4L),
                                 repeat_id = seq_len(lowpass_declaration$scaling_repeats))
scale_key <- function(x) sprintf("%s-%d-%d", x$axis, x$factor, x$repeat_id)
stopifnot(nrow(scaling) == nrow(expected_scaling), !anyDuplicated(scale_key(scaling)),
          setequal(scale_key(scaling), scale_key(expected_scaling)))
write_result(scaling, "scaling.csv")
scale_summary <- do.call(rbind, lapply(split(scaling, interaction(scaling$axis, scaling$factor)), function(x) {
  data.frame(axis = x$axis[[1L]], factor = x$factor[[1L]], input_rows = x$input_rows[[1L]],
    output_rows = x$output_rows[[1L]], median_seconds = median(x$seconds),
    min_seconds = min(x$seconds), max_seconds = max(x$seconds),
    max_rss_MiB = max(x$rss_bytes)/2^20, max_buffer_MiB = max(x$buffer_bytes)/2^20,
    ceiling_MiB = x$ceiling_bytes[[1L]]/2^20, spill_MiB = max(x$spill_bytes)/2^20)
}))
scale_summary <- scale_summary[order(scale_summary$axis, scale_summary$factor), ]
scale_summary$exponent <- NA_real_
scale_summary$exponent_verdict <- "not evaluated: 1x < 5 seconds"
for (axis in unique(scale_summary$axis)) {
  rows <- which(scale_summary$axis == axis)
  times <- scale_summary$median_seconds[rows]
  scale_summary$exponent[rows[-1L]] <- log2(times[-1L]/times[-length(times)])
  if (times[[1L]] >= lowpass_declaration$scaling_seconds) {
    scale_summary$exponent_verdict[rows[-1L]] <- ifelse(
      scale_summary$exponent[rows[-1L]] <= lowpass_declaration$exponent_limit, "met", "fail")
  }
}
write_result(scale_summary, "scaling_summary.csv")

memory_files <- list.files(probe, pattern = "^(ordered|unordered|materialized)-[0-9]+\\.rds$", full.names = TRUE)
memory <- do.call(rbind, lapply(memory_files, function(file) {
  x <- readRDS(file)
  data.frame(implementation = sub("-[0-9]+\\.rds$", "", basename(file)), sites = x$sites,
    threads = x$threads, seconds = x$seconds,
    rss_MiB = as.numeric(readLines(sub("\\.rds$", ".rss", file)))/1024,
    buffer_MiB = x$buffer_bytes/2^20, spill_bytes = x$temp_bytes)
}))
write_result(memory, "exploratory_memory.csv")

checks <- campaign$checks
criteria <- data.frame(criterion = c("R moment oracle", "Exact integer counters", "1/4-thread statuses",
  "1/4-thread shares", "1/4-thread moments", "Scaling RSS ceiling", "Scaling DuckDB buffer ceiling", "No scaling spill",
  "Timing exponent <= 1.25", "Genome-wide release recovery/interval targets"),
  denominator = c(rep(nrow(checks), 5), rep(nrow(scaling), 3), 2, 108),
  verdict = c(if (all(checks$moment_error <= lowpass_declaration$moment_tolerance)) "met" else "fail",
    if (all(checks$counters_exact)) "met" else "fail",
    if (all(checks$status_equal & checks$share_missing_equal)) "met" else "fail",
    if (all(checks$share_error <= lowpass_declaration$moment_tolerance)) "met" else "fail",
    if (all(checks$thread_moment_error <= lowpass_declaration$moment_tolerance)) "met" else "fail",
    if (all(scaling$rss_ok)) "met" else "fail", if (all(scaling$buffer_ok)) "met" else "fail",
    if (all(scaling$spill_ok)) "met" else "fail", "not evaluated: 1x < 5 seconds",
    "not evaluated: chr20, 50 < 200 replicates"))
write_result(criteria, "summary.csv")
empty_files <- list.files(scratch, pattern = "^empty-[0-9]+\\.rds$", full.names = TRUE)
if (length(empty_files)) {
  empty <- do.call(rbind, lapply(empty_files, function(file) {
    x <- readRDS(file)
    data.frame(rss_MiB = as.numeric(readLines(sub("\\.rds$", ".rss", file)))/1024,
               buffer_bytes = x$buffer_bytes, seconds = x$seconds)
  }))
  write_result(empty, "empty_process.csv")
}
metadata <- readRDS(file.path(scratch, "metadata.rds"))
creation <- system2("git", c("log", "--follow", "--diff-filter=A", "--format=%H", "--",
  "benchmarks/count_mixture_lowpass/declaration.R"), stdout = TRUE)
stopifnot(length(creation) == 1L)
original <- system2("git", c("rev-parse", metadata$declaration_commit), stdout = TRUE)
difference <- system2("git", c("diff", original, creation, "--",
  "benchmarks/count_mixture_lowpass/declaration.R"), stdout = TRUE)
stopifnot(length(difference) == 0L)
declared <- new.env(parent = baseenv())
connection <- textConnection(system2("git", c("show", paste0(creation,
  ":benchmarks/count_mixture_lowpass/declaration.R")), stdout = TRUE))
source(connection, local = declared)
close(connection)
stopifnot(identical(lowpass_declaration, declared$lowpass_declaration))
metadata$declaration_recorded_commit <- original
metadata$declaration_commit <- creation
metadata_lines <- trimws(capture.output(dput(metadata)), which = "right")
writeLines(metadata_lines, file.path(output, "metadata.txt"))
cat("Results bytes:", sum(file.info(list.files(output, full.names = TRUE))$size), "\n")
cat("Primary observed threshold counts:\n")
for (name in c("bias_observed", "accuracy_observed", "coverage_observed", "null_observed")) print(table(summary[[name]]))
