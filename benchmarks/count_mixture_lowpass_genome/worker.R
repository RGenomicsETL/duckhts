source("benchmarks/count_mixture_lowpass_genome/common.R")
args <- commandArgs(TRUE)
extension <- args[[1L]]
output <- args[[2L]]
setting <- list(
  depth = as.numeric(args[[3L]]), share = as.numeric(args[[4L]]),
  relationship = args[[5L]], error = as.numeric(args[[6L]])
)
seed <- as.integer(args[[7L]])
check <- as.logical(args[[8L]])
threads <- as.integer(args[[9L]])
stride <- as.integer(args[[10L]])
manifest <- readRDS(args[[11L]])
stopifnot(identical(manifest$input_identity, genome_input_identity()))
con <- genome_connect(extension, threads)
dbExecute(con, paste(
  "CREATE TABLE counts(sample_id VARCHAR, candidate VARCHAR, chrom VARCHAR,",
  "pos BIGINT, ref_count INTEGER, alt_count INTEGER, other_count INTEGER,",
  "sample_af DOUBLE, candidate_af DOUBLE)"
))
set.seed(seed)
first <- readRDS(duckhts_bench_artifact_path(genome_id(1L)))
identity <- genome_draw_identity(first$selected, setting$relationship)
rm(first)
retained <- list()
started <- proc.time()[["elapsed"]]
live_bytes <- 0
for (chrom in genome_settings$chromosomes) {
  path <- duckhts_bench_artifact_path(genome_id(chrom))
  stopifnot(identical(digest::digest(file = path, algo = "sha256"), manifest$artifacts$sha256[[chrom]]))
  panel <- readRDS(path)
  stopifnot(identical(panel$settings, genome_settings))
  if (stride > 1L) {
    keep <- seq.int(1L, nrow(panel$sites), by = stride)
    panel$sites <- panel$sites[keep, ]
    panel$genotypes <- lapply(panel$genotypes, `[`, keep)
  }
  counts <- genome_simulate(panel, setting, identity)
  counts$sample_id <- identity$sample
  counts$candidate <- if (setting$relationship == "other") "AFR" else "matched"
  live_bytes <- max(live_bytes, as.numeric(object.size(panel)) + as.numeric(object.size(counts)))
  dbAppendTable(con, "counts", counts)
  if (check) retained[[as.character(chrom)]] <- counts
  rm(panel, counts)
}
simulation_seconds <- proc.time()[["elapsed"]] - started
memory <- dbGetQuery(con, "SELECT sum(memory_usage_bytes) AS bytes FROM duckdb_memory()")$bytes
profile <- tempfile(fileext = ".json")
dbExecute(con, "PRAGMA enable_profiling='json'")
dbExecute(con, paste("SET profiling_output =", dbQuoteString(con, profile)))
started <- proc.time()[["elapsed"]]
fit <- genome_fit(con, setting$error)
fit_seconds <- proc.time()[["elapsed"]] - started
dbExecute(con, "PRAGMA disable_profiling")
metrics <- jsonlite::fromJSON(profile, simplifyVector = FALSE)
unlink(profile)
oracle <- NULL
thread_check <- NULL
if (check) {
  counts <- do.call(rbind, retained)
  lowpass_moment_macro(con)
  moments <- lowpass_moments(con, genome_settings$block_bases)
  oracle <- lowpass_oracle(counts, moments, genome_settings$block_bases)
  if (threads == 1L) {
    dbExecute(con, "SET threads=4")
    dbExecute(con, "PRAGMA enable_profiling='json'")
    dbExecute(con, paste("SET profiling_output =", dbQuoteString(con, profile)))
    four <- genome_fit(con, setting$error)
    four_metrics <- jsonlite::fromJSON(profile, simplifyVector = FALSE)
    four_moments <- lowpass_moments(con, genome_settings$block_bases)
    four_moment_metrics <- jsonlite::fromJSON(profile, simplifyVector = FALSE)
    dbExecute(con, "PRAGMA disable_profiling")
    unlink(profile)
    thread_check <- lowpass_thread_check(fit, four, moments, four_moments)
    thread_check$interval_missing_equal <- identical(is.na(fit$share_lower), is.na(four$share_lower)) &&
      identical(is.na(fit$share_upper), is.na(four$share_upper))
    endpoints <- c("share_lower", "share_upper")
    a <- unlist(fit[endpoints], use.names = FALSE)
    b <- unlist(four[endpoints], use.names = FALSE)
    finite <- is.finite(a) & is.finite(b)
    thread_check$interval_error <- if (any(finite)) lowpass_relative_error(a[finite], b[finite]) else 0
    thread_check$diagnostic_equal <- identical(fit$relation_diagnostic, four$relation_diagnostic)
    thread_check$peak_buffer_bytes <- max(
      four_metrics$system_peak_buffer_memory,
      four_moment_metrics$system_peak_buffer_memory
    )
    thread_check$spill_bytes <- max(
      four_metrics$system_peak_temp_dir_size,
      four_moment_metrics$system_peak_temp_dir_size
    )
  }
}
dbDisconnect(con, shutdown = TRUE)
result <- list(
  setting = setting, seed = seed, identity = identity, fit = fit,
  simulation_seconds = simulation_seconds, fit_seconds = fit_seconds,
  buffer_before_fit = memory, peak_buffer_bytes = metrics$system_peak_buffer_memory,
  spill_bytes = metrics$system_peak_temp_dir_size, live_bytes = live_bytes,
  threads = threads, stride = stride, oracle = oracle, thread_check = thread_check,
  manifest_sha256 = digest::digest(file = args[[11L]], algo = "sha256")
)
saveRDS(result, paste0(output, ".tmp"))
stopifnot(file.rename(paste0(output, ".tmp"), output))
