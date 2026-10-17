source("benchmarks/count_mixture_lowpass_genome/batch.R")
args <- commandArgs(TRUE)
extension <- args[[1L]]
manifest_path <- args[[2L]]
root <- args[[3L]]
ids <- as.integer(strsplit(args[[4L]], ",", fixed = TRUE)[[1L]])
key <- digest::digest(file = manifest_path, algo = "sha256")
manifest <- readRDS(manifest_path)
stopifnot(identical(manifest$input_identity, genome_input_identity()))
strata <- genome_strata()
panels <- genome_load_panels()
con <- genome_connect(extension, genome_declaration$threads)
dbExecute(con, "SET max_temp_directory_size = '0B'")
lowpass_moment_macro(con)

process_rss <- function() {
  line <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  as.numeric(strsplit(trimws(sub("VmHWM:", "", line)), "[[:space:]]+")[[1L]][[1L]]) * 1024
}

check_replicate <- function(counts, batched, error) {
  solo <- genome_batch_fit(con, counts, error, genome_declaration$threads)
  batch_exact <- identical(genome_order_fit(batched), genome_order_fit(solo$fit))
  stopifnot(batch_exact)
  one <- genome_batch_fit(con, counts, error, 1L)
  four <- genome_batch_fit(con, counts, error, 4L)
  genome_bulk_counts(con, counts)
  on.exit(dbExecute(con, "DROP TABLE counts"))
  dbExecute(con, "SET threads=1")
  one_moments <- lowpass_moments(con, genome_settings$block_bases)
  oracle <- lowpass_oracle(counts, one_moments, genome_settings$block_bases)
  dbExecute(con, "SET threads=4")
  four_moments <- lowpass_moments(con, genome_settings$block_bases)
  thread <- lowpass_thread_check(one$fit, four$fit, one_moments, four_moments)
  thread$interval_missing_equal <- identical(is.na(one$fit$share_lower), is.na(four$fit$share_lower)) &&
    identical(is.na(one$fit$share_upper), is.na(four$fit$share_upper))
  a <- unlist(one$fit[c("share_lower", "share_upper")], use.names = FALSE)
  b <- unlist(four$fit[c("share_lower", "share_upper")], use.names = FALSE)
  finite <- is.finite(a) & is.finite(b)
  thread$interval_error <- if (any(finite)) lowpass_relative_error(a[finite], b[finite]) else 0
  thread$diagnostic_equal <- identical(one$fit$relation_diagnostic, four$fit$relation_diagnostic)
  thread$bit_identical <- identical(one$fit[genome_exact_fields], four$fit[genome_exact_fields])
  thread$batch_exact <- batch_exact
  thread$peak_buffer_bytes <- max(solo$peak_buffer_bytes, one$peak_buffer_bytes, four$peak_buffer_bytes)
  thread$spill_bytes <- max(solo$spill_bytes, one$spill_bytes, four$spill_bytes)
  stopifnot(oracle$counters_exact, oracle$moment_error <= genome_declaration$oracle_relative_tolerance,
    thread$status_equal, thread$share_missing_equal, thread$interval_missing_equal,
    thread$diagnostic_equal, thread$share_error <= genome_declaration$thread_share_tolerance,
    thread$interval_error <= genome_declaration$thread_share_tolerance,
    thread$thread_moment_error <= genome_declaration$oracle_relative_tolerance)
  list(oracle = oracle, thread = thread)
}

for (id in ids) {
  setting <- strata[id, ]
  batch_size <- genome_declaration$microbatch_replicates[[as.character(setting$depth)]]
  started <- proc.time()[["elapsed"]]
  for (first in seq.int(1L, setting$replicates, by = batch_size)) {
    replicates <- seq.int(first, min(first + batch_size - 1L, setting$replicates))
    seeds <- genome_declaration$seed + id * 1000L + replicates
    paths <- file.path(root, sprintf("stratum%02d-seed%d.rds", id, seeds))
    present <- file.exists(paths)
    if (any(present)) {
      for (path in paths[present]) {
        saved <- readRDS(path)
        stopifnot(identical(saved$manifest_sha256, key))
      }
    }
    if (all(present)) next
    simulation_start <- proc.time()[["elapsed"]]
    simulated <- lapply(replicates, function(replicate) genome_counts(panels, setting, replicate))
    counts <- do.call(rbind, lapply(simulated, `[[`, "counts"))
    simulation_seconds <- proc.time()[["elapsed"]] - simulation_start
    live_bytes <- as.numeric(object.size(panels)) + as.numeric(object.size(simulated)) +
      as.numeric(object.size(counts))
    batch_start <- proc.time()[["elapsed"]]
    batch <- genome_batch_fit(con, counts, setting$error, genome_declaration$threads)
    stopifnot(batch$peak_buffer_bytes < 6e9, batch$spill_bytes == 0)
    checked <- NULL
    if (first == 1L) {
      sample_id <- simulated[[1L]]$counts$sample_id[[1L]]
      checked <- check_replicate(simulated[[1L]]$counts,
        batch$fit[batch$fit$sample == sample_id, ], setting$error)
    }
    batch_seconds <- proc.time()[["elapsed"]] - batch_start + simulation_seconds
    for (i in seq_along(replicates)) {
      x <- simulated[[i]]
      fit <- batch$fit[batch$fit$sample == x$counts$sample_id[[1L]], ]
      rownames(fit) <- NULL
      stopifnot(all(fit$sites == x$qc[["sites"]] - x$qc[["zero_depth_sites"]]),
        all(fit$zero_depth_sites == 0), all(fit$excluded_sites == x$qc[["excluded_sites"]]),
        all(fit$unavailable_sites == x$qc[["unavailable_sites"]]))
      result <- list(
        setting = as.list(setting[c("depth", "share", "relationship", "error")]),
        seed = x$seed, identity = x$identity, fit = fit, qc = x$qc,
        simulation_seconds = simulation_seconds / length(replicates),
        loading_seconds = batch$loading_seconds / length(replicates),
        fit_seconds = batch$seconds / length(replicates),
        batch_wall_seconds = batch_seconds, batch_replicates = length(replicates),
        batch_input_rows = nrow(counts), batch_first_seed = seeds[[1L]],
        process_rss_bytes = process_rss(), peak_buffer_bytes = batch$peak_buffer_bytes,
        spill_bytes = batch$spill_bytes, live_bytes = live_bytes,
        threads = genome_declaration$threads, stride = 1L,
        oracle = if (replicates[[i]] == 1L) checked$oracle else NULL,
        thread_check = if (replicates[[i]] == 1L) checked$thread else NULL,
        manifest_sha256 = key
      )
      if (!present[[i]]) {
        saveRDS(result, paste0(paths[[i]], ".tmp"), compress = "gzip")
        stopifnot(file.rename(paste0(paths[[i]], ".tmp"), paths[[i]]))
      }
    }
    rm(simulated, counts, batch, checked, result, x, fit)
    gc()
  }
  line <- sprintf("Completed genome stratum %d/%d: depth=%g share=%g relationship=%s error=%g; 200 replicates, %.1f seconds.",
    id, nrow(strata), setting$depth, setting$share, setting$relationship, setting$error,
    proc.time()[["elapsed"]] - started)
  cat(line, "\n")
  progress <- Sys.getenv("SOL_D_PROGRESS")
  if (nzchar(progress)) cat(line, "\n", file = progress, append = TRUE)
  flush.console()
}
dbDisconnect(con, shutdown = TRUE)
