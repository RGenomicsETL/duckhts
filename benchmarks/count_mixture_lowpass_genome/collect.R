source("benchmarks/count_mixture_lowpass_genome/declaration.R")

genome_summarize <- function(x, setting) {
  n <- nrow(x)
  conditional <- !is.na(x$status) & x$status == "ok_conditional"
  finite <- is.finite(x$share)
  interval <- is.finite(x$share_lower) & is.finite(x$share_upper) & x$share_lower <= x$share_upper
  covered <- interval & x$share_lower <= setting$share & x$share_upper >= setting$share
  detected <- interval & x$share_lower > 0
  error <- x$share - setting$share
  width <- x$share_upper - x$share_lower
  depth <- as.character(setting$depth)
  limit <- genome_declaration$accuracy_limits[[setting$relationship]][[depth]]
  accuracy <- finite & abs(error) <= limit
  bias <- if (any(finite)) mean(error[finite]) else NA_real_
  coverage_ci <- binom.test(sum(covered), n)$conf.int
  null_upper <- if (sum(detected) == n) 1 else qbeta(.95, sum(detected) + 1L, n - sum(detected))
  complete <- n == genome_declaration$replicates
  interior <- setting$share >= .05
  recovery <- complete && all(conditional & finite & interval)
  bias_pass <- recovery && abs(bias) <= genome_declaration$bias_limits[[depth]]
  accuracy_pass <- recovery && mean(accuracy) >= genome_declaration$minimum_accuracy
  coverage_pass <- complete && mean(covered) >= genome_declaration$minimum_coverage &&
    coverage_ci[[1L]] <= genome_declaration$nominal_coverage &&
    coverage_ci[[2L]] >= genome_declaration$nominal_coverage
  null_pass <- complete && all(interval) && null_upper <= genome_declaration$null_false_positive_upper
  diagnostic_pass <- setting$relationship == "other" ||
    all(!is.na(x$relation_diagnostic) & x$relation_diagnostic == genome_declaration$matched_diagnostic)
  quantile_width <- function(p) if (any(interval)) unname(quantile(width[interval], p)) else NA_real_
  data.frame(
    stratum = setting$stratum, primary = setting$primary,
    relationship = setting$relationship, depth = setting$depth, truth = setting$share,
    error = setting$error, replicates = n, evaluated = complete,
    mean_share = if (any(finite)) mean(x$share[finite]) else NA_real_, bias = bias,
    bias_mc_se = if (sum(finite) > 1L) sd(error[finite]) / sqrt(sum(finite)) else NA_real_,
    absolute_error_90 = if (all(finite)) unname(quantile(abs(error), .9)) else NA_real_,
    accuracy_limit = limit, accuracy = mean(accuracy), coverage = mean(covered),
    coverage_lower = coverage_ci[[1L]], coverage_upper = coverage_ci[[2L]],
    detection = mean(detected), null_upper = if (setting$share == 0) null_upper else NA_real_,
    failed_status = sum(!conditional), missing_share = sum(!finite), missing_interval = sum(!interval),
    interval_mean = if (any(interval)) mean(width[interval]) else NA_real_,
    interval_median = quantile_width(.5), interval_90 = quantile_width(.9),
    interval_max = if (any(interval)) max(width[interval]) else NA_real_,
    full_range_intervals = sum(interval & width >= .5),
    bias_pass = if (interior && complete) bias_pass else NA,
    accuracy_pass = if (interior && complete) accuracy_pass else NA,
    coverage_pass = if (complete) coverage_pass else NA,
    null_pass = if (setting$share == 0 && complete) null_pass else NA,
    diagnostic_pass = diagnostic_pass,
    release_pass = if (complete) {
      recovery && diagnostic_pass && coverage_pass &&
        (!interior || (bias_pass && accuracy_pass)) && (setting$share != 0 || null_pass)
    } else {
      NA
    }
  )
}

genome_log_metrics <- function(path) {
  lines <- readLines(path)
  rss <- as.numeric(sub(".*: ", "", grep("Maximum resident set size", lines, value = TRUE))) * 1024
  elapsed <- sub(".*\\): ", "", grep("Elapsed \\(wall clock\\)", lines, value = TRUE))
  parts <- as.numeric(strsplit(elapsed, ":", fixed = TRUE)[[1L]])
  wall <- sum(rev(parts) * 60^(seq_along(parts) - 1L))
  stopifnot(length(rss) == 1L, is.finite(rss), is.finite(wall))
  c(rss_bytes = rss, wall_seconds = wall)
}

genome_write_csv <- function(x, path) {
  con <- gzfile(path, "wt", compression = 9L)
  on.exit(close(con))
  write.csv(x, con, row.names = FALSE)
}

genome_collect <- function(root, destination) {
  manifest_path <- file.path(root, "manifest.rds")
  manifest <- readRDS(manifest_path)
  key <- digest::digest(file = manifest_path, algo = "sha256")
  run_root <- file.path(root, key)
  files <- list.files(run_root, "\\.rds$", full.names = TRUE)
  strata <- genome_strata()
  fits <- performance <- checks <- list()
  for (i in seq_along(files)) {
    saved <- readRDS(files[[i]])
    stopifnot(identical(saved$manifest_sha256, key), saved$threads == genome_declaration$threads, saved$stride == 1L)
    stratum <- strata$stratum[vapply(seq_len(nrow(strata)), function(j) {
      identical(as.list(strata[j, c("depth", "share", "relationship", "error")]), saved$setting)
    }, logical(1L))]
    stopifnot(length(stratum) == 1L)
    stopifnot(
      saved$seed > genome_declaration$seed + stratum * 1000L,
      saved$seed <= genome_declaration$seed + stratum * 1000L + genome_declaration$replicates,
      all(saved$fit$seq_error == saved$setting$error),
      all(saved$fit$artefact_weight == genome_declaration$artefact_weight),
      all(saved$fit$homozygosity_excess == 0)
    )
    x <- saved$fit
    fields <- c(
      "relationship_model", "share", "share_lower", "share_upper", "share_se",
      "singles_share", "pair_residual", "status",
      "relation_diagnostic", "jackknife_failed_refits", "at_bound", "interval_hull",
      "reads", "positive_pair_sites", "blocks", "pair_blocks"
    )
    fits[[i]] <- cbind(data.frame(
      stratum = stratum, seed = saved$seed,
      sample = saved$identity$sample, donor = saved$identity$donor,
      sample_haplotype = saved$identity$sample_haplotype,
      donor_haplotype = saved$identity$donor_haplotype
    ), x[fields])
    ceiling <- genome_declaration$memory_overhead_bytes + genome_declaration$live_memory_multiplier *
      (genome_declaration$sql_live_bytes_per_site * saved$batch_input_rows + saved$live_bytes +
        genome_declaration$max_state_bytes)
    performance[[i]] <- data.frame(
      stratum = stratum, seed = saved$seed,
      sites = saved$qc[["sites"]], covered_sites = max(x$sites),
      zero_depth_sites = saved$qc[["zero_depth_sites"]], excluded_sites = saved$qc[["excluded_sites"]],
      unavailable_sites = saved$qc[["unavailable_sites"]],
      output_rows = nrow(x), simulation_seconds = saved$simulation_seconds,
      loading_seconds = saved$loading_seconds, fit_seconds = saved$fit_seconds,
      wall_seconds = saved$batch_wall_seconds / saved$batch_replicates,
      batch_first_seed = saved$batch_first_seed, batch_replicates = saved$batch_replicates,
      batch_input_rows = saved$batch_input_rows,
      oracle_included = !is.null(saved$oracle), rss_bytes = saved$process_rss_bytes,
      peak_buffer_bytes = max(saved$peak_buffer_bytes, saved$thread_check$peak_buffer_bytes),
      spill_bytes = max(saved$spill_bytes, saved$thread_check$spill_bytes),
      live_R_bytes = saved$live_bytes, memory_ceiling_bytes = ceiling
    )
    if (!is.null(saved$oracle)) {
      checks[[length(checks) + 1L]] <- cbind(
        data.frame(stratum = stratum, seed = saved$seed),
        as.data.frame(c(saved$oracle, saved$thread_check))
      )
    }
  }
  fits <- do.call(rbind, fits)
  performance <- do.call(rbind, performance)
  checks <- do.call(rbind, checks)
  stopifnot(
    !anyDuplicated(fits[c("stratum", "seed", "relationship_model")]),
    !anyDuplicated(performance[c("stratum", "seed")])
  )
  summaries <- lapply(seq_len(nrow(strata)), function(i) {
    setting <- strata[i, ]
    model <- if (setting$relationship == "ibd1") "ibd1" else "unrelated"
    x <- fits[fits$stratum == setting$stratum & fits$relationship_model == model, ]
    if (!nrow(x)) {
      return(NULL)
    }
    genome_summarize(x, setting)
  })
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  genome_write_csv(fits, file.path(destination, "replicates.csv.gz"))
  genome_write_csv(performance, file.path(destination, "performance.csv.gz"))
  write.csv(do.call(rbind, summaries), file.path(destination, "strata.csv"), row.names = FALSE)
  write.csv(checks, file.path(destination, "oracle_and_threads.csv"), row.names = FALSE)
  status_fields <- fits[c("stratum", "relationship_model", "status", "relation_diagnostic")]
  status_fields$status[is.na(status_fields$status)] <- "missing_status"
  status_fields$relation_diagnostic[is.na(status_fields$relation_diagnostic)] <- "missing_diagnostic"
  status <- aggregate(list(count = rep(1L, nrow(fits))), status_fields, sum)
  write.csv(status, file.path(destination, "statuses.csv"), row.names = FALSE)
  saveRDS(manifest, file.path(destination, "manifest.rds"))
  metadata <- list(
    declaration_commit = system2("git", c("rev-parse", "1735221a"), stdout = TRUE),
    batching_amendment_commit = system2("git", c("rev-parse", "957d03c9"), stdout = TRUE),
    execution_amendment_commit = system2("git", c("log", "-1", "--format=%H", "--",
      "benchmarks/count_mixture_lowpass_genome/execution_declaration.md"), stdout = TRUE),
    collection_revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    manifest_sha256 = key, sources = manifest,
    collection_driver_sha256 = digest::digest(file = "benchmarks/count_mixture_lowpass_genome/collect.R", algo = "sha256"),
    cpu = system2("lscpu", stdout = TRUE), date = format(Sys.time(), tz = "UTC")
  )
  jsonlite::write_json(metadata, file.path(destination, "metadata.json"), auto_unbox = TRUE, pretty = TRUE)
  invisible(do.call(rbind, summaries))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  genome_collect(args[[1L]], args[[2L]])
}
