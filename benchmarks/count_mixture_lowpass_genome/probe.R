source("benchmarks/count_mixture_lowpass_genome/batch.R")
args <- commandArgs(TRUE)
panels <- genome_load_panels()
strata <- genome_strata()
checks <- list()
measurements <- list()
for (depth in genome_declaration$depths) {
  setting <- strata[strata$depth == depth & strata$share == .1 &
    strata$relationship == "same" & strata$error == .0001, ]
  for (replicate in 1:3) {
    con <- genome_connect(args[[1L]], 7L)
    dbExecute(con, "SET max_temp_directory_size = '0B'")
    started <- proc.time()[["elapsed"]]
    full <- genome_counts(panels, setting, replicate, covered_only = FALSE)
    simulation <- proc.time()[["elapsed"]] - started
    all_fit <- genome_batch_fit(con, full$counts, setting$error, 7L, materialize = FALSE)
    depth_count <- full$counts$ref_count + full$counts$alt_count + full$counts$other_count
    covered <- full$counts[depth_count > 0, ]
    sparse <- genome_counts(panels, setting, replicate, covered_only = TRUE)
    rownames(covered) <- NULL
    stopifnot(identical(covered, sparse$counts), identical(full$qc, sparse$qc))
    rm(sparse)
    dbDisconnect(con, shutdown = TRUE)
    con <- genome_connect(args[[1L]], 1L)
    dbExecute(con, "SET max_temp_directory_size = '0B'")
    covered_fit <- genome_batch_fit(con, covered, setting$error, 1L)
    a <- genome_order_fit(all_fit$fit)[genome_exact_fields]
    b <- genome_order_fit(covered_fit$fit)[genome_exact_fields]
    exact <- identical(a, b)
    full_ceiling <- genome_declaration$memory_overhead_bytes + genome_declaration$live_memory_multiplier *
      (genome_declaration$sql_live_bytes_per_site * nrow(full$counts) +
        as.numeric(object.size(panels)) + as.numeric(object.size(full)) +
        genome_declaration$max_state_bytes)
    stopifnot(all_fit$peak_buffer_bytes <= full_ceiling, covered_fit$peak_buffer_bytes < 6e9)
    checks[[length(checks) + 1L]] <- data.frame(
      depth = depth, replicate = replicate, seed = full$seed, exact = exact,
      sites = full$qc[["sites"]], zero_depth_sites = full$qc[["zero_depth_sites"]],
      excluded_sites = full$qc[["excluded_sites"]], unavailable_sites = full$qc[["unavailable_sites"]],
      simulation_seconds = simulation, all_seconds = all_fit$seconds,
      covered_seconds = covered_fit$seconds,
      all_buffer_bytes = all_fit$peak_buffer_bytes, covered_buffer_bytes = covered_fit$peak_buffer_bytes
    )
    saveRDS(list(all = all_fit, covered = covered_fit, qc = full$qc),
      file.path(args[[2L]], sprintf("equivalence-depth%g-rep%d.rds", depth, replicate)))
    cat("Equivalence", depth, replicate, exact, "simulation", simulation,
        "all", all_fit$seconds, "covered", covered_fit$seconds, "\n")
    flush.console()
    rm(full, covered, all_fit, covered_fit)
    dbDisconnect(con, shutdown = TRUE)
    gc()
  }
}
write.csv(do.call(rbind, checks), file.path(args[[2L]], "equivalence.csv"), row.names = FALSE)
stopifnot(all(vapply(checks, function(x) x$exact, logical(1L))))
con <- genome_connect(args[[1L]], 1L)
dbExecute(con, "SET max_temp_directory_size = '0B'")
for (depth in c(.1, .2, .4)) {
  batch_size <- c(`0.1` = 4L, `0.2` = 2L, `0.4` = 1L)[[as.character(depth)]]
  setting <- strata[strata$depth == depth & strata$share == .1 &
    strata$relationship == "same" & strata$error == .0001, ]
  started <- proc.time()[["elapsed"]]
  simulated <- lapply(seq_len(batch_size), function(replicate) genome_counts(panels, setting, replicate))
  counts <- do.call(rbind, lapply(simulated, `[[`, "counts"))
  seconds <- proc.time()[["elapsed"]] - started
  fit <- tryCatch(genome_batch_fit(con, counts, setting$error, 1L), error = identity)
  if (inherits(fit, "error")) {
    cat("Batch", batch_size, "failed:", conditionMessage(fit), "\n")
    break
  }
  measurements[[length(measurements) + 1L]] <- data.frame(
    depth = depth, batch_size = batch_size, input_rows = nrow(counts), simulation_seconds = seconds,
    loading_seconds = fit$loading_seconds, fit_seconds = fit$seconds, buffer_bytes = fit$peak_buffer_bytes,
    spill_bytes = fit$spill_bytes, live_R_bytes = as.numeric(object.size(panels)) +
      as.numeric(object.size(simulated)) + as.numeric(object.size(counts)))
  cat("Batch", batch_size, "rows", nrow(counts), "simulation", seconds, "fit", fit$seconds,
      "buffer", fit$peak_buffer_bytes, "\n")
  flush.console()
  saveRDS(fit, file.path(args[[2L]], paste0("batch-", batch_size, ".rds")))
  rm(simulated, counts, fit)
  gc()
}
dbDisconnect(con, shutdown = TRUE)
write.csv(do.call(rbind, measurements), file.path(args[[2L]], "microbatch.csv"), row.names = FALSE)
