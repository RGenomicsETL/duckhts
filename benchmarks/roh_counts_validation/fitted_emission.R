# Decode the fitted read-count emission check and write every observation as
# CSV. benchmark_roh_fitted_emission.Rmd states the design; the values are in
# fitted_declaration.R and declaration.R.
#
# Every query runs in a fresh DuckDB CLI process under /usr/bin/time, with the
# JSON profiler and max_temp_directory_size = 0, and is gated on a memory
# ceiling declared in the report. Inputs come from the duckhtsbench registry
# and are checked against their registered identities first.
#
# Run from the repository root.

source("benchmarks/roh_counts_validation/declaration.R")
source("benchmarks/roh_counts_validation/fitted_declaration.R")
source("benchmarks/roh_counts_validation/metrics.R")
source("benchmarks/roh_counts_validation/inputs.R")

roh_fitted_quote <- function(value) paste0("'", gsub("'", "''", value, fixed = TRUE), "'")

# Run setup statements and one query in a fresh CLI process. The query's rows
# go to `output` (CSV with a header, or Parquet). Returns the latency, peak
# RSS, peak buffer and temporary bytes. A failed process stops the render.
roh_fitted_run_fresh <- function(cli, extension, setup, query, threads, output,
                                 format = c("csv", "parquet")) {
  format <- match.arg(format)
  scratch <- tempfile("roh-fitted-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  sql_path <- file.path(scratch, "query.sql")
  profile_path <- file.path(scratch, "profile.json")
  time_path <- file.path(scratch, "time.tsv")
  error_path <- file.path(scratch, "stderr.txt")
  copy_options <- if (format == "csv") "(HEADER)" else "(FORMAT parquet)"
  writeLines(c(
    sprintf("SET threads = %d;", threads),
    "SET max_temp_directory_size = '0KiB';",
    sprintf("LOAD %s;", roh_fitted_quote(extension)),
    setup,
    "PRAGMA enable_profiling = 'json';",
    sprintf("PRAGMA profiling_output = %s;", roh_fitted_quote(profile_path)),
    sprintf("COPY (%s) TO %s %s;", query, roh_fitted_quote(output), copy_options),
    "PRAGMA disable_profiling;"
  ), sql_path)
  status <- system2("/usr/bin/time", c("-f", shQuote("%M\t%e"), "-o", shQuote(time_path),
                                       shQuote(cli), "-unsigned", "-batch", "-bail"),
                    stdin = sql_path, stdout = FALSE, stderr = error_path)
  if (status != 0L) {
    stop("DuckDB CLI failed (status ", status, "): ",
         paste(readLines(error_path, warn = FALSE), collapse = " "), call. = FALSE)
  }
  profile <- jsonlite::fromJSON(profile_path, simplifyVector = TRUE)
  timing <- utils::read.delim(time_path, header = FALSE,
                              col.names = c("max_rss_kib", "process_seconds"))
  data.frame(seconds = as.numeric(profile$latency),
             max_rss_mib = timing$max_rss_kib[[1L]] / 1024,
             peak_buffer_mib = as.numeric(profile$system_peak_buffer_memory) / 2^20,
             peak_temp_bytes = as.numeric(profile$system_peak_temp_dir_size))
}

# The four numbers of one count decode for each fitted sample. `spec` names
# the source of each number ("default" or "fit"); `perturbation` is NULL or
# one row of roh_fitted_declaration$perturbations. Returns the values and
# whether a value was set to the end of its guard range.
roh_fitted_parameters <- function(fit, spec, perturbation = NULL,
                                  declaration = roh_fitted_declaration) {
  parameters <- c("seq_error", "contamination", "allele_balance", "artefact_weight")
  values <- data.frame(sample = fit$sample, stringsAsFactors = FALSE)
  for (parameter in parameters) {
    values[[parameter]] <- if (spec[[parameter]] == "fit") {
      as.numeric(fit[[parameter]])
    } else {
      rep(declaration$defaults[[parameter]], nrow(fit))
    }
  }
  values$at_guard <- FALSE
  if (!is.null(perturbation)) {
    parameter <- perturbation$parameter
    current <- values[[parameter]]
    changed <- switch(perturbation$operation,
      add = current + perturbation$value,
      scale = current * perturbation$value,
      add_sd = current + perturbation$value *
        ifelse(is.na(fit$contamination_sd), declaration$sd_fallback, fit$contamination_sd),
      stop("unknown perturbation operation: ", perturbation$operation, call. = FALSE))
    range <- declaration$guard_ranges[[parameter]]
    values$at_guard <- changed < range[[1L]] | changed > range[[2L]]
    values[[parameter]] <- pmin(pmax(changed, range[[1L]]), range[[2L]])
  }
  values
}

roh_fitted_runs_of <- function(runs, sample) {
  runs[runs$sample == sample, , drop = FALSE]
}

roh_fitted_run <- function(extension, cli, output_dir, samples,
                           validation = roh_validation_declaration,
                           declaration = roh_fitted_declaration) {
  stopifnot(length(samples) >= 1L, all(samples %in% validation$samples))
  threads <- declaration$threads
  extension <- normalizePath(extension, mustWork = TRUE)
  work <- tempfile("roh-fitted-work-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Inputs: every declared input is checked against its registry identity,
  # also in a reduced run, before anything is decoded.
  count_paths <- vapply(validation$count_artifacts[samples],
                        duckhtsbench::duckhts_bench_artifact_path, character(1L))
  genotype_path <- duckhtsbench::duckhts_bench_artifact_path(validation$genotype_artifact)
  con <- DBI::dbConnect(duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  roh_validation_check_inputs(con, validation)
  source_sql <- sprintf("read_parquet([%s])",
                        paste(vapply(count_paths, roh_fitted_quote, character(1L)), collapse = ", "))
  shape <- DBI::dbGetQuery(con, sprintf(paste(
    "SELECT count(*) AS rows, count(DISTINCT sample_id) AS samples,",
    "count(DISTINCT (chrom, pos, ref, alt, af)) AS sites, count(DISTINCT chrom) AS chroms,",
    "min(pos) AS first_pos, max(pos) AS last_pos FROM %s"), source_sql))
  if (shape$samples != length(samples) || shape$chroms != 1 ||
      shape$rows != shape$samples * shape$sites) {
    stop("the count artifacts do not share one chromosome's sites and frequencies", call. = FALSE)
  }
  sites <- as.numeric(shape$sites)
  span_bases <- as.numeric(shape$last_pos) - as.numeric(shape$first_pos) + 1
  # The fit's denominators, as benchmark_count_error_fit.Rmd counts them.
  cells <- DBI::dbGetQuery(con, sprintf(paste(
    "WITH r AS (SELECT sample_id, chrom, (pos - 1) // 10000000 AS block_index,",
    "floor(af * 64) AS bin, ref_count + alt_count AS depth, alt_count AS alt FROM %s",
    "WHERE ref_count IS NOT NULL AND alt_count IS NOT NULL AND af > 0 AND af < 1",
    "AND ref_count + alt_count BETWEEN 1 AND 1000)",
    "SELECT (SELECT count(*) FROM (SELECT DISTINCT sample_id, chrom, block_index, bin,",
    "depth, alt FROM r)) AS cells,",
    "(SELECT count(*) FROM (SELECT DISTINCT sample_id, depth, alt FROM r)) AS pairs"),
    source_sql))

  setup <- c(
    sprintf(paste("CREATE TEMP VIEW counts AS SELECT sample_id, chrom, pos, ref, alt,",
                  "ref_count, alt_count, af FROM %s;"), source_sql),
    sprintf("CREATE TEMP VIEW af_sites AS SELECT chrom, pos, ref, alt, af FROM read_parquet(%s);",
            roh_fitted_quote(count_paths[[1L]])))

  # Ceilings in MiB: overhead + 3 x the live state declared in the report.
  observations <- list()
  measure <- function(name, kind, query, output, format = "csv", ceiling_mib = NA_real_,
                      extra_setup = character(0)) {
    measured <- roh_fitted_run_fresh(cli, extension, c(setup, extra_setup), query, threads,
                                     output, format)
    observations[[length(observations) + 1L]] <<- data.frame(
      observation = name, kind = kind, measured, ceiling_mib = ceiling_mib)
    invisible(measured)
  }
  for (repetition in 1:3) {
    measure(sprintf("empty_%d", repetition), "empty", "SELECT 1 AS n",
            file.path(work, "empty.csv"))
  }
  empty <- do.call(rbind, observations)
  overhead_mib <- stats::median(empty$max_rss_mib)
  mib <- function(bytes) bytes / 2^20
  fit_ceiling <- overhead_mib + 3 * mib(cells$cells * 312 + cells$pairs * 88 +
                                         24 * 1001 * length(samples))
  genotype_ceiling <- overhead_mib + 3 * mib(length(samples) * sites * 16 + sites * 128 +
                                              threads * sites * 36.25)
  counts_ceiling <- function(decoded_samples) {
    overhead_mib + 3 * mib(decoded_samples * sites * 24 + threads * sites * 36.25)
  }

  # The fit, written to Parquet so that the decodes read its numbers exactly.
  fit_path <- file.path(work, "fit.parquet")
  measure("fit", "fit", "SELECT * FROM duckhts_count_error_fit('counts') ORDER BY sample",
          fit_path, "parquet", fit_ceiling)
  fit <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s) ORDER BY sample",
                                      roh_fitted_quote(fit_path)))
  if (!setequal(fit$sample, samples) || nrow(fit) != length(samples)) {
    stop("the fit does not have one row per decoded sample", call. = FALSE)
  }
  fitted <- fit[fit$status %in% declaration$fit_status, , drop = FALSE]
  decoded <- fitted$sample
  if (!length(decoded)) stop("no sample has a fit with status ok", call. = FALSE)
  samples_sql <- roh_fitted_quote(paste(decoded, collapse = ","))
  fitted_setup <- sprintf(
    "CREATE TEMP VIEW fitted_counts AS SELECT * FROM counts WHERE sample_id IN (%s);",
    paste(vapply(decoded, roh_fitted_quote, character(1L)), collapse = ", "))
  runs_columns <- "sample, chrom, start, \"end\", length, n_markers, quality"
  read_runs <- function(path) {
    utils::read.csv(path, stringsAsFactors = FALSE,
                    colClasses = c(sample = "character", chrom = "character"))
  }

  # Genotype decodes.
  genotype_runs <- list()
  for (gt_error in declaration$gt_errors) {
    output <- file.path(work, sprintf("gt_%s.csv", format(gt_error)))
    measure(sprintf("genotype_gt%s", format(gt_error)), "genotype_decode", sprintf(paste(
      "SELECT %s FROM duckhts_roh_af_table(%s, 'af_sites', gt_error := %s, samples := %s)",
      "ORDER BY sample, chrom, start"), runs_columns, roh_fitted_quote(genotype_path),
      format(gt_error), samples_sql), output, ceiling_mib = genotype_ceiling)
    genotype_runs[[format(gt_error)]] <- read_runs(output)
  }

  # Count decodes. Each fit relation is a Parquet file of exactly the numbers
  # in `parameters`.
  decode_counts <- function(name, parameters) {
    fit_file <- file.path(work, sprintf("fit_%s.parquet", name))
    DBI::dbWriteTable(con, "decode_fit", parameters[c("sample", "seq_error", "contamination",
                                                       "allele_balance", "artefact_weight")],
                      overwrite = TRUE)
    DBI::dbExecute(con, sprintf("COPY decode_fit TO %s (FORMAT parquet)",
                                roh_fitted_quote(fit_file)))
    output <- file.path(work, sprintf("counts_%s.csv", name))
    measure(name, "count_decode", sprintf(paste(
      "SELECT %s FROM duckhts_roh_counts_fitted('fitted_counts', 'decode_fit')",
      "ORDER BY sample, chrom, start"), runs_columns), output,
      ceiling_mib = counts_ceiling(length(decoded)), extra_setup = c(
        fitted_setup,
        sprintf("CREATE TEMP VIEW decode_fit AS SELECT * FROM read_parquet(%s);",
                roh_fitted_quote(fit_file))))
    read_runs(output)
  }
  decodes <- declaration$decodes
  decode_parameters <- list()
  count_runs <- list()
  for (index in seq_len(nrow(decodes))) {
    name <- decodes$decode[[index]]
    decode_parameters[[name]] <- data.frame(
      decode = name, roh_fitted_parameters(fitted, decodes[index, ], declaration = declaration))
    count_runs[[name]] <- decode_counts(paste0("decode_", name), decode_parameters[[name]])
  }
  # Decode a goes through the fit relation with the scalar defaults; it must
  # equal duckhts_roh_counts with its scalar defaults, run for run.
  scalar_output <- file.path(work, "counts_scalar.csv")
  measure("scalar_defaults", "count_decode", sprintf(paste(
    "SELECT %s FROM duckhts_roh_counts('fitted_counts') ORDER BY sample, chrom, start"),
    runs_columns), scalar_output, ceiling_mib = counts_ceiling(length(decoded)),
    extra_setup = fitted_setup)
  scalar_runs <- read_runs(scalar_output)
  defaults_identical <- identical(scalar_runs, count_runs[["a"]])

  # Agreement with the genotype decodes.
  comparisons <- list()
  for (gt_error in names(genotype_runs)) {
    for (name in names(count_runs)) {
      for (sample in decoded) {
        comparisons[[length(comparisons) + 1L]] <- data.frame(
          gt_error = as.numeric(gt_error), decode = name, sample = sample,
          roh_run_verdict(roh_run_comparison_by_class(
            roh_fitted_runs_of(count_runs[[name]], sample),
            roh_fitted_runs_of(genotype_runs[[gt_error]], sample),
            span_bases, validation$long_run_bases), validation))
      }
    }
  }
  comparisons <- do.call(rbind, comparisons)
  criterion_rows <- comparisons[comparisons$gt_error == declaration$criterion_gt_error &
                                comparisons$run_class == "all", ]
  base_gap <- function(name) {
    rows <- criterion_rows[criterion_rows$decode == name, ]
    stats::setNames(abs(rows$test_bases - rows$reference_bases), rows$sample)[decoded]
  }
  criterion <- data.frame(sample = decoded, gap_a = base_gap("a"), gap_b = base_gap("b"),
                          gap_c = base_gap("c"), row.names = NULL)
  criterion$pass <- criterion$gap_c <= pmin(criterion$gap_a, criterion$gap_b)

  # Sensitivity of decode c to each number on its own.
  full <- decodes[decodes$decode == "c", ]
  perturbations <- declaration$perturbations
  perturbations$perturbation <- sprintf("%s_%s_%s", perturbations$parameter,
                                        perturbations$operation, format(perturbations$value))
  sensitivity <- list()
  perturbed_parameters <- list()
  perturbed_runs <- list()
  for (index in seq_len(nrow(perturbations))) {
    perturbation <- perturbations[index, ]
    parameters <- roh_fitted_parameters(fitted, full, perturbation, declaration)
    perturbed_parameters[[index]] <- data.frame(decode = perturbation$perturbation, parameters)
    runs <- decode_counts(paste0("perturbed_", index), parameters)
    perturbed_runs[[index]] <- if (nrow(runs)) data.frame(decode = perturbation$perturbation, runs)
    for (sample in decoded) {
      metrics <- roh_run_comparison_by_class(roh_fitted_runs_of(runs, sample),
                                             roh_fitted_runs_of(count_runs[["c"]], sample),
                                             span_bases, validation$long_run_bases)
      metrics <- metrics[metrics$run_class == "all", ]
      changed_bases <- metrics$either_bases - metrics$shared_bases
      changed_share <- if (metrics$reference_bases > 0) {
        changed_bases / metrics$reference_bases
      } else if (changed_bases == 0) 0 else Inf
      sample_parameters <- parameters[parameters$sample == sample, ]
      sensitivity[[length(sensitivity) + 1L]] <- data.frame(
        perturbation[c("perturbation", "parameter", "operation", "value")], sample = sample,
        fitted_value = as.numeric(fitted[fitted$sample == sample, perturbation$parameter]),
        perturbed_value = sample_parameters[[perturbation$parameter]],
        at_guard = sample_parameters$at_guard,
        metrics[c("reference_runs", "test_runs", "reference_bases", "test_bases",
                  "shared_bases", "either_bases", "froh_difference", "jaccard")],
        bases_change = metrics$test_bases - metrics$reference_bases,
        changed_bases = changed_bases, changed_share = changed_share, row.names = NULL)
    }
  }
  sensitivity <- do.call(rbind, sensitivity)
  sensitivity$froh_pass <- abs(sensitivity$froh_difference) < declaration$max_froh_change
  sensitivity$share_pass <- sensitivity$changed_share < declaration$max_changed_share
  sensitivity$pass <- sensitivity$froh_pass & sensitivity$share_pass

  observations <- do.call(rbind, observations)
  observations$within_ceiling <- is.na(observations$ceiling_mib) |
    observations$max_rss_mib <= observations$ceiling_mib
  observations$no_spill <- observations$peak_temp_bytes == 0

  revision <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
  src_tree <- system2("git", c("rev-parse", "HEAD:src"), stdout = TRUE)
  metadata <- data.frame(
    key = c("revision", "src_tree", "input_identities", "extension_sha256", "duckdb_cli",
            "r_version", "host_load", "cpu_affinity", "threads", "samples", "decoded_samples", "sites",
            "span_bases", "fit_cells", "fit_pairs", "overhead_mib",
            "defaults_identical_to_scalar"),
    value = c(revision, src_tree, "match the registry",
              strsplit(system2("sha256sum", shQuote(extension), stdout = TRUE), " ",
                       fixed = TRUE)[[1L]][[1L]],
              system2(cli, "--version", stdout = TRUE)[[1L]], R.version.string,
              sub(".*load average: ", "", system2("uptime", stdout = TRUE)),
              sub("^.*: *", "", system2("taskset", c("-pc", Sys.getpid()), stdout = TRUE)),
              threads,
              paste(samples, collapse = ","), paste(decoded, collapse = ","), format(sites),
              format(span_bases, scientific = FALSE), format(cells$cells), format(cells$pairs),
              sprintf("%.1f", overhead_mib), defaults_identical),
    stringsAsFactors = FALSE)

  all_parameters <- rbind(do.call(rbind, decode_parameters), do.call(rbind, perturbed_parameters))
  all_runs <- rbind(
    do.call(rbind, lapply(names(genotype_runs), function(gt_error) {
      if (nrow(genotype_runs[[gt_error]])) {
        data.frame(decode = paste0("genotype_gt", gt_error), genotype_runs[[gt_error]])
      }
    })),
    do.call(rbind, lapply(names(count_runs), function(name) {
      if (nrow(count_runs[[name]])) data.frame(decode = name, count_runs[[name]])
    })),
    do.call(rbind, perturbed_runs))
  write <- function(data, name) {
    utils::write.csv(data, file.path(output_dir, name), row.names = FALSE)
  }
  write(metadata, "run_metadata.csv")
  write(fit, "fit.csv")
  write(all_parameters, "decode_parameters.csv")
  write(all_runs, "runs.csv")
  write(comparisons, "agreement_metrics.csv")
  write(criterion, "criterion.csv")
  write(sensitivity, "sensitivity.csv")
  write(observations, "observations.csv")
  list(metadata = metadata, fit = fit, parameters = all_parameters, comparisons = comparisons,
       criterion = criterion, sensitivity = sensitivity, observations = observations,
       defaults_identical = defaults_identical, span_bases = span_bases)
}
