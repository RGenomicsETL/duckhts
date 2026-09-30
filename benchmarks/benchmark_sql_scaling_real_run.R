# Run from the repository root. Each timed query uses a fresh DuckDB process.
args <- commandArgs(trailingOnly = TRUE)
source("r/duckhtsbench/R/registry.R")
source("r/duckhtsbench/R/stage.R")
Sys.setenv(DUCKHTSBENCH_REGISTRY = normalizePath(
  "r/duckhtsbench/inst/benchmark_registry.tsv", mustWork = TRUE))

vcf_counts_budget <- list(
  memory_limit = "256MB",
  max_temp_directory_size = "0B",
  rss_mib = 256,
  peak_buffer_mib = 256,
  spill_mib = 0,
  file_doubling_delta_mib = 16
)

sha256_file <- function(path) {
  result <- system2("sha256sum", shQuote(path), stdout = TRUE)
  if (length(result) != 1L || !is.null(attr(result, "status"))) {
    stop("could not calculate SHA-256 for ", path)
  }
  sub("[[:space:]].*$", "", result[[1L]])
}

git_diff_sha256 <- function() {
  diff_file <- tempfile("duckhts-scaling-diff-")
  on.exit(unlink(diff_file))
  status <- system2("git", c("diff", "--binary", "HEAD"), stdout = diff_file)
  if (status != 0L) stop("could not identify the working-tree diff")
  sha256_file(diff_file)
}

package_version_or_na <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(package))
}

count_vcf_records <- function(path) {
  if (grepl("\\.vcf$", path, ignore.case = TRUE)) {
    count <- system2("grep", c("-vc", shQuote("^#"), shQuote(path)), stdout = TRUE)
  } else {
    count <- system2("bcftools", c("index", "-n", shQuote(path)), stdout = TRUE)
  }
  if (length(count) != 1L || !grepl("^[0-9]+$", count[[1L]])) {
    stop("could not count VCF records: ", path)
  }
  as.integer(count[[1L]])
}

profile_value <- function(profile, field) {
  if (is.null(profile) || is.null(profile[[field]])) return(NA_real_)
  as.numeric(profile[[field]]) / 1048576
}

private_temp_file_bytes <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  if (length(files) == 0L) return(0)
  info <- file.info(files)
  sum(info$size[!info$isdir], na.rm = TRUE)
}

prepare_vcf_count_panels <- function(extension, source_path, duplicate_path) {
  panel_dir <- tempfile("duckhts-vcf-count-panels-")
  dir.create(panel_dir)
  connection <- DBI::dbConnect(
    duckdb::duckdb(config = list(allow_unsigned_extensions = "true"),
                   shared_home = FALSE))
  on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
  quote_string <- function(value) as.character(DBI::dbQuoteString(connection, value))
  DBI::dbExecute(connection, paste("LOAD", quote_string(extension)))
  query <- paste0(
    "WITH candidates AS (SELECT CHROM AS region, POS AS position, ",
    "REF AS ref, ALT[1] AS alt FROM read_bcf(", quote_string(source_path),
    ", samples := '', scan_mode := 'sequential') ",
    "WHERE len(REF) = 1 AND len(ALT) = 1 AND len(ALT[1]) = 1 ",
    "AND REF IN ('A','C','G','T') AND ALT[1] IN ('A','C','G','T') ",
    "AND REF != ALT[1] QUALIFY row_number() OVER (PARTITION BY CHROM, POS ",
    "ORDER BY REF, ALT[1]) = 1) ",
    "SELECT region, position, ref, alt FROM candidates ",
    "ORDER BY region, position, ref, alt LIMIT 8000")
  sites <- DBI::dbGetQuery(connection, query)
  if (nrow(sites) != 8000L) stop("the 1x VCF yielded fewer than 8,000 panel sites")
  panels <- list()
  for (size in c(2000L, 4000L, 8000L)) {
    selected <- sites[seq_len(size), , drop = FALSE]
    panel <- data.frame(
      assembly = "GRCh38",
      site_index = seq_len(size) - 1L,
      region = as.character(selected$region),
      position = as.integer(selected$position),
      allele_a = pmin(as.character(selected$ref), as.character(selected$alt)),
      allele_b = pmax(as.character(selected$ref), as.character(selected$alt)),
      stringsAsFactors = FALSE
    )
    path <- file.path(panel_dir, sprintf("panel-%dx.tsv", size / 2000L))
    utils::write.table(panel, path, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = TRUE)
    panels[[as.character(size)]] <- list(
      path = path, sites = size, sha256 = sha256_file(path))
  }

  seed <- system2("grep", c("-m", "1", "-v", shQuote("^#"),
                            shQuote(duplicate_path)), stdout = TRUE)
  if (length(seed) != 1L || !is.null(attr(seed, "status"))) {
    stop("could not read the duplicate-error VCF seed record")
  }
  fields <- strsplit(seed[[1L]], "\t", fixed = TRUE)[[1L]]
  if (length(fields) < 5L) stop("duplicate-error VCF seed record is malformed")
  duplicate_panel <- data.frame(
    assembly = "GRCh38", site_index = 0L, region = fields[[1L]],
    position = as.integer(fields[[2L]]),
    allele_a = sort(fields[c(4L, 5L)])[[1L]],
    allele_b = sort(fields[c(4L, 5L)])[[2L]],
    stringsAsFactors = FALSE
  )
  duplicate_panel_path <- file.path(panel_dir, "panel-duplicate.tsv")
  utils::write.table(duplicate_panel, duplicate_panel_path, sep = "\t",
                     quote = FALSE, row.names = FALSE, col.names = TRUE)
  panels[["duplicate"]] <- list(
    path = duplicate_panel_path, sites = 1L,
    sha256 = sha256_file(duplicate_panel_path))
  list(directory = panel_dir, panels = panels)
}

baseline_git_sha <- function() {
  sha <- Sys.getenv("DUCKHTS_BASELINE_SHA")
  if (!grepl("^[0-9a-f]{40}$", sha)) {
    stop("DUCKHTS_BASELINE_SHA must name the full commit the baseline extension was built from")
  }
  sha
}

matrix_build_identity <- function(implementation, branch_extension,
                                  baseline_extension, run_dir) {
  if (implementation == "prechange") {
    extension <- file.path(run_dir, "duckhts.duckdb_extension")
    if (!file.copy(baseline_extension, extension)) {
      stop("could not stage frozen extension")
    }
    return(list(
      extension = extension,
      source_git_sha = baseline_git_sha(),
      source_dirty_diff_sha256 = paste0(
        "e3b0c44298fc1c149afbf4c8996fb924",
        "27ae41e4649b934ca495991b7852b855")
    ))
  }
  list(
    extension = normalizePath(branch_extension, mustWork = TRUE),
    source_git_sha = system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[[1L]],
    source_dirty_diff_sha256 = git_diff_sha256()
  )
}

vcf_counts_resource_verdicts <- function(peak_rss_mib, peak_buffer_mib,
                                         peak_spill_mib) {
  rss <- if (peak_rss_mib <= vcf_counts_budget$rss_mib) {
    "within_budget"
  } else {
    "over_budget"
  }
  buffer <- if (is.na(peak_buffer_mib)) {
    "unavailable"
  } else if (peak_buffer_mib <= vcf_counts_budget$peak_buffer_mib) {
    "within_budget"
  } else {
    "over_budget"
  }
  spill <- if (is.na(peak_spill_mib)) {
    "unavailable"
  } else if (peak_spill_mib == vcf_counts_budget$spill_mib) {
    "within_budget"
  } else {
    "over_budget"
  }
  overall <- if (rss == "within_budget" && buffer == "within_budget" &&
                 spill == "within_budget") {
    "within_budget"
  } else {
    "over_or_unmeasured"
  }
  list(rss = rss, buffer = buffer, spill = spill, overall = overall)
}

hash_vcf_counts_result <- function(result) {
  if (is.null(result)) return(NA_character_)
  comparison_keys <- c("sample_id", "site_index")
  if (!all(comparison_keys %in% names(result))) {
    stop("vcf_counts output lacks its sample/site comparison key")
  }
  result <- result[do.call(order, result[comparison_keys]), , drop = FALSE]
  row.names(result) <- NULL
  result_file <- tempfile("duckhts-vcf-count-result-")
  on.exit(unlink(result_file))
  saveRDS(result, result_file, version = 3, compress = FALSE)
  sha256_file(result_file)
}

measure_vcf_counts_query <- function(connection, query, profile_file,
                                     private_temp) {
  DBI::dbExecute(connection, "PRAGMA enable_profiling='json'")
  quote_string <- function(value) as.character(DBI::dbQuoteString(connection, value))
  DBI::dbExecute(connection, paste("PRAGMA profiling_output =",
                                   quote_string(profile_file)))
  result <- NULL
  error_message <- NA_character_
  started <- proc.time()[["elapsed"]]
  tryCatch({
    result <- DBI::dbGetQuery(connection, query)
  }, error = function(error) {
    error_message <<- gsub("[[:space:]]+", " ", conditionMessage(error))
  })
  seconds <- proc.time()[["elapsed"]] - started
  execution_status <- if (is.null(result)) {
    if (grepl("out of memory|memory limit", error_message, ignore.case = TRUE)) {
      "out_of_memory"
    } else if (grepl("multiple source records occur at a panel site",
                     error_message, fixed = TRUE)) {
      "duplicate_error"
    } else {
      "query_error"
    }
  } else {
    "success"
  }
  profile <- if (file.exists(profile_file)) jsonlite::fromJSON(profile_file) else NULL
  peak_buffer_mib <- profile_value(profile, "system_peak_buffer_memory")
  if (is.null(profile)) {
    peak_spill_mib <- private_temp_file_bytes(private_temp) / 1048576
    spill_source <- "private_temp_directory_after_query"
  } else {
    peak_spill_mib <- profile_value(profile, "system_peak_temp_dir_size")
    spill_source <- "duckdb_query_json_profile"
  }
  rss <- readLines("/proc/self/status", warn = FALSE)
  rss_line <- grep("^VmHWM:", rss, value = TRUE)
  peak_rss_mib <- as.numeric(gsub("[^0-9]", "", rss_line[[1L]])) / 1024
  verdicts <- vcf_counts_resource_verdicts(
    peak_rss_mib, peak_buffer_mib, peak_spill_mib)
  list(
    seconds = seconds, execution_status = execution_status,
    error_message = error_message,
    output_rows = if (is.null(result)) 0L else nrow(result),
    result_sha256 = hash_vcf_counts_result(result),
    peak_rss_mib = peak_rss_mib, peak_buffer_mib = peak_buffer_mib,
    peak_spill_mib = peak_spill_mib, spill_source = spill_source,
    rss_verdict = verdicts$rss, buffer_verdict = verdicts$buffer,
    spill_verdict = verdicts$spill, resource_verdict = verdicts$overall,
    profile_available = !is.null(profile)
  )
}

run_vcf_count_matrix_child <- function(plan_path, plan_row, branch_extension,
                                       baseline_extension) {
  row <- utils::read.delim(plan_path, stringsAsFactors = FALSE,
                           check.names = FALSE)[plan_row, , drop = FALSE]
  implementation <- row$implementation[[1L]]
  run_dir <- tempfile("duckhts-vcf-count-run-")
  dir.create(run_dir)
  private_temp <- file.path(run_dir, "duckdb-temp")
  dir.create(private_temp)
  profile_file <- tempfile("duckhts-vcf-count-profile-", fileext = ".json")
  build <- matrix_build_identity(implementation, branch_extension,
                                 baseline_extension, run_dir)
  extension <- build$extension
  extension_sha256 <- sha256_file(extension)
  connection <- DBI::dbConnect(
    duckdb::duckdb(config = list(allow_unsigned_extensions = "true"),
                   shared_home = FALSE))
  on.exit({
    DBI::dbDisconnect(connection, shutdown = TRUE)
    unlink(c(run_dir, profile_file), recursive = TRUE)
  })
  quote_string <- function(value) as.character(DBI::dbQuoteString(connection, value))
  DBI::dbExecute(connection, paste("LOAD", quote_string(extension)))
  DBI::dbExecute(connection, paste0(
    "SET memory_limit = ", quote_string(vcf_counts_budget$memory_limit)))
  DBI::dbExecute(connection, paste0(
    "SET max_temp_directory_size = ",
    quote_string(vcf_counts_budget$max_temp_directory_size)))
  DBI::dbExecute(connection, paste("SET temp_directory =", quote_string(private_temp)))
  DBI::dbExecute(connection, sprintf("SET threads = %d", row$threads[[1L]]))

  macro_text <- DBI::dbGetQuery(connection, paste0(
    "SELECT macro_definition FROM duckdb_functions() ",
    "WHERE function_name = 'duckhts_somalier_vcf_counts'"))$macro_definition
  if (length(macro_text) != 1L || is.na(macro_text[[1L]])) {
    stop("loaded extension did not expose the vcf_counts macro definition")
  }
  macro_file <- tempfile("duckhts-vcf-count-macro-")
  writeBin(charToRaw(macro_text[[1L]]), macro_file)
  macro_text_sha256 <- sha256_file(macro_file)
  unlink(macro_file)
  version <- DBI::dbGetQuery(connection, "PRAGMA version")
  htslib_version <- DBI::dbGetQuery(
    connection, "SELECT duckhts_htslib_version()")[[1L]][[1L]]
  DBI::dbExecute(connection, paste0(
    "CREATE TEMP TABLE panel AS SELECT assembly::VARCHAR AS assembly, ",
    "site_index::UBIGINT AS site_index, region::VARCHAR AS region, ",
    "position::UBIGINT AS position, allele_a::VARCHAR AS allele_a, ",
    "allele_b::VARCHAR AS allele_b FROM read_csv(",
    quote_string(row$panel_path[[1L]]),
    ", delim := '\\t', header := true)"))

  input_path <- row$input_path[[1L]]
  input_records <- count_vcf_records(input_path)
  input_sha256 <- sha256_file(input_path)
  panel_source_sha256 <- sha256_file(row$panel_source_path[[1L]])
  input_bytes <- as.numeric(file.info(input_path)$size)
  panel_source_bytes <- as.numeric(file.info(row$panel_source_path[[1L]])$size)
  started_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  query <- paste0("SELECT * FROM duckhts_somalier_vcf_counts(",
                  quote_string(input_path), ", 'panel')")
  measurement <- measure_vcf_counts_query(
    connection, query, profile_file, private_temp)

  output <- data.frame(
    run_id = row$run_id[[1L]], implementation = implementation,
    cell = row$cell[[1L]], axis = row$axis[[1L]],
    file_scale = row$file_scale[[1L]], panel_scale = row$panel_scale[[1L]],
    input_id = row$input_id[[1L]], input_sha256 = input_sha256,
    input_bytes = input_bytes, input_records = input_records,
    panel_source_id = row$panel_source_id[[1L]],
    panel_source_sha256 = panel_source_sha256,
    panel_source_bytes = panel_source_bytes,
    panel_sites = row$panel_sites[[1L]], panel_sha256 = row$panel_sha256[[1L]],
    threads = row$threads[[1L]], replicate = row$replicate[[1L]],
    started_utc = started_utc, seconds = measurement$seconds,
    execution_status = measurement$execution_status,
    error_message = measurement$error_message,
    output_rows = measurement$output_rows,
    result_sha256 = measurement$result_sha256,
    peak_rss_mib = measurement$peak_rss_mib,
    peak_buffer_mib = measurement$peak_buffer_mib,
    peak_spill_mib = measurement$peak_spill_mib,
    spill_source = measurement$spill_source,
    rss_verdict = measurement$rss_verdict,
    buffer_verdict = measurement$buffer_verdict,
    spill_verdict = measurement$spill_verdict,
    resource_verdict = measurement$resource_verdict,
    memory_limit = vcf_counts_budget$memory_limit,
    max_temp_directory_size = vcf_counts_budget$max_temp_directory_size,
    private_temp_directory = private_temp,
    source_git_sha = build$source_git_sha,
    source_dirty_diff_sha256 = build$source_dirty_diff_sha256,
    extension_sha256 = extension_sha256,
    macro_text_sha256 = macro_text_sha256,
    duckdb_version = version$library_version[[1L]],
    duckdb_source_id = version$source_id[[1L]],
    R_version = R.version.string,
    DBI_version = package_version_or_na("DBI"),
    duckdb_R_version = package_version_or_na("duckdb"),
    Rduckhts_installed_version = package_version_or_na("Rduckhts"),
    Rduckhts_source_version = read.dcf("r/Rduckhts/DESCRIPTION", "Version")[[1L]],
    duckhtsbench_source_version = read.dcf("r/duckhtsbench/DESCRIPTION", "Version")[[1L]],
    htslib_version = htslib_version,
    profile_available = measurement$profile_available,
    stringsAsFactors = FALSE
  )
  utils::write.table(output, row$result_path[[1L]], sep = "\t",
                     row.names = FALSE, col.names = TRUE, quote = TRUE,
                     na = "NA")
  invisible(output)
}

validate_vcf_count_matrix_inputs <- function(baseline_extension) {
  baseline_sha256 <- sha256_file(baseline_extension)
  if (!identical(baseline_sha256,
                 "9a23e378a4dd22b83526d05aef14c239ba144be604da75cbc44f07d3bc893f58")) {
    stop("baseline extension SHA-256 does not match the frozen lab artifact")
  }
  source_id <- "sql_scaling_giab_1x"
  duplicate_id <- "sql_scaling_vcf_counts_duplicate_million"
  input_ids <- c(source_id, "sql_scaling_giab_2x", "sql_scaling_giab_4x",
                 duplicate_id)
  input_paths <- setNames(vapply(input_ids, duckhts_bench_artifact_path,
                                 character(1)), input_ids)
  for (id in input_ids) {
    if (!file.exists(input_paths[[id]])) stop("staged input missing: ", id)
    duckhts_bench_validate_identity(id)
    needs_index <- grepl("\\.vcf\\.gz$", input_paths[[id]])
    index_exists <- file.exists(paste0(input_paths[[id]], ".tbi"))
    if (needs_index && !index_exists) stop("tabix index missing for ", id)
  }
  list(source_id = source_id, duplicate_id = duplicate_id,
       source_path = input_paths[[source_id]],
       duplicate_path = input_paths[[duplicate_id]], input_paths = input_paths)
}

vcf_count_matrix_cells <- function(source_id, duplicate_id) {
  data.frame(
    cell = c("panel_1x_fixed_file_4x", "panel_2x_fixed_file_4x",
             "panel_4x_fixed_file_4x", "file_1x_fixed_panel_1x",
             "file_2x_fixed_panel_1x", "file_4x_fixed_panel_1x",
             "joint_file_2x_panel_2x", "duplicate_error_1m"),
    axis = c("panel", "panel", "panel", "file", "file", "file",
             "joint", "duplicate_error"),
    input_id = c("sql_scaling_giab_4x", "sql_scaling_giab_4x",
                 "sql_scaling_giab_4x", source_id,
                 "sql_scaling_giab_2x", "sql_scaling_giab_4x",
                 "sql_scaling_giab_2x", duplicate_id),
    panel_key = c("2000", "4000", "8000", "2000", "2000", "2000",
                  "4000", "duplicate"),
    file_scale = c(4L, 4L, 4L, 1L, 2L, 4L, 2L, NA_integer_),
    panel_scale = c(1L, 2L, 4L, 1L, 1L, 1L, 2L, NA_integer_),
    stringsAsFactors = FALSE
  )
}

build_vcf_count_matrix_plan <- function(cells, inputs, source_id, source_path,
                                        panels) {
  schedule <- expand.grid(
    implementation_order = 1:2, replicate = 1:3, threads = c(1L, 4L),
    cell_index = seq_len(nrow(cells)), KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE)
  schedule$implementation <- ifelse(
    schedule$replicate %% 2L == 1L,
    ifelse(schedule$implementation_order == 1L, "prechange", "branch"),
    ifelse(schedule$implementation_order == 1L, "branch", "prechange"))
  selected <- cells[schedule$cell_index, , drop = FALSE]
  selected_panels <- lapply(selected$panel_key, function(key) panels[[key]])
  plan <- data.frame(
    implementation = schedule$implementation,
    cell = selected$cell, axis = selected$axis,
    file_scale = selected$file_scale, panel_scale = selected$panel_scale,
    input_id = selected$input_id,
    input_path = unname(inputs[selected$input_id]),
    panel_source_id = source_id,
    panel_source_path = source_path,
    panel_path = vapply(selected_panels, `[[`, character(1), "path"),
    panel_sites = vapply(selected_panels, `[[`, integer(1), "sites"),
    panel_sha256 = vapply(selected_panels, `[[`, character(1), "sha256"),
    threads = schedule$threads, replicate = schedule$replicate,
    stringsAsFactors = FALSE
  )
  plan$run_id <- sprintf("%s-%s-t%d-r%d", plan$implementation, plan$cell,
                         plan$threads, plan$replicate)
  plan[c("run_id", setdiff(names(plan), "run_id"))]
}

run_vcf_count_matrix_processes <- function(plan, plan_path, branch_extension,
                                           baseline_extension, output_path) {
  if (file.exists(output_path)) unlink(output_path)
  for (plan_row in seq_len(nrow(plan))) {
    stderr_file <- tempfile("duckhts-vcf-count-stderr-")
    child_args <- c("benchmarks/benchmark_sql_scaling_real_run.R",
                    "--vcf-counts-matrix-child", plan_path, plan_row,
                    branch_extension, baseline_extension)
    child_output <- suppressWarnings(system2(
      "Rscript", shQuote(child_args), stdout = TRUE, stderr = stderr_file))
    child_status <- attr(child_output, "status")
    child_failed <- !is.null(child_status) && child_status != 0L
    result_path <- plan$result_path[[plan_row]]
    if (child_failed || !file.exists(result_path)) {
      status_text <- if (is.null(child_status)) "zero" else as.character(child_status)
      details <- c(paste0("child exit status: ", status_text),
                   paste0("stdout: ", paste(child_output, collapse = " | ")),
                   paste0("stderr: ", paste(readLines(stderr_file, warn = FALSE),
                                            collapse = " | ")))
      unlink(stderr_file)
      stop("matrix child failed for ", plan$run_id[[plan_row]], ": ",
           paste(details, collapse = "\n"))
    }
    result <- utils::read.delim(result_path, stringsAsFactors = FALSE,
                                check.names = FALSE, na.strings = "NA")
    if (nrow(result) != 1L) stop("matrix child returned an invalid row")
    unlink(result_path)
    append_output <- file.exists(output_path)
    utils::write.table(result, output_path, sep = "\t", quote = TRUE,
                       row.names = FALSE, col.names = !append_output,
                       append = append_output, na = "NA")
    message(result$run_id[[1L]], ": ", result$execution_status[[1L]],
            "; RSS ", sprintf("%.1f", result$peak_rss_mib[[1L]]), " MiB; ",
            "buffer ", ifelse(is.na(result$peak_buffer_mib[[1L]]), "NA",
                              sprintf("%.1f", result$peak_buffer_mib[[1L]])),
            " MiB; spill ", sprintf("%.1f", result$peak_spill_mib[[1L]]), " MiB")
    unlink(stderr_file)
  }
  invisible(output_path)
}

run_vcf_count_matrix <- function(branch_extension, baseline_extension) {
  branch_extension <- normalizePath(branch_extension, mustWork = TRUE)
  baseline_extension <- normalizePath(baseline_extension, mustWork = TRUE)
  inputs <- validate_vcf_count_matrix_inputs(baseline_extension)
  panels <- prepare_vcf_count_panels(
    branch_extension, inputs$source_path, inputs$duplicate_path)
  on.exit(unlink(panels$directory, recursive = TRUE), add = TRUE)
  cells <- vcf_count_matrix_cells(inputs$source_id, inputs$duplicate_id)
  plan <- build_vcf_count_matrix_plan(
    cells, inputs$input_paths, inputs$source_id, inputs$source_path,
    panels$panels)
  results_dir <- tempfile("duckhts-vcf-count-results-")
  dir.create(results_dir)
  on.exit(unlink(results_dir, recursive = TRUE), add = TRUE)
  plan$result_path <- file.path(results_dir, paste0(plan$run_id, ".tsv"))
  plan_path <- tempfile("duckhts-vcf-count-plan-")
  on.exit(unlink(plan_path), add = TRUE)
  utils::write.table(plan, plan_path, sep = "\t", quote = TRUE,
                     row.names = FALSE, col.names = TRUE)
  output_path <- "benchmarks/sql_scaling_vcf_counts_runs.tsv"
  run_vcf_count_matrix_processes(
    plan, plan_path, branch_extension, baseline_extension, output_path)
  invisible(output_path)
}

if (length(args) == 5L && args[[1L]] == "--vcf-counts-matrix-child") {
  run_vcf_count_matrix_child(args[[2L]], as.integer(args[[3L]]),
                             args[[4L]], args[[5L]])
} else if (length(args) == 3L && args[[2L]] == "--vcf-counts-matrix") {
  run_vcf_count_matrix(args[[1L]], args[[3L]])
} else if (length(args) == 6L && args[[1L]] == "--child") {
  library(DBI)
  case <- args[[2L]]
  threads <- as.integer(args[[3L]])
  multiplier <- as.integer(args[[4L]])
  extension <- args[[5L]]
  input <- args[[6L]]
  con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true"),
                                  shared_home = FALSE))
  quote_string <- function(value) as.character(dbQuoteString(con, value))
  dbExecute(con, paste("LOAD", quote_string(extension)))
  dbExecute(con, sprintf("SET threads=%d", threads))
  if (case == "vcf_counts") {
    dbExecute(con, paste0(
      "CREATE TEMP TABLE panel AS WITH sites AS (",
      "SELECT CHROM AS region, POS AS position, REF AS ref, ALT[1] AS alt ",
      "FROM read_bcf(", quote_string(input), ", samples := '', ",
      "scan_mode := 'sequential') WHERE len(REF) = 1 AND len(ALT) = 1 ",
      "AND len(ALT[1]) = 1 AND REF IN ('A','C','G','T') ",
      "AND ALT[1] IN ('A','C','G','T') AND REF != ALT[1] ",
      "ORDER BY POS LIMIT 2000) ",
      "SELECT 'GRCh38' AS assembly, ",
      "CAST(row_number() OVER (ORDER BY position) - 1 AS UBIGINT) AS site_index, ",
      "region, CAST(position AS UBIGINT) AS position, ",
      "least(ref,alt) AS allele_a, greatest(ref,alt) AS allele_b FROM sites"))
    query <- paste0("CREATE TEMP TABLE result AS SELECT * FROM ",
                    "duckhts_somalier_vcf_counts(", quote_string(input), ", 'panel')")
  } else if (case == "import_sites") {
    query <- paste0("CREATE TEMP TABLE result AS SELECT * FROM ",
                    "duckhts_somalier_import_sites(", quote_string(input),
                    ", 'GRCh37', max_sites := 2000000)")
  } else if (case %in% c("munge", "munge_metal", "r_munge", "r_munge_metal")) {
    limit <- 200000L * multiplier
    dbExecute(con, paste0(
      "CREATE TEMP VIEW source AS SELECT CAST(CHR AS VARCHAR) AS CHR, ",
      "* EXCLUDE CHR FROM read_csv(", quote_string(input),
      ", delim := '\t', header := true) LIMIT ", limit))
    column_map <- paste0(
      "map(['CHR','BP','SNP','A1','A2','P','Z','FRQ','NEFF','HET_I2',",
      "'HET_P','DIRE'], ['CHR','BP','MarkerName','Allele1','Allele2',",
      "'P-value','Zscore','Freq1','Weight','HetISq','HetPVal','Direction'])")
    macro <- if (case %in% c("munge", "r_munge")) "duckdb_munge" else "duckdb_munge_metal"
    fasta <- duckhts_bench_artifact_path("liftover_grch37_fasta")
    stopifnot(file.exists(fasta))
    query <- paste0("CREATE TEMP TABLE result AS SELECT * FROM ", macro,
                    "('source', column_map := ", column_map,
                    ", fasta_ref := ", quote_string(fasta), ")")
  } else if (case == "liftover") {
    limit <- 250000L * multiplier
    chain <- duckhts_bench_artifact_path("liftover_grch37_grch38_chain")
    dst <- duckhts_bench_artifact_path("liftover_grch38_fasta")
    src <- duckhts_bench_artifact_path("liftover_grch37_fasta")
    stopifnot(file.exists(chain), file.exists(dst), file.exists(src))
    dbExecute(con, paste0(
      "CREATE TEMP VIEW source AS SELECT CHROM AS chrom, POS AS pos, ",
      "REF AS ref, ALT[1] AS alt FROM read_bcf(", quote_string(input),
      ", samples := '', scan_mode := 'sequential') LIMIT ", limit))
    query <- paste0(
      "CREATE TEMP TABLE result AS SELECT * FROM duckdb_liftover(",
      "'source', 'chrom', 'pos', ref_col := 'ref', alt_col := 'alt', ",
      "chain_path := ", quote_string(chain), ", dst_fasta_ref := ",
      quote_string(dst), ", src_fasta_ref := ", quote_string(src), ")")
  } else {
    stop("unknown case: ", case)
  }
  profile_file <- tempfile(fileext = ".json")
  dbExecute(con, "PRAGMA enable_profiling='json'")
  dbExecute(con, paste("PRAGMA profiling_output=", quote_string(profile_file)))
  started <- proc.time()[["elapsed"]]
  if (startsWith(case, "r_")) {
    mapping <- c(CHR = "CHR", BP = "BP", SNP = "MarkerName",
                 A1 = "Allele1", A2 = "Allele2", P = "P-value",
                 Z = "Zscore", FRQ = "Freq1", NEFF = "Weight")
    if (case == "r_munge_metal") {
      mapping <- c(mapping, HET_I2 = "HetISq", HET_P = "HetPVal",
                   DIRE = "Direction")
    }
    result <- Rduckhts::rduckhts_munge(con, "source", column_map = mapping,
                                      fasta_ref = fasta)
    seconds <- proc.time()[["elapsed"]] - started
    rows <- nrow(result)
    rm(result)
  } else {
    dbExecute(con, query)
    seconds <- proc.time()[["elapsed"]] - started
    profile <- jsonlite::fromJSON(profile_file)
    saved_profile <- Sys.getenv("DUCKHTS_SCALING_PROFILE", "")
    if (nzchar(saved_profile)) file.copy(profile_file, saved_profile, overwrite = TRUE)
    rows <- dbGetQuery(con, "SELECT count(*) AS n FROM result")$n[[1L]]
  }
  dbExecute(con, "PRAGMA disable_profiling")
  if (startsWith(case, "r_")) {
    profile <- NULL
  }
  saved_result <- Sys.getenv("DUCKHTS_SCALING_RESULT", "")
  if (nzchar(saved_result)) {
    stopifnot(case %in% c("vcf_counts", "import_sites"))
    keys <- if (case == "vcf_counts") "sample_id, site_index" else "site_index"
    dbExecute(con, paste0("COPY (SELECT * FROM result ORDER BY ", keys, ") TO ",
                          quote_string(saved_result), " (FORMAT PARQUET)"))
  }
  rss <- readLines("/proc/self/status", warn = FALSE)
  peak_rss_kib <- as.numeric(gsub("[^0-9]", "", rss[startsWith(rss, "VmHWM:")]))
  input_records <- if (case %in% c("munge", "munge_metal", "r_munge",
                                  "r_munge_metal", "liftover")) {
    limit
  } else {
    as.integer(system2("bcftools", c("index", "-n", shQuote(input)),
                       stdout = TRUE))
  }
  write.table(data.frame(case, multiplier, threads, input_records, seconds,
                         peak_rss_kib,
                         peak_buffer_mib = if (is.null(profile)) NA_real_ else {
                           profile$system_peak_buffer_memory / 1048576
                         },
                         peak_spill_mib = if (is.null(profile)) NA_real_ else {
                           profile$system_peak_temp_dir_size / 1048576
                         },
                         result_rows = rows), stdout(), sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  dbDisconnect(con, shutdown = TRUE)
  unlink(profile_file)
} else {
  stopifnot(length(args) %in% c(1L, 2L), file.exists(args[[1L]]))
  wrapper_run <- length(args) == 2L && identical(args[[2L]], "--wrappers")
  stopifnot(length(args) == 1L || wrapper_run)
  extension <- normalizePath(args[[1L]], mustWork = TRUE)
  output <- if (wrapper_run) "benchmarks/sql_scaling_audit_wrappers.tsv" else {
    "benchmarks/sql_scaling_audit_real.tsv"
  }
  cache <- file.path(duckhts_bench_cache_dir(), "benchmarks/sql-scaling-audit")
  inputs <- if (wrapper_run) {
    c(r_munge = "epilepsy", r_munge_metal = "epilepsy")
  } else {
    c(vcf_counts = "giab", import_sites = "phase3", munge = "epilepsy",
      munge_metal = "epilepsy", liftover = "phase3_source")
  }
  header <- paste(c("case", "multiplier", "threads", "input_records", "seconds",
                    "peak_rss_kib", "peak_buffer_mib", "peak_spill_mib",
                    "result_rows"), collapse = "\t")
  writeLines(header, output)
  for (case in names(inputs)) {
    for (multiplier in c(1L, 2L, 4L)) {
      if (case %in% c("munge", "munge_metal", "r_munge", "r_munge_metal")) {
        input <- duckhts_bench_artifact_path("sql_scaling_epilepsy")
      } else if (case == "liftover") {
        input <- duckhts_bench_artifact_path("sql_scaling_1000g_phase3_chr22")
      } else {
        input <- file.path(cache, sprintf("%s-%dx.vcf.gz", inputs[[case]], multiplier))
        stopifnot(file.exists(paste0(input, ".tbi")))
      }
      stopifnot(file.exists(input))
      for (threads in c(1L, 4L)) {
        command <- c("benchmarks/benchmark_sql_scaling_real_run.R", "--child",
                     case, threads, multiplier, extension, input)
        stderr_file <- tempfile("sql-scaling-stderr-")
        result <- suppressWarnings(system2("Rscript", shQuote(command),
                                           stdout = TRUE, stderr = stderr_file))
        if (length(result) != 1L || !is.null(attr(result, "status"))) {
          stop(paste("Failed:", case, multiplier, threads,
                     paste(readLines(stderr_file, warn = FALSE), collapse = "\n")))
        }
        write(result, file = output, append = TRUE)
        message(result)
        unlink(stderr_file)
      }
    }
  }
}
