#!/usr/bin/env Rscript
# Measure setup, first query and warm query in an independent R process.
script_start <- proc.time()[["elapsed"]]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L || !args[[1L]] %in% c("epilepsy", "genotypes", "joint")) {
  stop("Usage: benchmark_ancestry_memory.R epilepsy|genotypes|joint input.parquet reference.parquet threads scale")
}
workload <- args[[1L]]
input_path <- normalizePath(args[[2L]], mustWork = TRUE)
reference_path <- normalizePath(args[[3L]], mustWork = TRUE)
threads <- as.integer(args[[4L]])
scale <- as.integer(args[[5L]])
if (!threads %in% c(1L, 4L) || !scale %in% c(1L, 2L, 4L)) {
  stop("threads must be 1 or 4; scale must be 1, 2 or 4")
}
library(DBI)
library(Rduckhts)
library(duckhtsbench)
library_seconds <- proc.time()[["elapsed"]] - script_start
con <- rduckhts_connect()
extension_seconds <- proc.time()[["elapsed"]] - script_start - library_seconds
dbExecute(con, paste0("SET threads=", threads))
dbExecute(con, "SET memory_limit='1GB'")
dbExecute(con, "SET preserve_insertion_order=false")
spill <- tempfile("ancestry-spill-")
quoted <- function(x) as.character(dbQuoteString(con, x))
dbExecute(con, paste0("SET temp_directory=", quoted(spill)))
dbExecute(con, paste0("CREATE VIEW real_reference AS SELECT * FROM read_parquet(",
                      quoted(reference_path), ")"))
reference_open_seconds <- proc.time()[["elapsed"]] - script_start -
  library_seconds - extension_seconds
if (workload == "epilepsy") {
  input_scan <- if (scale == 4L) {
    paste0("SELECT * FROM read_parquet(", quoted(input_path), ") LIMIT 3343235")
  } else {
    paste0("SELECT * EXCLUDE (file_row_number) FROM read_parquet(",
           quoted(input_path), ", file_row_number=true) WHERE file_row_number < ",
           835809L * scale)
  }
  dbExecute(con, paste0("CREATE VIEW real_input AS ", input_scan))
} else {
  site_limit <- if (workload == "joint") paste0(" AND position IN (SELECT position ",
    "FROM read_parquet(", quoted(input_path), ") GROUP BY position ",
    "ORDER BY position LIMIT ", 5000L * scale, ")") else ""
  dbExecute(con, paste0("CREATE VIEW real_input AS SELECT * FROM read_parquet(",
                        quoted(input_path), ") WHERE sample_id IN (",
                        "SELECT DISTINCT sample_id FROM read_parquet(", quoted(input_path),
                        ") ORDER BY sample_id LIMIT ", 11L * scale, ")", site_limit))
}
correction_path <- duckhts_bench_stage_repository_fixtures(
  normalizePath("."), "ancestry-projection")[["ancestry_correction"]]
correction <- read.delim(correction_path)
dbWriteTable(con, "real_correction", correction)
input_rows <- dbGetQuery(con, "SELECT count(*) AS n FROM real_input")$n
profile <- tempfile("ancestry-profile-", fileext = ".json")
dbExecute(con, "PRAGMA enable_profiling='json'")
dbExecute(con, paste0("PRAGMA profiling_output=", quoted(profile)))
setup_seconds <- proc.time()[["elapsed"]] - script_start
profile_steps <- Sys.getenv("DUCKHTS_ANCESTRY_PROFILE_STEPS", unset = "")
if (nzchar(profile_steps)) {
  step_times <- numeric()
  options(duckhts.ancestry_step_hook = function(stage) {
    if (!is.character(stage) || length(stage) != 1L || is.na(stage) ||
        !nzchar(stage) || stage %in% names(step_times)) {
      stop("invalid or repeated ancestry profile stage")
    }
    step_times <<- c(step_times, stats::setNames(proc.time()[["elapsed"]], stage))
  })
}
start <- proc.time()[["elapsed"]]
first <- rduckhts_ancestry_proportions(
  con, "real_input", "real_reference", "real_correction")
first_query_seconds <- proc.time()[["elapsed"]] - start
if (nzchar(profile_steps)) {
  options(duckhts.ancestry_step_hook = NULL)
  stopifnot(length(step_times) > 0L, all(is.finite(step_times)),
            !is.unsorted(step_times))
  step_table <- data.frame(stage = c(names(step_times), "returned"),
                           seconds_from_query = c(step_times,
                             proc.time()[["elapsed"]]) - start)
  utils::write.table(step_table, profile_steps, sep = "\t",
                     row.names = FALSE, quote = FALSE)
}
start <- proc.time()[["elapsed"]]
result <- rduckhts_ancestry_proportions(
  con, "real_input", "real_reference", "real_correction")
query_seconds <- proc.time()[["elapsed"]] - start
stopifnot(identical(first[c("sample_id", "group_id", "status", "proportion")],
                    result[c("sample_id", "group_id", "status", "proportion")]))
repeat_difference <- max(abs(first$proportion - result$proportion))
stopifnot(nrow(result) == 21L * if (workload == "epilepsy") 1L else 11L * scale,
          all(result$status == "ok"),
          length(unique(result$sample_id)) == if (workload == "epilepsy") 1L else 11L * scale)
proportion_path <- Sys.getenv("DUCKHTS_ANCESTRY_PROPORTIONS")
if (nzchar(proportion_path)) {
  utils::write.table(result[c("sample_id", "group_id", "proportion")], proportion_path,
                     sep = "\t", row.names = FALSE, quote = FALSE)
}
metrics <- jsonlite::fromJSON(profile, simplifyVector = FALSE)
rss <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
rss_mib <- as.numeric(strsplit(trimws(rss), "[[:space:]]+")[[1L]][[2L]]) / 1024
output <- data.frame(workload = workload, threads = threads, scale = scale,
                     samples = if (workload == "epilepsy") 1L else 11L * scale,
                     input_rows = input_rows, used_variants = sum(unique(
                       result[c("sample_id", "used_variants")])$used_variants),
                     output_rows = nrow(result), query_seconds = query_seconds,
                     first_query_seconds = first_query_seconds,
                     repeat_max_difference = repeat_difference,
                     setup_seconds = setup_seconds,
                     library_seconds = library_seconds,
                     extension_seconds = extension_seconds,
                     reference_open_seconds = reference_open_seconds,
                     other_setup_seconds = setup_seconds - library_seconds -
                       extension_seconds - reference_open_seconds,
                     peak_rss_mib = rss_mib,
                     aligned_file_mib = as.numeric(attr(result, "aligned_bytes")) / 1048576,
                     peak_buffer_mib = metrics$system_peak_buffer_memory / 1048576,
                     peak_temp_mib = metrics$system_peak_temp_dir_size / 1048576)
dbDisconnect(con, shutdown = TRUE)
unlink(c(profile, spill), recursive = TRUE)
output$cleanup_seconds <- proc.time()[["elapsed"]] - script_start -
  setup_seconds - first_query_seconds - query_seconds
output$script_seconds <- proc.time()[["elapsed"]] - script_start
write.table(output, stdout(), sep = "\t", row.names = FALSE, quote = FALSE)
