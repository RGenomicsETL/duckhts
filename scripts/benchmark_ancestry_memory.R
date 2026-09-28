#!/usr/bin/env Rscript
# Measure the keyed Parquet ancestry query in an independent R process.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L || !args[[1L]] %in% c("epilepsy", "genotypes")) {
  stop("Usage: benchmark_ancestry_memory.R epilepsy|genotypes input.parquet reference.parquet threads scale")
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
con <- rduckhts_connect()
dbExecute(con, paste0("SET threads=", threads))
dbExecute(con, "SET memory_limit='1GB'")
dbExecute(con, "SET preserve_insertion_order=false")
spill <- tempfile("ancestry-spill-")
quoted <- function(x) as.character(dbQuoteString(con, x))
dbExecute(con, paste0("SET temp_directory=", quoted(spill)))
dbExecute(con, paste0("CREATE VIEW real_reference AS SELECT * FROM read_parquet(",
                      quoted(reference_path), ")"))
if (workload == "epilepsy") {
  n <- min(3343235L, 835809L * scale)
  dbExecute(con, paste0("CREATE VIEW real_input AS SELECT * FROM read_parquet(",
                        quoted(input_path), ") LIMIT ", n))
} else {
  dbExecute(con, paste0("CREATE VIEW real_input AS SELECT * FROM read_parquet(",
                        quoted(input_path), ") WHERE sample_id IN (",
                        "SELECT DISTINCT sample_id FROM read_parquet(", quoted(input_path),
                        ") ORDER BY sample_id LIMIT ", 11L * scale, ")"))
}
correction_path <- duckhts_bench_stage_repository_fixtures(
  normalizePath("."), "ancestry-projection")[["ancestry_correction"]]
correction <- read.delim(correction_path)
dbWriteTable(con, "real_correction", correction)
input_rows <- dbGetQuery(con, "SELECT count(*) AS n FROM real_input")$n
profile <- tempfile("ancestry-profile-", fileext = ".json")
dbExecute(con, "PRAGMA enable_profiling='json'")
dbExecute(con, paste0("PRAGMA profiling_output=", quoted(profile)))
start <- proc.time()[["elapsed"]]
result <- rduckhts_ancestry_proportions_wide(
  con, "real_input", "real_reference", "real_correction")
seconds <- proc.time()[["elapsed"]] - start
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
                     output_rows = nrow(result), seconds = seconds,
                     peak_rss_mib = rss_mib,
                     peak_buffer_mib = metrics$system_peak_buffer_memory / 1048576,
                     peak_temp_mib = metrics$system_peak_temp_dir_size / 1048576)
write.table(output, stdout(), sep = "\t", row.names = FALSE, quote = FALSE)
dbDisconnect(con, shutdown = TRUE)
unlink(c(profile, spill), recursive = TRUE)
