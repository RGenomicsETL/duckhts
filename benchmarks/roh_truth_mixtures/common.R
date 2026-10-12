# Shared pieces of the planted-truth and read-mixture evaluation
# (benchmark_roh_truth_and_mixtures.Rmd): registered inputs, their identity
# checks, fresh-process queries and the memory gate.
#
# Run from the repository root.

source("benchmarks/roh_truth_mixtures/declaration.R")
source("benchmarks/roh_counts_validation/fitted_emission.R")
source("benchmarks/bam_mismatch_stage.R")

roh_truth_quote <- roh_fitted_quote

roh_truth_registered <- function(id) {
  registry <- duckhtsbench::duckhts_bench_registry()
  row <- registry[registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("unknown or non-unique benchmark artifact: ", id, call. = FALSE)
  duckhtsbench:::duckhts_bench_identity_fields(row$supplier_identity[[1L]])
}

roh_truth_sha256 <- function(path) {
  value <- system2("sha256sum", shQuote(path), stdout = TRUE)
  if (!is.null(attr(value, "status")) || length(value) != 1L) {
    stop("could not hash ", path, call. = FALSE)
  }
  sub(" .*$", "", value)
}

# A staged file whose registered identity holds a sha256, bytes or both.
roh_truth_check_file <- function(id) {
  path <- duckhtsbench::duckhts_bench_artifact_path(id)
  if (!file.exists(path)) stop("artifact ", id, " is not staged: ", path, call. = FALSE)
  registered <- roh_truth_registered(id)
  checked <- intersect(c("sha256", "bytes"), names(registered))
  if (!length(checked)) stop("artifact ", id, " registers neither sha256 nor bytes", call. = FALSE)
  if ("bytes" %in% checked && file.size(path) != as.numeric(registered[["bytes"]])) {
    stop("artifact ", id, " differs from its registered bytes", call. = FALSE)
  }
  if ("sha256" %in% checked && !identical(roh_truth_sha256(path), registered[["sha256"]])) {
    stop("artifact ", id, " differs from its registered sha256", call. = FALSE)
  }
  path
}

roh_truth_connect <- function(extension, threads = 4L) {
  driver <- duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE, allow_extensions = TRUE,
    config = list(threads = as.character(threads), allow_unsigned_extensions = "true",
                  autoinstall_known_extensions = "false", autoload_known_extensions = "false"))
  con <- DBI::dbConnect(driver)
  DBI::dbExecute(con, sprintf("LOAD %s", roh_truth_quote(extension)))
  con
}

# A recorder of fresh-process queries: each call runs one query in a new CLI
# process (roh_fitted_run_fresh) and keeps its observation with its ceiling.
# The ceiling is the overhead, the median peak RSS of three empty runs, plus
# 3 x the declared live state in bytes.
roh_truth_recorder <- function(cli, extension, work) {
  observations <- list()
  record <- function(name, kind, setup, query, threads, output, format = "csv",
                     live_bytes = NA_real_, overhead_mib = NA_real_) {
    measured <- roh_fitted_run_fresh(cli, extension, setup, query, threads, output, format)
    ceiling_mib <- overhead_mib + 3 * live_bytes / 2^20
    observations[[length(observations) + 1L]] <<- data.frame(
      observation = name, kind = kind, threads = threads, measured,
      live_mib = live_bytes / 2^20, ceiling_mib = ceiling_mib)
    invisible(measured)
  }
  overhead <- function(threads) {
    for (repetition in 1:3) {
      record(sprintf("empty_%d_threads_%d", repetition, threads), "empty", character(0),
             "SELECT 1 AS n", threads, file.path(work, "empty.csv"))
    }
    empty <- do.call(rbind, observations)
    stats::median(empty$max_rss_mib[empty$kind == "empty" & empty$threads == threads])
  }
  table <- function() {
    if (!length(observations)) return(NULL)
    result <- do.call(rbind, observations)
    result$within_ceiling <- is.na(result$ceiling_mib) | result$max_rss_mib <= result$ceiling_mib
    result$no_spill <- result$peak_temp_bytes == 0
    result
  }
  list(record = record, overhead = overhead, table = table)
}

# Stop when any recorded query spilled or passed its ceiling.
roh_truth_gate <- function(observations) {
  bad <- observations[!observations$no_spill | !observations$within_ceiling, , drop = FALSE]
  if (nrow(bad)) {
    print(bad[c("observation", "max_rss_mib", "ceiling_mib", "peak_temp_bytes")])
    stop("a query spilled or passed its memory ceiling", call. = FALSE)
  }
  invisible(TRUE)
}

roh_truth_read_runs <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE,
                  colClasses = c(sample = "character", chrom = "character"))
}
