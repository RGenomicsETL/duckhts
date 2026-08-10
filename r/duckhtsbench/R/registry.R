duckhts_bench_cache_dir <- function() {
  configured <- Sys.getenv("DUCKHTS_CACHE_DIR", unset = "")
  if (nzchar(configured)) return(path.expand(configured))
  xdg <- Sys.getenv("XDG_CACHE_HOME", unset = "")
  if (nzchar(xdg)) return(file.path(path.expand(xdg), "duckhts"))
  file.path(path.expand("~"), ".cache", "duckhts")
}

#' Return a path below the configured DuckHTS cache.
#'
#' @param relative_path Relative cache path.
#' @return A cache path.
#' @export
duckhts_bench_cache_path <- function(relative_path = "") {
  if (length(relative_path) != 1L || is.na(relative_path) || grepl("^(/|[A-Za-z]:[/\\\\])", relative_path)) {
    stop("relative_path must be one relative cache path", call. = FALSE)
  }
  file.path(duckhts_bench_cache_dir(), relative_path)
}

bench_registry_path <- function() {
  path <- Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = "")
  if (nzchar(path) && file.exists(path)) return(path)
  installed <- system.file("benchmark_registry.tsv", package = "duckhtsbench")
  if (nzchar(installed)) return(installed)
  stop("duckhtsbench must be installed or DUCKHTSBENCH_REGISTRY must name benchmark_registry.tsv", call. = FALSE)
}

#' Read the DuckHTS benchmark input authority.
#'
#' @return A data frame with one row per source or derived artifact.
#' @export
duckhts_bench_registry <- function() {
  utils::read.delim(bench_registry_path(), stringsAsFactors = FALSE, check.names = FALSE)
}

#' Return the cache path for one benchmark artifact.
#'
#' @param id Artifact identifier from [duckhts_bench_registry()].
#' @return A cache path.
#' @export
duckhts_bench_artifact_path <- function(id) {
  registry <- duckhts_bench_registry()
  row <- registry[registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("unknown or non-unique benchmark artifact: ", id, call. = FALSE)
  file.path(duckhts_bench_cache_dir(), row$cache_relpath[[1L]])
}

#' Return a reproducible staging plan for a benchmark workload.
#'
#' @param workload Registry workload name.
#' @return Rows from the benchmark registry in staging order.
#' @export
duckhts_bench_stage_plan <- function(workload) {
  registry <- duckhts_bench_registry()
  rows <- registry[registry$workload == workload, , drop = FALSE]
  if (!nrow(rows)) stop("unknown benchmark workload: ", workload, call. = FALSE)
  rows[order(rows$stage_order), , drop = FALSE]
}

#' Reject machine-local paths in executable benchmark chunks.
#'
#' @param benchmark_dir Directory containing benchmark R Markdown sources.
#' @return Invisibly `TRUE`; errors name every active chunk with a `/root/` path.
#' @export
duckhts_bench_check_portability <- function(benchmark_dir = "benchmarks") {
  paths <- sort(list.files(benchmark_dir, pattern = "\\.Rmd$", full.names = TRUE))
  errors <- character()
  for (path in paths) {
    lines <- readLines(path, warn = FALSE)
    active <- FALSE
    for (line_number in seq_along(lines)) {
      line <- lines[[line_number]]
      chunk <- regmatches(line, regexec("^```\\{([^}]*)}\\s*$", line))[[1L]]
      if (length(chunk)) {
        header <- chunk[[2L]]
        engine <- sub("^\\s*([^,[:space:]]+).*$", "\\1", header)
        active <- engine %in% c("r", "bash", "sh", "python", "perl") &&
          !grepl("\\beval\\s*=\\s*FALSE\\b", header, ignore.case = TRUE, perl = TRUE)
        next
      }
      if (startsWith(line, "```")) {
        active <- FALSE
        next
      }
      if (active && grepl("/root/", line, fixed = TRUE)) {
        errors <- c(errors, sprintf("%s:%d: machine-local path in an active chunk", path, line_number))
      }
    }
  }
  if (length(errors)) stop(paste(errors, collapse = "\n"), call. = FALSE)
  message("benchmark portability: OK")
  invisible(TRUE)
}
