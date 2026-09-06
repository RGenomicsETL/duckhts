#' Stage a registry workload of committed fixtures without network access.
#'
#' @param repo DuckHTS checkout containing the registered fixture bytes.
#' @param workload Registry staging workload.
#' @return Named cache paths, indexed by artifact ID.
#' @export
duckhts_bench_stage_repository_fixtures <- function(repo, workload) {
  repo <- normalizePath(repo, winslash = "/", mustWork = TRUE)
  plan <- duckhts_bench_stage_plan(workload)
  stopifnot(all(plan$transform == "copy_committed_fixture"),
            all(grepl("^repo:test/data/[A-Za-z0-9._/-]+$", plan$locator)),
            !any(grepl("(^|/)\\.\\.?(/|$)", plan$locator)))
  paths <- stats::setNames(vapply(plan$id, duckhts_bench_artifact_path, character(1)), plan$id)
  for (i in seq_len(nrow(plan))) {
    source <- file.path(repo, sub("^repo:", "", plan$locator[[i]]))
    if (!file.exists(source)) stop("missing registered fixture: ", source, call. = FALSE)
    if (!startsWith(normalizePath(source, winslash = "/"), paste0(repo, "/test/data/"))) {
      stop("registered fixture resolves outside test/data: ", source, call. = FALSE)
    }
    duckhts_bench_validate_identity(plan$id[[i]], source)
    dir.create(dirname(paths[[i]]), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(source, paths[[i]], overwrite = TRUE))
    duckhts_bench_validate_identity(plan$id[[i]], paths[[i]])
    duckhts_bench_write_provenance(plan$id[[i]], paths[[i]])
  }
  paths
}

duckhts_bench_identity_fields <- function(identity) {
  identity <- as.character(identity)
  if (!length(identity) || is.na(identity) || !nzchar(identity)) return(character())
  fields <- strsplit(identity, ";", fixed = TRUE)[[1L]]
  if (any(!grepl("^[^=]+=", fields))) {
    stop("invalid semicolon-delimited artifact identity", call. = FALSE)
  }
  keys <- sub("=.*$", "", fields)
  if (anyDuplicated(keys)) stop("artifact identity fields must be unique", call. = FALSE)
  stats::setNames(sub("^[^=]+=", "", fields), keys)
}

#' Download one registry artifact into the DuckHTS cache.
#'
#' This handles only registry rows whose transform is `direct_download`; derived
#' artifacts must be produced by a named workload stage function.
#'
#' @param id Artifact identifier.
#' @param overwrite Whether to replace an existing cache artifact.
#' @param output Optional destination; defaults to the registry cache path.
#' @return The staged path, invisibly.
#' @export
duckhts_bench_fetch <- function(id, overwrite = FALSE, output = duckhts_bench_artifact_path(id)) {
  registry <- duckhts_bench_registry()
  row <- registry[registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("unknown or non-unique benchmark artifact: ", id, call. = FALSE)
  if (row$transform != "direct_download") {
    stop("artifact requires a named derivation, not direct download: ", id, call. = FALSE)
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  ready <- file.exists(output) && file.info(output)$size > 0L
  if (ready && !overwrite) {
    valid <- tryCatch({
      duckhts_bench_validate_identity(id, output)
      TRUE
    }, error = function(error) FALSE)
    if (!valid) {
      unlink(c(output, paste0(output, ".provenance.tsv")), force = TRUE)
      ready <- FALSE
    }
  }
  if (overwrite || !ready) {
    values <- duckhts_bench_identity_fields(row$supplier_identity)
    if ("etag" %in% names(values)) {
      curl <- Sys.which("curl")
      if (!nzchar(curl)) stop("curl is required to validate supplier ETag: ", id, call. = FALSE)
      headers <- system2(curl, c("-fsSIL", row$locator), stdout = TRUE)
      etags <- sub("^[[:space:]]*[Ee][Tt][Aa][Gg]:[[:space:]]*\\\"?([^\\\"[:space:]\\r]+).*", "\\1", headers)
      etags <- etags[grepl("^[[:space:]]*[Ee][Tt][Aa][Gg]:", headers)]
      if (!length(etags) || utils::tail(etags, 1L) != values[["etag"]]) {
        stop("supplier ETag does not match registered artifact: ", id, call. = FALSE)
      }
    }
    temporary <- paste0(output, ".partial-", Sys.getpid())
    unlink(temporary, force = TRUE)
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    old_timeout <- getOption("timeout")
    options(timeout = max(3600L, old_timeout))
    on.exit(options(timeout = old_timeout), add = TRUE)
    utils::download.file(row$locator, temporary, mode = "wb", quiet = FALSE)
    duckhts_bench_validate_identity(id, temporary)
    if (file.exists(output) && unlink(output, force = TRUE) != 0L) stop("could not replace cached artifact: ", id, call. = FALSE)
    if (!file.rename(temporary, output)) stop("could not publish downloaded artifact: ", id, call. = FALSE)
  }
  duckhts_bench_write_provenance(id, output)
  invisible(output)
}

#' Write the registry receipt adjacent to one staged artifact.
#'
#' @param id Artifact identifier.
#' @param output Staged file path.
#' @return The receipt path, invisibly.
#' @export
#' Validate a staged file against supplier-published identity metadata.
#'
#' @param id Artifact identifier.
#' @param output Staged file path.
#' @return Invisibly `TRUE`.
#' @export
duckhts_bench_validate_identity <- function(id, output = duckhts_bench_artifact_path(id)) {
  row <- duckhts_bench_registry()
  row <- row[row$id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("unknown or non-unique benchmark artifact: ", id, call. = FALSE)
  values <- duckhts_bench_identity_fields(row$supplier_identity)
  if (!length(values)) return(invisible(TRUE))
  if ("bytes" %in% names(values) && file.info(output)$size != as.numeric(values[["bytes"]])) {
    stop("supplier byte identity does not match cached artifact: ", id, call. = FALSE)
  }
  if ("md5" %in% names(values) && unname(tools::md5sum(output)) != values[["md5"]]) {
    stop("supplier MD5 identity does not match cached artifact: ", id, call. = FALSE)
  }
  if ("sha256" %in% names(values)) {
    sha256_bin <- Sys.which("sha256sum")
    sha256_args <- shQuote(output)
    if (!nzchar(sha256_bin)) {
      sha256_bin <- Sys.which("shasum")
      sha256_args <- c("-a", "256", shQuote(output))
    }
    if (!nzchar(sha256_bin)) stop("sha256sum or shasum is required to validate supplier checksum: ", id, call. = FALSE)
    actual_sha256 <- strsplit(trimws(system2(sha256_bin, sha256_args, stdout = TRUE)), "[[:space:]]+")[[1L]][[1L]]
    if (actual_sha256 != values[["sha256"]]) stop("supplier SHA-256 identity does not match cached artifact: ", id, call. = FALSE)
  }
  if ("sum" %in% names(values) && "blocks" %in% names(values)) {
    sum_bin <- Sys.which("sum")
    if (!nzchar(sum_bin)) stop("sum is required to validate supplier checksum: ", id, call. = FALSE)
    fields <- strsplit(trimws(system2(sum_bin, c("-r", shQuote(output)), stdout = TRUE)), "[[:space:]]+")[[1L]]
    if (length(fields) < 2L || as.integer(fields[[1L]]) != as.integer(values[["sum"]]) ||
        as.integer(fields[[2L]]) != as.integer(values[["blocks"]])) {
      stop("supplier sum identity does not match cached artifact: ", id, call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Write the registry receipt adjacent to one staged artifact.
#'
#' @param id Artifact identifier.
#' @param output Staged file path.
#' @return The receipt path, invisibly.
#' @export
duckhts_bench_write_provenance <- function(id, output = duckhts_bench_artifact_path(id)) {
  registry <- duckhts_bench_registry()
  row <- registry[registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("unknown or non-unique benchmark artifact: ", id, call. = FALSE)
  receipt <- paste0(output, ".provenance.tsv")
  fields <- data.frame(
    field = c("artifact_id", "workload", "release", "source_locator", "access", "transform", "supplier_identity", "cached_output", "consumer"),
    value = c(id, row$workload, row$release, row$locator, row$access, row$transform, row$supplier_identity, output, row$consumer),
    stringsAsFactors = FALSE
  )
  dir.create(dirname(receipt), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(fields, receipt, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(receipt)
}
