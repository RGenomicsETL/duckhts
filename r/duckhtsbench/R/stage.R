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
    temporary <- paste0(output, ".partial-", Sys.getpid())
    unlink(temporary, force = TRUE)
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
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
  identity <- as.character(row$supplier_identity)
  if (!length(identity) || is.na(identity) || !nzchar(identity)) return(invisible(TRUE))
  fields <- strsplit(identity, ";", fixed = TRUE)[[1L]]
  values <- stats::setNames(sub("^[^=]+=", "", fields), sub("=.*$", "", fields))
  if ("bytes" %in% names(values) && file.info(output)$size != as.numeric(values[["bytes"]])) {
    stop("supplier byte identity does not match cached artifact: ", id, call. = FALSE)
  }
  if ("md5" %in% names(values) && unname(tools::md5sum(output)) != values[["md5"]]) {
    stop("supplier MD5 identity does not match cached artifact: ", id, call. = FALSE)
  }
  if ("sum" %in% names(values) && "blocks" %in% names(values)) {
    sum_bin <- Sys.which("sum")
    if (!nzchar(sum_bin)) stop("sum is required to validate supplier checksum: ", id, call. = FALSE)
    fields <- strsplit(trimws(system2(sum_bin, c("-r", output), stdout = TRUE)), "[[:space:]]+")[[1L]]
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
    field = c("artifact_id", "workload", "release", "source_locator", "access", "transform", "cached_output", "consumer"),
    value = c(id, row$workload, row$release, row$locator, row$access, row$transform, output, row$consumer),
    stringsAsFactors = FALSE
  )
  dir.create(dirname(receipt), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(fields, receipt, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(receipt)
}
