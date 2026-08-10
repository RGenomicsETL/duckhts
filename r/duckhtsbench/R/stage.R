#' Download one registry artifact into the DuckHTS cache.
#'
#' This handles only registry rows whose transform is `direct_download`; derived
#' artifacts must be produced by a named workload stage function.
#'
#' @param id Artifact identifier.
#' @param overwrite Whether to replace an existing cache artifact.
#' @return The staged path, invisibly.
#' @export
duckhts_bench_fetch <- function(id, overwrite = FALSE) {
  registry <- duckhts_bench_registry()
  row <- registry[registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("unknown or non-unique benchmark artifact: ", id, call. = FALSE)
  if (row$transform != "direct_download") {
    stop("artifact requires a named derivation, not direct download: ", id, call. = FALSE)
  }
  output <- duckhts_bench_artifact_path(id)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  if (overwrite || !file.exists(output) || !file.info(output)$size) {
    utils::download.file(row$locator, output, mode = "wb", quiet = FALSE)
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
