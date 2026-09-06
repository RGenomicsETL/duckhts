#' Stage the committed small BCF scan benchmark inputs without network access.
#'
#' The registry pins the encoded files, not tool-version-dependent re-encodings.
#' Fixture reconstruction remains the repository's maintainer operation.
#' @param repo DuckHTS checkout containing the registered fixture bytes.
#' @return Named cache paths, indexed by artifact ID.
#' @export
duckhts_bench_stage_bcf_scan <- function(repo) {
  repo <- normalizePath(repo, mustWork = TRUE)
  plan <- duckhts_bench_stage_plan("bcf-scan-init")
  stopifnot(all(plan$transform == "copy_committed_fixture"),
            all(grepl("^repo:test/data/[^/]+$", plan$locator)))
  paths <- stats::setNames(vapply(plan$id, duckhts_bench_artifact_path, character(1)), plan$id)
  for (i in seq_len(nrow(plan))) {
    source <- file.path(repo, sub("^repo:", "", plan$locator[[i]]))
    if (!file.exists(source)) stop("missing registered fixture: ", source, call. = FALSE)
    duckhts_bench_validate_identity(plan$id[[i]], source)
    dir.create(dirname(paths[[i]]), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(source, paths[[i]], overwrite = TRUE))
    duckhts_bench_validate_identity(plan$id[[i]], paths[[i]])
    duckhts_bench_write_provenance(plan$id[[i]], paths[[i]])
  }
  paths
}
