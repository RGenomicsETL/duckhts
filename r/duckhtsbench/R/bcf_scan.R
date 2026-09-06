#' Stage the committed small BCF scan benchmark inputs without network access.
#'
#' The registry pins the encoded files, not tool-version-dependent re-encodings.
#' Fixture reconstruction remains the repository's maintainer operation.
#' @param repo DuckHTS checkout containing the registered fixture bytes.
#' @return Named cache paths, indexed by artifact ID.
#' @export
duckhts_bench_stage_bcf_scan <- function(repo) {
  duckhts_bench_stage_repository_fixtures(repo, "bcf-scan-init")
}
