#' Stage selected GIAB NIST v4.2.1 GRCh37 benchmark VCFs.
#'
#' @param samples Character sample identifiers, currently HG001, HG002, or HG006.
#' @param bcftools Path to `bcftools`.
#' @return Named character vector of cache VCF paths, invisibly.
#' @export
duckhts_bench_stage_giab <- function(samples = c("HG001", "HG002", "HG006"), bcftools = Sys.which("bcftools")) {
  if (!nzchar(bcftools)) stop("bcftools is required to index staged GIAB VCFs", call. = FALSE)
  ids <- paste0("giab_", tolower(samples), "_grch37_v421")
  registry <- duckhts_bench_registry()
  if (any(!ids %in% registry$id)) {
    stop("unsupported GIAB sample(s): ", paste(samples[!ids %in% registry$id], collapse = ", "), call. = FALSE)
  }
  paths <- vapply(ids, duckhts_bench_fetch, character(1L))
  names(paths) <- samples
  for (path in paths) {
    index <- paste0(path, ".tbi")
    if (file.exists(index) && unlink(index, force = TRUE) != 0L) stop("could not replace GIAB VCF index: ", path, call. = FALSE)
    status <- system2(bcftools, c("index", "-f", "-t", shQuote(path)))
    if (status != 0L || !file.exists(index) || !file.info(index)$size) stop("could not index staged GIAB VCF: ", path, call. = FALSE)
  }
  invisible(paths)
}
