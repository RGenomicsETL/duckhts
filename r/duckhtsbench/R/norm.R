#' Stage the pinned 1000 Genomes DRAGEN gVCF normalization slice.
#'
#' @param output Destination VCF, defaulting to the registered cache path.
#' @param bcftools Path to `bcftools`.
#' @param tabix Path to `tabix`.
#' @param threads bcftools worker threads.
#' @return The staged VCF path, invisibly.
#' @export
duckhts_bench_stage_norm <- function(
    output = duckhts_bench_artifact_path("norm_hg00096_chr22_20m_30m"),
    bcftools = Sys.which("bcftools"), tabix = Sys.which("tabix"), threads = 2L) {
  if (!nzchar(bcftools) || !nzchar(tabix)) stop("bcftools and tabix are required to stage normalization input", call. = FALSE)
  row <- duckhts_bench_registry()
  row <- row[row$id == "norm_hg00096_chr22_20m_30m", , drop = FALSE]
  if (nrow(row) != 1L) stop("normalization registry entry is missing", call. = FALSE)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(output) || !file.info(output)$size) {
    status <- system2(bcftools, c("view", "--threads", threads, "-r", "chr22:20000000-30000000", "-s", "HG00096", "-Oz", "-o", output, row$locator))
    if (status != 0L) stop("could not extract registered 1000G DRAGEN gVCF slice", call. = FALSE)
  }
  if (!file.exists(paste0(output, ".tbi")) || !file.info(paste0(output, ".tbi"))$size) {
    status <- system2(tabix, c("-f", "-p", "vcf", output))
    if (status != 0L) stop("could not index registered normalization gVCF slice", call. = FALSE)
  }
  duckhts_bench_write_provenance("norm_hg00096_chr22_20m_30m", output)
  invisible(output)
}
