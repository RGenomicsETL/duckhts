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
  registry <- duckhts_bench_registry()
  row <- registry[registry$id == "norm_hg00096_chr22_20m_30m", , drop = FALSE]
  if (nrow(row) != 1L || row$locator != "artifact:norm_hg00096_raw_gvcf") stop("normalization slice registry entry is missing or has an unregistered source", call. = FALSE)
  raw <- duckhts_bench_fetch("norm_hg00096_raw_gvcf")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(output, ".partial-", Sys.getpid())
  unlink(temporary, force = TRUE)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  status <- system2(bcftools, c("view", "--threads", threads, "-r", "chr22:20000000-30000000", "-s", "HG00096", "-Oz", "-o", shQuote(temporary), shQuote(raw)))
  if (status != 0L || !file.exists(temporary) || !file.info(temporary)$size) stop("could not extract registered 1000G DRAGEN gVCF slice", call. = FALSE)
  if (file.exists(output) && unlink(output, force = TRUE) != 0L) stop("could not replace normalization gVCF slice", call. = FALSE)
  if (!file.rename(temporary, output)) stop("could not publish registered normalization gVCF slice", call. = FALSE)
  index <- paste0(output, ".tbi")
  if (file.exists(index) && unlink(index, force = TRUE) != 0L) stop("could not replace normalization gVCF index", call. = FALSE)
  status <- system2(tabix, c("-f", "-p", "vcf", shQuote(output)))
  if (status != 0L || !file.exists(index) || !file.info(index)$size) stop("could not index registered normalization gVCF slice", call. = FALSE)
  duckhts_bench_write_provenance("norm_hg00096_chr22_20m_30m", output)
  invisible(output)
}
