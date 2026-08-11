#' Stage the GRCh37-to-GRCh38 liftover reference bundle.
#'
#' @param base_dir Destination directory, defaulting to the registered cache path.
#' @param samtools Path to `samtools`.
#' @param gzip Path to `gzip`.
#' @return A named character vector containing source FASTA, destination FASTA,
#'   and chain paths, invisibly.
#' @export
duckhts_bench_stage_liftover <- function(
    base_dir = duckhts_bench_cache_path("references/liftover/grch37-to-grch38"),
    samtools = Sys.which("samtools"), gzip = Sys.which("gzip")) {
  if (!nzchar(samtools) || !nzchar(gzip)) stop("samtools and gzip are required for liftover staging", call. = FALSE)
  source_gz <- duckhts_bench_fetch("liftover_grch37_fasta_gz", output = file.path(base_dir, "downloads", "human_g1k_v37.fasta.gz"))
  source_fasta <- file.path(base_dir, "human_g1k_v37.fasta")
  temporary <- paste0(source_fasta, ".partial-", Sys.getpid())
  unlink(temporary, force = TRUE)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  status <- system2(gzip, c("--decompress", "--keep", "--stdout", source_gz), stdout = temporary)
  if (status != 0L) stop("could not decompress the GRCh37 liftover FASTA", call. = FALSE)
  if (file.exists(source_fasta) && unlink(source_fasta, force = TRUE) != 0L) {
    stop("could not replace the decompressed GRCh37 liftover FASTA", call. = FALSE)
  }
  if (!file.rename(temporary, source_fasta)) stop("could not publish the decompressed GRCh37 FASTA", call. = FALSE)
  destination_fasta <- duckhts_bench_fetch("liftover_grch38_fasta", output = file.path(base_dir, "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"))
  chain <- duckhts_bench_fetch("liftover_grch37_grch38_chain", output = file.path(base_dir, "GRCh37_to_GRCh38.chain.gz"))
  for (fasta in c(source_fasta, destination_fasta)) {
    index <- paste0(fasta, ".fai")
    if (file.exists(index) && unlink(index, force = TRUE) != 0L) {
      stop("could not replace liftover FASTA index: ", index, call. = FALSE)
    }
    status <- system2(samtools, c("faidx", shQuote(fasta)))
    if (status != 0L || !file.exists(index) || !file.info(index)$size) {
      stop("could not index liftover FASTA: ", fasta, call. = FALSE)
    }
  }
  paths <- c(source_fasta = source_fasta, destination_fasta = destination_fasta, chain = chain)
  writeLines(c("field\tvalue", "workload\tliftover_reference_bundle", paste(names(paths), paths, sep = "\t")), file.path(base_dir, "provenance.tsv"))
  invisible(paths)
}
