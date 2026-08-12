#' Stage the GRCh37-to-GRCh38 liftover reference bundle.
#'
#' @param base_dir Destination directory, defaulting to the registered cache path.
#' @param samtools Path to `samtools`.
#' @param gzip Path to `gzip`.
#' @return A named character vector containing source FASTA, destination FASTA,
#'   and chain paths, invisibly.
#' @export
duckhts_bench_stage_liftover <- function(
    base_dir = NULL, samtools = Sys.which("samtools"), gzip = Sys.which("gzip")) {
  if (!nzchar(samtools) || !nzchar(gzip)) stop("samtools and gzip are required for liftover staging", call. = FALSE)
  use_registry_paths <- is.null(base_dir)
  source_gz <- if (use_registry_paths) {
    duckhts_bench_fetch("liftover_grch37_fasta_gz")
  } else {
    duckhts_bench_fetch("liftover_grch37_fasta_gz", output = file.path(base_dir, "downloads", "human_g1k_v37.fasta.gz"))
  }
  source_fasta <- if (use_registry_paths) duckhts_bench_artifact_path("liftover_grch37_fasta") else file.path(base_dir, "human_g1k_v37.fasta")
  if (use_registry_paths) base_dir <- dirname(source_fasta)
  registry <- duckhts_bench_registry()
  derived <- registry[registry$id == "liftover_grch37_fasta", , drop = FALSE]
  if (nrow(derived) != 1L || derived$locator != "artifact:liftover_grch37_fasta_gz") {
    stop("liftover derived-FASTA registry entry is missing or has an unregistered source", call. = FALSE)
  }
  temporary <- paste0(source_fasta, ".partial-", Sys.getpid())
  unlink(temporary, force = TRUE)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  status <- system2(gzip, c("--decompress", "--keep", "--stdout", shQuote(source_gz)), stdout = temporary)
  identity <- as.character(derived$supplier_identity[[1L]])
  expected_bytes <- if (!is.na(identity) && nzchar(identity) && grepl("(^|;)bytes=", identity)) {
    as.numeric(sub(".*(^|;)bytes=([0-9]+).*", "\\2", identity))
  } else NA_real_
  valid <- file.exists(temporary) && file.info(temporary)$size > 0L &&
    (is.na(expected_bytes) || file.info(temporary)$size == expected_bytes)
  if ((status != 0L && status != 2L) || !valid) {
    stop("could not derive the verified GRCh37 liftover FASTA", call. = FALSE)
  }
  if (file.exists(source_fasta) && unlink(source_fasta, force = TRUE) != 0L) {
    stop("could not replace the decompressed GRCh37 liftover FASTA", call. = FALSE)
  }
  if (!file.rename(temporary, source_fasta)) stop("could not publish the decompressed GRCh37 FASTA", call. = FALSE)
  destination_fasta <- if (use_registry_paths) {
    duckhts_bench_fetch("liftover_grch38_fasta")
  } else {
    duckhts_bench_fetch("liftover_grch38_fasta", output = file.path(base_dir, "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"))
  }
  chain <- if (use_registry_paths) {
    duckhts_bench_fetch("liftover_grch37_grch38_chain")
  } else {
    duckhts_bench_fetch("liftover_grch37_grch38_chain", output = file.path(base_dir, "GRCh37_to_GRCh38.chain.gz"))
  }
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
