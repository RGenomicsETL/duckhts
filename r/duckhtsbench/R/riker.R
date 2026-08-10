#' Stage the Riker-versus-DuckHTS WGS benchmark inputs.
#'
#' @param base_dir Destination directory, defaulting to the registered cache path.
#' @param threads Samtools decompression threads.
#' @param samtools Path to `samtools`.
#' @return Path to the derived BAM, invisibly.
#' @export
duckhts_bench_stage_riker <- function(base_dir = duckhts_bench_cache_path("benchmarks/riker-wgs"),
                                      threads = 8L, samtools = Sys.which("samtools")) {
  if (!nzchar(samtools)) stop("samtools is required to stage the Riker WGS benchmark", call. = FALSE)
  ref <- duckhts_bench_fetch("riker_hg00188_reference", output = file.path(base_dir, "reference", "GRCh38_full_analysis_set_plus_decoy_hla.fa"))
  cram <- duckhts_bench_fetch("riker_hg00188_cram", output = file.path(base_dir, "downloads", "HG00188.final.cram"))
  bam <- file.path(base_dir, "stage", "HG00188_30x", "input.bam")
  dir.create(dirname(bam), recursive = TRUE, showWarnings = FALSE)
  run <- function(args, error) {
    status <- system2(samtools, args)
    if (status != 0L) stop(error, call. = FALSE)
  }
  if (!file.exists(paste0(ref, ".fai")) || !file.info(paste0(ref, ".fai"))$size) run(c("faidx", ref), "could not index Riker reference")
  run(c("quickcheck", "-v", cram), "staged Riker CRAM failed samtools quickcheck")
  if (!file.exists(bam) || !file.info(bam)$size) run(c("view", "-@", threads, "-b", "-T", ref, "-o", bam, cram), "could not transcode Riker CRAM")
  run(c("quickcheck", "-v", bam), "derived Riker BAM failed samtools quickcheck")
  if (!file.exists(paste0(bam, ".bai")) || !file.info(paste0(bam, ".bai"))$size) run(c("index", "-@", threads, bam), "could not index Riker BAM")
  receipt <- file.path(base_dir, "provenance.tsv")
  writeLines(c(
    "field\tvalue",
    "workload\triker-wgs",
    "sample\tHG00188",
    "run_accession\tERR3240174",
    paste0("reference\t", ref), paste0("cram\t", cram), paste0("bam\t", bam), paste0("threads\t", threads)
  ), receipt)
  invisible(bam)
}
