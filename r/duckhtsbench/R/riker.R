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
  ref_index <- paste0(ref, ".fai")
  if (file.exists(ref_index) && unlink(ref_index, force = TRUE) != 0L) stop("could not replace Riker reference index", call. = FALSE)
  run(c("faidx", ref), "could not index Riker reference")
  if (!file.exists(ref_index) || !file.info(ref_index)$size) stop("could not publish Riker reference index", call. = FALSE)
  run(c("quickcheck", "-v", cram), "staged Riker CRAM failed samtools quickcheck")
  if (file.exists(bam) && unlink(bam, force = TRUE) != 0L) stop("could not replace Riker BAM", call. = FALSE)
  run(c("view", "-@", threads, "-b", "-T", ref, "-o", bam, cram), "could not transcode Riker CRAM")
  if (!file.exists(bam) || !file.info(bam)$size) stop("could not publish Riker BAM", call. = FALSE)
  run(c("quickcheck", "-v", bam), "derived Riker BAM failed samtools quickcheck")
  bam_index <- paste0(bam, ".bai")
  if (file.exists(bam_index) && unlink(bam_index, force = TRUE) != 0L) stop("could not replace Riker BAM index", call. = FALSE)
  run(c("index", "-@", threads, bam), "could not index Riker BAM")
  if (!file.exists(bam_index) || !file.info(bam_index)$size) stop("could not publish Riker BAM index", call. = FALSE)
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
