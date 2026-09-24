#' Stage the Oxford Nanopore E. coli K-12 MG1655 Alignment Benchmark Input
#'
#' Downloads ENA run `ERR14686255` (25,950 MinION reads of an E. coli K-12
#' MG1655 derivative, study PRJEB86481) and the NCBI RefSeq assembly FASTA
#' `GCF_000005845.2_ASM584v2_genomic.fna.gz`, verifies each against the
#' identity the registry pins (ENA's published MD5 and byte size for the reads,
#' NCBI's published MD5 for the archive, SHA-256 and byte size for the
#' uncompressed FASTA), and derives the coordinate-sorted, indexed BAM that
#' `benchmark_cigar_aligned_blocks.Rmd` reads by aligning with
#' `minimap2 -x map-ont` and sorting with `samtools`. Network access occurs only
#' in this explicit staging step; with `fetch = FALSE` both sources must already
#' be cached, which is how the report renders. A cached BAM that passes
#' `samtools quickcheck` is reused when its receipt's reference and read SHA-256
#' identities match the verified sources. Missing identities or changed inputs
#' require derivation; tool-version changes alone do not. `minimap2` is needed
#' only to derive the BAM.
#' @param fetch Whether to download a missing or invalid source.
#' @param threads Aligner, sort and index threads.
#' @param minimap2 Path to `minimap2`.
#' @param samtools Path to `samtools`. By default, use the executable on `PATH`,
#'   falling back to the optional RBCFTools package's bundled executable.
#' @return Named cache paths `reference` (FASTA), `reads` (FASTQ archive) and
#'   `bam` (the input).
#' @export
duckhts_bench_stage_ont_ecoli <- function(fetch = TRUE, threads = 8L,
                                          minimap2 = Sys.which("minimap2"),
                                          samtools = duckhts_bench_samtools()) {
  ids <- c("ont_ecoli_k12_reference_fna_gz", "ont_ecoli_k12_reference_fna",
           "ont_ecoli_k12_reads_fastq_gz", "ont_ecoli_k12_bam")
  plan <- duckhts_bench_stage_plan("ont-ecoli-k12")
  if (!identical(plan$id, ids)) stop("ont-ecoli-k12 registry plan is incomplete", call. = FALSE)
  transforms <- c("direct_download", "gunzip", "direct_download",
                  "minimap2_map_ont;samtools_sort;samtools_index")
  if (!identical(plan$transform, transforms) ||
      plan$locator[[2L]] != paste0("artifact:", ids[[1L]]) ||
      plan$locator[[4L]] != paste0("artifact:", ids[[2L]], ";artifact:", ids[[3L]])) {
    stop("ont-ecoli-k12 registry rows must be two downloads, a gunzip and a minimap2 derivation", call. = FALSE)
  }
  if (!nzchar(samtools)) stop("samtools is required to stage the ont-ecoli-k12 benchmark", call. = FALSE)
  paths <- stats::setNames(vapply(ids, duckhts_bench_artifact_path, character(1L)),
                           c("reference_gz", "reference", "reads", "bam"))

  for (source in c("reference_gz", "reads")) {
    id <- ids[[match(source, names(paths))]]
    if (fetch) {
      duckhts_bench_fetch(id)
    } else {
      if (!file.exists(paths[[source]])) {
        stop("ont-ecoli-k12 source is not staged; run duckhts_bench_stage_ont_ecoli(): ",
             paths[[source]], call. = FALSE)
      }
      duckhts_bench_validate_identity(id, paths[[source]])
    }
  }
  duckhts_bench_stage_gunzip(ids[[2L]], paths[["reference_gz"]], paths[["reference"]])

  source_hashes <- vapply(paths[c("reference", "reads")], function(path) {
    digest::digest(file = path, algo = "sha256")
  }, character(1L))
  names(source_hashes) <- paste0(names(source_hashes), "_sha256")
  bam <- paths[["bam"]]
  index <- paste0(bam, ".bai")
  receipt <- paste0(bam, ".provenance.tsv")
  run <- function(command, args, error) {
    status <- system2(command, args)
    if (status != 0L) stop(error, call. = FALSE)
  }
  if (file.exists(bam) && file.exists(index) && file.exists(receipt) &&
      system2(samtools, c("quickcheck", shQuote(bam))) == 0L) {
    fields <- utils::read.delim(receipt, colClasses = "character", quote = "", comment.char = "")
    recorded <- fields$value[match(names(source_hashes), fields$field)]
    if (identical(recorded, unname(source_hashes))) {
      return(invisible(paths[c("reference", "reads", "bam")]))
    }
  }
  if (!nzchar(minimap2)) stop("minimap2 is required to derive the ont-ecoli-k12 BAM", call. = FALSE)

  dir.create(dirname(bam), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(bam, ".partial-", Sys.getpid())
  alignment <- paste0(partial, ".sam")
  unlink(c(partial, paste0(partial, ".bai"), alignment), force = TRUE)
  on.exit(unlink(c(partial, paste0(partial, ".bai"), alignment), force = TRUE), add = TRUE)
  run(minimap2, c("-a", "-x", "map-ont", "-t", threads, "-o", shQuote(alignment),
                  shQuote(paths[["reference"]]), shQuote(paths[["reads"]])),
      "could not align the ont-ecoli-k12 reads")
  run(samtools, c("sort", "-@", threads, "-o", shQuote(partial), shQuote(alignment)),
      "could not sort the ont-ecoli-k12 alignment")
  run(samtools, c("quickcheck", "-v", shQuote(partial)), "derived ont-ecoli-k12 BAM failed samtools quickcheck")
  run(samtools, c("index", "-@", threads, shQuote(partial)), "could not index the ont-ecoli-k12 BAM")
  unlink(c(bam, index, receipt), force = TRUE)
  if (!file.rename(partial, bam) || !file.rename(paste0(partial, ".bai"), index)) {
    stop("could not publish the ont-ecoli-k12 BAM: ", bam, call. = FALSE)
  }
  tool_version <- function(command) {
    paste(system2(command, "--version", stdout = TRUE, stderr = TRUE)[[1L]], collapse = " ")
  }
  fields <- rbind(
    duckhts_bench_provenance_fields(ids[[4L]], bam),
    data.frame(
      field = c("run_accession", "study_accession", "sample_accession", "reference", "reads",
                "aligner", "aligner_preset", "sorter", "threads", names(source_hashes)),
      value = c("ERR14686255", "PRJEB86481", "SAMEA117787661", paths[["reference"]], paths[["reads"]],
                paste("minimap2", tool_version(minimap2)), "map-ont", tool_version(samtools),
                as.character(threads), unname(source_hashes)),
      stringsAsFactors = FALSE
    )
  )
  utils::write.table(fields, receipt, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(paths[c("reference", "reads", "bam")])
}

# Resolve the executable once for staging and its dependency checks.
duckhts_bench_samtools <- function() {
  samtools <- Sys.which("samtools")
  if (!nzchar(samtools) && requireNamespace("RBCFTools", quietly = TRUE) &&
      package_version(getNamespaceVersion("RBCFTools")) >= "1.24-1.1.0") {
    samtools <- RBCFTools::samtools_path()
  }
  unname(samtools)
}
