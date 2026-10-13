# Network-free BAM derivation from the registered, already staged chr20 CRAM.
# Shared executable selection and encoding-independent identity checks belong
# to bam_mismatch_stage.R. This stage never downloads its raw inputs.
source(file.path("benchmarks", "bam_mismatch_stage.R"))

stage_insert_size_bam_from_registry <- function(extension,
  output = duckhtsbench::duckhts_bench_artifact_path("insert_size_hg00403_chr20_bam")
) {
  id <- "insert_size_hg00403_chr20_bam"
  registry <- duckhtsbench::duckhts_bench_registry()
  row <- registry[registry$id == id, ]
  stopifnot(nrow(row) == 1L, row$transform == "stage_insert_size_bam_from_registry")
  inputs <- sub("^artifact:", "", strsplit(row$locator, ";", fixed = TRUE)[[1L]])
  stopifnot(length(inputs) == 2L)
  cram <- duckhtsbench::duckhts_bench_artifact_path(inputs[[1L]])
  reference <- duckhtsbench::duckhts_bench_artifact_path(inputs[[2L]])
  stopifnot(file.exists(cram), file.exists(reference))
  for (i in seq_along(inputs)) duckhtsbench::duckhts_bench_validate_identity(
    inputs[[i]], c(cram, reference)[[i]])
  expected <- getFromNamespace("duckhts_bench_identity_fields", "duckhtsbench")(row$supplier_identity)
  validate <- function(path) {
    observed <- bam_mismatch_slice_identity(path, reference, extension)
    if (!bam_mismatch_same_identity(observed, expected)) {
      stop("derived insert-size BAM alignment identity mismatch", call. = FALSE)
    }
  }
  if (file.exists(output)) {
    validate(output)
    stopifnot(file.exists(paste0(output, ".bai")),
              file.exists(paste0(output, ".provenance.tsv")))
    return(invisible(output))
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  available <- system2("df", c("--output=avail", "-B1", shQuote(dirname(output))), stdout = TRUE)
  free_bytes <- as.numeric(trimws(tail(available, 1L)))
  if (!is.finite(free_bytes) || free_bytes < (8 + 6) * 1024^3) {
    stop("staging requires 8 GiB free plus a 6 GiB transcode reserve", call. = FALSE)
  }
  temporary <- paste0(output, ".partial-", Sys.getpid(), ".bam")
  on.exit(unlink(c(temporary, paste0(temporary, ".bai"))), add = TRUE)
  bam_mismatch_run("samtools", c("view", "--no-PG", "-b", "-T", reference,
                                 "-o", temporary, cram))
  validate(temporary)
  bam_mismatch_run("samtools", c("index", temporary))
  stopifnot(file.rename(temporary, output), file.rename(paste0(temporary, ".bai"), paste0(output, ".bai")))
  duckhtsbench::duckhts_bench_write_provenance(id, output)
  invisible(output)
}
