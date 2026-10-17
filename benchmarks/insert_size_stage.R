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

# Coordinate-merge repeated copies into one larger file. samtools -c/-p keeps
# RG/PG identifiers unchanged, so scaling does not alter insert-size semantics.
stage_insert_size_scale_from_registry <- function(id, extension) {
  row <- bam_mismatch_registry_row(id)
  registry <- row$registry
  target <- registry[registry$id == id, , drop = FALSE]
  stopifnot(target$transform == "stage_insert_size_scale_from_registry",
            length(row$inputs) == 2L)
  bam_mismatch_require_fields(id, row$fields,
    c("copies", "alignments", "position_sum", "base_sum"))
  copies <- as.numeric(row$fields[["copies"]])
  if (!is.finite(copies) || !copies %in% c(2, 4)) {
    stop("insert-size scale must contain 2 or 4 copies", call. = FALSE)
  }
  source <- duckhtsbench::duckhts_bench_artifact_path(row$inputs[[1L]])
  reference <- duckhtsbench::duckhts_bench_artifact_path(row$inputs[[2L]])
  stopifnot(file.exists(source), file.exists(reference))
  duckhtsbench::duckhts_bench_validate_identity(row$inputs[[1L]], source)
  duckhtsbench::duckhts_bench_validate_identity(row$inputs[[2L]], reference)
  output <- duckhtsbench::duckhts_bench_artifact_path(id)
  suffix <- if (endsWith(output, ".bam")) ".bai" else if (endsWith(output, ".cram")) ".crai" else
    stop("scale output must be BAM or CRAM", call. = FALSE)
  sha256 <- function(path) unname(digest::digest(file = path, algo = "sha256"))
  receipt <- paste0(output, ".sources.tsv")
  expected <- data.frame(
    field = c("artifact_id", "supplier_identity", "source_sha256", "reference_sha256",
              "samtools_version", "derivation"),
    value = c(id, target$supplier_identity, sha256(source), sha256(reference),
      system2(bam_mismatch_tool("samtools"), "--version", stdout = TRUE,
              env = bam_mismatch_tool_environment())[[1L]],
      "samtools merge -c -p --no-PG; samtools index"))
  validate <- function(path) {
    observed <- bam_mismatch_slice_identity(path, reference, extension)
    if (!bam_mismatch_same_identity(observed, row$fields)) {
      stop("scaled insert-size input differs from its registered identity", call. = FALSE)
    }
  }
  if (file.exists(output)) {
    stored <- if (file.exists(receipt)) utils::read.delim(receipt, colClasses = "character") else NULL
    current <- !is.null(stored) && nrow(stored) == nrow(expected) + 2L &&
      identical(stored[seq_len(nrow(expected)), ], expected) &&
      file.exists(paste0(output, suffix)) &&
      identical(stored$value[stored$field == "output_sha256"], sha256(output)) &&
      identical(stored$value[stored$field == "index_sha256"], sha256(paste0(output, suffix)))
    if (!current) stop("cached scale input has a missing or mismatched receipt: ", id, call. = FALSE)
    validate(output)
    return(invisible(output))
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  available <- system2("df", c("--output=avail", "-B1", shQuote(dirname(output))), stdout = TRUE)
  free_bytes <- as.numeric(trimws(tail(available, 1L)))
  reserve <- 1.2 * copies * file.size(source)
  if (!is.finite(free_bytes) || free_bytes < 8 * 1024^3 + reserve) {
    stop("scale staging requires 8 GiB free plus 1.2 times the expanded input bytes", call. = FALSE)
  }
  temporary <- paste0(output, ".partial-", Sys.getpid(), if (suffix == ".bai") ".bam" else ".cram")
  temporary_receipt <- paste0(receipt, ".partial-", Sys.getpid())
  on.exit(unlink(c(temporary, paste0(temporary, suffix), temporary_receipt)), add = TRUE)
  bam_mismatch_run("samtools", c("merge", "-c", "-p", "--no-PG", "--reference", reference,
                                 "-o", temporary, rep(source, copies)))
  bam_mismatch_run("samtools", c("index", temporary))
  validate(temporary)
  hashes <- data.frame(field = c("output_sha256", "index_sha256"),
    value = c(sha256(temporary), sha256(paste0(temporary, suffix))))
  utils::write.table(rbind(expected, hashes), temporary_receipt,
                     sep = "\t", row.names = FALSE, quote = FALSE)
  duckhtsbench:::duckhts_bench_publish_files(
    c(temporary, paste0(temporary, suffix), temporary_receipt),
    c(output, paste0(output, suffix), receipt))
  duckhtsbench::duckhts_bench_write_provenance(id, output)
  invisible(output)
}
