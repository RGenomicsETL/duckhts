#!/usr/bin/env Rscript
# Network-free staging contract: committed smoke fixtures and a real CRAM-to-BAM
# transcode under an isolated registry/cache. No public source is downloaded.
root <- normalizePath(commandArgs(trailingOnly = TRUE)[[1L]], mustWork = TRUE)
setwd(root)
pkgload::load_all("r/duckhtsbench", quiet = TRUE)
source("benchmarks/insert_size_stage.R")
local({
  scratch <- tempfile("insert-size-stage-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE))
  old_cache <- Sys.getenv("DUCKHTS_CACHE_DIR", unset = NA_character_)
  Sys.setenv(DUCKHTS_CACHE_DIR = scratch)
  on.exit(if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else
    Sys.setenv(DUCKHTS_CACHE_DIR = old_cache), add = TRUE)
  paths <- duckhtsbench::duckhts_bench_stage_repository_fixtures(root, "insert-size-smoke")
  stopifnot(length(paths) == 8L, all(file.exists(paths)))
  extension <- normalizePath("build/release/duckhts.duckdb_extension")
  original <- getFromNamespace("duckhts_bench_registry", "duckhtsbench")
  registry <- original()
  target <- registry$id == "insert_size_hg00403_chr20_bam"
  registry$locator[target] <- "artifact:insert_size_smoke_edges_cram;artifact:insert_size_smoke_reference_fa"
  identity <- bam_mismatch_slice_identity(paths[["insert_size_smoke_edges_cram"]],
                                         paths[["insert_size_smoke_reference_fa"]], extension)
  registry$supplier_identity[target] <- paste(paste(names(identity), unlist(identity), sep = "="), collapse = ";")
  for (format in c("bam", "cram")) {
    for (scale in c(2L, 4L)) {
      id <- paste0("insert_size_hg00403_chr20_", format, "_", scale, "x")
      selected <- registry$id == id
      registry$locator[selected] <- paste0("artifact:insert_size_smoke_edges_", format,
                                          ";artifact:insert_size_smoke_reference_fa")
      scaled <- lapply(identity, function(value) value * scale)
      registry$supplier_identity[selected] <- paste0("copies=", scale, ";",
        paste(paste(names(scaled), unlist(scaled), sep = "="), collapse = ";"))
    }
  }
  assignInNamespace("duckhts_bench_registry", function() registry, "duckhtsbench")
  on.exit(assignInNamespace("duckhts_bench_registry", original, "duckhtsbench"), add = TRUE)
  output <- stage_insert_size_bam_from_registry(extension)
  stopifnot(file.exists(output), file.exists(paste0(output, ".bai")),
            file.exists(paste0(output, ".provenance.tsv")))
  stopifnot(identical(stage_insert_size_bam_from_registry(extension), output))
  observed <- bam_mismatch_slice_identity(output, paths[["insert_size_smoke_reference_fa"]], extension)
  stopifnot(identical(observed, identity))
  for (format in c("bam", "cram")) {
    for (scale in c(2L, 4L)) {
      id <- paste0("insert_size_hg00403_chr20_", format, "_", scale, "x")
      path <- stage_insert_size_scale_from_registry(id, extension)
      observed <- bam_mismatch_slice_identity(path, paths[["insert_size_smoke_reference_fa"]], extension)
      stopifnot(identical(observed, lapply(identity, function(value) value * scale)),
                identical(stage_insert_size_scale_from_registry(id, extension), path))
      # An unchanged-size damaged input must not be reused under its receipt.
      bytes <- readBin(path, "raw", n = file.size(path))
      damaged <- bytes
      damaged[[1L]] <- as.raw(bitwXor(as.integer(damaged[[1L]]), 1L))
      writeBin(damaged, path)
      failed <- tryCatch({
        stage_insert_size_scale_from_registry(id, extension)
        FALSE
      }, error = function(error) grepl("mismatched receipt", conditionMessage(error)))
      stopifnot(failed)
      writeBin(bytes, path)
    }
  }
  cat("PASS: 8 fixtures; BAM transcode; BAM/CRAM 2x/4x inputs, identities, receipts, reuse and corruption\n")
})
