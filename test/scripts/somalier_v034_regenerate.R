#!/usr/bin/env Rscript

# Rebuild the retained helper output directly from the pinned Somalier source.
somalier_v034_run_witness <- function(command, args = character()) {
  observed <- suppressWarnings(system2(command, args, stdout = TRUE))
  status <- attr(observed, "status")
  if (!is.null(status) && status != 0L) {
    stop("Pinned Somalier helper exited ", status)
  }
  observed
}

main <- function(args) {
  if (length(args) > 2L) {
    stop("usage: Rscript somalier_v034_regenerate.R [source-archive] [output.tsv]")
  }
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
  root <- normalizePath(file.path(dirname(script), "../.."))
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required")

  archive <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else {
    if (!requireNamespace("duckhtsbench", quietly = TRUE)) {
      stop("duckhtsbench is required when source-archive is omitted")
    }
    duckhtsbench::duckhts_bench_fetch("somalier_v034_source_archive")
  }
  expected_archive <- "acd2dc11be6051c80d15628703a9419965cea3fe7a0563b47e14743e7ac6339e"
  if (!identical(digest::digest(archive, file = TRUE, algo = "sha256"), expected_archive) ||
      file.info(archive)$size != 1208978) {
    stop("Somalier v0.3.4 source archive identity mismatch")
  }

  driver <- file.path(root, "test/scripts/somalier_v034_witness.nim")
  expected_driver <- "0f3bbd0c7088e98859221172c73d1399eefb2db33a7afc07052277726d7f0f32"
  if (!identical(digest::digest(driver, file = TRUE, algo = "sha256"), expected_driver)) {
    stop("Somalier witness driver identity mismatch")
  }

  directory <- tempfile("somalier-v034-source-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  prefix <- "somalier-ff58fdade8f4f8293d904f10e0a4a13f1fac808d/src/somalierpkg/"
  utils::untar(archive, files = paste0(prefix, c("charr.nim", "contamination.nim")),
    exdir = directory)
  source_dir <- file.path(directory,
    "somalier-ff58fdade8f4f8293d904f10e0a4a13f1fac808d/src/somalierpkg")
  charr_path <- file.path(source_dir, "charr.nim")
  pair_path <- file.path(source_dir, "contamination.nim")
  blob <- function(path) system2("git", c("hash-object", path), stdout = TRUE)
  if (!identical(blob(charr_path), "fb8563c5a5e1c4fa9c9edcf0956bf679a53b3bd6") ||
      !identical(blob(pair_path), "c1ef08d887b08889c9346bea770c3e20f5d518a4")) {
    stop("Pinned Somalier source blob identity mismatch")
  }

  charr <- readLines(charr_path, warn = FALSE)
  pair <- readLines(pair_path, warn = FALSE)
  if (length(charr) < 100L || length(pair) < 513L) stop("Pinned source is truncated")
  writeLines(charr[5:100], file.path(directory, "upstream_charr_helpers.nim"))
  ranges <- list(25:27, 51:61, 68:87, 268:288, 319:321, 336:433, 452:513)
  pair_helpers <- unlist(lapply(ranges, function(range) c(pair[range], "")), use.names = FALSE)
  writeLines(pair_helpers[-length(pair_helpers)],
    file.path(directory, "upstream_pair_helpers.nim"))
  if (!file.copy(driver, file.path(directory, "witness.nim"), overwrite = TRUE)) {
    stop("Could not copy the checked-in Somalier witness driver")
  }

  nim <- Sys.which("nim")
  if (!nzchar(nim)) stop("Nim is required to regenerate the pinned helper output")
  old <- setwd(directory)
  on.exit(setwd(old), add = TRUE)
  status <- system2(nim, c("c", "-d:release", "--out:witness", "witness.nim"),
    stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(status, "status")) && attr(status, "status") != 0L) {
    stop("Pinned Somalier helper compilation failed:\n", paste(status, collapse = "\n"))
  }
  observed <- somalier_v034_run_witness(file.path(directory, "witness"))
  bytes <- charToRaw(paste0(paste(observed, collapse = "\n"), "\n"))
  expected_output <- "34e8dff0151c4d556aaaf3c4b3c03bfabdc55e00861f60b058d74303ea3c4aeb"
  if (!identical(digest::digest(bytes, algo = "sha256", serialize = FALSE), expected_output)) {
    stop("Regenerated helper output differs from the retained upstream receipt")
  }
  if (length(args) == 2L) writeBin(bytes, args[[2L]])
  cat("Somalier v0.3.4 source-derived helper: OK (14 rows; SHA-256 ",
    expected_output, ")\n", sep = "")
}

if (sys.nframe() == 0L) main(commandArgs(trailingOnly = TRUE))
