#!/usr/bin/env Rscript

# Exercise the Somalier source-artifact cache and witness process contract
# without a network request or a Nim installation.
main <- function() {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
  root <- normalizePath(file.path(dirname(script), "../.."), winslash = "/")
  source(file.path(root, "r/duckhtsbench/R/registry.R"), local = FALSE)
  source(file.path(root, "r/duckhtsbench/R/stage.R"), local = FALSE)
  source(file.path(root, "test/scripts/somalier_v034_regenerate.R"), local = FALSE)

  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"),
    unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, as.list(previous[name]))
    }
  }, add = TRUE)
  directory <- tempfile("somalier-upstream-staging-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  registry <- utils::read.delim(
    file.path(root, "r/duckhtsbench/inst/benchmark_registry.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  row <- registry[registry$id == "somalier_v034_source_archive", , drop = FALSE]
  stopifnot(nrow(row) == 1L, row$transform == "direct_download")
  payload <- charToRaw("network-free-somalier-staging-v1\n")
  row$locator <- "https://invalid.invalid/somalier-source-fixture.tar.gz"
  row$cache_relpath <- "upstream/somalier/test/source-fixture.tar.gz"
  row$supplier_identity <- paste0(
    "sha256=027a320cd32bfc15d9d38825e03678de601db3aa0c35c1eeb6ca714a174b09a2;",
    "bytes=33;commit=test-fixture;license=MIT")
  registry_path <- file.path(directory, "registry.tsv")
  utils::write.table(row, registry_path, sep = "\t", row.names = FALSE,
    quote = FALSE)
  cache <- file.path(directory, "cache")
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path, DUCKHTS_CACHE_DIR = cache)
  destination <- duckhts_bench_artifact_path("somalier_v034_source_archive")
  dir.create(dirname(destination), recursive = TRUE)
  writeBin(payload, destination)
  staged <- duckhts_bench_fetch("somalier_v034_source_archive")
  stopifnot(identical(normalizePath(staged), normalizePath(destination)),
    identical(readBin(staged, "raw", n = length(payload)), payload),
    file.exists(paste0(staged, ".provenance.tsv")))

  good <- file.path(directory, "witness-good.R")
  bad <- file.path(directory, "witness-bad.R")
  writeLines("cat('metric\\tinput\\tvalue\\n')", good)
  writeLines(c("cat('metric\\tinput\\tvalue\\n')",
    "quit(save = 'no', status = 7L)"), bad)
  rscript <- file.path(R.home("bin"), "Rscript")
  observed <- somalier_v034_run_witness(rscript,
    c("--vanilla", shQuote(good)))
  stopifnot(identical(observed, "metric\tinput\tvalue"))
  error <- tryCatch({
    somalier_v034_run_witness(rscript, c("--vanilla", shQuote(bad)))
    NULL
  }, error = function(condition) conditionMessage(condition))
  stopifnot(identical(error, "Pinned Somalier helper exited 7"))

  cat("Somalier upstream network-free staging and witness exit: OK\n")
}

main()
