#!/usr/bin/env Rscript
# Network-free acquisition tests use real GNU tar and the canonical receipt validator.
repo <- normalizePath(".")
stopifnot(file.exists(file.path(repo, "scripts/duckvep_evidence.R")))
for (file in c("registry.R", "stage.R", "vep_cache.R")) {
  source(file.path(repo, "r/duckhtsbench/R", file))
}
Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(repo, "r/duckhtsbench/inst/benchmark_registry.tsv"))
plan <- duckhts_bench_stage_plan("duckvep-vep-cache")
stopifnot(nrow(plan) == 2L)
local({
  temporary <- tempfile("vep-cache-stage-")
  dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  source_root <- file.path(temporary, "archive source")
  prefix <- "homo_sapiens/116_GRCh38"
  members <- c("info.txt", "chr_synonyms.txt", ".metadata", "21/1-1000000.gz",
    "21/1000001-2000000.gz", "21/nested/index", "1/1-1000000.gz")
  for (member in members) {
    path <- file.path(source_root, prefix, member)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(paste("contents", member), path)
  }
  archive <- file.path(temporary, "test cache.tar.gz")
  stopifnot(system2("tar", shQuote(c("-czf", archive, "-C", source_root, "homo_sapiens"))) == 0L)
  registry <- file.path(temporary, "registry.tsv")
  plan$supplier_identity[[1L]] <- paste0("md5=", unname(tools::md5sum(archive)),
    ";bytes=", file.info(archive)$size)
  save_registry <- function() utils::write.table(plan, registry,
    sep = "\t", row.names = FALSE, quote = FALSE)
  save_registry()
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry, DUCKHTS_CACHE_DIR = file.path(temporary, "local cache"))
  paths <- duckhts_bench_stage_vep_cache(repo, archive = archive)
  kept <- members[!startsWith(members, "1/")]
  actual <- list.files(file.path(paths[["cache"]], prefix), recursive = TRUE, all.files = TRUE)
  stopifnot(identical(sort(actual), sort(kept)), file.exists(paths[["receipt"]]),
    identical(unname(tools::md5sum(file.path(paths[["cache"]], prefix, kept))),
      unname(tools::md5sum(file.path(source_root, prefix, kept)))),
    identical(duckhts_bench_stage_vep_cache(repo, archive = archive), paths))
  expect_error <- function(expr, pattern) {
    message <- tryCatch({ force(expr); "" }, error = conditionMessage)
    stopifnot(nzchar(message), grepl(pattern, message, fixed = TRUE))
  }
  writeLines("mutated shard", file.path(paths[["cache"]], prefix, "21/1-1000000.gz"))
  expect_error(duckhts_bench_stage_vep_cache(repo, archive = archive), "differs from its acquisition receipt")

  # Exercise pipefail even when tar receives a complete valid archive before curl fails.
  fake_curl <- file.path(temporary, "fake curl")
  writeLines(c("#!/usr/bin/env bash", "set -eu", "headers=''", "conditional=''",
    "while (($#)); do",
    "  case \"$1\" in",
    "    --dump-header) headers=$2; shift 2;;",
    "    --header) conditional=$2; shift 2;;",
    "    *) shift;;",
    "  esac",
    "done",
    "test \"$conditional\" = 'If-Match: \"abcd-1234\"'",
    "printf 'HTTP/1.1 200 OK\\r\\nETag: \"%s\"\\r\\n\\r\\n' \"${VEP_TEST_ETAG:-abcd-1234}\" > \"$headers\"",
    paste("cat", shQuote(archive)),
    paste0("printf '\\nduckhts_cache_bytes=%s\\n' \"${VEP_TEST_BYTES:-", file.info(archive)$size, "}\" >&2"),
    "exit \"${VEP_TEST_EXIT:-0}\""), fake_curl)
  Sys.chmod(fake_curl, "0755")
  plan$supplier_identity[[1L]] <- paste0("etag=abcd-1234;bytes=", file.info(archive)$size)
  save_registry()
  Sys.setenv(DUCKHTS_CACHE_DIR = file.path(temporary, "stream cache"))
  streamed <- duckhts_bench_stage_vep_cache(repo, curl = fake_curl)
  stopifnot(identical(sort(list.files(file.path(streamed[["cache"]], prefix), recursive = TRUE, all.files = TRUE)),
    sort(kept)))
  for (failure in c("exit", "etag", "bytes")) {
    Sys.setenv(DUCKHTS_CACHE_DIR = file.path(temporary, paste0("failed-", failure)))
    if (failure == "exit") Sys.setenv(VEP_TEST_EXIT = "22")
    if (failure == "etag") Sys.setenv(VEP_TEST_ETAG = "ffff-ffff")
    if (failure == "bytes") Sys.setenv(VEP_TEST_BYTES = "1")
    expect_error(duckhts_bench_stage_vep_cache(repo, curl = fake_curl),
      if (failure == "exit") "acquisition/extraction failed" else "ETag/byte identity")
    stopifnot(!file.exists(duckhts_bench_artifact_path("vep116_grch38_cache_chr21")))
    Sys.unsetenv(c("VEP_TEST_EXIT", "VEP_TEST_ETAG", "VEP_TEST_BYTES"))
  }

  plan$supplier_identity[[1L]] <- paste0(plan$supplier_identity[[1L]],
    ";md5=", unname(tools::md5sum(archive)))
  save_registry()
  Sys.setenv(DUCKHTS_CACHE_DIR = file.path(temporary, "unverified checksum"))
  expect_error(duckhts_bench_stage_vep_cache(repo, curl = fake_curl),
    "supply a checksum-pinned archive locally")

  # A valid checksum does not excuse an invalid/truncated gzip archive.
  truncated <- file.path(temporary, "truncated.tar.gz")
  input <- file(archive, "rb")
  bytes <- readBin(input, "raw", n = floor(file.info(archive)$size / 2))
  close(input)
  writeBin(bytes, truncated)
  plan$supplier_identity[[1L]] <- paste0("md5=", unname(tools::md5sum(truncated)),
    ";bytes=", length(bytes))
  save_registry()
  Sys.setenv(DUCKHTS_CACHE_DIR = file.path(temporary, "truncated cache"))
  expect_error(duckhts_bench_stage_vep_cache(repo, archive = truncated), "acquisition/extraction failed")
  stopifnot(!file.exists(duckhts_bench_artifact_path("vep116_grch38_cache_chr21")))

  # Refuse to replace an unreceipted destination or accept a changed selection.
  output <- duckhts_bench_artifact_path("vep116_grch38_cache_chr21")
  dir.create(output, recursive = TRUE)
  writeLines("preserve", file.path(output, "user-data"))
  expect_error(duckhts_bench_stage_vep_cache(repo, archive = archive), "acquisition.tsv")
  stopifnot(identical(readLines(file.path(output, "user-data")), "preserve"))
  Sys.setenv(DUCKHTS_CACHE_DIR = file.path(temporary, "stream cache"))
  plan$supplier_identity[[1L]] <- paste0("etag=abcd-1234;bytes=", file.info(archive)$size)
  plan$supplier_identity[[2L]] <- sub("regions=21", "regions=22", plan$supplier_identity[[2L]])
  save_registry()
  expect_error(duckhts_bench_stage_vep_cache(repo, curl = fake_curl), "registered source or scope")

  # Run the actual CLI with the same local, checksummed fixture source.
  Sys.setenv(DUCKHTS_CACHE_DIR = file.path(temporary, "CLI cache"))
  plan$supplier_identity[[1L]] <- paste0("md5=", unname(tools::md5sum(archive)),
    ";bytes=", file.info(archive)$size)
  plan$supplier_identity[[2L]] <- sub("regions=22", "regions=21", plan$supplier_identity[[2L]])
  save_registry()
  stopifnot(system2(file.path(R.home("bin"), "Rscript"), shQuote(c(
    file.path(repo, "r/duckhtsbench/scripts/stage_vep_cache.R"), "--archive", archive
  )), stdout = file.path(temporary, "cli.log"), stderr = file.path(temporary, "cli.err")) == 0L)
  output <- duckhts_bench_artifact_path("vep116_grch38_cache_chr21")
  stopifnot(identical(sort(list.files(file.path(output, prefix), recursive = TRUE, all.files = TRUE)),
    sort(kept)))
  cat("VEP regional cache staging: exact complete region/metadata, reuse, mutation, transfer failures, truncation and preserved destinations: OK\n")
})
