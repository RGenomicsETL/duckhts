library(tinytest)

plan <- duckhts_bench_stage_plan("tabix-single-split")
expect_equal(plan$id, c("tabix_split_gff3", "tabix_split_gtf", "tabix_split_bed"))
expect_equal(plan$transform[[3L]], "first_100000_gff3_records_to_bed16")
expect_equal(plan$locator[[3L]], "artifact:tabix_split_gff3")
expect_true(all(grepl("benchmark_tabix_split.Rmd", plan$consumer, fixed = TRUE)))

# Run the registered derivation using only synthetic inputs in a private cache.
test_tabix_split_staging <- function() {
  root <- Sys.getenv("DUCKHTS_REPO", unset = "")
  if (!nzchar(root)) return(invisible(NULL))
  script <- file.path(root, "scripts/stage_tabix_split.R")
  expect_true(file.exists(script))

  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  }, add = TRUE)
  work <- tempfile("tabix-split-stage-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  registry <- plan
  registry$cache_relpath <- c("raw/input.gff3.gz", "raw/input.gtf.gz", "derived/sample.bed")
  registry_path <- file.path(work, "registry.tsv")
  cache <- file.path(work, "cache")
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path, DUCKHTS_CACHE_DIR = cache)

  utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
  gff_line <- "chr1\tsrc\tgene\t1\t1\t.\t+\t.\tID=g"
  bed_line <- "chr1\t0\t1\tgene\t.\t+\tsrc\t.\tID=g\tchr1\tgene\t1\t1\t+\t.\tID=g"
  expected <- file.path(work, "expected.bed")
  writeLines(rep(bed_line, 100000L), expected, useBytes = TRUE)
  expected_bytes <- file.info(expected)$size
  expected_sha <- digest::digest(file = expected, algo = "sha256")

  raw_gff <- duckhts_bench_artifact_path("tabix_split_gff3")
  raw_gtf <- duckhts_bench_artifact_path("tabix_split_gtf")
  dir.create(dirname(raw_gff), recursive = TRUE)
  con <- gzfile(raw_gff, "wt")
  writeLines(c("##gff-version 3", rep(gff_line, 100000L)), con)
  close(con)
  con <- gzfile(raw_gtf, "wt")
  writeLines(rep("chr1\tsrc\tgene\t1\t1\t.\t+\t.\tgene_id \"g\";", 2L), con)
  close(con)
  registry$supplier_identity <- c(
    paste0("bytes=", file.info(raw_gff)$size, ";sha256=", digest::digest(file = raw_gff, algo = "sha256")),
    paste0("bytes=", file.info(raw_gtf)$size, ";sha256=", digest::digest(file = raw_gtf, algo = "sha256")),
    paste0("bytes=", expected_bytes, ";sha256=", expected_sha)
  )
  utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
  output <- duckhts_bench_artifact_path("tabix_split_bed")
  run_stage <- function() system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), "--offline"), stdout = TRUE, stderr = TRUE)
  stage <- function() {
    result <- run_stage()
    expect_equal(if (is.null(attr(result, "status"))) 0L else attr(result, "status"), 0L)
    expect_true(any(grepl("Staged 3 tabix split inputs", result, fixed = TRUE)))
  }

  registry$supplier_identity[[3L]] <- paste0("bytes=", expected_bytes, ";sha256=", strrep("0", 64L))
  utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
  failed <- suppressWarnings(run_stage())
  expect_true(!is.null(attr(failed, "status")) && attr(failed, "status") != 0L)
  expect_false(file.exists(output))
  expect_false(any(grepl("partial-", list.files(dirname(output)), fixed = TRUE)))
  registry$supplier_identity[[3L]] <- paste0("bytes=", expected_bytes, ";sha256=", expected_sha)
  utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)

  stage()
  expect_equal(file.info(output)$size, expected_bytes)
  expect_equal(digest::digest(file = output, algo = "sha256"), expected_sha)
  expect_identical(readLines(output, n = 1L), bed_line)
  expect_identical(tail(readLines(output), 1L), bed_line)
  for (id in plan$id) {
    path <- duckhts_bench_artifact_path(id)
    expect_true(duckhts_bench_validate_identity(id, path))
    receipt <- utils::read.delim(paste0(path, ".provenance.tsv"), colClasses = "character")
    expect_equal(receipt$value[receipt$field == "artifact_id"], id)
    expect_equal(receipt$value[receipt$field == "supplier_identity"],
      registry$supplier_identity[registry$id == id])
  }
  expect_false(any(grepl("partial-", list.files(dirname(output)), fixed = TRUE)))

  Sys.setFileTime(output, Sys.time() - 120)
  cached_time <- file.info(output)$mtime
  stage()
  expect_identical(file.info(output)$mtime, cached_time)

  writeLines("partial", output)
  stage()
  expect_equal(digest::digest(file = output, algo = "sha256"), expected_sha)
  expect_false(any(grepl("partial-", list.files(dirname(output)), fixed = TRUE)))

  writeLines(rep(sub("^chr1", "chr2", bed_line), 100000L), output, useBytes = TRUE)
  expect_equal(file.info(output)$size, expected_bytes)
  stage()
  expect_equal(digest::digest(file = output, algo = "sha256"), expected_sha)
  receipt <- utils::read.delim(paste0(output, ".provenance.tsv"), colClasses = "character")
  expect_equal(receipt$value[receipt$field == "supplier_identity"], registry$supplier_identity[[3L]])
  expect_false(any(grepl("partial-", list.files(dirname(output)), fixed = TRUE)))
}

test_tabix_split_staging()
