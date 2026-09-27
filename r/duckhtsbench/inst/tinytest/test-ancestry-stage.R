library(tinytest)

local({
  plan <- duckhts_bench_stage_plan("ancestry-projection")
  expect_equal(plan$id, "ancestry_correction")
  expect_equal(plan$transform, "copy_committed_fixture")
  root <- tempfile("ancestry-stage-")
  dir.create(file.path(root, "test/data"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE))
  source <- file.path(root, "test/data/correction.tsv")
  writeLines(c("pc\tcoefficient", "1\t1"), source)
  plan$locator <- "repo:test/data/correction.tsv"
  plan$supplier_identity <- paste0("bytes=", file.info(source)$size,
                                   ";md5=", unname(tools::md5sum(source)))
  registry <- file.path(root, "registry.tsv")
  utils::write.table(plan, registry, sep = "\t", row.names = FALSE, quote = FALSE)
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  }, add = TRUE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry, DUCKHTS_CACHE_DIR = file.path(root, "cache"))
  staged <- duckhts_bench_stage_repository_fixtures(root, "ancestry-projection")
  expect_equal(unname(tools::md5sum(staged)), unname(tools::md5sum(source)))
  expect_equal(duckhts_bench_stage_repository_fixtures(root, "ancestry-projection"), staged)
  writeLines("corrupted", staged)
  expect_error(duckhts_bench_stage_repository_fixtures(root, "ancestry-projection"),
               pattern = "identity does not match")
})

local({
  plan <- rbind(duckhts_bench_stage_plan("ancestry-reference"),
                duckhts_bench_stage_plan("ancestry-1000g-chr22"),
                duckhts_bench_stage_plan("ancestry-30x-cram"))
  expect_equal(plan$id, c("ancestry_ref_freqs", "ancestry_projection",
                          "ancestry_epilepsy", "ancestry_1000g_chr22",
                          "ancestry_1000g_panel", "ancestry_30x_na18507",
                          "ancestry_30x_hg00403"))
  expect_true(all(plan$transform == "direct_download"))
  expect_true(all(grepl("^https://", plan$locator)))
  root <- tempfile("ancestry-reference-stage-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"),
                         unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, as.list(previous[name]))
    }
  }, add = TRUE)
  plan$cache_relpath <- file.path("benchmarks/ancestry-reference", basename(plan$cache_relpath))
  plan$supplier_identity <- ""
  for (i in seq_len(nrow(plan))) {
    path <- file.path(root, plan$cache_relpath[[i]])
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(paste("reference", i), path)
    plan$supplier_identity[[i]] <- paste0(
      "bytes=", file.info(path)$size, ";sha256=",
      digest::digest(file = path, algo = "sha256")
    )
  }
  registry <- file.path(root, "registry.tsv")
  utils::write.table(plan, registry, sep = "\t", row.names = FALSE, quote = FALSE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry, DUCKHTS_CACHE_DIR = root)
  for (id in plan$id) {
    path <- duckhts_bench_artifact_path(id)
    expect_identical(duckhts_bench_fetch(id), path)
    expect_true(file.exists(paste0(path, ".provenance.tsv")))
    expect_equal(duckhts_bench_validate_identity(id, path), TRUE)
  }
  writeLines("changed", duckhts_bench_artifact_path(plan$id[[1L]]))
  expect_error(duckhts_bench_validate_identity(plan$id[[1L]]),
               pattern = "identity does not match")
})
