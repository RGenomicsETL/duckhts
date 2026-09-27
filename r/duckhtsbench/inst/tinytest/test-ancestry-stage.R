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
