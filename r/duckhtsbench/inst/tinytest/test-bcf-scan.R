library(tinytest)

test_bcf_scan_stage <- function() {
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  })
  directory <- tempfile("bcf-scan-stage-")
  dir.create(file.path(directory, "test/data"), recursive = TRUE)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  plan <- duckhts_bench_stage_plan("bcf-scan-init")
  expect_equal(nrow(plan), 12L)
  expect_true(all(grepl("bytes=[0-9]+;md5=[a-f0-9]{32}", plan$supplier_identity)))
  plan <- plan[1L, ]
  source <- file.path(directory, sub("^repo:", "", plan$locator))
  writeLines("fixture bytes", source)
  plan$supplier_identity <- paste0("bytes=", file.info(source)$size, ";md5=", tools::md5sum(source))
  registry <- file.path(directory, "registry.tsv")
  utils::write.table(plan, registry, sep = "\t", row.names = FALSE, quote = FALSE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry, DUCKHTS_CACHE_DIR = file.path(directory, "cache"))
  paths <- duckhts_bench_stage_bcf_scan(directory)
  expect_equal(names(paths), plan$id)
  expect_equal(unname(tools::md5sum(paths)), unname(tools::md5sum(source)))
  expect_true(all(file.exists(paste0(paths, ".provenance.tsv"))))
  writeLines("poisoned cached copy", paths[[1]])
  expect_equal(duckhts_bench_stage_bcf_scan(directory), paths)
  expect_equal(readLines(paths[[1]]), "fixture bytes")
  writeLines("wrong source revision", source)
  expect_error(duckhts_bench_stage_bcf_scan(directory), pattern = "identity does not match")
  unlink(source)
  expect_error(duckhts_bench_stage_bcf_scan(directory), pattern = "missing registered fixture")
}

test_bcf_scan_stage()
