library(tinytest)

local({
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  })
  directory <- tempfile("projection-stage-")
  dir.create(file.path(directory, "test/data/duckvep"), recursive = TRUE)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  plan <- duckhts_bench_stage_plan("duckvep-projection")
  expect_equal(nrow(plan), 3L)
  expect_equal(sort(plan$id), sort(c("projection_reference", "projection_reference_fai", "projection_model_gff")))
  sources <- file.path(directory, sub("^repo:", "", plan$locator))
  for (i in seq_along(sources)) {
    writeLines(paste("fixture", i), sources[[i]])
    plan$supplier_identity[[i]] <- paste0("bytes=", file.info(sources[[i]])$size,
      ";md5=", tools::md5sum(sources[[i]]))
  }
  registry <- file.path(directory, "registry.tsv")
  utils::write.table(plan, registry, sep = "\t", row.names = FALSE, quote = FALSE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry, DUCKHTS_CACHE_DIR = file.path(directory, "cache"))
  paths <- duckhts_bench_stage_repository_fixtures(directory, "duckvep-projection")
  expect_equal(unname(tools::md5sum(paths)), unname(tools::md5sum(sources)))
  expect_true(all(file.exists(paste0(paths, ".provenance.tsv"))))
  writeLines("bad source", sources[[1]])
  expect_error(duckhts_bench_stage_repository_fixtures(directory, "duckvep-projection"),
    pattern = "identity does not match")
  plan$locator[[1]] <- "repo:test/data/../outside"
  utils::write.table(plan, registry, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(duckhts_bench_stage_repository_fixtures(directory, "duckvep-projection"))
})
