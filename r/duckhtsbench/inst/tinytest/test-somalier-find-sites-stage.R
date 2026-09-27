library(tinytest)

local({
  plan <- duckhts_bench_stage_plan("somalier-find-sites")
  expect_equal(sort(plan$id), sort(c("somalier_find_sites_1000g",
                                     "somalier_find_sites_gnomad")))
  repo <- tempfile("find-sites-repo-")
  cache <- tempfile("find-sites-cache-")
  old_cache <- Sys.getenv("DUCKHTS_CACHE_DIR", unset = NA_character_)
  Sys.setenv(DUCKHTS_CACHE_DIR = cache)
  on.exit({
    if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else
      Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)
    unlink(c(repo, cache), recursive = TRUE)
  }, add = TRUE)
  for (i in seq_len(nrow(plan))) {
    source <- system.file("extdata", basename(plan$locator[[i]]),
                          package = "duckhtsbench", mustWork = TRUE)
    destination <- file.path(repo, sub("^repo:", "", plan$locator[[i]]))
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    expect_true(file.copy(source, destination))
  }
  paths <- duckhts_bench_stage_repository_fixtures(repo, "somalier-find-sites")
  expect_equal(length(paths), 2L)
  expect_true(all(file.exists(paths)))
  for (i in seq_len(nrow(plan))) {
    source <- file.path(repo, sub("^repo:", "", plan$locator[[i]]))
    expect_equal(unname(tools::md5sum(paths[[plan$id[[i]]]])),
                 unname(tools::md5sum(source)))
  }
  expect_identical(duckhts_bench_stage_repository_fixtures(
    repo, "somalier-find-sites"), paths)
})
