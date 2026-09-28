library(tinytest)

local({
  plan <- duckhts_bench_stage_plan("somalier-find-sites")
  expect_equal(sort(plan$id), sort(c("somalier_find_sites_1000g",
                                     "somalier_find_sites_gnomad")))
  repo <- tempfile("find-sites-repo-")
  cache <- tempfile("find-sites-cache-")
  registry_file <- tempfile(fileext = ".tsv")
  old_cache <- Sys.getenv("DUCKHTS_CACHE_DIR", unset = NA_character_)
  old_registry <- Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = NA_character_)
  on.exit({
    if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else
      Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)
    if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else
      Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
    unlink(c(repo, cache), recursive = TRUE)
    unlink(registry_file)
  }, add = TRUE)
  for (i in seq_len(nrow(plan))) {
    source <- file.path(repo, sub("^repo:", "", plan$locator[[i]]))
    dir.create(dirname(source), recursive = TRUE, showWarnings = FALSE)
    writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                 paste0("chr1\t", i, "\t.\tA\tG\t.\tPASS\tAF=0.5;AN=6")), source)
    plan$supplier_identity[[i]] <- paste0(
      "md5=", unname(tools::md5sum(source)),
      ";bytes=", file.info(source)$size)
  }
  utils::write.table(plan, registry_file, sep = "\t", row.names = FALSE,
                     quote = FALSE)
  Sys.setenv(DUCKHTS_CACHE_DIR = cache, DUCKHTSBENCH_REGISTRY = registry_file)
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
