library(tinytest)

local({
  id <- "somalier_find_sites_chr22_source"
  plan <- duckhts_bench_stage_plan("somalier-find-sites-chr22")
  expect_equal(plan$id, id)
  expect_equal(plan$transform, "direct_download")
  expect_true(grepl("chr22.filtered.shapeit2-duohmm-phased.vcf.gz$",
                    plan$locator))

  cache <- tempfile("find-sites-chr22-cache-")
  source <- tempfile(fileext = ".vcf.gz")
  registry_file <- tempfile(fileext = ".tsv")
  old_cache <- Sys.getenv("DUCKHTS_CACHE_DIR", unset = NA_character_)
  old_registry <- Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = NA_character_)
  on.exit({
    if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else
      Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)
    if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else
      Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
    unlink(c(cache, source, registry_file), recursive = TRUE)
  }, add = TRUE)
  writeLines("offline chromosome staging probe", source)
  plan$locator <- paste0("file://", source)
  plan$supplier_identity <- paste0("md5=", unname(tools::md5sum(source)),
                                   ";bytes=", file.info(source)$size)
  utils::write.table(plan, registry_file, sep = "\t", row.names = FALSE,
                     quote = FALSE)
  Sys.setenv(DUCKHTS_CACHE_DIR = cache, DUCKHTSBENCH_REGISTRY = registry_file)
  staged <- duckhts_bench_fetch(id)
  expect_true(file.exists(staged))
  expect_true(file.exists(paste0(staged, ".provenance.tsv")))
  expect_equal(unname(tools::md5sum(staged)), unname(tools::md5sum(source)))
  expect_identical(duckhts_bench_fetch(id), staged)
  writeLines("corrupt", staged)
  expect_identical(duckhts_bench_fetch(id), staged)
  expect_equal(unname(tools::md5sum(staged)), unname(tools::md5sum(source)))
})
