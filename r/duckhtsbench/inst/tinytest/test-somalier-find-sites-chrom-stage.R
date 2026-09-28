library(tinytest)

local({
  workloads <- c("somalier-find-sites-chr22", "somalier-find-sites-chrx")
  plans <- lapply(workloads, duckhts_bench_stage_plan)
  plan <- do.call(rbind, plans)
  expect_equal(plan$id, c("somalier_find_sites_chr22_source",
                          "somalier_find_sites_chrx_source"))
  expect_true(all(plan$transform == "direct_download"))
  expect_true(grepl("chr22.filtered.shapeit2-duohmm-phased.vcf.gz$",
                    plan$locator[[1L]]))
  expect_true(grepl("chrX.filtered.eagle2-phased.v2.vcf.gz$",
                    plan$locator[[2L]]))

  cache <- tempfile("find-sites-chrom-cache-")
  source_dir <- tempfile("find-sites-chrom-input-")
  registry_file <- tempfile(fileext = ".tsv")
  old_cache <- Sys.getenv("DUCKHTS_CACHE_DIR", unset = NA_character_)
  old_registry <- Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = NA_character_)
  on.exit({
    if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else
      Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)
    if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else
      Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
    unlink(c(cache, source_dir, registry_file), recursive = TRUE)
  }, add = TRUE)
  dir.create(source_dir)
  source <- file.path(source_dir, paste0(plan$id, ".vcf.gz"))
  for (i in seq_along(source)) {
    writeLines(paste("offline chromosome staging probe", i), source[[i]])
    plan$locator[[i]] <- paste0("file://", source[[i]])
    plan$supplier_identity[[i]] <- paste0(
      "md5=", unname(tools::md5sum(source[[i]])),
      ";bytes=", file.info(source[[i]])$size)
  }
  utils::write.table(plan, registry_file, sep = "\t", row.names = FALSE,
                     quote = FALSE)
  Sys.setenv(DUCKHTS_CACHE_DIR = cache, DUCKHTSBENCH_REGISTRY = registry_file)
  for (i in seq_along(source)) {
    staged <- duckhts_bench_fetch(plan$id[[i]])
    expect_true(file.exists(staged))
    expect_true(file.exists(paste0(staged, ".provenance.tsv")))
    expect_equal(unname(tools::md5sum(staged)),
                 unname(tools::md5sum(source[[i]])))
    expect_identical(duckhts_bench_fetch(plan$id[[i]]), staged)
    writeLines("corrupt", staged)
    expect_identical(duckhts_bench_fetch(plan$id[[i]]), staged)
    expect_equal(unname(tools::md5sum(staged)),
                 unname(tools::md5sum(source[[i]])))
  }
})
