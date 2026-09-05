library(tinytest)

test_multi_region_stage <- function() {
  previous <- Sys.getenv(c("DUCKHTSBENCH_REGISTRY", "DUCKHTS_CACHE_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  })
  directory <- tempfile("multi-region-stage-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  registry <- duckhts_bench_stage_plan("multi-region-readers")
  expect_equal(nrow(registry), 5L)
  expect_match(registry$supplier_identity[[1]], "records=2000000", fixed = TRUE)
  expect_equal(registry$locator[-1], paste0("artifact:multi_region_", c("vcf", "vcf_gz", "vcf_gz", "bcf")))
  registry$supplier_identity[[1]] <- "records=17"
  registry_file <- file.path(directory, "registry.tsv")
  utils::write.table(registry, registry_file, sep = "\t", quote = FALSE, row.names = FALSE)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_file, DUCKHTS_CACHE_DIR = directory)
  paths <- duckhts_bench_stage_multi_region()
  expect_true(all(file.exists(paths)))
  expect_true(all(file.exists(paste0(paths, ".provenance.tsv"))))
  for (path in paths[c("vcf_gz", "bcf")]) {
    positions <- system2("bcftools", c("query", "-f", shQuote("%POS\n"), shQuote(path)), stdout = TRUE)
    expect_equal(as.integer(positions), seq_len(17L))
  }
  first <- tools::md5sum(paths)
  expect_equal(tools::md5sum(duckhts_bench_stage_multi_region()), first)
}

test_multi_region_stage()
