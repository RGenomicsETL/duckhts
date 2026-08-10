library(tinytest)

registry_path <- Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = "")
if (!nzchar(registry_path)) {
  registry_path <- system.file("benchmark_registry.tsv", package = "duckhtsbench")
}
expect_true(file.exists(registry_path))
Sys.setenv(DUCKHTSBENCH_REGISTRY = registry_path)

registry <- duckhts_bench_registry()
expect_true(all(c("id", "workload", "release", "locator", "access", "cache_relpath", "transform", "consumer") %in% names(registry)))
expect_true(all(nzchar(registry$id)))
expect_equal(length(unique(registry$id)), nrow(registry))
expect_true(all(!grepl("^/", registry$cache_relpath)))
expect_true(all(!grepl("\\.\\.", registry$cache_relpath)))
expect_true(any(registry$id == "revel_v13_grch37"))
expect_match(duckhts_bench_artifact_path("revel_v13_grch37"), "revel_grch37\\.parquet$")
expect_equal(nrow(duckhts_bench_stage_plan("variantkey-providers")), nrow(registry))

old_registry <- Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = NA_character_)
old_cache <- Sys.getenv("DUCKHTS_CACHE_DIR", unset = NA_character_)
tmp <- tempfile("duckhtsbench-fetch-")
dir.create(tmp)
source_path <- file.path(tmp, "source.txt")
writeLines("registry fetch", source_path)
mini_registry <- file.path(tmp, "registry.tsv")
writeLines(c(
  paste(names(registry), collapse = "\t"),
  paste(c("fixture", "fixture", "fixture", "fixture-1", paste0("file://", source_path), "public", "fixture/output.txt", "direct_download", "tinytest", "1"), collapse = "\t")
), mini_registry)
Sys.setenv(DUCKHTSBENCH_REGISTRY = mini_registry, DUCKHTS_CACHE_DIR = file.path(tmp, "cache"))
output <- duckhts_bench_fetch("fixture")
expect_equal(readLines(output), "registry fetch")
expect_true(file.exists(paste0(output, ".provenance.tsv")))
if (is.na(old_registry)) Sys.unsetenv("DUCKHTSBENCH_REGISTRY") else Sys.setenv(DUCKHTSBENCH_REGISTRY = old_registry)
if (is.na(old_cache)) Sys.unsetenv("DUCKHTS_CACHE_DIR") else Sys.setenv(DUCKHTS_CACHE_DIR = old_cache)
