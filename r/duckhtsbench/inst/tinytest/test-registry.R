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
