#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("must be invoked as an R script", call. = FALSE)
package_dir <- normalizePath(
  file.path(dirname(sub("^--file=", "", script_arg)), ".."),
  mustWork = TRUE
)
if (!nzchar(Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = ""))) {
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(package_dir, "inst", "benchmark_registry.tsv"))
}
source(file.path(package_dir, "R", "registry.R"))
source(file.path(package_dir, "R", "stage.R"))

workload <- "giab-benchmark-vcf"
if (identical(args, "--plan")) {
  print(duckhts_bench_stage_plan(workload), row.names = FALSE)
  quit(status = 0L)
}
if (length(args)) stop("usage: stage_giab_benchmark_vcf.R [--plan]", call. = FALSE)
plan <- duckhts_bench_stage_plan(workload)
if (nrow(plan) != 1L || !identical(plan$transform, "direct_download")) {
  stop("GIAB benchmark VCF must be one direct-download artifact", call. = FALSE)
}
duckhts_bench_fetch(plan$id[[1L]])
message("GIAB benchmark VCF staged under ", duckhts_bench_cache_dir())
