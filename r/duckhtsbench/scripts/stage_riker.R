#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(Sys.getenv("DUCKHTS_ROOT", unset = "."), mustWork = TRUE)
package_root <- file.path(root, "r", "duckhtsbench")
Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(package_root, "inst", "benchmark_registry.tsv"))
for (file in c("registry.R", "stage.R", "riker.R")) source(file.path(package_root, "R", file))
base_dir <- if (length(args)) args[[1L]] else duckhts_bench_cache_path("benchmarks/riker-wgs")
threads <- if (length(args) > 1L) as.integer(args[[2L]]) else 8L
cat(duckhts_bench_stage_riker(base_dir, threads, Sys.getenv("SAMTOOLS_BIN", unset = Sys.which("samtools"))), "\n", sep = "")
