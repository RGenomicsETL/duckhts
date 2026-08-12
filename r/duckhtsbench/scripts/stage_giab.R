#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) && args[[1L]] == "--root") args[[2L]] else ".", mustWork = TRUE)
package_root <- file.path(root, "r", "duckhtsbench")
Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(package_root, "inst", "benchmark_registry.tsv"))
source(file.path(package_root, "R", "registry.R"))
source(file.path(package_root, "R", "stage.R"))
source(file.path(package_root, "R", "giab.R"))
samples <- strsplit(Sys.getenv("GIAB_SAMPLES", unset = "HG001,HG002,HG006"), ",", fixed = TRUE)[[1L]]
paths <- duckhts_bench_stage_giab(trimws(samples), Sys.getenv("BCFTOOLS_BIN", unset = Sys.which("bcftools")))
cat(paste(paths, collapse = "\n"), "\n", sep = "")
