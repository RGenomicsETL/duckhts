#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(Sys.getenv("DUCKHTS_ROOT", unset = "."), mustWork = TRUE)
package_root <- file.path(root, "r", "duckhtsbench")
if (!nzchar(Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = ""))) Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(package_root, "inst", "benchmark_registry.tsv"))
for (file in c("registry.R", "stage.R", "duckbedqc.R")) source(file.path(package_root, "R", file))
destination <- if (length(args)) args[[1L]] else duckhts_bench_artifact_path("duckbedqc_118fc21")
destination <- duckhts_bench_stage_duckbedqc(destination, Sys.getenv("GIT_BIN", unset = Sys.which("git")))
cat("DUCKBEDQC_DIR=", destination, "\nDUCKBEDQC_PROVENANCE=", file.path(destination, "provenance.tsv"), "\n", sep = "")
