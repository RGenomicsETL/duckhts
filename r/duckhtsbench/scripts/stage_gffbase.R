#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(Sys.getenv("DUCKHTS_ROOT", unset = "."), mustWork = TRUE)
package_root <- file.path(root, "r", "duckhtsbench")
if (!nzchar(Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = ""))) Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(package_root, "inst", "benchmark_registry.tsv"))
for (file in c("registry.R", "stage.R", "gffbase.R")) source(file.path(package_root, "R", file))
site_dir <- if (length(args)) args[[1L]] else duckhts_bench_artifact_path("gffbase_010")
site_dir <- duckhts_bench_stage_gffbase(site_dir, Sys.getenv("PYTHON_BIN", unset = Sys.which("python3")))
cat("GFFBASE_SITE=", site_dir, "\nGFFBASE_PROVENANCE=", file.path(site_dir, "provenance.tsv"), "\n", sep = "")
