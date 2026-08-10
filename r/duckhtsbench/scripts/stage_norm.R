#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(Sys.getenv("DUCKHTS_ROOT", unset = "."), mustWork = TRUE)
package_root <- file.path(root, "r", "duckhtsbench")
if (!nzchar(Sys.getenv("DUCKHTSBENCH_REGISTRY", unset = ""))) Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(package_root, "inst", "benchmark_registry.tsv"))
for (file in c("registry.R", "stage.R", "norm.R")) source(file.path(package_root, "R", file))
output <- if (length(args)) args[[1L]] else duckhts_bench_artifact_path("norm_hg00096_chr22_20m_30m")
output <- duckhts_bench_stage_norm(output, Sys.getenv("BCFTOOLS_BIN", unset = Sys.which("bcftools")), Sys.getenv("TABIX_BIN", unset = Sys.which("tabix")), as.integer(Sys.getenv("NORM_GVCF_THREADS", unset = "2")))
cat("NORM_GVCF_VCF=", output, "\nNORM_GVCF_PROVENANCE=", paste0(output, ".provenance.tsv"), "\n", sep = "")
