#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE)
package_root <- file.path(root, "r", "duckhtsbench")
source(file.path(package_root, "R", "registry.R"))
duckhts_bench_check_portability(file.path(root, "benchmarks"))
