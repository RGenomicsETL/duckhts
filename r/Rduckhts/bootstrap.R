#!/usr/bin/env Rscript
# Bootstrap script: copies extension sources into inst/duckhts_extension/
# Run from r/duckhts/:  Rscript bootstrap.R
# Or with explicit root: Rscript bootstrap.R /path/to/duckhts

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript bootstrap.R /path/to/duckhts", call. = FALSE)
}
repo_root <- args[1]

# Source the function directly (package may not be installed yet)
source(file.path("R", "bootstrap.R"))

duckhts_bootstrap(repo_root = repo_root)
