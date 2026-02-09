#!/usr/bin/env Rscript
# Bootstrap script: copies extension sources into inst/duckhts_extension/
# Run from r/Rduckhts/:  Rscript bootstrap.R
# Or with explicit root: Rscript bootstrap.R /path/to/duckhts

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) args[1] else NULL

# Source the function directly (package may not be installed yet)
source(file.path("R", "bootstrap.R"))

duckhts_bootstrap(repo_root = repo_root)
