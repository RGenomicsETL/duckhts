#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
run_conformance <- "--conformance" %in% args

suppressPackageStartupMessages({
    library(duckhts)
})

duckhts_bootstrap(vendor_conformance = run_conformance)
