#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[[1L]])) {
  stop("usage: audit_rduckhts_source_tarball.R <Rduckhts source tarball>",
       call. = FALSE)
}

tarball <- normalizePath(args[[1L]], mustWork = TRUE)
files <- utils::untar(tarball, list = TRUE)
compiled <- grep("[.](o|a|dll|so|dylib)$", files, value = TRUE)
if (length(compiled)) {
  stop(
    "source tarball contains compiled files: ",
    paste(compiled, collapse = ", "),
    call. = FALSE
  )
}

description_path <- "Rduckhts/DESCRIPTION"
if (!description_path %in% files) {
  stop("source tarball does not contain ", description_path, call. = FALSE)
}

extract_dir <- tempfile("rduckhts_source_audit_")
dir.create(extract_dir)
on.exit(unlink(extract_dir, recursive = TRUE, force = TRUE), add = TRUE)
utils::untar(tarball, files = description_path, exdir = extract_dir)
description <- read.dcf(file.path(extract_dir, description_path))
needs_compilation <- unname(description[1L, "NeedsCompilation"])
if (!identical(needs_compilation, "yes")) {
  stop(
    "Rduckhts source tarball must declare NeedsCompilation: yes; got ",
    dQuote(needs_compilation),
    call. = FALSE
  )
}

cat("Source tarball is clean and declares NeedsCompilation: yes\n")
