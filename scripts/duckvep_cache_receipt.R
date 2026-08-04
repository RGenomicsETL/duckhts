#!/usr/bin/env Rscript

suppressMessages(library(optparse))

root <- normalizePath(
  system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)[[1L]],
  mustWork = TRUE
)
source(file.path(root, "scripts", "duckvep_evidence.R"), local = TRUE)

`%||%` <- function(left, right) {
  if (is.null(left) || !length(left)) right else left
}

parser <- OptionParser(
  usage = "%prog --cache-dir DIR --output FILE [options]"
)
parser <- add_option(parser, "--cache-dir", dest = "cache_dir")
parser <- add_option(parser, "--species", default = "homo_sapiens")
parser <- add_option(parser, "--cache-version", dest = "cache_version")
parser <- add_option(parser, "--assembly")
parser <- add_option(parser, "--source-url", dest = "source_url")
parser <- add_option(parser, "--source-identity", dest = "source_identity")
parser <- add_option(parser, "--output")
parser <- add_option(
  parser,
  "--overwrite",
  action = "store_true",
  default = FALSE
)
options <- parse_args(parser)

required <- c(
  "cache_dir",
  "species",
  "cache_version",
  "assembly",
  "source_url",
  "source_identity",
  "output"
)
missing <- required[
  !nzchar(vapply(
    required,
    function(field) options[[field]] %||% "",
    character(1L)
  ))
]
if (length(missing)) {
  stop(
    "missing option(s): ",
    paste(paste0("--", gsub("_", "-", missing)), collapse = ", "),
    call. = FALSE
  )
}

path <- duckvep_evidence_write_cache_receipt(
  options$output,
  options$cache_dir,
  options$species,
  options$cache_version,
  options$assembly,
  options$source_url,
  options$source_identity,
  options$overwrite
)
validated <- duckvep_evidence_validate_cache_receipt(
  path,
  options$cache_dir,
  options$species,
  options$cache_version,
  options$assembly
)
cat(
  "receipt=",
  path,
  "\n",
  "entries=",
  validated$entries,
  "\n",
  "bytes=",
  validated$bytes,
  "\n",
  "inventory_sha256=",
  validated$inventory_sha256,
  "\n",
  sep = ""
)
