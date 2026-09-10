#!/usr/bin/env Rscript
# Run against an actual GitHub-issued bundle; no successful verifier is mocked.
source('scripts/duckvep_evidence.R')

main <- function(args) {
  stopifnot(length(args) == 2L)
  receipt <- jsonlite::fromJSON(args[1L], simplifyVector = FALSE)
  verify <- function(path, revision = receipt$source_revision, bundle = args[2L])
    duckvep_evidence_verify_ci_receipt(path, revision, bundle)
  verify(args[1L])
  verify(file.path(dirname(args[1L]), 'duckhts.duckdb_extension'))
  directory <- tempfile('duckvep-attestation-test-')
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  changed <- file.path(directory, 'receipt.json')
  rejected <- 0L
  fails <- function(expr) {
    result <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(result, 'error'))
    rejected <<- rejected + 1L
  }
  for (field in c('source_revision', 'extension_build_binding', 'tracked_changes', 'sha256')) {
    forged <- receipt
    forged[[field]] <- switch(field, source_revision = strrep('0', 40L),
      extension_build_binding = 'diagnostic_unbound', tracked_changes = list('src/geno_reader.c'),
      sha256 = list(forged = strrep('0', 64L)))
    jsonlite::write_json(forged, changed, auto_unbox = TRUE)
    fails(verify(changed))
  }
  fails(verify(args[1L], strrep('0', 40L)))
  fails(verify(args[1L], bundle = file.path(directory, 'missing.jsonl')))
  invalid <- file.path(directory, 'invalid.jsonl')
  writeLines('{}', invalid)
  fails(verify(args[1L], bundle = invalid))
  verify(args[1L])
  cat('GitHub execution provenance: authentic receipt accepted;', rejected,
    'altered receipts/identities/bundles rejected.\n')
}

main(commandArgs(trailingOnly = TRUE))
