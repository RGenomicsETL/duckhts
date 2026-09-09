#!/usr/bin/env Rscript
# Maintainer operation: retain two independently observed transcript records.
source('test/duckvep/conformance/haplotype_model_differential.R')
source('scripts/duckvep_evidence.R')

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
artifact <- duckvep_evidence_read_artifact(args[1L],
  c('inputs.rds', 'actual.rds', 'oracle.stdout', 'phase.jsonl', 'comparisons.rds'))
directory <- artifact$directory
inputs <- readRDS(file.path(directory, 'inputs.rds'))
inputs$models <- inputs$models[1:2, ]
inputs$cases <- inputs$cases[1L, , drop = FALSE]
inputs$exons <- inputs$exons[inputs$exons$transcript_index < 2L, ]
inputs$records <- inputs$records[inputs$records$seq_region == 0L, ]
inputs$geometry <- NULL
actual <- readRDS(file.path(directory, 'actual.rds'))
actual <- actual[actual$transcript_index < 2L, ]
observations <- function(file) {
  rows <- lapply(readLines(file.path(directory, file), n = 2L),
    jsonlite::fromJSON, simplifyVector = FALSE)
  names(rows) <- vapply(rows, `[[`, '', 'transcript')
  stopifnot(setequal(names(rows), inputs$models$transcript))
  rows
}
fixture <- list(inputs = inputs, actual = actual, oracle = observations('oracle.stdout'),
  phase = observations('phase.jsonl'), source_revision = artifact$receipt$source_revision,
  source_receipt = artifact$relative,
  source_receipt_sha256 = duckvep_evidence_sha256(artifact$receipt_path),
  oracle_revisions = artifact$receipt$oracle_revisions)
checked <- with(fixture, model_comparisons(inputs, actual, oracle, phase))
stopifnot(identical(checked$comparisons,
  readRDS(file.path(directory, 'comparisons.rds'))[1:2]), checked$metrics$failures == 0L)
saveRDS(fixture, 'test/duckvep/conformance/data/model_publication_fixture.rds', version = 2)
