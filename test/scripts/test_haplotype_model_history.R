#!/usr/bin/env Rscript
# Publication integrity, not a synthetic biological oracle.
source('test/duckvep/conformance/haplotype_model_differential.R')
source('scripts/duckvep_evidence.R')

main <- function() {
  fixture <- readRDS('test/duckvep/conformance/data/model_publication_fixture.rds')
  stopifnot(fixture$source_revision == '0857ec1faf5736299a4341f73ffa557d8a9691ee',
    identical(unlist(fixture$oracle_revisions), model_oracle_revisions))
  directory <- tempfile('model-history-test-', 'test/duckvep/conformance/results')
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  save <- function(actual = fixture$actual) {
    checked <- model_comparisons(fixture$inputs, actual, fixture$oracle, fixture$phase)
    saveRDS(actual, file.path(directory, 'actual.rds'))
    saveRDS(checked$comparisons, file.path(directory, 'comparisons.rds'))
    write.csv(checked$summary, file.path(directory, 'summary.csv'), row.names = FALSE)
    for (name in c('controls', 'metadata_controls')) write.csv(
      data.frame(control = names(checked[[name]]), rejected = checked[[name]]),
      file.path(directory, paste0(name, '.csv')), row.names = FALSE)
    checked
  }
  checked <- save()
  receipt <- checked$metrics
  required <- c('actual.rds', 'comparisons.rds', 'summary.csv', 'controls.csv', 'metadata_controls.csv')
  rehash <- function(value = receipt) {
    files <- file.path(directory, required)
    value$sha256 <- as.list(vapply(files, duckvep_evidence_sha256, ''))
    jsonlite::write_json(value, file.path(directory, 'receipt.json'), auto_unbox = TRUE)
  }
  verify <- function() {
    artifact <- duckvep_evidence_read_artifact(directory, required)
    with(fixture, model_check_saved(directory, inputs,
      readRDS(file.path(directory, 'actual.rds')), oracle, phase, artifact$receipt))
  }
  rejected <- 0L
  fails <- function(expr) {
    stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), 'error'))
    rejected <<- rejected + 1L
  }
  rehash()
  stopifnot(verify()$metrics$failures == 0L, checked$metrics$carriers == 12L,
    all(checked$controls), length(checked$controls) == 42L,
    all(checked$metadata_controls), length(checked$metadata_controls) == 18L)
  bad <- checked$comparisons
  bad[[1L]]$expected <- bad[[1L]]$observed <- list(list(cds = 'forged', count = 6))
  saveRDS(bad, file.path(directory, 'comparisons.rds'))
  fails(verify())
  rehash()
  fails(verify())
  save()
  # Coherent group-preserving swaps still contradict the upstream carrier lanes.
  swapped <- fixture$actual
  pair <- which(swapped$transcript_index == 1L)[1:2]
  for (field in c('cds', 'protein', 'coding_blocks'))
    swapped[[field]][pair] <- swapped[[field]][rev(pair)]
  saveRDS(swapped, file.path(directory, 'actual.rds'))
  rehash()
  fails(verify())
  genuine <- save(swapped)
  stopifnot(genuine$metrics$failures == 1L, genuine$metrics$replay_lane_failures == 1L)
  rehash(genuine$metrics)
  stopifnot(verify()$metrics$failures == 1L)
  save()
  for (name in names(receipt)[vapply(receipt, is.numeric, TRUE)]) {
    altered <- receipt
    altered[[name]] <- altered[[name]] + 1
    rehash(altered)
    fails(verify())
    altered[[name]] <- NULL
    rehash(altered)
    fails(verify())
  }
  rehash()
  summary <- checked$summary
  summary$transcript <- rev(summary$transcript)
  write.csv(summary, file.path(directory, 'summary.csv'), row.names = FALSE)
  rehash()
  fails(verify())
  save()
  for (name in c('controls', 'metadata_controls')) {
    write.csv(data.frame(control = 'invented', rejected = TRUE),
      file.path(directory, paste0(name, '.csv')), row.names = FALSE)
    rehash()
    fails(verify())
    save()
  }
  for (field in c('sequence_flags', 'nominal_length_diff')) {
    bad <- fixture$actual
    i <- which(bad$transcript_index == 1L)[1L]
    bad[[field]][i] <- bad[[field]][i] + 1
    saveRDS(bad, file.path(directory, 'actual.rds'))
    rehash()
    fails(verify())
    save()
  }
  unknown <- fixture$actual[1L, ]
  unknown$transcript_index <- 2L
  saveRDS(rbind(fixture$actual, unknown), file.path(directory, 'actual.rds'))
  rehash()
  fails(verify())
  save()
  for (name in c('oracle', 'phase')) {
    bad <- fixture[[name]][1L]
    fails(if (name == 'oracle') model_comparisons(fixture$inputs, fixture$actual, bad, fixture$phase)
      else model_comparisons(fixture$inputs, fixture$actual, fixture$oracle, bad))
    bad <- c(fixture[[name]], fixture[[name]][1L])
    fails(if (name == 'oracle') model_comparisons(fixture$inputs, fixture$actual, bad, fixture$phase)
      else model_comparisons(fixture$inputs, fixture$actual, fixture$oracle, bad))
  }
  stopifnot(identical(checked, model_comparisons(fixture$inputs, fixture$actual,
    rev(fixture$oracle), rev(fixture$phase))))
  # A passing subset cannot establish the full seeded model/input denominator.
  options <- list(rare_per_stratum = 0L, geometry_per_stratum = 0L,
    interaction_per_stratum = 0L, length_per_stratum = 0L)
  reference <- readLines('test/data/duckvep/haplotype_benchmark.fa')[2L]
  generated <- model_inputs(173L, options, reference)
  stopifnot(nrow(generated$inputs$models) == 1536L, nrow(generated$coverage) == 264L,
    !identical(generated$inputs, fixture$inputs))
  for (value in c(NA_real_, -1, 0.5, 1e8)) {
    invalid <- options
    invalid$rare_per_stratum <- value
    fails(model_inputs(173L, invalid, reference))
  }
  cat('Model publication: two pinned observations, all 60 controls, retained lane-swap failure;',
    rejected, 'forgeries/invalid inputs rejected.\n')
}

main()
