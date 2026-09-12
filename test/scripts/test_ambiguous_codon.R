#!/usr/bin/env Rscript
# Network-free comparator controls, independent of production annotations.
source('test/duckvep/conformance/ambiguous_codon_differential.R')

# Check object keys before array-to-data-frame conversion can discard duplicates.
read_unique_json <- function(text) {
  check_object_names <- function(value) {
    if (!is.list(value)) return(invisible(TRUE))
    stopifnot(!anyDuplicated(names(value)))
    for (child in value) check_object_names(child)
    invisible(TRUE)
  }
  check_object_names(jsonlite::fromJSON(text, simplifyVector = FALSE))
  jsonlite::fromJSON(text)
}

# Independently enumerate the declared axes and source ordinals. Do not infer
# eligibility or the expected substitutions from the retained case contents.
check_codon_matrix <- function(cases) {
  bases <- c('A', 'C', 'G', 'T', 'N')
  tables <- c(1:6, 9:14, 16, 21:31)
  stopifnot(length(cases) == 24L * 5L^3L)
  case_index <- 0L
  event_index <- 0L
  for (third in bases) for (second in bases) for (first in bases) for (table in tables) {
    case_index <- case_index + 1L
    case <- cases[[case_index]]
    codon <- paste0(first, second, third)
    stopifnot(!anyDuplicated(names(case)),
      setequal(names(case), c('id', 'cds', 'table', 'edits', 'variants')),
      identical(case$id, paste(table, codon, sep = '/')),
      identical(case$cds, paste0('ATG', codon, 'TAA')), is.numeric(case$table),
      identical(as.numeric(case$table), table),
      is.list(case$edits), length(case$edits) == 0L, is.data.frame(case$variants),
      !anyDuplicated(names(case$variants)),
      setequal(names(case$variants), c('id', 'position1', 'reference', 'alternate')))
    variant_index <- 0L
    for (position in 1:3) for (alternate in bases[1:4]) {
      reference <- substr(codon, position, position)
      if (reference == alternate) next
      variant_index <- variant_index + 1L
      event_index <- event_index + 1L
      stopifnot(variant_index <= nrow(case$variants))
      variant <- case$variants[variant_index, ]
      stopifnot(identical(variant$id, as.character(event_index)),
        is.numeric(variant$position1),
        identical(as.numeric(variant$position1), as.numeric(position + 3L)),
        identical(variant$reference, reference), identical(variant$alternate, alternate))
    }
    stopifnot(nrow(case$variants) == variant_index)
  }
  stopifnot(case_index == 3000L, event_index == 28800L)
  invisible(TRUE)
}

codon_matrix_controls <- function(cases) {
  rejected <- function(x) !isTRUE(tryCatch(check_codon_matrix(x), error = function(e) FALSE))
  rejected_json <- function(text) {
    tryCatch({
      changed <- cases
      changed[[1L]] <- read_unique_json(text)
      rejected(changed)
    }, error = function(e) TRUE)
  }
  stopifnot(!rejected(cases))
  controls <- c(missing_case = rejected(cases[-1L]))
  # Tables 1 and 11 share these internal AAA peptide observations. Keep all
  # event ordinals and totals while replacing the table-11 case with a copy.
  changed <- cases
  case_ids <- vapply(cases, `[[`, '', 'id')
  original <- match('1/AAA', case_ids)
  replaced <- match('11/AAA', case_ids)
  changed[[replaced]] <- cases[[original]]
  changed[[replaced]]$id <- 'duplicate/1/AAA'
  changed[[replaced]]$variants$id <- cases[[replaced]]$variants$id
  ids <- unlist(lapply(changed, function(x) x$variants$id))
  stopifnot(length(changed) == 3000L, length(ids) == 28800L, !anyDuplicated(ids),
    !anyDuplicated(vapply(changed, `[[`, '', 'id')))
  controls['count_preserving_case_swap'] <- rejected(changed)
  changed <- cases
  changed[[1L]]$variants[2L, c('position1', 'reference', 'alternate')] <-
    changed[[1L]]$variants[1L, c('position1', 'reference', 'alternate')]
  controls['count_preserving_variant_swap'] <- rejected(changed)
  for (field in c('id', 'cds', 'table', 'edits')) {
    changed <- cases
    changed[[1L]][[field]] <- switch(field, id = 'wrong', cds = 'ATGAACTAA',
      table = 2L, edits = list(list(position1 = 2L, alternate = 'A')))
    controls[paste0('case_', field)] <- rejected(changed)
  }
  for (field in c('id', 'position1', 'reference', 'alternate')) {
    changed <- cases
    value <- changed[[1L]]$variants[[field]][1L]
    changed[[1L]]$variants[[field]][1L] <- if (is.numeric(value)) value + 1L else paste0(value, 'X')
    controls[paste0('variant_', field)] <- rejected(changed)
  }
  # Append raw JSON properties: toJSON renames duplicate R list names, which
  # would test a different input from a source object with repeated properties.
  encoded <- jsonlite::toJSON(cases[[1L]], auto_unbox = TRUE, dataframe = 'rows')
  stopifnot(!rejected_json(encoded))
  for (field in c('id', 'cds', 'table', 'edits', 'variants')) {
    duplicate <- jsonlite::toJSON(cases[[1L]][field], auto_unbox = TRUE, dataframe = 'rows')
    text <- paste0(substr(encoded, 1L, nchar(encoded) - 1L), ',', substring(duplicate, 2L))
    changed <- cases
    changed[[1L]] <- jsonlite::fromJSON(text)
    stopifnot(anyDuplicated(names(changed[[1L]])) > 0L)
    controls[paste0('duplicate_case_', field)] <- rejected(changed) && rejected_json(text)
  }
  variants <- cases[[1L]]$variants
  encoded_variants <- vapply(seq_len(nrow(variants)), function(i)
    jsonlite::toJSON(as.list(variants[i, ]), auto_unbox = TRUE), '')
  case_fields <- cases[[1L]][setdiff(names(cases[[1L]]), 'variants')]
  prefix <- jsonlite::toJSON(case_fields, auto_unbox = TRUE)
  for (field in c('id', 'position1', 'reference', 'alternate')) {
    value <- variants[1L, field]
    duplicate <- setNames(list(if (is.numeric(value)) value + 1L else
      paste0(value, '_duplicate')), field)
    property <- jsonlite::toJSON(duplicate, auto_unbox = TRUE)
    changed_variants <- encoded_variants
    first <- encoded_variants[1L]
    changed_variants[1L] <- paste0(substr(first, 1L, nchar(first) - 1L), ',',
      substring(property, 2L))
    text <- paste0(substr(prefix, 1L, nchar(prefix) - 1L), ',"variants":[',
      paste(changed_variants, collapse = ','), ']}')
    raw <- jsonlite::fromJSON(text, simplifyVector = FALSE)
    stopifnot(anyDuplicated(names(raw$variants[[1L]])) > 0L)
    controls[paste0('duplicate_variant_', field)] <- rejected_json(text)
  }
  stopifnot(all(controls))
  controls
}

expected <- data.frame(event_index = 1:2, allele = c('A', 'C'),
  hgvsp = c('p.Ala2Thr', NA_character_), so = c('missense_variant', 'coding_sequence_variant'))
controls <- rbind(codon_controls(expected), codon_provenance_controls(FALSE),
  codon_provenance_controls(TRUE))
stopifnot(!anyDuplicated(controls$control), all(controls$rejected))
message('Ambiguous-codon comparator: ', nrow(controls), ' corruptions rejected')

codon_check_environment <- function(lines) {
  stopifnot(identical(duckvep_evidence_explicit_packages(lines),
    duckvep_evidence_explicit_packages(readLines(
      'test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt'))))
  invisible(TRUE)
}

codon_check_bundle_identity <- function(directory, receipt_sha256, source_revision, extension_sha256) {
  path <- file.path(directory, 'receipt.json')
  duckvep_evidence_check_receipt_pin(path, receipt_sha256)
  manifest <- read_unique_json(paste(readLines(path), collapse = '\n'))
  stopifnot(setequal(names(manifest), c('source_revision', 'source_binding', 'extension_path',
    'extension_sha256', 'oracle_revisions', 'models', 'source_snvs', 'oracle_pairs', 'decoded_gt',
    'source_gt', 'scope', 'sha256', 'local_receipt_sha256', 'source_sha256')),
    identical(manifest$source_binding, 'diagnostic_unbound'),
    identical(manifest$source_revision, source_revision),
    identical(manifest$extension_sha256, extension_sha256),
    identical(manifest$oracle_revisions,
      c('57ea5c52340acc1f156267f810ad162e26597082', '2fb834b987ede3824e200197a838ce11e91aeb4b')),
    identical(manifest$models, 3000L), identical(manifest$source_snvs, 28800L),
    identical(manifest$oracle_pairs, 28800L), identical(manifest$decoded_gt, 'haploid ALT'),
    identical(manifest$source_gt, '1|1'), identical(manifest$scope, paste0(
      'independent_SO_HGVSp_and_singleton_phased_HGVSp_internal_codons;',
      'complete_pairs_and_oracle_not_all_native_output_fields')))
  # The receipt pin fixes the complete historical source map and binary identity,
  # without substituting whichever checkout or executable is installed today.
  hashes <- unlist(manifest$sha256)
  sources <- unlist(manifest$source_sha256)
  stopifnot(length(sources) == 63L, is.character(sources), !anyDuplicated(names(sources)),
    !anyNA(sources), all(grepl('^[0-9a-f]{64}$', sources)),
    setequal(names(hashes), c('pairs.parquet', 'cases.jsonl.gz', 'oracle.stdout.gz',
      'summary.csv', 'controls.csv', 'environment.stdout')),
    setequal(list.files(directory, all.files = TRUE, no.. = TRUE), c(names(hashes), 'receipt.json')))
  for (name in names(hashes)) stopifnot(identical(unname(hashes[name]),
    duckvep_evidence_sha256(file.path(directory, name))))
  codon_check_environment(readLines(file.path(directory, 'environment.stdout')))
  manifest
}

codon_bundle_identity_controls <- function(directory, receipt_sha256, source_revision,
                                          extension_sha256, pairs, oracle_records) {
  temporary <- tempfile('duckvep-codon-identity-')
  stopifnot(dir.create(temporary))
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  stopifnot(all(file.copy(list.files(directory, full.names = TRUE), temporary)))
  path <- file.path(temporary, 'receipt.json')
  original <- readLines(path)
  receipt <- read_unique_json(paste(original, collapse = '\n'))
  rejected <- function(expression) tryCatch({ force(expression); FALSE }, error = function(e) TRUE)
  check <- function() codon_check_bundle_identity(temporary, receipt_sha256,
    source_revision, extension_sha256)
  stopifnot(!rejected(check()))
  jsonlite::write_json(receipt, path, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  stopifnot(!rejected(check()))
  corruptions <- list(source_binding = 'certified', source_revision = strrep('0', 40L),
    extension_sha256 = strrep('0', 64L), oracle_revisions = rep(strrep('0', 40L), 2L),
    source_sha256 = lapply(receipt$source_sha256, function(x) strrep('0', 64L)),
    local_receipt_sha256 = strrep('0', 64L), models = 3000.0000001,
    scope = 'authenticated_public_annotation')
  checks <- logical()
  for (field in names(corruptions)) {
    changed <- receipt
    changed[[field]] <- corruptions[[field]]
    jsonlite::write_json(changed, path, pretty = TRUE, auto_unbox = TRUE, digits = NA)
    checks[paste0('changed_', field)] <- rejected(check())
    changed[[field]] <- NULL
    jsonlite::write_json(changed, path, pretty = TRUE, auto_unbox = TRUE, digits = NA)
    checks[paste0('missing_', field)] <- rejected(check())
  }
  environment_path <- file.path(temporary, 'environment.stdout')
  environment <- readLines(environment_path)
  changed_environment <- sub('ensembl-vep-116.0-', 'ensembl-vep-117.0-', environment, fixed = TRUE)
  stopifnot(!identical(changed_environment, environment))
  checks['parsed_environment'] <- rejected(codon_check_environment(changed_environment))
  writeLines(changed_environment, environment_path)
  changed <- receipt
  changed$sha256[['environment.stdout']] <- duckvep_evidence_sha256(environment_path)
  jsonlite::write_json(changed, path, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  checks['rehashed_environment'] <- rejected(check())
  writeLines(environment, environment_path)

  # Replace one oracle answer and all eight matching native answers, preserving
  # row counts, equality flags, failure totals and otherwise valid file hashes.
  changed_oracle <- read_unique_json(oracle_records[1L])
  stopifnot(identical(changed_oracle$independent_hgvs$id[1L], '1'))
  changed_oracle$independent_hgvs$hgvsp[1L] <- '1/AAA_protein.1:p.Lys2Ala'
  oracle_records[1L] <- jsonlite::toJSON(changed_oracle, auto_unbox = TRUE, dataframe = 'rows')
  oracle_path <- file.path(temporary, 'oracle.stdout.gz')
  compressed <- gzfile(oracle_path, 'wt')
  tryCatch(writeLines(oracle_records, compressed), finally = close(compressed))
  at <- pairs$event_index == 1L
  stopifnot(sum(at) == 8L, all(pairs$hgvsp_equal[at]),
    all(pairs$hgvsp_actual[at] == 'p.Lys2Gln'), all(pairs$hgvsp_expected[at] == 'p.Lys2Gln'))
  pairs$hgvsp_actual[at] <- 'p.Lys2Ala'
  pairs$hgvsp_expected[at] <- 'p.Lys2Ala'
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(con, 'pairs', pairs)
  pairs_path <- file.path(temporary, 'pairs.parquet')
  DBI::dbExecute(con, paste('COPY pairs TO', DBI::dbQuoteString(con, pairs_path), '(FORMAT PARQUET)'))
  changed <- receipt
  for (name in c('oracle.stdout.gz', 'pairs.parquet')) changed$sha256[[name]] <-
    duckvep_evidence_sha256(file.path(temporary, name))
  jsonlite::write_json(changed, path, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  for (name in names(changed$sha256)) stopifnot(identical(changed$sha256[[name]],
    duckvep_evidence_sha256(file.path(temporary, name))))
  checks['coordinated_rehashed_oracle_and_pairs'] <- rejected(check())
  if (!all(checks)) stop('accepted identity corruption: ', paste(names(checks)[!checks], collapse = ', '))
  stopifnot(!anyDuplicated(names(checks)))
  message(basename(directory), ': ', length(checks), ' retained-identity corruptions rejected')
}

# Reconstruct every bundle from raw VEP observations, not stored equality flags.
# The failed baseline and the corrected result retain the same finite experiment.
check_codon_bundle <- function(directory, hgvsp_failures, so_failures,
                               receipt_sha256, source_revision, extension_sha256) {
  codon_check_bundle_identity(directory, receipt_sha256, source_revision, extension_sha256)
  read_records <- function(name) {
    connection <- gzfile(file.path(directory, name), 'rt')
    on.exit(close(connection))
    readLines(connection)
  }
  source_records <- read_records('cases.jsonl.gz')
  oracle_records <- read_records('oracle.stdout.gz')
  cases <- lapply(source_records, read_unique_json)
  oracle <- lapply(oracle_records, read_unique_json)
  stopifnot(check_codon_matrix(cases),
    identical(vapply(cases, `[[`, '', 'id'), vapply(oracle, `[[`, '', 'id')))
  matrix_controls <- codon_matrix_controls(cases)
  stopifnot(nrow(controls) == 64L, length(matrix_controls) == 20L)
  message(basename(directory), ': all 84 comparator, provenance and matrix corruptions rejected')
  for (i in seq_along(cases)) {
    variants <- cases[[i]]$variants
    observed <- oracle[[i]]$independent_hgvs
    stopifnot(identical(oracle[[i]]$prepared_cds, cases[[i]]$cds),
      nrow(observed) == nrow(variants), !anyDuplicated(observed$id),
      setequal(observed$id, variants$id),
      identical(observed$allele, variants$alternate[match(observed$id, variants$id)]))
  }
  events <- do.call(rbind, lapply(cases, function(x) data.frame(
    event_index = as.integer(x$variants$id), position = x$variants$position1,
    reference = x$variants$reference, allele = x$variants$alternate, cds = x$cds, table = x$table,
    case_id = x$id, codon = substr(x$cds, 4L, 6L))))
  expected <- do.call(rbind, lapply(oracle, function(x) {
    rows <- x$independent_hgvs
    data.frame(event_index = as.integer(rows$id), allele = rows$allele,
      hgvsp = sub('^.*:p\\.', 'p.', rows$hgvsp),
      so = vapply(rows$consequences, function(terms) paste(sort(terms), collapse = '&'), ''))
  }))
  stopifnot(nrow(events) == 28800L, nrow(expected) == 28800L,
    !anyDuplicated(events$event_index), !anyDuplicated(expected$event_index))
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  pairs <- DBI::dbGetQuery(con, paste('SELECT * FROM read_parquet(',
    DBI::dbQuoteString(con, file.path(directory, 'pairs.parquet')), ')'))
  stopifnot(nrow(pairs) == 230400L, all(pairs$actual_present), all(pairs$expected_present),
    setequal(unique(pairs$route), paste0(rep(c('independent', 'strict', 'vep116_compat', 'source_records'),
      each = 2L), '_', c(1L, 4L))),
    all(pairs$seq_region == pairs$event_index - 1L),
    all(pairs$transcript_index == pairs$event_index - 1L))
  input_at <- match(pairs$event_index, events$event_index)
  for (field in c('position', 'reference', 'allele', 'cds', 'table', 'case_id', 'codon'))
    stopifnot(all(pairs[[field]] == events[[field]][input_at]))
  stopifnot(all(pairs$source_n == ifelse(pairs$reference == 'N', 'yes', 'no')),
    all(pairs$codon_n == ifelse(grepl('N', pairs$codon, fixed = TRUE), 'yes', 'no')))
  for (route in unique(pairs$route)) {
    part <- pairs[pairs$route == route, ]
    actual <- data.frame(event_index = part$event_index, allele = part$allele,
      hgvsp = part$hgvsp_actual, so = part$so_actual)
    target <- expected
    phased <- !startsWith(route, 'independent_')
    if (phased) target$so <- NA_character_
    checked <- codon_equal(actual, target)
    if (phased) checked$so_equal <- NA
    at <- match(checked$event_index, part$event_index)
    for (field in names(checked)) stopifnot(identical(checked[[field]], part[[field]][at]))
  }
  summary <- aggregate(cbind(pairs = rep(1L, nrow(pairs)), hgvsp_failures = !pairs$hgvsp_equal,
    so_compared = !is.na(pairs$so_equal), so_failures = !pairs$so_equal),
    pairs[c('route', 'source_n', 'codon_n')], sum, na.rm = TRUE)
  retained_summary <- read.csv(file.path(directory, 'summary.csv'))
  stopifnot(identical(names(summary), names(retained_summary)),
    all(summary == retained_summary), sum(!pairs$hgvsp_equal) == hgvsp_failures,
    sum(!pairs$so_equal, na.rm = TRUE) == so_failures,
    sum(!is.na(pairs$so_equal)) == 57600L)
  retained_controls <- read.csv(file.path(directory, 'controls.csv'))
  stopifnot(identical(controls$control, retained_controls$control),
    identical(controls$rejected, retained_controls$rejected))
  message(basename(directory), ': 230,400 HGVSp / 57,600 SO comparisons reconstructed; ',
    hgvsp_failures, ' HGVSp / ', so_failures, ' SO failures')
  codon_bundle_identity_controls(directory, receipt_sha256, source_revision,
    extension_sha256, pairs, oracle_records)
  list(source_records = source_records, oracle_records = oracle_records)
}

baseline <- check_codon_bundle('test/duckvep/conformance/data/ambiguous_codon_baseline',
  hgvsp_failures = 12992L, so_failures = 3248L,
  receipt_sha256 = 'a518fb720ed903da955d43f1fc871ca3619b9d970dd4b065bfa9aaa8d6792c6f',
  source_revision = '69abd2b98c31ea3ef39706cf185b7fef7dacecc9',
  extension_sha256 = '14220a70249129773ce12da2c4ebc87f93d234b461b72453d5f8d21509f1b7f6')
consensus <- check_codon_bundle('test/duckvep/conformance/data/ambiguous_codon_consensus',
  hgvsp_failures = 0L, so_failures = 0L,
  receipt_sha256 = '9a905193758657a7cf210352996c2f533e55a0d5d264caa2583bd01291b34bf1',
  source_revision = 'e3ec6d231cc7e5c769685761a4739c8fcffc30ff',
  extension_sha256 = '6230d58e6701a71f51127224f9282e0ca69785dad283644e8703496b659de70e')
stopifnot(identical(baseline$source_records, consensus$source_records),
  identical(baseline$oracle_records, consensus$oracle_records))
message('Baseline and consensus retain identical raw source cases and oracle observations')
