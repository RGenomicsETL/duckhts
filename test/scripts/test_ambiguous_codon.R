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

# Reconstruct every bundle from raw VEP observations, not stored equality flags.
# The failed baseline and the corrected result retain the same finite experiment.
check_codon_bundle <- function(directory, hgvsp_failures, so_failures) {
  manifest <- read_unique_json(paste(readLines(file.path(directory, 'receipt.json')), collapse = '\n'))
  hashes <- unlist(manifest$sha256)
  stopifnot(setequal(names(hashes), c('pairs.parquet', 'cases.jsonl.gz', 'oracle.stdout.gz',
    'summary.csv', 'controls.csv', 'environment.stdout')))
  for (name in names(hashes)) stopifnot(identical(unname(hashes[name]),
    duckvep_evidence_sha256(file.path(directory, name))))
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
  list(source_records = source_records, oracle_records = oracle_records)
}

baseline <- check_codon_bundle('test/duckvep/conformance/data/ambiguous_codon_baseline',
  hgvsp_failures = 12992L, so_failures = 3248L)
consensus <- check_codon_bundle('test/duckvep/conformance/data/ambiguous_codon_consensus',
  hgvsp_failures = 0L, so_failures = 0L)
stopifnot(identical(baseline$source_records, consensus$source_records),
  identical(baseline$oracle_records, consensus$oracle_records))
message('Baseline and consensus retain identical raw source cases and oracle observations')
