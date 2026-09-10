#!/usr/bin/env Rscript
# Network-free matrix and transport checks; synthetic rows are not VEP evidence.
source('test/duckvep/conformance/ambiguous_codon_differential.R')

read_codon_records <- function(directory, name) {
  connection <- gzfile(file.path(directory, name), 'rt')
  on.exit(close(connection))
  readLines(connection)
}

# These axes specify source records, not an annotation or normalization algorithm.
indel_source_matrix <- function() {
  bases <- c('A', 'C', 'G', 'T', 'N')
  codons <- character()
  for (third in bases) for (second in bases) for (first in bases)
    codons <- c(codons, paste0(first, second, third))
  tables <- c(1:6, 9:14, 16, 21:31)
  model_codons <- rep(codons, each = length(tables))
  model_tables <- rep(tables, times = length(codons))
  operations <- data.frame(position1 = rep(3:6, each = 14L),
    removed = rep(c(rep(0L, 7L), 1:3, 2L, 2L, 3L, 3L), 4L),
    inserted = rep(c('A', 'C', 'G', 'T', 'AC', 'GCC', 'ACGT',
      '', '', '', 'A', 'ACGT', 'A', 'ACGT'), 4L))
  cases <- vector('list', length(model_codons))
  events <- vector('list', length(model_codons))
  for (i in seq_along(model_codons)) {
    cds <- paste0('ATG', model_codons[i], 'GCCTAA')
    ids <- (i - 1L) * 56L + seq_len(56L)
    reference <- substring(cds, operations$position1,
      operations$position1 + operations$removed)
    alternate <- paste0(substring(cds, operations$position1, operations$position1),
      operations$inserted)
    case_id <- paste(model_tables[i], model_codons[i], sep = '/')
    cases[[i]] <- list(id = case_id, cds = cds, table = model_tables[i],
      variants = data.frame(id = as.character(ids), position1 = operations$position1,
        reference = reference, alternate = alternate))
    events[[i]] <- data.frame(event_index = ids, seq_region = ids - 1L,
      transcript_index = ids - 1L, position = operations$position1 + 10L,
      reference = reference, alternate = alternate, cds = cds,
      table = model_tables[i], codon = model_codons[i], case_id = case_id)
  }
  list(cases = cases, events = do.call(rbind, events))
}

check_indel_inputs <- function(actual, expected) {
  stopifnot(setequal(names(actual), c('cases', 'events')),
    length(actual$cases) == 3000L, nrow(actual$events) == 168000L,
    !anyDuplicated(actual$events$event_index),
    identical(actual$events, expected$events))
  for (i in seq_along(expected$cases)) {
    observed <- actual$cases[[i]]
    target <- expected$cases[[i]]
    stopifnot(!anyDuplicated(names(observed)),
      setequal(names(observed), c('id', 'cds', 'table', 'edits', 'variants',
        'variant_format', 'genomic_sequence', 'cds_start1')),
      identical(observed$id, target$id), identical(observed$cds, target$cds),
      identical(observed$table, target$table), identical(observed$edits, list()),
      identical(observed$variant_format, 'vcf'), identical(observed$cds_start1, 11L),
      identical(observed$genomic_sequence, paste0(strrep('A', 10L), target$cds,
        strrep('A', 10L))), length(observed$variants) == 56L)
    for (j in seq_len(56L)) {
      stopifnot(!anyDuplicated(names(observed$variants[[j]])),
        identical(observed$variants[[j]], as.list(target$variants[j, ])))
    }
  }
  invisible(TRUE)
}

indel_inputs <- codon_inputs('indel')
declared <- indel_source_matrix()
check_indel_inputs(indel_inputs, declared)
events <- indel_inputs$events
stopifnot(all(nchar(events$reference) != nchar(events$alternate)),
  all(substring(events$reference, 1L, 1L) == substring(events$alternate, 1L, 1L)),
  all(nchar(events$cds) == 12L),
  all(table(events$table, events$codon) == 56L),
  sum(startsWith(events$reference, 'N')) == 25200L,
  sum(grepl('N', substring(events$reference, 2L), fixed = TRUE)) == 19512L)
message('Indel source matrix: 3,000 models / 168,000 original events; ',
  '25,200 N-anchor and 19,512 N-containing removed spans retained')

rejected <- function(expression) tryCatch({ force(expression); FALSE }, error = function(e) TRUE)
matrix_controls <- c(missing_case = rejected(check_indel_inputs(
  list(cases = indel_inputs$cases[-1L], events = events), declared)),
  missing_event = rejected(check_indel_inputs(
    list(cases = indel_inputs$cases, events = events[-1L, ]), declared)))
for (field in c('id', 'cds', 'table', 'variant_format', 'genomic_sequence', 'cds_start1')) {
  changed <- indel_inputs
  value <- changed$cases[[1L]][[field]]
  changed$cases[[1L]][[field]] <- if (is.numeric(value)) value + 1L else paste0(value, 'X')
  matrix_controls[paste0('case_', field)] <- rejected(check_indel_inputs(changed, declared))
}
for (field in c('id', 'position1', 'reference', 'alternate')) {
  changed <- indel_inputs
  value <- changed$cases[[1L]]$variants[[1L]][[field]]
  changed$cases[[1L]]$variants[[1L]][[field]] <-
    if (is.numeric(value)) value + 1L else paste0(value, 'X')
  matrix_controls[paste0('variant_', field)] <- rejected(check_indel_inputs(changed, declared))
}
for (field in names(events)) {
  changed <- indel_inputs
  value <- changed$events[[field]][1L]
  changed$events[[field]][1L] <- if (is.numeric(value)) value + 1L else paste0(value, 'X')
  matrix_controls[paste0('event_', field)] <- rejected(check_indel_inputs(changed, declared))
}
changed <- indel_inputs
changed$cases[[2L]] <- changed$cases[[1L]]
matrix_controls['duplicate_case'] <- rejected(check_indel_inputs(changed, declared))
changed <- indel_inputs
changed$cases[[2L]]$table <- changed$cases[[1L]]$table
matrix_controls['count_preserving_model_swap'] <- rejected(check_indel_inputs(changed, declared))
changed <- indel_inputs
changed$cases[[1L]]$variants[[2L]] <- changed$cases[[1L]]$variants[[1L]]
matrix_controls['duplicate_variant'] <- rejected(check_indel_inputs(changed, declared))
changed$cases[[1L]]$variants[[2L]]$id <- indel_inputs$cases[[1L]]$variants[[2L]]$id
matrix_controls['count_preserving_variant_swap'] <- rejected(check_indel_inputs(changed, declared))
changed <- indel_inputs
changed$events[2L, ] <- changed$events[1L, ]
matrix_controls['duplicate_event'] <- rejected(check_indel_inputs(changed, declared))
stopifnot(all(matrix_controls), !anyDuplicated(names(matrix_controls)))
message('Indel source matrix: ', length(matrix_controls), ' corruptions rejected')

partitions <- codon_partitions(events, TRUE)
stopifnot(identical(partitions, data.frame(model_partition = 1:3,
  event_offset = c(0L, 56000L, 112000L), first_event_index = c(1L, 56001L, 112001L),
  last_event_index = c(56000L, 112000L, 168000L), source_events = rep(56000L, 3L))))
restored <- vector('list', 3L)
for (i in 1:3) {
  local <- codon_partition_events(events, partitions[i, , drop = FALSE])
  stopifnot(identical(local$seq_region, 0:55999),
    identical(local$transcript_index, 0:55999))
  local$seq_region <- local$seq_region + partitions$event_offset[i]
  local$transcript_index <- local$transcript_index + partitions$event_offset[i]
  restored[[i]] <- local
}
stopifnot(identical(do.call(rbind, restored), events))
partition_controls <- c(missing_event = rejected(codon_partitions(events[-1L, ], TRUE)),
  duplicate_event = rejected(codon_partitions(rbind(events, events[1L, ]), TRUE)))
for (field in c('event_index', 'seq_region', 'transcript_index')) {
  changed <- events
  changed[[field]][1L] <- changed[[field]][1L] + 1L
  partition_controls[field] <- rejected(codon_partitions(changed, TRUE))
}
for (field in c('event_offset', 'first_event_index', 'last_event_index', 'source_events')) {
  changed <- partitions[2L, , drop = FALSE]
  changed[[field]] <- changed[[field]] + 1L
  partition_controls[field] <- rejected(codon_partition_events(events, changed))
}
stopifnot(all(partition_controls))
message('Native model partitions: all 168,000 events reconstruct exactly; ',
  length(partition_controls), ' ordinal/partition corruptions rejected')

# Literal DNA "NA" is a string; a missing HGVSp is JSON null. Annotation values
# below are synthetic transport sentinels, not expected biological observations.
# GAA>GA at genomic position 13 and A/- at position 15 reconstruct the same DNA.
transport <- jsonlite::fromJSON('[
  {"id":"transport_a","prepared_cds":"ATGNAAGCCTAA","independent_hgvs":[
    {"id":"1","source_reference":"N","source_alternate":"NA",
     "parser_start":15,"parser_end":14,"parser_allele_string":"-/A","allele":"A",
     "hgvsp":"fixture:p.Test1","consequences":["test_z","test_a"]},
    {"id":"2","source_reference":"NA","source_alternate":"N",
     "parser_start":15,"parser_end":15,"parser_allele_string":"A/-","allele":"-",
     "hgvsp":null,"consequences":["test_missing"]}]},
  {"id":"transport_b","prepared_cds":"ATGAAAGCCTAA","independent_hgvs":[
    {"id":"3","source_reference":"A","source_alternate":"AC",
     "parser_start":15,"parser_end":14,"parser_allele_string":"-/C","allele":"C",
     "hgvsp":"fixture:p.Test3","consequences":["test_b"]},
    {"id":"4","source_reference":"GAA","source_alternate":"GA",
     "parser_start":15,"parser_end":15,"parser_allele_string":"A/-","allele":"-",
     "hgvsp":null,"consequences":["test_minimized"]}]}]', simplifyVector = FALSE)
transport_events <- data.frame(event_index = 1:4, position = c(14L, 14L, 14L, 13L),
  reference = c('N', 'NA', 'A', 'GAA'), alternate = c('NA', 'N', 'AC', 'GA'),
  case_id = c('transport_a', 'transport_a', 'transport_b', 'transport_b'),
  cds = c('ATGNAAGCCTAA', 'ATGNAAGCCTAA', 'ATGAAAGCCTAA', 'ATGAAAGCCTAA'))
transport_expected <- data.frame(event_index = 1:4, allele = c('NA', 'N', 'AC', 'GA'),
  hgvsp = c('p.Test1', NA_character_, 'p.Test3', NA_character_),
  so = c('test_a&test_z', 'test_missing', 'test_b', 'test_minimized'))
stopifnot(identical(codon_expected(transport, transport_events, TRUE), transport_expected))
transport_controls <- c(missing_case = rejected(codon_expected(transport[-1L], transport_events, TRUE)),
  duplicate_case = rejected(codon_expected(c(transport, transport[1L]), transport_events, TRUE)))
changed <- transport
changed[[1L]]$independent_hgvs <- changed[[1L]]$independent_hgvs[-1L]
transport_controls['missing_observation'] <- rejected(codon_expected(changed, transport_events, TRUE))
changed <- transport
changed[[1L]]$independent_hgvs[[2L]] <- changed[[1L]]$independent_hgvs[[1L]]
transport_controls['duplicate_observation'] <- rejected(codon_expected(changed, transport_events, TRUE))
changed[[1L]]$independent_hgvs[[2L]]$id <- '2'
transport_controls['count_preserving_observation_swap'] <-
  rejected(codon_expected(changed, transport_events, TRUE))
for (field in c('id', 'source_reference', 'source_alternate', 'parser_start', 'parser_end',
    'parser_allele_string', 'allele')) {
  changed <- transport
  value <- changed[[1L]]$independent_hgvs[[1L]][[field]]
  changed[[1L]]$independent_hgvs[[1L]][[field]] <-
    if (is.numeric(value)) value + 1L else paste0(value, 'X')
  transport_controls[paste0('observation_', field)] <-
    rejected(codon_expected(changed, transport_events, TRUE))
  changed[[1L]]$independent_hgvs[[1L]][[field]] <- NULL
  transport_controls[paste0('missing_', field)] <-
    rejected(codon_expected(changed, transport_events, TRUE))
}
changed <- transport
changed[[1L]]$id <- 'transport_b'
transport_controls['wrong_case_identity'] <- rejected(codon_expected(changed, transport_events, TRUE))
changed <- transport
changed[[2L]]$independent_hgvs[[2L]]$parser_allele_string <- 'A/A'
changed[[2L]]$independent_hgvs[[2L]]$allele <- 'A'
transport_controls['complex_minimized_sequence_change'] <-
  rejected(codon_expected(changed, transport_events, TRUE))
for (value in c('01', '1.0', '+1')) {
  changed <- transport
  changed[[1L]]$independent_hgvs[[1L]]$id <- value
  transport_controls[paste0('noncanonical_id_', value)] <-
    rejected(codon_expected(changed, transport_events, TRUE))
}
changed <- transport
changed[[1L]]$prepared_cds <- 'ATGAAAGCCTAA'
transport_controls['changed_prepared_cds'] <- rejected(codon_expected(changed, transport_events, TRUE))
changed[[1L]]$prepared_cds <- NULL
transport_controls['missing_prepared_cds'] <- rejected(codon_expected(changed, transport_events, TRUE))
for (field in names(transport_events)) {
  changed <- transport_events
  value <- changed[[field]][1L]
  changed[[field]][1L] <- if (is.numeric(value)) value + 1L else paste0(value, 'X')
  transport_controls[paste0('source_', field)] <- rejected(codon_expected(transport, changed, TRUE))
}
for (field in c('hgvsp', 'consequences')) {
  changed <- transport
  changed[[1L]]$independent_hgvs[[1L]][[field]] <- NULL
  transport_controls[paste0('missing_', field)] <-
    rejected(codon_expected(changed, transport_events, TRUE))
}
for (field in names(transport[[1L]]$independent_hgvs[[1L]])) {
  changed <- transport
  row <- changed[[1L]]$independent_hgvs[[1L]]
  changed[[1L]]$independent_hgvs[[1L]] <- c(row, row[field])
  stopifnot(anyDuplicated(names(changed[[1L]]$independent_hgvs[[1L]])) > 0L)
  transport_controls[paste0('duplicate_observation_key_', field)] <-
    rejected(codon_expected(changed, transport_events, TRUE))
}
stopifnot(!anyDuplicated(names(transport_controls)))
if (!all(transport_controls)) stop('Accepted observation corruption: ',
  paste(names(transport_controls)[!transport_controls], collapse = ', '))
comparison_controls <- rbind(codon_controls(transport_expected), codon_provenance_controls(FALSE),
  codon_provenance_controls(TRUE))
stopifnot(all(comparison_controls$rejected))
message('Indel observation transport: ', length(transport_controls),
  ' corruptions rejected; literal NA and JSON null remain distinct')
message('Shared pair/provenance comparator: ', nrow(comparison_controls), ' corruptions rejected')

snv <- codon_inputs('snv')
stopifnot(identical(codon_partitions(snv$events, FALSE), data.frame(model_partition = 1L,
  event_offset = 0L, first_event_index = 1L, last_event_index = 28800L, source_events = 28800L)))
source_records <- vapply(snv$cases, jsonlite::toJSON, '', auto_unbox = TRUE)
for (bundle in c('ambiguous_codon_baseline', 'ambiguous_codon_consensus')) {
  directory <- file.path('test/duckvep/conformance/data', bundle)
  stopifnot(identical(source_records, read_codon_records(directory, 'cases.jsonl.gz')))
  raw <- lapply(read_codon_records(directory, 'oracle.stdout.gz'), jsonlite::fromJSON,
    simplifyVector = FALSE)
  records <- unlist(lapply(raw, `[[`, 'independent_hgvs'), recursive = FALSE)
  expected <- data.frame(event_index = vapply(records, function(x) as.integer(x$id), 0L),
    allele = vapply(records, `[[`, '', 'allele'),
    hgvsp = vapply(records, function(x) if (is.null(x$hgvsp)) NA_character_ else
      sub('^.*:p\\.', 'p.', x$hgvsp), ''),
    so = vapply(records, function(x) paste(sort(unlist(x$consequences)), collapse = '&'), ''))
  stopifnot(nrow(expected) == 28800L,
    identical(codon_expected(raw, snv$events, FALSE), expected))
  message(bundle, ': original SNV cases byte-identical; 28,800 raw expectations reconstructed')
}
