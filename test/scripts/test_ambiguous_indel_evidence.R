#!/usr/bin/env Rscript
# Read-only audit: DIRECTORY RECEIPT_SHA256 pairs. Supplied digests identify
# diagnostic bundles; they do not authenticate execution or establish agreement.
# Adjacent bundles in the default pinned chain must retain source and oracle
# expectations with no new failure. Explicit DIRECTORY SHA pairs are independent audits.
# --incomplete-model-diagnostic permits explicitly labelled old bundles with
# unspecified native transcript flanks; the default requires complete flanks.
source('test/duckvep/conformance/ambiguous_codon_differential.R')

indel_audit_json <- function(text, simplify = FALSE) {
  value <- jsonlite::fromJSON(text, simplifyVector = FALSE)
  check <- function(x) {
    if (!is.list(x)) return(invisible(TRUE))
    stopifnot(!anyDuplicated(names(x)))
    for (child in x) check(child)
    invisible(TRUE)
  }
  check(value)
  if (simplify) jsonlite::fromJSON(text) else value
}

indel_audit_records <- function(directory, name) {
  connection <- gzfile(file.path(directory, name), 'rt')
  on.exit(close(connection))
  readLines(connection)
}

# Reconstruct source geometry from declared axes, independently of codon_inputs.
indel_audit_sources <- function(cases) {
  bases <- c('A', 'C', 'G', 'T', 'N')
  tables <- c(1:6, 9:14, 16, 21:31)
  stopifnot(length(cases) == 3000L)
  position1 <- rep(3:6, each = 14L)
  removed <- rep(c(rep(0L, 7L), 1:3, 2L, 2L, 3L, 3L), 4L)
  inserted <- rep(c('A', 'C', 'G', 'T', 'AC', 'GCC', 'ACGT',
    '', '', '', 'A', 'ACGT', 'A', 'ACGT'), 4L)
  events <- vector('list', 3000L)
  for (i in seq_len(3000L)) {
    codon_index <- (i - 1L) %/% 24L
    codon <- paste0(bases[codon_index %% 5L + 1L],
      bases[codon_index %/% 5L %% 5L + 1L], bases[codon_index %/% 25L + 1L])
    table <- tables[(i - 1L) %% 24L + 1L]
    cds <- paste0('ATG', codon, 'GCCTAA')
    id <- paste(table, codon, sep = '/')
    case <- cases[[i]]
    stopifnot(!anyDuplicated(names(case)),
      setequal(names(case), c('id', 'cds', 'table', 'edits', 'variants',
        'variant_format', 'genomic_sequence', 'cds_start1')),
      identical(case$id, id), identical(case$cds, cds), is.numeric(case$table),
      length(case$table) == 1L, !is.na(case$table), case$table == table,
      identical(case$edits, list()), identical(case$variant_format, 'vcf'),
      is.numeric(case$cds_start1), length(case$cds_start1) == 1L,
      !is.na(case$cds_start1), case$cds_start1 == 11L,
      identical(case$genomic_sequence, paste0(strrep('A', 10L), cds, strrep('A', 10L))),
      length(case$variants) == 56L)
    indices <- (i - 1L) * 56L + seq_len(56L)
    reference <- substring(cds, position1, position1 + removed)
    alternate <- paste0(substring(cds, position1, position1), inserted)
    for (j in seq_len(56L)) {
      row <- case$variants[[j]]
      stopifnot(!anyDuplicated(names(row)),
        setequal(names(row), c('id', 'position1', 'reference', 'alternate')),
        identical(row$id, as.character(indices[j])), is.numeric(row$position1),
        length(row$position1) == 1L, !is.na(row$position1),
        row$position1 == position1[j], identical(row$reference, reference[j]),
        identical(row$alternate, alternate[j]))
    }
    events[[i]] <- data.frame(event_index = indices, position = position1 + 10L,
      reference = reference, alternate = alternate, cds = cds, table = table,
      codon = codon, case_id = id)
  }
  do.call(rbind, events)
}

indel_audit_partitions <- function(partitions) {
  expected <- data.frame(model_partition = 1:3, event_offset = c(0L, 56000L, 112000L),
    first_event_index = c(1L, 56001L, 112001L), last_event_index = c(56000L, 112000L, 168000L),
    source_events = rep(56000L, 3L))
  stopifnot(is.data.frame(partitions), identical(partitions, expected))
  invisible(TRUE)
}

indel_audit_source_controls <- function(cases) {
  rejected <- function(x) tryCatch({ indel_audit_sources(x); FALSE }, error = function(e) TRUE)
  checks <- c(missing_case = rejected(cases[-1L]))
  changed <- cases
  changed[[2L]] <- changed[[1L]]
  checks['duplicate_case'] <- rejected(changed)
  for (field in c('id', 'cds', 'table', 'variant_format', 'genomic_sequence', 'cds_start1')) {
    changed <- cases
    value <- changed[[1L]][[field]]
    changed[[1L]][[field]] <- if (is.numeric(value)) value + 1L else paste0(value, 'X')
    checks[paste0('case_', field)] <- rejected(changed)
  }
  for (field in c('id', 'position1', 'reference', 'alternate')) {
    changed <- cases
    value <- changed[[1L]]$variants[[1L]][[field]]
    changed[[1L]]$variants[[1L]][[field]] <-
      if (is.numeric(value)) value + 1L else paste0(value, 'X')
    checks[paste0('variant_', field)] <- rejected(changed)
  }
  changed <- cases
  changed[[1L]]$variants[[2L]] <- changed[[1L]]$variants[[1L]]
  changed[[1L]]$variants[[2L]]$id <- cases[[1L]]$variants[[2L]]$id
  checks['count_preserving_variant_swap'] <- rejected(changed)
  changed <- cases
  changed[[1L]]$table <- rep(changed[[1L]]$table, 2L)
  checks['vector_table'] <- rejected(changed)
  changed <- cases
  changed[[1L]]$cds_start1 <- NULL
  checks['missing_cds_start'] <- rejected(changed)
  stopifnot(all(checks))
  message('Indel source audit: ', length(checks), ' source-geometry corruptions rejected')
  invisible(TRUE)
}

indel_audit_pairs <- function(pairs, events, expected) {
  routes <- paste0(rep(c('independent', 'strict', 'vep116_compat', 'source_records'), each = 2L),
    '_', c(1L, 4L))
  numeric_fields <- c('event_index', 'seq_region', 'transcript_index', 'position', 'table',
    'model_partition', 'event_offset')
  character_fields <- c('route', 'allele', 'reference', 'codon', 'cds', 'case_id', 'source_n', 'codon_n')
  nullable_text <- c('hgvsp_actual', 'hgvsp_expected', 'so_actual', 'so_expected')
  flags <- c('actual_present', 'expected_present', 'hgvsp_equal', 'so_equal')
  stopifnot(!anyDuplicated(names(pairs)),
    setequal(names(pairs), c(numeric_fields, character_fields, nullable_text, flags)),
    nrow(pairs) == 8L * nrow(events), setequal(pairs$route, routes),
    !anyDuplicated(events$event_index), !anyDuplicated(expected$event_index),
    setequal(events$event_index, expected$event_index))
  for (field in numeric_fields) stopifnot(is.numeric(pairs[[field]]), !anyNA(pairs[[field]]),
    all(is.finite(pairs[[field]])), all(pairs[[field]] == floor(pairs[[field]])))
  for (field in character_fields) stopifnot(is.character(pairs[[field]]), !anyNA(pairs[[field]]))
  for (field in nullable_text) stopifnot(is.character(pairs[[field]]))
  for (field in flags) stopifnot(is.logical(pairs[[field]]))
  stopifnot(all(pairs$expected_present), !anyNA(pairs$hgvsp_equal),
    all(is.na(pairs$actual_present) | pairs$actual_present))
  at <- match(pairs$event_index, events$event_index)
  stopifnot(!anyNA(at), identical(pairs$allele, events$alternate[at]))
  for (field in c('position', 'reference', 'codon', 'table', 'cds', 'case_id'))
    stopifnot(all(pairs[[field]] == events[[field]][at]))
  partition <- (pairs$event_index - 1L) %/% 56000L + 1L
  offset <- (partition - 1L) * 56000L
  stopifnot(all(pairs$model_partition == partition), all(pairs$event_offset == offset),
    all(pairs$seq_region == pairs$event_index - 1L - offset),
    all(pairs$transcript_index == pairs$seq_region),
    identical(pairs$source_n, ifelse(grepl('N', pairs$reference, fixed = TRUE), 'yes', 'no')),
    identical(pairs$codon_n, ifelse(grepl('N', pairs$codon, fixed = TRUE), 'yes', 'no')))
  for (route in routes) {
    part <- pairs[pairs$route == route, , drop = FALSE]
    stopifnot(nrow(part) == nrow(events), !anyDuplicated(part$event_index),
      setequal(part$event_index, events$event_index))
    present <- !is.na(part$actual_present)
    stopifnot(all(is.na(part$hgvsp_actual[!present])), all(is.na(part$so_actual[!present])))
    target <- expected
    phased <- !startsWith(route, 'independent_')
    if (phased) {
      stopifnot(all(present), all(is.na(part$so_actual)))
      target$so <- NA_character_
    }
    actual <- data.frame(event_index = part$event_index[present], allele = part$allele[present],
      hgvsp = part$hgvsp_actual[present], so = part$so_actual[present])
    checked <- codon_equal(actual, target)
    if (phased) checked$so_equal <- NA
    ordered <- part[match(checked$event_index, part$event_index), , drop = FALSE]
    for (field in names(checked)) stopifnot(identical(checked[[field]], ordered[[field]]))
  }
  aggregate(cbind(pairs = rep(1L, nrow(pairs)), hgvsp_failures = !pairs$hgvsp_equal,
    so_compared = !is.na(pairs$so_equal), so_failures = !pairs$so_equal),
    pairs[c('route', 'source_n', 'codon_n')], sum, na.rm = TRUE)
}

indel_audit_summary <- function(observed, expected) {
  keys <- c('route', 'source_n', 'codon_n')
  stopifnot(identical(names(observed), names(expected)), nrow(observed) == nrow(expected))
  ordered <- function(x) {
    stopifnot(!anyNA(x), !anyDuplicated(x[keys]))
    x <- x[do.call(order, x[keys]), , drop = FALSE]
    rownames(x) <- NULL
    x
  }
  stopifnot(identical(ordered(observed), ordered(expected)))
  invisible(TRUE)
}

indel_audit_comparison <- function(pairs, events, expected) {
  list(events = events, expected = expected,
    comparisons = pairs[c('route', 'event_index', 'allele', 'hgvsp_equal', 'so_equal')])
}

indel_audit_compare <- function(baseline, consensus) {
  ordered <- function(relation, keys) {
    stopifnot(is.data.frame(relation), !anyDuplicated(names(relation)),
      !anyNA(relation[keys]), !anyDuplicated(relation[keys]))
    relation <- relation[do.call(order, relation[keys]), , drop = FALSE]
    rownames(relation) <- NULL
    relation
  }
  stopifnot(identical(ordered(baseline$events, 'event_index'),
      ordered(consensus$events, 'event_index')),
    identical(ordered(baseline$expected, c('event_index', 'allele')),
      ordered(consensus$expected, c('event_index', 'allele'))))
  keys <- c('route', 'event_index', 'allele')
  before <- ordered(baseline$comparisons, keys)
  after <- ordered(consensus$comparisons, keys)
  stopifnot(identical(names(before), c(keys, 'hgvsp_equal', 'so_equal')),
    identical(names(after), names(before)), identical(before[keys], after[keys]),
    !anyNA(before$hgvsp_equal), !anyNA(after$hgvsp_equal),
    identical(is.na(before$so_equal), is.na(after$so_equal)))
  report <- lapply(c('hgvsp', 'so'), function(metric) {
    field <- paste0(metric, '_equal')
    stopifnot(is.logical(before[[field]]), is.logical(after[[field]]))
    compared <- !is.na(before[[field]])
    prior_failures <- compared & !before[[field]]
    current_failures <- compared & !after[[field]]
    newly_failed <- which(current_failures & !prior_failures)
    if (length(newly_failed)) {
      first <- after[newly_failed[1L], keys]
      stop('New ', metric, ' failures: ', length(newly_failed), '; first route/event/ALT: ',
        paste(first, collapse = '/'), call. = FALSE)
    }
    data.frame(metric = metric, comparisons = sum(compared),
      baseline_failures = sum(prior_failures), consensus_failures = sum(current_failures),
      resolved_failures = sum(prior_failures & !current_failures), new_failures = 0L)
  })
  report <- do.call(rbind, report)
  message('Indel cross-bundle audit: identical source events, oracle expectations and route/event/ALT keys; ',
    paste(report$metric, report$comparisons, 'comparisons,', report$resolved_failures,
      'resolved failures, 0 new failures', collapse = '; '))
  invisible(report)
}

indel_audit_receipt <- function(receipt, allow_incomplete_model = FALSE) {
  fields <- c('source_revision', 'source_binding', 'extension_path', 'extension_sha256',
    'oracle_revisions', 'oracle_pairs', 'decoded_gt', 'source_gt', 'scope', 'sha256',
    'oracle_models', 'native_transcript_instances', 'native_partitions', 'native_ordinal_contract',
    'hgvsp_comparisons', 'so_comparisons', 'source_indels', 'variant_family',
    'local_receipt_sha256', 'source_sha256')
  has_flank_contract <- 'native_transcript_flanks' %in% names(receipt)
  stopifnot(!anyDuplicated(names(receipt)),
    setequal(names(receipt), c(fields, if (has_flank_contract) 'native_transcript_flanks')),
    identical(receipt$source_binding, 'diagnostic_unbound'), identical(receipt$variant_family, 'indel'),
    identical(receipt$oracle_models, 3000L), identical(receipt$source_indels, 168000L),
    identical(receipt$native_transcript_instances, 168000L), identical(receipt$oracle_pairs, 168000L),
    identical(receipt$hgvsp_comparisons, 1344000L), identical(receipt$so_comparisons, 336000L),
    identical(receipt$decoded_gt, 'haploid ALT'), identical(receipt$source_gt, '1|1'),
    identical(receipt$oracle_revisions, c('57ea5c52340acc1f156267f810ad162e26597082',
      '2fb834b987ede3824e200197a838ce11e91aeb4b')),
    identical(receipt$scope, paste0('independent_SO_HGVSp_and_singleton_phased_HGVSp_',
      'original_VCF_indels_forward_single_exon;complete_pairs_and_oracle_not_all_native_output_fields')),
    identical(receipt$native_ordinal_contract, paste(
      'seq_region and transcript_index are local to model_partition;',
      'event_index equals local ordinal plus event_offset plus one')))
  if (has_flank_contract) {
    stopifnot(identical(receipt$native_transcript_flanks,
      'complete_empty_pre_CDS_and_post_CDS_matching_oracle_transcript'))
  } else if (!allow_incomplete_model) {
    stop('Native transcript flanks are unspecified; use --incomplete-model-diagnostic only for retained diagnostics')
  }
  digest <- function(x, width) is.character(x) && length(x) == 1L && !is.na(x) &&
    grepl(paste0('^[0-9a-f]{', width, '}$'), x)
  stopifnot(digest(receipt$source_revision, 40L), digest(receipt$extension_sha256, 64L),
    digest(receipt$local_receipt_sha256, 64L), is.character(receipt$extension_path),
    length(receipt$extension_path) == 1L, !is.na(receipt$extension_path), nzchar(receipt$extension_path))
  indel_audit_partitions(receipt$native_partitions)
  for (field in c('sha256', 'source_sha256')) {
    hashes <- unlist(receipt[[field]])
    stopifnot(is.character(hashes), length(hashes) > 0L, !is.null(names(hashes)),
      !anyDuplicated(names(hashes)), !anyNA(hashes), all(grepl('^[0-9a-f]{64}$', hashes)),
      all(vapply(receipt[[field]], digest, TRUE, width = 64L)))
  }
  stopifnot(setequal(names(receipt$sha256), c('cases.jsonl.gz', 'oracle.stdout.gz', 'oracle.stderr.gz',
    'pairs.parquet', 'controls.csv', 'summary.csv', 'environment.stdout')),
    all(c('test/duckvep/conformance/ambiguous_codon_differential.R',
      'test/duckvep/conformance/reference_translation_oracle.pl', 'scripts/duckvep_evidence.R',
      'src/duckvep/kernel/src/duckvep_delta.c') %in% names(receipt$source_sha256)))
  invisible(TRUE)
}

indel_audit_receipt_controls <- function(receipt) {
  # This copy tests metadata validation only; it is not a new execution receipt.
  typed <- receipt
  typed$native_transcript_flanks <- 'complete_empty_pre_CDS_and_post_CDS_matching_oracle_transcript'
  rejected <- function(x) tryCatch({ indel_audit_receipt(x); FALSE }, error = function(e) TRUE)
  stopifnot(!rejected(typed))
  checks <- logical()
  for (field in names(typed)) {
    changed <- typed
    changed[[field]] <- NULL
    checks[paste0('missing_', field)] <- rejected(changed)
  }
  for (field in c('oracle_models', 'source_indels', 'native_transcript_instances',
      'oracle_pairs', 'hgvsp_comparisons', 'so_comparisons')) {
    changed <- typed
    changed[[field]] <- changed[[field]] + 0.5
    checks[paste0('fractional_', field)] <- rejected(changed)
  }
  for (field in names(typed$native_partitions)) {
    changed <- typed
    changed$native_partitions[[field]][2L] <- changed$native_partitions[[field]][2L] + 1L
    checks[paste0('partition_', field)] <- rejected(changed)
  }
  changed <- typed
  changed$native_transcript_flanks <- 'unknown'
  checks['invalid_flank_contract'] <- rejected(changed)
  changed$native_transcript_flanks <- NULL
  stopifnot(isTRUE(indel_audit_receipt(changed, allow_incomplete_model = TRUE)))
  stopifnot(all(checks))
  message('Indel receipt audit: ', length(checks), ' metadata corruptions rejected; ',
    'unspecified flanks require explicit diagnostic mode')
  invisible(TRUE)
}

indel_audit_bundle <- function(directory, receipt_sha256, allow_incomplete_model = FALSE) {
  directory <- normalizePath(directory, mustWork = TRUE)
  receipt_path <- file.path(directory, 'receipt.json')
  duckvep_evidence_check_receipt_pin(receipt_path, receipt_sha256)
  receipt <- indel_audit_json(paste(readLines(receipt_path), collapse = '\n'), simplify = TRUE)
  indel_audit_receipt(receipt, allow_incomplete_model)
  indel_audit_receipt_controls(receipt)
  stopifnot(setequal(list.files(directory, all.files = TRUE, no.. = TRUE),
    c(names(receipt$sha256), 'receipt.json')))
  for (name in names(receipt$sha256)) stopifnot(identical(receipt$sha256[[name]],
    duckvep_evidence_sha256(file.path(directory, name))))
  stopifnot(identical(duckvep_evidence_explicit_packages(readLines(file.path(directory, 'environment.stdout'))),
    duckvep_evidence_explicit_packages(readLines(
      'test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt'))))
  cases <- lapply(indel_audit_records(directory, 'cases.jsonl.gz'), indel_audit_json)
  events <- indel_audit_sources(cases)
  indel_audit_source_controls(cases)
  oracle <- lapply(indel_audit_records(directory, 'oracle.stdout.gz'), indel_audit_json)
  stopifnot(identical(vapply(cases, `[[`, '', 'id'), vapply(oracle, `[[`, '', 'id')))
  expected <- codon_expected(oracle, events, TRUE)
  stopifnot(nrow(expected) == 168000L)
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  pairs <- DBI::dbGetQuery(con, paste('SELECT * FROM read_parquet(',
    DBI::dbQuoteString(con, file.path(directory, 'pairs.parquet')), ')'))
  summary <- indel_audit_pairs(pairs, events, expected)
  indel_audit_summary(read.csv(file.path(directory, 'summary.csv')), summary)
  seed <- data.frame(event_index = 1:2, allele = c('A', 'C'),
    hgvsp = c('p.Test1', NA_character_), so = c('test_one', 'test_two'))
  controls <- rbind(codon_controls(seed), codon_provenance_controls(FALSE), codon_provenance_controls(TRUE))
  stopifnot(identical(read.csv(file.path(directory, 'controls.csv')), controls))
  warnings <- indel_audit_records(directory, 'oracle.stderr.gz')
  model_scope <- if ('native_transcript_flanks' %in% names(receipt)) 'complete native flank metadata' else
    'INCOMPLETE-MODEL DIAGNOSTIC: native transcript flanks unspecified'
  message(directory, ': ', model_scope, '; diagnostic receipt ', receipt_sha256, '; ',
    '1,344,000 HGVSp / 336,000 SO comparisons reconstructed; ',
    sum(summary$hgvsp_failures), ' HGVSp / ', sum(summary$so_failures), ' SO failures; ',
    length(warnings), ' retained oracle diagnostic lines')
  invisible(indel_audit_comparison(pairs, events, expected))
}

indel_audit_comparison_controls <- function(pairs, events, expected) {
  rejected <- function(expression) tryCatch({ force(expression); FALSE }, error = function(e) TRUE)
  third <- pairs$event_index == events$event_index[3L]
  pairs$actual_present[third] <- TRUE
  pairs$hgvsp_actual[third] <- expected$hgvsp[3L]
  pairs$hgvsp_equal[third] <- TRUE
  checked_bundle <- function(pair_rows, target = expected) {
    indel_audit_pairs(pair_rows, events, target)
    indel_audit_comparison(pair_rows, events, target)
  }
  baseline <- checked_bundle(pairs)
  baseline_summary <- indel_audit_pairs(pairs, events, expected)
  unchanged <- indel_audit_compare(baseline, baseline)
  stopifnot(all(unchanged$baseline_failures > 0L), all(unchanged$resolved_failures == 0L))
  reordered <- lapply(baseline, function(x) x[rev(seq_len(nrow(x))), , drop = FALSE])
  stopifnot(identical(unchanged, indel_audit_compare(baseline, reordered)))
  changed <- baseline
  changed$comparisons <- changed$comparisons[-1L, ]
  checks <- c(missing_key = rejected(indel_audit_compare(baseline, changed)))
  changed <- baseline
  changed$comparisons <- rbind(changed$comparisons, changed$comparisons[1L, ])
  checks['duplicate_key'] <- rejected(indel_audit_compare(baseline, changed))
  changed <- baseline
  changed$comparisons$event_index[1:2] <- rev(changed$comparisons$event_index[1:2])
  checks['count_preserving_key_swap'] <- rejected(indel_audit_compare(baseline, changed))
  for (metric in c('hgvsp', 'so')) {
    changed <- pairs
    if (metric == 'hgvsp') {
      changed$hgvsp_actual[2:3] <- c(NA_character_, 'p.ExchangedFailure')
      changed$hgvsp_equal[2:3] <- c(TRUE, FALSE)
    } else {
      changed$so_actual[2:3] <- c('exchanged_failure', expected$so[3L])
      changed$so_equal[2:3] <- c(FALSE, TRUE)
    }
    candidate <- checked_bundle(changed)
    indel_audit_summary(indel_audit_pairs(changed, events, expected), baseline_summary)
    field <- paste0(metric, '_equal')
    stopifnot(sum(!baseline$comparisons[[field]], na.rm = TRUE) ==
      sum(!candidate$comparisons[[field]], na.rm = TRUE))
    checks[paste0('count_preserving_', metric, '_failure_exchange')] <-
      rejected(indel_audit_compare(baseline, candidate))
    changed <- pairs
    target <- expected
    event <- if (metric == 'hgvsp') 2L else 3L
    target[[metric]][event] <- 'changed_oracle_expectation'
    selected <- changed$event_index == events$event_index[event]
    if (metric == 'so') selected <- selected & startsWith(changed$route, 'independent_')
    changed[[paste0(metric, '_expected')]][selected] <- target[[metric]][event]
    candidate <- checked_bundle(changed, target)
    checks[paste0('coordinated_', metric, '_expectation_change')] <-
      rejected(indel_audit_compare(baseline, candidate))
  }
  changed <- baseline
  changed$comparisons$so_equal[1L] <- NA
  checks['lost_so_comparison'] <- rejected(indel_audit_compare(baseline, changed))
  changed <- baseline
  changed$events$reference[1L] <- 'C'
  checks['changed_source_geometry'] <- rejected(indel_audit_compare(baseline, changed))
  changed <- pairs
  changed$hgvsp_actual[2L] <- NA_character_
  changed$hgvsp_equal[2L] <- TRUE
  improved <- indel_audit_compare(baseline, checked_bundle(changed))
  stopifnot(improved$resolved_failures[improved$metric == 'hgvsp'] == 1L)
  if (!all(checks)) stop('Accepted cross-bundle corruption: ', paste(names(checks)[!checks], collapse = ', '))
  message('Indel cross-bundle controls: ', length(checks), ' corruptions rejected; ',
    'unchanged nonzero failures, reordered rows and a genuine resolved failure accepted')
  invisible(TRUE)
}

indel_audit_self_test <- function() {
  rejected <- function(expression) tryCatch({ force(expression); FALSE }, error = function(e) TRUE)
  events <- data.frame(event_index = c(1L, 56001L, 112001L), position = 14L,
    reference = c('N', 'A', 'A'), alternate = c('NA', 'AC', 'AT'),
    cds = c('ATGNAAGCCTAA', 'ATGAAAGCCTAA', 'ATGAAAGCCTAA'),
    table = 1, codon = c('NAA', 'AAA', 'AAA'), case_id = c('a', 'b', 'c'))
  expected <- data.frame(event_index = events$event_index, allele = events$alternate,
    hgvsp = c('p.Test1', NA_character_, 'p.Test3'), so = c('test_a', 'test_b', 'test_c'))
  routes <- paste0(rep(c('independent', 'strict', 'vep116_compat', 'source_records'), each = 2L),
    '_', c(1L, 4L))
  pairs <- do.call(rbind, lapply(routes, function(route) {
    phased <- !startsWith(route, 'independent_')
    data.frame(event_index = events$event_index, route = route, allele = events$alternate,
      hgvsp_actual = c('p.Test1', 'p.Wrong2', NA_character_),
      so_actual = if (phased) rep(NA_character_, 3L) else c('test_a', 'test_b', NA_character_),
      actual_present = if (phased) rep(TRUE, 3L) else c(TRUE, TRUE, NA),
      hgvsp_expected = expected$hgvsp,
      so_expected = if (phased) rep(NA_character_, 3L) else expected$so,
      expected_present = TRUE, hgvsp_equal = c(TRUE, FALSE, FALSE),
      so_equal = if (phased) rep(NA, 3L) else c(TRUE, TRUE, FALSE),
      seq_region = 0L, transcript_index = 0L, position = events$position,
      reference = events$reference, codon = events$codon, table = events$table,
      cds = events$cds, case_id = events$case_id, model_partition = 1:3,
      event_offset = c(0L, 56000L, 112000L), source_n = c('yes', 'no', 'no'),
      codon_n = c('yes', 'no', 'no'))
  }))
  summary <- indel_audit_pairs(pairs, events, expected)
  stopifnot(sum(summary$hgvsp_failures) == 16L, sum(summary$so_failures) == 2L)
  checks <- c(missing_pair = rejected(indel_audit_pairs(pairs[-1L, ], events, expected)),
    extra_pair = rejected(indel_audit_pairs(rbind(pairs, pairs[1L, ]), events, expected)))
  for (field in names(pairs)) {
    changed <- pairs
    value <- changed[[field]][1L]
    changed[[field]][1L] <- if (is.logical(value)) !value else
      if (is.numeric(value)) value + 1L else paste0(value, 'X')
    checks[paste0('changed_', field)] <- rejected(indel_audit_pairs(changed, events, expected))
    changed <- pairs
    changed[[field]][1L] <- NA
    checks[paste0('missing_', field)] <- rejected(indel_audit_pairs(changed, events, expected))
  }
  changed <- pairs
  changed[2L, ] <- changed[1L, ]
  checks['count_preserving_duplicate'] <- rejected(indel_audit_pairs(changed, events, expected))
  changed <- pairs
  changed$so_equal[startsWith(changed$route, 'strict_')] <- TRUE
  checks['invented_phased_so'] <- rejected(indel_audit_pairs(changed, events, expected))
  changed <- summary
  changed$hgvsp_failures[1L] <- changed$hgvsp_failures[1L] + 1L
  checks['forged_summary'] <- rejected(indel_audit_summary(changed, summary))
  checks['duplicate_json_key'] <- rejected(indel_audit_json('{"id":"a","id":"b"}'))
  if (!all(checks)) stop('Accepted synthetic corruption: ', paste(names(checks)[!checks], collapse = ', '))
  message('Indel evidence audit: ', length(checks), ' in-memory corruptions rejected; ',
    'nonzero disagreements and explicit independent-native absence retained')
  indel_audit_comparison_controls(pairs, events, expected)
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  indel_audit_self_test()
  args <- commandArgs(trailingOnly = TRUE)
  allow_incomplete_model <- length(args) > 0L && args[1L] == '--incomplete-model-diagnostic'
  if (allow_incomplete_model) args <- args[-1L]
  stopifnot(length(args) %% 2L == 0L)
  compare_bundles <- !length(args)
  if (compare_bundles) {
    stopifnot(!allow_incomplete_model)
    args <- c(
      'test/duckvep/conformance/data/ambiguous_indel_baseline',
      '089860f4a90e07b33d73373f0683c9dbb80d8a3d5899f689517f092842037b03',
      'test/duckvep/conformance/data/ambiguous_indel_consensus',
      'a1d0f0cd56717f20cfa663de72efa22ab16a416a460d65fa6ffa7d6122039385',
      'test/duckvep/conformance/data/ambiguous_indel_translation',
      '263fe1773a382fdb76b27e034a3a126a069c0282f5afdc992d95e5a4762cc053')
  }
  baseline <- NULL
  for (i in seq.int(1L, length(args), by = 2L)) {
    audited <- indel_audit_bundle(args[i], args[i + 1L], allow_incomplete_model)
    if (compare_bundles) {
      if (!is.null(baseline)) indel_audit_compare(baseline, audited)
      baseline <- audited
    }
  }
}
