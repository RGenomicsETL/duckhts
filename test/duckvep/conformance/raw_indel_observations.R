# Fixed GT 1|1, forward single-exon container observations. These are not a
# raw-parser or HGVS oracle. Original records and complete containers stay on disk.
raw_indel_object <- function(x, required) {
  stopifnot(is.list(x), !anyDuplicated(names(x)), all(required %in% names(x)))
  invisible(TRUE)
}

raw_indel_scalar <- function(x, name, type = 'character') {
  raw_indel_object(x, name)
  value <- x[[name]]
  if (is.null(value)) return(if (type == 'character') NA_character_ else NA_real_)
  stopifnot(length(value) == 1L, !is.na(value), typeof(value) %in%
    if (type == 'character') 'character' else c('integer', 'double'))
  value
}

raw_indel_ids <- function(x) {
  if (is.null(x)) return(NA_character_)
  values <- unlist(x, use.names = FALSE)
  if (is.numeric(values)) {
    stopifnot(all(is.finite(values)), all(values > 0), all(values == floor(values)),
      all(values <= 2^53 - 1))
    values <- format(values, scientific = FALSE, trim = TRUE, digits = 22L)
  } else values <- as.character(values)
  stopifnot(!anyNA(values), !anyDuplicated(values), all(grepl('^[1-9][0-9]*$', values)))
  as.character(jsonlite::toJSON(sort(values), auto_unbox = FALSE))
}

raw_indel_masks <- function() {
  lines <- readLines('src/duckvep/kernel/src/duckvep_haplotype.h')
  vapply(c(indel = 'INDEL', frameshift = 'FRAMESHIFT'), function(name) {
    pattern <- paste0('^ *DUCKVEP_HAPLOTYPE_FLAG_', name, ' += 1u << ([0-9]+), *$')
    definition <- grep(pattern, lines, value = TRUE)
    stopifnot(length(definition) == 1L)
    bitwShiftL(1L, as.integer(sub(pattern, '\\1', definition)))
  }, 1L)
}

raw_indel_case <- function(case, events) {
  raw_indel_object(case, c('id', 'raw_haplotypes'))
  stopifnot(identical(case$id, unique(events$case_id)),
    length(case$raw_haplotypes) == nrow(events))
  ids <- vapply(case$raw_haplotypes, raw_indel_scalar, '', name = 'id')
  stopifnot(!anyNA(ids), !anyDuplicated(ids), setequal(ids, as.character(events$event_index)))
  events <- events[match(ids, as.character(events$event_index)), , drop = FALSE]
  rows <- vector('list', length(ids))
  for (i in seq_along(ids)) {
    source <- events[i, , drop = FALSE]
    observed <- case$raw_haplotypes[[i]]
    raw_indel_object(observed, c('id', 'source_reference', 'source_alternate',
      'source_position1', 'raw_gt', 'mapping_start', 'mapping_end', 'reference_cds',
      'reference_protein', 'lanes', 'container'))
    stopifnot(identical(raw_indel_scalar(observed, 'source_reference'), source$reference),
      identical(raw_indel_scalar(observed, 'source_alternate'), source$alternate),
      raw_indel_scalar(observed, 'source_position1', 'numeric') == source$position,
      identical(observed$raw_gt, '1|1'), identical(observed$reference_cds, source$cds),
      isTRUE(raw_indel_scalar(observed, 'mapping_start', 'numeric') == source$position - 10L),
      isTRUE(raw_indel_scalar(observed, 'mapping_end', 'numeric') ==
        source$position - 11L + nchar(source$reference)),
      length(observed$lanes) == 2L)
    raw_indel_scalar(observed, 'reference_protein')
    lanes <- observed$lanes
    stopifnot(identical(sort(vapply(lanes, raw_indel_scalar, 0,
      name = 'lane1', type = 'numeric')), c(1, 2)))
    container <- observed$container
    raw_indel_object(container, c('transcript_id', 'total_haplotype_count',
      'total_population_counts', 'cds_haplotypes', 'protein_haplotypes'))
    stopifnot(identical(container$transcript_id, case$id),
      isTRUE(raw_indel_scalar(container, 'total_haplotype_count', 'numeric') == 2L),
      identical(names(container$total_population_counts), '_all'),
      isTRUE(raw_indel_scalar(container$total_population_counts, '_all', 'numeric') == 2L))
    lane_rows <- lapply(lanes, function(lane) {
      raw_indel_object(lane, c('lane1', 'sample', 'cds', 'protein', 'flags', 'applied_sources'))
      stopifnot(identical(lane$sample, 'translation_sample'))
      raw_indel_object(lane$flags, character())
      stopifnot(all(names(lane$flags) %in% c('indel', 'frameshift', 'length_diff')))
      # TranscriptHaplotypeContainer constructs these fields with `$flag || 0`.
      flag <- function(name) if (name %in% names(lane$flags))
        raw_indel_scalar(lane$flags, name, 'numeric') else 0
      applied <- raw_indel_ids(lane$applied_sources)
      stopifnot(!is.na(applied), all(unlist(lane$applied_sources) %in% ids[i]))
      edit_count <- length(lane$applied_sources)
      for (axis in c('cds', 'protein')) {
        sequence <- raw_indel_scalar(lane, axis)
        groups <- container[[paste0(axis, '_haplotypes')]]
        stopifnot(length(groups) == 1L)
        matching <- vapply(groups, function(group) identical(group$seq, sequence), TRUE)
        stopifnot(sum(matching) == 1L)
        group <- groups[[which(matching)]]
        raw_indel_object(group, c('seq', 'count', 'samples', 'has_indel', 'contributing_variants'))
        stopifnot(isTRUE(raw_indel_scalar(group, 'count', 'numeric') == 2L),
          identical(names(group$samples), 'translation_sample'),
          isTRUE(raw_indel_scalar(group$samples, 'translation_sample', 'numeric') == 2L))
        # ProteinHaplotype filters contributors through TVA peptide effects.
        # Physical edit provenance is the observed CDS lane's applied_sources.
        if (axis == 'cds') stopifnot(
          isTRUE(raw_indel_scalar(group, 'has_indel', 'numeric') == flag('indel')),
          identical(raw_indel_ids(group$contributing_variants), applied))
      }
      data.frame(event_index = source$event_index, allele = source$alternate,
        lane = lane$lane1, model_partition = source$model_partition,
        event_offset = source$event_offset, transcript_index = source$transcript_index,
        seq_region = source$seq_region, position = observed$source_position1,
        reference = observed$source_reference, alt_index = 1L, sample_index = 0L,
        phase_set = NA_real_, ploidy = length(lanes), carrier_count = container$total_haplotype_count,
        source_ids = raw_indel_ids(list(observed$id)), applied_source_ids = applied,
        edit_count = edit_count, block_count = edit_count,
        indel = flag('indel'), frameshift = flag('frameshift'),
        nominal_length_diff = flag('length_diff'),
        block_cds_start = if (edit_count) observed$mapping_start else NA_real_,
        block_reference = if (edit_count) observed$source_reference else NA_character_,
        block_alternate = if (edit_count) observed$source_alternate else NA_character_,
        block_length_diff = flag('length_diff'), block_indel = flag('indel'),
        block_frameshift = flag('frameshift'),
        cds = raw_indel_scalar(lane, 'cds'), protein = raw_indel_scalar(lane, 'protein'))
    })
    rows[[i]] <- do.call(rbind, lane_rows)
  }
  do.call(rbind, rows)
}

raw_indel_expected <- function(path, events) {
  input <- file(path, 'rt')
  on.exit(close(input))
  case_ids <- unique(events$case_id)
  by_case <- split(events, events$case_id)
  rows <- vector('list', length(case_ids))
  seen <- character()
  repeat {
    line <- readLines(input, n = 1L)
    if (!length(line)) break
    case <- jsonlite::fromJSON(line, simplifyVector = FALSE)
    id <- raw_indel_scalar(case, 'id')
    stopifnot(!is.na(id), id %in% case_ids, !id %in% seen)
    seen <- c(seen, id)
    rows[[length(seen)]] <- raw_indel_case(case, by_case[[id]])
  }
  stopifnot(setequal(seen, case_ids))
  do.call(rbind, rows)
}

raw_indel_native <- function(actual, partition, masks) {
  required <- c('transcript_index', 'cds', 'protein', 'sequence_flags', 'edit_count',
    'carrier_count', 'carriers', 'contributors', 'coding_blocks', 'nominal_length_diff')
  stopifnot(all(required %in% names(actual)))
  rows <- vector('list', nrow(actual))
  for (i in seq_len(nrow(actual))) {
    source <- actual$contributors[[i]]
    stopifnot(is.data.frame(source), nrow(source) == 1L,
      all(c('event_index', 'seq_region', 'position', 'reference', 'alternate', 'alt_index') %in%
        names(source)))
    carriers <- actual$carriers[[i]]
    stopifnot(is.data.frame(carriers),
      all(c('sample_index', 'phase_set', 'haplotype_lane', 'ploidy') %in% names(carriers)))
    if (!nrow(carriers)) carriers[1L, ] <- NA
    blocks <- actual$coding_blocks[[i]]
    block_count <- if (is.null(blocks)) NA_integer_ else nrow(blocks)
    if (!is.na(block_count)) stopifnot(all(c('cds_start', 'reference', 'alternate',
      'length_change', 'sequence_flags', 'event_indices') %in% names(blocks)))
    # The declared matrix has one physical replacement or an observed skipped source.
    # Multiple blocks are retained as a failed block-count comparison.
    single <- !is.na(block_count) && block_count == 1L
    block_flags <- if (is.na(block_count)) NA_integer_ else
      Reduce(bitwOr, as.integer(blocks$sequence_flags), init = 0L)
    rows[[i]] <- data.frame(event_index = source$event_index, allele = source$alternate,
      lane = carriers$haplotype_lane, model_partition = partition$model_partition,
      event_offset = partition$event_offset, transcript_index = actual$transcript_index[i],
      seq_region = source$seq_region, position = source$position, reference = source$reference,
      alt_index = source$alt_index, sample_index = carriers$sample_index,
      phase_set = carriers$phase_set, ploidy = carriers$ploidy,
      carrier_count = actual$carrier_count[i], source_ids = raw_indel_ids(as.list(source$event_index)),
      applied_source_ids = if (is.na(block_count)) NA_character_ else
        raw_indel_ids(as.list(unlist(blocks$event_indices))), edit_count = actual$edit_count[i],
      block_count = block_count,
      indel = as.integer(bitwAnd(as.integer(actual$sequence_flags[i]), masks[['indel']]) != 0L),
      frameshift = as.integer(bitwAnd(as.integer(actual$sequence_flags[i]), masks[['frameshift']]) != 0L),
      nominal_length_diff = actual$nominal_length_diff[i],
      block_cds_start = if (single) blocks$cds_start else NA_real_,
      block_reference = if (single) blocks$reference else NA_character_,
      block_alternate = if (single) blocks$alternate else NA_character_,
      block_length_diff = if (is.na(block_count)) NA_real_ else sum(blocks$length_change),
      block_indel = as.integer(bitwAnd(block_flags, masks[['indel']]) != 0L),
      block_frameshift = as.integer(bitwAnd(block_flags, masks[['frameshift']]) != 0L),
      cds = actual$cds[i], protein = actual$protein[i])
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

raw_indel_equal <- function(actual, expected) {
  keys <- c('event_index', 'allele', 'lane')
  stopifnot(identical(names(actual), names(expected)))
  key <- function(x) do.call(paste, c(x[keys], sep = '/'))
  stopifnot(!anyDuplicated(key(actual)), !anyDuplicated(key(expected)),
    !anyNA(actual[c('event_index', 'allele')]), !anyNA(expected[keys]))
  fields <- setdiff(names(expected), keys)
  actual$actual_present <- TRUE
  expected$expected_present <- TRUE
  pairs <- merge(actual, expected, by = keys, all = TRUE,
    suffixes = c('_actual', '_expected'), sort = TRUE)
  present <- !is.na(pairs$actual_present) & !is.na(pairs$expected_present)
  for (field in fields) {
    a <- pairs[[paste0(field, '_actual')]]
    b <- pairs[[paste0(field, '_expected')]]
    pairs[[paste0(field, '_equal')]] <- present &
      ((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b))
  }
  pairs$all_equal <- present & Reduce(`&`, pairs[paste0(fields, '_equal')])
  pairs
}

raw_indel_controls <- function(expected) {
  # Transport/comparator controls use copied rows, not new biological observations.
  expected <- head(expected, 4L)
  accepted <- function(actual) tryCatch(all(raw_indel_equal(actual, expected)$all_equal),
    error = function(e) FALSE)
  stopifnot(accepted(expected), accepted(expected[rev(seq_len(nrow(expected))), ]))
  controls <- c(missing = !accepted(expected[-1L, ]),
    duplicate = !accepted(rbind(expected, expected[1L, ])))
  extra <- expected[1L, ]
  extra$event_index <- max(expected$event_index) + 1L
  controls['extra'] <- !accepted(rbind(expected, extra))
  for (field in names(expected)) {
    changed <- expected
    value <- changed[[field]][1L]
    changed[[field]][1L] <- if (is.character(value)) paste0(value, '_corrupt') else
      if (is.na(value)) 1 else value + 1
    controls[field] <- !accepted(changed)
    changed <- expected
    changed[[field]] <- NULL
    controls[paste0('missing_', field)] <- !accepted(changed)
  }
  text <- expected
  text$cds[1L] <- 'NA'
  null <- text
  null$cds[1L] <- NA_character_
  controls['literal_NA_not_null'] <- !all(raw_indel_equal(null, text)$all_equal)
  ids <- c(1, 100000, 168000)
  spelling <- c('1', '100000', '168000')
  stopifnot(identical(raw_indel_ids(as.list(ids)), raw_indel_ids(as.list(spelling))),
    identical(raw_indel_ids(list(100000)), raw_indel_ids(list('100000'))),
    identical(raw_indel_ids(list(integer())), '[]'), is.na(raw_indel_ids(NULL)))
  invalid_ids <- list(fractional = 1.5, zero = 0, negative = -1, missing = NA_real_,
    infinite = Inf, inexact = 2^53, exponent_text = '1e+05',
    leading_zero = '0100000', duplicate = c(100000, 100000))
  for (name in names(invalid_ids)) controls[paste0('source_id_', name)] <-
    tryCatch({ raw_indel_ids(as.list(invalid_ids[[name]])); FALSE }, error = function(e) TRUE)
  stopifnot(all(controls))
  data.frame(control = names(controls), rejected = unname(controls))
}

raw_indel_oracle_controls <- function(case, events) {
  rejected <- function(changed) !isTRUE(tryCatch({
    raw_indel_case(changed, events)
    TRUE
  }, error = function(e) FALSE))
  stopifnot(!rejected(case))
  changed <- case
  changed$raw_haplotypes <- changed$raw_haplotypes[-1L]
  controls <- c(missing_source = rejected(changed))
  changed$raw_haplotypes <- c(case$raw_haplotypes, case$raw_haplotypes[1L])
  controls['extra_source'] <- rejected(changed)
  changed <- case
  changed$raw_haplotypes[[2L]] <- changed$raw_haplotypes[[1L]]
  controls['duplicate_source'] <- rejected(changed)
  for (field in c('id', 'source_reference', 'source_alternate', 'source_position1',
    'raw_gt', 'mapping_start', 'mapping_end', 'reference_cds')) {
    changed <- case
    value <- changed$raw_haplotypes[[1L]][[field]]
    changed$raw_haplotypes[[1L]][[field]] <- if (is.character(value))
      paste0(value, '_corrupt') else value + 1L
    controls[field] <- rejected(changed)
  }
  for (field in names(case$raw_haplotypes[[1L]])) {
    changed <- case
    changed$raw_haplotypes[[1L]][[field]] <- NULL
    controls[paste0('missing_', field)] <- rejected(changed)
  }
  changed <- case
  changed$raw_haplotypes[[1L]] <- c(changed$raw_haplotypes[[1L]], list(id = '1'))
  controls['duplicate_object_key'] <- rejected(changed)
  changed <- case
  changed$raw_haplotypes[[1L]]$mapping_start <- list(NULL)
  controls['invalid_mapping_type'] <- rejected(changed)
  changed <- case
  changed$raw_haplotypes[[1L]]$lanes[[2L]]$lane1 <- 1L
  controls['duplicate_lane'] <- rejected(changed)
  changed <- case
  changed$raw_haplotypes[[1L]]$lanes[[1L]]$sample <- 'different_sample'
  controls['changed_sample'] <- rejected(changed)
  stopifnot(all(controls))
  data.frame(control = paste0('oracle_', names(controls)), rejected = unname(controls))
}

raw_indel_empty_block_control <- function(native, partition, masks) {
  empty <- vapply(native$coding_blocks, function(x) !is.null(x) && nrow(x) == 0L, TRUE)
  at <- which(native$edit_count == 0L & empty)[1L]
  stopifnot(!is.na(at))
  original <- native[at, , drop = FALSE]
  expected <- raw_indel_native(original, partition, masks)
  changed <- original
  blocks <- original$coding_blocks[[1L]][rep(NA_integer_, 2L), , drop = FALSE]
  blocks$cds_start <- c(1, 1)
  blocks$reference <- blocks$alternate <- c('', '')
  blocks$length_change <- blocks$sequence_flags <- c(0L, 0L)
  blocks$event_indices <- list(integer(), integer())
  changed$coding_blocks[[1L]] <- blocks
  pairs <- raw_indel_equal(raw_indel_native(changed, partition, masks), expected)
  stopifnot(all(!pairs$block_count_equal), all(!pairs$all_equal))
  data.frame(control = 'native_extra_empty_blocks', rejected = TRUE)
}

raw_indel_compare <- function(con, events, partitions, out) {
  source <- do.call(rbind, lapply(seq_len(nrow(partitions)), function(i) {
    local <- codon_partition_events(events, partitions[i, , drop = FALSE])
    local$model_partition <- partitions$model_partition[i]
    local$event_offset <- partitions$event_offset[i]
    local
  }))
  expected <- raw_indel_expected(file.path(out, 'raw_oracle.stdout'), source)
  raw_indel_compare_observations(con, source, partitions, out, expected)
}

raw_indel_compare_observations <- function(con, source, partitions, out, expected) {
  stopifnot(nrow(expected) == 2L * nrow(source))
  first_line <- readLines(file.path(out, 'raw_oracle.stdout'), n = 1L)
  first_case <- jsonlite::fromJSON(first_line, simplifyVector = FALSE)
  controls <- rbind(raw_indel_controls(expected),
    raw_indel_oracle_controls(first_case, source[source$case_id == first_case$id, , drop = FALSE]))
  write.csv(controls, file.path(out, 'raw_controls.csv'), row.names = FALSE)
  DBI::dbWriteTable(con, 'raw_expected', expected)
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbExecute(con, paste('COPY raw_expected TO', q(file.path(out, 'raw_expected.parquet')),
    '(FORMAT PARQUET)'))
  DBI::dbRemoveTable(con, 'raw_expected')
  masks <- raw_indel_masks()
  summaries <- list()
  for (threads in c(1L, 4L)) {
    route <- paste0('source_records_', threads)
    actual <- do.call(rbind, lapply(seq_len(nrow(partitions)), function(i) {
      file <- file.path(out, paste0(route, '_partition_', i, '.parquet'))
      native <- DBI::dbGetQuery(con, paste('SELECT * FROM read_parquet(', q(file), ')'))
      if (threads == 1L && i == 1L) {
        controls <<- rbind(controls, raw_indel_empty_block_control(
          native, partitions[i, , drop = FALSE], masks))
        write.csv(controls, file.path(out, 'raw_controls.csv'), row.names = FALSE)
      }
      raw_indel_native(native, partitions[i, , drop = FALSE], masks)
    }))
    if (is.null(actual)) actual <- expected[FALSE, ]
    pairs <- raw_indel_equal(actual, expected)
    pairs$route <- route
    DBI::dbWriteTable(con, 'raw_pairs', pairs)
    DBI::dbExecute(con, paste('COPY raw_pairs TO',
      q(file.path(out, paste0('raw_pairs_', threads, '.parquet'))), '(FORMAT PARQUET)'))
    DBI::dbRemoveTable(con, 'raw_pairs')
    checks <- grep('_equal$', names(pairs), value = TRUE)
    summaries[[route]] <- data.frame(route = route, check = checks,
      expected_lanes = nrow(expected), actual_lanes = nrow(actual), pairs = nrow(pairs),
      failures = vapply(pairs[checks], function(x) sum(!x), 0L))
  }
  summary <- do.call(rbind, summaries)
  write.csv(summary, file.path(out, 'raw_summary.csv'), row.names = FALSE)
  print(summary)
  list(scope = 'constructed_GT_1|1_container_mechanics_not_raw_parser_or_HGVS_oracle',
    sample_mapping = list(oracle = 'translation_sample', native_sample_index = 0L),
    source_events = nrow(source), oracle_models = length(unique(source$case_id)),
    routes = c('source_records_1', 'source_records_4'),
    expected_lane_comparisons = 4L * nrow(source), summary = summary,
    all_equal = all(summary$failures == 0L))
}
