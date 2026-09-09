# Complete sequence groups, applied-source identity sets and optional sample counts.
# Source sets do not prove physical-edit multiplicity or file-lane association.
canonical <- function(rows, provenance = TRUE, samples = FALSE) {
  if (!length(rows)) return(list())
  if (provenance) stopifnot(all(vapply(rows, function(row) {
    sources <- row[['contributors']]
    'contributors' %in% names(row) && !is.null(sources) &&
      (is.character(sources) || is.list(sources)) &&
      all(vapply(sources, function(x) is.character(x) && length(x) == 1L && !is.na(x), TRUE))
  }, TRUE)))
  keys <- vapply(rows, function(x) jsonlite::toJSON(list(cds=x$cds, protein=x$protein),
    auto_unbox=TRUE, na='null'), '')
  lapply(split(rows, keys), function(group) {
    value <- list(cds=group[[1L]]$cds, protein=group[[1L]]$protein,
      count=sum(vapply(group, function(x) as.numeric(x$count), 0)))
    if (provenance) value$contributors <- sort(unique(as.character(unlist(
      lapply(group, `[[`, 'contributors'), use.names=FALSE))))
    if (samples) {
      sample_names <- sort(unique(unlist(lapply(group, function(x) names(x$samples)), use.names=FALSE)))
      value$samples <- setNames(lapply(sample_names, function(name) sum(vapply(group, function(x) {
        count <- x$samples[[name]]
        if (is.null(count)) 0 else as.numeric(count)
      }, 0))), sample_names)
    }
    value
  })
}

haplotype_group_controls <- function() {
  oracle <- jsonlite::fromJSON('[
    {"cds":"ATG","protein":"M","count":1,"contributors":[],"samples":{"s0":1}},
    {"cds":"ATC","protein":"I","count":1,"contributors":["a","b"],"samples":{"s1":1}}
  ]', simplifyVector = FALSE)
  native <- oracle
  native[[1L]]$contributors <- character()
  native[[2L]]$contributors <- c('a', 'b')
  expected <- canonical(oracle, samples = TRUE)
  stopifnot(identical(expected, canonical(native, samples = TRUE)),
    length(expected) == 2L, all(vapply(expected, function(x) 'contributors' %in% names(x), TRUE)))
  merged <- native
  merged[[1L]]$cds <- merged[[2L]]$cds
  merged[[1L]]$protein <- merged[[2L]]$protein
  group <- canonical(merged, samples = TRUE)
  stopifnot(length(group) == 1L, identical(group[[1L]]$contributors, c('a', 'b')),
    group[[1L]]$count == 2, identical(group[[1L]]$samples, list(s0 = 1, s1 = 1)))
  corrupt <- list(missing_field = native, null_field = native, numeric_source = native,
    missing_source = native, extra_source = native, na_source = native,
    dropped_row = native[-1L], duplicate_row = c(native, native[1L]))
  corrupt$missing_field[[1L]]$contributors <- NULL
  corrupt$null_field[[1L]]['contributors'] <- list(NULL)
  corrupt$numeric_source[[1L]]$contributors <- 1L
  corrupt$missing_source[[2L]]$contributors <- 'a'
  corrupt$extra_source[[1L]]$contributors <- 'extra'
  corrupt$na_source[[1L]]$contributors <- NA_character_
  rejected <- vapply(corrupt, function(x) {
    observed <- try(canonical(x, samples = TRUE), silent = TRUE)
    inherits(observed, 'try-error') || !identical(expected, observed)
  }, TRUE)
  stopifnot(all(rejected))
  setNames(rejected, paste0('group_', names(rejected)))
}

# Every native row belongs to one declared transcript and at least one carrier.
native_haplotype_rows <- function(actual, transcript_indices) {
  stopifnot(is.data.frame(actual), is.numeric(transcript_indices),
    !anyNA(transcript_indices), !anyDuplicated(transcript_indices),
    all(c('transcript_index', 'carrier_count', 'carriers') %in% names(actual)),
    is.numeric(actual$transcript_index), !anyNA(actual$transcript_index),
    setequal(actual$transcript_index, transcript_indices),
    is.numeric(actual$carrier_count), all(actual$carrier_count > 0),
    all(actual$carrier_count == vapply(actual$carriers, nrow, 1L)))
  rows <- split(seq_len(nrow(actual)), actual$transcript_index)
  rows[as.character(transcript_indices)]
}

replay_comparisons_passed <- function(summary) {
  stopifnot(all(c('equal', 'counts_equal', 'replay_lanes_equal') %in% names(summary)))
  summary$equal & summary$counts_equal & summary$replay_lanes_equal
}

haplotype_output_controls <- function(witness) {
  stopifnot(nrow(witness) > 0L)
  indices <- unique(witness$transcript_index)
  native_haplotype_rows(witness, indices)
  rejects <- function(rows, domain = indices)
    inherits(try(native_haplotype_rows(rows, domain), silent = TRUE), 'try-error')
  unknown <- null_tx <- zero <- witness[1L, , drop = FALSE]
  unknown$transcript_index <- max(indices) + 1
  null_tx$transcript_index <- NA
  zero$carriers[[1L]] <- zero$carriers[[1L]][FALSE, , drop = FALSE]
  zero$carrier_count <- 0
  extra <- list(unknown_transcript = unknown, null_transcript = null_tx, zero_carriers = zero)
  rejected <- vapply(extra, function(row) rejects(rbind(witness, row)), TRUE)
  count <- witness
  count$carrier_count[1L] <- count$carrier_count[1L] + 1
  rejected <- c(rejected,
    missing_transcript = rejects(witness[witness$transcript_index != indices[1L], , drop = FALSE]),
    carrier_count = rejects(count),
    missing_count_field = rejects(witness[setdiff(names(witness), 'carrier_count')]),
    duplicate_transcript_domain = rejects(witness, c(indices, indices[1L])))
  comparison <- data.frame(equal = TRUE, counts_equal = TRUE, replay_lanes_equal = TRUE)
  stopifnot(isTRUE(replay_comparisons_passed(comparison)))
  comparison$counts_equal <- FALSE
  rejected <- c(rejected, total_count = !replay_comparisons_passed(comparison))
  stopifnot(all(rejected))
  c(setNames(rejected, paste0('output_', names(rejected))), haplotype_group_controls())
}

# One upstream mutator result per sample/file lane. Applied sources are identity
# sets, not physical-edit counts: a single source can contain several edit islands.
canonical_replay_lanes <- function(rows) {
  string <- function(x) is.character(x) && length(x) == 1L && !is.na(x)
  nullable <- function(x) is.null(x) || string(x)
  fields <- c('sample', 'lane1', 'cds', 'protein', 'applied_sources')
  source_fields <- c('allele_key', 'source_id', 'source_key')
  if (!is.list(rows)) return(NULL)
  keys <- character(length(rows))
  result <- vector('list', length(rows))
  for (i in seq_along(rows)) {
    x <- rows[[i]]
    if (!is.list(x) || anyDuplicated(names(x)) || !all(fields %in% names(x)) ||
        !string(x$sample) || !is.numeric(x$lane1) || length(x$lane1) != 1L ||
        !is.finite(x$lane1) || x$lane1 < 1 || x$lane1 > 65535 ||
        x$lane1 != floor(x$lane1) || !nullable(x$cds) || !nullable(x$protein) ||
        !is.list(x$applied_sources)) return(NULL)
    sources <- x$applied_sources
    if (any(!vapply(sources, function(s) is.list(s) && !anyDuplicated(names(s)) &&
        all(source_fields %in% names(s)) && all(vapply(s[source_fields], string, TRUE)), TRUE)))
      return(NULL)
    sources <- lapply(sources, function(s) s[source_fields])
    source_keys <- vapply(sources, jsonlite::toJSON, '', auto_unbox = TRUE)
    if (anyDuplicated(source_keys)) return(NULL)
    keys[i] <- jsonlite::toJSON(list(sample = x$sample, lane1 = as.integer(x$lane1)), auto_unbox = TRUE)
    result[[i]] <- list(sample = x$sample, lane1 = as.integer(x$lane1), cds = x$cds,
      protein = x$protein, applied_sources = sources[order(source_keys)])
  }
  if (anyDuplicated(keys)) return(NULL)
  result[order(keys)]
}

replay_lanes_equal <- function(expected, observed) {
  e <- canonical_replay_lanes(expected)
  o <- canonical_replay_lanes(observed)
  !is.null(e) && !is.null(o) && identical(e, o)
}

replay_lane_controls <- function(witness, require_shared = TRUE) {
  stopifnot(replay_lanes_equal(witness, witness), length(witness) >= 2L,
    length(witness[[1L]]$applied_sources) > 0L)
  if (require_shared) stopifnot(length(witness) == 6L)
  corrupt <- list(missing = witness[-1L], duplicate = c(witness, witness[1L]),
    cds = witness, protein = witness, allele_key = witness, source_key = witness,
    source_id = witness, sample = witness, absent_field = witness, invalid_source = witness,
    swapped_lanes = witness)
  corrupt$cds[[1L]]$cds <- paste0(witness[[1L]]$cds, 'A')
  corrupt$protein[[1L]]$protein <- paste0(witness[[1L]]$protein, 'X')
  for (field in c('allele_key', 'source_key', 'source_id'))
    corrupt[[field]][[1L]]$applied_sources[[1L]][[field]] <- 'wrong'
  corrupt$sample[[1L]]$sample <- 'unknown'
  corrupt$absent_field[[1L]]$cds <- NULL
  corrupt$invalid_source[[1L]]$applied_sources[[1L]] <- 'invalid'
  pair <- which(vapply(witness, function(x) x$sample == witness[[1L]]$sample, TRUE))
  stopifnot(length(pair) == 2L, !identical(witness[[pair[1L]]]$cds, witness[[pair[2L]]]$cds))
  corrupt$swapped_lanes[[pair[1L]]]$lane1 <- witness[[pair[2L]]]$lane1
  corrupt$swapped_lanes[[pair[2L]]]$lane1 <- witness[[pair[1L]]]$lane1
  if (require_shared) {
    shared <- which(vapply(witness, function(x) identical(x$cds, witness[[pair[2L]]]$cds), TRUE))
    stopifnot(length(shared) > 1L, length(witness[[shared[1L]]]$applied_sources) > 0L)
    corrupt$missing_lane_source <- witness
    corrupt$missing_lane_source[[shared[1L]]]$applied_sources <-
      witness[[shared[1L]]]$applied_sources[-1L]
  }
  groups <- function(rows) canonical(lapply(rows, function(x) list(cds = x$cds,
    protein = x$protein, count = 1L,
    contributors = vapply(x$applied_sources, `[[`, '', 'source_id'),
    samples = setNames(list(1L), x$sample))), samples = TRUE)
  for (name in intersect(c('allele_key', 'source_key', 'swapped_lanes', 'missing_lane_source'), names(corrupt)))
    stopifnot(identical(groups(witness), groups(corrupt[[name]])))
  rejected <- vapply(corrupt, function(x) !replay_lanes_equal(witness, x), TRUE)
  setNames(rejected, paste0('lane_', names(rejected)))
}

# Join physical edit IDs to the retained source relation. The allele key uses
# the full source ALT in transcript orientation, including its VCF anchor.
native_replay_lanes <- function(actual, records, strand, samples) {
  result <- list()
  for (i in seq_len(nrow(actual))) {
    ids <- unique(as.numeric(unlist(actual$coding_blocks[[i]]$event_indices, use.names = FALSE)))
    contributors <- actual$contributors[[i]]
    j <- match(ids, contributors$event_index)
    k <- match(ids, records$event_index)
    if (anyNA(j) || anyNA(k)) return(NULL)
    sources <- lapply(seq_along(ids), function(n) {
      r <- records[k[n], ]
      source_key <- paste(r$position, r$position + nchar(r$reference) - 1L,
        paste(c(r$reference, strsplit(r$alt, ',', fixed = TRUE)[[1L]]), collapse = '/'), sep = '_')
      alt <- contributors$alternate[j[n]]
      if (strand < 0L) alt <- paste(rev(strsplit(chartr('ACGT', 'TGCA', alt), '',
        fixed = TRUE)[[1L]]), collapse = '')
      list(allele_key = paste0(if (nchar(alt)) alt else '-', '|', source_key),
        source_id = r$source_id, source_key = source_key)
    })
    calls <- actual$carriers[[i]]
    if (anyNA(calls$sample_index) || any(calls$sample_index < 0L | calls$sample_index >= length(samples)))
      return(NULL)
    for (n in seq_len(nrow(calls))) result[[length(result) + 1L]] <- list(
      sample = samples[calls$sample_index[n] + 1L], lane1 = as.integer(calls$haplotype_lane[n]),
      cds = if (is.na(actual$cds[i])) NULL else actual$cds[i],
      protein = if (is.na(actual$protein[i])) NULL else actual$protein[i], applied_sources = sources)
  }
  result
}

# Raw mutation flags are distinct from sequence-group flags and from the final
# sequence length difference. Overlapping replacements can give different sums.
haplotype_raw_flags <- function(flags) {
  fields <- c('indel', 'frameshift', 'length_diff')
  if (!is.list(flags) || (length(flags) && is.null(names(flags))) || anyDuplicated(names(flags)) ||
      any(!names(flags) %in% fields)) return(NULL)
  values <- setNames(numeric(3L), fields)
  for (name in names(flags)) {
    x <- flags[[name]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x != floor(x) || abs(x) >= 2^53 ||
        (name != 'length_diff' && !x %in% 0:1)) return(NULL)
    values[name] <- x
  }
  values
}

haplotype_raw_flag_bits <- function(flags) {
  values <- haplotype_raw_flags(flags)
  if (is.null(values)) return(NA_integer_)
  as.integer(values['indel'] + if (values['frameshift'])
    if (values['length_diff'] %% 3) 2L else 4L else 0L)
}

haplotype_flag_categories <- function(bits) {
  c('frameshift', 'indel', 'resolved_frameshift')[bitwAnd(bits, c(2L, 1L, 4L)) != 0L]
}

# Validate the complete observed mutator-lane domain. Reference-only samples added
# by _add_reference_haplotypes are not synthesized as mutation observations.
haplotype_lane_metadata <- function(lanes) {
  scalar <- function(x) is.character(x) && length(x) == 1L && !is.na(x)
  integer_value <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == floor(x)
  fields <- c('sample', 'lane1', 'cds', 'protein', 'traversal_ordinal', 'flags')
  if (!is.list(lanes)) return(NULL)
  keys <- character(length(lanes)); ordinals <- numeric(length(lanes))
  for (i in seq_along(lanes)) {
    x <- lanes[[i]]
    if (!is.list(x) || anyDuplicated(names(x)) || !all(fields %in% names(x)) ||
        !scalar(x$sample) || !integer_value(x$lane1) || x$lane1 < 1L || x$lane1 > 65535L ||
        (!is.null(x$cds) && !scalar(x$cds)) || (!is.null(x$protein) && !scalar(x$protein)) ||
        !integer_value(x$traversal_ordinal) || is.null(haplotype_raw_flags(x$flags))) return(NULL)
    keys[i] <- jsonlite::toJSON(list(sample = x$sample, lane1 = x$lane1), auto_unbox = TRUE)
    ordinals[i] <- x$traversal_ordinal
  }
  if (anyDuplicated(keys) || !identical(sort(ordinals), as.numeric(seq_along(lanes)))) return(NULL)
  lanes[order(ordinals)]
}

native_lane_flags_equal <- function(lanes, actual, samples) {
  lanes <- haplotype_lane_metadata(lanes)
  if (is.null(lanes) || !is.data.frame(actual) ||
      !all(c('sequence_flags', 'carriers') %in% names(actual)) ||
      !is.character(samples) || anyNA(samples) || anyDuplicated(samples)) return(FALSE)
  expected <- vapply(lanes, function(x) haplotype_raw_flag_bits(x$flags), 1L)
  keys <- vapply(lanes, function(x) jsonlite::toJSON(list(sample = x$sample,
    lane1 = x$lane1), auto_unbox = TRUE), '')
  observed_keys <- character(); observed <- integer()
  for (i in seq_len(nrow(actual))) {
    flags <- actual$sequence_flags[i]
    carriers <- actual$carriers[[i]]
    if (!is.numeric(flags) || is.na(flags) || flags != floor(flags) || flags < 0 || flags > 15 ||
        !is.data.frame(carriers) || !all(c('sample_index', 'haplotype_lane') %in% names(carriers)) ||
        !all(vapply(carriers[c('sample_index', 'haplotype_lane')], is.numeric, TRUE)) ||
        anyNA(carriers[c('sample_index', 'haplotype_lane')]) ||
        any(carriers$sample_index < 0 | carriers$sample_index >= length(samples) |
          carriers$sample_index != floor(carriers$sample_index) |
          carriers$haplotype_lane < 1 | carriers$haplotype_lane > 65535 |
          carriers$haplotype_lane != floor(carriers$haplotype_lane))) return(FALSE)
    for (j in seq_len(nrow(carriers))) {
      observed_keys <- c(observed_keys, jsonlite::toJSON(list(
        sample = samples[carriers$sample_index[j] + 1L], lane1 = carriers$haplotype_lane[j]), auto_unbox = TRUE))
      observed <- c(observed, bitwAnd(as.integer(flags), 7L))
    }
  }
  !anyDuplicated(observed_keys) && setequal(keys, observed_keys) &&
    identical(expected, observed[match(keys, observed_keys)])
}

# Check the actual upstream group owner, not a union or a selected flag value.
# Output is the untouched custom oracle record or original Runner JSON. A CDS
# group with no mutator lane must be the zero-flag reference-only group.
haplotype_group_metadata_equal <- function(observation, output, container_json = FALSE) {
  lanes <- haplotype_lane_metadata(observation$replay_lanes)
  groups <- observation$cds_group_metadata
  published <- if (container_json) output$cds_haplotypes else output$haplotypes
  if (is.null(lanes) || !is.list(groups) || !is.list(published) ||
      !is.character(observation$reference_cds) || length(observation$reference_cds) != 1L ||
      is.na(observation$reference_cds)) return(FALSE)
  key <- function(rows, field) {
    if (any(!vapply(rows, function(x) is.list(x) && !anyDuplicated(names(x)) &&
        is.character(x[[field]]) && length(x[[field]]) == 1L && !is.na(x[[field]]), TRUE))) return(NULL)
    vapply(rows, `[[`, '', field)
  }
  keys <- key(groups, 'cds'); output_keys <- key(published, if (container_json) 'seq' else 'cds')
  if (is.null(keys) || is.null(output_keys) || anyDuplicated(keys) || anyDuplicated(output_keys) ||
      !setequal(keys, output_keys)) return(FALSE)
  eligible <- Filter(function(x) !is.null(x$protein) && nzchar(x$protein), lanes)
  if (any(!vapply(eligible, function(x) !is.null(x$cds) && x$cds %in% keys, TRUE))) return(FALSE)
  categories <- function(x) {
    if (!is.list(x) && !is.character(x)) return(NULL)
    if (any(!vapply(x, function(y) is.character(y) && length(y) == 1L && !is.na(y), TRUE))) return(NULL)
    x <- as.character(unlist(x, use.names = FALSE))
    if (anyDuplicated(x) || any(!x %in% c('indel', 'frameshift', 'resolved_frameshift'))) return(NULL)
    sort(x)
  }
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    candidates <- Filter(function(x) identical(x$cds, group$cds), eligible)
    if (!length(candidates) && !identical(group$cds, observation$reference_cds)) return(FALSE)
    first <- if (length(candidates)) candidates[[1L]]$flags else list()
    raw <- haplotype_raw_flags(group$flags)
    flags <- categories(group$categories)
    expected <- haplotype_flag_categories(haplotype_raw_flag_bits(first))
    if (!identical(sort(names(group$flags)), sort(c('indel', 'frameshift', 'length_diff'))) ||
        is.null(raw) || is.null(flags) || !identical(raw, haplotype_raw_flags(first)) ||
        !identical(flags, expected)) return(FALSE)
    row <- published[[match(keys[i], output_keys)]]
    if (container_json) {
      if (!is.numeric(row$has_indel) || length(row$has_indel) != 1L || is.na(row$has_indel) ||
          row$has_indel != raw['indel']) return(FALSE)
    } else {
      flags <- categories(row$flags)
      if (is.null(flags) || !identical(flags, expected)) return(FALSE)
    }
  }
  TRUE
}

haplotype_metadata_controls <- function(observation, output, container_json = FALSE) {
  check <- function(x, y = output) haplotype_group_metadata_equal(x, y, container_json)
  stopifnot(check(observation), length(observation$replay_lanes) > 0L)
  corrupt <- list(indel = observation, frameshift = observation, length_diff = observation,
    missing_flags = observation, invalid_flags = observation, duplicate_rank = observation,
    missing_rank = observation, invalid_rank = observation, group_flags = observation,
    group_categories = observation, missing_group_flag = observation,
    missing_group = observation, duplicate_group = observation)
  first <- which.min(vapply(observation$replay_lanes, `[[`, 1L, 'traversal_ordinal'))
  original <- haplotype_raw_flags(observation$replay_lanes[[first]]$flags)
  for (name in c('indel', 'frameshift'))
    corrupt[[name]]$replay_lanes[[first]]$flags[[name]] <- 1 - original[name]
  corrupt$length_diff$replay_lanes[[first]]$flags$length_diff <- original['length_diff'] + 1
  corrupt$missing_flags$replay_lanes[[first]]$flags <- NULL
  corrupt$invalid_flags$replay_lanes[[first]]$flags$indel <- 2L
  stopifnot(length(observation$replay_lanes) > 1L)
  corrupt$duplicate_rank$replay_lanes[[2L]]$traversal_ordinal <-
    corrupt$duplicate_rank$replay_lanes[[1L]]$traversal_ordinal
  corrupt$missing_rank$replay_lanes[[first]]$traversal_ordinal <- NULL
  corrupt$invalid_rank$replay_lanes[[first]]$traversal_ordinal <- 0L
  corrupt$group_flags$cds_group_metadata[[1L]]$flags$length_diff <-
    corrupt$group_flags$cds_group_metadata[[1L]]$flags$length_diff + 1
  corrupt$group_categories$cds_group_metadata[[1L]]$categories <- list('invalid')
  corrupt$missing_group_flag$cds_group_metadata[[1L]]$flags$length_diff <- NULL
  corrupt$missing_group$cds_group_metadata <- corrupt$missing_group$cds_group_metadata[-1L]
  corrupt$duplicate_group$cds_group_metadata <- c(corrupt$duplicate_group$cds_group_metadata,
    corrupt$duplicate_group$cds_group_metadata[1L])
  rejected <- vapply(corrupt, function(x) !check(x), TRUE)
  bad <- output
  if (container_json) bad$cds_haplotypes[[1L]]$has_indel <- 2L else bad$haplotypes[[1L]]$flags <- list('invalid')
  rejected <- c(rejected, published_flags = !check(observation, bad))
  setNames(rejected, paste0('metadata_', names(rejected)))
}
