# Complete sequence groups, applied-source identity sets and optional sample counts.
# Source sets do not prove physical-edit multiplicity or file-lane association.
canonical <- function(rows, provenance = TRUE, samples = FALSE) {
  if (!length(rows)) return(list())
  keys <- vapply(rows, function(x) jsonlite::toJSON(list(cds=x$cds, protein=x$protein),
    auto_unbox=TRUE, na='null'), '')
  lapply(split(rows, keys), function(group) {
    value <- list(cds=group[[1L]]$cds, protein=group[[1L]]$protein,
      count=sum(vapply(group, function(x) as.numeric(x$count), 0)))
    if (provenance) value$contributors <- sort(unique(unlist(lapply(group, `[[`, 'contributors'),
      use.names=FALSE)))
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
  setNames(rejected, paste0('output_', names(rejected)))
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
