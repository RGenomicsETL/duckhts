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
