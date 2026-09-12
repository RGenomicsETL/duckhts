# Source-event identity is independent of replay/projection success.
duckvep_check_contributors <- function(observed, source) {
  fields <- c('event_index', 'seq_region', 'position', 'reference', 'alternate')
  if ('alt_index' %in% names(source)) fields <- c(fields, 'alt_index')
  stopifnot(all(fields %in% names(observed)), all(fields %in% names(source)),
    nrow(observed) == nrow(source), !anyNA(observed$event_index),
    !anyDuplicated(observed$event_index), !anyDuplicated(source$event_index),
    setequal(observed$event_index, source$event_index))
  source <- source[match(observed$event_index, source$event_index), , drop = FALSE]
  for (field in fields) stopifnot(all(observed[[field]] == source[[field]]))
  invisible(TRUE)
}
