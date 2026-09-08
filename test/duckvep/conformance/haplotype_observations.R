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
