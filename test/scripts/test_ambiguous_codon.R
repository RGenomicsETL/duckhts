#!/usr/bin/env Rscript
# Network-free comparator controls, independent of production annotations.
source('test/duckvep/conformance/ambiguous_codon_differential.R')
expected <- data.frame(event_index = 1:2, allele = c('A', 'C'),
  hgvsp = c('p.Ala2Thr', NA_character_), so = c('missense_variant', 'coding_sequence_variant'))
controls <- rbind(codon_controls(expected), codon_provenance_controls(FALSE),
  codon_provenance_controls(TRUE))
stopifnot(!anyDuplicated(controls$control), all(controls$rejected))
message('Ambiguous-codon comparator: ', nrow(controls), ' corruptions rejected')

# Reconstruct the checked-in failed baseline from raw VEP observations, not
# from stored equality flags. This validates retention, not biological agreement.
directory <- 'test/duckvep/conformance/data/ambiguous_codon_baseline'
manifest <- jsonlite::read_json(file.path(directory, 'receipt.json'), simplifyVector = TRUE)
hashes <- unlist(manifest$sha256)
stopifnot(setequal(names(hashes), c('pairs.parquet', 'cases.jsonl.gz', 'oracle.stdout.gz',
  'summary.csv', 'controls.csv', 'environment.stdout')))
for (name in names(hashes)) stopifnot(identical(unname(hashes[name]),
  duckvep_evidence_sha256(file.path(directory, name))))
read_records <- function(name) {
  connection <- gzfile(file.path(directory, name), 'rt')
  on.exit(close(connection))
  lapply(readLines(connection), jsonlite::fromJSON)
}
cases <- read_records('cases.jsonl.gz')
oracle <- read_records('oracle.stdout.gz')
stopifnot(length(cases) == 3000L,
  identical(vapply(cases, `[[`, '', 'id'), vapply(oracle, `[[`, '', 'id')))
events <- do.call(rbind, lapply(cases, function(x) data.frame(
  event_index = as.integer(x$variants$id), position = x$variants$position1,
  reference = x$variants$reference, allele = x$variants$alternate, cds = x$cds, table = x$table)))
expected <- do.call(rbind, lapply(oracle, function(x) {
  rows <- x$independent_hgvs
  data.frame(event_index = as.integer(rows$id), allele = rows$allele,
    hgvsp = sub('^.*:p\\.', 'p.', rows$hgvsp),
    so = vapply(rows$consequences, function(terms) paste(sort(terms), collapse = '&'), ''))
}))
stopifnot(nrow(events) == 28800L, nrow(expected) == 28800L,
  !anyDuplicated(events$event_index), !anyDuplicated(expected$event_index))
con <- DBI::dbConnect(duckdb::duckdb())
pairs <- DBI::dbGetQuery(con, paste('SELECT * FROM read_parquet(',
  DBI::dbQuoteString(con, file.path(directory, 'pairs.parquet')), ')'))
DBI::dbDisconnect(con, shutdown = TRUE)
stopifnot(nrow(pairs) == 230400L, all(pairs$actual_present), all(pairs$expected_present),
  all(pairs$seq_region == pairs$event_index - 1L),
  all(pairs$transcript_index == pairs$event_index - 1L))
input_at <- match(pairs$event_index, events$event_index)
for (field in c('position', 'reference', 'allele', 'cds', 'table'))
  stopifnot(all(pairs[[field]] == events[[field]][input_at]))
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
  all(summary == retained_summary), sum(!pairs$hgvsp_equal) == 12992L,
  sum(!pairs$so_equal, na.rm = TRUE) == 3248L)
retained_controls <- read.csv(file.path(directory, 'controls.csv'))
stopifnot(identical(controls$control, retained_controls$control),
  identical(controls$rejected, retained_controls$rejected))
message('Complete retained baseline reconstructed: 230,400 HGVSp / 57,600 SO comparisons; failures preserved')
