#!/usr/bin/env Rscript
# Network-free tests of the matrix and its comparator; these are not VEP observations.
source('test/duckvep/conformance/hgvs_cis_codon_differential.R')

inputs <- cis_codon_inputs()
stopifnot(nrow(inputs$models) == 6588L, nrow(inputs$mnv) == 6588L,
  nrow(inputs$snv) == 16470L)
small <- lapply(inputs, function(x) x[x$transcript_index %in% 0:1, ])
models <- small$models
records <- small$snv
actual <- data.frame(transcript_index = models$transcript_index, carrier_count = 1L,
  sequence_status = 'ok', projection_status = 'ok')
actual$contributors <- lapply(split(records, records$transcript_index), function(x) {
  x <- x[c('event_index', 'seq_region', 'position', 'reference', 'alternate')]
  x$evidence_flags <- 1L
  x$projection_status <- 'ok'
  x
})
actual$carriers <- rep(list(data.frame(sample_index = 0L, phase_set = NA_real_,
  haplotype_lane = 1L, ploidy = 1L)), 2L)
expected <- data.frame(id = models$id, hgvsp = c('p.Pro2Lys', 'p.Ala2Lys'))
controls <- cis_codon_controls(actual, records, models, expected)
stopifnot(length(controls) == 34L, all(controls),
  all(cis_codon_equal(expected[2:1, ], expected)))
unknown <- expected
unknown$hgvsp[1L] <- NA_character_
stopifnot(all(cis_codon_equal(unknown, unknown)),
  identical(cis_codon_equal(unknown, expected), c(FALSE, TRUE)))
reversed <- actual
reversed$contributors <- lapply(actual$contributors, function(x) x[nrow(x):1L, ])
stopifnot(cis_codon_check_native(reversed, records, models))
diploid <- actual
diploid$carrier_count <- 2L
diploid$carriers <- rep(list(data.frame(sample_index = 0L, phase_set = NA_real_,
  haplotype_lane = 1:2, ploidy = 2L)), 2L)
stopifnot(cis_codon_check_native(diploid, records, models, ploidy = 2L))

header <- c('##VEP="v116.0" API="v116"',
  '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Feature|HGVSp|Amino_acids">')
body <- with(small$mnv, paste(models$chrom, position, models$id, reference, alternate,
  '.', 'PASS', paste0('CSQ=', models$id, '|', models$id,
    ':', c('p.Pro2%3D', 'p.Ala2Ter'), '|', c('P', 'A/*')), sep = '\t'))
oracle <- cis_codon_oracle(c(header, body), small)
stopifnot(identical(oracle$hgvsp, c('p.Pro2=', 'p.Ala2Ter')),
  identical(oracle$protein, c('MPA*', 'M*')))
absent <- sub(':p.Pro2%3D', '', body, fixed = TRUE)
absent <- sub(paste0('|', models$id[1L], '|'), '||', absent, fixed = TRUE)
stopifnot(is.na(cis_codon_oracle(c(header, absent), small)$hgvsp[1L]))
rejected <- 0L
fails <- function(lines) {
  stopifnot(inherits(tryCatch(cis_codon_oracle(lines, small), error = identity), 'error'))
  rejected <<- rejected + 1L
}
fails(c(header, body[-1L]))
fails(c(header, body[2:1]))
fails(c(header, body, body[1L]))
fails(c(header, body[1L], body[1L]))
fails(c(sub('v116.0', 'v115.0', header, fixed = TRUE), body))
fails(c(header, sub('CSQ=CIS00000|', 'CSQ=CIS99999|', body, fixed = TRUE)))
fails(c(header, sub('CSQ=', 'CSQ=extra,', body, fixed = TRUE)))
fields <- strsplit(body[1L], '\t', fixed = TRUE)[[1L]]
for (field in 1:7) {
  changed <- fields
  changed[field] <- 'corrupt'
  fails(c(header, paste(changed, collapse = '\t'), body[2L]))
}
cat('Cis-codon matrix: 6,588 models, 16,470 SNVs;', sum(controls),
  'native/HGVS and', rejected, 'oracle corruption controls rejected\n')
