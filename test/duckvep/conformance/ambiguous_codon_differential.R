#!/usr/bin/env Rscript
# Exact original SNVs: uploaded-allele ambiguity and surrounding codon ambiguity.
source('scripts/duckvep_evidence.R')
source('test/duckvep/conformance/contributor_identity.R')

codon_check_provenance <- function(actual, events, raw) {
  ploidy <- if (raw) 2L else 1L
  stopifnot(all(c('transcript_index', 'carrier_count', 'contributors', 'carriers') %in% names(actual)),
    nrow(actual) == nrow(events), !anyDuplicated(actual$transcript_index),
    setequal(actual$transcript_index, events$transcript_index), all(actual$carrier_count == ploidy))
  events <- events[match(actual$transcript_index, events$transcript_index), , drop = FALSE]
  if (raw) events$alt_index <- 1L
  for (i in seq_len(nrow(actual))) {
    duckvep_check_contributors(actual$contributors[[i]], events[i, , drop = FALSE])
    carriers <- actual$carriers[[i]]
    stopifnot(all(c('sample_index', 'phase_set', 'haplotype_lane', 'ploidy') %in% names(carriers)),
      nrow(carriers) == ploidy, all(carriers$sample_index == 0L), all(is.na(carriers$phase_set)),
      identical(sort(as.integer(carriers$haplotype_lane)), seq_len(ploidy)), all(carriers$ploidy == ploidy))
  }
  invisible(TRUE)
}

codon_provenance_controls <- function(raw) {
  ploidy <- if (raw) 2L else 1L
  events <- data.frame(event_index = 1L, transcript_index = 0L, seq_region = 0L,
    position = 4L, reference = 'G', alternate = 'A')
  actual <- data.frame(transcript_index = 0L, carrier_count = ploidy)
  actual$contributors <- list(events[setdiff(names(events), 'transcript_index')])
  if (raw) actual$contributors[[1L]]$alt_index <- 1L
  actual$carriers <- list(data.frame(sample_index = 0L, phase_set = NA_integer_,
    haplotype_lane = seq_len(ploidy), ploidy = ploidy))
  rejected <- function(x) !isTRUE(tryCatch(codon_check_provenance(x, events, raw),
    error = function(e) FALSE))
  stopifnot(!rejected(actual))
  controls <- c(missing_leaf = rejected(actual[FALSE, ]),
    duplicate_leaf = rejected(rbind(actual, actual)))
  for (field in c('transcript_index', 'carrier_count')) {
    changed <- actual
    changed[[field]][1L] <- changed[[field]][1L] + 1L
    controls[field] <- rejected(changed)
  }
  for (part in c('contributors', 'carriers')) {
    changed <- actual
    changed[[part]][[1L]] <- changed[[part]][[1L]][FALSE, ]
    controls[paste0('missing_', part)] <- rejected(changed)
    changed[[part]][[1L]] <- rbind(actual[[part]][[1L]], actual[[part]][[1L]])
    controls[paste0('duplicate_', part)] <- rejected(changed)
    for (field in names(actual[[part]][[1L]])) {
      changed <- actual
      value <- changed[[part]][[1L]][[field]][1L]
      changed[[part]][[1L]][[field]][1L] <- if (is.na(value)) 1L else
        if (is.character(value)) paste0(value, '_corrupt') else value + 1L
      controls[paste(part, field, sep = '_')] <- rejected(changed)
      changed[[part]][[1L]][[field]] <- NULL
      controls[paste('missing', part, field, sep = '_')] <- rejected(changed)
    }
  }
  stopifnot(all(controls))
  data.frame(control = paste(if (raw) 'raw' else 'decoded', names(controls), sep = '_'),
    rejected = unname(controls))
}

codon_equal <- function(actual, expected) {
  required <- c('event_index', 'allele', 'hgvsp', 'so')
  stopifnot(identical(names(actual), required), identical(names(expected), required))
  key <- function(x) paste(x$event_index, x$allele, sep = '/')
  if (anyDuplicated(key(actual)) || anyDuplicated(key(expected)))
    stop('duplicate event/ALT comparison key', call. = FALSE)
  actual$actual_present <- TRUE
  expected$expected_present <- TRUE
  pairs <- merge(actual, expected, by = c('event_index', 'allele'), all = TRUE,
    suffixes = c('_actual', '_expected'), sort = TRUE)
  present <- !is.na(pairs$actual_present) & !is.na(pairs$expected_present)
  equal <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
  pairs$hgvsp_equal <- present & equal(pairs$hgvsp_actual, pairs$hgvsp_expected)
  pairs$so_equal <- present & equal(pairs$so_actual, pairs$so_expected)
  pairs
}

codon_controls <- function(expected) {
  accepted <- function(x) {
    tryCatch({
      pairs <- codon_equal(x, expected)
      all(pairs$hgvsp_equal & pairs$so_equal)
    }, error = function(e) FALSE)
  }
  stopifnot(accepted(expected), accepted(expected[rev(seq_len(nrow(expected))), ]))
  controls <- c(missing = !accepted(expected[-1L, ]),
    duplicate = !accepted(rbind(expected, expected[1L, ])))
  extra <- expected[1L, ]
  extra$event_index <- max(expected$event_index) + 1L
  controls['extra'] <- !accepted(rbind(expected, extra))
  for (field in names(expected)) {
    changed <- expected
    at <- which(!is.na(changed[[field]]))[1L]
    changed[[field]][at] <- if (is.numeric(changed[[field]])) -1 else
      paste0(changed[[field]][at], '_corrupt')
    controls[field] <- !accepted(changed)
  }
  for (field in c('hgvsp', 'so')) {
    changed <- expected
    changed[[field]][which(!is.na(changed[[field]]))[1L]] <- NA_character_
    controls[paste0(field, '_missing')] <- !accepted(changed)
  }
  changed <- expected
  changed$hgvsp[which(is.na(changed$hgvsp))[1L]] <- 'p.Ala2Thr'
  controls['invented_hgvsp'] <- !accepted(changed)
  stopifnot(all(controls))
  data.frame(control = names(controls), rejected = unname(controls))
}

main <- function() {
  suppressPackageStartupMessages(library(DBI))
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option('--extension', default = 'build/release/duckhts.duckdb_extension'),
    optparse::make_option('--evidence-out', dest = 'evidence_out', default = '',
      help = 'Retain a compact complete comparison bundle in a new directory'),
    optparse::make_option('--vep-prefix', dest = 'vep_prefix',
      default = Sys.getenv('VEP_PREFIX', '/root/miniconda3/envs/vep'))
  )))
  revision <- duckvep_evidence_revision('.')
  extension <- normalizePath(opt$extension, mustWork = TRUE)
  sources <- c('test/duckvep/conformance/ambiguous_codon_differential.R',
    'test/duckvep/conformance/contributor_identity.R',
    'test/duckvep/conformance/reference_translation_oracle.pl', 'scripts/duckvep_evidence.R',
    list.files('src/duckvep', recursive = TRUE, full.names = TRUE, pattern = '\\.[ch]$'))
  source_hashes <- vapply(sources, duckvep_evidence_sha256, '')
  extension_hash <- duckvep_evidence_sha256(extension)
  pins <- c(vep = '57ea5c52340acc1f156267f810ad162e26597082',
    variation = '2fb834b987ede3824e200197a838ce11e91aeb4b')
  mirrors <- setNames(normalizePath(paste0('.sync/ensembl-', names(pins))), names(pins))
  for (name in names(pins)) {
    stopifnot(identical(duckvep_evidence_command('git',
      c('-C', mirrors[[name]], 'rev-parse', 'HEAD'), 'oracle revision'), unname(pins[name])),
      !length(duckvep_evidence_command('git', c('-C', mirrors[[name]], 'status', '--porcelain'),
        'oracle checkout')))
  }
  prefix <- normalizePath(opt$vep_prefix, mustWork = TRUE)
  out <- tempfile('ambiguous_codon_', 'test/duckvep/conformance/results')
  stopifnot(dir.create(out))
  out <- normalizePath(out)
  message('Ambiguous-codon artifacts: ', out)
  command <- function(args, label) {
    status <- system2('micromamba', shQuote(args),
      stdout = file.path(out, paste0(label, '.stdout')),
      stderr = file.path(out, paste0(label, '.stderr')))
    stopifnot(status == 0L)
  }
  command(c('list', '--explicit', '-p', prefix), 'environment')
  stopifnot(identical(duckvep_evidence_explicit_packages(readLines(file.path(out, 'environment.stdout'))),
    duckvep_evidence_explicit_packages(readLines(
      'test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt'))))
  modules <- c(file.path(prefix, 'share/ensembl-vep-116.0-0/Bio/EnsEMBL',
    c('Transcript.pm', 'Translation.pm', 'SeqEdit.pm')),
    file.path(prefix, 'lib/perl5/site_perl/Bio', c('Tools/CodonTable.pm', 'PrimarySeqI.pm')),
    file.path(mirrors[['variation']], 'modules/Bio/EnsEMBL/Variation',
      c('TranscriptHaplotypeContainer.pm', 'TranscriptVariation.pm', 'TranscriptVariationAllele.pm',
        'VariationFeatureOverlapAllele.pm', 'Utils/VariationEffect.pm')))
  module_hashes <- vapply(modules, duckvep_evidence_sha256, '')
  bases <- strsplit('ACGTN', '', fixed = TRUE)[[1L]]
  triplets <- do.call(paste0, expand.grid(rep(list(bases), 3L), stringsAsFactors = FALSE))
  grid <- expand.grid(table = c(1:6, 9:14, 16, 21:31), codon = triplets,
    stringsAsFactors = FALSE)
  events <- list()
  cases <- vector('list', nrow(grid))
  for (i in seq_len(nrow(grid))) {
    variants <- list()
    for (offset in 1:3) for (alt in bases[1:4]) {
      ref <- substr(grid$codon[i], offset, offset)
      if (ref == alt) next
      index <- length(events) + 1L
      variants[[length(variants) + 1L]] <- list(id = as.character(index),
        position1 = offset + 3L, reference = ref, alternate = alt)
      events[[index]] <- data.frame(event_index = index, seq_region = index - 1L,
        transcript_index = index - 1L, position = offset + 3L, reference = ref,
        alternate = alt, cds = paste0('ATG', grid$codon[i], 'TAA'), table = grid$table[i],
        codon = grid$codon[i], case_id = paste(grid$table[i], grid$codon[i], sep = '/'))
    }
    cases[[i]] <- list(id = paste(grid$table[i], grid$codon[i], sep = '/'), cds = paste0('ATG', grid$codon[i], 'TAA'),
      table = grid$table[i], edits = list(), variants = variants)
  }
  events <- do.call(rbind, events)
  stopifnot(nrow(grid) == 3000L, nrow(events) == 28800L, !anyDuplicated(events$event_index))
  saveRDS(events, file.path(out, 'events.rds'))
  input <- file.path(out, 'cases.jsonl')
  writeLines(vapply(cases, jsonlite::toJSON, '', auto_unbox = TRUE), input)
  libs <- paste(c(file.path(mirrors, 'modules'),
    file.path(prefix, 'share/ensembl-vep-116.0-0')), collapse = ':')
  command(c('run', '--clean-env', '--env', paste0('PERL5LIB=', libs), '-p', prefix, 'perl',
    normalizePath('test/duckvep/conformance/reference_translation_oracle.pl'), input), 'oracle')
  oracle <- lapply(readLines(file.path(out, 'oracle.stdout')), jsonlite::fromJSON)
  stopifnot(identical(vapply(oracle, `[[`, '', 'id'), vapply(cases, `[[`, '', 'id')))
  expected <- do.call(rbind, lapply(oracle, function(x) {
    rows <- x$independent_hgvs
    stopifnot(all(events$case_id[match(as.integer(rows$id), events$event_index)] == x$id))
    data.frame(event_index = as.integer(rows$id), allele = rows$allele,
      hgvsp = sub('^.*:p\\.', 'p.', rows$hgvsp),
      so = vapply(rows$consequences, function(terms) paste(sort(terms), collapse = '&'), ''))
  }))
  stopifnot(setequal(expected$event_index, events$event_index))
  controls <- rbind(codon_controls(expected), codon_provenance_controls(FALSE),
    codon_provenance_controls(TRUE))
  write.csv(expected, file.path(out, 'expected.csv'), row.names = FALSE)
  write.csv(controls, file.path(out, 'controls.csv'), row.names = FALSE)
  con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = 'true')))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  q <- function(x) as.character(dbQuoteString(con, x))
  dbExecute(con, paste('LOAD', q(extension)))
  dbWriteTable(con, 'inputs', events)
  tx <- paste('SELECT transcript_index::UINTEGER transcript_index,seq_region::UINTEGER seq_region,',
    '1::UBIGINT transcript_start,9::UBIGINT transcript_end,1::TINYINT strand,',
    'transcript_index::UINTEGER gene_index,3::UBIGINT transcript_flags,',
    'transcript_start cds_start,transcript_end cds_end,cds::BLOB cds_sequence,',
    '"table"::UTINYINT codon_table FROM inputs ORDER BY transcript_index')
  ex <- paste('SELECT transcript_index::UINTEGER transcript_index,1::UBIGINT exon_start,',
    '9::UBIGINT exon_end,1::UBIGINT exon_cdna_start,9::UBIGINT exon_cdna_end,',
    '0::TINYINT phase,0::TINYINT end_phase FROM inputs ORDER BY transcript_index')
  stopifnot(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('codons',",
    q('SELECT seq_region::UINTEGER seq_region FROM inputs ORDER BY seq_region'), ',', q(tx), ',', q(ex), ')'))$loaded)
  dbExecute(con, paste('CREATE TABLE events AS SELECT event_index,seq_region,position,reference,alternate,',
    'NULL::UBIGINT end_position,NULL::VARCHAR structural_type,NULL::VARCHAR copy_change,',
    'NULL::UINTEGER mate_seq_region,NULL::UBIGINT mate_position FROM inputs ORDER BY seq_region,position'))
  comparisons <- list()
  for (threads in c(1L, 4L)) {
    dbExecute(con, paste('SET threads=', threads))
    label <- paste0('independent_', threads)
    dbExecute(con, paste0('CREATE TABLE ', label, " AS SELECT * FROM duckvep_annotate('events',",
      "'codons',hgvs:=true,upstream_distance:=0,downstream_distance:=0)"))
    dbExecute(con, paste('COPY', label, 'TO', q(file.path(out, paste0(label, '.parquet'))), '(FORMAT PARQUET)'))
    geometry <- dbGetQuery(con, paste('SELECT a.event_index,a.transcript_index FROM', label, 'a'))
    stopifnot(!anyNA(geometry$event_index),
      all(geometry$transcript_index == geometry$event_index - 1L))
    actual <- dbGetQuery(con, paste('SELECT a.event_index::INTEGER event_index,i.alternate allele,',
      'a.protein_hgvs hgvsp,(SELECT string_agg(t.consequence,\'&\' ORDER BY t.consequence)',
      'FROM duckvep_so_terms() t WHERE (a.consequence_mask & t.consequence_mask)<>0) so',
      'FROM', label, 'a LEFT JOIN inputs i USING(event_index)'))
    comparisons[[label]] <- codon_equal(actual, expected)
    for (route in c('strict', 'vep116_compat', 'source_records')) {
      label <- paste0(route, '_', threads)
      raw <- route == 'source_records'
      policy <- if (raw) 'vep116_compat' else route
      calls <- paste('SELECT event_index,seq_region,position,reference,alternate,transcript_index,',
        '1 alt_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set FROM inputs')
      if (raw) calls <- paste('SELECT event_index,seq_region,position,reference,transcript_index,',
        "0 sample_index,[alternate] alternates,'1|1' gt FROM inputs")
      dbExecute(con, paste0('CREATE TABLE ', label, ' AS SELECT * FROM duckvep_haplotypes(', q(calls),
        ",'codons',hgvs:=true,phase_policy:=", q(policy), ',input_mode:=',
        q(if (raw) 'source_records' else 'alt_events'), ')'))
      dbExecute(con, paste('COPY', label, 'TO', q(file.path(out, paste0(label, '.parquet'))), '(FORMAT PARQUET)'))
      provenance <- dbGetQuery(con, paste('SELECT transcript_index,contributors,carriers,carrier_count FROM', label))
      codon_check_provenance(provenance, events, raw)
      actual <- dbGetQuery(con, paste('SELECT h.contributors[1].event_index::INTEGER event_index,',
        'h.contributors[1].alternate allele,h.hgvsp,NULL::VARCHAR so FROM', label, 'h'))
      # Only the documented outer prediction wrapper is removed. No equality,
      # unknown-residue, absent-result or compound-expression normalization.
      actual$hgvsp <- sub('^p\\.\\((.*)\\)$', 'p.\\1', actual$hgvsp)
      phased_expected <- expected
      phased_expected$so <- NA_character_
      comparisons[[label]] <- codon_equal(actual, phased_expected)
      comparisons[[label]]$so_equal <- NA
    }
  }
  pairs <- do.call(rbind, lapply(names(comparisons), function(route)
    cbind(route = route, comparisons[[route]])))
  pairs <- merge(pairs, events[c('event_index', 'seq_region', 'transcript_index', 'position',
    'reference', 'codon', 'table', 'cds', 'case_id')],
    by = 'event_index', all.x = TRUE, sort = FALSE)
  pairs$source_n <- ifelse(is.na(pairs$reference), 'unmatched',
    ifelse(pairs$reference == 'N', 'yes', 'no'))
  pairs$codon_n <- ifelse(is.na(pairs$codon), 'unmatched',
    ifelse(grepl('N', pairs$codon, fixed = TRUE), 'yes', 'no'))
  write.csv(pairs, file.path(out, 'pairs.csv'), row.names = FALSE)
  summary <- aggregate(cbind(pairs = rep(1L, nrow(pairs)), hgvsp_failures = !pairs$hgvsp_equal,
    so_compared = !is.na(pairs$so_equal), so_failures = !pairs$so_equal),
    pairs[c('route', 'source_n', 'codon_n')], sum, na.rm = TRUE)
  write.csv(summary, file.path(out, 'summary.csv'), row.names = FALSE)
  # Arbitrary supplied binaries are explicitly diagnostic; hashes do not prove
  # that the checkout built the extension. Retain both identities independently.
  files <- list.files(out, full.names = TRUE)
  jsonlite::write_json(list(source_revision = revision, source_binding = 'diagnostic_unbound',
    extension_path = extension, extension_sha256 = extension_hash, oracle_revisions = pins,
    models = length(cases), source_snvs = nrow(events), oracle_pairs = nrow(expected),
    decoded_gt = 'haploid ALT', source_gt = '1|1',
    scope = 'independent_SO_HGVSp_and_singleton_phased_HGVSp_internal_codons',
    sha256 = as.list(c(source_hashes, module_hashes, vapply(files, duckvep_evidence_sha256, '')))),
    file.path(out, 'receipt.json'), pretty = TRUE, auto_unbox = TRUE)
  print(summary)
  stopifnot(identical(source_hashes, vapply(sources, duckvep_evidence_sha256, '')),
    identical(module_hashes, vapply(modules, duckvep_evidence_sha256, '')),
    identical(extension_hash, duckvep_evidence_sha256(extension)),
    identical(revision, duckvep_evidence_revision('.')))
  if (nzchar(opt$evidence_out)) {
    destination <- opt$evidence_out
    stopifnot(!dir.exists(destination), dir.create(destination))
    dbWriteTable(con, 'comparison_pairs', pairs)
    dbExecute(con, paste('COPY comparison_pairs TO', q(file.path(destination, 'pairs.parquet')),
      '(FORMAT PARQUET)'))
    for (name in c('cases.jsonl', 'oracle.stdout')) {
      compressed <- gzfile(file.path(destination, paste0(name, '.gz')), 'wt')
      writeLines(readLines(file.path(out, name)), compressed)
      close(compressed)
    }
    stopifnot(all(file.copy(file.path(out, c('summary.csv', 'controls.csv', 'environment.stdout')),
      destination)))
    manifest <- jsonlite::read_json(file.path(out, 'receipt.json'), simplifyVector = TRUE)
    manifest$local_receipt_sha256 <- duckvep_evidence_sha256(file.path(out, 'receipt.json'))
    manifest$scope <- paste(manifest$scope, 'complete_pairs_and_oracle_not_all_native_output_fields', sep = ';')
    manifest$source_sha256 <- as.list(c(source_hashes, module_hashes))
    published <- list.files(destination, full.names = TRUE)
    manifest$sha256 <- as.list(setNames(vapply(published, duckvep_evidence_sha256, ''), basename(published)))
    jsonlite::write_json(manifest, file.path(destination, 'receipt.json'), pretty = TRUE, auto_unbox = TRUE)
  }
  stopifnot(all(pairs$hgvsp_equal), all(pairs$so_equal, na.rm = TRUE))
}
if (sys.nframe() == 0L) main()
