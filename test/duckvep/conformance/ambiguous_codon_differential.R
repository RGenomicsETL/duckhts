#!/usr/bin/env Rscript
# Original VCF alleles: uploaded ambiguity and surrounding codon ambiguity.
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

codon_inputs <- function(family) {
  stopifnot(family %in% c('snv', 'indel'))
  bases <- strsplit('ACGTN', '', fixed = TRUE)[[1L]]
  triplets <- do.call(paste0, expand.grid(rep(list(bases), 3L), stringsAsFactors = FALSE))
  grid <- expand.grid(table = c(1:6, 9:14, 16, 21:31), codon = triplets,
    stringsAsFactors = FALSE)
  events <- vector('list', nrow(grid))
  cases <- vector('list', nrow(grid))
  event_count <- 0L
  for (i in seq_len(nrow(grid))) {
    cds <- paste0('ATG', grid$codon[i], if (family == 'snv') 'TAA' else 'GCCTAA')
    variants <- list()
    add <- function(position, reference, alternate) {
      variants[[length(variants) + 1L]] <<- list(
        id = as.character(event_count + length(variants) + 1L),
        position1 = position, reference = reference, alternate = alternate)
    }
    if (family == 'snv') {
      for (offset in 1:3) for (alt in bases[1:4]) {
        ref <- substr(grid$codon[i], offset, offset)
        if (ref != alt) add(offset + 3L, ref, alt)
      }
    } else {
      for (position in 3:6) {
        anchor <- substr(cds, position, position)
        for (inserted in c('A', 'C', 'G', 'T', 'AC', 'GCC', 'ACGT'))
          add(position, anchor, paste0(anchor, inserted))
        for (deleted in 1:3)
          add(position, substr(cds, position, position + deleted), anchor)
        for (replaced in 2:3) for (inserted in c('A', 'ACGT'))
          add(position, substr(cds, position, position + replaced), paste0(anchor, inserted))
      }
    }
    id <- paste(grid$table[i], grid$codon[i], sep = '/')
    cases[[i]] <- list(id = id, cds = cds, table = grid$table[i], edits = list(), variants = variants)
    indices <- event_count + seq_along(variants)
    events[[i]] <- data.frame(event_index = indices, seq_region = indices - 1L,
      transcript_index = indices - 1L, position = vapply(variants, `[[`, 0L, 'position1'),
      reference = vapply(variants, `[[`, '', 'reference'),
      alternate = vapply(variants, `[[`, '', 'alternate'), cds = cds,
      table = grid$table[i], codon = grid$codon[i], case_id = id)
    if (family == 'indel') {
      cases[[i]]$variant_format <- 'vcf'
      cases[[i]]$genomic_sequence <- paste0(strrep('A', 10L), cds, strrep('A', 10L))
      cases[[i]]$cds_start1 <- 11L
      events[[i]]$position <- events[[i]]$position + 10L
    }
    event_count <- event_count + length(variants)
  }
  events <- do.call(rbind, events)
  stopifnot(length(cases) == 3000L, !anyDuplicated(events$event_index),
    nrow(events) == if (family == 'snv') 28800L else 168000L)
  list(events = events, cases = cases)
}

codon_expected <- function(oracle, events, indel) {
  case_ids <- vapply(oracle, `[[`, '', 'id')
  stopifnot(!anyDuplicated(case_ids), setequal(case_ids, events$case_id))
  do.call(rbind, lapply(oracle, function(x) {
    # Scalar JSON decoding preserves literal DNA "NA" separately from JSON null.
    rows <- x$independent_hgvs
    required <- c('id', 'allele', 'hgvsp', 'consequences')
    if (indel) required <- c(required, 'source_reference', 'source_alternate',
      'parser_start', 'parser_end', 'parser_allele_string')
    stopifnot(!anyDuplicated(names(x)), is.list(rows), length(rows) > 0L)
    for (row in rows) {
      stopifnot(!anyDuplicated(names(row)), setequal(names(row), required),
        is.null(row$hgvsp) || (is.character(row$hgvsp) && length(row$hgvsp) == 1L &&
          !is.na(row$hgvsp)), is.list(row$consequences), length(row$consequences) > 0L,
        all(vapply(row$consequences, function(term)
          is.character(term) && length(term) == 1L && !is.na(term), TRUE)))
    }
    strings <- function(field) vapply(rows, function(row) row[[field]], '')
    stopifnot(all(grepl('^[1-9][0-9]*$', strings('id'))))
    indices <- as.integer(strings('id'))
    selected <- match(indices, events$event_index)
    stopifnot(!anyNA(selected), !anyDuplicated(indices),
      identical(strings('id'), as.character(indices)),
      setequal(indices, events$event_index[events$case_id == x$id]),
      all(events$case_id[selected] == x$id),
      identical(x$prepared_cds, unique(events$cds[selected])))
    if (indel) {
      reference <- strings('source_reference')
      alternate <- strings('source_alternate')
      stopifnot(identical(reference, events$reference[selected]),
        identical(alternate, events$alternate[selected]))
      for (j in seq_along(rows)) {
        row <- rows[[j]]
        source <- events[selected[j], ]
        parsed <- strsplit(row$parser_allele_string, '/', fixed = TRUE)[[1L]]
        stopifnot(length(parsed) == 2L, all(grepl('^([ACGTN]+|-)$', parsed)),
          identical(row$allele, parsed[2L]))
        parsed[parsed == '-'] <- ''
        first <- row$parser_start
        last <- row$parser_end
        stopifnot(is.numeric(first), is.numeric(last), length(first) == 1L, length(last) == 1L,
          is.finite(first), is.finite(last), first == floor(first), last == floor(last),
          first >= source$position + 1L,
          last <= source$position + nchar(source$reference) - 1L,
          first <= last + 1L, nchar(parsed[1L]) == last - first + 1L)
        # VEP also minimises some complex indels without --minimal. Check the
        # parsed span against the source sequence without duplicating its rules.
        genome <- paste0(strrep('A', 10L), source$cds, strrep('A', 10L))
        stopifnot(identical(substr(genome, first, last), parsed[1L]))
        original <- paste0(substr(genome, 1L, source$position - 1L), source$alternate,
          substring(genome, source$position + nchar(source$reference)))
        observed <- paste0(substr(genome, 1L, first - 1L), parsed[2L], substring(genome, last + 1L))
        stopifnot(identical(original, observed))
      }
    }
    data.frame(event_index = indices, allele = if (indel) alternate else strings('allele'),
      hgvsp = sub('^.*:p\\.', 'p.', vapply(rows, function(row)
        if (is.null(row$hgvsp)) NA_character_ else row$hgvsp, '')),
      so = vapply(rows, function(row) paste(sort(unlist(row$consequences)), collapse = '&'), ''))
  }))
}

codon_partitions <- function(events, indel) {
  stopifnot(nrow(events) > 0L, identical(events$event_index, seq_len(nrow(events))),
    identical(events$seq_region, events$event_index - 1L),
    identical(events$transcript_index, events$event_index - 1L))
  # Each original event has its own native transcript/region. The model stores
  # region ordinals in 16 bits; partitioning changes ownership, not eligibility.
  capacity <- if (indel) 56000L else nrow(events)
  offsets <- seq.int(0L, nrow(events) - 1L, by = capacity)
  data.frame(model_partition = seq_along(offsets), event_offset = offsets,
    first_event_index = offsets + 1L,
    last_event_index = pmin(offsets + capacity, nrow(events)),
    source_events = pmin(capacity, nrow(events) - offsets))
}

codon_partition_events <- function(events, partition) {
  stopifnot(nrow(partition) == 1L,
    partition$first_event_index == partition$event_offset + 1L,
    partition$last_event_index == partition$event_offset + partition$source_events)
  selected <- events[events$event_index >= partition$first_event_index &
    events$event_index <= partition$last_event_index, , drop = FALSE]
  stopifnot(nrow(selected) == partition$source_events,
    identical(selected$event_index, seq.int(partition$first_event_index,
      partition$last_event_index)))
  selected$seq_region <- selected$event_index - 1L - partition$event_offset
  selected$transcript_index <- selected$seq_region
  selected
}

codon_native_partition <- function(con, events, expected, partition, indel, fasta, out) {
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbWriteTable(con, 'inputs', events)
  start <- if (indel) 11L else 1L
  tx <- paste('SELECT transcript_index::UINTEGER transcript_index,seq_region::UINTEGER seq_region,',
    start, '::UBIGINT transcript_start,', start - 1L,
    '+length(cds)::UBIGINT transcript_end,1::TINYINT strand,',
    'transcript_index::UINTEGER gene_index,3::UBIGINT transcript_flags,',
    'transcript_start cds_start,transcript_end cds_end,cds::BLOB cds_sequence,',
    '"table"::UTINYINT codon_table',
    if (indel) ",''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence" else '',
    'FROM inputs ORDER BY transcript_index')
  ex <- paste('SELECT transcript_index::UINTEGER transcript_index,', start, '::UBIGINT exon_start,',
    start - 1L, '+length(cds)::UBIGINT exon_end,1::UBIGINT exon_cdna_start,',
    'length(cds)::UBIGINT exon_cdna_end,',
    '0::TINYINT phase,0::TINYINT end_phase FROM inputs ORDER BY transcript_index')
  regions <- 'SELECT seq_region::UINTEGER seq_region FROM inputs ORDER BY seq_region'
  reference_option <- ''
  if (indel) {
    regions <- paste('SELECT seq_region::UINTEGER seq_region,',
      '(length(cds)+20)::UBIGINT sequence_length,event_index::VARCHAR seq_region_name',
      'FROM inputs ORDER BY seq_region')
    reference_option <- paste0(',reference_fasta:=', q(fasta))
  }
  stopifnot(DBI::dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('codons',",
    q(regions), ',', q(tx), ',', q(ex), reference_option, ')'))$loaded)
  DBI::dbExecute(con, paste('CREATE TABLE events AS SELECT event_index,seq_region,position,reference,alternate,',
    'NULL::UBIGINT end_position,NULL::VARCHAR structural_type,NULL::VARCHAR copy_change,',
    'NULL::UINTEGER mate_seq_region,NULL::UBIGINT mate_position FROM inputs ORDER BY seq_region,position'))
  comparisons <- list()
  artifact <- function(label) file.path(out, paste0(label,
    if (indel) paste0('_partition_', partition$model_partition) else '', '.parquet'))
  for (threads in c(1L, 4L)) {
    DBI::dbExecute(con, paste('SET threads=', threads))
    label <- paste0('independent_', threads)
    message('Native partition ', partition$model_partition, ' route: ', label)
    DBI::dbExecute(con, paste0('CREATE TABLE ', label, " AS SELECT * FROM duckvep_annotate('events',",
      "'codons',hgvs:=true,upstream_distance:=0,downstream_distance:=0)"))
    DBI::dbExecute(con, paste('COPY', label, 'TO', q(artifact(label)), '(FORMAT PARQUET)'))
    geometry <- DBI::dbGetQuery(con, paste('SELECT a.event_index,a.transcript_index FROM', label, 'a'))
    stopifnot(!anyNA(geometry$event_index),
      all(geometry$transcript_index == geometry$event_index - 1L - partition$event_offset))
    actual <- DBI::dbGetQuery(con, paste('SELECT a.event_index::INTEGER event_index,i.alternate allele,',
      'a.protein_hgvs hgvsp,(SELECT string_agg(t.consequence,\'&\' ORDER BY t.consequence)',
      'FROM duckvep_so_terms() t WHERE (a.consequence_mask & t.consequence_mask)<>0) so',
      'FROM', label, 'a LEFT JOIN inputs i USING(event_index)'))
    comparisons[[label]] <- codon_equal(actual, expected)
    DBI::dbRemoveTable(con, label)
    for (route in c('strict', 'vep116_compat', 'source_records')) {
      label <- paste0(route, '_', threads)
      message('Native partition ', partition$model_partition, ' route: ', label)
      raw <- route == 'source_records'
      policy <- if (raw) 'vep116_compat' else route
      calls <- paste('SELECT event_index,seq_region,position,reference,alternate,transcript_index,',
        '1 alt_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set FROM inputs')
      if (raw) calls <- paste('SELECT event_index,seq_region,position,reference,transcript_index,',
        "0 sample_index,[alternate] alternates,'1|1' gt FROM inputs")
      DBI::dbExecute(con, paste0('CREATE TABLE ', label, ' AS SELECT * FROM duckvep_haplotypes(', q(calls),
        ",'codons',hgvs:=true,phase_policy:=", q(policy), ',input_mode:=',
        q(if (raw) 'source_records' else 'alt_events'), ')'))
      DBI::dbExecute(con, paste('COPY', label, 'TO', q(artifact(label)), '(FORMAT PARQUET)'))
      provenance <- DBI::dbGetQuery(con, paste('SELECT transcript_index,contributors,carriers,carrier_count FROM', label))
      codon_check_provenance(provenance, events, raw)
      actual <- DBI::dbGetQuery(con, paste('SELECT h.contributors[1].event_index::INTEGER event_index,',
        'h.contributors[1].alternate allele,h.hgvsp,NULL::VARCHAR so FROM', label, 'h'))
      # Only the documented outer prediction wrapper is removed.
      actual$hgvsp <- sub('^p\\.\\((.*)\\)$', 'p.\\1', actual$hgvsp)
      phased_expected <- expected
      phased_expected$so <- NA_character_
      comparisons[[label]] <- codon_equal(actual, phased_expected)
      comparisons[[label]]$so_equal <- NA
      DBI::dbRemoveTable(con, label)
    }
  }
  # Each CREATE TABLE above has consumed its scan to EOF. No live result borrows
  # the model when it is dropped; failure cleanup belongs to the owned connection.
  stopifnot(DBI::dbGetQuery(con, "SELECT duckvep_model_drop('codons') dropped")$dropped)
  DBI::dbRemoveTable(con, 'events')
  DBI::dbRemoveTable(con, 'inputs')
  comparisons
}

main <- function() {
  suppressPackageStartupMessages(library(DBI))
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option('--extension', default = 'build/release/duckhts.duckdb_extension'),
    optparse::make_option('--evidence-out', dest = 'evidence_out', default = '',
      help = 'Retain a compact complete comparison bundle in a new directory'),
    optparse::make_option('--variant-family', dest = 'variant_family', default = 'snv',
      help = 'snv or indel: complete declared original-allele matrix'),
    optparse::make_option('--vep-prefix', dest = 'vep_prefix',
      default = Sys.getenv('VEP_PREFIX', '/root/miniconda3/envs/vep'))
  )))
  stopifnot(opt$variant_family %in% c('snv', 'indel'))
  indel <- opt$variant_family == 'indel'
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
  out <- tempfile(if (indel) 'ambiguous_indel_' else 'ambiguous_codon_',
    'test/duckvep/conformance/results')
  stopifnot(dir.create(out))
  out <- normalizePath(out)
  message('Ambiguous-codon artifacts: ', out)
  execution_complete <- FALSE
  stage <- 'oracle_environment'
  if (indel) on.exit({
    if (!execution_complete) tryCatch({
      retained <- list.files(out, full.names = TRUE)
      failure <- list(execution_status = 'failed', last_started_stage = stage,
        source_binding = 'diagnostic_unbound', source_revision = revision,
        extension_path = extension, extension_sha256 = extension_hash,
        oracle_revisions = pins, declared_oracle_models = 3000L,
        declared_source_indels = 168000L, source_sha256 = as.list(source_hashes),
        sha256 = as.list(setNames(vapply(retained, duckvep_evidence_sha256, ''), basename(retained))))
      jsonlite::write_json(failure, file.path(out, 'failure_receipt.json'),
        pretty = TRUE, auto_unbox = TRUE)
    }, error = function(e) warning('Could not retain failure receipt: ', conditionMessage(e)))
  }, add = TRUE)
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
  if (indel) modules <- c(modules,
    file.path(mirrors[['vep']], 'modules/Bio/EnsEMBL/VEP', c('Config.pm', 'Parser.pm', 'Parser/VCF.pm')),
    file.path(mirrors[['variation']], 'modules/Bio/EnsEMBL/Variation/Utils/FastaSequence.pm'),
    file.path(prefix, 'share/ensembl-vep-116.0-0/Bio/EnsEMBL/Slice.pm'))
  module_hashes <- vapply(modules, duckvep_evidence_sha256, '')
  inputs <- codon_inputs(opt$variant_family)
  events <- inputs$events
  cases <- inputs$cases
  saveRDS(events, file.path(out, 'events.rds'))
  input <- file.path(out, 'cases.jsonl')
  writeLines(vapply(cases, jsonlite::toJSON, '', auto_unbox = TRUE), input)
  libs <- paste(c(file.path(mirrors, 'modules'),
    file.path(prefix, 'share/ensembl-vep-116.0-0')), collapse = ':')
  message('VEP oracle: ', length(cases), ' models / ', nrow(events), ' original alleles')
  stage <- 'oracle_execution'
  command(c('run', '--clean-env', '--env', paste0('PERL5LIB=', libs), '-p', prefix, 'perl',
    normalizePath('test/duckvep/conformance/reference_translation_oracle.pl'), input), 'oracle')
  stage <- 'oracle_validation'
  oracle <- lapply(readLines(file.path(out, 'oracle.stdout')), jsonlite::fromJSON,
    simplifyVector = FALSE)
  stopifnot(identical(vapply(oracle, `[[`, '', 'id'), vapply(cases, `[[`, '', 'id')))
  expected <- codon_expected(oracle, events, indel)
  stopifnot(setequal(expected$event_index, events$event_index))
  controls <- rbind(codon_controls(expected), codon_provenance_controls(FALSE),
    codon_provenance_controls(TRUE))
  write.csv(expected, file.path(out, 'expected.csv'), row.names = FALSE)
  write.csv(controls, file.path(out, 'controls.csv'), row.names = FALSE)
  con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = 'true')))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  q <- function(x) as.character(dbQuoteString(con, x))
  dbExecute(con, paste('LOAD', q(extension)))
  stage <- 'native_reference_preparation'
  fasta <- NULL
  if (indel) {
    fasta <- file.path(out, 'reference.fa')
    writeLines(as.vector(rbind(paste0('>', events$event_index),
      paste0(strrep('A', 10L), events$cds, strrep('A', 10L)))), fasta)
    stopifnot(system2('samtools', c('faidx', shQuote(fasta)),
      stdout = file.path(out, 'faidx.stdout'), stderr = file.path(out, 'faidx.stderr')) == 0L)
  }
  partitions <- codon_partitions(events, indel)
  if (indel) write.csv(partitions, file.path(out, 'native_partitions.csv'), row.names = FALSE)
  comparisons <- vector('list', nrow(partitions))
  source_facts <- vector('list', nrow(partitions))
  for (i in seq_len(nrow(partitions))) {
    stage <- paste0('native_partition_', i)
    partition <- partitions[i, , drop = FALSE]
    local_events <- codon_partition_events(events, partition)
    local_expected <- expected[expected$event_index %in% local_events$event_index, , drop = FALSE]
    comparisons[[i]] <- codon_native_partition(con, local_events, local_expected,
      partition, indel, fasta, out)
    source_facts[[i]] <- local_events[c('event_index', 'seq_region', 'transcript_index', 'position',
      'reference', 'codon', 'table', 'cds', 'case_id')]
    if (indel) {
      source_facts[[i]]$model_partition <- partition$model_partition
      source_facts[[i]]$event_offset <- partition$event_offset
    }
  }
  pairs <- do.call(rbind, lapply(comparisons, function(part)
    do.call(rbind, lapply(names(part), function(route) cbind(route = route, part[[route]])))))
  pairs <- merge(pairs, do.call(rbind, source_facts), by = 'event_index', all.x = TRUE, sort = FALSE)
  pairs$source_n <- ifelse(is.na(pairs$reference), 'unmatched',
    ifelse(grepl('N', pairs$reference, fixed = TRUE), 'yes', 'no'))
  pairs$codon_n <- ifelse(is.na(pairs$codon), 'unmatched',
    ifelse(grepl('N', pairs$codon, fixed = TRUE), 'yes', 'no'))
  write.csv(pairs, file.path(out, 'pairs.csv'), row.names = FALSE)
  summary <- aggregate(cbind(pairs = rep(1L, nrow(pairs)), hgvsp_failures = !pairs$hgvsp_equal,
    so_compared = !is.na(pairs$so_equal), so_failures = !pairs$so_equal),
    pairs[c('route', 'source_n', 'codon_n')], sum, na.rm = TRUE)
  write.csv(summary, file.path(out, 'summary.csv'), row.names = FALSE)
  stage <- 'complete_evidence_retention'
  # Arbitrary supplied binaries are explicitly diagnostic; hashes do not prove
  # that the checkout built the extension. Retain both identities independently.
  files <- list.files(out, full.names = TRUE)
  manifest <- list(source_revision = revision, source_binding = 'diagnostic_unbound',
    extension_path = extension, extension_sha256 = extension_hash, oracle_revisions = pins,
    models = length(cases), source_snvs = nrow(events), oracle_pairs = nrow(expected),
    decoded_gt = 'haploid ALT', source_gt = '1|1',
    scope = 'independent_SO_HGVSp_and_singleton_phased_HGVSp_internal_codons',
    sha256 = as.list(c(source_hashes, module_hashes, vapply(files, duckvep_evidence_sha256, ''))))
  if (indel) {
    manifest$models <- NULL
    manifest$oracle_models <- length(cases)
    manifest$native_transcript_instances <- nrow(events)
    manifest$native_transcript_flanks <- 'complete_empty_pre_CDS_and_post_CDS_matching_oracle_transcript'
    manifest$native_partitions <- partitions
    manifest$native_ordinal_contract <- paste(
      'seq_region and transcript_index are local to model_partition;',
      'event_index equals local ordinal plus event_offset plus one')
    manifest$hgvsp_comparisons <- nrow(pairs)
    manifest$so_comparisons <- sum(!is.na(pairs$so_equal))
    manifest$source_snvs <- NULL
    manifest$source_indels <- nrow(events)
    manifest$variant_family <- 'indel'
    manifest$scope <- 'independent_SO_HGVSp_and_singleton_phased_HGVSp_original_VCF_indels_forward_single_exon'
  }
  jsonlite::write_json(manifest,
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
    compressed_files <- c('cases.jsonl', 'oracle.stdout', if (indel) 'oracle.stderr')
    for (name in compressed_files) {
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
  execution_complete <- TRUE
  stopifnot(all(pairs$hgvsp_equal), all(pairs$so_equal, na.rm = TRUE))
}
if (sys.nframe() == 0L) main()
