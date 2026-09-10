#!/usr/bin/env Rscript
# Exact MNV compatibility and separate cis-SNV/MNV representation comparisons.
source('scripts/duckvep_evidence.R')

cis_codon_reverse_complement <- function(x) {
  paste0(rev(strsplit(chartr('ACGT', 'TGCA', x), '', fixed = TRUE)[[1L]]), collapse = '')
}

cis_codon_inputs <- function() {
  rc <- cis_codon_reverse_complement
  codons <- do.call(paste0, expand.grid(rep(list(c('A', 'C', 'G', 'T')), 3L)))
  models <- expand.grid(ref = setdiff(codons, c('TAA', 'TAG', 'TGA')), alt = codons,
    strand = c(1L, -1L), stringsAsFactors = FALSE)
  models$edits <- mapply(function(r, a)
    sum(strsplit(r, '')[[1L]] != strsplit(a, '')[[1L]]), models$ref, models$alt)
  models <- models[models$edits >= 2L, ]
  models$transcript_index <- models$seq_region <- seq_len(nrow(models)) - 1L
  models$id <- sprintf('CIS%05d', models$transcript_index)
  models$chrom <- paste0('chr', models$id)
  models$cds <- paste0('ATG', models$ref, 'GCCTAA')
  models$expected_cds <- paste0('ATG', models$alt, 'GCCTAA')
  models$genome <- paste0(strrep('A', 10L), ifelse(models$strand == 1L, models$cds,
    vapply(models$cds, rc, '')), strrep('A', 10L))
  mnv <- data.frame(transcript_index = models$transcript_index, seq_region = models$seq_region,
    position = ifelse(models$strand == 1L, 14L, 17L),
    reference = ifelse(models$strand == 1L, models$ref, vapply(models$ref, rc, '')),
    alternate = ifelse(models$strand == 1L, models$alt, vapply(models$alt, rc, '')),
    event_index = seq_len(nrow(models)))
  snv <- do.call(rbind, lapply(seq_len(nrow(models)), function(i) {
    m <- models[i, ]
    ref <- strsplit(m$ref, '')[[1L]]
    alt <- strsplit(m$alt, '')[[1L]]
    changed <- which(ref != alt)
    cds_position1 <- changed + 3L
    data.frame(transcript_index = m$transcript_index, seq_region = m$seq_region,
      position = if (m$strand == 1L) 10L + cds_position1 else 23L - cds_position1,
      reference = if (m$strand == 1L) ref[changed] else vapply(ref[changed], rc, ''),
      alternate = if (m$strand == 1L) alt[changed] else vapply(alt[changed], rc, ''))
  }))
  snv$event_index <- seq_len(nrow(snv))
  # Rebuild both physical input representations on genomic DNA independently
  # of DuckVEP's projection, edit grouping and transcript replay.
  for (records in list(mnv, snv)) {
    by_transcript <- split(records, records$transcript_index)
    for (i in seq_len(nrow(models))) {
      genome <- strsplit(models$genome[i], '', fixed = TRUE)[[1L]]
      events <- by_transcript[[as.character(models$transcript_index[i])]]
      for (j in seq_len(nrow(events))) {
        event <- events[j, ]
        at <- event$position + seq_len(nchar(event$reference)) - 1L
        stopifnot(paste0(genome[at], collapse = '') == event$reference)
        genome[at] <- strsplit(event$alternate, '', fixed = TRUE)[[1L]]
      }
      cds <- paste0(genome[11:22], collapse = '')
      if (models$strand[i] == -1L) cds <- rc(cds)
      stopifnot(cds == models$expected_cds[i])
    }
  }
  stopifnot(nrow(models) == 6588L, nrow(snv) == 16470L,
    all(table(models$edits, models$strand) == 1647L))
  list(models = models, mnv = mnv, snv = snv)
}

cis_codon_oracle <- function(lines, inputs) {
  models <- inputs$models
  mnv <- inputs$mnv
  stopifnot(any(startsWith(lines, '##VEP="v116.0" API="v116"')))
  header <- lines[startsWith(lines, '##INFO=<ID=CSQ,')]
  stopifnot(length(header) == 1L)
  columns <- strsplit(sub('.*Format: ([^"]+).*', '\\1', header), '|', fixed = TRUE)[[1L]]
  stopifnot(!anyDuplicated(columns), all(c('Feature', 'HGVSp', 'Amino_acids') %in% columns))
  rows <- strsplit(lines[!startsWith(lines, '#')], '\t', fixed = TRUE)
  stopifnot(length(rows) == nrow(models), all(lengths(rows) == 8L))
  fields <- do.call(rbind, lapply(rows, function(x) x[1:7]))
  expected <- cbind(models$chrom, as.character(mnv$position), models$id,
    mnv$reference, mnv$alternate, '.', 'PASS')
  stopifnot(identical(unname(fields), unname(expected)))
  csq <- do.call(rbind, lapply(rows, function(x) {
    stopifnot(startsWith(x[8L], 'CSQ='), !grepl(',', x[8L], fixed = TRUE))
    values <- strsplit(sub('^CSQ=', '', x[8L]), '|', fixed = TRUE)[[1L]]
    stopifnot(length(values) <= length(columns))
    length(values) <- length(columns)
    values[!is.na(values) & values == ''] <- NA_character_
    values[!is.na(values)] <- vapply(values[!is.na(values)], utils::URLdecode, '')
    setNames(values, columns)
  }))
  stopifnot(identical(unname(csq[, 'Feature']), models$id))
  # Only the internal second codon differs; VEP's final amino-acid allele
  # determines the displayed protein, including a new stop at that codon.
  alternate_aa <- sub('^.*/', '', csq[, 'Amino_acids'])
  stopifnot(!anyNA(alternate_aa), all(nchar(alternate_aa) == 1L))
  data.frame(id = models$id, hgvsp = sub('^.*:p\\.', 'p.', csq[, 'HGVSp']),
    protein = ifelse(alternate_aa == '*', 'M*', paste0('M', alternate_aa, 'A*')))
}

cis_codon_equal <- function(actual, expected) {
  stopifnot(all(c('id', 'hgvsp') %in% names(actual)),
    all(c('id', 'hgvsp') %in% names(expected)),
    nrow(actual) == nrow(expected), !anyNA(actual$id), !anyNA(expected$id),
    !anyDuplicated(actual$id), !anyDuplicated(expected$id), setequal(actual$id, expected$id))
  wanted <- expected$hgvsp[match(actual$id, expected$id)]
  equal <- (is.na(actual$hgvsp) & is.na(wanted)) |
    (!is.na(actual$hgvsp) & !is.na(wanted) & actual$hgvsp == wanted)
  equal[is.na(equal)] <- FALSE
  equal
}

cis_codon_sequence_equal <- function(actual, expected) {
  stopifnot(length(actual) == length(expected), !anyNA(expected))
  !is.na(actual) & actual == expected
}

cis_codon_check_native <- function(actual, records, models, ploidy = 1L) {
  stopifnot(nrow(actual) == nrow(models),
    all(c('transcript_index', 'carrier_count', 'sequence_status', 'projection_status',
      'contributors', 'carriers') %in% names(actual)),
    identical(as.numeric(actual$transcript_index), as.numeric(models$transcript_index)),
    all(actual$carrier_count == ploidy), all(actual$sequence_status == 'ok'),
    all(actual$projection_status == 'ok'))
  grouped <- split(records, records$transcript_index)
  for (i in seq_len(nrow(models))) {
    source <- grouped[[as.character(models$transcript_index[i])]]
    observed <- actual$contributors[[i]]
    fields <- c('event_index', 'seq_region', 'position', 'reference', 'alternate')
    stopifnot(all(c(fields, 'evidence_flags', 'projection_status') %in% names(observed)),
      nrow(observed) == nrow(source), !anyDuplicated(observed$event_index),
      setequal(observed$event_index, source$event_index))
    source <- source[match(observed$event_index, source$event_index), ]
    for (field in fields)
      stopifnot(all(observed[[field]] == source[[field]]))
    carrier <- actual$carriers[[i]]
    stopifnot(all(c('sample_index', 'phase_set', 'haplotype_lane', 'ploidy') %in% names(carrier)),
      all(observed$evidence_flags == 1L), all(observed$projection_status == 'ok'),
      nrow(carrier) == ploidy, all(carrier$sample_index == 0L), all(is.na(carrier$phase_set)),
      identical(sort(as.integer(carrier$haplotype_lane)), seq_len(ploidy)),
      all(carrier$ploidy == ploidy))
  }
  invisible(TRUE)
}

cis_codon_controls <- function(actual, records, models, expected) {
  rejected <- function(expr) {
    !isTRUE(tryCatch(force(expr), error = function(e) FALSE))
  }
  stopifnot(cis_codon_check_native(actual, records, models),
    all(cis_codon_equal(expected, expected)))
  controls <- c(dropped_hgvs = rejected(all(cis_codon_equal(expected[-1L, ], expected))))
  for (field in c('id', 'hgvsp')) {
    changed <- expected
    changed[[field]][1L] <- 'deliberate_corruption'
    controls[paste0('changed_', field)] <- rejected(all(cis_codon_equal(changed, expected)))
    changed[[field]] <- NULL
    controls[paste0('absent_', field)] <- rejected(all(cis_codon_equal(changed, expected)))
  }
  changed <- expected
  changed$id[2L] <- changed$id[1L]
  controls['duplicate_hgvs_id'] <- rejected(all(cis_codon_equal(changed, expected)))
  at <- which(!is.na(expected$hgvsp))[1L]
  stopifnot(!is.na(at))
  changed <- expected
  changed$hgvsp[at] <- NA_character_
  controls['missing_hgvs'] <- rejected(all(cis_codon_equal(changed, expected)))
  absent <- changed
  controls['invented_hgvs'] <- rejected(all(cis_codon_equal(expected, absent)))
  controls['dropped_native'] <- rejected(cis_codon_check_native(actual[-1L, ], records, models))
  changed <- actual
  changed$contributors[[1L]] <- changed$contributors[[1L]][-1L, ]
  controls['dropped_contributor'] <- rejected(cis_codon_check_native(changed, records, models))
  for (field in c('event_index', 'seq_region', 'position', 'reference', 'alternate',
                  'evidence_flags', 'projection_status')) {
    changed <- actual
    value <- changed$contributors[[1L]][[field]]
    value[1L] <- if (is.character(value)) 'corrupt' else value[1L] + 1L
    changed$contributors[[1L]][[field]] <- value
    controls[paste0('contributor_', field)] <- rejected(
      cis_codon_check_native(changed, records, models))
    changed$contributors[[1L]][[field]] <- NULL
    controls[paste0('absent_contributor_', field)] <- rejected(
      cis_codon_check_native(changed, records, models))
  }
  for (field in c('sample_index', 'phase_set', 'haplotype_lane', 'ploidy')) {
    changed <- actual
    changed$carriers[[1L]][[field]][1L] <- 2L
    controls[paste0('carrier_', field)] <- rejected(
      cis_codon_check_native(changed, records, models))
    changed$carriers[[1L]][[field]] <- NULL
    controls[paste0('absent_carrier_', field)] <- rejected(
      cis_codon_check_native(changed, records, models))
  }
  controls['missing_sequence'] <- rejected(all(cis_codon_sequence_equal(NA_character_, 'ATG')))
  controls['changed_sequence'] <- rejected(all(cis_codon_sequence_equal('ATA', 'ATG')))
  stopifnot(all(controls))
  controls
}

cis_codon_main <- function() {
  suppressPackageStartupMessages({ library(DBI); library(duckdb) })
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option('--extension', default = 'build/release/duckhts.duckdb_extension'),
    optparse::make_option('--vep-prefix', dest = 'vep_prefix', default = Sys.getenv('VEP_PREFIX')),
    optparse::make_option('--out', default = '')
  )))
  extension <- normalizePath(opt$extension, mustWork = TRUE)
  stopifnot(nzchar(opt$vep_prefix))
  prefix <- normalizePath(opt$vep_prefix, mustWork = TRUE)
  pins <- c(vep = '57ea5c52340acc1f156267f810ad162e26597082',
    variation = '2fb834b987ede3824e200197a838ce11e91aeb4b')
  for (name in names(pins)) {
    mirror <- file.path('.sync', paste0('ensembl-', name))
    stopifnot(identical(duckvep_evidence_command('git', c('-C', mirror, 'rev-parse', 'HEAD'),
      'oracle revision'), unname(pins[name])), !length(duckvep_evidence_command('git',
      c('-C', mirror, 'status', '--porcelain'), 'oracle state')))
  }
  out <- if (nzchar(opt$out)) opt$out else tempfile('hgvs_cis_codon_',
    'test/duckvep/conformance/results')
  stopifnot(!dir.exists(out), dir.create(out, recursive = TRUE))
  out <- normalizePath(out)
  message('Cis-codon artifacts: ', out)
  sha <- duckvep_evidence_sha256
  revision <- duckvep_evidence_revision('.')
  sources <- unique(c('test/duckvep/conformance/hgvs_cis_codon_differential.R',
    'scripts/duckvep_evidence.R', list.files('src/duckvep', recursive = TRUE,
      full.names = TRUE, pattern = '\\.[ch]$')))
  hashes <- vapply(sources, sha, '')
  extension_hash <- sha(extension)
  stopifnot(file.copy(extension, file.path(out, 'duckhts.duckdb_extension')))
  command <- function(exe, args, label) {
    status <- system2(exe, shQuote(args), stdout = file.path(out, paste0(label, '.stdout')),
      stderr = file.path(out, paste0(label, '.stderr')))
    stopifnot(status == 0L)
  }
  command('micromamba', c('list', '--explicit', '-p', prefix), 'environment')
  stopifnot(identical(duckvep_evidence_explicit_packages(readLines(file.path(out,
    'environment.stdout'))), duckvep_evidence_explicit_packages(readLines(
      'test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt'))))
  inputs <- cis_codon_inputs()
  models <- inputs$models
  saveRDS(inputs, file.path(out, 'inputs.rds'))
  writeLines(as.vector(rbind(paste0('>', models$chrom), models$genome)),
    file.path(out, 'reference.fa'))
  gff <- lapply(seq_len(nrow(models)), function(i) {
    id <- models$id[i]
    attrs <- c(paste0('ID=gene:', id, ';biotype=protein_coding'),
      paste0('ID=transcript:', id, ';Parent=gene:', id, ';biotype=protein_coding'),
      paste0('ID=exon:', id, ';Parent=transcript:', id), paste0('Parent=transcript:', id))
    paste(models$chrom[i], 'cis', c('gene', 'mRNA', 'exon', 'CDS'), 11, 22, '.',
      if (models$strand[i] == 1L) '+' else '-', c('.', '.', '.', '0'), attrs, sep = '\t')
  })
  writeLines(c('##gff-version 3', unlist(gff)), file.path(out, 'model.gff3'))
  writeLines(c('##fileformat=VCFv4.4', paste0('##contig=<ID=', models$chrom, ',length=32>'),
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    with(inputs$mnv, paste(models$chrom, position, models$id, reference, alternate,
      '.', 'PASS', '.', sep = '\t'))), file.path(out, 'input.vcf'))
  command('samtools', c('faidx', file.path(out, 'reference.fa')), 'faidx')
  command('bgzip', file.path(out, 'model.gff3'), 'bgzip')
  command('tabix', c('-p', 'gff', file.path(out, 'model.gff3.gz')), 'tabix')
  libs <- paste(normalizePath(c('.sync/ensembl-vep/modules', '.sync/ensembl-variation/modules',
    file.path(prefix, 'share/ensembl-vep-116.0-0'))), collapse = ':')
  oracle <- list()
  for (buffer in c(1L, 5000L)) {
    label <- paste0('vep', buffer)
    command('micromamba', c('run', '--clean-env', '--env', paste0('PERL5LIB=', libs), '-p', prefix,
      'perl', normalizePath('.sync/ensembl-vep/vep'), '--gff', file.path(out, 'model.gff3.gz'),
      '--fasta', file.path(out, 'reference.fa'), '--format', 'vcf', '--vcf', '--hgvs',
      '--buffer_size', buffer, '--no_stats', '--force_overwrite',
      '--input_file', file.path(out, 'input.vcf'),
      '--output_file', file.path(out, paste0(label, '.vcf'))), label)
    oracle[[label]] <- cis_codon_oracle(readLines(file.path(out, paste0(label, '.vcf'))), inputs)
  }
  stopifnot(identical(oracle[[1L]], oracle[[2L]]))
  expected <- oracle[[1L]]
  saveRDS(oracle, file.path(out, 'oracle.rds'))
  con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = 'true')))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  q <- function(x) as.character(dbQuoteString(con, x))
  dbExecute(con, paste('LOAD', q(file.path(out, 'duckhts.duckdb_extension'))))
  dbWriteTable(con, 'models', models)
  for (name in c('mnv', 'snv')) dbWriteTable(con, name, inputs[[name]])
  tx <- paste('SELECT transcript_index::UINTEGER transcript_index,seq_region::UINTEGER seq_region,',
    '11::UBIGINT transcript_start,22::UBIGINT transcript_end,strand::TINYINT strand,',
    'transcript_index::UINTEGER gene_index,3::UBIGINT transcript_flags,',
    '11::UBIGINT cds_start,22::UBIGINT cds_end,cds::BLOB cds_sequence,',
    "1::UTINYINT codon_table,''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence",
    'FROM models ORDER BY transcript_index')
  ex <- paste('SELECT transcript_index::UINTEGER transcript_index,11::UBIGINT exon_start,',
    '22::UBIGINT exon_end,1::UBIGINT exon_cdna_start,12::UBIGINT exon_cdna_end,',
    '0::TINYINT phase,0::TINYINT end_phase FROM models ORDER BY transcript_index')
  stopifnot(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('cis',",
    q('SELECT seq_region::UINTEGER seq_region,32::UBIGINT sequence_length,chrom seq_region_name FROM models ORDER BY seq_region'),
    ',', q(tx), ',', q(ex), ',reference_fasta:=', q(file.path(out, 'reference.fa')), ')'))$loaded)
  dbExecute(con, paste('CREATE TABLE independent_input AS SELECT event_index,seq_region,position,',
    'reference,alternate,NULL::UBIGINT end_position,NULL::VARCHAR structural_type,',
    'NULL::VARCHAR copy_change,NULL::UINTEGER mate_seq_region,NULL::UBIGINT mate_position FROM mnv'))
  comparisons <- list()
  controls <- NULL
  for (threads in c(1L, 4L)) {
    dbExecute(con, paste('SET threads=', threads))
    independent <- paste0('independent_', threads)
    dbExecute(con, paste('CREATE TABLE', independent, 'AS SELECT * FROM duckvep_annotate(',
      "'independent_input','cis',hgvs:=true,upstream_distance:=0,downstream_distance:=0)"))
    actual <- dbGetQuery(con, paste('SELECT event_index,transcript_index,protein_hgvs,',
      'protein_hgvs_status,status_code FROM', independent,
      'ORDER BY event_index'))
    stopifnot(identical(as.numeric(actual$event_index), as.numeric(inputs$mnv$event_index)),
      identical(as.numeric(actual$transcript_index), as.numeric(models$transcript_index)),
      all(actual$protein_hgvs_status == 'supported'), all(actual$status_code == 0L))
    comparison <- cis_codon_equal(data.frame(id = models$id, hgvsp = actual$protein_hgvs), expected)
    comparisons[[independent]] <- data.frame(route = 'independent_mnv', threads = threads,
      id = models$id, expected = expected$hgvsp, actual = actual$protein_hgvs,
      hgvsp_equal = comparison, cds_equal = NA, protein_equal = NA)
    dbExecute(con, paste('COPY', independent, 'TO', q(file.path(out,
      paste0(independent, '.parquet'))), '(FORMAT PARQUET)'))
    for (representation in c('mnv', 'snv')) for (mode in c('alt_events', 'source_records')) {
      table <- paste(representation, mode, threads, sep = '_')
      calls <- paste('SELECT event_index,seq_region,position,reference,transcript_index,',
        '0 sample_index,', if (mode == 'alt_events')
          'alternate,1 alt_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set' else
          "[alternate] alternates,'1|1' gt", 'FROM', representation)
      dbExecute(con, paste0('CREATE TABLE ', table, ' AS SELECT * FROM duckvep_haplotypes(',
        q(calls), ",'cis',hgvs:=true,input_mode:=", q(mode), ',phase_policy:=',
        q(if (mode == 'alt_events') 'strict' else 'vep116_compat'), ')'))
      actual <- dbGetQuery(con, paste('SELECT * FROM', table, 'ORDER BY transcript_index'))
      saveRDS(actual, file.path(out, paste0(table, '.rds')))
      cis_codon_check_native(actual, inputs[[representation]], models,
        ploidy = if (mode == 'alt_events') 1L else 2L)
      if (mode == 'source_records') stopifnot(all(vapply(actual$contributors, function(x)
        !is.null(x$alt_index) && all(x$alt_index == 1L), TRUE)))
      stopifnot(all(actual$hgvsp_status == 'ok'))
      if (is.null(controls)) controls <- cis_codon_controls(actual, inputs[[representation]],
        models, expected)
      # Compare suffixes; remove only one outer prediction wrapper. Compound
      # brackets, equality positions and all other HGVS text remain intact.
      observed <- sub('^p\\.\\((.*)\\)$', 'p.\\1', actual$hgvsp)
      comparison <- cis_codon_equal(data.frame(id = models$id, hgvsp = observed), expected)
      comparisons[[table]] <- data.frame(route = paste(representation, mode, sep = '_'),
        threads = threads, id = models$id, expected = expected$hgvsp, actual = observed,
        hgvsp_equal = comparison, cds_equal = cis_codon_sequence_equal(actual$cds, models$expected_cds),
        protein_equal = cis_codon_sequence_equal(actual$protein, expected$protein))
    }
  }
  routes <- c('independent', 'mnv_alt_events', 'mnv_source_records',
    'snv_alt_events', 'snv_source_records')
  thread_differences <- setNames(vapply(routes, function(route) {
    first <- paste0('SELECT * FROM ', route, '_1')
    last <- paste0('SELECT * FROM ', route, '_4')
    as.numeric(dbGetQuery(con, paste('SELECT count(*) n FROM ((', first, 'EXCEPT ALL', last,
      ') UNION ALL (', last, 'EXCEPT ALL', first, '))'))$n)
  }, 0), routes)
  comparisons <- do.call(rbind, comparisons)
  write.csv(comparisons, file.path(out, 'comparisons.csv'), row.names = FALSE, na = '')
  summary <- aggregate(cbind(comparisons = rep(1L, nrow(comparisons)),
    hgvsp_mismatches = as.integer(!comparisons$hgvsp_equal)),
    comparisons[c('route', 'threads')], sum)
  write.csv(summary, file.path(out, 'summary.csv'), row.names = FALSE)
  stopifnot(identical(hashes, vapply(sources, sha, '')), sha(extension) == extension_hash,
    sha(file.path(out, 'duckhts.duckdb_extension')) == extension_hash,
    revision == duckvep_evidence_revision('.'))
  for (path in sources) {
    destination <- file.path(out, 'source', path)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(path, destination), sha(destination) == hashes[[path]])
  }
  files <- list.files(out, recursive = TRUE, full.names = TRUE)
  jsonlite::write_json(list(scope = 'internal-codon exact MNV compatibility and cis-SNV/MNV contrast',
    oracle_revisions = as.list(pins), source_revision = revision, build_binding = 'diagnostic_unbound',
    extension_sha256 = extension_hash, models = nrow(models), mnv_records = nrow(inputs$mnv),
    snv_records = nrow(inputs$snv), comparisons = nrow(comparisons),
    hgvsp_mismatches = sum(!comparisons$hgvsp_equal),
    cds_mismatches = sum(!comparisons$cds_equal, na.rm = TRUE),
    protein_mismatches = sum(!comparisons$protein_equal, na.rm = TRUE),
    failures_waived = 0L, compound_hgvs_certified = FALSE,
    corruption_controls = as.list(controls),
    thread_differences = as.list(thread_differences),
    source_sha256 = as.list(hashes),
    sha256 = as.list(setNames(vapply(files, sha, ''), substring(files, nchar(out) + 2L)))),
    file.path(out, 'receipt.json'), pretty = TRUE, auto_unbox = TRUE)
  print(summary)
  stopifnot(all(comparisons$hgvsp_equal), all(comparisons$cds_equal, na.rm = TRUE),
    all(comparisons$protein_equal, na.rm = TRUE), all(thread_differences == 0))
}

if (sys.nframe() == 0L) cis_codon_main()
