#!/usr/bin/env Rscript
# Original container JSON keeps CDS and protein groups separate: one CDS can
# belong to both a curated reference protein and a raw mutation protein.
source('scripts/duckvep_evidence.R')
source('test/duckvep/conformance/haplotype_observations.R')

reference_axes <- function(rows) {
  if (!is.data.frame(rows) || !identical(names(rows), c('axis', 'sample', 'sequence', 'count')) ||
      !is.character(rows$axis) || !is.character(rows$sample) || !is.character(rows$sequence) ||
      !is.numeric(rows$count) ||
      anyNA(rows) || any(!rows$axis %in% c('cds', 'protein')) ||
      any(!is.finite(rows$count) | rows$count <= 0 | rows$count != floor(rows$count))) return(NULL)
  result <- aggregate(count ~ axis + sample + sequence, rows, sum)
  result <- result[order(result$axis, result$sample, result$sequence), ]
  rownames(result) <- NULL
  result
}

reference_carriers_valid <- function(carriers) {
  is.data.frame(carriers) && nrow(carriers) == 4L &&
    all(c('sample_index', 'haplotype_lane', 'ploidy', 'phase_set') %in% names(carriers)) &&
    all(vapply(carriers[c('sample_index', 'haplotype_lane', 'ploidy', 'phase_set')], is.numeric, TRUE)) &&
    !anyNA(carriers[c('sample_index', 'haplotype_lane', 'ploidy')]) &&
    !anyDuplicated(carriers[c('sample_index', 'haplotype_lane')]) &&
    all(sort(carriers$sample_index) == c(0, 0, 2, 2)) &&
    all(carriers$haplotype_lane %in% 1:2) && all(carriers$ploidy == 2L) && all(is.na(carriers$phase_set))
}

reference_sources_valid <- function(actual, records) {
  for (i in seq_len(nrow(actual))) {
    c <- actual$contributors[[i]]
    c <- c[order(c$event_index), ]
    carriers <- actual$carriers[[i]]
    for (j in seq_len(nrow(carriers))) {
      gt <- records[[paste0('s', carriers$sample_index[j])]]
      missing <- gt %in% c('.', './.', '.|.')
      alt <- gt == '1|1' | (gt == '0|1' & carriers$haplotype_lane[j] == 2L)
      r <- records[missing | alt, ]
      if (nrow(c) != nrow(r) || anyNA(c[c('event_index', 'seq_region', 'position', 'reference',
          'alternate', 'alt_index', 'evidence_flags', 'projection_status')]) ||
          any(c$event_index != r$event_index | c$seq_region != r$seq_region |
            c$position != r$position | c$reference != r$reference |
            c$alternate != ifelse(alt[missing | alt], r$alt, r$reference) |
            c$alt_index != as.integer(alt[missing | alt]) |
            c$evidence_flags != ifelse(missing[missing | alt], 10L, 1L) |
            c$projection_status != ifelse(r$intronic, 'outside_cds', 'ok'))) return(FALSE)
    }
  }
  TRUE
}

# Missing-only carriers remain in the full native checks. The upstream mutator
# observes only samples with retained calls; its lane flags use that exact domain.
reference_mutator_rows <- function(actual, retained, samples) {
  actual$carriers <- lapply(actual$carriers, function(x)
    x[samples[x$sample_index + 1L] %in% retained, , drop = FALSE])
  actual$carrier_count <- vapply(actual$carriers, nrow, 1L)
  actual[actual$carrier_count > 0L, , drop = FALSE]
}

main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option('--seed', type = 'integer', default = 173L),
    optparse::make_option('--rare-per-stratum', dest = 'quota', type = 'integer', default = 1L),
    optparse::make_option('--extension-receipt', dest = 'extension_receipt', default = NULL),
    optparse::make_option('--vep-prefix', dest = 'prefix', default = '/root/miniconda3/envs/vep'))))
  stopifnot(!is.na(opt$seed), !is.na(opt$quota), opt$quota >= 1L, opt$quota <= 32768L %/% 360L)
  root <- normalizePath('.')
  revision <- duckvep_evidence_revision(root)
  extension <- normalizePath('build/release/duckhts.duckdb_extension')
  binding <- 'diagnostic_unbound'
  if (!is.null(opt$extension_receipt)) {
    duckvep_evidence_assert_checkout(root, revision)
    binding <- duckvep_evidence_read_extension_receipt(opt$extension_receipt, root, extension, revision)$binding
  }
  out <- tempfile(paste0('haplotype_reference_seed', opt$seed, '_'), tmpdir = 'test/duckvep/conformance/results')
  dir.create(out)
  message('Artifacts: ', out)
  pins <- c(vep = '57ea5c52340acc1f156267f810ad162e26597082', variation = '2fb834b987ede3824e200197a838ce11e91aeb4b')
  mirrors <- normalizePath(c('.sync/ensembl-vep', '.sync/ensembl-variation'))
  for (i in seq_along(pins)) stopifnot(
    identical(duckvep_evidence_command('git', c('-C', mirrors[i], 'rev-parse', 'HEAD'), 'revision'), unname(pins[i])),
    !length(duckvep_evidence_command('git', c('-C', mirrors[i], 'status', '--porcelain'), 'clean oracle')))
  prefix <- normalizePath(opt$prefix)
  environment <- duckvep_evidence_command('micromamba', c('list', '-p', prefix, '--explicit'), 'environment')
  stopifnot(identical(duckvep_evidence_explicit_packages(environment),
    duckvep_evidence_explicit_packages(readLines('test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt'))))
  writeLines(environment, file.path(out, 'environment.txt'))
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(root, 'r/duckhtsbench/inst/benchmark_registry.tsv'))
  paths <- duckhtsbench::duckhts_bench_stage_repository_fixtures(root, 'duckvep-haplotypes')
  reference <- readLines(paths[['haplotype_benchmark_reference']])[2L]
  strata <- expand.grid(start = c('ATG', 'CTG', 'TTG'), internal_stop = c(FALSE, TRUE),
    terminal_stop = c(FALSE, TRUE), strand = c(1L, -1L), missing_gt = c('.', './.', '.|.'),
    route = c('missing_only', 'retained_ref_lane', 'later_retained',
      'intronic_short', 'intronic_long'), stringsAsFactors = FALSE)
  cases <- strata[rep(seq_len(nrow(strata)), each = opt$quota), ]
  cases$draw <- rep(seq_len(opt$quota), nrow(strata))
  stopifnot(nrow(strata) == 360L)
  set.seed(opt$seed)
  rc <- function(x) paste(rev(strsplit(chartr('ACGT', 'TGCA', x), '', fixed = TRUE)[[1L]]), collapse = '')
  models <- exons <- records <- fasta <- gff <- vector('list', nrow(cases))
  for (i in seq_len(nrow(cases))) {
    x <- cases[i, ]
    cds <- reference
    substr(cds, 1L, 3L) <- x$start
    if (!x$terminal_stop) substr(cds, 178L, 180L) <- 'GCT'
    if (x$internal_stop) {
      stop_start <- 3L * sample(1:54, 1L) + 1L
      substr(cds, stop_start, stop_start + 2L) <- sample(c('TAA', 'TAG', 'TGA'), 1L)
    }
    pos <- sort(sample(14:187, 2L))
    intronic <- x$route %in% c('intronic_short', 'intronic_long')
    intron_length <- if (x$route == 'intronic_short') sample(1:12, 1L) else
      if (x$route == 'intronic_long') sample(13:120, 1L) else 0L
    genomic_cds <- if (x$strand > 0L) cds else rc(cds)
    genome <- paste0(strrep('A', 10L), substring(genomic_cds, 1L, 90L),
      strrep('C', intron_length), substring(genomic_cds, 91L, 180L), strrep('A', 10L))
    if (intronic) pos <- c(sample(14:97, 1L), 100L + sample.int(intron_length, 1L))
    refs <- substring(genome, pos, pos)
    r <- data.frame(position = pos, reference = refs, alt = chartr('ACGT', 'CGTA', refs),
      source_id = c('first', 'second'), intronic = c(FALSE, intronic),
      s0 = c(x$missing_gt, if (x$route == 'later_retained' || intronic) '0|1' else x$missing_gt),
      s1 = '0|0', s2 = c(if (x$route == 'retained_ref_lane') '0|1' else '1|1', '0|0'))
    chr <- sprintf('chr%05d', i)
    tx <- paste0('T', i)
    r$chrom <- chr; r$seq_region <- i - 1L; r$event_index <- 2L * i - 2L + seq_len(2L)
    records[[i]] <- r
    models[[i]] <- data.frame(transcript_index = i - 1L, seq_region = i - 1L, transcript = tx,
      strand = x$strand, cds = cds, transcript_end = 190L + intron_length,
      intron_length = intron_length)
    fasta[[i]] <- c(paste0('>', chr), genome)
    e <- if (intronic) data.frame(exon_start = c(11L, 101L + intron_length),
      exon_end = c(100L, 190L + intron_length),
      exon_cdna_start = if (x$strand > 0L) c(1L, 91L) else c(91L, 1L),
      exon_cdna_end = if (x$strand > 0L) c(90L, 180L) else c(180L, 90L)) else
      data.frame(exon_start = 11L, exon_end = 190L, exon_cdna_start = 1L, exon_cdna_end = 180L)
    e$transcript_index <- i - 1L
    exons[[i]] <- e
    gff[[i]] <- c(paste(chr, 'reference_route', c('gene', 'mRNA'), 11, 190L + intron_length, '.',
      if (x$strand > 0L) '+' else '-', '.',
      c(paste0('ID=gene:', tx, ';biotype=protein_coding'),
        paste0('ID=transcript:', tx, ';Parent=gene:', tx, ';biotype=protein_coding')), sep = '\t'),
      unlist(lapply(seq_len(nrow(e)), function(j)
        paste(chr, 'reference_route', c('exon', 'CDS'), e$exon_start[j], e$exon_end[j], '.',
          if (x$strand > 0L) '+' else '-', c('.', '0'),
          c(paste0('ID=exon:', tx, ':', j, ';Parent=transcript:', tx),
            paste0('Parent=transcript:', tx)), sep = '\t'))))
  }
  models <- do.call(rbind, models); records <- do.call(rbind, records); exons <- do.call(rbind, exons)
  saveRDS(list(cases = cases, models = models, records = records, exons = exons), file.path(out, 'inputs.rds'))
  key <- function(x) do.call(paste, c(x[names(strata)], sep = ':'))
  coverage <- strata; coverage$required <- opt$quota
  coverage$observed <- tabulate(match(key(cases), key(strata)), nrow(strata))
  stopifnot(all(coverage$observed == opt$quota))
  write.csv(coverage, file.path(out, 'coverage.csv'), row.names = FALSE)
  writeLines(unlist(fasta), file.path(out, 'reference.fa'))
  writeLines(c('##gff-version 3', unlist(gff)), file.path(out, 'model.gff3'))
  samples <- c('missing', 'reference', 'called')
  writeLines(c('##fileformat=VCFv4.4', paste0('##contig=<ID=', unique(records$chrom), ',length=', models$transcript_end + 10L, '>'),
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste(c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', samples), collapse = '\t'),
    with(records, paste(chrom, position, source_id, reference, alt, '.', 'PASS', '.', 'GT', s0, s1, s2, sep = '\t'))),
    file.path(out, 'calls.vcf'))
  run <- function(cmd, args, name) stopifnot(system2(cmd, shQuote(args),
    stdout = file.path(out, paste0(name, '.stdout')), stderr = file.path(out, paste0(name, '.stderr'))) == 0L)
  run('samtools', c('faidx', file.path(out, 'reference.fa')), 'faidx')
  run('bcftools', c('norm', '-c', 'e', '-f', file.path(out, 'reference.fa'), '-o', file.path(out, 'checked.vcf'),
    file.path(out, 'calls.vcf')), 'refcheck')
  run('bgzip', file.path(out, 'model.gff3'), 'bgzip')
  run('tabix', c('-p', 'gff', file.path(out, 'model.gff3.gz')), 'tabix')
  libs <- paste(c(file.path(mirrors, 'modules'), file.path(prefix, 'share/ensembl-vep-116.0-0')), collapse = ':')
  run('micromamba', c('run', '--clean-env', '--env', paste0('PERL5LIB=', libs), '-p', prefix,
    'perl', 'test/duckvep/conformance/haplotype_oracle.pl', '--container-json', file.path(out, 'calls.vcf'),
    file.path(out, 'reference.fa'), file.path(out, 'model.gff3.gz'), file.path(out, 'phase.jsonl')), 'oracle')
  run('micromamba', c('run', '--clean-env', '--env', paste0('PERL5LIB=', libs), '-p', prefix,
    'perl', 'test/duckvep/conformance/haplotype_oracle.pl', '--container-json',
    file.path(out, c('calls.vcf', 'reference.fa', 'model.gff3.gz'))), 'unobserved')
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = 'true')))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbExecute(con, paste('LOAD', q(extension))); DBI::dbExecute(con, 'SET threads=4')
  DBI::dbWriteTable(con, 'models', models); DBI::dbWriteTable(con, 'records', records)
  DBI::dbWriteTable(con, 'exons', exons)
  queries <- c('SELECT seq_region::UINTEGER seq_region FROM models ORDER BY seq_region',
    'SELECT transcript_index::UINTEGER transcript_index,seq_region::UINTEGER seq_region,11::UBIGINT transcript_start,
     transcript_end::UBIGINT transcript_end,strand::TINYINT strand,transcript_index::UINTEGER gene_index,3::UBIGINT transcript_flags,
     11::UBIGINT cds_start,transcript_end::UBIGINT cds_end,cds::BLOB cds_sequence,1::UTINYINT codon_table FROM models ORDER BY transcript_index',
    'SELECT transcript_index::UINTEGER transcript_index,exon_start::UBIGINT exon_start,exon_end::UBIGINT exon_end,
     exon_cdna_start::UBIGINT exon_cdna_start,exon_cdna_end::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase
     FROM exons ORDER BY transcript_index,exon_cdna_start')
  stopifnot(DBI::dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('reference_route',", paste(q(queries), collapse = ','), ')'))$loaded)
  calls <- 'SELECT event_index,r.seq_region,position,reference,[alt] alternates,m.transcript_index,s.i sample_index,
    CASE s.i WHEN 0 THEN s0 WHEN 1 THEN s1 ELSE s2 END gt FROM records r JOIN models m USING(seq_region),range(3) s(i)'
  actual <- DBI::dbGetQuery(con, paste0('SELECT * FROM duckvep_haplotypes(', q(calls),
    ",'reference_route',input_mode:='source_records',phase_policy:='vep116_compat')"))
  saveRDS(actual, file.path(out, 'actual.rds'))
  oracle_lines <- readLines(file.path(out, 'oracle.stdout'))
  oracle <- lapply(oracle_lines, jsonlite::fromJSON, simplifyVector = FALSE)
  phase <- lapply(readLines(file.path(out, 'phase.jsonl')), jsonlite::fromJSON, simplifyVector = FALSE)
  unobserved_lines <- readLines(file.path(out, 'unobserved.stdout'))
  unobserved <- lapply(unobserved_lines, jsonlite::fromJSON, simplifyVector = FALSE)
  names(unobserved) <- vapply(unobserved, `[[`, '', 'transcript_id')
  names(oracle) <- vapply(oracle, `[[`, '', 'transcript_id'); names(phase) <- vapply(phase, `[[`, '', 'transcript')
  names(oracle_lines) <- names(oracle); names(unobserved_lines) <- names(unobserved)
  stopifnot(!anyDuplicated(names(oracle)), !anyDuplicated(names(phase)),
    !anyDuplicated(names(unobserved)), setequal(names(unobserved), models$transcript),
    setequal(names(oracle), models$transcript), setequal(names(phase), models$transcript))
  by_tx <- native_haplotype_rows(actual, models$transcript_index)
  comparisons <- lapply(seq_len(nrow(models)), function(i) {
    model <- models[i, ]; a <- actual[by_tx[[as.character(model$transcript_index)]], ]
    r <- records[records$seq_region == model$seq_region, ]
    upstream <- oracle[[model$transcript]]; p <- phase[[model$transcript]]
    expected <- do.call(rbind, lapply(c('cds', 'protein'), function(axis)
      do.call(rbind, lapply(upstream[[paste0(axis, '_haplotypes')]], function(h)
        data.frame(axis = axis, sample = names(h$samples), sequence = h$seq,
          count = as.numeric(unlist(h$samples)))))))
    observed <- do.call(rbind, lapply(seq_len(nrow(a)), function(j)
      do.call(rbind, lapply(c('cds', 'protein'), function(axis)
        data.frame(axis = axis, sample = samples[a$carriers[[j]]$sample_index + 1L],
          sequence = a[[axis]][j], count = 1)))))
    carriers <- do.call(rbind, a$carriers)
    retained <- unique(vapply(p$calls, `[[`, '', 'sample'))
    lanes <- native_replay_lanes(a, r, model$strand, samples)
    lanes <- Filter(function(x) x$sample %in% retained, lanes)
    e <- reference_axes(expected[expected$sample != 'reference', ])
    o <- reference_axes(observed)
    list(expected = expected, observed = observed, expected_lanes = p$replay_lanes, observed_lanes = lanes,
      observer_equal = identical(oracle_lines[[model$transcript]], unobserved_lines[[model$transcript]]),
      axes_equal = !is.null(e) && !is.null(o) && identical(e, o),
      lanes_equal = replay_lanes_equal(p$replay_lanes, lanes),
      lane_flags_equal = native_lane_flags_equal(p$replay_lanes,
        reference_mutator_rows(a, retained, samples), samples),
      group_metadata_equal = haplotype_group_metadata_equal(p, upstream, container_json = TRUE),
      model_equal = identical(model$cds, p$reference_cds),
      sources_equal = reference_sources_valid(a, r),
      counts_equal = all(vapply(split(expected$count, expected$axis), sum, 0) == 6) &&
        sum(a$carrier_count) == 4L && nrow(carriers) == 4L &&
        all(a$carrier_count == vapply(a$carriers, nrow, 1L)),
      carriers_equal = reference_carriers_valid(carriers),
      implicit_reference_equal = sum(expected$count[expected$sample == 'reference']) == 4L &&
        all(expected$sequence[expected$sample == 'reference' & expected$axis == 'cds'] == model$cds) &&
        all(expected$sequence[expected$sample == 'reference' & expected$axis == 'protein'] == p$reference_protein))
  })
  saveRDS(comparisons, file.path(out, 'comparisons.rds'))
  summary <- cbind(models[setdiff(names(models), 'cds')], cases[setdiff(names(cases), 'strand')])
  for (name in c('axes_equal', 'lanes_equal', 'model_equal', 'sources_equal', 'counts_equal', 'carriers_equal', 'implicit_reference_equal', 'observer_equal'))
    summary[[name]] <- vapply(comparisons, `[[`, TRUE, name)
  summary$passed <- with(summary, axes_equal & lanes_equal & model_equal & sources_equal & counts_equal & carriers_equal & implicit_reference_equal & observer_equal)
  for (name in c('lane_flags_equal', 'group_metadata_equal'))
    summary[[name]] <- vapply(comparisons, `[[`, TRUE, name)
  summary$metadata_passed <- summary$lane_flags_equal & summary$group_metadata_equal
  write.csv(summary, file.path(out, 'summary.csv'), row.names = FALSE)
  witness <- comparisons[[1L]]$expected
  equal <- function(x) identical(reference_axes(x), reference_axes(witness))
  controls <- c(missing_group = !equal(witness[-1L, ]), duplicate_group = !equal(rbind(witness, witness[1L, ])))
  for (field in c('axis', 'sample', 'sequence', 'count')) {
    bad <- witness
    bad[[field]][1L] <- if (field == 'count') bad$count[1L] + 1 else paste0(bad[[field]][1L], 'X')
    controls[field] <- !equal(bad)
  }
  lane_witness <- comparisons[[which(cases$route == 'retained_ref_lane')[1L]]]$expected_lanes
  controls['missing_lane'] <- !replay_lanes_equal(lane_witness, lane_witness[-1L])
  bad <- lane_witness; bad[[1L]]$protein <- paste0(bad[[1L]]$protein, 'X')
  controls['lane_protein'] <- !replay_lanes_equal(lane_witness, bad)
  carrier_witness <- do.call(rbind, actual$carriers[actual$transcript_index == 0L])
  stopifnot(reference_carriers_valid(carrier_witness))
  integer_keys <- carrier_witness; integer_keys$sample_index <- as.integer(integer_keys$sample_index)
  stopifnot(reference_carriers_valid(integer_keys))
  controls['missing_carrier'] <- !reference_carriers_valid(carrier_witness[-1L, ])
  controls['duplicate_carrier'] <- !reference_carriers_valid(rbind(carrier_witness, carrier_witness[1L, ]))
  for (field in c('sample_index', 'phase_set', 'haplotype_lane', 'ploidy')) {
    bad <- carrier_witness; bad[[field]][1L] <- 3
    controls[paste0('carrier_', field)] <- !reference_carriers_valid(bad)
  }
  bad <- carrier_witness; bad$sample_index[1L] <- 0.5
  controls['fractional_sample'] <- !reference_carriers_valid(bad)
  bad$sample_index[1L] <- NA_real_
  controls['missing_sample'] <- !reference_carriers_valid(bad)
  container_witness <- oracle_lines[[models$transcript[1L]]]
  container_equal <- function(x) identical(x, container_witness)
  controls['container_transcript'] <- !container_equal(sub('"transcript_id":"T1"',
    '"transcript_id":"wrong"', container_witness, fixed = TRUE))
  controls['container_sequence'] <- !container_equal(sub(models$cds[1L],
    paste0('X', substring(models$cds[1L], 2L)), container_witness, fixed = TRUE))
  controls['container_total'] <- !container_equal(sub('"total_haplotype_count":6',
    '"total_haplotype_count":7', container_witness, fixed = TRUE))
  bad_cds <- models$cds[1L]; substr(bad_cds, 1L, 1L) <- 'C'
  controls['model_sequence'] <- !identical(bad_cds, phase[[models$transcript[1L]]]$reference_cds)
  source_witness <- actual[actual$transcript_index == 0L, ]
  source_records <- records[records$seq_region == 0L, ]
  stopifnot(reference_sources_valid(source_witness, source_records))
  bad <- source_witness; bad$contributors[[1L]] <- bad$contributors[[1L]][-1L, ]
  controls['missing_source'] <- !reference_sources_valid(bad, source_records)
  for (field in c('event_index', 'seq_region', 'position', 'reference', 'alternate',
      'alt_index', 'evidence_flags', 'projection_status')) {
    bad <- source_witness
    bad$contributors[[1L]][[field]][1L] <- if (field %in% c('reference', 'alternate', 'projection_status'))
      'X' else bad$contributors[[1L]][[field]][1L] + 1L
    controls[paste0('source_', field)] <- !reference_sources_valid(bad, source_records)
  }
  write.csv(data.frame(control = names(controls), rejected = controls), file.path(out, 'controls.csv'), row.names = FALSE)
  phase_witness <- phase[[models$transcript[1L]]]
  upstream_witness <- oracle[[models$transcript[1L]]]
  metadata_controls <- haplotype_metadata_controls(phase_witness, upstream_witness, container_json = TRUE)
  flag_witness <- reference_mutator_rows(source_witness,
    unique(vapply(phase_witness$calls, `[[`, '', 'sample')), samples)
  flags_equal <- function(x) native_lane_flags_equal(phase_witness$replay_lanes, x, samples)
  stopifnot(flags_equal(flag_witness))
  for (bit in c(1L, 2L, 4L)) {
    bad <- flag_witness
    bad$sequence_flags[1L] <- bitwXor(bad$sequence_flags[1L], bit)
    metadata_controls[paste0('metadata_native_bit_', bit)] <- !flags_equal(bad)
  }
  bad <- flag_witness
  bad$carriers[[1L]] <- bad$carriers[[1L]][-1L, , drop = FALSE]
  metadata_controls['metadata_missing_mutator_carrier'] <- !flags_equal(bad)
  bad <- flag_witness
  extra <- bad$carriers[[1L]][1L, , drop = FALSE]
  extra$haplotype_lane <- 3L
  bad$carriers[[1L]] <- rbind(bad$carriers[[1L]], extra)
  metadata_controls['metadata_extra_mutator_carrier'] <- !flags_equal(bad)
  reference_group <- which(vapply(phase_witness$cds_group_metadata,
    function(x) identical(x$cds, phase_witness$reference_cds), TRUE))
  stopifnot(length(reference_group) == 1L,
    !any(vapply(phase_witness$replay_lanes,
      function(x) identical(x$cds, phase_witness$reference_cds), TRUE)))
  for (field in c('indel', 'frameshift', 'length_diff')) {
    bad <- phase_witness
    bad$cds_group_metadata[[reference_group]]$flags[[field]] <- 1L
    metadata_controls[paste0('metadata_reference_', field)] <-
      !haplotype_group_metadata_equal(bad, upstream_witness, container_json = TRUE)
  }
  write.csv(data.frame(control = names(metadata_controls), rejected = metadata_controls),
    file.path(out, 'metadata_controls.csv'), row.names = FALSE)
  output_controls <- haplotype_output_controls(actual)
  write.csv(data.frame(control = names(output_controls), rejected = output_controls),
    file.path(out, 'output_controls.csv'), row.names = FALSE)
  identities <- c('test/duckvep/conformance/haplotype_reference_differential.R',
    'test/duckvep/conformance/haplotype_oracle.pl', 'test/duckvep/conformance/haplotype_observations.R',
    'scripts/duckvep_evidence.R', 'r/duckhtsbench/inst/benchmark_registry.tsv',
    paths[['haplotype_benchmark_reference']], extension, list.files(out, full.names = TRUE))
  jsonlite::write_json(list(source_revision = revision, extension_build_binding = binding,
    tracked_changes = duckvep_evidence_tracked_changes(root), oracle_revisions = as.list(pins),
    seed = opt$seed, rare_per_stratum = opt$quota, strata = nrow(strata), profiles = nrow(summary),
    source_records = nrow(records), source_calls = nrow(records) * 3L, oracle_carriers = nrow(summary) * 6L,
    native_carriers = sum(actual$carrier_count), implicit_reference_carriers = nrow(summary) * 2L,
    oracle_mutation_lanes = sum(vapply(comparisons, function(x) length(x$expected_lanes), 1L)),
    failures = sum(!summary$passed), axes_failures = sum(!summary$axes_equal),
    observer_failures = sum(!summary$observer_equal),
    source_failures = sum(!summary$sources_equal),
    lane_failures = sum(!summary$lanes_equal), controls_rejected = sum(controls), threads = 4L,
    lane_flag_failures = sum(!summary$lane_flags_equal),
    group_metadata_failures = sum(!summary$group_metadata_equal),
    cds_groups_compared = sum(vapply(phase, function(x) length(x$cds_group_metadata), 1L)),
    metadata_controls_rejected = sum(metadata_controls),
    output_controls_rejected = sum(output_controls),
    scope = 'one_or_two_exon_standard_code_source_retention_and_reference_protein_routes_with_implicit_REF_samples',
    sha256 = as.list(vapply(identities, duckvep_evidence_sha256, ''))), file.path(out, 'receipt.json'), pretty = TRUE, auto_unbox = TRUE)
  if (!is.null(opt$extension_receipt)) duckvep_evidence_assert_checkout(root, revision)
  print(aggregate(cbind(profiles = rep(1L, nrow(summary)), failures = !summary$passed,
    axes_failures = !summary$axes_equal, lane_failures = !summary$lanes_equal) ~ route, summary, sum), row.names = FALSE)
  stopifnot(all(controls), all(metadata_controls), all(output_controls))
  if (any(!summary$passed)) stop('Reference-route disagreements retained: ', out, call. = FALSE)
  if (any(!summary$metadata_passed)) stop('Reference-route lane/group metadata disagreements retained: ', out, call. = FALSE)
}
main()
