#!/usr/bin/env Rscript
# Seeded shared-source replay over full and two-exon transcript models.
source('test/duckvep/conformance/haplotype_observations.R')


# This generated grammar uses pipe-separated diploid calls, including .|1.
# The pinned parser compacts called slots before consuming two file lanes.
file_slots <- function(gt) {
  slots <- strsplit(gt, '|', fixed = TRUE)[[1L]]
  slots <- as.integer(slots[slots != '.'])
  slots[1:2]
}

input_provenance_ok <- function(actual, records) {
  carriers <- character()
  for (j in seq_len(nrow(actual))) {
    calls <- actual$carriers[[j]]
    sources <- actual$contributors[[j]]
    matched <- match(sources$event_index, records$event_index)
    if (anyNA(matched) || anyDuplicated(sources$event_index) ||
        any(sources$seq_region != records$seq_region[matched]) ||
        any(sources$position != records$position[matched]) ||
        any(sources$reference != records$reference[matched])) return(FALSE)
    for (k in seq_len(nrow(calls))) {
      sample <- calls$sample_index[k]
      lane <- calls$haplotype_lane[k]
      if (is.na(sample) || !sample %in% 0:2 || !lane %in% 1:2 ||
          calls$ploidy[k] != 2L || !is.na(calls$phase_set[k])) return(FALSE)
      carriers <- c(carriers, paste(sample, lane))
      gt <- records[[paste0('s', sample)]]
      retained <- gt != '0|0'
      slots <- vapply(gt, function(x) file_slots(x)[lane], 1L, USE.NAMES = FALSE)
      required <- which(retained & (is.na(slots) | slots != 0L))
      observed <- matched[is.na(sources$alt_index) | sources$alt_index != 0L]
      if (!identical(sort(required), sort(observed))) return(FALSE)
      for (n in seq_len(nrow(sources))) {
        index <- matched[n]
        if (!retained[index] ||
            !identical(as.numeric(sources$alt_index[n]), as.numeric(slots[index]))) return(FALSE)
        alleles <- c(records$reference[index], strsplit(records$alt[index], ',', fixed = TRUE)[[1L]])
        alt <- if (is.na(slots[index])) '' else alleles[slots[index] + 1L]
        if (!identical(sources$alternate[n], alt)) return(FALSE)
      }
    }
  }
  identical(sort(carriers), sort(as.vector(outer(0:2, 1:2, paste))))
}

mappings_ok <- function(observation, model, spans, records) {
  if (is.null(observation) || !identical(observation$transcript, model$transcript)) return(FALSE)
  buffer <- observation$source_buffer
  if (length(buffer) != nrow(records)) return(FALSE)
  ids <- vapply(buffer, function(x) paste(x$ids, collapse = ','), '')
  matched <- match(ids, records$source_id)
  if (anyNA(matched) || anyDuplicated(matched)) return(FALSE)
  for (i in seq_along(buffer)) {
    record <- records[matched[i], ]
    b <- buffer[[i]]
    alleles <- c(record$reference, strsplit(record$alt, ',', fixed = TRUE)[[1L]])
    if (!identical(b$chrom, record$chrom) || b$start != record$position ||
        b$end != record$position + nchar(record$reference) - 1L ||
        !identical(b$alleles, paste(alleles, collapse = ','))) return(FALSE)
  }
  calls <- observation$calls
  expected_keys <- character()
  for (i in seq_len(nrow(records))) {
    end <- records$position[i] + nchar(records$reference[i]) - 1L
    copies <- sum(spans$start <= end & spans$end >= records$position[i])
    for (sample in 0:2) if (records[[paste0('s', sample)]][i] != '0|0')
      expected_keys <- c(expected_keys, rep(paste(records$source_id[i], paste0('s', sample)), copies))
  }
  keys <- vapply(calls, function(x) paste(x$source_id, x$sample), '')
  if (!identical(sort(keys), sort(expected_keys))) return(FALSE)
  identifiers <- vapply(calls, `[[`, '', 'source_key')
  selected <- vapply(split(calls, identifiers), function(group) tail(group, 1L)[[1L]]$source_id, '')
  for (call in calls) {
    record <- records[match(call$source_id, records$source_id), ]
    end <- record$position + nchar(record$reference) - 1L
    alleles <- c(record$reference, strsplit(record$alt, ',', fixed = TRUE)[[1L]])
    identifier <- paste(record$position, end, paste(alleles, collapse = '/'), sep = '_')
    if (!identical(call$source_key, identifier)) return(FALSE)
    gt <- file_slots(record[[call$sample]])
    expected_gt <- unname(alleles[gt[!is.na(gt)] + 1L])
    if (!identical(unlist(call$genotype, use.names = FALSE), expected_gt)) return(FALSE)
    span <- spans[spans$start <= record$position & spans$end >= end, ]
    if (nrow(span) && call$source_id == selected[[identifier]]) {
      start_cds <- if (model$strand == 1L) span$cs + record$position - span$start else span$cs + span$end - end
      end_cds <- if (model$strand == 1L) span$cs + end - span$start else span$cs + span$end - record$position
      if (nrow(span) != 1L || !identical(as.numeric(call$mapping_start), as.numeric(start_cds)) ||
          !identical(as.numeric(call$mapping_end), as.numeric(end_cds))) return(FALSE)
    } else if (!is.null(call$mapping_start) || !is.null(call$mapping_end)) return(FALSE)
  }
  TRUE
}

observation_controls <- function(witness, actual, phase, model, exons, records) {
  corrupt <- list(missing_row = witness[-1L], duplicate_row = c(witness, witness[1L]),
    cds = witness, protein = witness, contributor = witness, sample = witness)
  corrupt$cds[[1L]]$cds <- paste0(witness[[1L]]$cds, 'A')
  corrupt$protein[[1L]]$protein <- paste0(witness[[1L]]$protein, 'X')
  corrupt$contributor[[1L]]$contributors <- c(witness[[1L]]$contributors, 'extra_source')
  for (i in seq_along(witness)) {
    names <- names(corrupt$sample[[i]]$samples)
    names(corrupt$sample[[i]]$samples) <- ifelse(names == 's0', 's1', ifelse(names == 's1', 's0', names))
  }
  rejected <- vapply(corrupt, function(x)
    !identical(canonical(x, samples = TRUE), canonical(witness, samples = TRUE)), TRUE)
  stopifnot(identical(canonical(corrupt$sample, FALSE), canonical(witness, FALSE)))
  stopifnot(input_provenance_ok(actual, records), mappings_ok(phase, model, exons, records))
  wrong <- missing <- extra <- actual
  wrong$contributors[[1L]]$position[1L] <- wrong$contributors[[1L]]$position[1L] + 1L
  missing$contributors[[1L]] <- missing$contributors[[1L]][-1L, ]
  extra$carriers[[1L]] <- rbind(extra$carriers[[1L]], extra$carriers[[1L]][1L, ])
  rejected <- c(rejected, wrong_input_position = !input_provenance_ok(wrong, records),
    missing_source = !input_provenance_ok(missing, records),
    duplicate_carrier = !input_provenance_ok(extra, records))
  mapped <- which(vapply(phase$calls, function(x) !is.null(x$mapping_start), TRUE))[1L]
  wrong <- missing <- dropped <- buffer <- phase
  wrong$calls[[mapped]]$mapping_start <- wrong$calls[[mapped]]$mapping_start + 1L
  missing$calls[[mapped]]$mapping_start <- missing$calls[[mapped]]$mapping_end <- NULL
  dropped$calls <- dropped$calls[-1L]
  buffer$source_buffer <- buffer$source_buffer[-1L]
  c(rejected, wrong_mapping = !mappings_ok(wrong, model, exons, records),
    missing_mapping = !mappings_ok(missing, model, exons, records),
    missing_genotype = !mappings_ok(dropped, model, exons, records),
    missing_buffer_record = !mappings_ok(buffer, model, exons, records))
}

main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option('--seed', type = 'integer', default = 173L),
    optparse::make_option('--rare-per-stratum', dest = 'rare_per_stratum', type = 'integer', default = 1L),
    optparse::make_option('--extension-receipt', dest = 'extension_receipt', default = NULL),
    optparse::make_option('--vep-prefix', dest = 'vep_prefix', default = '/root/miniconda3/envs/vep')
  )))
  source('scripts/duckvep_evidence.R', local = TRUE)
  root <- normalizePath('.')
  seed <- opt$seed
  stopifnot(!is.na(seed), !is.na(opt$rare_per_stratum), opt$rare_per_stratum >= 0L,
    opt$rare_per_stratum <= (32768L - 768L) %/% 264L)
  revision <- duckvep_evidence_revision(root)
  extension <- normalizePath('build/release/duckhts.duckdb_extension')
  binding <- 'diagnostic_unbound'
  if (!is.null(opt$extension_receipt)) {
    duckvep_evidence_assert_checkout(root, revision)
    binding <- duckvep_evidence_read_extension_receipt(opt$extension_receipt, root, extension, revision)$binding
  }
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(root, 'r/duckhtsbench/inst/benchmark_registry.tsv'))
  paths <- duckhtsbench::duckhts_bench_stage_repository_fixtures(root, 'duckvep-haplotypes')
  cds <- readLines(paths[['haplotype_benchmark_reference']])[2L]
  rc <- function(x) paste(rev(strsplit(chartr('ACGT', 'TGCA', x), '', fixed = TRUE)[[1L]]), collapse = '')
  out <- tempfile(paste0('haplotype_models_seed', seed, '_'),
    tmpdir = 'test/duckvep/conformance/results')
  dir.create(out)
  message('Artifacts: ', out)
  pins <- c(vep = '57ea5c52340acc1f156267f810ad162e26597082',
    variation = '2fb834b987ede3824e200197a838ce11e91aeb4b')
  mirrors <- normalizePath(c('.sync/ensembl-vep', '.sync/ensembl-variation'))
  for (i in seq_along(pins)) stopifnot(
    identical(duckvep_evidence_command('git', c('-C', mirrors[i], 'rev-parse', 'HEAD'), 'revision'), unname(pins[i])),
    !length(duckvep_evidence_command('git', c('-C', mirrors[i], 'status', '--porcelain'), 'clean oracle')))
  prefix <- normalizePath(opt$vep_prefix)
  environment <- duckvep_evidence_command('micromamba', c('list', '-p', prefix, '--explicit'), 'environment')
  stopifnot(identical(duckvep_evidence_explicit_packages(environment),
    duckvep_evidence_explicit_packages(readLines('test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt'))))
  writeLines(environment, file.path(out, 'environment.txt'))
  set.seed(seed)
  cases <- expand.grid(shape = c('disjoint', 'duplicate', 'different_alt_lists', 'exon_spanning'),
    gt_pair = c('cis', 'trans', 'mixed_missing'), strand = c(1L, -1L), draw = 1:32,
    stringsAsFactors = FALSE)
  cases$neutral_count <- c(0L, 1L, 4L, 5L, 6L, 14L, 15L, 16L, 34L, 35L, 36L)[(cases$draw - 1L) %% 11L + 1L]
  cases$cohort <- 'baseline'
  strata <- expand.grid(shape = unique(cases$shape), gt_pair = unique(cases$gt_pair),
    strand = c(1L, -1L), neutral_count = c(0L, 1L, 4L, 5L, 6L, 14L, 15L, 16L, 34L, 35L, 36L),
    stringsAsFactors = FALSE)
  if (opt$rare_per_stratum) {
    rare <- strata[rep(seq_len(nrow(strata)), each = opt$rare_per_stratum), ]
    rare$draw <- rep(seq_len(opt$rare_per_stratum), nrow(strata))
    rare$cohort <- 'quota'
    cases <- rbind(cases, rare[names(cases)])
  }
  models <- exons <- records <- fasta <- gff <- vector('list', nrow(cases))
  gts <- list(cis = c('1|1', '1|0', '0|1'), trans = c('1|0', '0|1', '1|1'),
    mixed_missing = c('.|1', '0|0', '1|1'))
  for (i in seq_len(nrow(cases))) {
    chr <- sprintf('chr%04d', i)
    strand <- cases$strand[i]
    genome <- paste0(strrep('A', 10L), if (strand == 1L) cds else rc(cds), strrep('A', 10L))
    spliced <- paste0(substr(genome, 11L, 70L), substr(genome, 101L, 190L))
    if (strand == -1L) spliced <- rc(spliced)
    fasta[[i]] <- c(paste0('>', chr), genome)
    models[[i]] <- data.frame(transcript_index = 2L * (i - 1L) + 0:1, seq_region = i - 1L,
      transcript = paste0('T', i, c('full', 'spliced')), strand = strand,
      kind = c('full', 'spliced'), cds = c(cds, spliced))
    for (j in 0:1) {
      tx <- models[[i]]$transcript[j + 1L]
      span <- if (j) data.frame(start = c(11L, 101L), end = c(70L, 190L)) else
        data.frame(start = 11L, end = 190L)
      if (strand == -1L) span <- span[rev(seq_len(nrow(span))), ]
      span$ce <- cumsum(span$end - span$start + 1L)
      span$cs <- c(1L, head(span$ce, -1L) + 1L)
      exons[[i]] <- rbind(exons[[i]], data.frame(transcript_index = 2L * (i - 1L) + j, span))
      gf <- data.frame(feature = c('gene', 'mRNA'), start = 11L, end = 190L, phase = '.',
        attr = c(paste0('ID=gene:', tx, ';biotype=protein_coding'),
          paste0('ID=transcript:', tx, ';Parent=gene:', tx, ';biotype=protein_coding')))
      for (k in seq_len(nrow(span))) gf <- rbind(gf, data.frame(feature = c('exon', 'CDS'),
        start = span$start[k], end = span$end[k], phase = c('.', '0'),
        attr = c(paste0('ID=exon:', tx, k, ';Parent=transcript:', tx), paste0('Parent=transcript:', tx))))
      gff[[i]] <- with(gf, paste(chr, 'probe', feature, start, end, '.',
        if (strand == 1L) '+' else '-', phase, attr, sep = '\t')) |>
        c(gff[[i]])
    }
    pos <- c(sample(35:60, 1L), sample(105:140, 1L), 165L)
    len <- rep(1L, 3L)
    if (cases$shape[i] %in% c('duplicate', 'different_alt_lists')) pos[2L] <- pos[1L]
    if (cases$shape[i] == 'exon_spanning') {
      pos[1:2] <- c(sample(61:70, 1L), sample(101:115, 1L))
      len[1L] <- sample(116:140, 1L) - pos[1L] + 1L
    }
    ref <- substring(genome, pos, pos + len - 1L)
    alt <- chartr('ACGT', 'CGTA', ref)
    if (cases$shape[i] == 'different_alt_lists') {
      extra <- setdiff(c('A', 'C', 'G', 'T'), c(ref[1L], alt[1L]))
      alt[1:2] <- paste0(alt[1:2], ',', extra)
    }
    gt <- rbind(gts[[cases$gt_pair[i]]], rev(gts[[cases$gt_pair[i]]]), rep('1|1', 3L))
    r <- data.frame(chrom = chr, seq_region = i - 1L, source_id = c('a', 'b', 'anchor'),
      position = pos, reference = ref, alt = alt, s0 = gt[, 1L], s1 = gt[, 2L], s2 = gt[, 3L])
    if (cases$neutral_count[i]) {
      neutral <- r[rep(1L, cases$neutral_count[i]), ]
      neutral$source_id <- paste0('neutral', seq_len(nrow(neutral)))
      neutral$position <- 12L + (seq_len(nrow(neutral)) - 1L) %% 16L
      neutral$reference <- substring(genome, neutral$position, neutral$position)
      neutral$alt <- chartr('ACGT', 'CGTA', neutral$reference)
      neutral$s0 <- neutral$s1 <- neutral$s2 <- '0|0'
      r <- rbind(neutral, r)
    }
    records[[i]] <- r
  }
  models <- do.call(rbind, models); exons <- do.call(rbind, exons); records <- do.call(rbind, records)
  records <- records[order(records$seq_region, records$position, seq_len(nrow(records))), ]
  records$event_index <- seq_len(nrow(records))
  stopifnot(nrow(cases) == 768L + 264L * opt$rare_per_stratum, nrow(models) == 2L * nrow(cases),
    all(table(interaction(cases[cases$cohort == 'baseline', c('shape', 'gt_pair', 'strand')])) == 32L),
    all(records$reference != records$alt))
  key <- function(x) do.call(paste, c(x[c('shape', 'gt_pair', 'strand', 'neutral_count')], sep = ':'))
  coverage <- strata
  coverage$required <- opt$rare_per_stratum
  coverage$observed <- tabulate(match(key(cases[cases$cohort == 'quota', ]), key(strata)), nrow(strata))
  stopifnot(nrow(strata) == 264L, all(coverage$observed == coverage$required))
  write.csv(coverage, file.path(out, 'coverage.csv'), row.names = FALSE)
  saveRDS(list(cases = cases, models = models, exons = exons, records = records), file.path(out, 'inputs.rds'))
  writeLines(unlist(fasta), file.path(out, 'reference.fa'))
  gff_lines <- unlist(gff)
  parts <- strsplit(gff_lines, '\t')
  gff_order <- order(vapply(parts, `[[`, '', 1L), as.integer(vapply(parts, `[[`, '', 4L)))
  writeLines(c('##gff-version 3', gff_lines[gff_order]), file.path(out, 'model.gff3'))
  writeLines(c('##fileformat=VCFv4.4', paste0('##contig=<ID=', unique(records$chrom), ',length=200>'),
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2',
    with(records, paste(chrom, position, source_id, reference, alt, '.', 'PASS', '.', 'GT', s0, s1, s2,
      sep = '\t'))), file.path(out, 'calls.vcf'))
  run <- function(cmd, args, name) stopifnot(system2(cmd, shQuote(args),
    stdout = file.path(out, paste0(name, '.stdout')), stderr = file.path(out, paste0(name, '.stderr'))) == 0L)
  run('samtools', c('faidx', file.path(out, 'reference.fa')), 'faidx')
  run('bcftools', c('norm', '-c', 'e', '-f', file.path(out, 'reference.fa'), '-o',
    file.path(out, 'checked.vcf'), file.path(out, 'calls.vcf')), 'refcheck')
  run('bgzip', file.path(out, 'model.gff3'), 'bgzip')
  run('tabix', c('-p', 'gff', file.path(out, 'model.gff3.gz')), 'tabix')
  libs <- paste(c(file.path(mirrors, 'modules'), file.path(prefix, 'share/ensembl-vep-116.0-0')), collapse = ':')
  run('micromamba', c('run', '--clean-env', '--env', paste0('PERL5LIB=', libs), '-p', prefix,
    'perl', 'test/duckvep/conformance/haplotype_oracle.pl', file.path(out, 'calls.vcf'),
    file.path(out, 'reference.fa'), file.path(out, 'model.gff3.gz'), file.path(out, 'phase.jsonl')), 'oracle')
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = 'true')))
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbExecute(con, paste('LOAD', q(extension)))
  DBI::dbExecute(con, 'SET threads=4')
  DBI::dbWriteTable(con, 'models', models); DBI::dbWriteTable(con, 'exons', exons)
  DBI::dbWriteTable(con, 'records', records)
  queries <- c('SELECT DISTINCT seq_region::UINTEGER seq_region FROM models ORDER BY seq_region',
    'SELECT transcript_index::UINTEGER transcript_index,seq_region::UINTEGER seq_region,11::UBIGINT transcript_start,
     190::UBIGINT transcript_end,strand::TINYINT strand,transcript_index::UINTEGER gene_index,3::UBIGINT transcript_flags,
     11::UBIGINT cds_start,190::UBIGINT cds_end,cds::BLOB cds_sequence,1::UTINYINT codon_table FROM models ORDER BY transcript_index',
    'SELECT transcript_index::UINTEGER transcript_index,start::UBIGINT exon_start,"end"::UBIGINT exon_end,
     cs::UBIGINT exon_cdna_start,ce::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase FROM exons ORDER BY transcript_index,cs')
  stopifnot(DBI::dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('probe',",
    paste(q(queries), collapse = ','), ')'))$loaded)
  calls <- "SELECT event_index,r.seq_region,position,reference,string_split(alt,',') alternates,m.transcript_index,
    s.i sample_index,CASE s.i WHEN 0 THEN s0 WHEN 1 THEN s1 ELSE s2 END gt
    FROM records r JOIN models m USING(seq_region),range(3) s(i)"
  actual <- DBI::dbGetQuery(con, paste0('SELECT * FROM duckvep_haplotypes(', q(calls),
    ",'probe',input_mode:='source_records',phase_policy:='vep116_compat')"))
  saveRDS(actual, file.path(out, 'actual.rds'))
  oracle <- lapply(readLines(file.path(out, 'oracle.stdout')), jsonlite::fromJSON, simplifyVector = FALSE)
  names(oracle) <- vapply(oracle, `[[`, '', 'transcript')
  stopifnot(!anyDuplicated(names(oracle)), setequal(names(oracle), models$transcript))
  phase <- lapply(readLines(file.path(out, 'phase.jsonl')), jsonlite::fromJSON, simplifyVector = FALSE)
  names(phase) <- vapply(phase, `[[`, '', 'transcript')
  stopifnot(!anyDuplicated(names(phase)), setequal(names(phase), models$transcript))
  rows_by_tx <- split(seq_len(nrow(actual)), actual$transcript_index)
  records_by_region <- split(seq_len(nrow(records)), records$seq_region)
  exons_by_tx <- split(seq_len(nrow(exons)), exons$transcript_index)
  stopifnot(all(actual$carrier_count == vapply(actual$carriers, nrow, 1L)),
    setequal(names(rows_by_tx), as.character(models$transcript_index)))
  summary <- cbind(models[setdiff(names(models), 'cds')],
    cases[rep(seq_len(nrow(cases)), each = 2L), setdiff(names(cases), 'strand')])
  comparisons <- lapply(seq_len(nrow(models)), function(i) {
    model <- models[i, ]
    a <- actual[rows_by_tx[[as.character(model$transcript_index)]], , drop = FALSE]
    r <- records[records_by_region[[as.character(model$seq_region)]], ]
    spans <- exons[exons_by_tx[[as.character(model$transcript_index)]], ]
    expected <- oracle[[model$transcript]]$haplotypes
    observed <- lapply(seq_len(nrow(a)), function(j) {
      ids <- suppressWarnings(as.numeric(unlist(a$coding_blocks[[j]]$event_indices, use.names = FALSE)))
      matched <- match(ids, r$event_index)
      stopifnot(!anyNA(matched))
      list(cds = a$cds[j], protein = a$protein[j], count = a$carrier_count[j],
        contributors = r$source_id[matched],
        samples = as.list(table(paste0('s', a$carriers[[j]]$sample_index))))
    })
    list(expected = canonical(expected, samples = TRUE), observed = canonical(observed, samples = TRUE),
      equal = identical(canonical(expected, samples = TRUE), canonical(observed, samples = TRUE)),
      sequences_equal = identical(canonical(expected, FALSE), canonical(observed, FALSE)),
      counts_equal = sum(a$carrier_count) == oracle[[model$transcript]]$total_haplotype_count,
      unavailable_carriers = sum(a$carrier_count[is.na(a$cds)]),
      input_provenance_equal = input_provenance_ok(a, r),
      mappings_equal = mappings_ok(phase[[model$transcript]], model, spans, r))
  })
  for (name in c('equal', 'sequences_equal', 'counts_equal', 'input_provenance_equal',
      'mappings_equal', 'unavailable_carriers'))
    summary[[name]] <- vapply(comparisons, `[[`, if (name == 'unavailable_carriers') 0 else TRUE, name)
  summary$passed <- with(summary, equal & counts_equal & input_provenance_equal & mappings_equal)
  saveRDS(comparisons, file.path(out, 'comparisons.rds'))
  write.csv(summary, file.path(out, 'summary.csv'), row.names = FALSE)
  controls <- observation_controls(oracle[['T1full']]$haplotypes, actual[rows_by_tx[['0']], ],
    phase[['T1full']], models[1L, ], exons[exons_by_tx[['0']], ], records[records_by_region[['0']], ])
  write.csv(data.frame(control = names(controls), rejected = controls), file.path(out, 'controls.csv'), row.names = FALSE)
  DBI::dbDisconnect(con, shutdown = TRUE)
  identities <- unique(c('test/duckvep/conformance/haplotype_model_differential.R',
    'test/duckvep/conformance/haplotype_observations.R', 'scripts/duckvep_evidence.R', extension,
    'test/duckvep/conformance/haplotype_oracle.pl', paths[['haplotype_benchmark_reference']],
    'r/duckhtsbench/inst/benchmark_registry.tsv', list.files(out, full.names = TRUE)))
  jsonlite::write_json(list(source_revision = revision, extension_build_binding = binding,
    tracked_changes = duckvep_evidence_tracked_changes(root),
    scope = 'shared_transcript_source_replay_sequence_sample_counts_and_source_identity_sets',
    seed = seed, rare_per_stratum = opt$rare_per_stratum, required_strata = nrow(coverage),
    minimum_stratum_draws = min(coverage$observed), controls_rejected = sum(controls),
    oracle_revisions = as.list(pins), source_artifact = 'haplotype_benchmark_reference',
    profiles = nrow(summary), records = nrow(records), leaves = nrow(actual),
    carriers = sum(actual$carrier_count), threads = 4L, failures = sum(!summary$passed),
    sequence_failures = sum(!summary$sequences_equal), count_failures = sum(!summary$counts_equal),
    input_provenance_failures = sum(!summary$input_provenance_equal),
    mapping_failures = sum(!summary$mappings_equal),
    sha256 = as.list(vapply(identities, duckvep_evidence_sha256, ''))),
    file.path(out, 'receipt.json'), auto_unbox = TRUE, pretty = TRUE)
  if (!is.null(opt$extension_receipt)) duckvep_evidence_assert_checkout(root, revision)
  print(aggregate(cbind(profiles = rep(1L, nrow(summary)), failures = as.integer(!summary$passed),
    sequence_failures = as.integer(!summary$sequences_equal),
    input_provenance_failures = as.integer(!summary$input_provenance_equal),
    mapping_failures = as.integer(!summary$mappings_equal)) ~ cohort + shape + kind,
    summary, sum), row.names = FALSE)
  stopifnot(all(controls))
  if (any(!summary$passed)) stop('Shared-transcript replay differences retained: ', out, call. = FALSE)
}
main()
