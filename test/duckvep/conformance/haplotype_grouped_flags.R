#!/usr/bin/env Rscript
# Reproduce sequence-group metadata under controlled Perl sample iteration.
source('scripts/duckvep_evidence.R')
source('test/duckvep/conformance/haplotype_observations.R')
main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option('--vep-prefix', dest = 'prefix', default = '/root/miniconda3/envs/vep'),
    optparse::make_option('--oracle', default = 'test/duckvep/conformance/haplotype_grouped_flags.pl'))))
  root <- normalizePath('.')
  pins <- c(vep = '57ea5c52340acc1f156267f810ad162e26597082', variation = '2fb834b987ede3824e200197a838ce11e91aeb4b')
  mirrors <- normalizePath(c('.sync/ensembl-vep', '.sync/ensembl-variation'))
  for (i in seq_along(pins)) stopifnot(
    identical(duckvep_evidence_command('git', c('-C', mirrors[i], 'rev-parse', 'HEAD'), 'oracle revision'), unname(pins[i])),
    !length(duckvep_evidence_command('git', c('-C', mirrors[i], 'status', '--porcelain'), 'oracle checkout')))
  prefix <- normalizePath(opt$prefix)
  environment <- duckvep_evidence_command('micromamba', c('list', '-p', prefix, '--explicit'), 'oracle environment')
  stopifnot(identical(duckvep_evidence_explicit_packages(environment),
    duckvep_evidence_explicit_packages(readLines('test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt'))))
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(root, 'r/duckhtsbench/inst/benchmark_registry.tsv'))
  paths <- duckhtsbench::duckhts_bench_stage_repository_fixtures(root, 'duckvep-haplotypes')
  reference <- paths[['haplotype_benchmark_reference']]
  cds <- readLines(reference)[2L]
  stopifnot(nchar(cds) == 180L, substring(cds, 178L) == 'TAA')
  spans <- data.frame(start = c(11L, 46L, 116L, 125L, 165L, 208L, 256L),
    end = c(37L, 105L, 122L, 141L, 204L, 235L, 256L))
  lengths <- with(spans, end - start + 1L)
  stopifnot(sum(lengths) == 180L)
  offsets <- c(0L, head(cumsum(lengths), -1L))
  genome <- strrep('A', 266L)
  for (i in seq_len(nrow(spans))) substr(genome, spans$start[i], spans$end[i]) <-
    substring(cds, offsets[i] + 1L, offsets[i] + lengths[i])
  out <- tempfile('haplotype_grouped_flags_', tmpdir = 'test/duckvep/conformance/results')
  dir.create(out)
  message('Grouped-flag artifacts: ', out)
  writeLines(environment, file.path(out, 'environment.txt'))
  writeLines(c('>chrG13940', genome), file.path(out, 'reference.fa'))
  gff <- c('##gff-version 3', paste('chrG13940', 'geometry', c('gene', 'mRNA'), 11, 256, '.', '+', '.',
    c('ID=gene:TG13940;biotype=protein_coding',
      'ID=transcript:TG13940;Parent=gene:TG13940;biotype=protein_coding'), sep = '\t'))
  for (i in seq_len(nrow(spans))) gff <- c(gff, paste('chrG13940', 'geometry', c('exon', 'CDS'),
    spans$start[i], spans$end[i], '.', '+', c('.', (3L - offsets[i] %% 3L) %% 3L),
    c(paste0('ID=exon:TG13940_', i, ';Parent=transcript:TG13940'), 'Parent=transcript:TG13940'), sep = '\t'))
  writeLines(gff, file.path(out, 'model.gff3'))
  writeLines(c('##fileformat=VCFv4.4', '##contig=<ID=chrG13940,length=266>',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2',
    'chrG13940\t256\tedge\tA\tAGC\t.\tPASS\t.\tGT\t1|0\t0|1\t.|1',
    'chrG13940\t256\tanchor\tA\tC\t.\tPASS\t.\tGT\t1|1\t1|1\t1|1'), file.path(out, 'calls.vcf'))
  run <- function(command, args, name) stopifnot(system2(command, shQuote(args),
    stdout = file.path(out, paste0(name, '.stdout')), stderr = file.path(out, paste0(name, '.stderr'))) == 0L)
  run('samtools', c('faidx', file.path(out, 'reference.fa')), 'faidx')
  run('bcftools', c('norm', '-c', 'e', '-f', file.path(out, 'reference.fa'), '-o',
    file.path(out, 'checked.vcf'), file.path(out, 'calls.vcf')), 'reference_check')
  run('bgzip', file.path(out, 'model.gff3'), 'bgzip')
  run('tabix', c('-p', 'gff', file.path(out, 'model.gff3.gz')), 'tabix')
  libs <- paste(c(file.path(mirrors, 'modules'), file.path(prefix, 'share/ensembl-vep-116.0-0')), collapse = ':')
  extension <- normalizePath('build/release/duckhts.duckdb_extension')
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = 'true')))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbExecute(con, paste('LOAD', q(extension)))
  DBI::dbExecute(con, 'SET threads=1')
  DBI::dbWriteTable(con, 'grouped_exons', data.frame(spans,
    cs = offsets + 1L, ce = offsets + lengths,
    phase = offsets %% 3L, end_phase = (offsets + lengths) %% 3L))
  queries <- c('SELECT 0::UINTEGER seq_region', paste0(
    'SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,11::UBIGINT transcript_start,
     256::UBIGINT transcript_end,1::TINYINT strand,0::UINTEGER gene_index,3::UBIGINT transcript_flags,
     11::UBIGINT cds_start,256::UBIGINT cds_end,', q(cds), '::BLOB cds_sequence,1::UTINYINT codon_table'),
    'SELECT 0::UINTEGER transcript_index,start::UBIGINT exon_start,"end"::UBIGINT exon_end,
     cs::UBIGINT exon_cdna_start,ce::UBIGINT exon_cdna_end,phase::TINYINT phase,end_phase::TINYINT end_phase
     FROM grouped_exons ORDER BY cs')
  stopifnot(DBI::dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('grouped',",
    paste(q(queries), collapse = ','), ')'))$loaded)
  calls <- paste0('SELECT record_index event_index,0::UINTEGER seq_region,POS AS position,REF AS reference,
    ALT alternates,0::UINTEGER transcript_index,c.sample_index,c.raw_gt gt
    FROM (SELECT record_index,POS,REF,ALT,unnest(calls) c FROM read_geno(',
    q(normalizePath(file.path(out, 'calls.vcf'))), ',raw_gt:=true))')
  native_query <- paste0('SELECT * FROM duckvep_haplotypes(', q(calls),
    ",'grouped',input_mode:='source_records',phase_policy:='vep116_compat')")
  writeLines(native_query, file.path(out, 'native_query.sql'))
  actual <- DBI::dbGetQuery(con, native_query)
  saveRDS(actual, file.path(out, 'actual.rds'))
  stopifnot(nrow(actual) == 4L, sum(actual$carrier_count) == 6L)
  rows <- list()
  canonical_lanes <- function(x) {
    x <- lapply(x, function(lane) lane[setdiff(names(lane), 'traversal_ordinal')])
    x[order(vapply(x, `[[`, '', 'sample'), vapply(x, `[[`, 1L, 'lane1'))]
  }
  first_lanes <- NULL
  metadata_controls <- NULL
  for (seed in 0:31) for (repeat_index in 1:2) {
    name <- sprintf('seed%02d_repeat%d', seed, repeat_index)
    lane_path <- file.path(out, paste0(name, '.lanes.jsonl'))
    for (mode in c('plain', 'observed')) run('micromamba', c('run', '--clean-env',
      '--env', paste0('PERL5LIB=', libs), '--env', paste0('PERL_HASH_SEED=', seed),
      '--env', 'PERL_PERTURB_KEYS=0', '-p', prefix, 'perl', opt$oracle,
      file.path(out, c('calls.vcf', 'reference.fa', 'model.gff3.gz')),
      if (mode == 'observed') lane_path), paste(name, mode, sep = '.'))
    plain <- readLines(file.path(out, paste0(name, '.plain.stdout')))
    observed <- readLines(file.path(out, paste0(name, '.observed.stdout')))
    stopifnot(length(plain) == 1L, identical(plain, observed))
    x <- jsonlite::fromJSON(plain, simplifyVector = FALSE)
    lanes <- lapply(readLines(lane_path), jsonlite::fromJSON, simplifyVector = FALSE)
    group_observation <- jsonlite::fromJSON(readLines(paste0(lane_path, '.groups.jsonl')),
      simplifyVector = FALSE)
    observation <- list(replay_lanes = lanes, cds_group_metadata = group_observation$groups,
      reference_cds = group_observation$reference_cds)
    stopifnot(group_observation$transcript == 'TG13940',
      haplotype_group_metadata_equal(observation, x, container_json = TRUE),
      native_lane_flags_equal(lanes, actual, paste0('s', 0:2)))
    stopifnot(length(lanes) == 6L, x$transcript_id == 'TG13940', x$total_haplotype_count == 6L,
      length(x$cds_haplotypes) == 2L, length(x$protein_haplotypes) == 1L)
    ordered <- canonical_lanes(lanes)
    if (is.null(first_lanes)) first_lanes <- ordered else stopifnot(identical(first_lanes, ordered))
    if (is.null(metadata_controls)) {
      metadata_controls <- haplotype_metadata_controls(observation, x, container_json = TRUE)
      ordered_metadata <- haplotype_lane_metadata(lanes)
      same_cds <- Filter(function(lane) nchar(lane$cds) == 180L, ordered_metadata)
      first <- same_cds[[1L]]
      other <- Filter(function(lane) !identical(haplotype_raw_flags(lane$flags),
        haplotype_raw_flags(first$flags)), same_cds)[[1L]]
      wrong_order <- observation
      ordinals <- vapply(lanes, `[[`, 1L, 'traversal_ordinal')
      wrong_order$replay_lanes[[match(first$traversal_ordinal, ordinals)]]$traversal_ordinal <- other$traversal_ordinal
      wrong_order$replay_lanes[[match(other$traversal_ordinal, ordinals)]]$traversal_ordinal <- first$traversal_ordinal
      stopifnot(!is.null(haplotype_lane_metadata(wrong_order$replay_lanes)))
      metadata_controls <- c(metadata_controls, metadata_wrong_first_lane =
        !haplotype_group_metadata_equal(wrong_order, x, container_json = TRUE))
      for (bit in c(1L, 2L, 4L)) {
        wrong <- actual
        wrong$sequence_flags[1L] <- bitwXor(wrong$sequence_flags[1L], bit)
        metadata_controls[paste0('metadata_native_bit_', bit)] <-
          !native_lane_flags_equal(lanes, wrong, paste0('s', 0:2))
      }
      stopifnot(all(metadata_controls))
    }
    rows[[name]] <- do.call(rbind, lapply(x$cds_haplotypes, function(h) {
      first <- Filter(function(lane) identical(lane$cds, h$seq), lanes)[[1L]]
      stopifnot(!is.null(h$has_indel), length(h$has_indel) == 1L, h$has_indel %in% 0:1,
        h$has_indel == as.integer(isTRUE(first$flags$indel == 1L)),
        h$count == 3L, identical(unlist(h$samples)[sort(names(h$samples))], c(s0 = 1L, s1 = 1L, s2 = 1L)))
      raw <- haplotype_raw_flags(first$flags)
      data.frame(hash_seed = seed, repeat_index = repeat_index, cds = h$seq, has_indel = h$has_indel,
        first_sample = first$sample, first_lane = first$lane1, count = h$count,
        first_traversal_ordinal = first$traversal_ordinal, raw_indel = unname(raw['indel']),
        raw_frameshift = unname(raw['frameshift']), nominal_length_diff = unname(raw['length_diff']),
        categories = paste(haplotype_flag_categories(haplotype_raw_flag_bits(first$flags)), collapse = ','))
    }))
  }
  rows <- do.call(rbind, rows)
  write.csv(rows, file.path(out, 'observations.csv'), row.names = FALSE)
  write.csv(data.frame(control = names(metadata_controls), rejected = metadata_controls),
    file.path(out, 'metadata_controls.csv'), row.names = FALSE)
  a <- subset(rows, repeat_index == 1L); b <- subset(rows, repeat_index == 2L)
  a <- a[order(a$hash_seed, a$cds), ]; b <- b[order(b$hash_seed, b$cds), ]
  stopifnot(all(vapply(setdiff(names(a), 'repeat_index'), function(k) identical(a[[k]], b[[k]]), TRUE)),
    setequal(a$has_indel[nchar(a$cds) == 180L], 0:1))
  identities <- c(reference, opt$oracle, 'test/duckvep/conformance/haplotype_grouped_flags.R',
    'test/duckvep/conformance/haplotype_observations.R', extension,
    file.path(mirrors[2L], 'modules/Bio/EnsEMBL/Variation/TranscriptHaplotypeContainer.pm'),
    list.files(out, full.names = TRUE))
  jsonlite::write_json(list(source_revision = duckvep_evidence_revision(root),
    tracked_changes = duckvep_evidence_tracked_changes(root), oracle_revisions = as.list(pins),
    source_artifact = 'haplotype_benchmark_reference', hash_seeds = 0:31, repeats = 2L,
    perturb_keys = 0L, original_runs = 64L, observed_runs = 64L, full_output_unchanged = TRUE,
    per_lane_flags_unchanged = TRUE, first_lane_flags_match = TRUE,
    lane_flag_comparisons = 384L, cds_group_metadata_comparisons = 128L,
    lane_flag_failures = 0L, group_metadata_failures = 0L,
    metadata_controls_rejected = sum(metadata_controls), extension_build_binding = 'diagnostic_unbound',
    sha256 = as.list(vapply(identities, duckvep_evidence_sha256, ''))),
    file.path(out, 'receipt.json'), pretty = TRUE, auto_unbox = TRUE)
  print(table(cds_length = nchar(rows$cds), has_indel = rows$has_indel))
}
main()
