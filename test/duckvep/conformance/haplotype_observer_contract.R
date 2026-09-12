#!/usr/bin/env Rscript
# Verify transcript-owned observations against a shared-source, two-exon fixture.

main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option('--vep-prefix', dest = 'vep_prefix', default = '/root/miniconda3/envs/vep'),
    optparse::make_option('--oracle', default = 'test/duckvep/conformance/haplotype_oracle.pl')
  )))
  source('scripts/duckvep_evidence.R', local = TRUE)
  source('test/duckvep/conformance/haplotype_observations.R', local = TRUE)
  root <- normalizePath('.')
  prefix <- normalizePath(opt$vep_prefix)
  pins <- c(vep = '57ea5c52340acc1f156267f810ad162e26597082',
    variation = '2fb834b987ede3824e200197a838ce11e91aeb4b')
  mirrors <- normalizePath(c('.sync/ensembl-vep', '.sync/ensembl-variation'))
  for (i in seq_along(pins)) stopifnot(
    identical(duckvep_evidence_command('git', c('-C', mirrors[i], 'rev-parse', 'HEAD'), 'oracle revision'), unname(pins[i])),
    !length(duckvep_evidence_command('git', c('-C', mirrors[i], 'status', '--porcelain'), 'oracle checkout')))
  environment <- duckvep_evidence_command('micromamba', c('list', '-p', prefix, '--explicit'), 'oracle environment')
  stopifnot(identical(duckvep_evidence_explicit_packages(environment),
    duckvep_evidence_explicit_packages(readLines('test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt'))))
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(root, 'r/duckhtsbench/inst/benchmark_registry.tsv'))
  paths <- duckhtsbench::duckhts_bench_stage_repository_fixtures(root, 'duckvep-haplotypes')
  reference <- paths[['haplotype_benchmark_reference']]
  cds <- readLines(reference)[2L]
  stopifnot(nchar(cds) == 180L)
  reverse <- paste(rev(strsplit(chartr('ACGT', 'TGCA', cds), '', fixed = TRUE)[[1L]]), collapse = '')
  genome <- paste0(strrep('A', 10L), c(cds, reverse), strrep('A', 10L))
  out <- tempfile('haplotype_observer_', tmpdir = 'test/duckvep/conformance/results')
  dir.create(out)
  message('Observer artifacts: ', out)
  writeLines(environment, file.path(out, 'environment.txt'))
  chr <- c('chrForward', 'chrReverse')
  positions <- c(105L, 40L)
  ref <- substring(genome, positions, positions)
  alt <- chartr('ACGT', 'CGTA', ref)
  writeLines(as.vector(rbind(paste0('>', chr), genome)), file.path(out, 'reference.fa'))
  gff <- character()
  expected <- data.frame(transcript = c('F_full', 'F_spliced', 'R_full', 'R_spliced'),
    position = c(95L, 65L, 151L, 121L), alternate = rep(alt, each = 2L))
  for (i in 1:2) for (kind in c('full', 'spliced')) {
    tx <- paste(if (i == 1L) 'F' else 'R', kind, sep = '_')
    spans <- if (kind == 'full') data.frame(start = 11L, end = 190L) else
      data.frame(start = c(11L, 101L), end = c(70L, 190L))
    features <- data.frame(type = c('gene', 'mRNA'), start = 11L, end = 190L, phase = '.',
      attributes = c(paste0('ID=gene:', tx, ';biotype=protein_coding'),
        paste0('ID=transcript:', tx, ';Parent=gene:', tx, ';biotype=protein_coding')))
    for (j in seq_len(nrow(spans))) features <- rbind(features, data.frame(type = c('exon', 'CDS'),
      start = spans$start[j], end = spans$end[j], phase = c('.', '0'),
      attributes = c(paste0('ID=exon:', tx, j, ';Parent=transcript:', tx), paste0('Parent=transcript:', tx))))
    gff <- c(gff, with(features, paste(chr[i], 'observer', type, start, end, '.',
      if (i == 1L) '+' else '-', phase, attributes, sep = '\t')))
  }
  fields <- strsplit(gff, '\t')
  gff <- gff[order(vapply(fields, `[[`, '', 1L), as.integer(vapply(fields, `[[`, '', 4L)))]
  writeLines(c('##gff-version 3', gff), file.path(out, 'model.gff3'))
  writeLines(c('##fileformat=VCFv4.4', paste0('##contig=<ID=', chr, ',length=200>'),
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample',
    paste(chr, positions, 'site', ref, alt, '.', 'PASS', '.', 'GT', '1|1', sep = '\t')),
    file.path(out, 'calls.vcf'))
  run <- function(command, args, name) {
    status <- system2(command, shQuote(args), stdout = file.path(out, paste0(name, '.stdout')),
      stderr = file.path(out, paste0(name, '.stderr')))
    if (status != 0L) stop(name, ' failed; artifacts retained: ', out, call. = FALSE)
  }
  run('samtools', c('faidx', file.path(out, 'reference.fa')), 'faidx')
  run('bcftools', c('norm', '-c', 'e', '-f', file.path(out, 'reference.fa'), '-o',
    file.path(out, 'checked.vcf'), file.path(out, 'calls.vcf')), 'reference_check')
  run('bgzip', file.path(out, 'model.gff3'), 'bgzip')
  run('tabix', c('-p', 'gff', file.path(out, 'model.gff3.gz')), 'tabix')
  libs <- paste(c(file.path(mirrors, 'modules'), file.path(prefix, 'share/ensembl-vep-116.0-0')), collapse = ':')
  oracle <- normalizePath(opt$oracle)
  for (mode in c('plain', 'observed')) run('micromamba', c('run', '--clean-env',
    '--env', paste0('PERL5LIB=', libs), '-p', prefix, 'perl', oracle,
    file.path(out, c('calls.vcf', 'reference.fa', 'model.gff3.gz')),
    if (mode == 'observed') file.path(out, 'phase.jsonl')), mode)
  plain <- readLines(file.path(out, 'plain.stdout'))
  observed <- readLines(file.path(out, 'observed.stdout'))
  stopifnot(length(plain) == 4L, identical(sort(plain), sort(observed)))
  phase <- lapply(readLines(file.path(out, 'phase.jsonl')), jsonlite::fromJSON, simplifyVector = FALSE)
  names(phase) <- vapply(phase, `[[`, '', 'transcript')
  check <- function(value) {
    if (length(value) != 4L || anyDuplicated(names(value)) ||
        !setequal(names(value), expected$transcript)) return(FALSE)
    for (i in seq_len(nrow(expected))) {
      calls <- value[[expected$transcript[i]]]$calls
      if (length(calls) != 1L) return(FALSE)
      call <- calls[[1L]]
      if (!identical(call$source_id, 'site') || !identical(call$sample, 'sample') ||
          !identical(as.numeric(call$mapping_start), as.numeric(expected$position[i])) ||
          !identical(as.numeric(call$mapping_end), as.numeric(expected$position[i])) ||
          !identical(call$genotype, rep(list(expected$alternate[i]), 2L))) return(FALSE)
    }
    TRUE
  }
  stopifnot(check(phase))
  stale <- missing <- dropped <- phase
  stale$F_spliced$calls[[1L]]$mapping_start <- 95L
  stale$F_spliced$calls[[1L]]$mapping_end <- 95L
  missing$F_full$calls[[1L]]$mapping_start <- NULL
  missing$F_full$calls[[1L]]$mapping_end <- NULL
  dropped$R_full$calls <- list()
  stopifnot(!check(stale), !check(missing), !check(dropped))
  output <- lapply(plain, jsonlite::fromJSON, simplifyVector = FALSE)
  names(output) <- vapply(output, `[[`, '', 'transcript')
  metadata_controls <- haplotype_metadata_controls(phase[['F_full']], output[['F_full']])
  for (i in seq_len(nrow(expected))) {
    tx <- expected$transcript[i]
    haplotypes <- output[[tx]]$haplotypes
    stopifnot(length(haplotypes) == 1L, haplotypes[[1L]]$count == 2L)
    region <- if (startsWith(tx, 'F')) 1L else 2L
    source_key <- paste(positions[region], positions[region],
      paste(ref[region], alt[region], sep = '/'), sep = '_')
    allele <- if (region == 1L) alt[region] else chartr('ACGT', 'TGCA', alt[region])
    lanes <- lapply(1:2, function(lane) list(sample = 'sample', lane1 = lane,
      cds = haplotypes[[1L]]$cds, protein = haplotypes[[1L]]$protein,
      applied_sources = list(list(allele_key = paste0(allele, '|', source_key),
        source_id = 'site', source_key = source_key))))
    stopifnot(replay_lanes_equal(lanes, phase[[tx]]$replay_lanes))
    stopifnot(haplotype_group_metadata_equal(phase[[tx]], output[[tx]]))
    wrong <- lanes
    wrong[[1L]]$applied_sources[[1L]]$allele_key <- 'wrong'
    stopifnot(!replay_lanes_equal(lanes, lanes[-1L]),
      !replay_lanes_equal(lanes, c(lanes, lanes[1L])), !replay_lanes_equal(lanes, wrong))
  }
  for (mode in c('plain', 'observed')) run('micromamba', c('run', '--clean-env',
    '--env', paste0('PERL5LIB=', libs), '--env', 'PERL_HASH_SEED=0',
    '--env', 'PERL_PERTURB_KEYS=0', '-p', prefix, 'perl', oracle, '--container-json',
    file.path(out, c('calls.vcf', 'reference.fa', 'model.gff3.gz')),
    if (mode == 'observed') file.path(out, 'container_phase.jsonl')), paste0('container_', mode))
  container_plain <- readLines(file.path(out, 'container_plain.stdout'))
  container_observed <- readLines(file.path(out, 'container_observed.stdout'))
  stopifnot(length(container_plain) == 4L, identical(container_plain, container_observed))
  container_phase <- lapply(readLines(file.path(out, 'container_phase.jsonl')),
    jsonlite::fromJSON, simplifyVector = FALSE)
  container_output <- lapply(container_plain, jsonlite::fromJSON, simplifyVector = FALSE)
  names(container_phase) <- vapply(container_phase, `[[`, '', 'transcript')
  for (x in container_output) stopifnot(haplotype_group_metadata_equal(
    container_phase[[x$transcript_id]], x, container_json = TRUE))
  first <- container_output[[1L]]
  metadata_controls <- c(metadata_controls, setNames(haplotype_metadata_controls(
    container_phase[[first$transcript_id]], first, container_json = TRUE),
    paste0('container_', names(metadata_controls))))
  reference_observation <- phase[['F_full']]
  reference_observation$replay_lanes <- list()
  reference_observation$reference_cds <- reference_observation$cds_group_metadata[[1L]]$cds
  stopifnot(haplotype_group_metadata_equal(reference_observation, output[['F_full']]))
  reference_observation$reference_cds <- paste0(reference_observation$reference_cds, 'A')
  metadata_controls <- c(metadata_controls, metadata_unobserved_nonreference_group =
    !haplotype_group_metadata_equal(reference_observation, output[['F_full']]))
  stopifnot(all(metadata_controls))
  write.csv(data.frame(control = names(metadata_controls), rejected = metadata_controls),
    file.path(out, 'metadata_controls.csv'), row.names = FALSE)
  identities <- c(reference, oracle, 'test/duckvep/conformance/haplotype_observer_contract.R',
    'test/duckvep/conformance/haplotype_observations.R',
    list.files(out, full.names = TRUE))
  jsonlite::write_json(list(source_revision = duckvep_evidence_revision(root),
    tracked_changes = duckvep_evidence_tracked_changes(root), oracle_revisions = as.list(pins),
    source_artifact = 'haplotype_benchmark_reference', transcripts = 4L, source_records = 2L,
    controls_rejected = 15L, replay_lanes = 8L, full_output_unchanged = TRUE,
    metadata_controls_rejected = sum(metadata_controls), container_output_unchanged = TRUE,
    sha256 = as.list(vapply(identities, duckvep_evidence_sha256, ''))),
    file.path(out, 'receipt.json'), pretty = TRUE, auto_unbox = TRUE)
  message('Observer ownership: 4 transcript mappings, 8 replay lanes, unchanged complete output, 15 corruptions rejected')
  message('Metadata: ', sum(metadata_controls), ' corruptions rejected; unchanged original Runner JSON')
}
main()
