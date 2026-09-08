#!/usr/bin/env Rscript
# Reproduce sequence-group metadata under controlled Perl sample iteration.
source('scripts/duckvep_evidence.R')
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
  rows <- list()
  canonical_lanes <- function(x) x[order(vapply(x, `[[`, '', 'sample'), vapply(x, `[[`, 1L, 'lane1'))]
  first_lanes <- NULL
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
    stopifnot(length(lanes) == 6L, x$transcript_id == 'TG13940', x$total_haplotype_count == 6L,
      length(x$cds_haplotypes) == 2L, length(x$protein_haplotypes) == 1L)
    ordered <- canonical_lanes(lanes)
    if (is.null(first_lanes)) first_lanes <- ordered else stopifnot(identical(first_lanes, ordered))
    rows[[name]] <- do.call(rbind, lapply(x$cds_haplotypes, function(h) {
      first <- Filter(function(lane) identical(lane$cds, h$seq), lanes)[[1L]]
      stopifnot(!is.null(h$has_indel), length(h$has_indel) == 1L, h$has_indel %in% 0:1,
        h$has_indel == as.integer(isTRUE(first$flags$indel == 1L)),
        h$count == 3L, identical(unlist(h$samples)[sort(names(h$samples))], c(s0 = 1L, s1 = 1L, s2 = 1L)))
      data.frame(hash_seed = seed, repeat_index = repeat_index, cds = h$seq, has_indel = h$has_indel,
        first_sample = first$sample, first_lane = first$lane1, count = h$count)
    }))
  }
  rows <- do.call(rbind, rows)
  write.csv(rows, file.path(out, 'observations.csv'), row.names = FALSE)
  a <- subset(rows, repeat_index == 1L); b <- subset(rows, repeat_index == 2L)
  a <- a[order(a$hash_seed, a$cds), ]; b <- b[order(b$hash_seed, b$cds), ]
  stopifnot(all(vapply(setdiff(names(a), 'repeat_index'), function(k) identical(a[[k]], b[[k]]), TRUE)),
    setequal(a$has_indel[nchar(a$cds) == 180L], 0:1))
  identities <- c(reference, opt$oracle, 'test/duckvep/conformance/haplotype_grouped_flags.R',
    file.path(mirrors[2L], 'modules/Bio/EnsEMBL/Variation/TranscriptHaplotypeContainer.pm'),
    list.files(out, full.names = TRUE))
  jsonlite::write_json(list(source_revision = duckvep_evidence_revision(root),
    tracked_changes = duckvep_evidence_tracked_changes(root), oracle_revisions = as.list(pins),
    source_artifact = 'haplotype_benchmark_reference', hash_seeds = 0:31, repeats = 2L,
    perturb_keys = 0L, original_runs = 64L, observed_runs = 64L, full_output_unchanged = TRUE,
    per_lane_flags_unchanged = TRUE, first_lane_flags_match = TRUE,
    sha256 = as.list(vapply(identities, duckvep_evidence_sha256, ''))),
    file.path(out, 'receipt.json'), pretty = TRUE, auto_unbox = TRUE)
  print(table(cds_length = nchar(rows$cds), has_indel = rows$has_indel))
}
main()
