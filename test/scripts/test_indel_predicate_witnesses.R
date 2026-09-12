#!/usr/bin/env Rscript
# Eight unsigned VEP-116 diagnostic witnesses, not general indel/compound conformance.
source('scripts/duckvep_evidence.R')
receipt_pin <- 'd782e33f3ead75901fe178b820ea6d5797982594911c6f160f2acc11b1971db7'
oracle_pins <- c('57ea5c52340acc1f156267f810ad162e26597082', '2fb834b987ede3824e200197a838ce11e91aeb4b')
read_json <- function(text) {
  check <- function(x) {
    if (is.list(x)) {
      stopifnot(!anyDuplicated(names(x)))
      for (child in x) check(child)
    }
    x
  }
  check(jsonlite::fromJSON(text, simplifyVector = FALSE))
}
records <- function(text) lapply(strsplit(sub('\n$', '', text), '\n', fixed = TRUE)[[1L]], read_json)
same <- function(a, b) {
  if (!is.list(b)) return(identical(a, b))
  if (!is.list(a) || length(a) != length(b) || !setequal(names(a), names(b))) return(FALSE)
  if (!is.null(names(b))) a <- a[names(b)]
  all(vapply(seq_along(b), function(i) same(a[[i]], b[[i]]), FALSE))
}
hash_text <- function(text) {
  path <- tempfile('indel-text-')
  on.exit(unlink(path))
  writeChar(text, path, eos = NULL, useBytes = TRUE)
  duckvep_evidence_sha256(path)
}
expected_models <- function() {
  frame <- lapply(1:2, function(i) {
    id <- c('FRAME_LOCAL_GX', 'FRAME_LOCAL_PX')[i]
    cds <- c('ATGAAACCCGGGTTTTAA', 'ATGACGTATGTAGTAGATCCTTCCGAATAT')[i]
    list(id = id, table = 1L, cds = cds, genomic_sequence = paste0(strrep('A', 10L), cds, strrep('A', 10L)),
      cds_start1 = 11L, edits = list(), variant_format = 'vcf', variants = list(list(id = id,
        position1 = c(10L, 19L)[i], reference = c('GG', 'CCT')[i], alternate = c('C', 'A')[i])))
  })
  phase_models <- list()
  original <- 'CGTACGTACGTACGTACGTACGTTACGTACGTACGTACTGGTAA'
  cdna <- paste0(strrep('A', 34L), 'N', substring(original, 2L), 'ACGTACGTAC')
  for (strand in c(-1L, 1L)) for (phase in 0:2) {
    starts <- if (strand == 1L) c(100L, 150L, 220L) else c(225L, 170L, 100L)
    ends <- if (strand == 1L) c(125L, 180L, 250L) else c(250L, 200L, 130L)
    phases <- c(-1L, phase, (phase + 23L) %% 3L)
    exons <- lapply(1:3, function(i) list(start1 = starts[i], end1 = ends[i],
      cdna_start1 = c(1L, 27L, 58L)[i], cdna_end1 = c(26L, 57L, 88L)[i],
      phase = phases[i], end_phase = c(-1L, phases[3L], phases[3L])[i]))
    positions <- unlist(Map(function(lo, hi) if (strand == 1L) lo:hi else hi:lo, starts, ends))
    genome <- rep('A', 350L)
    genome[positions] <- strsplit(if (strand == 1L) cdna else chartr('ACGT', 'TGCA', cdna), '', fixed = TRUE)[[1L]]
    phase_models[[length(phase_models) + 1L]] <- list(
      id = paste0('LATER_', if (strand == 1L) 'PLUS' else 'MINUS', '_P', phase), strand = strand,
      phase = phase, exons = exons, original_unpadded_cds = original,
      padded_cds = paste0(strrep('N', phase + 1L), substring(original, 2L)), cdna = cdna,
      genomic_sequence = paste(genome, collapse = ''), source = list(
        position1 = if (strand == 1L) 160L else 190L, reference = if (strand == 1L) 'T' else 'A',
        alternate = if (strand == 1L) 'AC' else 'GT'))
  }
  list(raw_missense = frame, phase_n = phase_models)
}
check_family <- function(bundle, expected) {
  frame <- bundle$family == 'raw_missense'
  files <- bundle$files
  models <- records(files[[if (frame) 'cases.jsonl' else 'models.jsonl']])
  tva <- records(files[['tva.stdout']])
  stopifnot(same(models, expected), length(tva) == length(expected))
  if (frame) stopifnot(length(records(files[['observer.stdout']])) == 2L)
  source <- strsplit(strsplit(files[['input.vcf']], '\n', fixed = TRUE)[[1L]], '\t', fixed = TRUE)
  source <- source[!vapply(source, function(x) startsWith(x[1L], '#'), FALSE)]
  cli_lines <- strsplit(files[['vep.vcf']], '\n', fixed = TRUE)[[1L]]
  header <- grep('^##INFO=<ID=CSQ,', cli_lines, value = TRUE)
  stopifnot(length(header) == 1L)
  fields <- strsplit(sub('.*Format: (.*)".*', '\\1', header), '|', fixed = TRUE)[[1L]]
  cli <- strsplit(cli_lines[!startsWith(cli_lines, '#')], '\t', fixed = TRUE)
  stopifnot(length(source) == length(expected), length(cli) == length(expected), !anyDuplicated(fields),
    identical(files[['reference.fa']], paste0(paste(unlist(lapply(expected, function(x)
      c(paste0('>', x$id), x$genomic_sequence))), collapse = '\n'), '\n')))
  for (i in seq_along(expected)) {
    model <- expected[[i]]
    allele <- if (frame) model$variants[[1L]] else model$source
    position <- allele$position1 + if (frame) 10L else 0L
    stopifnot(substr(model$genomic_sequence, position, position + nchar(allele$reference) - 1L) == allele$reference,
      identical(source[[i]], c(model$id, as.character(position), model$id, allele$reference, allele$alternate, '.', 'PASS', '.')),
      length(cli[[i]]) == 8L, identical(cli[[i]][1:7], source[[i]][1:7]))
    suffix <- if (frame) c('p.Gly4ArgfsTer?', 'p.Pro7IlefsTer?')[i] else if (model$phase == 0L) 'p.Thr2_?1' else 'p.Met1?'
    ref <- if (frame) c('GGG', 'CCT')[i] else c('NGT', 'NNG', 'NNN')[model$phase + 1L]
    alt <- if (frame) c('CG', 'A')[i] else c('NGAC', 'NNAC', 'NNAC')[model$phase + 1L]
    terms <- if (frame) list('frameshift_variant') else list('frameshift_variant', 'start_lost')
    row <- list(id = model$id, source = allele, parser_start = position,
      parser_end = position + nchar(allele$reference) - 1L,
      parser_allele_string = paste(allele$reference, allele$alternate, sep = '/'), codon_table = 1L,
      ref_codon = ref, alt_codon = alt, ref_peptide = if (frame) c('G', 'P')[i] else 'X',
      alt_peptide = if (frame) 'X' else 'XX', consequences = terms,
      hgvsp = paste0(model$id, '_protein.1:', suffix))
    if (frame) row <- c(row, list(frameshift = i, missense = 1L, coding_unknown = 0L)) else {
      predicates <- setNames(as.list(rep(0L, 13L)), c('start_lost', 'start_retained_variant', 'stop_gained',
        'stop_lost', 'stop_retained', 'frameshift', 'inframe_insertion', 'inframe_deletion',
        'missense_variant', 'synonymous_variant', 'protein_altering_variant', 'partial_codon', 'coding_unknown'))
      predicates$start_lost <- predicates$frameshift <- 1L
      row <- c(row, list(strand = model$strand, phase = model$phase, prepared_cds = model$padded_cds,
        cdna = model$cdna, cds_start = 3L, cds_end = 3L, translation_start = 1L,
        translation_end = 1L, raw_predicates = predicates))
    }
    stopifnot(same(tva[[i]], row))
    csq <- strsplit(paste0(sub('^CSQ=', '', cli[[i]][8L]), '|END'), '|', fixed = TRUE)[[1L]]
    stopifnot(length(csq) == length(fields) + 1L, tail(csq, 1L) == 'END')
    csq <- setNames(head(csq, -1L), fields)
    stopifnot(csq['Allele'] == allele$alternate, csq['Feature'] == model$id,
      csq['Consequence'] == paste(unlist(terms), collapse = '&'), csq['HGVSp'] == paste0('.1:', suffix),
      csq['Amino_acids'] == paste(row$ref_peptide, row$alt_peptide, sep = '/'),
      toupper(csq['Codons']) == paste(ref, alt, sep = '/'))
    if (frame) {
      observed <- records(files[['observer.stdout']])[[i]]
      protein <- c('MKPGF*', 'MTYVVDPSEY')[i]
      stopifnot(same(observed, list(id = model$id, prepared_cds = model$cds,
        core_reference = c('MKPGF', 'MTYVVDPSEY')[i], reference = protein,
        alternate_full = protein, alternate = protein, independent_hgvs = list(list(
          id = model$id, allele = allele$alternate, consequences = terms, hgvsp = row$hgvsp,
          source_reference = allele$reference, source_alternate = allele$alternate,
          parser_start = row$parser_start, parser_end = row$parser_end, parser_allele_string = row$parser_allele_string)))))
    }
  }
  invisible(TRUE)
}
directory <- 'test/duckvep/conformance/data/indel_predicate_witnesses'
receipt_path <- file.path(directory, 'receipt.json')
duckvep_evidence_check_receipt_pin(receipt_path, receipt_pin)
receipt <- read_json(paste(readLines(receipt_path), collapse = '\n'))
stopifnot(receipt$format == 'indel_predicate_witnesses', receipt$source_binding == 'unsigned_diagnostic',
  receipt$scope == 'raw_indel_predicates_vs_emitted_SO_and_HGVSp_only_not_general_indel_or_compound_conformance',
  same(receipt$cases, list(raw_missense = 2L, phase_n = 6L)), identical(unlist(receipt$oracle_revisions), oracle_pins),
  identical(names(receipt$sha256), 'evidence.jsonl.gz'),
  setequal(list.files(directory), c('receipt.json', 'evidence.jsonl.gz')),
  duckvep_evidence_sha256(file.path(directory, 'evidence.jsonl.gz')) == receipt$sha256[['evidence.jsonl.gz']])
connection <- gzfile(file.path(directory, 'evidence.jsonl.gz'), 'rt')
bundles <- lapply(readLines(connection), read_json)
close(connection)
expected <- expected_models()
stopifnot(identical(vapply(bundles, `[[`, '', 'family'), names(expected)))
for (bundle in bundles) {
  original <- read_json(bundle$files[['receipt.json']])
  retained <- c('reference.fa', 'input.vcf', 'model.gff3', 'tva.stdout', 'tva.stderr', 'vep.vcf',
    'vep.stderr', 'vep.vcf_warnings.txt', 'observe.pl', 'run.R', 'receipt.json',
    if (bundle$family == 'raw_missense') c('cases.jsonl', 'observer.stdout', 'observer.stderr') else 'models.jsonl')
  stopifnot(setequal(names(bundle$files), retained),
    hash_text(bundle$files[['receipt.json']]) == receipt$original_receipt_sha256[[bundle$family]],
    identical(unlist(original$oracle_revisions), oracle_pins))
  for (file in setdiff(names(bundle$files), 'receipt.json'))
    stopifnot(hash_text(bundle$files[[file]]) == original$sha256[[file]])
  check_family(bundle, expected[[bundle$family]])
}

# Semantic controls bypass hashes deliberately; each must fail the actual comparison.
rejects <- function(expression) tryCatch({ force(expression); FALSE }, error = function(e) TRUE)
controls <- logical()
mutations <- list(wrong_ref = function(x) { x$source$reference <- 'A'; x },
  false_missense = function(x) { x$missense <- 0L; x }, wrong_key = function(x) { x$id <- 'wrong'; x },
  false_so = function(x) { x$consequences <- c(x$consequences, list('missense_variant')); x },
  wrong_codon = function(x) { x$alt_codon <- 'CCC'; x }, fractional = function(x) { x$frameshift <- 1.00000001; x })
for (name in names(mutations)) {
  bad <- bundles[[1L]]
  rows <- records(bad$files[['tva.stdout']])
  rows[[1L]] <- mutations[[name]](rows[[1L]])
  bad$files[['tva.stdout']] <- paste(vapply(rows, jsonlite::toJSON, '', auto_unbox = TRUE, digits = NA), collapse = '\n')
  controls[name] <- rejects(check_family(bad, expected$raw_missense))
}
for (family in seq_along(bundles)) for (kind in c('missing', 'duplicate', 'source', 'cli', 'raw_fact')) {
  bad <- bundles[[family]]
  if (kind %in% c('missing', 'duplicate')) {
    lines <- strsplit(bad$files[['tva.stdout']], '\n', fixed = TRUE)[[1L]]
    bad$files[['tva.stdout']] <- paste(if (kind == 'missing') lines[-1L] else c(lines, lines[1L]), collapse = '\n')
  } else if (kind == 'source') bad$files[['reference.fa']] <- sub('A', 'C', bad$files[['reference.fa']], fixed = TRUE)
  else if (kind == 'raw_fact') bad$files[['tva.stdout']] <- sub('"frameshift":1', '"frameshift":0', bad$files[['tva.stdout']], fixed = TRUE)
  else bad$files[['vep.vcf']] <- sub('frameshift_variant', 'missense_variant', bad$files[['vep.vcf']], fixed = TRUE)
  controls[paste(family, kind)] <- rejects(check_family(bad, expected[[family]]))
}
controls['duplicate_json_key'] <- rejects(read_json('{"id":"a","id":"b"}'))
controls <- c(controls, local({
  temporary <- tempfile('indel-receipt-')
  payload <- tempfile('indel-payload-', fileext = '.gz')
  on.exit(unlink(c(temporary, payload)))
  lines <- readLines(receipt_path)
  writeLines(lines, temporary)
  duckvep_evidence_check_receipt_pin(temporary, receipt_pin)
  writeLines(sub(oracle_pins[1L], strrep('0', 40L), lines, fixed = TRUE), temporary)
  identity <- rejects(duckvep_evidence_check_receipt_pin(temporary, receipt_pin))
  changed <- bundles
  changed[[1L]]$files[['tva.stdout']] <- sub('"missense":1', '"missense":0', changed[[1L]]$files[['tva.stdout']], fixed = TRUE)
  stopifnot(changed[[1L]]$files[['tva.stdout']] != bundles[[1L]]$files[['tva.stdout']])
  connection <- gzfile(payload, 'wb')
  writeLines(vapply(changed, jsonlite::toJSON, '', auto_unbox = TRUE), connection)
  close(connection)
  writeLines(sub(receipt$sha256[['evidence.jsonl.gz']], duckvep_evidence_sha256(payload), lines, fixed = TRUE), temporary)
  c(identity_swap = identity, rehashed_payload = rejects(duckvep_evidence_check_receipt_pin(temporary, receipt_pin)))
}))
stopifnot(all(controls))
message('Indel predicate witnesses: 2 + 6 cases verified; ', length(controls), ' corruptions rejected')
