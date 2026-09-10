#!/usr/bin/env Rscript
# Reference-translation and curated-reference HGVS observations from pinned Ensembl.
main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--vep-prefix", dest = "vep_prefix",
      default = Sys.getenv("VEP_PREFIX", "/root/miniconda3/envs/vep")),
    optparse::make_option("--evidence-out", dest = "evidence_out", default = "",
      help = "Retain complete comparison pairs and raw oracle in a new directory")
  )))
  source("scripts/duckvep_evidence.R", local = TRUE)
  root <- normalizePath(".")
  revision <- duckvep_evidence_revision(root)
  mirrors <- setNames(normalizePath(c(".sync/ensembl-vep", ".sync/ensembl-variation")),
    c("vep", "variation"))
  revisions <- c(vep = "57ea5c52340acc1f156267f810ad162e26597082",
    variation = "2fb834b987ede3824e200197a838ce11e91aeb4b")
  for (name in names(revisions)) stopifnot(
    identical(system2("git", c("-C", mirrors[[name]], "rev-parse", "HEAD"), stdout = TRUE), revisions[[name]]),
    !length(system2("git", c("-C", mirrors[[name]], "status", "--porcelain"), stdout = TRUE)))
  prefix <- normalizePath(opt$vep_prefix, mustWork = TRUE)
  out <- tempfile("reference_translation_", tmpdir = "test/duckvep/conformance/results")
  dir.create(out)
  message("Reference-translation artifacts: ", out)
  conda <- function(...) do.call(blit::conda,
    c(as.list(duckvep_blit_quote(c(...))), list(conda = duckvep_blit_quote("micromamba"))))
  environment <- file.path(out, "environment.txt")
  stopifnot(blit::cmd_run(conda("list", "-p", prefix, "--explicit"),
    stdout = environment, stderr = "2>&1", verbose = FALSE) == 0L,
    identical(duckvep_evidence_explicit_packages(readLines(environment)),
      duckvep_evidence_explicit_packages(readLines(
        "test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt"))))
  sources <- c("test/duckvep/conformance/reference_translation_probe.c",
    "src/duckvep/kernel/src/duckvep_codon.c", "src/duckvep/kernel/src/duckvep_haplotype.c")
  code <- c(sources, "src/duckvep/kernel/src/duckvep_codon.h",
    "src/duckvep/kernel/src/duckvep_dna.h", "src/duckvep/kernel/src/duckvep_haplotype.h",
    "test/duckvep/conformance/reference_translation_oracle.pl",
    "test/duckvep/conformance/reference_translation_differential.R", "scripts/duckvep_evidence.R")
  code_hashes <- vapply(code, duckvep_evidence_sha256, "")
  # Explicit actual inputs are hashed; dirty/untracked code remains diagnostic.
  changed <- c(duckvep_evidence_tracked_changes(root), duckvep_evidence_untracked_build_inputs(root))
  binding <- if (length(changed)) "diagnostic_unbound" else "clean_checkout"
  shared <- file.path(out, paste0("reference_translation_probe", .Platform$dynlib.ext))
  stopifnot(system2(Sys.getenv("CC", "cc"), shQuote(c("-std=c99", "-O1", "-g", "-Wall", "-Wextra",
    "-fPIC", "-shared", "-I", "src/duckvep/kernel/src", sources, "-o", shared)),
    stdout = file.path(out, "compiler.log"), stderr = file.path(out, "compiler.log")) == 0L)
  dll <- dyn.load(shared)
  on.exit(dyn.unload(shared), add = TRUE)
  symbol <- getNativeSymbolInfo("duckhts_test_raw_translation", PACKAGE = dll)
  reference_symbol <- getNativeSymbolInfo("duckhts_test_reference_protein", PACKAGE = dll)
  translate <- function(cds, table) {
    x <- .C(symbol, as.character(cds), as.integer(table), peptide = raw(nchar(cds) %/% 3L + 1L),
      as.integer(nchar(cds) %/% 3L + 1L), status = integer(1L), facts = double(3L))
    stopifnot(x$status == 0L)
    full <- rawToChar(x$peptide[seq_len(x$facts[1L])])
    list(full = full, visible = substr(full, 1L, if (x$facts[2L]) x$facts[2L] else nchar(full)))
  }
  reference_protein <- function(cds, table, edits) {
    positions <- vapply(edits, `[[`, 1L, "position1")
    alternates <- charToRaw(paste0(vapply(edits, `[[`, "", "alternate"), collapse = ""))
    capacity <- nchar(cds) %/% 3L + 2L
    x <- .C(reference_symbol, as.character(cds), as.integer(table), positions, alternates,
      as.integer(length(edits)), peptide = raw(capacity), as.integer(capacity),
      status = integer(1L), length = double(1L), coding_peptide = raw(capacity),
      coding_facts = double(3L))
    stopifnot(x$status == 0L)
    list(protein = rawToChar(x$peptide[seq_len(x$length)]),
      coding = rawToChar(x$coding_peptide[seq_len(x$coding_facts[1L])]),
      coding_length = x$coding_facts[1L], first_stop = x$coding_facts[2L],
      unambiguous = x$coding_facts[3L])
  }
  # Exhaustive ACGTN triplets in each codon role, every supported table and
  # trailing-partial length. The matrix is declared independently of C tables.
  tables <- c(1:6, 9:14, 16, 21:31)
  triplets <- do.call(paste0, expand.grid(rep(list(strsplit("ACGTN", "")[[1L]]), 3L),
    stringsAsFactors = FALSE))
  grid <- expand.grid(table = tables, triplet = triplets,
    role = c("start", "internal", "terminal"), tail = c("", "A", "AA"), stringsAsFactors = FALSE)
  cases <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    cds <- switch(g$role, start = paste0(g$triplet, "GCCTAA"),
      internal = paste0("ATG", g$triplet, "GCCTAA"), terminal = paste0("ATGGCC", g$triplet))
    list(id = paste(g$table, g$triplet, g$role, nchar(g$tail), sep = "/"),
      family = g$role, cds = paste0(cds, g$tail), table = g$table, edits = list())
  })
  witnesses <- list(
    legitimate_start = list(cds = "CTGGCCTAA"),
    internal_stop = list(cds = "ATGTGAGCCTAA"),
    consensus_n = list(cds = "ATGGCNTAA"),
    initial_met = list(cds = "GTGGCCTAA", edits = list(
      list(code = "initial_met", position1 = 1L, alternate = "M"))),
    selenocysteine = list(cds = "ATGTGAGCCTAA", edits = list(
      list(code = "_selenocysteine", position1 = 2L, alternate = "U"))),
    substitution = list(cds = "ATGGCCTAA", edits = list(
      list(code = "amino_acid_sub", position1 = 2L, alternate = "V"))),
    internal_readthrough = list(cds = "ATGTGAGCCTAA", edits = list(
      list(code = "_stop_codon_rt", position1 = 2L, alternate = "W"))),
    terminal_readthrough = list(cds = "ATGGCCTAA", edits = list(
      list(code = "_stop_codon_rt", position1 = 3L, alternate = "W"))),
    mito_terminal_tga = list(cds = "ATGGCCTGA", table = 2L),
    mito_terminal_aga = list(cds = "ATGGCCAGA", table = 2L),
    partial_stop_spelling = list(cds = "ATGGCCTAAA"),
    lowercase_stop = list(cds = "atggcctaa"),
    lowercase_start = list(cds = "ctggcctaa"),
    multiple_edits = list(cds = "GTGTGAGCCTAA", edits = list(
      list(code = "initial_met", position1 = 1L, alternate = "M"),
      list(code = "_selenocysteine", position1 = 2L, alternate = "U"),
      list(code = "amino_acid_sub", position1 = 3L, alternate = "V")))
  )
  hgvs_witnesses <- c("legitimate_start", "selenocysteine", "terminal_readthrough",
    "mito_terminal_tga", "lowercase_start")
  for (name in names(witnesses)) {
    x <- witnesses[[name]]
    cases[[length(cases) + 1L]] <- list(id = paste0("witness/", name), family = "witness",
      cds = x$cds, table = if (is.null(x$table)) 1L else x$table,
      edits = if (is.null(x$edits)) list() else x$edits)
    if (name %in% hgvs_witnesses) cases[[length(cases)]]$variants <- list(list(
      id = "synonymous", position1 = if (name == "selenocysteine") 9L else 6L,
      reference = "C", alternate = "T"))
  }
  ids <- vapply(cases, `[[`, "", "id")
  stopifnot(!anyDuplicated(ids), nrow(grid) == 27000L, length(cases) == 27014L)
  input <- file.path(out, "cases.jsonl")
  writeLines(vapply(cases, jsonlite::toJSON, "", auto_unbox = TRUE), input)
  perl_lib <- paste(c(file.path(mirrors, "modules"),
    file.path(prefix, "share/ensembl-vep-116.0-0")), collapse = ":")
  oracle_file <- file.path(out, "oracle.jsonl")
  run <- conda("run", "--clean-env", "--env", paste0("PERL5LIB=", perl_lib), "-p", prefix,
    "perl", file.path(root, "test/duckvep/conformance/reference_translation_oracle.pl"), normalizePath(input))
  stopifnot(blit::cmd_run(run, stdout = oracle_file, stderr = file.path(out, "oracle.log"),
    stdin = NULL, verbose = FALSE) == 0L)
  oracle <- lapply(readLines(oracle_file), jsonlite::fromJSON)
  stopifnot(identical(vapply(oracle, `[[`, "", "id"), ids))
  independent <- do.call(rbind, lapply(hgvs_witnesses, function(name) {
    id <- paste0("witness/", name)
    observed <- oracle[[match(id, ids)]]$independent_hgvs
    stopifnot(is.data.frame(observed), all(observed$id == "synonymous"),
      all(unlist(observed$consequences) == "synonymous_variant"))
    data.frame(case_id = id, allele = observed$allele, hgvsp = observed$hgvsp)
  }))
  expected_hgvs <- do.call(rbind, lapply(hgvs_witnesses, function(name) data.frame(
    case_id = paste0("witness/", name), allele = if (name == "lowercase_start") c("C", "T") else "T",
    hgvsp = paste0("witness/", name, "_protein.1:p.Ala", if (name == "selenocysteine") 3 else 2, "="))))
  equal_hgvs <- function(x) identical(x, expected_hgvs)
  write.csv(independent, file.path(out, "independent_hgvs.csv"), row.names = FALSE)
  rows <- lapply(seq_along(cases), function(i) {
    case <- cases[[i]]; expected <- oracle[[i]]
    actual <- translate(expected$prepared_cds, case$table)
    reference <- reference_protein(expected$prepared_cds, case$table, case$edits)
    # Coding operands use the uncurated pinned consensus translation. Raw DNA
    # ambiguity is an independent fact, including the trailing partial codon.
    coding <- strsplit(expected$alternate_full, '', fixed = TRUE)[[1L]]
    stopifnot(length(coding) == nchar(expected$prepared_cds) %/% 3L)
    coding_stop <- which(coding == '*')
    expected_coding <- paste0(coding, collapse = '')
    expected_stop <- if (length(coding_stop)) coding_stop[1L] else 0
    coding_metadata_equal <- reference$coding_length == length(coding) &&
      reference$first_stop == expected_stop
    coding_metadata_equal <- coding_metadata_equal &&
      reference$unambiguous == !grepl('N', toupper(expected$prepared_cds), fixed = TRUE)
    data.frame(id = case$id, family = case$family, table = case$table,
      cds = case$cds, prepared_cds = expected$prepared_cds,
      core_reference = expected$core_reference, expected_reference = expected$reference,
      raw_reference = actual$full, actual_reference = reference$protein,
      expected_coding = expected_coding, actual_coding = reference$coding,
      expected_coding_length = length(coding), actual_coding_length = reference$coding_length,
      expected_coding_stop = expected_stop, actual_coding_stop = reference$first_stop,
      expected_coding_unambiguous = !grepl('N', toupper(expected$prepared_cds), fixed = TRUE),
      actual_coding_unambiguous = reference$unambiguous != 0,
      coding_failure = !identical(reference$coding, expected_coding) || !coding_metadata_equal,
      expected_alternate_full = expected$alternate_full,
      actual_alternate_full = actual$full, expected_alternate = expected$alternate,
      actual_alternate = actual$visible,
      reference_failure = !identical(reference$protein, expected$reference),
      raw_reference_hypothesis_failure = !identical(actual$full, expected$reference),
      alternate_full_failure = !identical(actual$full, expected$alternate_full),
      alternate_failure = !identical(actual$visible, expected$alternate))
  })
  pairs <- do.call(rbind, rows)
  write.csv(pairs, file.path(out, "pairs.csv"), row.names = FALSE)
  write.csv(pairs[pairs$family == "witness", ], file.path(out, "witnesses.csv"), row.names = FALSE)
  summary <- aggregate(cbind(cases = rep(1L, nrow(pairs)), reference_failures = pairs$reference_failure,
    coding_failures = pairs$coding_failure,
    raw_reference_hypothesis_failures = pairs$raw_reference_hypothesis_failure,
    alternate_full_failures = pairs$alternate_full_failure, alternate_failures = pairs$alternate_failure),
    pairs[c("family")], sum)
  write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
  # Deliberate corruption must fail the complete identity/sequence comparison.
  expected <- pairs[c("id", "expected_reference", "expected_alternate_full", "expected_alternate",
    "expected_coding", "expected_coding_length", "expected_coding_stop",
    "expected_coding_unambiguous")]
  equal <- function(x) identical(x, expected)
  controls <- c(missing = !equal(expected[-1L, ]), duplicate = !equal(rbind(expected, expected[1L, ])))
  for (field in names(expected)) {
    corrupt <- expected
    value <- corrupt[[field]][1L]
    corrupt[[field]][1L] <- if (is.logical(value)) !value else if (is.numeric(value))
      value + 1 else paste0(value, "X")
    controls[field] <- !equal(corrupt)
  }
  corrupt_hgvs <- independent
  corrupt_hgvs$hgvsp[1L] <- "p.Met1Leu"
  controls <- c(controls, hgvs_missing = !equal_hgvs(independent[-1L, ]),
    hgvs_duplicate = !equal_hgvs(rbind(independent, independent[1L, ])),
    hgvs_wrong_reference_contrast = !equal_hgvs(corrupt_hgvs))
  write.csv(data.frame(control = names(controls), rejected = controls), file.path(out, "controls.csv"), row.names = FALSE)
  perl_sources <- c(file.path(prefix, "share/ensembl-vep-116.0-0/Bio/EnsEMBL",
      c("Transcript.pm", "Translation.pm", "SeqEdit.pm")),
    file.path(prefix, "lib/perl5/site_perl/Bio", c("Tools/CodonTable.pm", "PrimarySeqI.pm")),
    file.path(mirrors[["variation"]], "modules/Bio/EnsEMBL/Variation",
      c("TranscriptHaplotypeContainer.pm", "TranscriptVariation.pm", "TranscriptVariationAllele.pm",
        "VariationFeatureOverlap.pm", "Utils/VariationEffect.pm")))
  identities <- c(code, perl_sources, shared,
    file.path(out, c("cases.jsonl", "oracle.jsonl", "pairs.csv", "witnesses.csv", "summary.csv", "controls.csv",
      "independent_hgvs.csv", "environment.txt")))
  jsonlite::write_json(list(source_revision = revision, source_binding = binding,
    oracle_revisions = revisions, cases = length(cases),
    independent_hgvs_observations = nrow(independent), independent_hgvs_equal = equal_hgvs(independent),
    scope = "native_reference_and_alternate_proteins_vs_actual_Ensembl_not_public_protein_differences",
    sha256 = as.list(vapply(identities, duckvep_evidence_sha256, ""))), file.path(out, "receipt.json"),
    pretty = TRUE, auto_unbox = TRUE)
  print(summary)
  print(pairs[pairs$family == "witness", c("id", "expected_reference", "actual_reference", "expected_alternate", "actual_alternate")])
  stopifnot(all(controls), identical(code_hashes, vapply(code, duckvep_evidence_sha256, "")),
    identical(revision, duckvep_evidence_revision(root)))
  if (nzchar(opt$evidence_out)) {
    destination <- opt$evidence_out
    stopifnot(!dir.exists(destination), dir.create(destination))
    con <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    DBI::dbWriteTable(con, "comparison_pairs", pairs)
    DBI::dbExecute(con, paste("COPY comparison_pairs TO",
      DBI::dbQuoteString(con, file.path(destination, "pairs.parquet")), "(FORMAT PARQUET)"))
    for (name in c("cases.jsonl", "oracle.jsonl")) {
      compressed <- gzfile(file.path(destination, paste0(name, ".gz")), "wt")
      writeLines(readLines(file.path(out, name)), compressed)
      close(compressed)
    }
    stopifnot(all(file.copy(file.path(out, c("summary.csv", "controls.csv",
      "independent_hgvs.csv", "environment.txt")), destination)))
    manifest <- jsonlite::read_json(file.path(out, "receipt.json"), simplifyVector = TRUE)
    manifest$local_receipt_sha256 <- duckvep_evidence_sha256(file.path(out, "receipt.json"))
    manifest$probe_sha256 <- manifest$sha256[[shared]]
    manifest$source_sha256 <- as.list(vapply(c(code, perl_sources), duckvep_evidence_sha256, ""))
    published <- list.files(destination, full.names = TRUE)
    manifest$sha256 <- as.list(setNames(vapply(published, duckvep_evidence_sha256, ""),
      basename(published)))
    jsonlite::write_json(manifest, file.path(destination, "receipt.json"),
      pretty = TRUE, auto_unbox = TRUE)
  }
  # Every sequence comparison and the raw translation metadata must agree.
  stopifnot(!any(pairs$reference_failure | pairs$alternate_full_failure | pairs$alternate_failure |
    pairs$coding_failure),
    equal_hgvs(independent))
}
main()
