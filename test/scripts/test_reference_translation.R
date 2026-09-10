#!/usr/bin/env Rscript
# Network-free reconstruction of the retained VEP-116 translation comparisons.
source("scripts/duckvep_evidence.R")

reference_read_json <- function(text) {
  check_names <- function(value) {
    if (!is.list(value)) return(invisible(TRUE))
    stopifnot(!anyDuplicated(names(value)))
    for (child in value) check_names(child)
    invisible(TRUE)
  }
  value <- jsonlite::fromJSON(text, simplifyVector = FALSE)
  check_names(value)
  value
}

reference_same_json <- function(actual, expected) {
  if (!is.list(expected)) return(identical(actual, expected))
  if (!is.list(actual) || length(actual) != length(expected)) return(FALSE)
  if (!is.null(names(expected))) {
    if (anyDuplicated(names(actual)) || !setequal(names(actual), names(expected))) return(FALSE)
    actual <- actual[names(expected)]
  } else if (!is.null(names(actual))) return(FALSE)
  all(vapply(seq_along(expected), function(i)
    reference_same_json(actual[[i]], expected[[i]]), FALSE))
}

reference_expected_cases <- function() {
  cases <- vector("list", 27014L)
  index <- 0L
  bases <- c("A", "C", "G", "T", "N")
  tables <- c(1:6, 9:14, 16L, 21:31)
  for (tail in c("", "A", "AA")) for (role in c("start", "internal", "terminal")) {
    for (third in bases) for (second in bases) for (first in bases) for (table in tables) {
      index <- index + 1L
      codon <- paste0(first, second, third)
      cds <- switch(role, start = paste0(codon, "GCCTAA"),
        internal = paste0("ATG", codon, "GCCTAA"), terminal = paste0("ATGGCC", codon))
      cases[[index]] <- list(id = paste(table, codon, role, nchar(tail), sep = "/"),
        family = role, cds = paste0(cds, tail), table = table, edits = list())
    }
  }
  stopifnot(index == 27000L)
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
  variant_positions <- c(legitimate_start = 6L, selenocysteine = 9L,
    terminal_readthrough = 6L, mito_terminal_tga = 6L, lowercase_start = 6L)
  for (name in names(witnesses)) {
    index <- index + 1L
    witness <- witnesses[[name]]
    cases[[index]] <- list(id = paste0("witness/", name), family = "witness",
      cds = witness$cds, table = if (is.null(witness$table)) 1L else witness$table,
      edits = if (is.null(witness$edits)) list() else witness$edits)
    if (name %in% names(variant_positions)) cases[[index]]$variants <- list(list(
      id = "synonymous", position1 = unname(variant_positions[name]), reference = "C", alternate = "T"))
  }
  stopifnot(index == length(cases), !anyDuplicated(vapply(cases, `[[`, "", "id")))
  cases
}

reference_check_cases <- function(cases, expected) {
  stopifnot(length(cases) == 27014L, reference_same_json(cases, expected))
  invisible(TRUE)
}

reference_expected_hgvs <- function() {
  data.frame(case_id = paste0("witness/", c("legitimate_start", "selenocysteine",
    "terminal_readthrough", "mito_terminal_tga", "lowercase_start", "lowercase_start")),
    allele = c("T", "T", "T", "T", "C", "T"),
    hgvsp = c("witness/legitimate_start_protein.1:p.Ala2=",
      "witness/selenocysteine_protein.1:p.Ala3=",
      "witness/terminal_readthrough_protein.1:p.Ala2=",
      "witness/mito_terminal_tga_protein.1:p.Ala2=",
      "witness/lowercase_start_protein.1:p.Ala2=",
      "witness/lowercase_start_protein.1:p.Ala2="))
}

reference_oracle_expectations <- function(cases, oracle) {
  stopifnot(length(oracle) == length(cases),
    identical(vapply(oracle, `[[`, "", "id"), vapply(cases, `[[`, "", "id")))
  hgvs <- list()
  rows <- lapply(seq_along(cases), function(i) {
    case <- cases[[i]]
    observed <- oracle[[i]]
    fields <- c("id", "prepared_cds", "core_reference", "reference", "alternate_full", "alternate")
    if (!is.null(case$variants)) fields <- c(fields, "independent_hgvs")
    stopifnot(setequal(names(observed), fields),
      all(vapply(observed[setdiff(fields, "independent_hgvs")], function(x)
        is.character(x) && length(x) == 1L && !is.na(x), FALSE)),
      identical(toupper(observed$prepared_cds), toupper(case$cds)))
    if (!is.null(case$variants)) {
      expected_alleles <- if (case$id == "witness/lowercase_start") c("C", "T") else "T"
      stopifnot(length(observed$independent_hgvs) == length(expected_alleles))
      for (j in seq_along(expected_alleles)) {
        row <- observed$independent_hgvs[[j]]
        stopifnot(setequal(names(row), c("id", "allele", "consequences", "hgvsp")),
          identical(row$id, "synonymous"), identical(row$allele, expected_alleles[j]),
          identical(row$consequences, list("synonymous_variant")))
        hgvs[[length(hgvs) + 1L]] <<- data.frame(case_id = case$id, allele = row$allele, hgvsp = row$hgvsp)
      }
    }
    coding_length <- nchar(observed$alternate_full)
    stopifnot(coding_length == nchar(observed$prepared_cds) %/% 3L)
    stop_position <- regexpr("*", observed$alternate_full, fixed = TRUE)[1L]
    data.frame(id = case$id, family = case$family, table = case$table, cds = case$cds,
      prepared_cds = observed$prepared_cds, core_reference = observed$core_reference,
      expected_reference = observed$reference, expected_coding = observed$alternate_full,
      expected_coding_length = coding_length, expected_coding_stop = max(0L, stop_position),
      expected_coding_unambiguous = !grepl("N", toupper(observed$prepared_cds), fixed = TRUE),
      expected_alternate_full = observed$alternate_full, expected_alternate = observed$alternate)
  })
  observed_hgvs <- do.call(rbind, hgvs)
  stopifnot(identical(observed_hgvs, reference_expected_hgvs()))
  list(pairs = do.call(rbind, rows), hgvs = observed_hgvs)
}

reference_check_pairs <- function(pairs, expected) {
  actual_fields <- c("raw_reference", "actual_reference", "actual_coding", "actual_coding_length",
    "actual_coding_stop", "actual_coding_unambiguous", "actual_alternate_full", "actual_alternate")
  flag_fields <- c("coding_failure", "reference_failure", "raw_reference_hypothesis_failure",
    "alternate_full_failure", "alternate_failure")
  stopifnot(nrow(pairs) == nrow(expected), !anyDuplicated(names(pairs)), !anyDuplicated(pairs$id),
    setequal(names(pairs), c(names(expected), actual_fields, flag_fields)), !anyNA(pairs),
    identical(pairs$id, expected$id))
  for (field in names(expected)) stopifnot(isTRUE(all.equal(pairs[[field]], expected[[field]], tolerance = 0)))
  for (field in c(flag_fields, "actual_coding_unambiguous")) stopifnot(is.logical(pairs[[field]]))
  stopifnot(identical(pairs$raw_reference, pairs$actual_alternate_full))
  checked <- pairs
  checked$coding_failure <- pairs$actual_coding != expected$expected_coding |
    pairs$actual_coding_length != expected$expected_coding_length |
    pairs$actual_coding_stop != expected$expected_coding_stop |
    pairs$actual_coding_unambiguous != expected$expected_coding_unambiguous
  checked$reference_failure <- pairs$actual_reference != expected$expected_reference
  checked$raw_reference_hypothesis_failure <- pairs$raw_reference != expected$expected_reference
  checked$alternate_full_failure <- pairs$actual_alternate_full != expected$expected_alternate_full
  checked$alternate_failure <- pairs$actual_alternate != expected$expected_alternate
  for (field in flag_fields) stopifnot(identical(pairs[[field]], checked[[field]]))
  # This retained run is a passing translation comparison. The rejected raw-reference
  # hypothesis remains a separate measured result, never a reason to drop a pair.
  stopifnot(!any(checked$coding_failure | checked$reference_failure |
    checked$alternate_full_failure | checked$alternate_failure))
  aggregate(cbind(cases = rep(1L, nrow(checked)), reference_failures = checked$reference_failure,
    coding_failures = checked$coding_failure,
    raw_reference_hypothesis_failures = checked$raw_reference_hypothesis_failure,
    alternate_full_failures = checked$alternate_full_failure, alternate_failures = checked$alternate_failure),
    checked["family"], sum)
}

reference_check_summary <- function(actual, expected) {
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = 0)))
  invisible(TRUE)
}

reference_check_receipt <- function(receipt) {
  stopifnot(all(c("source_revision", "source_binding", "oracle_revisions", "cases",
    "independent_hgvs_observations", "independent_hgvs_equal", "scope", "sha256",
    "source_sha256", "local_receipt_sha256", "probe_sha256") %in% names(receipt)),
    is.character(receipt$source_revision), length(receipt$source_revision) == 1L,
    grepl("^[0-9a-f]{40}$", receipt$source_revision),
    identical(receipt$source_binding, "diagnostic_unbound"), identical(receipt$cases, 27014L),
    identical(receipt$independent_hgvs_observations, 6L), isTRUE(receipt$independent_hgvs_equal),
    identical(unlist(receipt$oracle_revisions, use.names = FALSE),
      c("57ea5c52340acc1f156267f810ad162e26597082", "2fb834b987ede3824e200197a838ce11e91aeb4b")),
    identical(receipt$scope,
      "native_reference_and_alternate_proteins_vs_actual_Ensembl_not_public_protein_differences"),
    is.character(receipt$local_receipt_sha256), length(receipt$local_receipt_sha256) == 1L,
    grepl("^[0-9a-f]{64}$", receipt$local_receipt_sha256),
    is.character(receipt$probe_sha256), length(receipt$probe_sha256) == 1L,
    grepl("^[0-9a-f]{64}$", receipt$probe_sha256))
  for (field in c("sha256", "source_sha256")) {
    hashes <- unlist(receipt[[field]])
    stopifnot(is.character(hashes), length(hashes) > 0L, !is.null(names(hashes)),
      !anyDuplicated(names(hashes)), !anyNA(hashes), all(grepl("^[0-9a-f]{64}$", hashes)))
  }
  stopifnot(setequal(names(receipt$sha256), c("cases.jsonl.gz", "oracle.jsonl.gz", "pairs.parquet",
    "summary.csv", "controls.csv", "independent_hgvs.csv", "environment.txt")))
  invisible(TRUE)
}

reference_rejects <- function(expression) {
  tryCatch({ force(expression); FALSE }, error = function(e) TRUE)
}

reference_corruption_controls <- function(cases, oracle, expected, pairs, receipt) {
  controls <- c(missing = reference_rejects(reference_check_pairs(pairs[-1L, ], expected)),
    duplicate = reference_rejects(reference_check_pairs(rbind(pairs, pairs[1L, ]), expected)))
  retained_fields <- c("id", "expected_reference", "expected_alternate_full", "expected_alternate",
    "expected_coding", "expected_coding_length", "expected_coding_stop", "expected_coding_unambiguous")
  for (field in retained_fields) {
    changed <- pairs
    value <- changed[[field]][1L]
    changed[[field]][1L] <- if (is.logical(value)) !value else if (is.numeric(value))
      value + 1 else paste0(value, "X")
    controls[field] <- reference_rejects(reference_check_pairs(changed, expected))
  }
  hgvs <- reference_expected_hgvs()
  check_hgvs <- function(x) stopifnot(identical(x, hgvs))
  changed <- hgvs
  changed$hgvsp[1L] <- "p.Met1Leu"
  controls <- c(controls, hgvs_missing = reference_rejects(check_hgvs(hgvs[-1L, ])),
    hgvs_duplicate = reference_rejects(check_hgvs(rbind(hgvs, hgvs[1L, ]))),
    hgvs_wrong_reference_contrast = reference_rejects(check_hgvs(changed)))
  retained <- data.frame(control = names(controls), rejected = unname(controls))
  for (field in c("actual_reference", "actual_coding", "actual_alternate_full", "actual_alternate",
      "actual_coding_length", "actual_coding_stop", "actual_coding_unambiguous", "coding_failure")) {
    changed <- pairs
    value <- changed[[field]][1L]
    changed[[field]][1L] <- if (is.logical(value)) !value else if (is.numeric(value))
      value + 1 else paste0(value, "X")
    controls[paste0("native_", field)] <- reference_rejects(reference_check_pairs(changed, expected))
  }
  changed <- pairs
  changed$expected_coding_length[1L] <- changed$expected_coding_length[1L] + 1e-10
  controls["fractional_expected_metadata"] <- reference_rejects(reference_check_pairs(changed, expected))
  summary <- reference_check_pairs(pairs, expected)
  changed <- summary
  changed$cases[1L] <- changed$cases[1L] + 1e-7
  controls["fractional_summary_metadata"] <- reference_rejects(reference_check_summary(changed, summary))
  declared <- reference_expected_cases()
  controls["missing_case"] <- reference_rejects(reference_check_cases(cases[-1L], declared))
  controls["duplicate_case"] <- reference_rejects(reference_check_cases(c(cases, cases[1L]), declared))
  for (field in c("id", "family", "cds", "table", "edits")) {
    changed <- cases
    changed[[1L]][[field]] <- switch(field, id = "wrong", family = "terminal", cds = "AAAGCCTAAA",
      table = 2L, edits = list(list(code = "initial_met", position1 = 1L, alternate = "M")))
    controls[paste0("case_", field)] <- reference_rejects(reference_check_cases(changed, declared))
  }
  changed <- cases
  changed[[27004L]]$edits[[1L]]$alternate <- "V"
  controls["witness_edit"] <- reference_rejects(reference_check_cases(changed, declared))
  changed <- cases
  changed[[27001L]]$variants[[1L]]$position1 <- 7L
  controls["witness_variant"] <- reference_rejects(reference_check_cases(changed, declared))
  controls["missing_oracle"] <- reference_rejects(reference_oracle_expectations(cases, oracle[-1L]))
  controls["duplicate_oracle"] <- reference_rejects(reference_oracle_expectations(cases, c(oracle, oracle[1L])))
  changed <- oracle
  changed[[1L]]$prepared_cds <- "AAAGCCTAAA"
  controls["oracle_source"] <- reference_rejects(reference_oracle_expectations(cases, changed))
  for (field in c("cases", "source_binding", "scope", "oracle_revisions", "local_receipt_sha256")) {
    changed <- receipt
    changed[[field]] <- switch(field, cases = 27013L, source_binding = "certified", scope = "public_hgvsp",
      oracle_revisions = list("unversioned"), local_receipt_sha256 = "missing")
    controls[paste0("receipt_", field)] <- reference_rejects(reference_check_receipt(changed))
  }
  changed <- receipt
  changed$local_receipt_sha256 <- NULL
  controls["receipt_missing_hash"] <- reference_rejects(reference_check_receipt(changed))
  duplicate_json <- c('{"id":"a","id":"b"}',
    '{"edits":[{"position1":1,"position1":2}]}',
    '{"variants":[{"reference":"A","reference":"N"}]}',
    '{"independent_hgvs":[{"hgvsp":"p.Ala2=","hgvsp":null}]}',
    '{"sha256":{"cases.jsonl.gz":"a","cases.jsonl.gz":"b"}}')
  for (i in seq_along(duplicate_json)) controls[paste0("json_duplicate_", i)] <-
    reference_rejects(reference_read_json(duplicate_json[i]))
  stopifnot(!anyDuplicated(names(controls)), all(controls))
  list(retained = retained, count = length(controls))
}

main <- function() {
  directory <- "test/duckvep/conformance/data/reference_translation_consensus"
  receipt <- reference_read_json(paste(readLines(file.path(directory, "receipt.json")), collapse = "\n"))
  reference_check_receipt(receipt)
  stopifnot(setequal(list.files(directory), c(names(receipt$sha256), "receipt.json")))
  for (name in names(receipt$sha256)) stopifnot(identical(receipt$sha256[[name]],
    duckvep_evidence_sha256(file.path(directory, name))))
  stopifnot(identical(duckvep_evidence_explicit_packages(readLines(file.path(directory, "environment.txt"))),
    duckvep_evidence_explicit_packages(readLines(
      "test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt"))))
  read_records <- function(name) {
    connection <- gzfile(file.path(directory, name), "rt")
    on.exit(close(connection))
    lapply(readLines(connection), reference_read_json)
  }
  cases <- read_records("cases.jsonl.gz")
  reference_check_cases(cases, reference_expected_cases())
  oracle <- read_records("oracle.jsonl.gz")
  expected <- reference_oracle_expectations(cases, oracle)
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  pairs <- DBI::dbGetQuery(con, paste("SELECT * FROM read_parquet(",
    DBI::dbQuoteString(con, file.path(directory, "pairs.parquet")), ")"))
  summary <- reference_check_pairs(pairs, expected$pairs)
  reference_check_summary(summary, read.csv(file.path(directory, "summary.csv")))
  stopifnot(identical(expected$hgvs, read.csv(file.path(directory, "independent_hgvs.csv"))))
  controls <- reference_corruption_controls(cases, oracle, expected$pairs, pairs, receipt)
  stopifnot(identical(controls$retained, read.csv(file.path(directory, "controls.csv"))))
  message("Reference translation: all 27,014 comparisons reconstructed; ", controls$count,
    " corruption controls rejected. Six HGVS observations checked against literals; no native HGVSp claim.")
}

if (sys.nframe() == 0L) main()
