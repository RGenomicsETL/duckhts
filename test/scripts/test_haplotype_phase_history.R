#!/usr/bin/env Rscript
# Synthetic receipts exercise publication integrity without an oracle installation.
source("test/duckvep/conformance/haplotype_phase_differential.R")
source("scripts/duckvep_evidence.R")

main <- function() {
  directory <- tempfile("phase-history-test-", "test/duckvep/conformance/results")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  plain <- genotypes(2L)
  gt <- c(plain, paste0("|", plain), paste0("/", plain))
  n <- length(gt)
  ploidy <- lengths(strsplit(sub("^[|/]", "", gt), "[|/]"))
  groups <- canonical(list(list(cds = "ATG", protein = "M", count = 2, contributors = "a")))
  comparison <- list(expected = groups, observed = groups, equal = TRUE,
    oracle_lanes = 2L, native_lanes = 2, native_unknown = 0L, native_unavailable_carriers = 0)
  decoded <- rep(list(comparison), n)
  decoded[[1L]]$observed[[1L]]$cds <- "CTG"
  decoded[[1L]]$observed <- canonical(decoded[[1L]]$observed)
  decoded[[1L]]$equal <- FALSE
  summary <- data.frame(transcript_index = seq_len(n) - 1L, seq_region = seq_len(n) - 1L,
    chrom = sprintf("chrP%05d", seq_len(n)), transcript = sprintf("HP%05d", seq_len(n)), GT = gt, ploidy,
    prefix = grepl("^[|/]", gt), missing = grepl(".", gt, fixed = TRUE),
    mixed = grepl("|", gt, fixed = TRUE) & grepl("/", gt, fixed = TRUE),
    equal = seq_len(n) != 1L, oracle_lanes = 2, native_lanes = 2,
    native_unknown = 0, native_unavailable_carriers = 0)
  parser_rows <- data.frame(status = rep(0L, 2L * n), retained = TRUE, ploidy = 2L,
    missing = 0L, slots = 2L, first = 0L, second = 0L)
  semantics <- data.frame(source_indices = rep(0L, 4L * n), source_evidence = 0L,
    evidence = 0L, sequence_status = 0L)
  raw <- list(output = list(errors = integer(n)), comparisons = rep(list(comparison), n),
    expected_semantics = semantics, observed_semantics = semantics)
  native <- list(records = data.frame(record_index = seq_len(2L * n) - 1L, ID = rep(c("a", "b"), n)),
    calls = data.frame(event_index = seq_len(3L * n)),
    actual = data.frame(transcript_index = rep(seq_len(n) - 1L, each = 2L),
      carrier_count = 1, cds = "ATG", protein = "M"))
  native$actual$contributors <- lapply(native$actual$transcript_index,
    function(i) data.frame(event_index = 6L * i + 1L))
  raw$actual <- native$actual
  native$actual$cds[native$actual$transcript_index == 0L] <- "CTG"
  native$records$ALT <- rep(list(c("A", "T"), "C"), n)
  native$records$calls <- lapply(rep(ploidy, each = 2L), function(p) {
    x <- data.frame(sample_index = 0L)
    x$alleles <- list(integer(p))
    x
  })
  saveRDS(decoded, file.path(directory, "comparisons.rds"))
  parser <- list(keys = data.frame(transcript = rep(summary$transcript, each = 2L),
    source_id = rep(c("a", "b"), n), GT = as.vector(rbind(gt, vapply(ploidy,
      function(p) paste(rep("1", p), collapse = "|"), "")))),
    expected = parser_rows, observed = parser_rows)
  saveRDS(parser,
    file.path(directory, "phase_comparisons.rds"))
  saveRDS(raw, file.path(directory, "raw_replay.rds"))
  saveRDS(raw, file.path(directory, "public_raw_replay.rds"))
  saveRDS(native, file.path(directory, "native.rds"))
  saveRDS(list(), file.path(directory, "decoded_collisions.rds"))
  writeLines("synthetic source input", file.path(directory, "calls.vcf"))
  write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
  cases <- summary[c("transcript_index", "seq_region", "chrom", "transcript", "GT",
    "ploidy", "prefix", "missing", "mixed")]
  write.table(transform(cases, cds = "ATG"), file.path(directory, "cases.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  oracle <- lapply(cases$transcript, function(transcript) list(transcript = transcript,
    haplotypes = unname(groups), total_haplotype_count = 2))
  writeLines(vapply(oracle, jsonlite::toJSON, "", auto_unbox = TRUE),
    file.path(directory, "oracle.stdout"))
  for (kind in c("raw", "public_raw"))
    write.csv(transform(cases, equal = TRUE), file.path(directory, paste0(kind, "_replay_summary.csv")),
      row.names = FALSE)
  for (file in c("controls.csv", "phase_controls.csv", "raw_replay_controls.csv")) {
    control <- switch(file, "controls.csv" = c("duplicate", "cds", "protein", "contributor"),
      "phase_controls.csv" = names(parser_rows), "raw_replay_controls.csv" = names(semantics))
    write.csv(data.frame(control, rejected = TRUE), file.path(directory, file), row.names = FALSE)
  }
  receipt <- list(source_revision = strrep("a", 40L),
    extension_build_binding = "htslib_distclean_make_release",
    scope = "raw_GT_finite_phase_audit_not_conformance", max_ploidy = 2L,
    oracle_revisions = list(vep = strrep("b", 40L), variation = strrep("c", 40L)),
    cases = n, disagreements = 1L, oracle_lanes = 2L * n, native_lanes = 2L * n,
    native_unavailable_carriers = 0L, decoded_collision_groups = 0L, controls_rejected = 4L,
    raw_parser_calls = 2L * n, raw_parser_disagreements = 0L, raw_parser_controls_rejected = 7L,
    raw_replay_cases = n, raw_replay_disagreements = 0L, raw_replay_record_observations = 4L * n,
    raw_replay_record_disagreements = 0L, raw_replay_controls_rejected = 4L,
    public_raw_replay_cases = n, public_raw_replay_disagreements = 0L,
    public_raw_record_observations = 4L * n, public_raw_record_disagreements = 0L,
    input_records = 2L * n, source_alt_events = 3L * n, input_genotype_calls = 2L * n,
    input_allele_slots = 2L * sum(ploidy), candidate_alt_calls = 3L * n, native_leaves = 2L * n)
  update_receipt <- function(value = receipt) {
    files <- list.files(directory, full.names = TRUE)
    files <- files[basename(files) != "receipt.json"]
    value$sha256 <- as.list(c(vapply(files, duckvep_evidence_sha256, ""),
      "/synthetic/duckhts.duckdb_extension" = strrep("d", 64L)))
    jsonlite::write_json(value, file.path(directory, "receipt.json"), auto_unbox = TRUE)
  }
  fails <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), "error"))
  update_receipt()
  rows <- phase_history_rows(directory)
  stopifnot(sum(rows$cases) == 108L, sum(rows$disagreements) == 1L,
    sum(rows$raw_parser_calls) == 216L, sum(rows$public_raw_record_observations) == 432L)
  unlink(file.path(directory, "calls.vcf"))
  fails(phase_history_rows(directory))
  writeLines("synthetic source input", file.path(directory, "calls.vcf"))
  outside <- tempfile("phase-history-outside-")
  dir.create(outside)
  on.exit(unlink(outside, recursive = TRUE), add = TRUE)
  fails(phase_history_rows(outside))
  for (field in c("cases", "disagreements", "input_records", "input_allele_slots",
    "controls_rejected", "raw_parser_disagreements", "public_raw_record_disagreements")) {
    changed <- receipt
    changed[[field]] <- changed[[field]] + 1L
    update_receipt(changed)
    fails(phase_history_rows(directory))
  }
  update_receipt(modifyList(receipt, list(extension_build_binding = "diagnostic_unbound")))
  fails(phase_history_rows(directory))
  update_receipt()
  changed <- summary
  changed$equal[1L] <- TRUE
  write.csv(changed, file.path(directory, "summary.csv"), row.names = FALSE)
  fails(phase_history_rows(directory))  # Unreceipted byte change.
  update_receipt()
  fails(phase_history_rows(directory))  # Rehashed summary still contradicts retained comparisons.
  changed$GT[1L] <- changed$GT[2L]
  write.csv(changed, file.path(directory, "summary.csv"), row.names = FALSE)
  update_receipt()
  fails(phase_history_rows(directory))
  write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
  swapped <- decoded
  swapped[c(1L, n)] <- swapped[c(n, 1L)]
  changed <- summary
  changed$equal[c(1L, n)] <- changed$equal[c(n, 1L)]
  saveRDS(swapped, file.path(directory, "comparisons.rds"))
  write.csv(changed, file.path(directory, "summary.csv"), row.names = FALSE)
  update_receipt()
  fails(phase_history_rows(directory))
  saveRDS(decoded, file.path(directory, "comparisons.rds"))
  write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
  for (change in c("missing_payload", "null_payload", "empty_payload", "missing_provenance", "wrong_count")) {
    changed_raw <- raw
    if (change == "missing_payload") changed_raw$comparisons <- rep(list(list(equal = TRUE)), n)
    if (change == "null_payload") changed_raw$comparisons[[1L]][c("expected", "observed")] <- list(NULL, NULL)
    if (change == "empty_payload") changed_raw$comparisons[[1L]][c("expected", "observed")] <- list(list(), list())
    if (change == "missing_provenance") changed_raw$comparisons[[1L]]$observed[[1L]]$contributors <- NULL
    if (change == "wrong_count") {
      changed_raw$comparisons[[1L]]$expected[[1L]]$count <- 1
      changed_raw$comparisons[[1L]]$observed[[1L]]$count <- 1
    }
    saveRDS(changed_raw, file.path(directory, "raw_replay.rds"))
    update_receipt()
    fails(phase_history_rows(directory))
  }
  saveRDS(raw, file.path(directory, "raw_replay.rds"))
  changed_parser <- parser
  changed_parser$expected <- changed_parser$expected["status"]
  changed_parser$observed <- changed_parser$observed["status"]
  saveRDS(changed_parser, file.path(directory, "phase_comparisons.rds"))
  update_receipt()
  fails(phase_history_rows(directory))
  changed_parser <- parser
  changed_parser$keys$GT[c(1L, 4L)] <- changed_parser$keys$GT[c(4L, 1L)]
  saveRDS(changed_parser, file.path(directory, "phase_comparisons.rds"))
  update_receipt()
  fails(phase_history_rows(directory))
  saveRDS(parser, file.path(directory, "phase_comparisons.rds"))
  # A genuinely absent observed lane remains a disagreement, not an excluded case.
  absent <- decoded
  absent[[1L]]$observed <- list()
  absent[[1L]]$native_lanes <- 0
  absent_native <- native
  absent_native$actual <- native$actual[native$actual$transcript_index != 0L, ]
  absent_summary <- summary
  absent_summary$native_lanes[1L] <- 0
  saveRDS(absent, file.path(directory, "comparisons.rds"))
  saveRDS(absent_native, file.path(directory, "native.rds"))
  write.csv(absent_summary, file.path(directory, "summary.csv"), row.names = FALSE)
  update_receipt(modifyList(receipt, list(native_lanes = 2L * n - 2L, native_leaves = 2L * n - 2L)))
  absent_rows <- phase_history_rows(directory)
  stopifnot(sum(absent_rows$cases) == n, sum(absent_rows$disagreements) == 1L,
    sum(absent_rows$native_lanes) == 2L * n - 2L)
  saveRDS(decoded, file.path(directory, "comparisons.rds"))
  saveRDS(native, file.path(directory, "native.rds"))
  write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
  saveRDS(raw[-match("observed_semantics", names(raw))], file.path(directory, "raw_replay.rds"))
  update_receipt()
  fails(phase_history_rows(directory))
  saveRDS(raw, file.path(directory, "raw_replay.rds"))
  write.csv(data.frame(control = c("duplicate", "cds", "protein", "protein"), rejected = TRUE),
    file.path(directory, "controls.csv"), row.names = FALSE)
  update_receipt()
  fails(phase_history_rows(directory))
  write.csv(data.frame(control = c("duplicate", "cds", "protein", "contributor"), rejected = TRUE),
    file.path(directory, "controls.csv"), row.names = FALSE)
  update_receipt()
  history <- tempfile("phase-history-csv-", "test/duckvep/conformance/results", fileext = ".csv")
  on.exit(unlink(history), add = TRUE)
  write.csv(read.csv("test/duckvep/conformance/data/haplotype_phase_history.csv")[FALSE, ],
    history, row.names = FALSE)
  stopifnot(dir.create(paste0(history, ".lock")))
  before <- readBin(history, "raw", file.info(history)$size)
  fails(publish_phase_history(directory, history))
  stopifnot(identical(before, readBin(history, "raw", file.info(history)$size)),
    dir.exists(paste0(history, ".lock")))
  unlink(paste0(history, ".lock"), recursive = TRUE)
  publish_phase_history(directory, history)
  before <- readBin(history, "raw", file.info(history)$size)
  fails(publish_phase_history(directory, history))
  stopifnot(identical(before, readBin(history, "raw", file.info(history)$size)))
  stopifnot(!dir.exists(paste0(history, ".lock")))
  message("Raw-GT history publication: receipt, denominator, failure and duplicate controls pass")
}
main()
