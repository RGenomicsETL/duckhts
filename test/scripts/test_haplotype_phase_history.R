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
  parser_rows <- data.frame(status = rep(0L, 2L * n), retained = TRUE, ploidy = rep(ploidy, each = 2L),
    missing = as.vector(rbind(as.integer(summary$missing), 0L)), slots = 2L,
    first = rep(c(1L, 0L), n), second = rep(c(1L, 0L), n))
  lane_evidence <- rep(ifelse(summary$missing, 11L, 1L), each = 2L)
  lane_status <- rep(ifelse(summary$missing, 8L, 0L), each = 2L)
  semantics <- data.frame(source_indices = rep(c(1L, -2L), 2L * n),
    source_evidence = as.vector(rbind(lane_evidence, 0L)),
    evidence = rep(lane_evidence, each = 2L), sequence_status = rep(lane_status, each = 2L))
  native <- list(records = data.frame(record_index = seq_len(2L * n) - 1L,
      CHROM = rep(summary$chrom, each = 2L), POS = rep(c(41, 44), n),
      ID = rep(c("a", "b"), n), REF = "G"),
    calls = data.frame(event_index = seq_len(3L * n)),
    actual = data.frame(transcript_index = rep(seq_len(n) - 1L, each = 2L),
      carrier_count = 1, cds = "ATG", protein = "M"))
  native$actual$contributors <- lapply(native$actual$transcript_index,
    function(i) data.frame(event_index = 6L * i + 1L))
  public <- list(actual = native$actual,
    comparisons = rep(list(comparison[c("expected", "observed", "equal")]), n),
    expected_semantics = semantics, observed_semantics = semantics)
  public$actual$coding_blocks <- lapply(public$actual$transcript_index,
    function(i) data.frame(event_indices = I(list(2L * i))))
  public$actual$contributors <- lapply(seq_len(nrow(public$actual)), function(row) {
    i <- public$actual$transcript_index[row]
    data.frame(event_index = 2L * i, seq_region = i, position = 41,
      reference = "G", alt_index = 1L, alternate = "A", evidence_flags = lane_evidence[row])
  })
  public$actual$carriers <- lapply(seq_len(nrow(public$actual)), function(row)
    data.frame(sample_index = 0L, ploidy = 2L, phase_set = NA_integer_,
      haplotype_lane = (row - 1L) %% 2L + 1L))
  public$actual$evidence_flags <- lane_evidence
  public$actual$sequence_status <- ifelse(lane_status == 0L, "ok", "conditional")
  native$actual$cds[native$actual$transcript_index == 0L] <- "CTG"
  native$records$ALT <- rep(list(c("A", "T"), "C"), n)
  raw_gt <- as.vector(rbind(gt, vapply(ploidy,
    function(p) paste(rep("1", p), collapse = "|"), "")))
  capacity <- 512L
  replay <- c(list("ATG", 11L, raw_gt, c(41L, 44L), c("G", "G"), c("A", "T", "C"),
      c(2L, 1L), n, capacity),
    list(cds = rep(c(charToRaw("ATG"), raw(capacity - 3L)), 2L * n),
      protein = rep(c(charToRaw("M"), raw(capacity - 1L)), 2L * n),
      cds_lengths = rep(3L, 2L * n), protein_lengths = rep(1L, 2L * n),
      sequence_status = lane_status, evidence = lane_evidence, edit_masks = rep(1L, 2L * n),
      source_indices = semantics$source_indices, source_evidence = semantics$source_evidence,
      errors = integer(n)))
  raw <- list(output = replay, comparisons = public$comparisons,
    expected_semantics = semantics, observed_semantics = semantics)
  native$records$calls <- lapply(seq_along(raw_gt), function(i) {
    x <- data.frame(sample_index = 0L, raw_gt = raw_gt[i])
    x$alleles <- list(integer(ploidy[(i + 1L) %/% 2L]))
    x
  })
  saveRDS(decoded, file.path(directory, "comparisons.rds"))
  parser <- list(keys = data.frame(transcript = rep(summary$transcript, each = 2L),
    source_id = rep(c("a", "b"), n), GT = as.vector(rbind(gt, vapply(ploidy,
      function(p) paste(rep("1", p), collapse = "|"), "")))),
    expected = parser_rows, observed = parser_rows, disposition = rep(3L, 2L * n))
  saveRDS(parser,
    file.path(directory, "phase_comparisons.rds"))
  saveRDS(raw, file.path(directory, "raw_replay.rds"))
  saveRDS(public, file.path(directory, "public_raw_replay.rds"))
  saveRDS(native, file.path(directory, "native.rds"))
  saveRDS(list(), file.path(directory, "decoded_collisions.rds"))
  vcf <- c("##fileformat=VCFv4.4",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
    paste(native$records$CHROM, native$records$POS, native$records$ID, "G",
      rep(c("A,T", "C"), n), ".", "PASS", ".", "GT:PS",
      paste0(raw_gt, rep(c(":10", ":20"), n)), sep = "\t"))
  writeLines(vcf, file.path(directory, "calls.vcf"))
  write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
  cases <- summary[c("transcript_index", "seq_region", "chrom", "transcript", "GT",
    "ploidy", "prefix", "missing", "mixed")]
  write.table(transform(cases, cds = "ATG"), file.path(directory, "cases.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  oracle <- lapply(cases$transcript, function(transcript) list(transcript = transcript,
    haplotypes = unname(groups), total_haplotype_count = 2))
  names(oracle) <- cases$transcript
  writeLines(vapply(oracle, jsonlite::toJSON, "", auto_unbox = TRUE),
    file.path(directory, "oracle.stdout"))
  phase <- lapply(cases$transcript, function(transcript) list(transcript = transcript,
    default_ploidy = 2L, sample_ploidy = list(sample = 2L), calls = list(
      list(source_id = "a", sample = "sample", genotype = c("A", "A")),
      list(source_id = "b", sample = "sample", genotype = c("G", "G")))))
  writeLines(vapply(phase, jsonlite::toJSON, "", auto_unbox = TRUE), file.path(directory, "phase.jsonl"))
  write_replay_summaries <- function(profiles) {
    for (kind in c("raw", "public_raw"))
      write.csv(transform(profiles, equal = TRUE),
        file.path(directory, paste0(kind, "_replay_summary.csv")), row.names = FALSE)
  }
  write_replay_summaries(cases)
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
  for (field in c("chrom", "transcript", "seq_region")) {
    changed_cases <- cases
    changed_cases[c(1L, n), field] <- changed_cases[c(n, 1L), field]
    fails(phase_check_source_records(directory, changed_cases, native$records))
  }
  for (sample_index in c(0.5, NA_real_, 1)) {
    changed_records <- native$records
    changed_records$calls[[1L]]$sample_index <- sample_index
    fails(phase_check_source_records(directory, cases, changed_records))
  }
  update_receipt()
  rows <- phase_history_rows(directory)
  stopifnot(sum(rows$cases) == 108L, sum(rows$disagreements) == 1L,
    sum(rows$raw_parser_calls) == 216L, sum(rows$public_raw_record_observations) == 432L)
  unlink(file.path(directory, "calls.vcf"))
  fails(phase_history_rows(directory))
  writeLines(vcf, file.path(directory, "calls.vcf"))
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
  # Reassign the profile labels coherently while keeping the physical input and
  # transcript-keyed observations intact. Global denominators remain unchanged.
  changed_cases <- cases
  profile_fields <- c("GT", "ploidy", "prefix", "missing", "mixed")
  changed_cases[c(1L, n), profile_fields] <- changed_cases[c(n, 1L), profile_fields]
  changed <- summary
  changed[profile_fields] <- changed_cases[profile_fields]
  changed_parser <- parser
  changed_parser$keys$GT <- as.vector(rbind(changed_cases$GT,
    vapply(changed_cases$ploidy, function(p) paste(rep("1", p), collapse = "|"), "")))
  write.table(transform(changed_cases, cds = "ATG"), file.path(directory, "cases.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  write.csv(changed, file.path(directory, "summary.csv"), row.names = FALSE)
  write_replay_summaries(changed_cases)
  saveRDS(changed_parser, file.path(directory, "phase_comparisons.rds"))
  update_receipt()
  fails(phase_history_rows(directory))
  write.table(transform(cases, cds = "ATG"), file.path(directory, "cases.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
  write_replay_summaries(cases)
  saveRDS(parser, file.path(directory, "phase_comparisons.rds"))
  changed_native <- native
  changed_native$records$calls[[1L]]$raw_gt <- "1"
  saveRDS(changed_native, file.path(directory, "native.rds"))
  update_receipt()
  fails(phase_history_rows(directory))
  saveRDS(native, file.path(directory, "native.rds"))
  for (file in c("raw_replay.rds", "public_raw_replay.rds")) {
    retained <- if (file == "raw_replay.rds") raw else public
    changed_raw <- retained
    forged <- canonical(list(list(cds = "FORGED", protein = "FORGED", count = 2,
      contributors = "unobserved_source")))
    changed_raw$comparisons <- rep(list(list(expected = forged, observed = forged, equal = TRUE)), n)
    saveRDS(changed_raw, file.path(directory, file))
    update_receipt()
    fails(phase_history_rows(directory))
    saveRDS(retained, file.path(directory, file))
  }
  for (kind in c("raw", "public_raw")) {
    file <- paste0(kind, "_replay.rds")
    retained <- if (kind == "raw") raw else public
    changed <- retained
    if (kind == "raw") {
      changed$output$cds[1L] <- charToRaw("C")
      changed$comparisons <- phase_raw_comparisons(transform(cases, cds = "ATG"), oracle, changed$output)
    } else {
      changed$actual$cds[1L] <- "CTG"
      changed$comparisons <- phase_public_comparisons(cases, oracle, changed$actual, native$records)
    }
    # One genuinely different observed lane stays in its source-GT stratum.
    lane_summary <- transform(cases, equal = seq_len(n) != 1L)
    lane_receipt <- receipt
    lane_receipt[[paste0(kind, "_replay_disagreements")]] <- 1L
    saveRDS(changed, file.path(directory, file))
    write.csv(lane_summary, file.path(directory, paste0(kind, "_replay_summary.csv")), row.names = FALSE)
    update_receipt(lane_receipt)
    failed_rows <- phase_history_rows(directory)
    stopifnot(sum(failed_rows[[paste0(kind, "_replay_cases")]]) == n,
      sum(failed_rows[[paste0(kind, "_replay_disagreements")]]) == 1L)
    changed$comparisons[c(1L, n)] <- changed$comparisons[c(n, 1L)]
    lane_summary$equal[c(1L, n)] <- lane_summary$equal[c(n, 1L)]
    saveRDS(changed, file.path(directory, file))
    write.csv(lane_summary, file.path(directory, paste0(kind, "_replay_summary.csv")), row.names = FALSE)
    update_receipt(lane_receipt)
    fails(phase_history_rows(directory))
    saveRDS(retained, file.path(directory, file))
    write_replay_summaries(cases)
  }
  for (field in names(parser_rows)) {
    changed_parser <- parser
    value <- parser_rows[[field]][1L]
    value <- if (is.logical(value)) !value else value + 1L
    changed_parser$expected[[field]][1L] <- value
    changed_parser$observed[[field]][1L] <- value
    if (field == "retained") changed_parser$disposition[1L] <- 1L
    saveRDS(changed_parser, file.path(directory, "phase_comparisons.rds"))
    update_receipt()
    fails(phase_history_rows(directory))
  }
  saveRDS(parser, file.path(directory, "phase_comparisons.rds"))
  for (file in c("raw_replay.rds", "public_raw_replay.rds")) {
    retained <- if (file == "raw_replay.rds") raw else public
    for (field in names(semantics)) {
      changed <- retained
      changed$expected_semantics[[field]][1L] <- changed$expected_semantics[[field]][1L] + 1L
      changed$observed_semantics[[field]][1L] <- changed$expected_semantics[[field]][1L]
      saveRDS(changed, file.path(directory, file))
      update_receipt()
      fails(phase_history_rows(directory))
    }
    changed <- retained
    if (file == "raw_replay.rds") changed$output$source_indices[1L] <- 0L
    else changed$actual$evidence_flags[1L] <- 0L
    saveRDS(changed, file.path(directory, file))
    update_receipt()
    fails(phase_history_rows(directory))
    saveRDS(retained, file.path(directory, file))
  }
  for (change in c("gt_order", "short_buffer", "long_sequence", "negative_sequence", "edit_mask")) {
    changed <- raw
    if (change == "gt_order") changed$output[[3L]][c(1L, 2L * n - 1L)] <-
      changed$output[[3L]][c(2L * n - 1L, 1L)]
    if (change == "short_buffer") changed$output$cds <- changed$output$cds[-1L]
    if (change == "long_sequence") changed$output$cds_lengths[1L] <- capacity + 1L
    if (change == "negative_sequence") changed$output$protein_lengths[1L] <- -1L
    if (change == "edit_mask") changed$output$edit_masks[1L] <- 4L
    saveRDS(changed, file.path(directory, "raw_replay.rds"))
    update_receipt()
    fails(phase_history_rows(directory))
  }
  saveRDS(raw, file.path(directory, "raw_replay.rds"))
  changed <- public
  changed$actual <- public$actual[rev(seq_len(nrow(public$actual))), ]
  saveRDS(changed, file.path(directory, "public_raw_replay.rds"))
  update_receipt()
  reordered <- phase_history_rows(directory)
  stopifnot(sum(reordered$public_raw_replay_cases) == n,
    sum(reordered$public_raw_replay_disagreements) == 0L)
  saveRDS(public, file.path(directory, "public_raw_replay.rds"))
  for (missing_lanes in 1:2) {
    absent_public <- public
    absent_public$actual <- public$actual[-seq_len(missing_lanes), ]
    absent_public$comparisons <- phase_public_comparisons(
      cases, oracle, absent_public$actual, native$records)
    absent_public$observed_semantics <- phase_public_semantics(
      cases, absent_public$actual, native$records)
    stopifnot(all(is.na(absent_public$observed_semantics[seq_len(2L * missing_lanes), ])))
    saveRDS(absent_public, file.path(directory, "public_raw_replay.rds"))
    write.csv(transform(cases, equal = seq_len(n) != 1L),
      file.path(directory, "public_raw_replay_summary.csv"), row.names = FALSE)
    update_receipt(modifyList(receipt, list(public_raw_replay_disagreements = 1L,
      public_raw_record_disagreements = 2L * missing_lanes)))
    absent_rows <- phase_history_rows(directory)
    stopifnot(sum(absent_rows$public_raw_replay_cases) == n,
      sum(absent_rows$public_raw_replay_disagreements) == 1L,
      sum(absent_rows$public_raw_record_observations) == 4L * n,
      sum(absent_rows$public_raw_record_disagreements) == 2L * missing_lanes)
  }
  saveRDS(public, file.path(directory, "public_raw_replay.rds"))
  write_replay_summaries(cases)
  for (field in c("transcript_index", "haplotype_lane", "alt_index")) {
    invalid <- public$actual
    if (field == "transcript_index") invalid$transcript_index[1L] <- 0.5
    if (field == "haplotype_lane") invalid$carriers[[1L]]$haplotype_lane <- 1.5
    if (field == "alt_index") invalid$contributors[[1L]]$alt_index <- 0.5
    fails(phase_public_semantics(cases, invalid, native$records))
  }
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
