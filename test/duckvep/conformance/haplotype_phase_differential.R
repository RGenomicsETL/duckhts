#!/usr/bin/env Rscript
# Finite raw-GT audit against unmodified Haplosaurus. Disagreements are retained,
# not accepted as conformance. The existing seeded replay corpus is untouched.

source('test/duckvep/conformance/haplotype_observations.R')

genotypes <- function(max_ploidy) {
  unlist(lapply(seq_len(max_ploidy), function(n) {
    alleles <- as.matrix(expand.grid(rep(list(c("0", "1", "2", ".")), n),
      stringsAsFactors = FALSE))
    separators <- if (n == 1L) matrix(character(), 1L, 0L) else
      as.matrix(expand.grid(rep(list(c("|", "/")), n - 1L), stringsAsFactors = FALSE))
    unlist(lapply(seq_len(nrow(alleles)), function(i) {
      vapply(seq_len(nrow(separators)), function(j) {
        paste0(paste0(alleles[i, ], c(separators[j, ], "")), collapse = "")
      }, "")
    }), use.names = FALSE)
  }), use.names = FALSE)
}

phase_equal <- function(expected, observed) {
  stopifnot(identical(names(expected), names(observed)), nrow(expected) == nrow(observed))
  Reduce(`&`, Map(`==`, expected, observed))
}

phase_check_source_records <- function(directory, cases, records) {
  text <- read.delim(file.path(directory, "calls.vcf"), header = FALSE, comment.char = "#",
    col.names = c("chrom", "position", "id", "ref", "alt", "qual", "filter", "info", "format", "sample"),
    colClasses = "character", quote = "", fill = FALSE)
  n <- nrow(cases)
  gt <- as.vector(rbind(cases$GT, vapply(cases$ploidy,
    function(p) paste(rep("1", p), collapse = "|"), "")))
  stopifnot(identical(cases$transcript_index, seq_len(n) - 1L),
    identical(cases$seq_region, cases$transcript_index),
    identical(cases$chrom, sprintf("chrP%05d", seq_len(n))),
    identical(cases$transcript, sprintf("HP%05d", seq_len(n))),
    nrow(text) == 2L * n, nrow(records) == nrow(text),
    identical(text$chrom, rep(cases$chrom, each = 2L)),
    identical(text$position, rep(c("41", "44"), n)),
    identical(text$id, rep(c("a", "b"), n)), all(text$ref == "G"),
    identical(text$alt, rep(c("A,T", "C"), n)), all(text$format == "GT:PS"),
    identical(text$sample, paste0(gt, rep(c(":10", ":20"), n))),
    !anyDuplicated(records$record_index))
  at <- match(seq_len(nrow(text)) - 1L, records$record_index)
  stopifnot(!anyNA(at))
  records <- records[at, , drop = FALSE]
  reader_gt <- vapply(records$calls, function(call) {
    stopifnot(is.data.frame(call), nrow(call) == 1L,
      is.numeric(call$sample_index), length(call$sample_index) == 1L,
      !is.na(call$sample_index), call$sample_index == 0,
      is.character(call$raw_gt), length(call$raw_gt) == 1L, !is.na(call$raw_gt))
    call$raw_gt
  }, "")
  stopifnot(identical(records$CHROM, text$chrom),
    identical(as.character(records$POS), text$position), identical(records$ID, text$id),
    identical(records$REF, text$ref),
    identical(vapply(records$ALT, paste, "", collapse = ","), text$alt),
    identical(reader_gt, gt))
}

phase_decoded_comparisons <- function(cases, oracle, actual, records) {
  stopifnot(!anyDuplicated(names(oracle)), setequal(names(oracle), cases$transcript),
    all(actual$transcript_index %in% cases$transcript_index),
    !anyDuplicated(records$record_index))
  rows <- split(seq_len(nrow(actual)), actual$transcript_index)
  lapply(seq_len(nrow(cases)), function(i) {
    o <- oracle[[cases$transcript[i]]]
    a <- actual[rows[[as.character(cases$transcript_index[i])]], , drop = FALSE]
    observed <- lapply(seq_len(nrow(a)), function(j) list(cds = a$cds[j], protein = a$protein[j],
      count = a$carrier_count[j], contributors = records$ID[match(
        (a$contributors[[j]]$event_index - 1) %/% 3, records$record_index)]))
    expected <- canonical(o$haplotypes)
    observed <- canonical(observed)
    list(expected = expected, observed = observed, equal = identical(expected, observed),
      oracle_lanes = o$total_haplotype_count, native_lanes = sum(a$carrier_count),
      native_unknown = sum(is.na(a$cds)),
      native_unavailable_carriers = sum(a$carrier_count[is.na(a$cds)]))
  })
}

phase_raw_comparisons <- function(cases, oracle, replay) {
  n <- nrow(cases)
  nlanes <- 2L * n
  stopifnot(!anyDuplicated(names(oracle)), setequal(names(oracle), cases$transcript))
  # The retained .C result starts with the nine inputs to
  # duckhts_test_raw_phase_haplotypes; these positions bind each output lane.
  input_names <- c("reference", "genomic_start", "gt", "positions", "refs", "alts",
    "alt_counts", "nprofiles", "capacity")
  input <- setNames(replay[seq_along(input_names)], input_names)
  gt <- as.vector(rbind(cases$GT, vapply(cases$ploidy,
    function(p) paste(rep("1", p), collapse = "|"), "")))
  stopifnot(identical(input$reference, unique(cases$cds)), identical(input$genomic_start, 11L),
    identical(input$gt, gt), identical(input$positions, c(41L, 44L)),
    identical(input$refs, c("G", "G")), identical(input$alts, c("A", "T", "C")),
    identical(input$alt_counts, c(2L, 1L)), identical(input$nprofiles, n),
    is.integer(input$capacity), length(input$capacity) == 1L,
    !is.na(input$capacity), input$capacity > 0L)
  capacity <- input$capacity
  for (axis in c("cds", "protein")) {
    sizes <- replay[[paste0(axis, "_lengths")]]
    stopifnot(is.raw(replay[[axis]]), length(replay[[axis]]) == as.double(nlanes) * capacity,
      is.integer(sizes), length(sizes) == nlanes, !anyNA(sizes),
      all(sizes >= 0L & sizes <= capacity))
  }
  stopifnot(is.integer(replay$edit_masks), length(replay$edit_masks) == nlanes,
    !anyNA(replay$edit_masks), all(replay$edit_masks %in% 0:3),
    is.integer(replay$errors), length(replay$errors) == n, !anyNA(replay$errors))
  sequences <- lapply(seq_len(nlanes), function(i) list(
    cds = rawToChar(replay$cds[(i - 1L) * capacity + seq_len(replay$cds_lengths[i])]),
    protein = rawToChar(replay$protein[(i - 1L) * capacity + seq_len(replay$protein_lengths[i])]),
    count = 1, contributors = c("a", "b")[bitwAnd(replay$edit_masks[i], c(1L, 2L)) != 0L]))
  lapply(seq_len(n), function(i) {
    expected <- canonical(oracle[[cases$transcript[i]]]$haplotypes)
    observed <- canonical(sequences[2L * i - c(1L, 0L)])
    list(expected = expected, observed = observed,
      equal = replay$errors[i] == 0L && identical(expected, observed))
  })
}

phase_public_comparisons <- function(cases, oracle, actual, records) {
  stopifnot(!anyDuplicated(names(oracle)), setequal(names(oracle), cases$transcript),
    is.data.frame(actual), all(c("transcript_index", "cds", "protein", "carrier_count",
      "coding_blocks") %in% names(actual)),
    is.numeric(actual$transcript_index), all(is.finite(actual$transcript_index)),
    all(actual$transcript_index == floor(actual$transcript_index)),
    all(actual$transcript_index %in% cases$transcript_index),
    !anyDuplicated(records$record_index))
  rows <- split(seq_len(nrow(actual)), actual$transcript_index)
  lapply(seq_len(nrow(cases)), function(i) {
    a <- actual[rows[[as.character(cases$transcript_index[i])]], , drop = FALSE]
    observed <- lapply(seq_len(nrow(a)), function(j) {
      blocks <- a$coding_blocks[[j]]
      stopifnot(is.null(blocks) || (is.data.frame(blocks) && "event_indices" %in% names(blocks)))
      ids <- unlist(blocks$event_indices, use.names = FALSE)
      source_rows <- match(ids, records$record_index)
      stopifnot(!anyNA(source_rows), all(records$CHROM[source_rows] == cases$chrom[i]))
      list(cds = a$cds[j], protein = a$protein[j], count = a$carrier_count[j],
        contributors = records$ID[source_rows])
    })
    expected <- canonical(oracle[[cases$transcript[i]]]$haplotypes)
    observed <- canonical(observed)
    list(expected = expected, observed = observed, equal = identical(expected, observed))
  })
}

phase_parser_expected <- function(cases, phase) {
  stopifnot(!anyDuplicated(names(phase)), setequal(names(phase), cases$transcript))
  do.call(rbind, lapply(seq_len(nrow(cases)), function(i) {
    p <- phase[[cases$transcript[i]]]
    stopifnot(p$default_ploidy == 2L, identical(names(p$sample_ploidy), "sample"),
      p$sample_ploidy$sample == 2L)
    ids <- vapply(p$calls, `[[`, "", "source_id")
    stopifnot(!anyDuplicated(ids), all(ids %in% c("a", "b")),
      all(vapply(p$calls, function(x) identical(x$sample, "sample"), TRUE)))
    do.call(rbind, lapply(c("a", "b"), function(id) {
      call <- p$calls[ids == id]
      alleles <- if (id == "a") c("G", "A", "T") else c("G", "C")
      indices <- if (length(call)) match(unlist(call[[1L]]$genotype), alleles) - 1L else integer()
      stopifnot(!anyNA(indices))
      data.frame(status = 0L, retained = length(call) == 1L, ploidy = cases$ploidy[i],
        missing = as.integer(id == "a" && cases$missing[i]), slots = length(indices),
        first = if (length(indices)) indices[1L] else -1L,
        second = if (length(indices) > 1L) indices[2L] else -1L)
    }))
  }))
}

phase_expected_semantics <- function(expected_phase) {
  stopifnot(nrow(expected_phase) %% 2L == 0L)
  n <- nrow(expected_phase) %/% 2L
  nlanes <- 2L * n
  sources <- rep(-2L, 2L * nlanes)
  evidence <- integer(2L * nlanes)
  for (i in seq_len(n)) for (lane in 1:2) for (record in 1:2) {
    call <- expected_phase[2L * (i - 1L) + record, ]
    at <- (2L * (i - 1L) + lane - 1L) * 2L + record
    if (!call$retained) {
      if (call$missing) {
        sources[at] <- 0L
        evidence[at] <- 10L
      }
      next
    }
    allele <- if (lane == 1L) call$first else call$second
    if (allele == 0L && !call$missing) next
    sources[at] <- allele
    evidence[at] <- if (allele == -1L) 8L else if (allele > 0L) 1L else 0L
    if (call$missing) evidence[at] <- bitwOr(evidence[at], 10L)
  }
  lane_evidence <- bitwOr(evidence[seq(1L, 2L * nlanes, 2L)],
    evidence[seq(2L, 2L * nlanes, 2L)])
  data.frame(source_indices = sources, source_evidence = evidence,
    evidence = rep(lane_evidence, each = 2L),
    sequence_status = rep(ifelse(bitwAnd(lane_evidence, 8L) != 0L, 8L, 0L), each = 2L))
}

phase_raw_semantics <- function(replay) {
  nlanes <- 2L * length(replay$errors)
  for (field in c("source_indices", "source_evidence", "evidence", "sequence_status")) {
    count <- if (field %in% c("source_indices", "source_evidence")) 2L * nlanes else nlanes
    stopifnot(is.integer(replay[[field]]), length(replay[[field]]) == count, !anyNA(replay[[field]]))
  }
  data.frame(source_indices = replay$source_indices, source_evidence = replay$source_evidence,
    evidence = rep(replay$evidence, each = 2L), sequence_status = rep(replay$sequence_status, each = 2L))
}

phase_public_semantics <- function(cases, actual, records) {
  stopifnot(is.data.frame(actual), all(c("transcript_index", "carrier_count", "carriers",
    "contributors", "evidence_flags", "sequence_status") %in% names(actual)),
    is.numeric(actual$transcript_index), all(is.finite(actual$transcript_index)),
    all(actual$transcript_index == floor(actual$transcript_index)),
    all(actual$transcript_index %in% cases$transcript_index))
  semantics <- data.frame(source_indices = rep(NA_integer_, 4L * nrow(cases)),
    source_evidence = NA_integer_, evidence = NA_integer_, sequence_status = NA_integer_)
  seen <- logical(2L * nrow(cases))
  for (i in seq_len(nrow(actual))) {
    a <- actual[i, ]
    carriers <- a$carriers[[1L]]
    stopifnot(is.data.frame(carriers), all(c("sample_index", "ploidy", "phase_set",
      "haplotype_lane") %in% names(carriers)), is.numeric(a$carrier_count),
      is.finite(a$carrier_count), a$carrier_count > 0, a$carrier_count == nrow(carriers),
      all(carriers$sample_index == 0L), all(carriers$ploidy == 2L),
      all(is.na(carriers$phase_set)), is.numeric(carriers$haplotype_lane),
      all(carriers$haplotype_lane %in% 1:2))
    for (lane in carriers$haplotype_lane) {
      key <- 2L * a$transcript_index + lane
      stopifnot(!seen[key])
      seen[key] <- TRUE
      at <- 2L * (key - 1L) + 1:2
      semantics$source_indices[at] <- -2L
      semantics$source_evidence[at] <- 0L
      semantics$evidence[at] <- a$evidence_flags
      semantics$sequence_status[at] <- match(a$sequence_status, c("ok", "conditional")) * 8L - 8L
      contributors <- a$contributors[[1L]]
      stopifnot(is.data.frame(contributors), all(c("event_index", "seq_region", "position",
        "reference", "alternate", "alt_index", "evidence_flags") %in% names(contributors)),
        is.numeric(contributors$alt_index))
      source_rows <- match(contributors$event_index, records$record_index)
      stopifnot(!anyNA(source_rows),
        all(records$CHROM[source_rows] == cases$chrom[a$transcript_index + 1L]),
        all(contributors$seq_region == cases$seq_region[a$transcript_index + 1L]),
        all(contributors$position == records$POS[source_rows]),
        identical(contributors$reference, records$REF[source_rows]))
      for (j in seq_len(nrow(contributors))) {
        ordinal <- contributors$alt_index[j]
        stopifnot(is.na(ordinal) || (is.finite(ordinal) && ordinal == floor(ordinal) &&
          ordinal >= 0L && ordinal <= length(records$ALT[[source_rows[j]]])))
        expected_alt <- if (is.na(ordinal)) "" else if (ordinal == 0L) records$REF[source_rows[j]]
          else records$ALT[[source_rows[j]]][ordinal]
        stopifnot(identical(contributors$alternate[j], expected_alt))
      }
      record <- match(records$ID[source_rows], c("a", "b"))
      stopifnot(!anyNA(record), !anyDuplicated(record))
      semantics$source_indices[at[record]] <- ifelse(is.na(contributors$alt_index), -1L,
        as.integer(contributors$alt_index))
      semantics$source_evidence[at[record]] <- contributors$evidence_flags
    }
  }
  semantics
}

phase_history_rows <- function(directory) {
  source("scripts/duckvep_evidence.R", local = TRUE)
  directory <- normalizePath(directory, mustWork = TRUE)
  artifact_directory <- duckvep_evidence_repo_path(getwd(), directory)
  stopifnot(nzchar(artifact_directory))
  receipt_path <- file.path(directory, "receipt.json")
  receipt <- jsonlite::fromJSON(receipt_path)
  stopifnot(receipt$extension_build_binding == "htslib_distclean_make_release",
    receipt$scope == "raw_GT_finite_phase_audit_not_conformance",
    grepl("^[0-9a-f]{40}$", receipt$source_revision), receipt$max_ploidy %in% 2:4)
  hashes <- unlist(receipt$sha256)
  files <- c("cases.tsv", "summary.csv", "comparisons.rds", "phase_comparisons.rds", "raw_replay.rds",
    "public_raw_replay.rds", "controls.csv", "phase_controls.csv", "raw_replay_controls.csv",
    "decoded_collisions.rds", "native.rds", "calls.vcf", "oracle.stdout", "raw_replay_summary.csv",
    "public_raw_replay_summary.csv", "phase.jsonl")
  paths <- normalizePath(file.path(directory, files), mustWork = TRUE)
  absolute <- startsWith(names(hashes), "/") | grepl("^[A-Za-z]:", names(hashes))
  recorded <- normalizePath(ifelse(absolute, names(hashes), file.path(getwd(), names(hashes))),
    mustWork = FALSE)
  stopifnot(!anyDuplicated(recorded), all(paths %in% recorded))
  # Verify retained observations before reading them. Source and binary hashes
  # identify the execution receipt, not the checkout used to publish it.
  retained <- list.files(directory, recursive = TRUE, full.names = TRUE)
  retained <- retained[basename(retained) != "receipt.json"]
  in_artifact <- startsWith(recorded, paste0(directory, .Platform$file.sep))
  stopifnot(setequal(recorded[in_artifact], normalizePath(retained)))
  at <- match(normalizePath(retained), recorded)
  stopifnot(!anyNA(at), identical(unname(vapply(retained, duckvep_evidence_sha256, "")),
    unname(hashes[at])))
  read <- function(name) readRDS(file.path(directory, name))
  summary <- read.csv(file.path(directory, "summary.csv"), stringsAsFactors = FALSE)
  cases <- read.delim(file.path(directory, "cases.tsv"), stringsAsFactors = FALSE)
  keys <- setdiff(names(cases), "cds")
  stopifnot(all(c("transcript_index", "seq_region", "chrom", "transcript", "GT",
    "ploidy", "prefix", "missing", "mixed", "cds") %in% names(cases)),
    identical(summary[keys], cases[keys]))
  n <- nrow(summary)
  stopifnot(n == 3 * sum(4^(1:receipt$max_ploidy) * 2^(0:(receipt$max_ploidy - 1L))),
    !anyDuplicated(summary$GT), !anyNA(summary),
    all(grepl("^[|/]?[012.]([|/][012.]){0,3}$", summary$GT)),
    identical(summary$transcript_index, seq_len(n) - 1L),
    identical(summary$ploidy, lengths(strsplit(sub("^[|/]", "", summary$GT), "[|/]"))),
    max(summary$ploidy) == receipt$max_ploidy,
    identical(summary$prefix, grepl("^[|/]", summary$GT)),
    identical(summary$missing, grepl(".", summary$GT, fixed = TRUE)),
    identical(summary$mixed, grepl("|", summary$GT, fixed = TRUE) &
      grepl("/", summary$GT, fixed = TRUE)))
  decoded <- read("comparisons.rds")
  parser <- read("phase_comparisons.rds")
  raw <- read("raw_replay.rds")
  public <- read("public_raw_replay.rds")
  native <- read("native.rds")
  phase_check_source_records(directory, cases, native$records)
  oracle <- lapply(readLines(file.path(directory, "oracle.stdout")), jsonlite::fromJSON,
    simplifyVector = FALSE)
  names(oracle) <- vapply(oracle, `[[`, "", "transcript")
  phase <- lapply(readLines(file.path(directory, "phase.jsonl")), jsonlite::fromJSON,
    simplifyVector = FALSE)
  names(phase) <- vapply(phase, `[[`, "", "transcript")
  stopifnot(identical(decoded, phase_decoded_comparisons(cases, oracle, native$actual, native$records)),
    identical(raw$comparisons, phase_raw_comparisons(cases, oracle, raw$output)),
    identical(public$comparisons, phase_public_comparisons(cases, oracle, public$actual, native$records)))
  expected_phase <- phase_parser_expected(cases, phase)
  expected_semantics <- phase_expected_semantics(expected_phase)
  stopifnot(identical(parser$expected, expected_phase),
    identical(parser$observed$retained, parser$disposition == 3L),
    identical(raw$expected_semantics, expected_semantics),
    identical(public$expected_semantics, expected_semantics),
    identical(raw$observed_semantics, phase_raw_semantics(raw$output)),
    identical(public$observed_semantics, phase_public_semantics(cases, public$actual, native$records)))
  raw_gt <- as.vector(rbind(cases$GT, vapply(cases$ploidy,
    function(p) paste(rep("1", p), collapse = "|"), "")))
  stopifnot(identical(parser$keys, data.frame(transcript = rep(cases$transcript, each = 2L),
    source_id = rep(c("a", "b"), n), GT = raw_gt)))
  carrier_counts <- function(actual) {
    stopifnot(is.data.frame(actual), all(c("transcript_index", "carrier_count") %in% names(actual)),
      is.numeric(actual$transcript_index), is.numeric(actual$carrier_count),
      !anyNA(actual[c("transcript_index", "carrier_count")]),
      all(actual$transcript_index %in% summary$transcript_index),
      all(is.finite(actual$carrier_count) & actual$carrier_count > 0 &
        actual$carrier_count == floor(actual$carrier_count)))
    unname(vapply(split(actual$carrier_count,
      factor(actual$transcript_index, levels = summary$transcript_index)), sum, 0))
  }
  groups_valid <- function(groups) is.list(groups) && all(vapply(groups, function(group) {
    nullable_string <- function(x) is.null(x) || (is.character(x) && length(x) == 1L)
    is.list(group) && !anyDuplicated(names(group)) &&
      all(c("cds", "protein", "count", "contributors") %in% names(group)) &&
      nullable_string(group$cds) && nullable_string(group$protein) &&
      is.numeric(group$count) && length(group$count) == 1L && is.finite(group$count) &&
      group$count > 0 && group$count == floor(group$count) &&
      is.character(group$contributors) && !anyNA(group$contributors)
  }, TRUE))
  comparison_matches <- function(x, observed_counts, errors = integer(n)) {
    stopifnot(length(x) == n, length(errors) == n, !anyNA(errors), length(observed_counts) == n)
    stopifnot(all(vapply(x, function(row) is.list(row) && !anyDuplicated(names(row)) &&
      all(c("expected", "observed", "equal") %in% names(row)) &&
      is.logical(row$equal) && length(row$equal) == 1L && !is.na(row$equal) &&
      groups_valid(row$expected) && groups_valid(row$observed), TRUE)))
    counts <- function(groups) sum(vapply(groups, `[[`, 0, "count"))
    stopifnot(all(vapply(x, function(row) counts(row$expected), 0) == 2),
      identical(vapply(x, function(row) counts(row$observed), 0), as.numeric(observed_counts)))
    result <- vapply(x, function(row) identical(row$expected, row$observed), TRUE) & errors == 0L
    stopifnot(identical(result, vapply(x, `[[`, TRUE, "equal")))
    result
  }
  stopifnot(all(summary$oracle_lanes == 2),
    identical(as.numeric(summary$native_lanes), carrier_counts(native$actual)),
    identical(summary$equal, comparison_matches(decoded, summary$native_lanes)))
  for (field in c("oracle_lanes", "native_lanes", "native_unknown", "native_unavailable_carriers"))
    stopifnot(identical(as.numeric(summary[[field]]), vapply(decoded, `[[`, 0, field)))
  row_matches <- function(expected, observed, count, fields) {
    stopifnot(is.data.frame(expected), is.data.frame(observed), nrow(expected) == count,
      identical(names(expected), fields), identical(names(observed), fields))
    result <- phase_equal(expected, observed)
    result[is.na(result)] <- FALSE
    result
  }
  parser_fields <- c("status", "retained", "ploidy", "missing", "slots", "first", "second")
  record_fields <- c("source_indices", "source_evidence", "evidence", "sequence_status")
  parser_matches <- row_matches(parser$expected, parser$observed, 2L * n, parser_fields)
  raw_matches <- comparison_matches(raw$comparisons, rep(2, n), raw$output$errors)
  public_matches <- comparison_matches(public$comparisons, carrier_counts(public$actual))
  raw_records <- row_matches(raw$expected_semantics, raw$observed_semantics, 4L * n, record_fields)
  public_records <- row_matches(public$expected_semantics, public$observed_semantics, 4L * n, record_fields)
  stopifnot(identical(raw$expected_semantics, public$expected_semantics))
  for (kind in c("raw", "public_raw")) {
    lane <- read.csv(file.path(directory, paste0(kind, "_replay_summary.csv")), stringsAsFactors = FALSE)
    stopifnot(identical(lane[keys], cases[keys]),
      identical(lane$equal, if (kind == "raw") raw_matches else public_matches))
  }
  controls <- function(file, required) {
    x <- read.csv(file.path(directory, file), stringsAsFactors = FALSE)
    stopifnot(identical(names(x), c("control", "rejected")), !anyDuplicated(x$control),
      all(required %in% x$control), is.logical(x$rejected), all(x$rejected))
    nrow(x)
  }
  control_counts <- c(controls_rejected = controls("controls.csv",
      c("duplicate", "cds", "protein", "contributor")),
    raw_parser_controls_rejected = controls("phase_controls.csv", parser_fields),
    raw_replay_controls_rejected = controls("raw_replay_controls.csv", record_fields))
  values <- data.frame(cases = 1L, disagreements = as.integer(!summary$equal),
    summary[c("oracle_lanes", "native_lanes", "native_unavailable_carriers")],
    raw_parser_calls = 2L, raw_parser_disagreements = rowSums(matrix(!parser_matches, ncol = 2L, byrow = TRUE)),
    raw_replay_cases = 1L, raw_replay_disagreements = as.integer(!raw_matches),
    raw_replay_record_observations = 4L,
    raw_replay_record_disagreements = rowSums(matrix(!raw_records, ncol = 4L, byrow = TRUE)),
    public_raw_replay_cases = 1L, public_raw_replay_disagreements = as.integer(!public_matches),
    public_raw_record_observations = 4L,
    public_raw_record_disagreements = rowSums(matrix(!public_records, ncol = 4L, byrow = TRUE)))
  totals <- c(colSums(values), control_counts,
    decoded_collision_groups = length(read("decoded_collisions.rds")),
    input_records = nrow(native$records), source_alt_events = sum(lengths(native$records$ALT)),
    input_genotype_calls = sum(vapply(native$records$calls, nrow, 1L)),
    input_allele_slots = sum(vapply(native$records$calls,
      function(calls) sum(lengths(calls$alleles)), 0)),
    candidate_alt_calls = nrow(native$calls), native_leaves = nrow(native$actual))
  for (field in names(totals)) stopifnot(identical(as.numeric(receipt[[field]]), totals[[field]]))
  extension_hash <- hashes[endsWith(names(hashes), "/duckhts.duckdb_extension")]
  stopifnot(length(extension_hash) == 1L, grepl("^[0-9a-f]{64}$", extension_hash),
    all(grepl("^[0-9a-f]{40}$", unlist(receipt$oracle_revisions))))
  strata <- aggregate(values, summary[c("ploidy", "prefix", "missing", "mixed")], sum)
  metadata <- data.frame(source_revision = receipt$source_revision,
    extension_build_binding = receipt$extension_build_binding, extension_sha256 = unname(extension_hash),
    oracle_vep_revision = receipt$oracle_revisions$vep,
    oracle_variation_revision = receipt$oracle_revisions$variation,
    receipt_sha256 = duckvep_evidence_sha256(receipt_path),
    artifact_directory = artifact_directory, max_ploidy = receipt$max_ploidy,
    decoded_collision_groups = receipt$decoded_collision_groups)
  cbind(metadata[rep(1L, nrow(strata)), ], strata,
    as.data.frame(as.list(control_counts)), row.names = NULL)
}

publish_phase_history <- function(directory, history_path) {
  rows <- phase_history_rows(directory)
  stopifnot(file.exists(history_path), !dir.exists(history_path))
  lock <- paste0(history_path, ".lock")
  if (!dir.create(lock, showWarnings = FALSE))
    stop("Phase history publication is busy: ", lock, call. = FALSE)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  history <- read.csv(history_path, stringsAsFactors = FALSE)
  stopifnot(setequal(names(history), names(rows)),
    !any(rows$source_revision %in% history$source_revision))
  rows <- rows[names(history)]
  output <- tempfile("phase-history-", tmpdir = dirname(history_path))
  on.exit(unlink(output), add = TRUE)
  write.csv(rbind(history, rows), output, row.names = FALSE)
  stopifnot(file.rename(output, history_path))
  message("Published ", sum(rows$cases), " profiles with ", sum(rows$disagreements),
    " retained decoded/raw disagreements; publication does not change their verdict.")
}

main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--max-ploidy", dest = "max_ploidy", type = "integer", default = 4L),
    optparse::make_option("--vep-prefix", dest = "vep_prefix",
      default = Sys.getenv("VEP_PREFIX", "/root/miniconda3/envs/vep")),
    optparse::make_option("--extension-receipt", dest = "extension_receipt", default = NULL),
    optparse::make_option("--publish-artifact", dest = "publish_artifact", default = NULL),
    optparse::make_option("--history", default = "test/duckvep/conformance/data/haplotype_phase_history.csv")
  )))
  if (!is.null(opt$publish_artifact)) return(publish_phase_history(opt$publish_artifact, opt$history))
  stopifnot(opt$max_ploidy >= 2L, opt$max_ploidy <= 4L)
  source("scripts/duckvep_evidence.R", local = TRUE)
  root <- normalizePath(".")
  revision <- duckvep_evidence_revision(root)
  extension <- normalizePath("build/release/duckhts.duckdb_extension")
  binding <- "diagnostic_unbound"
  if (!is.null(opt$extension_receipt)) {
    duckvep_evidence_assert_checkout(root, revision)
    binding <- duckvep_evidence_read_extension_receipt(opt$extension_receipt,
      root, extension, revision)$binding
  }
  mirrors <- c(vep = normalizePath(".sync/ensembl-vep"),
    variation = normalizePath(".sync/ensembl-variation"))
  pins <- c(vep = "57ea5c52340acc1f156267f810ad162e26597082",
    variation = "2fb834b987ede3824e200197a838ce11e91aeb4b")
  for (name in names(pins)) {
    stopifnot(identical(duckvep_evidence_command("git", c("-C", mirrors[[name]],
      "rev-parse", "HEAD"), "oracle revision"), pins[[name]]),
      !length(duckvep_evidence_command("git", c("-C", mirrors[[name]], "status",
        "--porcelain"), "oracle worktree")))
  }
  prefix <- normalizePath(opt$vep_prefix)
  out <- tempfile("haplotype_phase_", tmpdir = "test/duckvep/conformance/results")
  dir.create(out)
  message("Phase artifacts: ", out)
  environment <- duckvep_evidence_command("micromamba", c("list", "-p", prefix, "--explicit"),
    "oracle environment")
  writeLines(environment, file.path(out, "environment.txt"))
  stopifnot(identical(duckvep_evidence_explicit_packages(environment),
    duckvep_evidence_explicit_packages(readLines(
      "test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt"))))
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(root, "r/duckhtsbench/inst/benchmark_registry.tsv"))
  paths <- duckhtsbench::duckhts_bench_stage_repository_fixtures(root, "duckvep-haplotypes")
  cds <- readLines(paths[["haplotype_benchmark_reference"]])[2L]
  stopifnot(nchar(cds) == 180L, substr(cds, 31L, 31L) == "G", substr(cds, 34L, 34L) == "G")
  plain <- genotypes(opt$max_ploidy)
  gt <- c(plain, paste0("|", plain), paste0("/", plain))
  stopifnot(!anyDuplicated(gt), length(gt) == 3 * sum(4^(1:opt$max_ploidy) * 2^(0:(opt$max_ploidy - 1L))))
  ploidy <- lengths(strsplit(sub("^[|/]", "", gt), "[|/]"))
  cases <- data.frame(transcript_index = seq_along(gt) - 1L, seq_region = seq_along(gt) - 1L,
    chrom = sprintf("chrP%05d", seq_along(gt)), transcript = sprintf("HP%05d", seq_along(gt)),
    GT = gt, ploidy, prefix = grepl("^[|/]", gt), missing = grepl(".", gt, fixed = TRUE),
    mixed = grepl("|", gt, fixed = TRUE) & grepl("/", gt, fixed = TRUE), cds = cds)
  write.table(cases, file.path(out, "cases.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(as.vector(rbind(paste0(">", cases$chrom), paste0(strrep("A", 10), cds, strrep("A", 10)))),
    file.path(out, "reference.fa"))
  gff <- unlist(lapply(seq_len(nrow(cases)), function(i) {
    id <- cases$transcript[i]
    attrs <- c(paste0("ID=gene:", id, ";biotype=protein_coding"),
      paste0("ID=transcript:", id, ";Parent=gene:", id, ";biotype=protein_coding"),
      paste0("ID=exon:", id, ";Parent=transcript:", id), paste0("Parent=transcript:", id))
    paste(cases$chrom[i], "phase", c("gene", "mRNA", "exon", "CDS"), 11, 190, ".", "+",
      c(".", ".", ".", "0"), attrs, sep = "\t")
  }))
  writeLines(c("##gff-version 3", gff), file.path(out, "model.gff3"))
  vcf <- unlist(lapply(seq_len(nrow(cases)), function(i) c(
    paste(cases$chrom[i], 41, "a", "G", "A,T", ".", "PASS", ".", "GT:PS", paste0(gt[i], ":10"), sep = "\t"),
    paste(cases$chrom[i], 44, "b", "G", "C", ".", "PASS", ".", "GT:PS",
      paste0(paste(rep("1", ploidy[i]), collapse = "|"), ":20"), sep = "\t"))))
  writeLines(c("##fileformat=VCFv4.4", paste0("##contig=<ID=", cases$chrom, ",length=200>"),
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample", vcf), file.path(out, "calls.vcf"))
  run <- function(command, args, name) {
    status <- system2(command, shQuote(args), stdout = file.path(out, paste0(name, ".stdout")),
      stderr = file.path(out, paste0(name, ".stderr")))
    if (status != 0L) stop(name, " failed; retained artifacts: ", out)
  }
  run("samtools", c("faidx", file.path(out, "reference.fa")), "faidx")
  run("bcftools", c("norm", "-c", "e", "-f", file.path(out, "reference.fa"),
    "-o", file.path(out, "ref_checked.vcf"), file.path(out, "calls.vcf")), "reference_check")
  run("bgzip", file.path(out, "model.gff3"), "bgzip")
  run("tabix", c("-p", "gff", file.path(out, "model.gff3.gz")), "tabix")
  perl_lib <- paste(c(file.path(mirrors, "modules"), file.path(prefix, "share/ensembl-vep-116.0-0")), collapse = ":")
  run("micromamba", c("run", "--clean-env", "--env", paste0("PERL5LIB=", perl_lib), "-p", prefix,
    "perl", "test/duckvep/conformance/haplotype_oracle.pl", file.path(out, "calls.vcf"),
    file.path(out, "reference.fa"), file.path(out, "model.gff3.gz"),
    file.path(out, "phase.jsonl")), "oracle")
  oracle <- lapply(readLines(file.path(out, "oracle.stdout")), jsonlite::fromJSON, simplifyVector = FALSE)
  names(oracle) <- vapply(oracle, `[[`, "", "transcript")
  stopifnot(!anyDuplicated(names(oracle)), setequal(names(oracle), cases$transcript))
  for (o in oracle) {
    stopifnot(sum(vapply(o$haplotypes, function(h) as.numeric(h$count), 0)) == o$total_haplotype_count,
      all(vapply(o$haplotypes, function(h)
        identical(names(h$samples), "sample") && h$samples$sample == h$count, TRUE)))
  }

  # Observe actual upstream objects, not a second implementation of its parser.
  # Source ploidy/missingness is independently known from the generated grammar;
  # effective file ploidy and retained allele slots come from Haplosaurus.
  phase <- lapply(readLines(file.path(out, "phase.jsonl")), jsonlite::fromJSON, simplifyVector = FALSE)
  names(phase) <- vapply(phase, `[[`, "", "transcript")
  stopifnot(!anyDuplicated(names(phase)), setequal(names(phase), cases$transcript))
  expected_phase <- phase_parser_expected(cases, phase)
  raw_gt <- as.vector(rbind(cases$GT, vapply(cases$ploidy,
    function(n) paste(rep("1", n), collapse = "|"), "")))
  probe_sources <- c("test/duckvep/conformance/phase_probe.c", paste0("src/duckvep/kernel/src/duckvep_",
    c("phase", "haplotype", "carriers", "haplotype_stream", "classify", "codon", "coding", "projection", "delta"), ".c"))
  probe <- file.path(out, paste0("phase_probe", .Platform$dynlib.ext))
  compiler <- Sys.getenv("CC", "cc")
  run(compiler, c("-std=c11", "-O2", "-Wall", "-Wextra", "-Werror", "-Wpedantic",
    "-pedantic-errors", "-fPIC", "-shared",
    "-Isrc/duckvep/kernel/src", "-Isrc/duckvep/kernel/include", probe_sources, "-o", probe), "phase_compile")
  writeLines(duckvep_evidence_command(compiler, "--version", "phase compiler"),
    file.path(out, "phase_compiler.txt"))
  dll <- dyn.load(probe)
  parsed <- .C(getNativeSymbolInfo("duckhts_test_vep116_raw_gt", dll), raw_gt,
    rep(c(2L, 1L), nrow(cases)), as.integer(length(raw_gt)), status = integer(length(raw_gt)),
    disposition = integer(length(raw_gt)), ploidy = integer(length(raw_gt)),
    missing = integer(length(raw_gt)), slots = integer(length(raw_gt)),
    first = integer(length(raw_gt)), second = integer(length(raw_gt)))
  observed_phase <- as.data.frame(parsed[c("status", "ploidy", "missing", "slots", "first", "second")])
  observed_phase$retained <- parsed$disposition == 3L
  observed_phase <- observed_phase[names(expected_phase)]
  phase_matches <- phase_equal(expected_phase, observed_phase)
  phase_keys <- data.frame(transcript = rep(cases$transcript, each = 2L),
    source_id = rep(c("a", "b"), nrow(cases)), GT = raw_gt)
  saveRDS(list(keys = phase_keys, expected = expected_phase, observed = observed_phase,
    disposition = parsed$disposition), file.path(out, "phase_comparisons.rds"))
  write.csv(cbind(phase_keys, equal = phase_matches), file.path(out, "phase_summary.csv"), row.names = FALSE)
  write.csv(cbind(phase_keys[!phase_matches, ], expected_phase[!phase_matches, ],
    observed_phase[!phase_matches, ]), file.path(out, "phase_mismatches.csv"), row.names = FALSE)
  phase_rejected <- vapply(names(expected_phase), function(field) {
    corrupt <- expected_phase
    corrupt[[field]][1L] <- if (is.logical(corrupt[[field]])) !corrupt[[field]][1L] else
      corrupt[[field]][1L] + 1L
    !all(phase_equal(expected_phase, corrupt))
  }, TRUE)
  stopifnot(all(phase_rejected))
  write.csv(data.frame(control = names(phase_rejected), rejected = phase_rejected),
    file.path(out, "phase_controls.csv"), row.names = FALSE)

  # Raw-record native replay and public decoded-call replay have independent
  # verdicts. The bridge binds the two source sites verbatim.
  nlanes <- 2L * nrow(cases)
  capacity <- 512L
  replay <- .C(getNativeSymbolInfo("duckhts_test_raw_phase_haplotypes", dll), as.character(cds),
    11L, raw_gt, c(41L, 44L), c("G", "G"), c("A", "T", "C"), c(2L, 1L),
    as.integer(nrow(cases)), capacity, cds = raw(nlanes * capacity), protein = raw(nlanes * capacity),
    cds_lengths = integer(nlanes), protein_lengths = integer(nlanes),
    sequence_status = integer(nlanes), evidence = integer(nlanes), edit_masks = integer(nlanes),
    source_indices = integer(2L * nlanes), source_evidence = integer(2L * nlanes), errors = integer(nrow(cases)))
  dyn.unload(dll[["path"]])
  raw_comparisons <- phase_raw_comparisons(cases, oracle, replay)
  raw_matches <- vapply(raw_comparisons, `[[`, TRUE, "equal")
  # Expected record observations come from the upstream object sidecar, not
  # the native parser. Missing omissions retain a conditional REF observation;
  # only physical edits participate in upstream's contributing-variant list.
  expected_semantics <- phase_expected_semantics(expected_phase)
  observed_semantics <- phase_raw_semantics(replay)
  semantics_matches <- phase_equal(expected_semantics, observed_semantics)
  raw_rejected <- vapply(names(expected_semantics), function(field) {
    corrupt <- expected_semantics
    corrupt[[field]][1L] <- corrupt[[field]][1L] + 1L
    !all(phase_equal(expected_semantics, corrupt))
  }, TRUE)
  stopifnot(all(raw_rejected))
  saveRDS(list(output = replay, comparisons = raw_comparisons,
    expected_semantics = expected_semantics, observed_semantics = observed_semantics),
    file.path(out, "raw_replay.rds"))
  write.csv(data.frame(control = names(raw_rejected), rejected = raw_rejected),
    file.path(out, "raw_replay_controls.csv"), row.names = FALSE)
  write.csv(cbind(cases[setdiff(names(cases), "cds")], equal = raw_matches, error = replay$errors),
    file.path(out, "raw_replay_summary.csv"), row.names = FALSE)

  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbExecute(con, paste("LOAD", q(extension)))
  DBI::dbExecute(con, "SET threads=4")
  DBI::dbWriteTable(con, "models", cases)
  queries <- c("SELECT seq_region::UINTEGER seq_region FROM models ORDER BY seq_region",
    "SELECT transcript_index::UINTEGER transcript_index,seq_region::UINTEGER seq_region,
     11::UBIGINT transcript_start,190::UBIGINT transcript_end,1::TINYINT strand,0::UINTEGER gene_index,
     3::UBIGINT transcript_flags,11::UBIGINT cds_start,190::UBIGINT cds_end,cds::BLOB cds_sequence,
     1::UTINYINT codon_table FROM models ORDER BY transcript_index",
    "SELECT transcript_index::UINTEGER transcript_index,11::UBIGINT exon_start,190::UBIGINT exon_end,
     1::UBIGINT exon_cdna_start,180::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase
     FROM models ORDER BY transcript_index")
  stopifnot(DBI::dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('phase',",
    paste(q(queries), collapse = ","), ")"))$loaded)
  DBI::dbExecute(con, paste0("CREATE TABLE records AS SELECT * FROM read_geno(",
    q(normalizePath(file.path(out, "calls.vcf"))), ",raw_gt:=true)"))
  DBI::dbExecute(con, "CREATE TABLE calls AS SELECT r.record_index*3+a.i event_index,
    m.seq_region,r.POS AS position,r.REF reference,r.ALT[a.i] alternate,a.i alt_index,
    m.transcript_index,c.sample_index,c.alleles,c.phase_before,c.phase_set
    FROM records r JOIN models m ON m.chrom=r.CHROM,unnest(r.calls) u(c),range(1,len(r.ALT)+1) a(i)")
  records <- DBI::dbGetQuery(con, "SELECT * FROM records")
  calls <- DBI::dbGetQuery(con, "SELECT * FROM calls")
  stopifnot(nrow(records) == 2L * nrow(cases), nrow(calls) == 3L * nrow(cases),
    all(lengths(calls$alleles) == cases$ploidy[calls$transcript_index + 1L]))
  actual <- DBI::dbGetQuery(con, "SELECT * FROM duckvep_haplotypes('SELECT * FROM calls',
    'phase',phase_policy:='vep116_compat')")
  stopifnot(all(actual$carrier_count == vapply(actual$carriers, nrow, 1L)),
    all(vapply(actual$carriers, function(c) all(c$sample_index == 0L), TRUE)))
  saveRDS(list(records = records, calls = calls, actual = actual), file.path(out, "native.rds"))
  comparisons <- phase_decoded_comparisons(cases, oracle, actual, records)
  saveRDS(comparisons, file.path(out, "comparisons.rds"))
  # The independent text read checks every physical record ordinal. Production
  # source calls consume only the original GT retained by the VCF reader.
  phase_check_source_records(out, cases, records)
  DBI::dbExecute(con, "CREATE TABLE source_calls AS SELECT r.record_index event_index,
    m.seq_region,r.POS AS position,r.REF reference,r.ALT alternates,m.transcript_index,
    c.sample_index,c.raw_gt AS gt FROM records r JOIN models m ON m.chrom=r.CHROM,
    unnest(r.calls) u(c)")
  stopifnot(DBI::dbGetQuery(con, "SELECT count(*) n FROM source_calls")$n == nrow(records))
  public_raw <- DBI::dbGetQuery(con, "SELECT * FROM duckvep_haplotypes('SELECT * FROM source_calls',
    'phase',phase_policy:='vep116_compat',input_mode:='source_records')")
  saveRDS(public_raw, file.path(out, "public_raw_output.rds"))
  public_raw_comparisons <- phase_public_comparisons(cases, oracle, public_raw, records)
  public_semantics <- phase_public_semantics(cases, public_raw, records)
  public_matches <- vapply(public_raw_comparisons, `[[`, TRUE, "equal")
  public_semantics_matches <- phase_equal(expected_semantics, public_semantics)
  public_semantics_matches[is.na(public_semantics_matches)] <- FALSE
  saveRDS(list(actual = public_raw, comparisons = public_raw_comparisons,
    expected_semantics = expected_semantics, observed_semantics = public_semantics),
    file.path(out, "public_raw_replay.rds"))
  write.csv(cbind(cases[setdiff(names(cases), "cds")], equal = public_matches),
    file.path(out, "public_raw_replay_summary.csv"), row.names = FALSE)
  summary <- cases[setdiff(names(cases), "cds")]
  for (name in c("equal", "oracle_lanes", "native_lanes", "native_unknown", "native_unavailable_carriers"))
    summary[[name]] <- vapply(comparisons, `[[`, if (name == "equal") TRUE else 0, name)
  write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
  # A fully called ordinary diploid GT is a positive control, not an exclusion:
  # all other profiles and their failures remain in the complete comparison table.
  control <- with(summary, ploidy == 2L & !prefix & !mixed & !missing)
  stopifnot(sum(control) == 18L, all(summary$equal[control]))
  witness <- oracle[[cases$transcript[which(control)[1L]]]]$haplotypes
  controls <- list(duplicate = c(witness, witness[1L]), cds = witness,
    protein = witness, contributor = witness)
  controls$cds[[1L]]$cds <- paste0("C", substring(witness[[1L]]$cds, 2L))
  controls$protein[[1L]]$protein <- paste0("X", substring(witness[[1L]]$protein, 2L))
  controls$contributor[[1L]]$contributors <- c(witness[[1L]]$contributors, "deliberate_extra_event")
  rejected <- vapply(controls, function(x) !identical(canonical(x), canonical(witness)), TRUE)
  rejected <- c(rejected, haplotype_group_controls())
  stopifnot(all(rejected))
  write.csv(data.frame(control = names(rejected), rejected),
    file.path(out, "controls.csv"), row.names = FALSE)
  decoded <- DBI::dbGetQuery(con, "SELECT m.transcript_index,to_json(struct_pack(
    alleles:=c.alleles,phase_before:=c.phase_before))::VARCHAR decoded_gt
    FROM records r JOIN models m ON m.chrom=r.CHROM,unnest(r.calls) u(c) WHERE r.ID='a'
    ORDER BY m.transcript_index")
  stopifnot(identical(as.integer(decoded$transcript_index), cases$transcript_index))
  summary$decoded_gt <- decoded$decoded_gt
  summary$oracle_signature <- vapply(comparisons, function(x)
    jsonlite::toJSON(x$expected, auto_unbox = TRUE, na = "null"), "")
  collisions <- split(summary, summary$decoded_gt)
  collisions <- collisions[vapply(collisions, function(x) length(unique(x$oracle_signature)) > 1L, TRUE)]
  saveRDS(collisions, file.path(out, "decoded_collisions.rds"))
  inputs <- c(paths, extension, probe_sources,
    list.files("src/duckvep/kernel", pattern = "\\.(h|inc|def)$", recursive = TRUE, full.names = TRUE),
    "test/duckvep/conformance/haplotype_phase_differential.R",
    "test/duckvep/conformance/haplotype_observations.R",
    "test/duckvep/conformance/haplotype_oracle.pl", file.path(prefix,
      "share/ensembl-vep-116.0-0/Bio/EnsEMBL/IO/Parser/BaseVCF4.pm"))
  jsonlite::write_json(list(source_revision = revision, extension_build_binding = binding,
    oracle_revisions = as.list(pins), scope = "raw_GT_finite_phase_audit_not_conformance",
    max_ploidy = opt$max_ploidy, cases = nrow(summary), disagreements = sum(!summary$equal),
    decoded_collision_groups = length(collisions), controls_rejected = sum(rejected),
    raw_parser_calls = length(raw_gt), raw_parser_disagreements = sum(!phase_matches),
    raw_parser_controls_rejected = sum(phase_rejected),
    raw_replay_cases = length(raw_matches), raw_replay_disagreements = sum(!raw_matches),
    raw_replay_record_observations = length(semantics_matches),
    raw_replay_record_disagreements = sum(!semantics_matches),
    raw_replay_controls_rejected = sum(raw_rejected),
    public_raw_replay_cases = length(public_matches),
    public_raw_replay_disagreements = sum(!public_matches),
    public_raw_record_observations = length(public_semantics_matches),
    public_raw_record_disagreements = sum(!public_semantics_matches),
    input_records = nrow(records), source_alt_events = 3L * nrow(cases),
    input_genotype_calls = nrow(records), input_allele_slots = 2L * sum(cases$ploidy),
    candidate_alt_calls = nrow(calls), native_leaves = nrow(actual),
    oracle_lanes = sum(summary$oracle_lanes), native_lanes = sum(summary$native_lanes),
    native_unavailable_carriers = sum(summary$native_unavailable_carriers),
    sha256 = as.list(vapply(unique(c(inputs, list.files(out, full.names = TRUE))),
      duckvep_evidence_sha256, ""))), file.path(out, "receipt.json"), pretty = TRUE, auto_unbox = TRUE)
  if (!is.null(opt$extension_receipt)) duckvep_evidence_assert_checkout(root, revision)
  print(aggregate(cbind(cases = rep(1L, nrow(summary)), disagreements = as.integer(!summary$equal)) ~
    ploidy + prefix + missing + mixed, data = summary, FUN = sum), row.names = FALSE)
  message("Decoded-equivalence collision groups: ", length(collisions))
  message("Public raw-record replay: ", sum(!public_matches), "/", length(public_matches),
    " disagreements; source-record observations: ", sum(!public_semantics_matches), "/",
    length(public_semantics_matches))
  message("Raw parser: ", sum(!phase_matches), " disagreements / ", length(raw_gt), " calls")
  message("Raw native replay: ", sum(!raw_matches), " disagreements / ", length(raw_matches), " profiles")
  message("Raw record observations: ", sum(!semantics_matches), " disagreements / ", length(semantics_matches))
  if (any(!phase_matches)) stop("Raw parser disagreements retained: ", out, call. = FALSE)
  if (any(!raw_matches)) stop("Raw native replay disagreements retained: ", out, call. = FALSE)
  if (any(!semantics_matches)) stop("Raw record observations disagree: ", out, call. = FALSE)
  if (any(!public_matches)) stop("Public raw replay disagreements retained: ", out, call. = FALSE)
  if (any(!public_semantics_matches)) stop("Public raw record observations disagree: ", out, call. = FALSE)
  if (any(!summary$equal)) stop("Raw-GT compatibility disagreements retained: ", out, call. = FALSE)
}
if (sys.nframe() == 0L) main()
