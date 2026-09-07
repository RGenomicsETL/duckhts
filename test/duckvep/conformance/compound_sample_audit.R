#!/usr/bin/env Rscript
# Additional diagnostic, never a replacement for the unmodified cohort audit.
# Leaf-local rendering must not depend on which other samples share an edit prefix.
main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--audit-artifacts", dest = "audit_artifacts", type = "character")
  )))
  source("scripts/duckvep_evidence.R", local = TRUE)
  root <- normalizePath(".")
  revision <- duckvep_evidence_revision(root)
  audit <- normalizePath(opt$audit_artifacts, mustWork = TRUE)
  receipt <- jsonlite::read_json(file.path(audit, "receipt.json"), simplifyVector = TRUE)
  stopifnot(identical(receipt$bcftools_source_revision, "9adeaf4cfcc3bff40efca6237749fefb53391678"),
    identical(receipt$bcftools_source_tree, "src/bcftools-1.23"))
  input <- dirname(receipt$source_receipt)
  source_paths <- file.path(input, c("inputs.rds", "native.rds", "reference.fa", "model.gff3.gz", "carriers.vcf"))
  unpatched <- file.path(audit, "src/bcftools-1.23/bcftools")
  archive <- file.path(audit, "bcftools-source.tar")
  original_paths <- c(source_paths, unpatched, archive)
  hashes <- vapply(original_paths, duckvep_evidence_sha256, "")
  for (p in original_paths) stopifnot(identical(hashes[[p]], receipt$sha256[[p]]))
  out <- tempfile("compound_sample_", tmpdir = dirname(audit))
  dir.create(out)
  message("Sample-separability audit: ", out)
  utils::untar(archive, exdir = out)
  build <- file.path(out, "src/bcftools-1.23")
  patch <- normalizePath("test/duckvep/conformance/patches/bcftools_csq_prefix_local.patch")
  code <- c(patch, "test/duckvep/conformance/compound_sample_audit.R", "scripts/duckvep_evidence.R")
  code_hashes <- vapply(code, duckvep_evidence_sha256, "")
  stopifnot(identical(duckvep_evidence_sha256(file.path(build, "csq.c")),
    "f407e55bc307029b1a82ca6c6d0470423818d0b6e8551464eee030b5d53affe6"))
  stopifnot(system2("patch", shQuote(c("--fuzz=0", "-p1", "-d", build, "-i", patch)),
    stdout = file.path(out, "patch.log"), stderr = file.path(out, "patch.log")) == 0L)
  stopifnot(system2("make", shQuote(c("-C", build, "-j2", "bcftools", "PLUGINS_ENABLED=no")),
    stdout = file.path(out, "build.log"), stderr = file.path(out, "build.log")) == 0L)
  patched <- file.path(build, "bcftools")
  cases <- readRDS(source_paths[1L])$cases
  native <- readRDS(source_paths[2L])
  samples <- unique(vapply(native[[1L]], `[[`, "", "sample"))
  occupied <- unlist(lapply(cases, function(case) {
    paths <- native[[case$transcript]]
    stopifnot(length(paths) == 6L)
    vapply(seq_along(paths), function(i) {
      p <- paths[[i]]
      if (!length(p$contributors)) return(NA_character_)
      paste(case$chrom, p$sample, (i - 1L) %% 2L + 1L, sep = "/")
    }, "")
  }), use.names = FALSE)
  occupied <- occupied[!is.na(occupied)] # Reference lanes have no consequence rows.
  stopifnot(!anyDuplicated(occupied))
  observations <- list()
  for (sample in c(samples, "cohort")) {
    command <- c("csq", "-f", source_paths[3L], "-g", source_paths[4L], "-p", "a", "-Ot",
      if (sample != "cohort") c("-s", sample), source_paths[5L])
    file <- file.path(out, paste0(sample, ".tsv"))
    status <- system2(if (sample == "cohort") patched else unpatched, shQuote(command),
      stdout = file, stderr = paste0(file, ".log"))
    lines <- readLines(file)
    observations[[sample]] <- list(status = status, command = command,
      rows = lines[startsWith(lines, "CSQ\t")])
  }
  expected <- unlist(lapply(observations[samples], `[[`, "rows"), use.names = FALSE)
  actual <- observations$cohort$rows
  original_counts <- table(expected)
  patched_counts <- table(actual)
  keys <- union(names(original_counts), names(patched_counts))
  pairs <- data.frame(row = keys, unpatched_individual = as.integer(original_counts[keys]),
    patched_cohort = as.integer(patched_counts[keys]))
  pairs$unpatched_individual[is.na(pairs$unpatched_individual)] <- 0L
  pairs$patched_cohort[is.na(pairs$patched_cohort)] <- 0L
  pairs$equal <- pairs$unpatched_individual == pairs$patched_cohort
  # Retain full multiset differences; a set comparison cannot hide extra rows.
  write.csv(pairs, file.path(out, "pairs.csv"), row.names = FALSE)
  saveRDS(observations, file.path(out, "observations.rds"))
  carrier_keys <- function(rows) unique(vapply(strsplit(rows, "\t", fixed = TRUE),
    function(x) if (length(x) == 6L) paste(x[4L], x[2L], x[3L], sep = "/") else "invalid", ""))
  individual_keys <- carrier_keys(expected)
  cohort_keys <- carrier_keys(actual)
  summary <- data.frame(transcripts = length(cases), lanes = 6L * length(cases),
    occupied = length(occupied), individual_rows = length(expected), cohort_rows = length(actual),
    individual_failures = sum(vapply(observations[samples], function(x) x$status != 0L, TRUE)),
    cohort_exit_status = observations$cohort$status,
    missing_individual = sum(!occupied %in% individual_keys), extra_individual = sum(!individual_keys %in% occupied),
    missing_cohort = sum(!occupied %in% cohort_keys), extra_cohort = sum(!cohort_keys %in% occupied),
    differing_distinct_rows = sum(pairs$unpatched_individual == 0L | pairs$patched_cohort == 0L),
    differing_multiset_rows = sum(!pairs$equal),
    extra_cohort_copies = sum(pmax(0L, pairs$patched_cohort - pairs$unpatched_individual)),
    missing_cohort_copies = sum(pmax(0L, pairs$unpatched_individual - pairs$patched_cohort)))
  write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
  files <- c(original_paths, file.path(audit, "receipt.json"), code, patched, file.path(build, "csq.c"),
    file.path(out, c("patch.log", "build.log", "pairs.csv", "observations.rds", "summary.csv")),
    list.files(out, "\\.tsv(\\.log)?$", full.names = TRUE))
  jsonlite::write_json(list(source_revision = revision, source_binding = "diagnostic_unbound",
    scope = "proposed_bcftools_prefix_fix_vs_unpatched_individual_samples_not_vep_conformance",
    source_audit_receipt = file.path(audit, "receipt.json"),
    bcftools_source_revision = receipt$bcftools_source_revision,
    sha256 = as.list(vapply(files, duckvep_evidence_sha256, ""))), file.path(out, "receipt.json"),
    pretty = TRUE, auto_unbox = TRUE)
  print(summary)
  stopifnot(identical(hashes, vapply(original_paths, duckvep_evidence_sha256, "")),
    identical(code_hashes, vapply(code, duckvep_evidence_sha256, "")),
    identical(revision, duckvep_evidence_revision(root)), summary$individual_failures == 0L,
    summary$cohort_exit_status == 0L, summary$missing_individual == 0L, summary$extra_individual == 0L,
    summary$missing_cohort == 0L, summary$extra_cohort == 0L, all(pairs$equal))
}
main()
