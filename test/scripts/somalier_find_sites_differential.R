#!/usr/bin/env Rscript

# Somalier v0.3.4 vs a caller-connection DuckHTS selection on a retained
# three-sample 1000 Genomes chromosome-1 VCF slice.
library(DBI)
library(Rduckhts)
library(duckhtsbench)

site_disagreements <- function(left, right, slice) {
  keys <- c("region", "position", "source_ref", "source_alt")
  left <- left[, keys, drop = FALSE]
  right <- right[, keys, drop = FALSE]
  if (anyDuplicated(left) || anyDuplicated(right)) {
    stop("differential keys must be unique in both outputs: ", slice)
  }
  flag <- right
  flag$found <- rep(TRUE, nrow(flag))
  left_only <- merge(left, flag, by = keys, all.x = TRUE)
  left_only <- left_only[is.na(left_only$found), keys, drop = FALSE]
  flag <- left
  flag$found <- rep(TRUE, nrow(flag))
  right_only <- merge(right, flag, by = keys, all.x = TRUE)
  right_only <- right_only[is.na(right_only$found), keys, drop = FALSE]
  rows <- rbind(left_only, right_only)
  rows$origin <- c(rep("DuckHTS only", nrow(left_only)),
    rep("Somalier only", nrow(right_only)))
  rows$slice <- rep(slice, nrow(rows))
  rows
}

run_differential <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  inputs <- duckhts_bench_stage_repository_fixtures(root, "somalier-find-sites")
  input <- inputs[["somalier_find_sites_1000g"]]
  executable <- duckhtsbench:::duckhts_bench_stage_somalier_v034()
  output <- tempfile(fileext = ".vcf.gz")
  on.exit(unlink(output), add = TRUE)
  log <- suppressWarnings(system2(executable, c(
    "find-sites", "--min-AN", "6", "--snp-dist", "100",
    "--output-vcf", output, input
  ), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(log, "status")) && attr(log, "status") != 0L) {
    stop("Somalier find-sites failed: ", paste(log, collapse = "\n"))
  }
  candidate_message <- grep("^[0-9]+ candidate variants$", log, value = TRUE)
  if (length(candidate_message) != 1L) stop("Somalier candidate count missing")
  upstream_candidates <- as.integer(sub(" candidate variants$", "", candidate_message))

  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  dbExecute(con, "SET threads = 1")
  dbExecute(con, sprintf(
    "CREATE TEMP VIEW input_variants AS SELECT * FROM read_bcf(%s, samples := '', scan_mode := 'sequential')",
    as.character(dbQuoteString(con, input))
  ))
  parameters <- list(
    con = con, source = "input_variants", assembly = "GRCh37",
    min_af = 0.15, min_an = 6, af_field = "AF", an_field = "AN",
    snp_dist = 100, target_af = 0.48,
    intervals = list(include = NULL, exclude = NULL, gnotate = NULL),
    mode = "somalier_v0.3.4", max_autosomal = 65535,
    max_x = 10001, max_y = 5001
  )
  counts <- dbGetQuery(con, do.call(Rduckhts:::.somalier_find_sites_query,
    c(parameters, list(diagnostics = TRUE))))
  times <- numeric(5L)
  for (i in seq_along(times)) {
    times[[i]] <- system.time({
      selected <- rduckhts_somalier_find_sites(con, source_table = "input_variants",
        assembly = "GRCh37", min_an = 6, snp_dist = 100)
    })[["elapsed"]]
  }
  upstream <- dbGetQuery(con, sprintf(
    "SELECT CHROM AS region, CAST(POS AS UBIGINT) AS position, REF AS source_ref, ALT[1] AS source_alt FROM read_bcf(%s, samples := '', scan_mode := 'sequential')",
    as.character(dbQuoteString(con, output))
  ))
  disagreements <- site_disagreements(selected, upstream, "1000G")
  first_candidate_mismatch <- counts$records[counts$gate == "gated"] !=
    upstream_candidates
  second <- inputs[["somalier_find_sites_gnomad"]]
  second_output <- tempfile(fileext = ".vcf.gz")
  on.exit(unlink(second_output), add = TRUE)
  second_log <- suppressWarnings(system2(executable, c(
    "find-sites", "--min-AN", "100", "--snp-dist", "100",
    "--output-vcf", second_output, second
  ), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(second_log, "status")) && attr(second_log, "status") != 0L) {
    stop("Somalier gnomAD find-sites failed: ", paste(second_log, collapse = "\n"))
  }
  second_message <- grep("^[0-9]+ candidate variants$", second_log, value = TRUE)
  if (length(second_message) != 1L) stop("Somalier gnomAD candidate count missing")
  second_candidates <- as.integer(sub(" candidate variants$", "", second_message))
  dbExecute(con, sprintf(
    "CREATE TEMP VIEW gnomad_variants AS SELECT * FROM read_bcf(%s, samples := '', scan_mode := 'sequential')",
    as.character(dbQuoteString(con, second))
  ))
  second_parameters <- parameters
  second_parameters$source <- "gnomad_variants"
  second_parameters$min_an <- 100L
  second_counts <- dbGetQuery(con, do.call(Rduckhts:::.somalier_find_sites_query,
    c(second_parameters, list(diagnostics = TRUE))))
  second_selected <- rduckhts_somalier_find_sites(con,
    source_table = "gnomad_variants", assembly = "GRCh37", min_an = 100,
    snp_dist = 100)
  second_upstream <- dbGetQuery(con, sprintf(
    "SELECT CHROM AS region, CAST(POS AS UBIGINT) AS position, REF AS source_ref, ALT[1] AS source_alt FROM read_bcf(%s, samples := '', scan_mode := 'sequential')",
    as.character(dbQuoteString(con, second_output))
  ))
  disagreements <- rbind(disagreements,
    site_disagreements(second_selected, second_upstream, "gnomAD"))
  second_candidate_mismatch <- second_counts$records[
    second_counts$gate == "gated"] != second_candidates
  dir.create(file.path(root, "benchmarks"), showWarnings = FALSE)
  utils::write.table(second_counts, file.path(root,
    "benchmarks/somalier_find_sites_gnomad_gates.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE)
  utils::write.table(counts, file.path(root, "benchmarks/somalier_find_sites_gates.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(disagreements, file.path(root,
    "benchmarks/somalier_find_sites_disagreements.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE)
  observation <- data.frame(revision = system2("git", "rev-parse HEAD", stdout = TRUE),
    source_records = counts$records[counts$gate == "src"],
    upstream_candidates = upstream_candidates, duckhts_candidates = counts$records[counts$gate == "gated"],
    upstream_sites = nrow(upstream), duckhts_sites = nrow(selected),
    threads = 1L, repetitions = length(times),
    median_seconds = stats::median(times), min_seconds = min(times))
  utils::write.table(observation, file.path(root,
    "benchmarks/somalier_find_sites_observations.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE)
  print(observation)
  print(counts)
  print(disagreements)
  print(data.frame(slice = "gnomAD", source_records = second_counts$records[
    second_counts$gate == "src"], upstream_candidates = second_candidates,
  duckhts_candidates = second_counts$records[second_counts$gate == "gated"],
  upstream_sites = nrow(second_upstream), duckhts_sites = nrow(second_selected)))
  if (nrow(disagreements) || first_candidate_mismatch || second_candidate_mismatch) {
    stop("Somalier differential has disagreements; inspect retained TSV evidence")
  }
}

root <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[[1L]] else getwd()
run_differential(root)
