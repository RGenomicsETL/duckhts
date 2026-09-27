#!/usr/bin/env Rscript

# Whole-chromosome comparison against the pinned Somalier v0.3.4 executable.
library(DBI)
library(Rduckhts)
library(duckhtsbench)

args <- commandArgs(TRUE)
if (length(args) != 3L || !args[[1L]] %in% c("upstream", "select", "compare") ||
    !args[[3L]] %in% c("chr22", "chrX")) {
  stop("usage: somalier_find_sites_chromosome.R upstream|select|compare REPO_ROOT chr22|chrX")
}
root <- normalizePath(args[[2L]], mustWork = TRUE)
chromosome <- args[[3L]]
prefix <- paste0("somalier_find_sites_", tolower(chromosome), "_")
id <- paste0(prefix, "source")
input <- duckhts_bench_fetch(id)
cache <- dirname(input)
upstream_file <- file.path(cache, paste0(chromosome, ".upstream.vcf.gz"))
selected_file <- file.path(cache, paste0(chromosome, ".duckhts.tsv"))

if (args[[1L]] == "upstream") {
  executable <- duckhtsbench:::duckhts_bench_stage_somalier_v034()
  log_file <- file.path(cache, paste0(chromosome, ".upstream.log"))
  time_file <- file.path(cache, paste0(chromosome, ".upstream.time"))
  status <- system2("/usr/bin/time", c(
    "-v", "-o", shQuote(time_file), shQuote(executable), "find-sites",
    "--min-AN", "6000", "--output-vcf", shQuote(upstream_file),
    shQuote(input)
  ), stdout = log_file, stderr = log_file)
  if (status != 0L) stop("Somalier failed; inspect ", log_file)

  tie_input <- file.path(cache, "tie-input.vcf")
  tie_output <- file.path(cache, "tie-output.vcf.gz")
  tie_log <- file.path(cache, "tie-probe.log")
  writeLines(c(
    "##fileformat=VCFv4.2", "##contig=<ID=chr1>",
    '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">',
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t200\t.\tA\tG\t.\tPASS\tAF=0.48;AN=120000",
    "chr1\t100\t.\tA\tG\t.\tPASS\tAF=0.48;AN=120000"
  ), tie_input)
  status <- system2(executable, c(
    "find-sites", "--snp-dist", "1000", "--output-vcf",
    shQuote(tie_output), shQuote(tie_input)
  ), stdout = tie_log, stderr = tie_log)
  if (status != 0L) stop("Somalier tie probe failed")
  tie_con <- rduckhts_connect()
  tie_result <- dbGetQuery(tie_con, sprintf(
    "SELECT POS FROM read_bcf(%s, samples := '', scan_mode := 'sequential')",
    as.character(dbQuoteString(tie_con, tie_output))
  ))
  dbDisconnect(tie_con, shutdown = TRUE)
  if (!identical(as.integer(tie_result$POS), 200L)) {
    stop("Somalier no longer retains input order on equal AF scores")
  }
}

if (args[[1L]] == "select") {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  dbExecute(con, "SET threads = 1")
  selected <- rduckhts_somalier_find_sites(
    con, source_vcf = input, assembly = "GRCh38", min_an = 6000
  )
  selected <- selected[, c("region", "position", "source_ref", "source_alt")]
  utils::write.table(selected, selected_file, sep = "\t", row.names = FALSE,
                     quote = FALSE, na = "")
}

if (args[[1L]] == "compare") {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  dbExecute(con, "SET threads = 1")
  dbExecute(con, sprintf(
    "CREATE TEMP VIEW population AS SELECT * FROM read_bcf(%s, samples := '', scan_mode := 'sequential')",
    as.character(dbQuoteString(con, input))
  ))
  parameters <- list(
    con = con, source = "population", assembly = "GRCh38",
    min_af = 0.15, min_an = 6000, af_field = "AF", an_field = "AN",
    snp_dist = 10000, target_af = 0.48,
    intervals = list(include = NULL, exclude = NULL, gnotate = NULL),
    mode = "somalier_v0.3.4", max_autosomal = 65535,
    max_x = 10001, max_y = 5001
  )
  counts <- dbGetQuery(con, do.call(Rduckhts:::.somalier_find_sites_query,
    c(parameters, list(diagnostics = TRUE))))
  selected <- utils::read.delim(selected_file, check.names = FALSE)
  upstream <- dbGetQuery(con, sprintf(
    "SELECT CHROM AS region, POS AS position, REF AS source_ref, ALT[1] AS source_alt FROM read_bcf(%s, samples := '', scan_mode := 'sequential')",
    as.character(dbQuoteString(con, upstream_file))
  ))
  keys <- c("region", "position", "source_ref", "source_alt")
  if (anyDuplicated(selected[, keys]) || anyDuplicated(upstream[, keys])) {
    stop("key collision in selected sites")
  }
  left <- merge(selected[, keys], upstream[, keys], by = keys)
  duck_membership <- merge(selected[, keys],
    transform(upstream[, keys], present = TRUE), by = keys, all.x = TRUE,
    sort = FALSE)
  upstream_membership <- merge(upstream[, keys],
    transform(selected[, keys], present = TRUE), by = keys, all.x = TRUE,
    sort = FALSE)
  duck_only <- duck_membership[is.na(duck_membership$present), keys]
  upstream_only <- upstream_membership[is.na(upstream_membership$present), keys]
  disagreements <- rbind(
    data.frame(duck_only, origin = rep("DuckHTS only", nrow(duck_only))),
    data.frame(upstream_only, origin = rep("Somalier only", nrow(upstream_only)))
  )
  lexical <- rduckhts_somalier_find_sites(
    con, source_table = "population", assembly = "GRCh38", min_an = 6000,
    tie_order = "lexical"
  )
  lexical_diff <- nrow(merge(selected[, keys], lexical[, keys], by = keys))
  log <- readLines(file.path(cache, paste0(chromosome, ".upstream.log")))
  candidates <- grep("^[0-9]+ candidate variants$", log, value = TRUE)
  if (length(candidates) != 1L) stop("missing upstream candidate denominator")
  candidates <- as.integer(sub(" candidate variants$", "", candidates))
  counts$somalier_records <- NA_integer_
  counts$somalier_records[counts$gate == "src"] <- counts$records[counts$gate == "src"]
  counts$somalier_records[counts$gate == "gated"] <- candidates
  counts$somalier_records[counts$gate == "selected"] <- nrow(upstream)
  time_value <- function(path, label) {
    line <- readLines(path)
    metric <- function(name) grep(name, line, value = TRUE, fixed = TRUE)[[1L]]
    wall <- sub(".*\\):[[:space:]]*", "", metric("Elapsed (wall clock) time"))
    parts <- rev(as.numeric(strsplit(wall, ":", fixed = TRUE)[[1L]]))
    rss <- sub("^[^:]*:[[:space:]]*", "", metric("Maximum resident set size"))
    data.frame(tool = label,
      wall_seconds = sum(parts * 60^(seq_along(parts) - 1L)),
      peak_rss_kib = as.integer(rss))
  }
  timing <- rbind(
    time_value(file.path(cache, paste0(chromosome, ".upstream.time")),
               "Somalier v0.3.4"),
    time_value(file.path(cache, paste0(chromosome, ".duckhts.time")),
               "DuckHTS")
  )
  dir.create(file.path(root, "benchmarks"), showWarnings = FALSE)
  write_tsv <- function(data, name) utils::write.table(data,
    file.path(root, "benchmarks", name), sep = "\t", row.names = FALSE,
    quote = FALSE)
  write_tsv(counts, paste0(prefix, "gates.tsv"))
  if (chromosome == "chrX") {
    par <- dbGetQuery(con, paste(
      "SELECT min(POS) AS first_position, max(POS) AS last_position,",
      "count_if(POS < 2781480) AS par1_records,",
      "count_if(POS > 154931045) AS par2_records FROM population"
    ))
    write_tsv(par, paste0(prefix, "par_records.tsv"))
  }
  write_tsv(utils::head(disagreements, 100L),
            paste0(prefix, "disagreements.tsv"))
  write_tsv(data.frame(class = c("DuckHTS only", "Somalier only"),
    records = c(nrow(duck_only), nrow(upstream_only))),
    paste0(prefix, "disagreement_summary.tsv"))
  summary <- data.frame(
    revision = system2("git", "rev-parse HEAD", stdout = TRUE),
    source_records = counts$records[counts$gate == "src"],
    upstream_candidates = candidates,
    duckhts_candidates = counts$records[counts$gate == "gated"],
    max_per_chromosome_candidates = counts$records[counts$gate == "neighbors"],
    upstream_sites = nrow(upstream), duckhts_sites = nrow(selected),
    common_sites = nrow(left), duckhts_only = nrow(duck_only),
    upstream_only = nrow(upstream_only),
    lexical_tie_dependent_sites = nrow(selected) - lexical_diff,
    threads_duckhts = 1L, threads_somalier_reader = 2L,
    threads_somalier_writer = 1L
  )
  write_tsv(summary, paste0(prefix, "observations.tsv"))
  write_tsv(timing, paste0(prefix, "timing.tsv"))
  print(summary)
  print(counts)
  print(timing)
  print(utils::head(disagreements, 10L))
  if (nrow(disagreements) || candidates != summary$duckhts_candidates) {
    stop(chromosome, " differential disagrees; inspect retained evidence")
  }
}
