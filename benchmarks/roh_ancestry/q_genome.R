# Genome-wide ancestry proportions of every child for the ROH ancestry
# evaluation, from dosages at the matched non-palindromic reference sites of
# the run's autosomes, in batches of children.
# Usage, from the repository root: Rscript benchmarks/roh_ancestry/q_genome.R [BATCH_SIZE]
# A batch is reused only when its marker records the same SHA-256 for the
# extension, this script, the shared helpers, the installed Rduckhts R code
# (it holds the solver), the ancestry reference, the correction table, the
# children list and every children BCF of the run, with the same chromosome
# set, batch members and package versions.
source("benchmarks/roh_ancestry/genome_common.R")
args <- commandArgs(trailingOnly = TRUE)
batch_size <- if (length(args)) as.integer(args[[1L]]) else 16L
if (length(batch_size) != 1L || is.na(batch_size) || batch_size < 1L || batch_size > 32L) {
  stop("batch size must be in 1..32", call. = FALSE)
}
chromosomes <- genome_chromosomes()
reference_path <- duckhtsbench::duckhts_bench_artifact_path("ancestry_reference_grch38_parquet")
correction_path <- normalizePath("test/data/ancestry_correction_bigsnpr_1.12.21.tsv",
                                 mustWork = TRUE)
bcf_paths <- stats::setNames(vapply(chromosomes, genome_children_bcf, character(1L)),
                             paste0("bcf_", chromosomes))
extension <- normalizePath("build/release/duckhts.duckdb_extension", mustWork = TRUE)
common_identity <- c(genome_set_line(chromosomes), genome_identity(
  c(extension = extension,
    driver = normalizePath("benchmarks/roh_ancestry/q_genome.R", mustWork = TRUE),
    helpers = normalizePath("benchmarks/roh_ancestry/genome_common.R", mustWork = TRUE),
    solver = system.file("R", "Rduckhts.rdb", package = "Rduckhts"),
    reference = reference_path, correction = correction_path,
    children = genome_children_list, bcf_paths),
  c(batch_size = batch_size, duckdb_r = as.character(utils::packageVersion("duckdb")),
    rduckhts = as.character(utils::packageVersion("Rduckhts")))))
connection <- Rduckhts::rduckhts_connect(
  extension_path = extension,
  config = list(memory_limit = "12GB", threads = "4",
                temp_directory = file.path(genome_cache, "duckdb-tmp")))
quote_string <- function(value) as.character(DBI::dbQuoteString(connection, value))
DBI::dbExecute(connection, sprintf(
  "CREATE VIEW reference_all AS SELECT * FROM read_parquet(%s)",
  quote_string(reference_path)))
DBI::dbExecute(connection, sprintf(
  "CREATE VIEW ancestry_correction AS SELECT * FROM read_csv(%s, delim='\\t', header=true)",
  quote_string(correction_path)))
all_samples <- readLines(genome_children_list)
batches <- split(all_samples, ceiling(seq_along(all_samples) / batch_size))
batch_valid <- function(batch_result, samples) {
  nrow(batch_result) == length(samples) * 21L &&
    length(unique(batch_result$group_id)) == 21L &&
    setequal(unique(batch_result$sample_id), samples) &&
    all(batch_result$status == "ok") && all(batch_result$cor_pred >= 0.4)
}
results <- vector("list", length(batches))
for (batch_index in seq_along(batches)) {
  samples <- batches[[batch_index]]
  batch_output <- file.path(genome_cache, sprintf("autosome.q_batch%02d_%03d.csv",
                                                  batch_size, batch_index))
  identity <- c(common_identity, paste0("batch_samples=", paste(samples, collapse = ",")))
  if (genome_current(batch_output, identity)) {
    batch_result <- utils::read.csv(batch_output, stringsAsFactors = FALSE)
    if (!batch_valid(batch_result, samples)) {
      stop("cached q batch failed validation", call. = FALSE)
    }
    results[[batch_index]] <- batch_result
    cat("reused batch", batch_index, "of", length(batches), "children", length(samples), "\n")
    next
  }
  unlink(c(batch_output, genome_marker(batch_output)))
  sample_selector <- quote_string(paste(samples, collapse = ","))
  reads <- vapply(bcf_paths, function(path) sprintf(paste0(
    "SELECT SAMPLE_ID, CHROM, POS, REF, ALT, FORMAT_GT ",
    "FROM read_bcf(%s, tidy_format := true, samples := %s)"),
    quote_string(path), sample_selector), character(1L))
  input_sql <- paste(reads, collapse = " UNION ALL BY NAME ")
  DBI::dbExecute(connection, "DROP TABLE IF EXISTS child_dosage")
  sql <- sprintf(paste0(
    "CREATE TEMP TABLE child_dosage AS ",
    "SELECT b.SAMPLE_ID AS sample_id, b.CHROM AS chromosome, ",
    "b.POS AS position, b.REF AS allele_a, b.ALT[1] AS allele_b, ",
    "CASE WHEN regexp_matches(b.FORMAT_GT, '^[01][/|][01]$') THEN ",
    "(length(b.FORMAT_GT) - length(replace(b.FORMAT_GT, '1', '')))::DOUBLE END AS dosage ",
    "FROM (%s) b JOIN reference_all r ON ",
    "try_cast(regexp_replace(b.CHROM, '^chr', '') AS INTEGER) = r.chromosome ",
    "AND b.POS = r.position AND ((b.REF = r.allele_a AND b.ALT[1] = r.allele_b) ",
    "OR (b.REF = r.allele_b AND b.ALT[1] = r.allele_a)) ",
    "WHERE len(b.ALT) = 1 AND NOT ((r.allele_a IN ('A','T') AND ",
    "r.allele_b IN ('A','T')) OR (r.allele_a IN ('C','G') AND r.allele_b IN ('C','G')))"),
    input_sql)
  DBI::dbExecute(connection, sql)
  dosage_rows <- DBI::dbGetQuery(connection,
    "SELECT count(*) AS n, count(DISTINCT sample_id) AS samples FROM child_dosage")
  if (dosage_rows$samples != length(samples)) {
    stop("dosage batch omitted a child", call. = FALSE)
  }
  cat("batch", batch_index, "of", length(batches), "children", length(samples),
      "dosage_rows", dosage_rows$n, "\n")
  flush.console()
  batch_result <- Rduckhts::rduckhts_ancestry_proportions(
    connection, "child_dosage", "reference_all", "ancestry_correction",
    input_kind = "dosage", min_cor = 0.4)
  if (!batch_valid(batch_result, samples)) {
    stop("genome-wide q batch validation failed", call. = FALSE)
  }
  utils::write.csv(batch_result, batch_output, row.names = FALSE)
  genome_record(batch_output, identity)
  results[[batch_index]] <- batch_result
  DBI::dbExecute(connection, "DROP TABLE child_dosage")
  invisible(gc())
}
proportions <- do.call(rbind, results)
if (!setequal(unique(proportions$sample_id), all_samples) ||
    any(proportions$status != "ok") || any(proportions$cor_pred < 0.4)) {
  stop("genome-wide q validation failed", call. = FALSE)
}
q_path <- file.path(genome_cache, "autosome.q.csv")
utils::write.csv(proportions, q_path, row.names = FALSE)
genome_record(q_path, common_identity)
cat("q children", length(unique(proportions$sample_id)),
    "groups", length(unique(proportions$group_id)), "chromosomes", length(chromosomes),
    "cor_pred range", range(proportions$cor_pred), "\n")
DBI::dbDisconnect(connection, shutdown = TRUE)
