args <- commandArgs(trailingOnly = TRUE)
batch_size <- if (length(args)) as.integer(args[[1L]]) else 16L
if (length(batch_size) != 1L || is.na(batch_size) || batch_size < 1L || batch_size > 32L) {
  stop("batch size must be in 1..32", call. = FALSE)
}
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
reference_path <- duckhtsbench::duckhts_bench_artifact_path("ancestry_reference_grch38_parquet")
correction_path <- normalizePath("test/data/ancestry_correction_bigsnpr_1.12.21.tsv")
bcf_paths <- file.path(cache, sprintf("chr%d.children.bcf", 1:22))
if (!all(file.exists(bcf_paths))) stop("one or more chromosome BCFs are missing", call. = FALSE)
connection <- Rduckhts::rduckhts_connect(
  extension_path = normalizePath("build/release/duckhts.duckdb_extension"),
  config = list(memory_limit = "12GB", threads = "4",
                temp_directory = file.path(cache, "duckdb-tmp")))
quote_string <- function(value) as.character(DBI::dbQuoteString(connection, value))
DBI::dbExecute(connection, sprintf(
  "CREATE VIEW reference_all AS SELECT * FROM read_parquet(%s)",
  quote_string(reference_path)))
DBI::dbExecute(connection, sprintf(
  "CREATE VIEW ancestry_correction AS SELECT * FROM read_csv(%s, delim='\\t', header=true)",
  quote_string(correction_path)))
all_samples <- readLines(file.path(cache, "chr20.children.present.txt"))
batches <- split(all_samples, ceiling(seq_along(all_samples) / batch_size))
results <- vector("list", length(batches))
for (batch_index in seq_along(batches)) {
  batch_output <- file.path(cache, sprintf("autosome.q_batch%02d_%03d.csv",
                                          batch_size, batch_index))
  if (file.exists(batch_output)) {
    batch_result <- utils::read.csv(batch_output, stringsAsFactors = FALSE)
    if (nrow(batch_result) != length(batches[[batch_index]]) * 21L ||
        length(unique(batch_result$group_id)) != 21L ||
        !setequal(unique(batch_result$sample_id), batches[[batch_index]]) ||
        any(batch_result$status != "ok") || any(batch_result$cor_pred < 0.4)) {
      stop("cached q batch failed validation", call. = FALSE)
    }
    results[[batch_index]] <- batch_result
    cat("reused batch", batch_index, "of", length(batches), "children",
        length(batches[[batch_index]]), "\n")
    next
  }
  sample_selector <- quote_string(paste(batches[[batch_index]], collapse = ","))
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
  if (dosage_rows$samples != length(batches[[batch_index]])) {
    stop("dosage batch omitted a child", call. = FALSE)
  }
  cat("batch", batch_index, "of", length(batches), "children",
      length(batches[[batch_index]]), "dosage_rows", dosage_rows$n, "\n")
  flush.console()
  results[[batch_index]] <- Rduckhts::rduckhts_ancestry_proportions(
    connection, "child_dosage", "reference_all", "ancestry_correction",
    input_kind = "dosage", min_cor = 0.4)
  batch_result <- results[[batch_index]]
  if (nrow(batch_result) != length(batches[[batch_index]]) * 21L ||
      length(unique(batch_result$group_id)) != 21L ||
      !setequal(unique(batch_result$sample_id), batches[[batch_index]]) ||
      any(batch_result$status != "ok") || any(batch_result$cor_pred < 0.4)) {
    stop("genome-wide q batch validation failed", call. = FALSE)
  }
  utils::write.csv(batch_result, batch_output, row.names = FALSE)
  DBI::dbExecute(connection, "DROP TABLE child_dosage")
  invisible(gc())
}
proportions <- do.call(rbind, results)
if (!setequal(unique(proportions$sample_id), all_samples) ||
    any(proportions$status != "ok") || any(proportions$cor_pred < 0.4)) {
  stop("genome-wide q validation failed", call. = FALSE)
}
utils::write.csv(proportions, file.path(cache, "autosome.q.csv"), row.names = FALSE)
utils::write.csv(proportions,
                 "benchmarks/results/roh-ancestry/autosome_q_by_child.csv",
                 row.names = FALSE)
cat("q children", length(unique(proportions$sample_id)),
    "groups", length(unique(proportions$group_id)),
    "cor_pred range", range(proportions$cor_pred), "\n")
DBI::dbDisconnect(connection, shutdown = TRUE)
