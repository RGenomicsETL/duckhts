cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
reference_path <- duckhtsbench::duckhts_bench_artifact_path("ancestry_reference_grch38_parquet")
connection <- Rduckhts::rduckhts_connect(
  extension_path = normalizePath("build/release/duckhts.duckdb_extension"),
  config = list(memory_limit = "12GB", threads = "4",
                temp_directory = file.path(cache, "duckdb-tmp")))
quote_string <- function(value) as.character(DBI::dbQuoteString(connection, value))
for (chromosome in 1:22) {
  long_output <- file.path(cache, sprintf("chr%d.reference_long.parquet", chromosome))
  af_output <- file.path(cache, sprintf("chr%d.af_reference_sites.parquet", chromosome))
  if (chromosome == 20L && file.exists(long_output) && file.exists(af_output)) {
    cat("reference chromosome 20 reused\n")
    next
  }
  unlink(c(long_output, af_output))
  DBI::dbExecute(connection, "DROP VIEW IF EXISTS ref_current")
  DBI::dbExecute(connection, sprintf(paste0(
    "CREATE TEMP VIEW ref_current AS SELECT * FROM read_parquet(%s) ",
    "WHERE try_cast(regexp_replace(CAST(chromosome AS VARCHAR), '^chr', '') AS INTEGER) = %d"),
    quote_string(reference_path), chromosome))
  columns <- DBI::dbGetQuery(connection, "DESCRIBE ref_current")$column_name
  pc_columns <- paste0("PC", seq_len(16L))
  group_columns <- setdiff(columns, c(
    "chromosome", "position", "allele_a", "allele_b", pc_columns))
  if (length(group_columns) != 21L) {
    stop("unexpected ancestry group columns", call. = FALSE)
  }
  unpivot_columns <- paste(as.character(DBI::dbQuoteIdentifier(
    connection, group_columns)), collapse = ", ")
  long_query <- sprintf(paste0(
    "UNPIVOT ref_current ON %s INTO NAME group_id VALUE frequency"),
    unpivot_columns)
  DBI::dbExecute(connection, sprintf(paste0(
    "COPY (SELECT chromosome::VARCHAR AS chromosome, position::BIGINT AS position, ",
    "allele_a, allele_b, group_id, frequency::DOUBLE AS frequency FROM (%s)) ",
    "TO %s (FORMAT PARQUET, COMPRESSION ZSTD)"), long_query, quote_string(long_output)))
  bcf_path <- file.path(cache, sprintf("chr%d.children.bcf", chromosome))
  DBI::dbExecute(connection, sprintf(paste0(
    "COPY (SELECT v.CHROM AS chrom, v.POS AS pos, v.REF AS ref, v.ALT[1] AS alt, ",
    "v.INFO_AF, v.INFO_AF_AFR, v.INFO_AF_AMR, v.INFO_AF_EAS, v.INFO_AF_EUR ",
    "FROM read_bcf(%s) v JOIN (SELECT DISTINCT position, allele_a, allele_b ",
    "FROM ref_current) r ON v.POS = r.position AND ",
    "((v.REF = r.allele_a AND v.ALT[1] = r.allele_b) OR ",
    "(v.REF = r.allele_b AND v.ALT[1] = r.allele_a)) ",
    "WHERE len(v.ALT) = 1 AND NOT ((r.allele_a IN ('A','T') AND ",
    "r.allele_b IN ('A','T')) OR (r.allele_a IN ('C','G') AND ",
    "r.allele_b IN ('C','G'))) AND v.INFO_AF IS NOT NULL ",
    "AND v.INFO_AF_AFR IS NOT NULL AND v.INFO_AF_AMR IS NOT NULL ",
    "AND v.INFO_AF_EAS IS NOT NULL AND v.INFO_AF_EUR IS NOT NULL) ",
    "TO %s (FORMAT PARQUET, COMPRESSION ZSTD)"),
    quote_string(bcf_path), quote_string(af_output)))
  long_rows <- DBI::dbGetQuery(connection, sprintf(
    "SELECT count(*) AS n FROM read_parquet(%s)", quote_string(long_output)))$n
  af_sites <- DBI::dbGetQuery(connection, sprintf(
    "SELECT count(*) AS n FROM read_parquet(%s)", quote_string(af_output)))$n
  cat("reference chromosome", chromosome, "long rows", long_rows,
      "AF sites", af_sites, "\n")
  flush.console()
}
DBI::dbDisconnect(connection, shutdown = TRUE)
