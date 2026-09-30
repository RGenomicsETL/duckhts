# Compare complete VCF-count outputs from two extension builds, keyed by sample/site.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, all(file.exists(args)))
con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
quote_string <- function(value) as.character(DBI::dbQuoteString(con, value))
before <- quote_string(normalizePath(args[[1L]], mustWork = TRUE))
after <- quote_string(normalizePath(args[[2L]], mustWork = TRUE))
invisible(DBI::dbExecute(con, paste0(
  "CREATE TEMP VIEW baseline AS FROM read_parquet(", before, ")")))
invisible(DBI::dbExecute(con, paste0(
  "CREATE TEMP VIEW candidate AS FROM read_parquet(", after, ")")))
query <- paste0(
  "SELECT (SELECT count(*) FROM baseline) AS rows_before, ",
  "(SELECT count(*) FROM candidate) AS rows_after, ",
  "(SELECT count(*) FROM (SELECT sample_id, site_index FROM baseline ",
  "GROUP BY 1, 2 HAVING count(*) > 1)) AS duplicate_keys_before, ",
  "(SELECT count(*) FROM (SELECT sample_id, site_index FROM candidate ",
  "GROUP BY 1, 2 HAVING count(*) > 1)) AS duplicate_keys_after, ",
  "(SELECT count(*) FROM (SELECT * FROM baseline EXCEPT ALL ",
  "SELECT * FROM candidate)) AS lost, ",
  "(SELECT count(*) FROM (SELECT * FROM candidate EXCEPT ALL ",
  "SELECT * FROM baseline)) AS gained")
result <- DBI::dbGetQuery(con, query)
DBI::dbDisconnect(con, shutdown = TRUE)
write.table(result, stdout(), sep = "\t", row.names = FALSE, quote = FALSE)
stopifnot(result$rows_before > 0, result$rows_before == result$rows_after,
          result$duplicate_keys_before == 0, result$duplicate_keys_after == 0,
          result$lost == 0, result$gained == 0)
