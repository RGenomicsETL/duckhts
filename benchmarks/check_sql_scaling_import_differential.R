# Compare complete site-import outputs from two builds on the same staged VCF.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, all(file.exists(args)))
con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
quote_string <- function(value) as.character(DBI::dbQuoteString(con, value))
invisible(DBI::dbExecute(con, paste0(
  "CREATE TEMP VIEW baseline AS FROM read_parquet(",
  quote_string(normalizePath(args[[1L]])), ")")))
invisible(DBI::dbExecute(con, paste0(
  "CREATE TEMP VIEW candidate AS FROM read_parquet(",
  quote_string(normalizePath(args[[2L]])), ")")))
query <- paste0(
  "SELECT (SELECT count(*) FROM baseline) AS rows_before, ",
  "(SELECT count(*) FROM candidate) AS rows_after, ",
  "(SELECT count(*) FROM (SELECT assembly, site_index FROM baseline ",
  "GROUP BY ALL HAVING count(*) > 1)) AS duplicate_keys_before, ",
  "(SELECT count(*) FROM (SELECT assembly, site_index FROM candidate ",
  "GROUP BY ALL HAVING count(*) > 1)) AS duplicate_keys_after, ",
  "(SELECT count(*) FROM baseline b FULL OUTER JOIN candidate c ",
  "USING (assembly, site_index) WHERE b.site_index IS NULL ",
  "OR c.site_index IS NULL OR b IS DISTINCT FROM c) AS keyed_disagreements, ",
  "(SELECT count(*) FROM (SELECT * FROM baseline EXCEPT ALL ",
  "SELECT * FROM candidate)) AS lost, ",
  "(SELECT count(*) FROM (SELECT * FROM candidate EXCEPT ALL ",
  "SELECT * FROM baseline)) AS gained")
result <- DBI::dbGetQuery(con, query)
DBI::dbDisconnect(con, shutdown = TRUE)
write.table(result, stdout(), sep = "\t", row.names = FALSE, quote = FALSE)
stopifnot(result$rows_before > 0, result$rows_before == result$rows_after,
          result$duplicate_keys_before == 0, result$duplicate_keys_after == 0,
          result$keyed_disagreements == 0, result$lost == 0, result$gained == 0)
