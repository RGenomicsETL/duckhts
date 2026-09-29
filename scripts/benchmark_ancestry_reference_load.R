#!/usr/bin/env Rscript
# Compare full numeric-column scans of the registered CSVs and derived Parquet.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: benchmark_ancestry_reference_load.R ref_freqs.csv.gz projection.csv.gz reference.parquet threads")
}
threads <- as.integer(args[[4L]])
if (!threads %in% c(1L, 4L)) stop("threads must be 1 or 4")
library(DBI)
con <- dbConnect(duckdb::duckdb(shared_home = FALSE), dbdir = ":memory:")
dbExecute(con, paste0("SET threads=", threads))
quote_str <- function(path) as.character(dbQuoteString(con, normalizePath(path, mustWork = TRUE)))
quote_id <- function(name) as.character(dbQuoteIdentifier(con, name))
paths <- vapply(args[1:3], quote_str, character(1L))
dbExecute(con, paste0("CREATE VIEW frequencies AS SELECT * FROM read_csv_auto(", paths[[1L]], ")"))
dbExecute(con, paste0("CREATE VIEW loadings AS SELECT * FROM read_csv_auto(", paths[[2L]], ")"))
dbExecute(con, paste0("CREATE VIEW wide AS SELECT * FROM read_parquet(", paths[[3L]], ")"))
measure <- function(view, columns) {
  columns <- paste0("sum(", quote_id(columns), ")", collapse = ", ")
  start <- proc.time()[["elapsed"]]
  values <- dbGetQuery(con, paste0("SELECT count(*) AS n, ", columns, " FROM ", view))
  list(seconds = proc.time()[["elapsed"]] - start, n = values$n,
       sums = unname(unlist(values[-1L])))
}
groups <- setdiff(dbListFields(con, "frequencies"), c("chr", "pos", "rsid", "a0", "a1"))
pcs <- paste0("PC", seq_len(16L))
frequency <- measure("frequencies", groups)
loading <- measure("loadings", pcs)
parquet <- measure("wide", c(groups, pcs))
stopifnot(frequency$n == parquet$n, loading$n == parquet$n,
          max(abs(c(frequency$sums, loading$sums) - parquet$sums)) < 1e-5)
result <- data.frame(threads = threads, reference_rows = parquet$n,
                     numeric_columns = length(groups) + length(pcs),
                     csv_full_column_seconds = round(frequency$seconds + loading$seconds, 3L),
                     parquet_full_column_seconds = round(parquet$seconds, 3L))
write.table(result, stdout(), sep = "\t", row.names = FALSE, quote = FALSE)
dbDisconnect(con, shutdown = TRUE)
