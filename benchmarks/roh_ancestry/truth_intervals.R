library(DBI)
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:", config = list(threads = 2L))
path <- file.path(cache, "chr20.truth_windows.parquet")
query <- sprintf(paste0(
  "WITH eligible AS (SELECT sample_id, window_id FROM read_parquet(%s) WHERE heterozygotes <= 1), ",
  "numbered AS (SELECT sample_id, window_id, window_id - row_number() OVER ",
  "(PARTITION BY sample_id ORDER BY window_id) AS island FROM eligible), ",
  "runs AS (SELECT sample_id, min(window_id) AS start_window, max(window_id) AS end_window, ",
  "count(*) AS windows FROM numbered GROUP BY sample_id, island HAVING count(*) >= 10) ",
  "SELECT sample_id, start_window, end_window, windows, start_window * 100000 + 1 AS start, ",
  "least((end_window + 1) * 100000, 64444167) AS end ",
  "FROM runs ORDER BY sample_id, start_window"), as.character(dbQuoteString(con, path)))
runs <- dbGetQuery(con, query)
write.csv(runs, file.path(cache, "chr20.truth_intervals.csv"), row.names = FALSE)
print(table(runs$windows))
cat(nrow(runs), "truth intervals across", length(unique(runs$sample_id)), "children\n")
dbDisconnect(con, shutdown = TRUE)
