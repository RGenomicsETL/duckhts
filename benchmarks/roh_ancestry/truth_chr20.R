# The installed Rduckhts; load_all() would compile r/Rduckhts in place.
library(Rduckhts)
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
path <- file.path(cache, "chr20.children.bcf")
con <- rduckhts_connect(extension_path = normalizePath("build/release/duckhts.duckdb_extension"),
                        config = list(memory_limit = "12GB", threads = "4",
                                      temp_directory = file.path(cache, "duckdb-tmp")))
DBI::dbExecute(con, sprintf(
  "CREATE TEMP TABLE observed AS SELECT SAMPLE_ID AS sample_id, floor((POS - 1) / 100000)::INTEGER AS window_id, count(*) FILTER (WHERE FORMAT_GT IN ('0/1','1/0','0|1','1|0'))::INTEGER AS heterozygotes FROM read_bcf(%s, tidy_format := true) GROUP BY SAMPLE_ID, window_id",
  as.character(DBI::dbQuoteString(con, path))))
children <- readLines(file.path(cache, "chr20.children.present.txt"))
length_bp <- 64444167L
windows <- data.frame(window_id = seq.int(0L, (length_bp - 1L) %/% 100000L))
DBI::dbWriteTable(con, "children", data.frame(sample_id = children), temporary = TRUE)
DBI::dbWriteTable(con, "windows", windows, temporary = TRUE)
DBI::dbExecute(con, sprintf(
  "COPY (SELECT c.sample_id, w.window_id, coalesce(o.heterozygotes, 0)::INTEGER AS heterozygotes FROM children c CROSS JOIN windows w LEFT JOIN observed o USING(sample_id, window_id) ORDER BY c.sample_id, w.window_id) TO %s (FORMAT PARQUET, COMPRESSION ZSTD)",
  as.character(DBI::dbQuoteString(con, file.path(cache, "chr20.truth_windows.parquet")))))
summary <- DBI::dbGetQuery(con, sprintf(
  "SELECT heterozygotes, count(*) AS child_windows FROM read_parquet(%s) GROUP BY heterozygotes ORDER BY heterozygotes",
  as.character(DBI::dbQuoteString(con, file.path(cache, "chr20.truth_windows.parquet")))))
utils::write.csv(summary, file.path(cache, "chr20.truth_distribution.csv"), row.names = FALSE)
print(summary)
truth <- DBI::dbGetQuery(con, sprintf(
  "WITH eligible AS (SELECT sample_id, window_id FROM read_parquet(%s) WHERE heterozygotes <= 1), numbered AS (SELECT sample_id, window_id, window_id - row_number() OVER (PARTITION BY sample_id ORDER BY window_id) AS island FROM eligible), runs AS (SELECT sample_id, min(window_id) AS start_window, max(window_id) AS end_window, count(*) AS windows FROM numbered GROUP BY sample_id, island HAVING count(*) >= 10) SELECT sample_id, count(*) AS truth_runs, sum(windows) AS truth_windows FROM runs GROUP BY sample_id ORDER BY sample_id",
  as.character(DBI::dbQuoteString(con, file.path(cache, "chr20.truth_windows.parquet")))))
utils::write.csv(truth, file.path(cache, "chr20.truth_runs.csv"), row.names = FALSE)
print(head(truth))
DBI::dbDisconnect(con, shutdown = TRUE)
