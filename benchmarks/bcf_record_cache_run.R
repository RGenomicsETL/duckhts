# One extension per fresh process; invoked by benchmark_bcf_record_cache.Rmd.
library(DBI)
library(duckdb)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 7L)
extension <- normalizePath(args[1], mustWork = TRUE)
input <- normalizePath(args[2], mustWork = TRUE)
reader <- match.arg(args[3], c("read_bcf", "read_bcf_v2"))
tidy <- match.arg(args[4], c("wide", "tidy")) == "tidy"
expected_rows <- as.numeric(args[5])
snapshot <- args[6]
output <- args[7]
con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = "true")))
dbExecute(con, paste("LOAD", dbQuoteString(con, extension)))
dbExecute(con, "SET threads=1")
dbExecute(con, "SET memory_limit='16GB'")
sql <- sprintf(paste0("CREATE TABLE result AS SELECT * FROM %s(%s,",
  "tidy_format:=%s,scan_mode:='sequential',decompression_threads:=0)"),
  reader, dbQuoteString(con, input), if (tidy) "true" else "false")
timing <- system.time(dbExecute(con, sql))
hwm <- readLines("/proc/self/status")
rss <- as.numeric(sub("^VmHWM:\\s+([0-9]+).*", "\\1", hwm[grepl("^VmHWM:", hwm)]))
rows <- as.numeric(dbGetQuery(con, "SELECT count(*) AS n FROM result")$n)
stopifnot(rows == expected_rows)
# Every materialized field contributes to the repetition checksum. The Rmd
# additionally compares complete first-run snapshots with EXCEPT ALL both ways.
fingerprint <- dbGetQuery(con, paste(
  "SELECT sum(hash(r)::HUGEINT)::VARCHAR AS checksum FROM result r"))$checksum
if (nzchar(snapshot)) dbExecute(con, sprintf(
  "COPY result TO %s (FORMAT PARQUET, COMPRESSION ZSTD)", dbQuoteString(con, snapshot)))
write.table(data.frame(rows = rows, columns = ncol(dbGetQuery(con, "SELECT * FROM result LIMIT 0")),
  checksum = fingerprint, elapsed = unname(timing["elapsed"]),
  cpu = sum(timing[c("user.self", "sys.self")]), peak_rss_kib = rss),
  output, sep = "\t", row.names = FALSE, quote = FALSE)
dbDisconnect(con, shutdown = TRUE)
