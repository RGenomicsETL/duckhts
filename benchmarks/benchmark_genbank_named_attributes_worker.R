# One fresh process, one workload: load the extension, run the query once, and
# print elapsed seconds, peak resident memory, and the result row.
# Usage: Rscript benchmark_genbank_named_attributes_worker.R <extension> <sql>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

peak_rss_mib <- function() {
  status <- readLines("/proc/self/status")
  hwm <- grep("^VmHWM:", status, value = TRUE)
  as.numeric(gsub("[^0-9]", "", hwm)) / 1024
}

con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = "true")))
invisible(dbExecute(con, "SET threads = 1"))
invisible(dbExecute(con, paste("LOAD", dbQuoteString(con, args[[1L]]))))
before <- peak_rss_mib()
start <- proc.time()[["elapsed"]]
result <- dbGetQuery(con, args[[2L]])
elapsed <- proc.time()[["elapsed"]] - start
cat(sprintf("%.6f\t%.1f\t%.1f\t%s\t%s\n", elapsed, before, peak_rss_mib(),
            result$n[[1L]], result$metric[[1L]]))
dbDisconnect(con, shutdown = TRUE)
