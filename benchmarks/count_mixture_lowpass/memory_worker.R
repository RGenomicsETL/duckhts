suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(duckdb))
args <- commandArgs(TRUE)
mode <- args[[1L]]
size <- as.integer(args[[2L]])
threads <- as.integer(args[[3L]])
output <- args[[4L]]
extension <- normalizePath(args[[5L]])
groups <- if (length(args) >= 6L) as.integer(args[[6L]]) else 1L
samples <- if (groups == 1L) 1L else 2L
candidates <- if (groups == 4L) 2L else 1L
temporary <- tempfile("lowpass-memory-")
dir.create(temporary)
con <- dbConnect(duckdb(shared_home = FALSE, config = list(
  allow_unsigned_extensions = "true", threads = as.character(threads), temp_directory = temporary)))
q <- function(x) as.character(dbQuoteString(con, x))
dbExecute(con, paste("LOAD", q(extension)))
dbExecute(con, "SET memory_limit='3GB'")
if (mode != "empty") {
  dbExecute(con, sprintf(paste(
    "CREATE VIEW counts AS SELECT 's'||s.k AS sample_id, 'c'||c.j AS candidate, '20' AS chrom,",
    "i+1 AS pos, (i%%3)::BIGINT AS ref_count, (i%%2)::BIGINT AS alt_count, 0::BIGINT AS other_count,",
    "(.05+.9*(i%%997)/997.0)::DOUBLE AS sample_af,",
    "(.05+.9*(i%%991)/991.0)::DOUBLE AS candidate_af FROM range(%d) t(i)",
    "CROSS JOIN range(%d) s(k) CROSS JOIN range(%d) c(j)"), size, samples, candidates))
}
profile <- paste0(output, ".json")
dbExecute(con, "PRAGMA enable_profiling='json'")
dbExecute(con, paste("SET profiling_output=", q(profile)))
query <- if (mode == "empty") "SELECT 1" else paste(
  "SELECT * FROM duckhts_count_mixture_lowpass('counts', .01, .02, 0, block_bases := 1000000)")
start <- proc.time()[["elapsed"]]
result <- dbGetQuery(con, query)
elapsed <- proc.time()[["elapsed"]] - start
dbExecute(con, "PRAGMA disable_profiling")
p <- jsonlite::fromJSON(profile, simplifyVector = FALSE)
saveRDS(list(result = result, seconds = elapsed, buffer_bytes = p$system_peak_buffer_memory,
             temp_bytes = p$system_peak_temp_dir_size, sites = size, groups = groups,
             samples = samples, candidates = candidates, output_rows = nrow(result),
             input_rows = size*groups, decoded_bytes = size*groups*56, threads = threads), output)
dbDisconnect(con, shutdown = TRUE)
