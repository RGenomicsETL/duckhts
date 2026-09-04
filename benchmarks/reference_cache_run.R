# One fresh process per revision/function/format/thread count. Invoked by Rmd.
library(DBI)
library(duckdb)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L)
operator <- match.arg(args[1], c("norm", "munge"))
threads <- as.integer(args[2])
reference <- normalizePath(args[3], mustWork = TRUE)
rows <- as.integer(args[4])
runs <- as.integer(args[5])
output <- args[6]
extension <- normalizePath(Sys.getenv("REFERENCE_CACHE_EXTENSION"), mustWork = TRUE)
probe <- Sys.getenv("REFERENCE_CACHE_METRICS", "")
if (nzchar(probe)) dyn.load(probe)
metrics <- function() {
  if (nzchar(probe)) .C("reference_cache_metrics", values = double(4))$values else rep(NA_real_, 4)
}
con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = "true")))
dbExecute(con, paste("LOAD", dbQuoteString(con, extension)))
dbExecute(con, paste("SET threads =", threads))
path <- dbQuoteString(con, reference)
if (!file.exists(paste0(reference, ".fai")))
  stopifnot(dbGetQuery(con, paste0("SELECT * FROM fasta_index(", path, ")"))$success[1])
# Synthetic SNVs at consecutive A/C/G/T positions on the registered full GRCh38
# reference. Staging and table construction are outside the measured query.
bases <- strsplit(toupper(dbGetQuery(con, sprintf(
  "SELECT SEQUENCE FROM read_fasta(%s, region := 'chr22:20000000-21999999')",
  path))$SEQUENCE), "", fixed = TRUE)[[1L]]
positions <- head(which(bases %in% c("A", "C", "G", "T")), rows)
stopifnot(length(positions) == rows)
ref <- bases[positions]
dbWriteTable(con, "inputs", data.frame(i = seq_len(rows) - 1L,
  pos = 19999999L + positions, ref = ref, alt = ifelse(ref == "A", "C", "A")))
rm(bases, positions, ref)
invisible(gc())
stopifnot(dbGetQuery(con, "SELECT count(*) AS n FROM inputs")$n == rows)
call <- if (operator == "norm") {
  sprintf("bcftools_norm_row('chr22',pos,ref,[alt],%s,NULL,NULL,NULL,NULL)", path)
} else {
  sprintf(paste0("bcftools_munge_row('chr22',pos,alt,ref,NULL,",
    "NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,",
    "%s,'IFFY','REF_MISMATCH',NULL,NULL,NULL)"), path)
}
predicate <- if (operator == "norm") {
  "value.pos_normed = pos AND value.ref_normed = ref AND value.alt_normed = [alt] AND value.norm_status = 'Unchanged' AND NOT value.normed"
} else {
  "value.pos = pos AND value.ref = ref AND value.alt = alt AND value.filter IS NULL AND NOT value.alleles_swapped"
}
sql <- sprintf(paste("SELECT count(*) AS n, count(*) FILTER (WHERE %s) AS correct,",
  "sum(i)::DOUBLE AS input_ids FROM (SELECT *, %s AS value FROM inputs)"), predicate, call)
measure <- function(run) {
  before <- metrics()
  timing <- system.time(result <- dbGetQuery(con, sql))
  after <- metrics()
  stopifnot(result$n == rows, result$correct == rows,
    result$input_ids == (rows - 1) * rows / 2)
  delta <- after - before
  if (nzchar(probe)) stopifnot(delta[1] > 0, delta[1] <= rows)
  hwm <- readLines("/proc/self/status")
  rss <- as.numeric(sub("^VmHWM:\\s+([0-9]+).*", "\\1", hwm[grepl("^VmHWM:", hwm)]))
  data.frame(run = run, rows = rows, reference_bases_requested = rows,
    output_rows = as.numeric(result$n), correct_rows = as.numeric(result$correct),
    elapsed = unname(timing["elapsed"]), cpu = sum(timing[c("user.self", "sys.self")]),
    peak_rss_kib = rss, fetches = delta[1], fetched_bases = delta[2],
    cache_hits = rows - delta[1], handle_loads = delta[3], handle_destroys = delta[4])
}
results <- do.call(rbind, lapply(seq_len(runs), measure))
write.table(results, output, sep = "\t", row.names = FALSE, quote = FALSE)
dbDisconnect(con, shutdown = TRUE)
