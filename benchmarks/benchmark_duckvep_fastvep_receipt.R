#!/usr/bin/env Rscript

# Produce one output-file receipt for the DuckVEP/FastVEP comparison. The
# fingerprint is independent of row order but remains sensitive to every
# projected column. Callers retain one row per timed or diagnostic run.

suppressMessages({
  library(DBI)
  library(duckdb)
  library(glue)
  library(optparse)
})

op <- OptionParser()
op <- add_option(op, "--input", default = "")
op <- add_option(op, "--tool", default = "")
op <- add_option(op, "--output-contract", dest = "output_contract", default = "")
op <- add_option(op, "--threads", type = "integer", default = 1L)
op <- add_option(op, "--run", type = "integer", default = 1L)
op <- add_option(op, "--timing-file", dest = "timing_file", default = "")
op <- add_option(op, "--skip-lines", dest = "skip_lines", type = "integer", default = 0L)
op <- add_option(op, "--output", default = "")
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
if (!nzchar(opt$input) || !file.exists(opt$input)) {
  die("--input must name an existing output file")
}
if (!nzchar(opt$tool) || !nzchar(opt$output_contract) || !nzchar(opt$output)) {
  die("--tool, --output-contract, and --output are required")
}
if (opt$threads < 1L || opt$run < 1L || opt$skip_lines < 0L) {
  die("--threads and --run must be positive; --skip-lines must be non-negative")
}

input <- normalizePath(opt$input)
drv <- duckdb(dbdir = ":memory:")
con <- dbConnect(drv)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
input_sql <- as.character(dbQuoteString(con, input))

fingerprint <- dbGetQuery(
  con,
  glue(
    "WITH rows AS (
       SELECT hash(COLUMNS(*))::UBIGINT AS h
       FROM read_csv(
         {input_sql}, delim = '\\t', header = true, skip = {opt$skip_lines},
         quote = '', escape = '', all_varchar = true
       )
     ), receipt AS (
       SELECT
         count(*)::UBIGINT AS row_count,
         bit_xor(h)::UBIGINT AS xor_hash,
         sum((h & 4294967295::UBIGINT)::HUGEINT)::HUGEINT AS low32_sum,
         sum((h >> 32)::HUGEINT)::HUGEINT AS high32_sum
       FROM rows
     )
     SELECT
       row_count::VARCHAR AS row_count,
       xor_hash::VARCHAR AS xor_hash,
       low32_sum::VARCHAR AS low32_sum,
       high32_sum::VARCHAR AS high32_sum
     FROM receipt"
  )
)

sha256 <- strsplit(system2("sha256sum", input, stdout = TRUE), " +")[[1L]][1L]
lines <- strsplit(system2("wc", c("-l", input), stdout = TRUE), " +")[[1L]]
lines <- lines[nzchar(lines)][[1L]]
timing_file <- if (nzchar(opt$timing_file)) basename(opt$timing_file) else ""

receipt <- data.frame(
  tool = opt$tool,
  output_contract = opt$output_contract,
  threads = opt$threads,
  run = opt$run,
  timing_file = timing_file,
  row_count = fingerprint$row_count,
  bytes = as.character(file.info(input)$size),
  lines = lines,
  sha256 = sha256,
  multiset_checked = TRUE,
  xor_hash = fingerprint$xor_hash,
  low32_sum = fingerprint$low32_sum,
  high32_sum = fingerprint$high32_sum,
  stringsAsFactors = FALSE
)

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(receipt, opt$output, row.names = FALSE, quote = TRUE)
