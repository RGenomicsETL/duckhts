#!/usr/bin/env Rscript

# Record the stable DuckDB C-API path over a large, sorted variant stream.
# The checked-in fixture has one transcript: this measures vector/adapter and
# hot-transcript kernel cost, not whole-genome model throughput.

suppressMessages({
  library(DBI)
  library(duckdb)
  library(glue)
  library(optparse)
})
options(rlang_backtrace_on_error = "none")

root <- tryCatch(
  system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE),
  error = function(e) "."
)

op <- OptionParser()
op <- add_option(
  op,
  "--extension",
  default = file.path(root, "build", "release", "duckhts.duckdb_extension")
)
op <- add_option(
  op,
  "--model-sql",
  dest = "model_sql",
  default = file.path(
    root,
    "test",
    "duckvep",
    "conformance",
    "minimal_model.sql"
  )
)
op <- add_option(op, "--variants", type = "double", default = 10000000)
op <- add_option(op, "--passes", type = "integer", default = 3L)
op <- add_option(op, "--warmup", type = "double", default = 100000)
op <- add_option(op, "--threads", type = "integer", default = 1L)
op <- add_option(
  op,
  "--history",
  default = file.path(root, "benchmarks", "data", "duckvep_throughput.csv")
)
op <- add_option(
  op,
  "--source-revision",
  dest = "source_revision",
  default = ""
)
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
whole_count <- function(x, name) {
  if (!is.finite(x) || x < 1 || x != floor(x)) {
    die("{name} must be a positive integer")
  }
  format(x, scientific = FALSE, trim = TRUE)
}
variant_sql <- whole_count(opt$variants, "--variants")
warmup_sql <- whole_count(min(opt$warmup, opt$variants), "--warmup")
if (opt$passes < 1L) {
  die("--passes must be positive")
}
if (opt$threads < 1L) {
  die("--threads must be positive")
}
missing <- c(opt$extension, opt$model_sql)[
  !file.exists(c(opt$extension, opt$model_sql))
]
if (length(missing) != 0L) {
  die("missing input(s):\n{paste(missing, collapse = '\n')}")
}

drv <- duckdb(config = list(allow_unsigned_extensions = "true"))
con <- dbConnect(drv)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
sql_q <- function(x) as.character(dbQuoteString(con, x))
invisible(dbExecute(con, glue("SET threads = {opt$threads}")))
invisible(dbExecute(con, glue("LOAD {sql_q(normalizePath(opt$extension))}")))
invisible(dbExecute(
  con,
  paste(readLines(opt$model_sql, warn = FALSE), collapse = "\n")
))

model_name <- "throughput"
load_queries <- c(
  "SELECT seq_region FROM duckvep_sequence_regions ORDER BY seq_region",
  paste(
    "SELECT transcript_index, seq_region, transcript_start, transcript_end,",
    "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table",
    "FROM duckvep_transcripts ORDER BY seq_region, transcript_start, transcript_index"
  ),
  paste(
    "SELECT transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end,",
    "phase, end_phase FROM duckvep_exons",
    "ORDER BY transcript_index, exon_cdna_start"
  )
)
loaded <- dbGetQuery(
  con,
  glue(
    "SELECT loaded FROM duckvep_model_load(
       {sql_q(model_name)}, {sql_q(load_queries[1])},
       {sql_q(load_queries[2])}, {sql_q(load_queries[3])})"
  )
)$loaded
if (!identical(loaded, TRUE)) {
  die("model load failed")
}

annotation_query <- function(n) {
  glue(
    "WITH annotated AS (
       SELECT unnest(duckvep_annotate(
         {sql_q(model_name)}, 1::UINTEGER, 124::UBIGINT, 'T',
         CASE WHEN i % 2 = 0 THEN 'C' ELSE 'G' END,
         0::UBIGINT
       )) AS annotation
       FROM range({n}) r(i)
     )
     SELECT
       count(*)::DOUBLE AS annotated_rows,
       sum(length(annotation.consequence))::DOUBLE AS term_bytes
     FROM annotated"
  )
}

invisible(dbGetQuery(con, annotation_query(warmup_sql)))
elapsed <- numeric(opt$passes)
checks <- vector("list", opt$passes)
for (i in seq_len(opt$passes)) {
  timing <- system.time({
    checks[[i]] <- dbGetQuery(con, annotation_query(variant_sql))
  })
  elapsed[[i]] <- unname(timing[["elapsed"]])
}
check <- do.call(rbind, checks)
if (any(check$annotated_rows != opt$variants)) {
  die("annotation cardinality changed across benchmark passes")
}
if (length(unique(check$term_bytes)) != 1L) {
  die("annotation checksum changed across benchmark passes")
}

if (!nzchar(opt$source_revision)) {
  revision <- suppressWarnings(system2(
    "git",
    c("-C", root, "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  ))
  revision_status <- attr(revision, "status")
  if (!is.null(revision_status) && revision_status != 0L) {
    die("cannot determine source revision")
  }
  opt$source_revision <- trimws(revision[[1L]])
}

declared_version <- sub(
  "^version:[[:space:]]*",
  "",
  grep(
    "^version:[[:space:]]*",
    readLines(file.path(root, "description.yml"), warn = FALSE),
    value = TRUE
  )[[1L]]
)
loaded_version <- dbGetQuery(
  con,
  "SELECT extension_version FROM duckdb_extensions()
   WHERE extension_name = 'duckhts'"
)$extension_version[[1L]]
if (!identical(loaded_version, declared_version)) {
  die(
    "loaded extension version {loaded_version} does not match description.yml ",
    "({declared_version}); run `make configure release` first"
  )
}
duckdb_version <- dbGetQuery(con, "SELECT version() AS version")$version[[1L]]
cpu <- if (file.exists("/proc/cpuinfo")) {
  line <- grep(
    "^model name[[:space:]]*:",
    readLines("/proc/cpuinfo"),
    value = TRUE
  )[[1L]]
  trimws(sub("^[^:]+:", "", line))
} else {
  Sys.info()[["machine"]]
}

median_seconds <- stats::median(elapsed)
row <- data.frame(
  run_date = as.character(Sys.Date()),
  source_revision = opt$source_revision,
  declared_extension_version = declared_version,
  loaded_extension_version = loaded_version,
  duckdb_version = duckdb_version,
  host = unname(Sys.info()[["nodename"]]),
  cpu = cpu,
  workload = "fixture_one_transcript_sorted",
  input_order = "nondecreasing_seq_region_position",
  threads = opt$threads,
  variants = as.integer(opt$variants),
  transcripts = 1,
  exons = 2,
  passes = opt$passes,
  warmup_variants = as.integer(min(opt$warmup, opt$variants)),
  min_seconds = min(elapsed),
  median_seconds = median_seconds,
  max_seconds = max(elapsed),
  variants_per_second = opt$variants / median_seconds,
  ns_per_variant = median_seconds * 1e9 / opt$variants,
  annotated_rows = as.integer(check$annotated_rows[[1L]]),
  term_bytes = as.integer(check$term_bytes[[1L]]),
  stringsAsFactors = FALSE
)

if (nzchar(opt$history)) {
  dir.create(dirname(opt$history), recursive = TRUE, showWarnings = FALSE)
  rows <- row
  if (file.exists(opt$history)) {
    old <- utils::read.csv(
      opt$history,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (!identical(names(old), names(row))) {
      die("history schema does not match: {opt$history}")
    }
    same <- old$source_revision == row$source_revision &
      old$workload == row$workload &
      old$threads == row$threads &
      old$variants == row$variants
    rows <- rbind(old[!same, , drop = FALSE], row)
  }
  rows <- rows[
    order(rows$run_date, rows$source_revision, rows$workload),
    ,
    drop = FALSE
  ]
  tmp <- tempfile("duckvep-throughput-", dirname(opt$history), ".csv")
  utils::write.csv(rows, tmp, row.names = FALSE)
  if (!file.rename(tmp, opt$history)) {
    unlink(tmp)
    die("cannot replace history: {opt$history}")
  }
}

invisible(dbGetQuery(
  con,
  glue("SELECT duckvep_model_drop({sql_q(model_name)})")
))
cat(
  glue(
    "{format(opt$variants, big.mark = ',', scientific = FALSE)} sorted variants; median ",
    "{sprintf('%.3f', median_seconds)} s; ",
    "{sprintf('%.1f', row$ns_per_variant)} ns/variant"
  ),
  "\n",
  sep = ""
)
if (nzchar(opt$history)) {
  cat(glue("history -> {opt$history}"), "\n", sep = "")
}
