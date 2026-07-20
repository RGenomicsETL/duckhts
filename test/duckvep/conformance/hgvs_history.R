#!/usr/bin/env Rscript

# Record exact independent-event HGVS agreement from the pair-level Parquet
# emitted by corpus_differential.R. VEP remains the sole executable oracle.

suppressMessages({
  library(DBI)
  library(duckdb)
  library(glue)
  library(optparse)
})
options(rlang_backtrace_on_error = "none")

root <- tryCatch(
  system2(
    "git",
    c("rev-parse", "--show-toplevel"),
    stdout = TRUE,
    stderr = FALSE
  ),
  error = function(e) "."
)

op <- OptionParser()
op <- add_option(op, "--pairs", default = "")
op <- add_option(
  op,
  "--history",
  default = file.path(
    root,
    "test",
    "duckvep",
    "conformance",
    "data",
    "hgvs_history.csv"
  )
)
op <- add_option(
  op,
  "--source-revision",
  dest = "source_revision",
  default = ""
)
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
if (!nzchar(opt$pairs)) {
  die("--pairs is required")
}
opt$pairs <- normalizePath(opt$pairs, mustWork = TRUE)

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
if (!grepl("^[0-9a-f]{40}$", opt$source_revision)) {
  die("--source-revision must be a full 40-character Git revision")
}

con <- dbConnect(duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
sql_q <- function(x) as.character(dbQuoteString(con, x))
required <- c(
  "corpus", "model", "oracle_version", "oracle_build",
  "hgvsc_comparison", "hgvsp_comparison"
)
columns <- dbGetQuery(
  con,
  glue(
    "DESCRIBE SELECT * FROM read_parquet({sql_q(opt$pairs)})"
  )
)$column_name
missing <- setdiff(required, columns)
if (length(missing)) {
  die("HGVS pair artifact lacks column(s): {paste(missing, collapse = ', ')}")
}

rows <- dbGetQuery(
  con,
  glue(
    "WITH pairs AS (
       SELECT * FROM read_parquet({sql_q(opt$pairs)})
     ), comparisons AS (
       SELECT
         corpus, model, oracle_version, oracle_build,
         'hgvsc'::VARCHAR AS metric,
         hgvsc_comparison AS comparison
       FROM pairs
       UNION ALL
       SELECT
         corpus, model, oracle_version, oracle_build,
         'hgvsp', hgvsp_comparison
       FROM pairs
     )
     SELECT
       corpus,
       model,
       oracle_version,
       oracle_build,
       metric,
       comparison,
       count(*)::DOUBLE AS n
     FROM comparisons
     GROUP BY ALL
     ORDER BY corpus, model, metric, comparison"
  )
)
if (!nrow(rows)) {
  die("HGVS pair artifact contains no comparisons")
}
if (anyNA(rows) || any(!nzchar(rows$comparison))) {
  die("HGVS pair artifact contains incomplete comparison metadata")
}

history_columns <- c(
  "run_date", "source_revision", "artifact_md5", "corpus", "model",
  "oracle_version", "oracle_build", "metric", "comparison", "n"
)
rows <- data.frame(
  run_date = rep(as.character(Sys.Date()), nrow(rows)),
  source_revision = rep(opt$source_revision, nrow(rows)),
  artifact_md5 = rep(unname(tools::md5sum(opt$pairs)), nrow(rows)),
  rows,
  stringsAsFactors = FALSE,
  check.names = FALSE
)[, history_columns]

dir.create(dirname(opt$history), recursive = TRUE, showWarnings = FALSE)
if (file.exists(opt$history)) {
  old <- utils::read.csv(
    opt$history,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c(source_revision = "character", oracle_version = "character")
  )
  if (!identical(names(old), history_columns)) {
    die("history schema does not match: {opt$history}")
  }
  new_keys <- unique(paste(
    rows$source_revision,
    rows$corpus,
    rows$model,
    sep = "\034"
  ))
  old_keys <- paste(
    old$source_revision,
    old$corpus,
    old$model,
    sep = "\034"
  )
  rows <- rbind(old[!(old_keys %in% new_keys), , drop = FALSE], rows)
}
rows <- rows[
  order(
    rows$run_date,
    rows$source_revision,
    rows$corpus,
    rows$model,
    rows$metric,
    rows$comparison
  ),
  ,
  drop = FALSE
]
tmp <- tempfile("duckvep-hgvs-", dirname(opt$history), ".csv")
utils::write.csv(rows, tmp, row.names = FALSE)
if (!file.rename(tmp, opt$history)) {
  unlink(tmp)
  die("cannot replace history: {opt$history}")
}

summary <- aggregate(
  n ~ corpus + model + metric,
  rows[rows$source_revision == opt$source_revision, , drop = FALSE],
  sum
)
for (i in seq_len(nrow(summary))) {
  cat(
    glue(
      "{summary$corpus[[i]]} {summary$metric[[i]]}: ",
      "{format(summary$n[[i]], big.mark = ',', scientific = FALSE)} pairs"
    ),
    "\n",
    sep = ""
  )
}
cat(glue("history -> {opt$history}"), "\n", sep = "")
