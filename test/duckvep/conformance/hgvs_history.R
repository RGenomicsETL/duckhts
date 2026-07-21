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
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
root <- normalizePath(root[[1L]], mustWork = TRUE)
source(file.path(root, "scripts", "duckvep_evidence.R"), local = TRUE)
if (!nzchar(opt$pairs)) {
  die("--pairs is required")
}
opt$pairs <- normalizePath(opt$pairs, mustWork = TRUE)
revision <- duckvep_evidence_assert_checkout(
  root,
  allowed_outputs = opt$history,
  context = "DuckVEP HGVS history"
)

con <- dbConnect(duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
sql_q <- function(x) as.character(dbQuoteString(con, x))
required <- c(
  "source_revision", "extension_build_binding", "extension_sha256",
  "model_artifact_kind", "model_artifact_sha256",
  "reference_fasta_sha256", "reference_fai_sha256", "source_vcf_sha256",
  "input_vcf_sha256",
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

receipt <- dbGetQuery(
  con,
  glue(
    "SELECT
       count(DISTINCT source_revision) AS source_revisions,
       any_value(source_revision) AS source_revision,
       count(DISTINCT extension_build_binding) AS extension_bindings,
       any_value(extension_build_binding) AS extension_build_binding,
       count(DISTINCT extension_sha256) AS extension_digests,
       any_value(extension_sha256) AS extension_sha256,
       count(DISTINCT model_artifact_kind) AS model_artifact_kinds,
       any_value(model_artifact_kind) AS model_artifact_kind,
       count(DISTINCT model_artifact_sha256) AS model_artifact_digests,
       any_value(model_artifact_sha256) AS model_artifact_sha256,
       count(DISTINCT reference_fasta_sha256) AS reference_fasta_digests,
       any_value(reference_fasta_sha256) AS reference_fasta_sha256,
       count(DISTINCT reference_fai_sha256) AS reference_fai_digests,
       any_value(reference_fai_sha256) AS reference_fai_sha256,
       count(DISTINCT source_vcf_sha256) AS source_vcf_digests,
       any_value(source_vcf_sha256) AS source_vcf_sha256,
       count(DISTINCT input_vcf_sha256) AS input_vcf_digests,
       any_value(input_vcf_sha256) AS input_vcf_sha256
     FROM read_parquet({sql_q(opt$pairs)})"
  )
)
count_columns <- grep("(revisions|bindings|digests|kinds)$", names(receipt), value = TRUE)
if (nrow(receipt) != 1L || anyNA(receipt) ||
    any(unlist(receipt[count_columns], use.names = FALSE) != 1)) {
  die("HGVS pair artifact has incomplete or non-constant build receipts")
}
if (!identical(receipt$source_revision[[1L]], revision)) {
  die(
    "HGVS pair artifact revision {receipt$source_revision[[1L]]} is not ",
    "the current clean checkout {revision}"
  )
}
if (!identical(
  receipt$extension_build_binding[[1L]],
  "htslib_distclean_make_release"
)) {
  die("HGVS pair artifact was not produced by a clean vendored release build")
}
digest_columns <- c(
  "extension_sha256", "model_artifact_sha256", "reference_fasta_sha256",
  "reference_fai_sha256", "source_vcf_sha256", "input_vcf_sha256"
)
if (any(!vapply(
  receipt[digest_columns],
  function(value) grepl("^[0-9a-f]{64}$", value[[1L]]),
  logical(1L)
))) {
  die("HGVS pair artifact contains an invalid SHA-256 receipt")
}
if (!nzchar(receipt$model_artifact_kind[[1L]])) {
  die("HGVS pair artifact omits the model artifact kind")
}

rows <- dbGetQuery(
  con,
  glue(
    "WITH pairs AS (
       SELECT * FROM read_parquet({sql_q(opt$pairs)})
     ), comparisons AS (
       SELECT
         source_revision, extension_build_binding, extension_sha256,
         model_artifact_kind, model_artifact_sha256,
         reference_fasta_sha256, reference_fai_sha256, source_vcf_sha256,
         input_vcf_sha256,
         corpus, model, oracle_version, oracle_build,
         'hgvsc'::VARCHAR AS metric,
         hgvsc_comparison AS comparison
       FROM pairs
       UNION ALL
       SELECT
         source_revision, extension_build_binding, extension_sha256,
         model_artifact_kind, model_artifact_sha256,
         reference_fasta_sha256, reference_fai_sha256, source_vcf_sha256,
         input_vcf_sha256,
         corpus, model, oracle_version, oracle_build,
         'hgvsp', hgvsp_comparison
       FROM pairs
     )
     SELECT
       source_revision,
       extension_build_binding,
       extension_sha256,
       model_artifact_kind,
       model_artifact_sha256,
       reference_fasta_sha256,
       reference_fai_sha256,
       source_vcf_sha256,
       input_vcf_sha256,
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
  "run_date", "source_revision", "artifact_md5", "artifact_sha256",
  "extension_build_binding", "extension_sha256", "model_artifact_kind",
  "model_artifact_sha256", "reference_fasta_sha256",
  "reference_fai_sha256", "source_vcf_sha256", "input_vcf_sha256",
  "corpus", "model",
  "oracle_version", "oracle_build", "metric", "comparison", "n"
)
rows <- data.frame(
  run_date = rep(as.character(Sys.Date()), nrow(rows)),
  artifact_md5 = rep(unname(tools::md5sum(opt$pairs)), nrow(rows)),
  artifact_sha256 = rep(duckvep_evidence_sha256(opt$pairs), nrow(rows)),
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
    colClasses = c(
      source_revision = "character",
      artifact_md5 = "character",
      artifact_sha256 = "character",
      extension_sha256 = "character",
      model_artifact_sha256 = "character",
      reference_fasta_sha256 = "character",
      reference_fai_sha256 = "character",
      source_vcf_sha256 = "character",
      input_vcf_sha256 = "character",
      oracle_version = "character"
    )
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
utils::write.csv(rows, tmp, row.names = FALSE, na = "")
if (!file.rename(tmp, opt$history)) {
  unlink(tmp)
  die("cannot replace history: {opt$history}")
}

duckvep_evidence_assert_checkout(
  root,
  revision,
  opt$history,
  context = "DuckVEP HGVS history"
)

summary <- aggregate(
  n ~ corpus + model + metric,
  rows[rows$source_revision == revision, , drop = FALSE],
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
