#!/usr/bin/env Rscript
# Exact statistical audit over a DuckVEP annotation dump.
#
# Reads a *_annotations.parquet dump emitted by corpus_differential.R and
# computes the VEP --gff statistical table: for each comparator
# engine and VEP consequence stratum, use the UNION of (variant, transcript) pairs
# that either VEP or the engine emitted, count exact SO-term-set discordance and
# emission misses/extras, and attach a Clopper-Pearson 95% upper bound via
# binom.test. This is a report, not a pass/fail gate. VEP remains the sole oracle.

suppressMessages({
  library(optparse)
  library(duckdb)
  library(DBI)
  library(glue)
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

op <- OptionParser(
  usage = "%prog [options]  (default: newest test/duckvep/conformance/results/*_annotations.parquet)"
)
op <- add_option(op, "--annotations", default = "")
op <- add_option(op, "--out", default = "")
op <- add_option(op, "--audit-out", dest = "audit_out", default = "")
op <- add_option(op, "--pairs-out", dest = "pairs_out", default = "")
opt <- parse_args(op)

data_dir <- file.path(root, "test", "duckvep", "conformance", "results")
if (!nzchar(opt$annotations)) {
  candidates <- list.files(
    data_dir,
    pattern = "_annotations[.]parquet$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(candidates) == 0L) {
    stop(
      "no *_annotations.parquet dump found; run corpus_differential.R first"
    )
  }
  info <- file.info(candidates)
  opt$annotations <- candidates[order(info$mtime, decreasing = TRUE)[1L]]
}
opt$annotations <- normalizePath(opt$annotations, mustWork = TRUE)
stem <- sub("_annotations[.]parquet$", "", basename(opt$annotations))
out_dir <- dirname(opt$annotations)
if (!nzchar(opt$out)) {
  opt$out <- file.path(out_dir, glue("{stem}_statistical_conformance.csv"))
}
if (!nzchar(opt$audit_out)) {
  opt$audit_out <- file.path(out_dir, glue("{stem}_methodology_audit.csv"))
}
if (!nzchar(opt$pairs_out)) {
  opt$pairs_out <- file.path(out_dir, glue("{stem}_pairs.parquet"))
}

upper95 <- function(k, n) {
  if (is.na(n) || n <= 0L) {
    return(NA_real_)
  }
  stats::binom.test(k, n, conf.level = 0.95)$conf.int[[2L]]
}

con <- dbConnect(duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
sql_q <- function(x) as.character(dbQuoteString(con, x))

invisible(dbExecute(
  con,
  glue(
    "CREATE TABLE ann AS SELECT * FROM read_parquet({sql_q(opt$annotations)})"
  )
))
ann_columns <- dbListFields(con, "ann")
if (!("status" %in% ann_columns)) {
  invisible(dbExecute(con, "ALTER TABLE ann ADD COLUMN status VARCHAR"))
}
if (!("reason" %in% ann_columns)) {
  invisible(dbExecute(con, "ALTER TABLE ann ADD COLUMN reason VARCHAR"))
}
invisible(dbExecute(
  con,
  "CREATE TABLE ann_pairs AS
   SELECT
     source,
     any_value(corpus) AS corpus,
     variant_id,
     coalesce(tx, '') AS tx,
     any_value(chrom) AS chrom,
     any_value(pos) AS pos,
     any_value(ref) AS ref,
     any_value(alt) AS alt,
     any_value(var_type) AS var_type,
     any_value(length_bin) AS length_bin,
     coalesce(string_agg(DISTINCT nullif(term, ''), '&' ORDER BY nullif(term, '')), '') AS consequence,
     any_value(impact) AS impact,
     CASE
       WHEN source = 'vep' THEN 'oracle'
       WHEN count(*) FILTER (WHERE status = 'unresolved') > 0 THEN 'unresolved'
       ELSE coalesce(any_value(status), 'unknown')
     END AS engine_status,
     string_agg(DISTINCT reason, ';' ORDER BY reason) AS engine_reason
   FROM ann
   LEFT JOIN LATERAL unnest(string_split(coalesce(consequence, ''), '&')) u(term) ON true
   GROUP BY source, variant_id, coalesce(tx, '')"
))
engines <- dbGetQuery(
  con,
  "SELECT DISTINCT source FROM ann_pairs WHERE source <> 'vep' ORDER BY source"
)$source
if (length(engines) == 0L) {
  stop("annotation dump has no comparator source rows")
}

pair_queries <- vapply(
  engines,
  function(engine) {
    glue(
      "WITH
       v AS (SELECT * FROM ann_pairs WHERE source = 'vep'),
       e AS (SELECT * FROM ann_pairs WHERE source = {sql_q(engine)}),
       u AS (
         SELECT variant_id, tx FROM v
         UNION
         SELECT variant_id, tx FROM e
       )
       SELECT
         {sql_q(engine)} AS engine,
         coalesce(v.corpus, e.corpus, '') AS corpus,
         u.variant_id,
         u.tx,
         coalesce(v.chrom, e.chrom, '') AS chrom,
         coalesce(v.pos, e.pos) AS pos,
         coalesce(v.ref, e.ref, '') AS ref,
         coalesce(v.alt, e.alt, '') AS alt,
         coalesce(v.var_type, e.var_type, '') AS var_type,
         coalesce(v.length_bin, e.length_bin, '') AS length_bin,
         v.consequence AS vep_consequence,
         e.consequence AS engine_consequence,
         v.impact AS vep_impact,
         e.impact AS engine_impact,
         coalesce(e.engine_status, 'missing') AS engine_status,
         e.engine_reason,
         CASE
           WHEN v.variant_id IS NULL THEN 'engine_extra'
           WHEN e.variant_id IS NULL THEN 'engine_missing'
           WHEN v.consequence = e.consequence THEN 'match'
           ELSE 'term_mismatch'
         END AS comparison
       FROM u
       LEFT JOIN v USING (variant_id, tx)
       LEFT JOIN e USING (variant_id, tx)"
    )
  },
  character(1)
)
invisible(dbExecute(
  con,
  glue("CREATE TABLE pairs AS {paste(pair_queries, collapse = ' UNION ALL ')}")
))
dir.create(dirname(opt$pairs_out), showWarnings = FALSE, recursive = TRUE)
invisible(dbExecute(
  con,
  glue("COPY pairs TO {sql_q(opt$pairs_out)} (FORMAT parquet, COMPRESSION zstd)")
))

stats <- dbGetQuery(
  con,
  "SELECT
     engine,
     corpus,
     CASE
       WHEN vep_consequence IS NULL THEN '(no_vep_emission)'
       WHEN vep_consequence = '' THEN '(empty_vep_terms)'
       ELSE vep_consequence
     END AS consequence_class,
     coalesce(nullif(vep_impact, ''), '(no_vep_emission)') AS impact,
     var_type,
     length_bin,
     engine_status,
     coalesce(engine_reason, '') AS engine_reason,
     count(*) AS n,
     count(*) FILTER (WHERE comparison = 'match') AS exact_agree,
     count(*) FILTER (WHERE comparison <> 'match') AS exact_discordant,
     count(*) FILTER (WHERE engine_status = 'unresolved') AS unresolved,
     count(*) FILTER (WHERE engine_status <> 'unresolved') AS resolved_n,
     count(*) FILTER (
       WHERE engine_status <> 'unresolved' AND comparison = 'match'
     ) AS resolved_agree,
     count(*) FILTER (
       WHERE engine_status <> 'unresolved' AND comparison <> 'match'
     ) AS resolved_discordant,
     count(*) FILTER (WHERE comparison = 'term_mismatch') AS term_mismatch,
     count(*) FILTER (WHERE comparison = 'engine_extra') AS engine_extra,
     count(*) FILTER (WHERE comparison = 'engine_missing') AS engine_missing
   FROM pairs
   GROUP BY ALL
   ORDER BY engine, resolved_discordant DESC, unresolved DESC, n DESC,
            consequence_class, var_type, length_bin"
)
stats$upper95 <- mapply(upper95, stats$resolved_discordant, stats$resolved_n)
stats <- stats[
  order(
    stats$engine,
    -stats$resolved_discordant,
    -stats$unresolved,
    -stats$n,
    stats$consequence_class,
    stats$var_type,
    stats$length_bin
  ),
]

audit <- do.call(
  rbind,
  lapply(split(stats, stats$engine), function(x) {
    n <- sum(x$n)
    resolved_n <- sum(x$resolved_n)
    resolved_discordant <- sum(x$resolved_discordant)
    data.frame(
      run_date = as.character(Sys.Date()),
      corpus = paste(sort(unique(x$corpus)), collapse = ","),
      engine = x$engine[[1L]],
      n = n,
      exact_agree = sum(x$exact_agree),
      exact_discordant = sum(x$exact_discordant),
      unresolved = sum(x$unresolved),
      resolved_n = resolved_n,
      resolved_agree = sum(x$resolved_agree),
      resolved_discordant = resolved_discordant,
      term_mismatch = sum(x$term_mismatch),
      engine_extra = sum(x$engine_extra),
      engine_missing = sum(x$engine_missing),
      upper95 = upper95(resolved_discordant, resolved_n),
      annotations = opt$annotations,
      stringsAsFactors = FALSE
    )
  })
)

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$audit_out), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$pairs_out), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(stats, opt$out, row.names = FALSE)
utils::write.csv(audit, opt$audit_out, row.names = FALSE)

cat("Statistical VEP --gff audit\n")
cat(glue("  annotations -> {opt$annotations}"), "\n", sep = "")
for (i in seq_len(nrow(audit))) {
  cat(
    glue(
      "  {audit$engine[i]}: exact {audit$exact_agree[i]}/{audit$n[i]}; resolved {audit$resolved_agree[i]}/{audit$resolved_n[i]}, discord {audit$resolved_discordant[i]} (<= {sprintf('%.2e', audit$upper95[i])} @95%); unresolved {audit$unresolved[i]}, extra {audit$engine_extra[i]}, missing {audit$engine_missing[i]}"
    ),
    "\n",
    sep = ""
  )
}
cat(glue("  strata -> {opt$out}"), "\n", sep = "")
cat(glue("  audit  -> {opt$audit_out}"), "\n", sep = "")
cat(glue("  pairs  -> {opt$pairs_out}"), "\n", sep = "")
