#!/usr/bin/env Rscript
# Exact statistical audit over a DuckVEP annotation dump.
#
# Reads a *_annotations.parquet dump emitted by corpus_differential.R and
# computes the VEP executable-oracle statistical table: for each comparator
# engine and VEP consequence stratum, use the UNION of (variant, transcript) pairs
# that either VEP or the engine emitted, count exact SO-term-set discordance and
# emission misses/extras, and attach a Clopper-Pearson 95% upper bound via
# binom.test. It writes every pair and summary before failing closed on any
# consequence discordance, missing/extra row, or DuckVEP unresolved state. VEP
# remains the sole oracle.

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
op <- add_option(op, "--nmd-out", dest = "nmd_out", default = "")
op <- add_option(
  op,
  "--history",
  default = "",
  help = "append or replace this revision/corpus/model in an audit history CSV"
)
op <- add_option(
  op,
  "--source-revision",
  dest = "source_revision",
  default = "",
  help = "source revision recorded with --history [current Git HEAD]"
)
op <- add_option(
  op,
  "--duckdb-memory-limit",
  dest = "duckdb_memory_limit",
  default = "8GB",
  help = "DuckDB memory cap for spillable pair construction [%default]"
)
op <- add_option(
  op,
  "--duckdb-threads",
  dest = "duckdb_threads",
  type = "integer",
  default = 4L,
  help = "DuckDB worker threads for evidence projection [%default]"
)
op <- add_option(
  op,
  "--pair-level-input",
  dest = "pair_level_input",
  action = "store_true",
  default = FALSE,
  help = paste(
    "assert that each source has at most one row per variant/transcript and",
    "skip legacy consequence-row normalization [%default]"
  )
)
opt <- parse_args(op)
if (opt$duckdb_threads < 1L) {
  stop("--duckdb-threads must be positive")
}

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
if (!nzchar(opt$nmd_out)) {
  opt$nmd_out <- file.path(out_dir, glue("{stem}_nmd_conformance.csv"))
}

upper95 <- function(k, n) {
  if (is.na(n) || n <= 0L) {
    return(NA_real_)
  }
  stats::binom.test(k, n, conf.level = 0.95)$conf.int[[2L]]
}

con <- dbConnect(duckdb())
spill_dir <- NULL
on.exit({
  try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  if (!is.null(spill_dir)) {
    unlink(spill_dir, recursive = TRUE, force = TRUE)
  }
}, add = TRUE)
sql_q <- function(x) as.character(dbQuoteString(con, x))
spill_dir <- tempfile("duckvep-statistical-spill-")
if (!dir.create(spill_dir)) {
  stop("cannot create DuckDB spill directory: ", spill_dir)
}
invisible(dbExecute(con, glue("SET temp_directory = {sql_q(spill_dir)}")))
invisible(dbExecute(
  con,
  glue("SET memory_limit = {sql_q(opt$duckdb_memory_limit)}")
))
invisible(dbExecute(con, glue("SET threads = {opt$duckdb_threads}")))
invisible(dbExecute(con, "SET preserve_insertion_order = false"))

so_spec_path <- file.path(
  root,
  "test",
  "duckvep",
  "conformance",
  "data",
  "so_consequences.tsv"
)
if (!file.exists(so_spec_path)) {
  stop("missing generated VEP SO table: ", so_spec_path)
}
so_spec <- utils::read.delim(
  so_spec_path,
  stringsAsFactors = FALSE,
  na.strings = "\\N"
)
invisible(dbWriteTable(con, "so_spec", so_spec, overwrite = TRUE))

ann_columns <- dbGetQuery(
  con,
  glue("DESCRIBE SELECT * FROM read_parquet({sql_q(opt$annotations)})")
)$column_name
compat_columns <- c(
  status = "NULL::VARCHAR",
  reason = "NULL::VARCHAR",
  model = "NULL::VARCHAR",
  oracle_version = "NULL::VARCHAR",
  oracle_build = "NULL::VARCHAR",
  nmd_prediction = "'not_measured'::VARCHAR"
)
missing_compat <- setdiff(names(compat_columns), ann_columns)
compat_projection <- if (length(missing_compat)) {
  paste0(
    ", ",
    paste(
      glue("{compat_columns[missing_compat]} AS {missing_compat}"),
      collapse = ", "
    )
  )
} else {
  ""
}
invisible(dbExecute(
  con,
  glue(
    "CREATE TEMP VIEW ann AS
     SELECT *{compat_projection}
     FROM read_parquet({sql_q(opt$annotations)})"
  )
))
if (isTRUE(opt$pair_level_input)) {
  duplicate_pair <- dbGetQuery(
    con,
    "SELECT source, variant_id, coalesce(tx, '') AS tx, count(*) AS n
     FROM ann
     GROUP BY source, variant_id, coalesce(tx, '')
     HAVING count(*) > 1
     LIMIT 1"
  )
  if (nrow(duplicate_pair) != 0L) {
    stop(glue(
      "pair-level input is not unique for source={duplicate_pair$source[[1L]]}, ",
      "variant_id={duplicate_pair$variant_id[[1L]]}, ",
      "tx={duplicate_pair$tx[[1L]]}: {duplicate_pair$n[[1L]]} rows"
    ))
  }
  invisible(dbExecute(
    con,
    "CREATE TEMP VIEW ann_pairs AS
     SELECT
       source,
       corpus,
       coalesce(model, '') AS model,
       coalesce(oracle_version, '') AS oracle_version,
       coalesce(oracle_build, '') AS oracle_build,
       variant_id,
       coalesce(tx, '') AS tx,
       chrom,
       pos,
       ref,
       alt,
       var_type,
       length_bin,
       coalesce(consequence, '') AS consequence,
       impact,
       CASE
         WHEN source = 'vep' THEN 'oracle'
         ELSE coalesce(status, 'unknown')
       END AS engine_status,
       reason AS engine_reason,
       coalesce(nmd_prediction, 'not_measured') AS nmd_prediction
     FROM ann"
  ))
} else {
  ann_pairs_path <- file.path(spill_dir, "ann_pairs.parquet")
  invisible(dbExecute(
    con,
    glue("COPY (
   SELECT
     source,
     any_value(corpus) AS corpus,
     any_value(coalesce(model, '')) AS model,
     any_value(coalesce(oracle_version, '')) AS oracle_version,
     any_value(coalesce(oracle_build, '')) AS oracle_build,
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
         WHEN count(*) FILTER (
           WHERE status IS NULL OR status NOT IN ('supported', 'unresolved')
         ) > 0 THEN 'invalid'
         WHEN count(*) FILTER (WHERE status = 'unresolved') > 0 THEN 'unresolved'
         ELSE 'supported'
       END AS engine_status,
     string_agg(DISTINCT reason, ';' ORDER BY reason) AS engine_reason,
     any_value(coalesce(nmd_prediction, 'not_measured')) AS nmd_prediction
   FROM ann
   LEFT JOIN LATERAL unnest(string_split(coalesce(consequence, ''), '&')) u(term) ON true
   GROUP BY source, variant_id, coalesce(tx, '')
   ) TO {sql_q(ann_pairs_path)} (FORMAT parquet, COMPRESSION zstd)")
  ))
  invisible(dbExecute(
    con,
    glue(
      "CREATE TEMP VIEW ann_pairs AS
       SELECT * FROM read_parquet({sql_q(ann_pairs_path)})"
    )
  ))
}
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
       e AS (SELECT * FROM ann_pairs WHERE source = {sql_q(engine)})
       SELECT
         {sql_q(engine)} AS engine,
         coalesce(v.corpus, e.corpus, '') AS corpus,
         coalesce(v.model, e.model, '') AS model,
         coalesce(v.oracle_version, e.oracle_version, '') AS oracle_version,
         coalesce(v.oracle_build, e.oracle_build, '') AS oracle_build,
         coalesce(v.variant_id, e.variant_id) AS variant_id,
         coalesce(v.tx, e.tx) AS tx,
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
         coalesce(v.nmd_prediction, 'not_measured') AS vep_nmd_prediction,
         coalesce(e.nmd_prediction, 'not_measured') AS engine_nmd_prediction,
         CASE
           -- A present row with no NMD observation remains unmeasured even if
           -- the other engine omitted its consequence row. Only measured
           -- emissions can establish that a missing peer is not comparable.
           WHEN (
             v.variant_id IS NOT NULL AND
             coalesce(v.nmd_prediction, 'not_measured') = 'not_measured'
           ) OR (
             e.variant_id IS NOT NULL AND
             coalesce(e.nmd_prediction, 'not_measured') = 'not_measured'
           )
             THEN 'not_measured'
           -- A missing peer does not turn an ineligible consequence into an
           -- NMD observation. Reserve not_comparable for an eligible NMD
           -- prediction whose peer consequence row is absent.
           WHEN (
             v.variant_id IS NULL AND e.nmd_prediction = 'not_applicable'
           ) OR (
             e.variant_id IS NULL AND v.nmd_prediction = 'not_applicable'
           )
             THEN 'not_measured'
           WHEN v.variant_id IS NULL OR e.variant_id IS NULL THEN 'not_comparable'
           WHEN coalesce(v.nmd_prediction, 'not_measured') =
                coalesce(e.nmd_prediction, 'not_measured') THEN 'match'
           ELSE 'mismatch'
         END AS nmd_comparison,
         CASE
           WHEN v.variant_id IS NULL THEN 'engine_extra'
           WHEN e.variant_id IS NULL THEN 'engine_missing'
           WHEN v.consequence = e.consequence THEN 'match'
           ELSE 'term_mismatch'
         END AS comparison
       FROM v
       FULL OUTER JOIN e USING (variant_id, tx)"
    )
  },
  character(1)
)
invisible(dbExecute(
  con,
  glue(
    "CREATE TEMP VIEW pairs_projection AS
     {paste(pair_queries, collapse = ' UNION ALL ')}"
  )
))
dir.create(dirname(opt$pairs_out), showWarnings = FALSE, recursive = TRUE)
invisible(dbExecute(
  con,
  glue(
    "COPY pairs_projection TO {sql_q(opt$pairs_out)}
     (FORMAT parquet, COMPRESSION zstd)"
  )
))
invisible(dbExecute(con, "DROP VIEW pairs_projection"))
invisible(dbExecute(
  con,
  glue(
    "CREATE TEMP VIEW pairs AS
     SELECT * FROM read_parquet({sql_q(opt$pairs_out)})"
  )
))

stats <- dbGetQuery(
  con,
  "WITH pair_terms AS (
     SELECT
       p.* EXCLUDE (comparison),
       CASE
         WHEN list_contains(
                string_split(coalesce(p.vep_consequence, ''), '&'),
                u.so_term
              )
          AND list_contains(
                string_split(coalesce(p.engine_consequence, ''), '&'),
                u.so_term
              ) THEN 'term_match'
         WHEN list_contains(
                string_split(coalesce(p.vep_consequence, ''), '&'),
                u.so_term
              ) THEN 'term_missing'
         ELSE 'term_extra'
       END AS comparison,
       u.so_term
     FROM pairs p
     CROSS JOIN UNNEST(list_distinct(list_concat(
       string_split(coalesce(p.vep_consequence, ''), '&'),
       string_split(coalesce(p.engine_consequence, ''), '&')
     ))) u(so_term)
     WHERE u.so_term <> ''
   ), strata AS (
     SELECT
       p.*,
       'consequence_set'::VARCHAR AS stratum_kind,
       CASE
         WHEN vep_consequence IS NULL THEN '(no_vep_emission)'
         WHEN vep_consequence = '' THEN '(empty_vep_terms)'
         ELSE vep_consequence
       END AS consequence_class,
       coalesce(nullif(vep_impact, ''), '(no_vep_emission)') AS stratum_impact
     FROM pairs p
     UNION ALL
     SELECT
       p.* EXCLUDE (so_term),
       'so_term'::VARCHAR AS stratum_kind,
       p.so_term AS consequence_class,
       coalesce(
         nullif(s.impact, ''),
         nullif(p.vep_impact, ''),
         nullif(p.engine_impact, ''),
         '(unknown)'
       ) AS stratum_impact
     FROM pair_terms p
     LEFT JOIN so_spec s ON s.SO_term = p.so_term
   )
   SELECT
     engine,
     corpus,
     model,
     oracle_version,
     oracle_build,
     stratum_kind,
     consequence_class,
     stratum_impact AS impact,
     var_type,
     length_bin,
     engine_status,
     coalesce(engine_reason, '') AS engine_reason,
     count(*) AS n,
     count(*) FILTER (WHERE comparison IN ('match', 'term_match')) AS exact_agree,
     count(*) FILTER (WHERE comparison NOT IN ('match', 'term_match')) AS exact_discordant,
       count(*) FILTER (WHERE engine_status = 'unresolved') AS unresolved,
       count(*) FILTER (
         WHERE engine_status NOT IN ('supported', 'unresolved')
       ) AS invalid_status,
       count(*) FILTER (WHERE engine_status = 'supported') AS resolved_n,
       count(*) FILTER (
         WHERE engine_status = 'supported'
           AND comparison IN ('match', 'term_match')
       ) AS resolved_agree,
       count(*) FILTER (
         WHERE engine_status = 'supported'
           AND comparison NOT IN ('match', 'term_match')
     ) AS resolved_discordant,
     count(*) FILTER (
       WHERE comparison IN ('term_mismatch', 'term_extra', 'term_missing')
     ) AS term_mismatch,
     count(*) FILTER (WHERE comparison IN ('engine_extra', 'term_extra')) AS engine_extra,
     count(*) FILTER (WHERE comparison IN ('engine_missing', 'term_missing')) AS engine_missing
   FROM strata
   GROUP BY ALL
   ORDER BY engine, stratum_kind, resolved_discordant DESC, unresolved DESC,
            n DESC, consequence_class, var_type, length_bin"
)
stats$upper95 <- mapply(upper95, stats$resolved_discordant, stats$resolved_n)
stats <- stats[
  order(
    stats$engine,
    stats$stratum_kind,
    -stats$resolved_discordant,
    -stats$unresolved,
    -stats$n,
    stats$consequence_class,
    stats$var_type,
    stats$length_bin
  ),
]

audit <- dbGetQuery(
  con,
  "SELECT
     engine,
     corpus,
     model,
     any_value(oracle_version) AS oracle_version,
     any_value(oracle_build) AS oracle_build,
     count(*) AS n,
     count(*) FILTER (WHERE comparison = 'match') AS exact_agree,
     count(*) FILTER (WHERE comparison <> 'match') AS exact_discordant,
     count(*) FILTER (WHERE engine_status = 'unresolved') AS unresolved,
     count(*) FILTER (
       WHERE engine_status NOT IN ('supported', 'unresolved')
     ) AS invalid_status,
     count(*) FILTER (WHERE engine_status = 'supported') AS resolved_n,
     count(*) FILTER (
       WHERE engine_status = 'supported' AND comparison = 'match'
     ) AS resolved_agree,
     count(*) FILTER (
       WHERE engine_status = 'supported' AND comparison <> 'match'
     ) AS resolved_discordant,
     count(*) FILTER (WHERE comparison = 'term_mismatch') AS term_mismatch,
     count(*) FILTER (WHERE comparison = 'engine_extra') AS engine_extra,
     count(*) FILTER (WHERE comparison = 'engine_missing') AS engine_missing
   FROM pairs
   GROUP BY engine, corpus, model
   ORDER BY engine, corpus, model"
)
audit$run_date <- as.character(Sys.Date())
audit$upper95 <- mapply(upper95, audit$resolved_discordant, audit$resolved_n)
audit$annotations <- opt$annotations

nmd_stats <- dbGetQuery(
  con,
  "SELECT
     engine,
     corpus,
     model,
     oracle_version,
     oracle_build,
     vep_nmd_prediction,
     engine_nmd_prediction,
     nmd_comparison,
     CASE
       WHEN vep_consequence IS NULL THEN '(no_vep_emission)'
       WHEN vep_consequence = '' THEN '(empty_vep_terms)'
       ELSE vep_consequence
     END AS consequence_class,
     var_type,
     length_bin,
     count(*) AS n
   FROM pairs
   WHERE nmd_comparison <> 'not_measured'
     AND (
       vep_nmd_prediction <> 'not_applicable' OR
       engine_nmd_prediction <> 'not_applicable'
     )
   GROUP BY ALL
   ORDER BY engine, nmd_comparison DESC, n DESC, consequence_class,
            var_type, length_bin"
)
nmd_audit <- dbGetQuery(
  con,
  "SELECT
     engine,
     corpus,
     model,
     count(*) AS n,
     count(*) FILTER (WHERE nmd_comparison = 'match') AS exact_agree,
     count(*) FILTER (WHERE nmd_comparison = 'mismatch') AS exact_discordant,
     count(*) FILTER (WHERE nmd_comparison = 'not_comparable') AS not_comparable,
     count(*) FILTER (
       WHERE vep_nmd_prediction = 'unresolved'
     ) AS oracle_unresolved,
     count(*) FILTER (
       WHERE engine_nmd_prediction = 'unresolved'
     ) AS engine_unresolved
   FROM pairs
   WHERE nmd_comparison <> 'not_measured'
     AND (
       vep_nmd_prediction <> 'not_applicable' OR
       engine_nmd_prediction <> 'not_applicable'
     )
   GROUP BY engine, corpus, model
   ORDER BY engine, corpus, model"
)
if (nrow(nmd_audit) != 0L) {
  comparable <- nmd_audit$exact_agree + nmd_audit$exact_discordant
  nmd_audit$upper95 <- mapply(
    upper95,
    nmd_audit$exact_discordant,
    comparable
  )
}
nmd_history_stats <- dbGetQuery(
  con,
  "SELECT
     engine,
     corpus,
     model,
     any_value(oracle_version) AS oracle_version,
     any_value(oracle_build) AS oracle_build,
     vep_nmd_prediction,
     engine_nmd_prediction,
     count(*) AS n,
     count(*) FILTER (WHERE nmd_comparison = 'match') AS exact_agree,
     count(*) FILTER (WHERE nmd_comparison <> 'match') AS exact_discordant,
     count(*) FILTER (
       WHERE engine_nmd_prediction = 'unresolved'
     ) AS unresolved,
     count(*) FILTER (
       WHERE nmd_comparison <> 'not_comparable'
         AND vep_nmd_prediction <> 'unresolved'
         AND engine_nmd_prediction <> 'unresolved'
     ) AS resolved_n,
     count(*) FILTER (
       WHERE nmd_comparison = 'match'
         AND vep_nmd_prediction <> 'unresolved'
         AND engine_nmd_prediction <> 'unresolved'
     ) AS resolved_agree,
     count(*) FILTER (
       WHERE nmd_comparison = 'mismatch'
         AND vep_nmd_prediction <> 'unresolved'
         AND engine_nmd_prediction <> 'unresolved'
     ) AS resolved_discordant,
     count(*) FILTER (WHERE nmd_comparison = 'mismatch') AS term_mismatch,
     count(*) FILTER (WHERE comparison = 'engine_extra') AS engine_extra,
     count(*) FILTER (WHERE comparison = 'engine_missing') AS engine_missing
   FROM pairs
   WHERE nmd_comparison <> 'not_measured'
     AND (
       vep_nmd_prediction <> 'not_applicable' OR
       engine_nmd_prediction <> 'not_applicable'
     )
   GROUP BY engine, corpus, model, vep_nmd_prediction,
            engine_nmd_prediction
   ORDER BY engine, corpus, model, vep_nmd_prediction,
            engine_nmd_prediction"
)
if (nrow(nmd_history_stats) != 0L) {
  nmd_history_stats$upper95 <- mapply(
    upper95,
    nmd_history_stats$resolved_discordant,
    nmd_history_stats$resolved_n
  )
}

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$audit_out), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$pairs_out), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$nmd_out), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(stats, opt$out, row.names = FALSE)
utils::write.csv(audit, opt$audit_out, row.names = FALSE)
utils::write.csv(nmd_stats, opt$nmd_out, row.names = FALSE)

core_failure <- audit[
  audit$exact_discordant != 0L |
    audit$unresolved != 0L |
    audit$invalid_status != 0L |
    audit$engine_extra != 0L |
    audit$engine_missing != 0L,
  ,
  drop = FALSE
]
nmd_failure_count <- dbGetQuery(
  con,
  "SELECT count(*) AS n
   FROM pairs
   WHERE nmd_comparison IN ('mismatch', 'not_comparable')
      OR (
        engine_nmd_prediction = 'unresolved' AND
        vep_nmd_prediction <> 'unresolved'
      )"
)$n[[1L]]
if (nrow(core_failure) != 0L || nmd_failure_count != 0L) {
  core_detail <- if (nrow(core_failure) == 0L) {
    "none"
  } else {
    paste(
      glue(
        "{core_failure$engine}: discord={core_failure$exact_discordant}, ",
        "unresolved={core_failure$unresolved}, ",
        "invalid_status={core_failure$invalid_status}, ",
        "extra={core_failure$engine_extra}, ",
        "missing={core_failure$engine_missing}"
      ),
      collapse = "; "
    )
  }
  stop(glue(
    "VEP executable audit failed closed ({core_detail}; ",
    "NMD failures={nmd_failure_count}). Every pair remains in {opt$pairs_out}; ",
    "summaries remain in {opt$out}, {opt$audit_out}, and {opt$nmd_out}."
  ))
}

if (nzchar(opt$history)) {
  if (!nzchar(opt$source_revision)) {
    revision <- suppressWarnings(system2(
      "git",
      c("-C", root, "rev-parse", "HEAD"),
      stdout = TRUE,
      stderr = FALSE
    ))
    revision_status <- attr(revision, "status")
    if (!is.null(revision_status) && revision_status != 0L) {
      stop("cannot determine source revision")
    }
    opt$source_revision <- trimws(revision[[1L]])
  }
  if (!nzchar(opt$source_revision)) {
    stop("empty source revision")
  }

  run_date <- as.character(Sys.Date())
  artifact_md5 <- unname(tools::md5sum(opt$annotations))
  history_columns <- c(
    "run_date",
    "source_revision",
    "artifact_md5",
    "engine",
    "corpus",
    "model",
    "oracle_version",
    "oracle_build",
    "stratum_kind",
    "consequence_class",
    "impact",
    "var_type",
    "length_bin",
    "engine_status",
    "engine_reason",
    "n",
    "exact_agree",
    "exact_discordant",
    "unresolved",
    "resolved_n",
    "resolved_agree",
    "resolved_discordant",
    "term_mismatch",
    "engine_extra",
    "engine_missing",
    "upper95"
  )
  history_stats <- data.frame(
    run_date = rep(run_date, nrow(stats)),
    source_revision = rep(opt$source_revision, nrow(stats)),
    artifact_md5 = rep(artifact_md5, nrow(stats)),
    stats,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )[, history_columns]
  history_all <- data.frame(
    run_date = rep(run_date, nrow(audit)),
    source_revision = rep(opt$source_revision, nrow(audit)),
    artifact_md5 = rep(artifact_md5, nrow(audit)),
    engine = audit$engine,
    corpus = audit$corpus,
    model = audit$model,
    oracle_version = audit$oracle_version,
    oracle_build = audit$oracle_build,
    stratum_kind = "all",
    consequence_class = "(all)",
    impact = "(all)",
    var_type = "(all)",
    length_bin = "(all)",
    engine_status = "(all)",
    engine_reason = "",
    n = audit$n,
    exact_agree = audit$exact_agree,
    exact_discordant = audit$exact_discordant,
    unresolved = audit$unresolved,
    resolved_n = audit$resolved_n,
    resolved_agree = audit$resolved_agree,
    resolved_discordant = audit$resolved_discordant,
    term_mismatch = audit$term_mismatch,
    engine_extra = audit$engine_extra,
    engine_missing = audit$engine_missing,
    upper95 = audit$upper95,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )[, history_columns]
  history_nmd <- data.frame()
  if (nrow(nmd_history_stats) != 0L) {
    history_nmd <- data.frame(
      run_date = rep(run_date, nrow(nmd_history_stats)),
      source_revision = rep(opt$source_revision, nrow(nmd_history_stats)),
      artifact_md5 = rep(artifact_md5, nrow(nmd_history_stats)),
      engine = nmd_history_stats$engine,
      corpus = nmd_history_stats$corpus,
      model = nmd_history_stats$model,
      oracle_version = nmd_history_stats$oracle_version,
      oracle_build = nmd_history_stats$oracle_build,
      stratum_kind = "nmd_prediction",
      consequence_class = nmd_history_stats$vep_nmd_prediction,
      impact = "(nmd)",
      var_type = "(all)",
      length_bin = "(all)",
      engine_status = nmd_history_stats$engine_nmd_prediction,
      engine_reason = "",
      n = nmd_history_stats$n,
      exact_agree = nmd_history_stats$exact_agree,
      exact_discordant = nmd_history_stats$exact_discordant,
      unresolved = nmd_history_stats$unresolved,
      resolved_n = nmd_history_stats$resolved_n,
      resolved_agree = nmd_history_stats$resolved_agree,
      resolved_discordant = nmd_history_stats$resolved_discordant,
      term_mismatch = nmd_history_stats$term_mismatch,
      engine_extra = nmd_history_stats$engine_extra,
      engine_missing = nmd_history_stats$engine_missing,
      upper95 = nmd_history_stats$upper95,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )[, history_columns]
  }
  history_rows <- rbind(history_all, history_stats, history_nmd)

  dir.create(dirname(opt$history), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(opt$history)) {
    old <- utils::read.csv(
      opt$history,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      colClasses = c(oracle_version = "character")
    )
    if (!identical(names(old), history_columns)) {
      stop("history schema does not match: ", opt$history)
    }
    new_keys <- unique(paste(
      history_rows$source_revision,
      history_rows$engine,
      history_rows$corpus,
      history_rows$model,
      sep = "\034"
    ))
    old_keys <- paste(
      old$source_revision,
      old$engine,
      old$corpus,
      old$model,
      sep = "\034"
    )
    old <- old[!(old_keys %in% new_keys), , drop = FALSE]
    history_rows <- rbind(old, history_rows)
  }
  history_rows <- history_rows[
    order(
      history_rows$run_date,
      history_rows$source_revision,
      history_rows$corpus,
      history_rows$model,
      history_rows$engine,
      history_rows$stratum_kind,
      history_rows$consequence_class,
      history_rows$impact,
      history_rows$var_type,
      history_rows$length_bin,
      history_rows$engine_status,
      history_rows$engine_reason
    ),
    ,
    drop = FALSE
  ]
  history_tmp <- tempfile(
    pattern = "duckvep-conformance-",
    tmpdir = dirname(opt$history),
    fileext = ".csv"
  )
  utils::write.csv(history_rows, history_tmp, row.names = FALSE)
  if (!file.rename(history_tmp, opt$history)) {
    unlink(history_tmp)
    stop("cannot replace history: ", opt$history)
  }
}

cat("Statistical VEP executable audit\n")
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
if (nrow(nmd_audit) != 0L) {
  for (i in seq_len(nrow(nmd_audit))) {
    comparable <- nmd_audit$exact_agree[i] + nmd_audit$exact_discordant[i]
    cat(
      glue(
        "  NMD {nmd_audit$engine[i]}: exact ",
        "{nmd_audit$exact_agree[i]}/{comparable}; ",
        "discord {nmd_audit$exact_discordant[i]} ",
        "(<= {sprintf('%.2e', nmd_audit$upper95[i])} @95%); ",
        "not comparable {nmd_audit$not_comparable[i]}, ",
        "oracle unresolved {nmd_audit$oracle_unresolved[i]}, ",
        "engine unresolved {nmd_audit$engine_unresolved[i]}"
      ),
      "\n",
      sep = ""
    )
  }
  cat(glue("  NMD    -> {opt$nmd_out}"), "\n", sep = "")
}
if (nzchar(opt$history)) {
  cat(glue("  history -> {opt$history}"), "\n", sep = "")
}
