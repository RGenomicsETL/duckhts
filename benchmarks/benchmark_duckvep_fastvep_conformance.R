#!/usr/bin/env Rscript

# Compare DuckVEP and FastVEP with the same Ensembl VEP executable oracle.
# The oracle Parquet is produced by corpus_differential.R and retains every
# transcript pair. FastVEP CSQ rows are read through DuckHTS's header-driven
# VEP parser. Missing and extra transcript keys remain in the denominator.

suppressMessages({
  library(DBI)
  library(duckdb)
  library(glue)
  library(optparse)
})

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
op <- add_option(op, "--fastvep-vcf", dest = "fastvep_vcf", default = "")
op <- add_option(op, "--input-vcf", dest = "input_vcf", default = "")
op <- add_option(op, "--oracle-parquet", dest = "oracle_parquet", default = "")
op <- add_option(op, "--output", default = "")
op <- add_option(op, "--examples", default = "")
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
root <- normalizePath(root[[1L]], mustWork = TRUE)
required <- c(opt$extension, opt$input_vcf, opt$fastvep_vcf, opt$oracle_parquet)
missing <- required[!nzchar(required) | !file.exists(required)]
if (length(missing) != 0L) {
  die("missing input(s):\n{paste(missing, collapse = '\n')}")
}
if (!nzchar(opt$output)) {
  die("--output is required")
}

drv <- duckdb(
  dbdir = ":memory:",
  config = list(allow_unsigned_extensions = "true")
)
con <- dbConnect(drv)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
sql_q <- function(x) as.character(dbQuoteString(con, normalizePath(x)))

invisible(dbExecute(con, glue("LOAD {sql_q(opt$extension)}")))
invisible(dbExecute(con, "PRAGMA threads=1"))

# CSQ uses percent encoding for separators and '=' in VCF transport. Decode
# those bytes and remove the accession prefix before comparing HGVS suffixes
# with the suffixes retained by the executable-VEP oracle artifact.
decode_hgvs <- function(expression) {
  glue(
    "replace(replace(replace(replace(replace({expression},
       '%25', '%'), '%3D', '='), '%2C', ','), '%3B', ';'), '%20', ' ')"
  )
}

fast_hgvsc <- decode_hgvs(
  "regexp_replace(VEP_HGVSc[u.i], '^[^:]+:', '')"
)
fast_hgvsp <- decode_hgvs(
  "regexp_replace(VEP_HGVSp[u.i], '^[^:]+:', '')"
)

invisible(dbExecute(
  con,
  glue(
    "CREATE TEMP TABLE fastvep_raw AS
     SELECT
       ID AS variant_id, u.tx,
       VEP_Consequence[u.i] AS consequence,
       {fast_hgvsc} AS hgvsc,
       {fast_hgvsp} AS hgvsp
     FROM read_bcf({sql_q(opt$fastvep_vcf)})
     CROSS JOIN LATERAL UNNEST(VEP_Feature)
       WITH ORDINALITY AS u(tx, i)
     WHERE u.tx IS NOT NULL"
  )
))

invisible(dbExecute(
  con,
  "CREATE TEMP TABLE fastvep_pairs AS
   SELECT
     variant_id, tx,
     list_sort(list_distinct(flatten(list(string_split(consequence, '&')))))
       AS terms,
     CASE count(DISTINCT hgvsc) FILTER (WHERE hgvsc IS NOT NULL)
       WHEN 0 THEN NULL
       WHEN 1 THEN min(hgvsc) FILTER (WHERE hgvsc IS NOT NULL)
       ELSE concat(
         '__multiple__:',
         string_agg(DISTINCT hgvsc, '|' ORDER BY hgvsc)
           FILTER (WHERE hgvsc IS NOT NULL)
       )
     END AS hgvsc,
     CASE count(DISTINCT hgvsp) FILTER (WHERE hgvsp IS NOT NULL)
       WHEN 0 THEN NULL
       WHEN 1 THEN min(hgvsp) FILTER (WHERE hgvsp IS NOT NULL)
       ELSE concat(
         '__multiple__:',
         string_agg(DISTINCT hgvsp, '|' ORDER BY hgvsp)
           FILTER (WHERE hgvsp IS NOT NULL)
       )
     END AS hgvsp,
     count(*) AS source_rows
   FROM fastvep_raw
   GROUP BY variant_id, tx"
))

oracle_path <- sql_q(opt$oracle_parquet)
for (source in c("vep", "duckvep")) {
  table <- paste0(source, "_pairs")
  invisible(dbExecute(
    con,
    glue(
      "CREATE TEMP TABLE {table} AS
       SELECT
         variant_id, tx,
         list_sort(list_distinct(flatten(list(string_split(consequence, '&')))))
           AS terms,
         CASE count(DISTINCT hgvsc) FILTER (WHERE hgvsc IS NOT NULL)
           WHEN 0 THEN NULL
           WHEN 1 THEN min(hgvsc) FILTER (WHERE hgvsc IS NOT NULL)
           ELSE concat(
             '__multiple__:',
             string_agg(DISTINCT hgvsc, '|' ORDER BY hgvsc)
               FILTER (WHERE hgvsc IS NOT NULL)
           )
         END AS hgvsc,
         CASE count(DISTINCT hgvsp) FILTER (WHERE hgvsp IS NOT NULL)
           WHEN 0 THEN NULL
           WHEN 1 THEN min(hgvsp) FILTER (WHERE hgvsp IS NOT NULL)
           ELSE concat(
             '__multiple__:',
             string_agg(DISTINCT hgvsp, '|' ORDER BY hgvsp)
               FILTER (WHERE hgvsp IS NOT NULL)
           )
         END AS hgvsp,
         count(*) AS source_rows
       FROM read_parquet({oracle_path})
       WHERE source = {as.character(dbQuoteString(con, source))}
         AND tx IS NOT NULL
       GROUP BY variant_id, tx"
    )
  ))
}

metric_sql <- function(engine, table) {
  glue(
    "WITH joined AS (
       SELECT
         coalesce(v.variant_id, e.variant_id) AS variant_id,
         coalesce(v.tx, e.tx) AS tx,
         v.terms AS oracle_terms, e.terms AS engine_terms,
         v.hgvsc AS oracle_hgvsc, e.hgvsc AS engine_hgvsc,
         v.hgvsp AS oracle_hgvsp, e.hgvsp AS engine_hgvsp
       FROM vep_pairs v
       FULL OUTER JOIN {table} e USING (variant_id, tx)
     ), counts AS (
       SELECT
         count(*) AS union_keys,
         count(*) FILTER (
           WHERE oracle_terms IS NOT NULL AND engine_terms IS NOT NULL
         ) AS shared_keys,
         count(*) FILTER (
           WHERE oracle_terms IS NOT NULL AND engine_terms IS NULL
         ) AS missing_keys,
         count(*) FILTER (
           WHERE oracle_terms IS NULL AND engine_terms IS NOT NULL
         ) AS extra_keys,
         count(*) FILTER (WHERE oracle_terms = engine_terms) AS consequence_exact,
         count(*) FILTER (
           WHERE oracle_terms IS NOT NULL AND engine_terms IS NOT NULL
             AND oracle_terms != engine_terms
         ) AS consequence_discordant,
         count(*) FILTER (
           WHERE oracle_terms IS NOT NULL AND engine_terms IS NOT NULL
             AND oracle_hgvsc IS NOT DISTINCT FROM engine_hgvsc
         ) AS hgvsc_exact_including_absent,
         count(*) FILTER (
           WHERE oracle_terms IS NOT NULL AND engine_terms IS NOT NULL
             AND oracle_hgvsc IS NULL AND engine_hgvsc IS NULL
         ) AS hgvsc_both_absent,
         count(*) FILTER (
           WHERE oracle_terms IS NOT NULL AND engine_terms IS NOT NULL
             AND oracle_hgvsc = engine_hgvsc
         ) AS hgvsc_exact_present,
         count(*) FILTER (
           WHERE oracle_terms IS NOT NULL AND engine_terms IS NOT NULL
             AND oracle_hgvsc IS DISTINCT FROM engine_hgvsc
         ) AS hgvsc_discordant,
         count(*) FILTER (
           WHERE oracle_terms IS NOT NULL AND engine_terms IS NOT NULL
             AND oracle_hgvsp IS NOT DISTINCT FROM engine_hgvsp
         ) AS hgvsp_exact_including_absent,
         count(*) FILTER (
           WHERE oracle_terms IS NOT NULL AND engine_terms IS NOT NULL
             AND oracle_hgvsp IS NULL AND engine_hgvsp IS NULL
         ) AS hgvsp_both_absent,
         count(*) FILTER (
           WHERE oracle_terms IS NOT NULL AND engine_terms IS NOT NULL
             AND oracle_hgvsp = engine_hgvsp
         ) AS hgvsp_exact_present,
         count(*) FILTER (
           WHERE oracle_terms IS NOT NULL AND engine_terms IS NOT NULL
             AND oracle_hgvsp IS DISTINCT FROM engine_hgvsp
         ) AS hgvsp_discordant
       FROM joined
     )
     SELECT {as.character(dbQuoteString(con, engine))} AS engine,
       metric, value::UBIGINT AS n, denominator::UBIGINT AS denominator
     FROM counts
     CROSS JOIN LATERAL (VALUES
       ('union_keys', union_keys, union_keys),
       ('shared_keys', shared_keys, union_keys),
       ('missing_keys', missing_keys, union_keys),
       ('extra_keys', extra_keys, union_keys),
       ('consequence_exact', consequence_exact, union_keys),
       ('consequence_discordant', consequence_discordant, shared_keys),
       ('hgvsc_suffix_exact_including_absent', hgvsc_exact_including_absent, shared_keys),
       ('hgvsc_suffix_both_absent', hgvsc_both_absent, shared_keys),
       ('hgvsc_suffix_exact_present', hgvsc_exact_present, shared_keys),
       ('hgvsc_suffix_discordant', hgvsc_discordant, shared_keys),
       ('hgvsp_suffix_exact_including_absent', hgvsp_exact_including_absent, shared_keys),
       ('hgvsp_suffix_both_absent', hgvsp_both_absent, shared_keys),
       ('hgvsp_suffix_exact_present', hgvsp_exact_present, shared_keys),
       ('hgvsp_suffix_discordant', hgvsp_discordant, shared_keys)
     ) metrics(metric, value, denominator)"
  )
}

metrics <- dbGetQuery(
  con,
  paste(
    paste0("(", metric_sql("duckvep", "duckvep_pairs"), ")"),
    paste0("(", metric_sql("fastvep", "fastvep_pairs"), ")"),
    sep = " UNION ALL "
  )
)
metrics$run_date <- as.character(Sys.Date())
metrics$oracle <- "Ensembl VEP 116.0 executable indexed cache"
metrics$input_vcf_sha256 <- strsplit(
  system2("sha256sum", opt$input_vcf, stdout = TRUE),
  " +"
)[[1L]][1L]
metrics$fastvep_output_sha256 <- strsplit(
  system2("sha256sum", opt$fastvep_vcf, stdout = TRUE),
  " +"
)[[1L]][1L]
metrics$oracle_parquet_sha256 <- strsplit(
  system2("sha256sum", opt$oracle_parquet, stdout = TRUE),
  " +"
)[[1L]][1L]

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(metrics, opt$output, row.names = FALSE)

if (nzchar(opt$examples)) {
  examples <- dbGetQuery(
    con,
    "SELECT
       coalesce(v.variant_id, f.variant_id) AS variant_id,
       coalesce(v.tx, f.tx) AS transcript_id,
       array_to_string(v.terms, '&') AS vep_terms,
       array_to_string(f.terms, '&') AS fastvep_terms,
       v.hgvsc AS vep_hgvsc, f.hgvsc AS fastvep_hgvsc,
       v.hgvsp AS vep_hgvsp, f.hgvsp AS fastvep_hgvsp
     FROM vep_pairs v
     FULL OUTER JOIN fastvep_pairs f USING (variant_id, tx)
     WHERE v.terms IS DISTINCT FROM f.terms
        OR v.hgvsc IS DISTINCT FROM f.hgvsc
        OR v.hgvsp IS DISTINCT FROM f.hgvsp
     ORDER BY variant_id, transcript_id
     LIMIT 100"
  )
  dir.create(dirname(opt$examples), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(examples, opt$examples, row.names = FALSE)
}
