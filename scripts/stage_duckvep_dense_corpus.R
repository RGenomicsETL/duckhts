#!/usr/bin/env Rscript

# Stage a deterministic, annotation-dense small-variant benchmark relation.
# Network access is deliberately absent: the model and VCF must already exist.

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
op <- add_option(op, "--model-database", dest = "model_database")
op <- add_option(op, "--source-vcf", dest = "source_vcf")
op <- add_option(op, "--output-database", dest = "output_database")
op <- add_option(op, "--tile-size", dest = "tile_size", type = "integer", default = 250000L)
op <- add_option(
  op,
  "--tiles-per-axis",
  dest = "tiles_per_axis",
  type = "integer",
  default = 64L
)
op <- add_option(op, "--threads", type = "integer", default = 4L)
op <- add_option(op, "--overwrite", action = "store_true", default = FALSE)
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
required <- c("model_database", "source_vcf", "output_database")
missing_options <- required[vapply(
  required,
  function(name) is.null(opt[[name]]) || !nzchar(opt[[name]]),
  logical(1L)
)]
if (length(missing_options)) {
  die("missing required option(s): {paste0('--', gsub('_', '-', missing_options), collapse = ', ')}")
}
if (opt$tile_size < 1L) die("--tile-size must be positive")
if (opt$tiles_per_axis < 1L) die("--tiles-per-axis must be positive")
if (opt$threads < 1L) die("--threads must be positive")

input_paths <- c(opt$extension, opt$model_database, opt$source_vcf)
missing_paths <- input_paths[!file.exists(input_paths)]
if (length(missing_paths)) {
  die("missing input(s):\n{paste(missing_paths, collapse = '\n')}")
}

extension <- normalizePath(opt$extension)
model_database <- normalizePath(opt$model_database)
source_vcf <- normalizePath(opt$source_vcf)
dir.create(
  dirname(opt$output_database),
  recursive = TRUE,
  showWarnings = FALSE
)
output_database <- file.path(
  normalizePath(dirname(opt$output_database)),
  basename(opt$output_database)
)
input_paths <- c(extension, model_database, source_vcf)
if (output_database %in% input_paths) {
  die("output database must not replace the extension, model database, or source VCF")
}
file_identity <- function(path) {
  if (!file.exists(path) || !identical(Sys.info()[["sysname"]], "Linux")) {
    return(NA_character_)
  }
  value <- system2(
    "stat",
    c("-Lc", "%d:%i", "--", path),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(value, "status")
  if (!is.null(status) && status != 0L) die("stat failed for {path}")
  trimws(value[[1L]])
}
if (file.exists(output_database)) {
  output_identity <- file_identity(output_database)
  input_identities <- vapply(input_paths, file_identity, character(1L))
  if (!is.na(output_identity) && output_identity %in% input_identities) {
    die("output database is a hard-link alias of an input")
  }
  if (!isTRUE(opt$overwrite)) {
    die("output already exists; pass --overwrite to replace it: {output_database}")
  }
}
staged_database <- paste0(
  output_database,
  ".partial.",
  Sys.getpid()
)
unlink(c(staged_database, paste0(staged_database, ".wal")))
on.exit(unlink(c(staged_database, paste0(staged_database, ".wal"))), add = TRUE)

drv <- duckdb(config = list(allow_unsigned_extensions = "true"))
con <- dbConnect(drv, dbdir = staged_database)
on.exit({
  if (!is.null(con)) dbDisconnect(con, shutdown = TRUE)
}, add = TRUE)
sql_q <- function(x) as.character(dbQuoteString(con, x))

invisible(dbExecute(con, glue("SET threads = {opt$threads}")))
invisible(dbExecute(con, glue("LOAD {sql_q(extension)}")))
invisible(dbExecute(
  con,
  glue("ATTACH {sql_q(model_database)} AS duckvep_model (READ_ONLY)")
))

primary_contigs <- paste(
  c(as.character(1:22), "X"),
  collapse = "', '"
)
small_ncrna <- paste(
  c(
    "misc_RNA", "snRNA", "miRNA", "snoRNA", "rRNA", "scaRNA",
    "Mt_tRNA", "Mt_rRNA", "ribozyme", "sRNA", "vault_RNA"
  ),
  collapse = "', '"
)

invisible(dbExecute(
  con,
  glue(
    "CREATE TEMP TABLE dense_all_tiles AS
     SELECT
       r.seq_region,
       r.seq_region_name AS chrom,
       tile_id::UINTEGER AS tile_id,
       (tile_id * {opt$tile_size} + 1)::UBIGINT AS tile_start,
       least(
         (tile_id + 1) * {opt$tile_size},
         r.sequence_length
       )::UBIGINT AS tile_end
     FROM duckvep_model.main.model_regions r,
          range(ceil(r.sequence_length::DOUBLE / {opt$tile_size})::BIGINT) q(tile_id)
     WHERE r.seq_region_name IN ('{primary_contigs}')"
  )
))

invisible(dbExecute(
  con,
  glue(
    "CREATE TEMP TABLE dense_transcript_metrics AS
     SELECT
       d.seq_region,
       d.tile_id,
       count(*)::UBIGINT AS transcript_count,
       count(*) FILTER (WHERE t.cds_start != 0)::UBIGINT AS coding_transcript_count,
       count(*) FILTER (WHERE t.transcript_biotype = 'lncRNA')::UBIGINT
         AS lncrna_transcript_count,
       count(*) FILTER (
         WHERE t.transcript_biotype IN ('{small_ncrna}')
       )::UBIGINT AS small_ncrna_transcript_count,
       count(DISTINCT t.gene_index)::UBIGINT AS gene_count
     FROM dense_all_tiles d
     JOIN duckvep_model.main.model_transcripts t
       ON t.seq_region = d.seq_region
      AND t.transcript_start <= d.tile_end
      AND t.transcript_end >= d.tile_start
     GROUP BY d.seq_region, d.tile_id"
  )
))

invisible(dbExecute(
  con,
  "CREATE TEMP TABLE dense_exon_metrics AS
   SELECT
     d.seq_region,
     d.tile_id,
     count(*)::UBIGINT AS exon_membership_count,
     count(*) FILTER (
       WHERE t.cds_start != 0
         AND e.exon_end >= t.cds_start
         AND e.exon_start <= t.cds_end
     )::UBIGINT AS coding_exon_membership_count
   FROM dense_all_tiles d
   JOIN duckvep_model.main.model_transcripts t
     ON t.seq_region = d.seq_region
    AND t.transcript_start <= d.tile_end
    AND t.transcript_end >= d.tile_start
   CROSS JOIN UNNEST(t.exons) x(e)
   WHERE e.exon_start <= d.tile_end
     AND e.exon_end >= d.tile_start
   GROUP BY d.seq_region, d.tile_id"
))

invisible(dbExecute(
  con,
  "CREATE TEMP TABLE dense_regulation_metrics AS
   SELECT
     d.seq_region,
     d.tile_id,
     count(*) FILTER (
       WHERE f.feature_kind = 1
     )::UBIGINT AS regulatory_region_count,
     count(*) FILTER (
       WHERE f.feature_kind = 2
     )::UBIGINT AS motif_feature_count
   FROM dense_all_tiles d
   JOIN duckvep_model.main.duckvep_regulation_features f
     ON f.seq_region = d.seq_region
    AND f.feature_start <= d.tile_end
    AND f.feature_end >= d.tile_start
   GROUP BY d.seq_region, d.tile_id"
))

invisible(dbExecute(
  con,
  "CREATE TABLE dense_tile_metrics AS
   SELECT
     d.*,
     coalesce(t.transcript_count, 0)::UBIGINT AS transcript_count,
     coalesce(t.coding_transcript_count, 0)::UBIGINT AS coding_transcript_count,
     coalesce(t.lncrna_transcript_count, 0)::UBIGINT AS lncrna_transcript_count,
     coalesce(t.small_ncrna_transcript_count, 0)::UBIGINT
       AS small_ncrna_transcript_count,
     coalesce(t.gene_count, 0)::UBIGINT AS gene_count,
     coalesce(e.exon_membership_count, 0)::UBIGINT AS exon_membership_count,
     coalesce(e.coding_exon_membership_count, 0)::UBIGINT
       AS coding_exon_membership_count,
     coalesce(f.regulatory_region_count, 0)::UBIGINT AS regulatory_region_count,
     coalesce(f.motif_feature_count, 0)::UBIGINT AS motif_feature_count
   FROM dense_all_tiles d
   LEFT JOIN dense_transcript_metrics t USING (seq_region, tile_id)
   LEFT JOIN dense_exon_metrics e USING (seq_region, tile_id)
   LEFT JOIN dense_regulation_metrics f USING (seq_region, tile_id)
   ORDER BY d.seq_region, d.tile_id"
))

invisible(dbExecute(
  con,
  glue(
    "CREATE TABLE dense_tiles AS
     WITH axes AS (
       SELECT seq_region, tile_id, 'transcript_overlap' AS axis,
              transcript_count AS metric FROM dense_tile_metrics
       UNION ALL
       SELECT seq_region, tile_id, 'coding_transcript_overlap',
              coding_transcript_count FROM dense_tile_metrics
       UNION ALL
       SELECT seq_region, tile_id, 'lncrna_overlap',
              lncrna_transcript_count FROM dense_tile_metrics
       UNION ALL
       SELECT seq_region, tile_id, 'small_ncrna_overlap',
              small_ncrna_transcript_count FROM dense_tile_metrics
       UNION ALL
       SELECT seq_region, tile_id, 'exon_membership',
              exon_membership_count FROM dense_tile_metrics
       UNION ALL
       SELECT seq_region, tile_id, 'coding_exon_membership',
              coding_exon_membership_count FROM dense_tile_metrics
       UNION ALL
       SELECT seq_region, tile_id, 'regulatory_region',
              regulatory_region_count FROM dense_tile_metrics
       UNION ALL
       SELECT seq_region, tile_id, 'motif_feature',
              motif_feature_count FROM dense_tile_metrics
     ), ranked AS (
       SELECT *, row_number() OVER (
         PARTITION BY axis
         ORDER BY metric DESC, seq_region, tile_id
       ) AS rank_in_axis
       FROM axes
     ), selected AS (
       SELECT
         seq_region,
         tile_id,
         list(axis ORDER BY axis) AS selection_categories
       FROM ranked
       WHERE rank_in_axis <= {opt$tiles_per_axis}
       GROUP BY seq_region, tile_id
     )
     SELECT m.*, s.selection_categories
     FROM selected s
     JOIN dense_tile_metrics m USING (seq_region, tile_id)
     ORDER BY m.seq_region, m.tile_id"
  )
))

source_sql <- glue(
  "WITH records AS (
     SELECT
       (row_number() OVER () - 1)::UBIGINT AS source_record_index,
       CHROM, POS, REF, ALT
     FROM read_bcf({sql_q(source_vcf)}, scan_mode := 'sequential')
   )
   SELECT
     source_record_index,
     generate_subscripts(ALT, 1)::UINTEGER AS source_alt_index,
     CHROM, POS, REF, unnest(ALT) AS alternate
   FROM records"
)
source_alt_count <- dbGetQuery(
  con,
  glue("SELECT count(*)::UBIGINT AS n FROM ({source_sql})")
)$n[[1L]]

invisible(dbExecute(
  con,
  glue(
    "CREATE TABLE dense_source_alleles AS
     WITH source AS ({source_sql}), selected AS (
       SELECT
         source.source_record_index,
         source.source_alt_index,
         r.seq_region,
         source.POS::UBIGINT AS position,
         upper(source.REF) AS reference,
         upper(source.alternate) AS alternate
       FROM source
       JOIN duckvep_model.main.model_regions r
         ON regexp_replace(source.CHROM, '^chr', '') = r.seq_region_name
       JOIN dense_tiles d
         ON d.seq_region = r.seq_region
        AND source.POS BETWEEN d.tile_start AND d.tile_end
     )
     SELECT
       source_record_index,
       source_alt_index,
       seq_region,
       position,
       reference,
       alternate,
       regexp_full_match(reference, '[ACGTN]+')
         AND regexp_full_match(alternate, '[ACGTN]+') AS supported_literal,
       CASE
         WHEN NOT regexp_full_match(reference, '[ACGTN]+')
           THEN 'non_literal_reference'
         WHEN NOT regexp_full_match(alternate, '[ACGTN]+')
           THEN 'non_literal_alternate'
         ELSE NULL
       END AS exclusion_reason,
       CASE
         WHEN length(reference) = 1 AND length(alternate) = 1
           THEN 'single_base_substitution'
         WHEN length(reference) = length(alternate)
           THEN 'length_preserving_multibase'
         WHEN length(reference) < length(alternate)
           THEN 'lengthening'
         ELSE 'shortening'
       END AS raw_allele_shape,
       length(reference)::UINTEGER AS reference_length,
       length(alternate)::UINTEGER AS alternate_length
     FROM selected
     ORDER BY seq_region, position, reference, alternate,
              source_record_index, source_alt_index"
  )
))

invisible(dbExecute(
  con,
  "CREATE TABLE bench_dense_variants AS
   SELECT source_record_index, source_alt_index,
          seq_region, position, reference, alternate
   FROM dense_source_alleles
   WHERE supported_literal
   ORDER BY seq_region, position, reference, alternate,
            source_record_index, source_alt_index"
))
invisible(dbExecute(
  con,
  "CREATE TABLE bench_dense_snv_variants AS
   SELECT source_record_index, source_alt_index,
          seq_region, position, reference, alternate
   FROM dense_source_alleles
   WHERE supported_literal AND raw_allele_shape = 'single_base_substitution'
   ORDER BY seq_region, position, reference, alternate,
            source_record_index, source_alt_index"
))
invisible(dbExecute(
  con,
  "CREATE TABLE bench_dense_nonsnv_variants AS
   SELECT source_record_index, source_alt_index,
          seq_region, position, reference, alternate
   FROM dense_source_alleles
   WHERE supported_literal AND raw_allele_shape != 'single_base_substitution'
   ORDER BY seq_region, position, reference, alternate,
            source_record_index, source_alt_index"
))

model_receipt <- dbGetQuery(
  con,
  "SELECT * FROM duckvep_model.main.model_receipt"
)
if (nrow(model_receipt) != 1L) die("model_receipt must contain exactly one row")
region_ordinal_sha256 <- dbGetQuery(
  con,
  "SELECT sha256(string_agg(
     seq_region::VARCHAR || ':' || seq_region_name || ':' ||
     sequence_length::VARCHAR,
     '|' ORDER BY seq_region
   )) AS digest
   FROM duckvep_model.main.model_regions"
)$digest[[1L]]
sha256 <- function(path) {
  value <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(value, "status")
  if (!is.null(status) && status != 0L) die("sha256sum failed for {path}")
  strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
source_sha256 <- sha256(source_vcf)
staged_counts <- dbGetQuery(
  con,
  "SELECT
     (SELECT count(*) FROM dense_tiles)::UBIGINT AS selected_tile_count,
     (SELECT sum(tile_end - tile_start + 1) FROM dense_tiles)::UBIGINT
       AS selected_tile_bases,
     count(*)::UBIGINT AS selected_source_alt_count,
     count(*) FILTER (WHERE supported_literal)::UBIGINT AS supported_alt_count,
     count(*) FILTER (WHERE NOT supported_literal)::UBIGINT AS excluded_alt_count,
     max(reference_length)::UINTEGER AS max_reference_length,
     max(alternate_length)::UINTEGER AS max_alternate_length
   FROM dense_source_alleles"
)

receipt <- data.frame(
  schema_version = 2L,
  source_name = "ClinVar",
  source_path = source_vcf,
  source_sha256 = source_sha256,
  source_alt_count = source_alt_count,
  model_source_version = as.character(model_receipt$source_version[[1L]]),
  assembly = as.character(model_receipt$assembly[[1L]]),
  model_sha256 = as.character(model_receipt$model_sha256[[1L]]),
  region_ordinal_sha256 = as.character(region_ordinal_sha256),
  tile_size = opt$tile_size,
  tiles_per_axis = opt$tiles_per_axis,
  selected_tile_count = staged_counts$selected_tile_count,
  selected_tile_bases = staged_counts$selected_tile_bases,
  selected_source_alt_count = staged_counts$selected_source_alt_count,
  supported_alt_count = staged_counts$supported_alt_count,
  excluded_alt_count = staged_counts$excluded_alt_count,
  max_reference_length = staged_counts$max_reference_length,
  max_alternate_length = staged_counts$max_alternate_length,
  stringsAsFactors = FALSE
)
dbWriteTable(con, "dense_corpus_receipt", receipt, overwrite = TRUE)

invisible(dbExecute(con, "CHECKPOINT"))
dbDisconnect(con, shutdown = TRUE)
con <- NULL

if (!file.rename(staged_database, output_database)) {
  die("cannot publish staged database: {output_database}")
}
cat(
  glue(
    "staged {format(receipt$supported_alt_count, big.mark = ',', scientific = FALSE)} ",
    "literal ALT alleles from {format(receipt$selected_tile_count, big.mark = ',')} ",
    "annotation-dense tiles -> {output_database}"
  ),
  "\n",
  sep = ""
)
