#!/usr/bin/env Rscript

# End-to-end DuckVEP worker for the FastVEP comparison. The outer benchmark
# process pins this worker and records GNU time. This worker deliberately
# includes extension/model loading, sequential VCF decoding, explicit sorting,
# annotation, text projection, and COPY to a real local file.

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
op <- add_option(op, "--model", default = "")
op <- add_option(op, "--input", default = "")
op <- add_option(op, "--output", default = "")
op <- add_option(op, "--threads", type = "integer", default = 1L)
op <- add_option(op, "--distance", type = "integer", default = 5000L)
op <- add_option(op, "--memory-limit", dest = "memory_limit", default = "4GB")
op <- add_option(
  op,
  "--profile-json",
  dest = "profile_json",
  default = "",
  help = "optional DuckDB JSON profile path for the measured COPY query"
)
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
root <- normalizePath(root[[1L]], mustWork = TRUE)

required <- c(opt$extension, opt$model, opt$input)
missing <- required[!nzchar(required) | !file.exists(required)]
if (length(missing) != 0L) {
  die("missing input(s):\n{paste(missing, collapse = '\n')}")
}
if (!nzchar(opt$output)) {
  die("--output is required")
}
if (opt$threads < 1L || opt$threads > 1024L) {
  die("--threads must be from 1 through 1024")
}
if (opt$distance < 0L || opt$distance > 2^32 - 1) {
  die("--distance must fit an unsigned 32-bit integer")
}

extension <- normalizePath(opt$extension)
model <- normalizePath(opt$model)
input <- normalizePath(opt$input)
output <- normalizePath(opt$output, mustWork = FALSE)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

drv <- duckdb(
  dbdir = ":memory:",
  config = list(allow_unsigned_extensions = "true")
)
con <- dbConnect(drv)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
sql_q <- function(x) as.character(dbQuoteString(con, x))

invisible(dbExecute(con, glue("LOAD {sql_q(extension)}")))
invisible(dbExecute(con, glue("PRAGMA threads={opt$threads}")))
invisible(dbExecute(con, glue("SET memory_limit = {sql_q(opt$memory_limit)}")))
invisible(dbExecute(con, "SET preserve_insertion_order = false"))
if (nzchar(opt$profile_json)) {
  profile_json <- normalizePath(opt$profile_json, mustWork = FALSE)
  dir.create(dirname(profile_json), recursive = TRUE, showWarnings = FALSE)
  invisible(dbExecute(con, "PRAGMA enable_profiling = 'json'"))
  invisible(dbExecute(
    con,
    glue("PRAGMA profiling_output = {sql_q(profile_json)}")
  ))
}
invisible(dbExecute(
  con,
  glue("ATTACH {sql_q(model)} AS duckvep_bench_model (READ_ONLY)")
))
invisible(dbExecute(
  con,
  "CREATE TEMP TABLE duckvep_bench_regions AS
   SELECT seq_region, name
   FROM duckvep_bench_model.duckvep_sequence_regions"
))
invisible(dbExecute(
  con,
  "CREATE TEMP TABLE duckvep_bench_transcript_labels AS
   SELECT transcript_index, gene_stable_id, transcript_stable_id, strand
   FROM duckvep_bench_model.model_transcripts"
))

region_query <- paste(
  "SELECT seq_region, sequence_length",
  "FROM duckvep_bench_model.duckvep_sequence_regions",
  "ORDER BY seq_region"
)
transcript_query <- paste(
  "SELECT transcript_index, seq_region, transcript_start, transcript_end,",
  "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence,",
  "codon_table, pre_cds_sequence, post_cds_sequence",
  "FROM duckvep_bench_model.duckvep_transcripts",
  "ORDER BY seq_region, transcript_start, transcript_index"
)
exon_query <- paste(
  "SELECT transcript_index, exon_start, exon_end, exon_cdna_start,",
  "exon_cdna_end, phase, end_phase",
  "FROM duckvep_bench_model.duckvep_exons",
  "ORDER BY transcript_index, exon_cdna_start"
)
mature_mirna_query <- paste(
  "SELECT transcript_index, mature_mirna_start, mature_mirna_end",
  "FROM duckvep_bench_model.duckvep_mature_mirna",
  "ORDER BY transcript_index, mature_mirna_start"
)
peptide_edit_query <- paste(
  "SELECT transcript_index, protein_position, alternate_amino_acid",
  "FROM duckvep_bench_model.duckvep_peptide_edits",
  "ORDER BY transcript_index, protein_position"
)

loaded <- dbGetQuery(
  con,
  glue(
    "SELECT loaded FROM duckvep_model_load(
       'fastvep_comparison',
       {sql_q(region_query)}, {sql_q(transcript_query)}, {sql_q(exon_query)},
       mature_mirna_query := {sql_q(mature_mirna_query)},
       peptide_edit_query := {sql_q(peptide_edit_query)},
       transcript_coverage_complete := TRUE)"
  )
)$loaded
if (length(loaded) != 1L || !isTRUE(loaded[[1L]])) {
  die("DuckVEP model load failed")
}
invisible(dbExecute(con, "DETACH duckvep_bench_model"))

# FastVEP's native tab writer has a fixed 17-column contract. DuckVEP emits the
# same columns here so real-file costs are comparable. Fields not owned by the
# current rich consequence row remain '-' rather than being fabricated. This
# projection intentionally does not load regulatory/motif intervals because
# FastVEP's transcript cache does not contain them.
copy_query <- glue(
  "COPY (
     WITH source AS (
       SELECT
         row_number() OVER ()::UBIGINT AS record_index,
         CHROM AS chrom, POS::UBIGINT AS position, ID AS variant_id,
         REF AS reference, ALT AS alternates
       FROM read_bcf(
         {sql_q(input)}, scan_mode := 'sequential', decompression_threads := 0
       )
     ),
     alleles AS (
       SELECT
         s.record_index, a.alt_index, s.chrom, s.position, s.variant_id,
         s.reference, a.alternate
       FROM source s
       CROSS JOIN UNNEST(s.alternates) WITH ORDINALITY AS a(alternate, alt_index)
       WHERE regexp_full_match(s.reference, '[ACGTNacgtn]+')
         AND regexp_full_match(a.alternate, '[ACGTNacgtn]+')
         AND upper(s.reference) <> upper(a.alternate)
     ),
     prepared AS (
       SELECT
         a.record_index, a.alt_index, r.seq_region, a.chrom, a.position,
         a.variant_id, upper(a.reference) AS reference,
         upper(a.alternate) AS alternate
       FROM alleles a
       JOIN duckvep_bench_regions r
         ON r.name = regexp_replace(a.chrom, '^chr', '')
     ),
     annotated AS (
       SELECT
         v.*,
         unnest(_duckvep_annotate_small_rich(
           'fastvep_comparison', v.seq_region, v.position,
           v.reference, v.alternate, {opt$distance}
         )) AS annotation
       FROM (
         SELECT *
         FROM prepared
         ORDER BY seq_region, position, record_index, alt_index
       ) v
     )
     SELECT
       coalesce(
         a.variant_id,
         concat(a.chrom, ':', a.position, ':', a.reference, ':', a.alternate)
       ) AS Uploaded_variation,
       concat(a.chrom, ':', a.position) AS Location,
       a.alternate AS Allele,
       coalesce(t.gene_stable_id, '-') AS Gene,
       coalesce(t.transcript_stable_id, '-') AS Feature,
       CASE WHEN a.annotation.transcript_index IS NULL THEN '-' ELSE 'Transcript' END
         AS Feature_type,
       replace(a.annotation.consequence, '&', ',') AS Consequence,
       coalesce(a.annotation.cdna_position::VARCHAR, '') AS cDNA_position,
       coalesce(a.annotation.cds_position::VARCHAR, '') AS CDS_position,
       coalesce(a.annotation.protein_position::VARCHAR, '') AS Protein_position,
       CASE
         WHEN a.annotation.reference_amino_acid IS NULL
           AND a.annotation.alternate_amino_acid IS NULL THEN '-'
         ELSE concat(
           coalesce(a.annotation.reference_amino_acid, ''), '/',
           coalesce(a.annotation.alternate_amino_acid, '')
         )
       END AS Amino_acids,
       '-' AS Codons,
       '-' AS Existing_variation,
       a.annotation.impact AS IMPACT,
       '-' AS DISTANCE,
       coalesce(t.strand::VARCHAR, '-') AS STRAND,
       '-' AS FLAGS
     FROM annotated a
     LEFT JOIN duckvep_bench_transcript_labels t
       ON t.transcript_index = a.annotation.transcript_index
   ) TO {sql_q(output)} (
     FORMAT CSV, DELIMITER E'\\t', HEADER TRUE, QUOTE '', ESCAPE ''
   )"
)

invisible(dbExecute(con, copy_query))
