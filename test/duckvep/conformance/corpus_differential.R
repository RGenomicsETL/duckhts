#!/usr/bin/env Rscript

# Compare DuckVEP with Ensembl VEP 116 over a deterministic sample of a VCF.
# VEP is the oracle. Every emitted transcript pair is retained, including
# DuckVEP rows marked unresolved; there is no hand-written supported-case filter.

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
results_dir <- file.path(root, "test", "duckvep", "conformance", "results")
ncores <- parallel::detectCores()
if (is.na(ncores)) {
  ncores <- 1L
}

op <- OptionParser(
  usage = "%prog [options]  (defaults run the in-tree witness differential)"
)
op <- add_option(op, "--corpus", default = "witnesses")
op <- add_option(op, "--vcf", default = "")
op <- add_option(
  op,
  "--event-mode",
  dest = "event_mode",
  default = "small",
  help = "event family: small or structural [default: %default]"
)
op <- add_option(
  op,
  "--gff",
  default = file.path(root, "test", "data", "duckvep", "minimal.gff3")
)
op <- add_option(
  op,
  "--cache-dir",
  dest = "cache_dir",
  default = "",
  help = paste(
    "VEP cache root; when set, use --cache --offline instead of --gff",
    "[default: use --gff]"
  )
)
op <- add_option(op, "--assembly", default = "GRCh38")
op <- add_option(op, "--species", default = "homo_sapiens")
op <- add_option(
  op,
  "--cache-version",
  dest = "cache_version",
  default = "",
  help = paste(
    "VEP cache release when it differs from the VEP executable release,",
    "for example 63 for Ensembl Genomes paired with VEP 116 [%default]"
  )
)
op <- add_option(
  op,
  "--fasta",
  default = file.path(root, "test", "data", "duckvep", "minimal.fa")
)
op <- add_option(op, "--database", default = ":memory:")
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
  ),
  help = "optional SQL that creates the four duckvep_* model relations [%default]"
)
op <- add_option(
  op,
  "--model-name",
  dest = "model_name",
  default = "differential"
)
op <- add_option(
  op,
  "--sample-per-shape",
  dest = "sample_per_shape",
  type = "integer",
  default = 10000L,
  help = "deterministic limit per allele type and length bin; 0 keeps all [%default]"
)
op <- add_option(op, "--seed", type = "integer", default = 1L)
op <- add_option(op, "--chrom", default = "")
op <- add_option(op, "--distance", type = "integer", default = 5000L)
op <- add_option(
  op,
  "--annotations-out",
  dest = "annotations_out",
  default = ""
)
op <- add_option(op, "--sample-vcf", dest = "sample_vcf", default = "")
op <- add_option(
  op,
  "--keep-sample-vcf",
  dest = "keep_sample_vcf",
  action = "store_true",
  default = FALSE
)
op <- add_option(
  op,
  "--extension",
  default = Sys.getenv(
    "DUCKHTS_EXT",
    file.path(root, "build", "release", "duckhts.duckdb_extension")
  )
)
op <- add_option(op, "--fork", default = as.character(max(1L, min(8L, ncores))))
op <- add_option(
  op,
  "--vep-prefix",
  dest = "vep_prefix",
  default = Sys.getenv("VEP_PREFIX", "/root/miniconda3/envs/vep"),
  help = "micromamba environment prefix containing VEP 116 [%default]"
)
op <- add_option(
  op,
  "--micromamba",
  default = Sys.getenv("MICROMAMBA", "micromamba"),
  help = "micromamba executable [%default]"
)
op <- add_option(
  op,
  "--nmd-plugin-dir",
  dest = "nmd_plugin_dir",
  default = Sys.getenv("VEP_PLUGIN_DIR", ""),
  help = paste(
    "directory containing the pinned release/116 NMD.pm; when set, compare",
    "DuckVEP variant-induced NMD predictions with the executable plugin"
  )
)
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
if (!(opt$event_mode %in% c("small", "structural"))) {
  die("--event-mode must be small or structural")
}
if (
  !requireNamespace("blit", quietly = TRUE) ||
    utils::packageVersion("blit") < "0.2.0.9000"
) {
  die(
    "the current WangLabCSU/blit checkout is required ",
    "(blit >= 0.2.0.9000)"
  )
}
oracle_mode <- if (nzchar(opt$cache_dir)) "cache" else "gff"
required_files <- c(opt$fasta, opt$extension)
external_model_database <- !nzchar(opt$model_sql) &&
  !identical(opt$database, ":memory:")
if (external_model_database) {
  required_files <- c(required_files, opt$database)
}
if (identical(oracle_mode, "gff")) {
  required_files <- c(opt$gff, required_files)
} else if (!dir.exists(opt$cache_dir)) {
  die("VEP cache root does not exist: {opt$cache_dir}")
}
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) != 0L) {
  die("missing input(s):\n{paste(missing_files, collapse = '\n')}")
}
if (opt$sample_per_shape < 0L) {
  die("--sample-per-shape must be non-negative")
}
if (opt$distance < 0L) {
  die("--distance must be non-negative")
}
nmd_plugin_sha256 <- ""
nmd_state_plugin_sha256 <- ""
nmd_state_plugin <- file.path(
  root,
  "test",
  "duckvep",
  "conformance",
  "DuckVEPNMDState.pm"
)
if (nzchar(opt$nmd_plugin_dir)) {
  if (identical(opt$event_mode, "structural")) {
    die("the structural differential does not yet compare NMD plugin output")
  }
  if (!dir.exists(opt$nmd_plugin_dir)) {
    die("VEP plugin directory does not exist: {opt$nmd_plugin_dir}")
  }
  opt$nmd_plugin_dir <- normalizePath(opt$nmd_plugin_dir)
  nmd_plugin <- file.path(opt$nmd_plugin_dir, "NMD.pm")
  if (!file.exists(nmd_plugin)) {
    die("NMD.pm does not exist in {opt$nmd_plugin_dir}")
  }
  sha_line <- system2("sha256sum", nmd_plugin, stdout = TRUE, stderr = FALSE)
  if (length(sha_line) != 1L) {
    die("cannot checksum {nmd_plugin}")
  }
  nmd_plugin_sha256 <- strsplit(trimws(sha_line), "[[:space:]]+")[[1L]][1L]
  expected_nmd_sha256 <-
    "1e38bd67783ff09bad2775d09235dd77f23a7e5ade50fa56d4777235092e0eeb"
  if (!identical(nmd_plugin_sha256, expected_nmd_sha256)) {
    die(
      "NMD.pm checksum mismatch: expected {expected_nmd_sha256}, ",
      "found {nmd_plugin_sha256}"
    )
  }
  if (!file.exists(nmd_state_plugin)) {
    die("missing NMD coordinate observer: {nmd_state_plugin}")
  }
  state_sha_line <- system2(
    "sha256sum",
    nmd_state_plugin,
    stdout = TRUE,
    stderr = FALSE
  )
  if (length(state_sha_line) != 1L) {
    die("cannot checksum {nmd_state_plugin}")
  }
  nmd_state_plugin_sha256 <-
    strsplit(trimws(state_sha_line), "[[:space:]]+")[[1L]][1L]
}
nmd_oracle_enabled <- nzchar(nmd_plugin_sha256)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
if (!nzchar(opt$annotations_out)) {
  opt$annotations_out <- file.path(
    results_dir,
    glue("{opt$corpus}_annotations.parquet")
  )
}

temporary_files <- character()
cleanup <- function() unlink(temporary_files, recursive = TRUE, force = TRUE)
on.exit(cleanup(), add = TRUE)
nmd_plugin_run_dir <- opt$nmd_plugin_dir
if (nmd_oracle_enabled) {
  nmd_plugin_run_dir <- tempfile(pattern = "duckvep-nmd-plugins-")
  if (!dir.create(nmd_plugin_run_dir)) {
    die("cannot create staged NMD plugin directory: {nmd_plugin_run_dir}")
  }
  temporary_files <- c(temporary_files, nmd_plugin_run_dir)
  copied <- file.copy(
    c(nmd_plugin, nmd_state_plugin),
    nmd_plugin_run_dir,
    overwrite = FALSE
  )
  if (!all(copied)) {
    die("cannot stage NMD oracle plugins in {nmd_plugin_run_dir}")
  }
}

drv <- duckdb(
  dbdir = if (external_model_database) ":memory:" else opt$database,
  config = list(allow_unsigned_extensions = "true")
)
con <- dbConnect(drv)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
sql_q <- function(x) as.character(dbQuoteString(con, x))
invisible(dbExecute(con, glue("LOAD {sql_q(normalizePath(opt$extension))}")))
invisible(tryCatch(
  dbExecute(con, "LOAD json"),
  error = function(e) {
    die("DuckDB JSON support is required: {conditionMessage(e)}")
  }
))

if (external_model_database) {
  invisible(dbExecute(
    con,
    glue(
      "ATTACH {sql_q(normalizePath(opt$database))} ",
      "AS duckvep_model_source (READ_ONLY)"
    )
  ))
  model_relations <- dbGetQuery(
    con,
    "SELECT table_name
     FROM information_schema.tables
     WHERE table_catalog = 'duckvep_model_source'
       AND table_schema = 'main'
       AND table_name LIKE 'duckvep_%'
     ORDER BY table_name"
  )$table_name
  for (relation in model_relations) {
    relation_id <- as.character(dbQuoteIdentifier(con, relation))
    invisible(dbExecute(
      con,
      glue(
        "CREATE TEMP VIEW {relation_id} AS SELECT * FROM ",
        "duckvep_model_source.main.{relation_id}"
      )
    ))
  }
}

if (nzchar(opt$model_sql)) {
  if (!file.exists(opt$model_sql)) {
    die("model SQL does not exist: {opt$model_sql}")
  }
  invisible(dbExecute(
    con,
    paste(readLines(opt$model_sql, warn = FALSE), collapse = "\n")
  ))
}

needed_relations <- c(
  "duckvep_sequence_regions",
  "duckvep_transcripts",
  "duckvep_exons",
  "duckvep_transcript_names"
)
present_relations <- dbGetQuery(
  con,
  paste(
    "SELECT DISTINCT table_name FROM information_schema.tables",
    "WHERE table_schema = 'main'"
  )
)$table_name
missing_relations <- setdiff(needed_relations, present_relations)
if (length(missing_relations) != 0L) {
  die(
    "model database is missing relation(s): {paste(missing_relations, collapse = ', ')}"
  )
}

region_columns <- dbGetQuery(
  con,
  "SELECT column_name FROM duckdb_columns()
   WHERE table_name = 'duckvep_sequence_regions'"
)$column_name
complete_coverage <- "sequence_length" %in% region_columns
model_query_relation <- function(relation) {
  relation_id <- as.character(dbQuoteIdentifier(con, relation))
  if (external_model_database) {
    return(glue("duckvep_model_source.main.{relation_id}"))
  }
  relation_id
}
load_queries <- c(
  paste(
    "SELECT seq_region",
    if (complete_coverage) ", sequence_length" else "",
    "FROM", model_query_relation("duckvep_sequence_regions"),
    "ORDER BY seq_region"
  ),
  paste(
    "SELECT transcript_index, seq_region, transcript_start, transcript_end,",
    "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table,",
    "pre_cds_sequence, post_cds_sequence",
    "FROM", model_query_relation("duckvep_transcripts"),
    "ORDER BY seq_region, transcript_start, transcript_index"
  ),
  paste(
    "SELECT transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end,",
    "phase, end_phase FROM", model_query_relation("duckvep_exons"),
    "ORDER BY transcript_index, exon_cdna_start"
  )
)
load_options <- character()
if ("duckvep_mature_mirna" %in% present_relations) {
  mature_mirna_query <- paste(
    "SELECT transcript_index, mature_mirna_start, mature_mirna_end",
    "FROM", model_query_relation("duckvep_mature_mirna"),
    "ORDER BY transcript_index, mature_mirna_start"
  )
  load_options <- c(
    load_options,
    glue("mature_mirna_query := {sql_q(mature_mirna_query)}")
  )
}
if ("duckvep_peptide_edits" %in% present_relations) {
  peptide_edit_query <- paste(
    "SELECT transcript_index, protein_position, alternate_amino_acid",
    "FROM", model_query_relation("duckvep_peptide_edits"),
    "ORDER BY transcript_index, protein_position"
  )
  load_options <- c(
    load_options,
    glue("peptide_edit_query := {sql_q(peptide_edit_query)}")
  )
}
if (complete_coverage) {
  load_options <- c(load_options, "transcript_coverage_complete := TRUE")
}
loaded <- dbGetQuery(
  con,
  glue(
    "SELECT loaded FROM duckvep_model_load(
       {sql_q(opt$model_name)}, {sql_q(load_queries[1])},
       {sql_q(load_queries[2])}, {sql_q(load_queries[3])}
       {if (length(load_options)) paste0(', ', paste(load_options, collapse = ', ')) else ''})"
  )
)$loaded
if (length(loaded) != 1L || !isTRUE(loaded[[1L]])) {
  die("DuckVEP model load failed")
}

source_vcf <- opt$vcf
generate_structural <- identical(opt$event_mode, "structural") &&
  !nzchar(source_vcf)
if (!nzchar(source_vcf) && !generate_structural) {
  source_vcf <- tempfile(fileext = ".vcf")
  temporary_files <- c(temporary_files, source_vcf)
  rc <- system2(
    "Rscript",
    c(
      file.path(root, "test", "duckvep", "conformance", "generate_witnesses.R"),
      "--gff",
      opt$gff,
      "--fasta",
      opt$fasta,
      "--tx",
      "DUCK1-201",
      "--out",
      source_vcf,
      "--ext",
      opt$extension
    )
  )
  if (rc != 0L || !file.exists(source_vcf)) die("witness generation failed")
} else if (nzchar(source_vcf) && !file.exists(source_vcf)) {
  die("VCF does not exist: {source_vcf}")
}

chrom_filter <- ""
chroms <- trimws(strsplit(opt$chrom, ",", fixed = TRUE)[[1L]])
chroms <- chroms[nzchar(chroms)]
if (length(chroms) != 0L) {
  chrom_filter <- glue(
    "AND CHROM IN ({paste(vapply(chroms, sql_q, character(1)), collapse = ', ')})"
  )
}
generated_chrom_filter <- ""
if (length(chroms) != 0L) {
  generated_chrom_filter <- glue(
    "AND r.name IN ({paste(vapply(chroms, sql_q, character(1)), collapse = ', ')})"
  )
}
sample_filter <- if (opt$sample_per_shape == 0L) {
  ""
} else {
  glue("WHERE sample_rank <= {opt$sample_per_shape}")
}
generated_filter <- if (opt$sample_per_shape == 0L) {
  ""
} else {
  glue("WHERE geometry_rank <= {opt$sample_per_shape}")
}

if (generate_structural) {
  invisible(dbExecute(
    con,
    glue(
      "CREATE OR REPLACE TEMP TABLE duckvep_generated_structural_source AS
       WITH transcripts AS MATERIALIZED (
         SELECT t.*, r.name AS chrom, r.sequence_length
         FROM duckvep_transcripts t
         JOIN duckvep_sequence_regions r USING (seq_region)
         WHERE t.transcript_start > 1
           AND t.transcript_end <= r.sequence_length
           {generated_chrom_filter}
       ), exons AS MATERIALIZED (
         SELECT
           e.*, t.chrom, t.sequence_length, t.transcript_start,
           t.transcript_end, t.strand, t.cds_start, t.cds_end,
           lead(e.exon_start) OVER (
             PARTITION BY e.transcript_index ORDER BY e.exon_start, e.exon_end
           ) AS next_exon_start
         FROM duckvep_exons e
         JOIN transcripts t USING (transcript_index)
       ), coding_segments AS MATERIALIZED (
         SELECT
           transcript_index, chrom, sequence_length, strand,
           greatest(exon_start, cds_start) AS segment_start,
           least(exon_end, cds_end) AS segment_end
         FROM exons
         WHERE cds_start IS NOT NULL AND cds_end IS NOT NULL
           AND greatest(exon_start, cds_start) <= least(exon_end, cds_end)
       ), introns AS MATERIALIZED (
         SELECT
           transcript_index, chrom, sequence_length,
           exon_end + 1 AS intron_start, next_exon_start - 1 AS intron_end
         FROM exons
         WHERE next_exon_start IS NOT NULL AND next_exon_start > exon_end + 1
       ), span_geometries AS (
         SELECT
           'transcript_exact' AS event_state, transcript_index, chrom,
           transcript_start AS event_start, transcript_end AS event_end
         FROM transcripts
         UNION ALL
         SELECT
           'transcript_containing', transcript_index, chrom,
           CASE WHEN transcript_start > 18 THEN transcript_start - 17 ELSE 2 END,
           least(sequence_length, transcript_end + 17)
         FROM transcripts
         UNION ALL
         SELECT
           'left_partial', transcript_index, chrom,
           CASE WHEN transcript_start > 18 THEN transcript_start - 17 ELSE 2 END,
           least(transcript_end, transcript_start + 17)
         FROM transcripts
         UNION ALL
         SELECT
           'right_partial', transcript_index, chrom,
           greatest(transcript_start,
                    CASE WHEN transcript_end > 17 THEN transcript_end - 17 ELSE 2 END),
           least(sequence_length, transcript_end + 17)
         FROM transcripts
         UNION ALL
         SELECT
           'start_codon', transcript_index, chrom,
           CASE WHEN strand = 1 THEN cds_start ELSE greatest(cds_start, cds_end - 2) END,
           CASE WHEN strand = 1 THEN least(cds_end, cds_start + 2) ELSE cds_end END
         FROM transcripts
         WHERE cds_start IS NOT NULL AND cds_end IS NOT NULL
         UNION ALL
         SELECT
           'stop_codon', transcript_index, chrom,
           CASE WHEN strand = 1 THEN greatest(cds_start, cds_end - 2) ELSE cds_start END,
           CASE WHEN strand = 1 THEN cds_end ELSE least(cds_end, cds_start + 2) END
         FROM transcripts
         WHERE cds_start IS NOT NULL AND cds_end IS NOT NULL
         UNION ALL
         SELECT
           concat('coding_len', CAST(length_value AS VARCHAR)),
           transcript_index, chrom,
           segment_start + (
             hash(transcript_index, length_value, {opt$seed}) %
             (segment_end - segment_start - length_value + 2)
           ) AS event_start,
           segment_start + (
             hash(transcript_index, length_value, {opt$seed}) %
             (segment_end - segment_start - length_value + 2)
           ) + length_value - 1 AS event_end
         FROM coding_segments
         CROSS JOIN (VALUES (1), (2), (3), (4), (9), (10)) lengths(length_value)
         WHERE segment_end - segment_start + 1 >= length_value
         UNION ALL
         SELECT
           concat('intron_len', CAST(length_value AS VARCHAR)),
           transcript_index, chrom,
           intron_start + (
             hash(transcript_index, length_value, {opt$seed}, 1) %
             (intron_end - intron_start - length_value + 2)
           ) AS event_start,
           intron_start + (
             hash(transcript_index, length_value, {opt$seed}, 1) %
             (intron_end - intron_start - length_value + 2)
           ) + length_value - 1 AS event_end
         FROM introns
         CROSS JOIN (VALUES (1), (3), (10)) lengths(length_value)
         WHERE intron_end - intron_start + 1 >= length_value
         UNION ALL
         SELECT
           'splice_cross', transcript_index, chrom,
           greatest(2, exon_end - 1), least(sequence_length, next_exon_start + 1)
         FROM exons
         WHERE next_exon_start IS NOT NULL AND next_exon_start > exon_end + 1
         UNION ALL
         SELECT
           CASE WHEN strand = 1 THEN 'five_prime_utr' ELSE 'three_prime_utr' END,
           transcript_index, chrom, exon_start, least(exon_end, cds_start - 1)
         FROM exons
         WHERE cds_start IS NOT NULL AND exon_start < cds_start
           AND exon_start <= least(exon_end, cds_start - 1)
         UNION ALL
         SELECT
           CASE WHEN strand = 1 THEN 'three_prime_utr' ELSE 'five_prime_utr' END,
           transcript_index, chrom, greatest(exon_start, cds_end + 1), exon_end
         FROM exons
         WHERE cds_end IS NOT NULL AND exon_end > cds_end
           AND greatest(exon_start, cds_end + 1) <= exon_end
       ), insertion_geometries AS (
         SELECT
           'insertion_coding' AS event_state, transcript_index, chrom,
           segment_start + (
             hash(transcript_index, {opt$seed}, 2) %
             (segment_end - segment_start + 1)
           ) AS event_start,
           segment_start + (
             hash(transcript_index, {opt$seed}, 2) %
             (segment_end - segment_start + 1)
           ) AS event_end
         FROM coding_segments
         UNION ALL
         SELECT
           'insertion_start_left', transcript_index, chrom,
           greatest(1, CASE WHEN strand = 1 THEN cds_start - 1 ELSE cds_end - 1 END),
           greatest(1, CASE WHEN strand = 1 THEN cds_start - 1 ELSE cds_end - 1 END)
         FROM transcripts
         WHERE cds_start IS NOT NULL AND cds_end IS NOT NULL
         UNION ALL
         SELECT
           'insertion_start_right', transcript_index, chrom,
           CASE WHEN strand = 1 THEN cds_start ELSE cds_end END,
           CASE WHEN strand = 1 THEN cds_start ELSE cds_end END
         FROM transcripts
         WHERE cds_start IS NOT NULL AND cds_end IS NOT NULL
         UNION ALL
         SELECT
           'insertion_exon_edge', transcript_index, chrom, exon_end, exon_end
         FROM exons
         WHERE next_exon_start IS NOT NULL
         UNION ALL
         SELECT
           'insertion_intron', transcript_index, chrom,
           intron_start + (
             hash(transcript_index, {opt$seed}, 3) %
             (intron_end - intron_start + 1)
           ) AS event_start,
           intron_start + (
             hash(transcript_index, {opt$seed}, 3) %
             (intron_end - intron_start + 1)
           ) AS event_end
         FROM introns
       ), geometries AS (
         SELECT false AS is_insertion, * FROM span_geometries
         UNION ALL
         SELECT true AS is_insertion, * FROM insertion_geometries
       ), ranked_geometries AS (
         SELECT *, row_number() OVER (
           PARTITION BY is_insertion, event_state
           ORDER BY hash(
             chrom, event_start, event_end, transcript_index, {opt$seed}
           )
         ) AS geometry_rank
         FROM geometries
         WHERE event_start > 1 AND event_end >= event_start
           AND event_end < 4294967295
       ), selected_geometries AS (
         SELECT * FROM ranked_geometries {generated_filter}
       ), span_events AS (
         SELECT
           concat(
             'duckvep-generated:', event_state, ':', transcript_index, ':',
             event_start, ':', event_end, ':', structural_type
           ) AS source_id,
           chrom, event_start - 1 AS raw_position, event_end AS raw_end,
           'N' AS reference, alternate, source_svtype, structural_type,
           event_state
         FROM selected_geometries
         CROSS JOIN (VALUES
           ('DEL', 'DEL', '<DEL>'),
           ('DUP', 'DUP', '<DUP>'),
           ('TDUP', 'DUP', '<DUP:TANDEM>'),
           ('INV', 'INV', '<INV>'),
           ('CNV', 'CNV', '<CNV>')
         ) operations(structural_type, source_svtype, alternate)
         WHERE NOT is_insertion
       ), insertion_events AS (
         SELECT
           concat(
             'duckvep-generated:', event_state, ':', transcript_index, ':',
             event_start, ':', event_end, ':INS'
           ) AS source_id,
           chrom, event_start AS raw_position, event_end AS raw_end,
           'N' AS reference, '<INS>' AS alternate, 'INS' AS source_svtype,
           'INS' AS structural_type, event_state
         FROM selected_geometries WHERE is_insertion
       )
       SELECT * FROM span_events
       UNION ALL BY NAME SELECT * FROM insertion_events"
    )
  ))
}

structural_source_sql <- if (generate_structural) {
  "SELECT * FROM duckvep_generated_structural_source"
} else if (identical(opt$event_mode, "structural")) {
  glue(
    "SELECT DISTINCT
       coalesce(nullif(ID, ''), '.') AS source_id,
       CHROM AS chrom,
       POS::UBIGINT AS raw_position,
       coalesce(INFO_END::UBIGINT, POS::UBIGINT) AS raw_end,
       REF AS reference,
       ALT[1] AS alternate,
       upper(INFO_SVTYPE) AS source_svtype,
       CASE
         WHEN upper(ALT[1]) = '<DUP:TANDEM>' THEN 'TDUP'
         WHEN upper(INFO_SVTYPE) = 'DEL' THEN 'DEL'
         WHEN upper(INFO_SVTYPE) = 'DUP' THEN 'DUP'
         WHEN upper(INFO_SVTYPE) = 'INV' THEN 'INV'
         WHEN upper(INFO_SVTYPE) = 'INS' THEN 'INS'
         WHEN upper(INFO_SVTYPE) = 'CNV' THEN 'CNV'
         ELSE NULL
       END AS structural_type,
       'source' AS event_state
     FROM read_bcf({sql_q(normalizePath(source_vcf))}, scan_mode := 'sequential')
     WHERE len(ALT) = 1
       AND regexp_full_match(ALT[1], '<[^>]+>')
       AND upper(INFO_SVTYPE) IN ('DEL', 'DUP', 'INV', 'INS', 'CNV')
       {chrom_filter}"
  )
} else {
  ""
}

sample_sql <- if (identical(opt$event_mode, "small")) {
  glue(
    "CREATE OR REPLACE TEMP TABLE duckvep_sample AS
     WITH alleles AS (
       SELECT DISTINCT
         CHROM AS chrom,
         POS::UBIGINT AS position,
         REF AS reference,
         ALT[1] AS alternate
       FROM read_bcf({sql_q(normalizePath(source_vcf))}, scan_mode := 'sequential')
       WHERE len(ALT) = 1
         AND REF <> ALT[1]
         AND length(REF) BETWEEN 1 AND 50
         AND length(ALT[1]) BETWEEN 1 AND 50
         AND regexp_full_match(REF, '[ACGT]+')
         AND regexp_full_match(ALT[1], '[ACGT]+')
         {chrom_filter}
     ), shaped AS (
       SELECT *,
         CASE
           WHEN length(reference) = 1 AND length(alternate) = 1 THEN 'snv'
           WHEN length(reference) = length(alternate) THEN 'mnv'
           WHEN length(reference) < length(alternate) THEN 'insertion_like'
           ELSE 'deletion_like'
         END AS var_type,
         length(alternate) - length(reference) AS length_change
       FROM alleles
     ), binned AS (
       SELECT *,
         CASE
           WHEN length_change = 0 THEN '0'
           WHEN abs(length_change) = 1 THEN concat(if(length_change > 0, '+', '-'), '1')
           WHEN abs(length_change) <= 3 THEN concat(if(length_change > 0, '+', '-'), '2..3')
           WHEN abs(length_change) <= 10 THEN concat(if(length_change > 0, '+', '-'), '4..10')
           ELSE concat(if(length_change > 0, '+', '-'), '11..50')
         END AS length_bin
       FROM shaped
     ), ranked AS (
       SELECT *, row_number() OVER (
         PARTITION BY var_type, length_bin
         ORDER BY hash(chrom, position, reference, alternate, {opt$seed})
       ) AS sample_rank
       FROM binned
     )
     SELECT
       concat('duckvep:', chrom, ':', position, ':', reference, ':', alternate) AS variant_id,
       r.seq_region,
       chrom,
       position,
       reference,
       alternate,
       var_type,
       length_bin
     FROM ranked
     JOIN duckvep_sequence_regions r ON r.name = chrom
     {sample_filter}
     ORDER BY r.seq_region, position, reference, alternate"
  )
} else {
  glue(
    "CREATE OR REPLACE TEMP TABLE duckvep_sample AS
     WITH source_events AS (
       {structural_source_sql}
     ), prepared AS (
       SELECT *,
         CASE WHEN structural_type = 'INS'
              THEN raw_position ELSE raw_position + 1 END AS event_start,
         CASE WHEN structural_type = 'INS'
              THEN raw_position ELSE raw_end END AS event_end,
         CASE
           WHEN structural_type = 'DEL' THEN 'LOSS'
           WHEN structural_type IN ('DUP', 'TDUP') THEN 'GAIN'
           WHEN structural_type = 'INV' THEN 'NEUTRAL'
           ELSE 'UNKNOWN'
         END AS copy_change
       FROM source_events
       WHERE raw_position < 4294967295
     ), size_binned AS (
       SELECT *, lower(structural_type) AS var_type,
         CASE
           WHEN structural_type = 'INS' THEN '0'
           WHEN event_end - event_start + 1 = 1 THEN '1'
           WHEN event_end - event_start + 1 <= 10 THEN '2..10'
           WHEN event_end - event_start + 1 <= 100 THEN '11..100'
           WHEN event_end - event_start + 1 <= 1000 THEN '101..1000'
           WHEN event_end - event_start + 1 <= 10000 THEN '1001..10000'
           ELSE '>10000'
         END AS length_class
       FROM prepared
       WHERE event_start > 0 AND event_end >= event_start
     ), binned AS (
       SELECT *, concat(event_state, '/', length_class) AS length_bin
       FROM size_binned
     ), ranked AS (
       SELECT *, row_number() OVER (
         PARTITION BY structural_type, copy_change, event_state, length_class
         ORDER BY hash(
           chrom, raw_position, raw_end, alternate, structural_type,
           source_id, {opt$seed}
         )
       ) AS sample_rank
       FROM binned
     )
     SELECT
       CASE
         WHEN source_id LIKE 'duckvep-generated:%' THEN source_id
         ELSE concat(
           'duckvep-sv:', chrom, ':', raw_position, ':', raw_end, ':',
           structural_type, ':', alternate
         )
       END AS variant_id,
       r.seq_region,
       chrom,
       CAST(event_start AS UBIGINT) AS position,
       CAST(event_start AS UBIGINT) AS event_start,
       CAST(event_end AS UBIGINT) AS event_end,
       raw_position,
       raw_end,
       reference,
       alternate,
       structural_type,
       copy_change,
       source_svtype,
       var_type,
       length_bin
     FROM ranked
     JOIN duckvep_sequence_regions r ON r.name = chrom
     {sample_filter}
     ORDER BY r.seq_region, event_start, event_end, structural_type"
  )
}
invisible(dbExecute(con, sample_sql))
sample_count <- dbGetQuery(con, "SELECT count(*) AS n FROM duckvep_sample")$n[[
  1L
]]
if (sample_count == 0) {
  die("the sampler found no supported {opt$event_mode} events in model regions")
}

sample_vcf <- opt$sample_vcf
if (!nzchar(sample_vcf)) {
  sample_vcf <- tempfile(fileext = ".vcf")
  if (!isTRUE(opt$keep_sample_vcf)) {
    temporary_files <- c(temporary_files, sample_vcf)
  }
}
dir.create(dirname(sample_vcf), recursive = TRUE, showWarnings = FALSE)
vc <- file(sample_vcf, open = "wt")
vcf_header <- "##fileformat=VCFv4.2"
if (identical(opt$event_mode, "structural")) {
  vcf_header <- c(
    vcf_header,
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">"
  )
}
writeLines(
  c(vcf_header, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"),
  vc
)
sample_vcf_query <- if (identical(opt$event_mode, "small")) {
  paste(
    "SELECT chrom, position, variant_id, reference, alternate",
    "FROM duckvep_sample ORDER BY seq_region, position, reference, alternate"
  )
} else {
  paste(
    "SELECT chrom, raw_position AS position, raw_end, source_svtype,",
    "variant_id, reference, alternate FROM duckvep_sample",
    "ORDER BY seq_region, event_start, event_end, structural_type"
  )
}
res <- dbSendQuery(
  con,
  sample_vcf_query
)
repeat {
  chunk <- dbFetch(res, n = 100000L)
  if (nrow(chunk) == 0L) {
    break
  }
  info <- if (identical(opt$event_mode, "small")) {
    rep(".", nrow(chunk))
  } else {
    paste0("END=", chunk$raw_end, ";SVTYPE=", chunk$source_svtype)
  }
  writeLines(
    paste(
      chunk$chrom, chunk$position, chunk$variant_id, chunk$reference,
      chunk$alternate, ".", "PASS", info, sep = "\t"
    ),
    vc
  )
}
dbClearResult(res)
close(vc)

stage_gff <- function(path) {
  if (grepl("[.]gz$", path) && file.exists(paste0(path, ".tbi"))) {
    return(path)
  }
  header <- tempfile(fileext = ".header")
  body <- tempfile(fileext = ".body")
  sorted_body <- tempfile(fileext = ".sorted-body")
  sorted_gff <- tempfile(fileext = ".gff")
  output <- tempfile(fileext = ".gff.gz")
  temporary_files <<- c(
    temporary_files,
    header,
    body,
    sorted_body,
    sorted_gff,
    output,
    paste0(output, ".tbi")
  )
  input <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
  header_con <- file(header, "wt")
  body_con <- file(body, "wt")
  repeat {
    lines <- readLines(input, n = 100000L, warn = FALSE)
    if (length(lines) == 0L) {
      break
    }
    is_header <- startsWith(lines, "#")
    writeLines(lines[is_header], header_con)
    writeLines(lines[!is_header], body_con)
  }
  close(input)
  close(header_con)
  close(body_con)
  rc <- system2(
    "sort",
    c("-k1,1", "-k4,4n", body),
    stdout = sorted_body,
    env = "LC_ALL=C"
  )
  if (rc != 0L) {
    die("sorting GFF failed")
  }
  out_con <- file(sorted_gff, "wt")
  for (part in c(header, sorted_body)) {
    in_con <- file(part, "rt")
    repeat {
      lines <- readLines(in_con, n = 100000L, warn = FALSE)
      if (length(lines) == 0L) {
        break
      }
      writeLines(lines, out_con)
    }
    close(in_con)
  }
  close(out_con)
  bgzip <- dbGetQuery(
    con,
    glue(
      "SELECT success FROM bgzip({sql_q(sorted_gff)}, output_path := {sql_q(output)},
       keep := TRUE, overwrite := TRUE)"
    )
  )$success
  if (length(bgzip) != 1L || !isTRUE(bgzip[[1L]])) {
    die("bgzip of GFF failed")
  }
  indexed <- dbGetQuery(
    con,
    glue(
      "SELECT success FROM tabix_index({sql_q(output)}, preset := 'gff',
       index_path := {sql_q(paste0(output, '.tbi'))}, threads := 1)"
    )
  )$success
  if (length(indexed) != 1L || !isTRUE(indexed[[1L]])) {
    die("tabix indexing of GFF failed")
  }
  output
}
gff_for_vep <- if (identical(oracle_mode, "gff")) {
  stage_gff(opt$gff)
} else {
  ""
}

engine_call_sql <- if (identical(opt$event_mode, "small")) {
  glue(
    "duckvep_annotate(
       {sql_q(opt$model_name)}, v.seq_region, v.position,
       v.reference, v.alternate, {opt$distance}::UBIGINT
     )"
  )
} else {
  glue(
    "duckvep_annotate_sv(
       {sql_q(opt$model_name)}, v.seq_region, v.event_start, v.event_end,
       v.structural_type, v.copy_change, {opt$distance}::UBIGINT
     )"
  )
}

engine_time <- system.time({
  engine_nmd_sql <- if (nmd_oracle_enabled) {
    "coalesce(v.annotation.nmd_prediction, 'not_applicable')"
  } else {
    "'not_measured'::VARCHAR"
  }
  dbExecute(
    con,
    glue(
      "CREATE OR REPLACE TEMP TABLE duckvep_annotation AS
       WITH annotated AS (
         SELECT
           v.variant_id,
           v.seq_region,
           v.position,
           unnest({engine_call_sql}) AS annotation
         FROM duckvep_sample v
       )
       SELECT
         v.variant_id,
         n.transcript_id AS tx,
         list_aggregate(
           list_sort(list_distinct(string_split(v.annotation.consequence, '&'))),
           'string_agg', '&'
         ) AS consequence,
         v.annotation.impact,
         v.annotation.status,
         v.annotation.reason,
         {engine_nmd_sql} AS nmd_prediction
       FROM annotated v
       JOIN duckvep_transcript_names n
         ON n.transcript_index = v.annotation.transcript_index
       ORDER BY v.seq_region, v.position, n.transcript_index"
    )
  )
})

micromamba <- if (file.exists(opt$micromamba)) {
  normalizePath(opt$micromamba)
} else {
  unname(Sys.which(opt$micromamba))
}
if (!nzchar(micromamba)) {
  die("micromamba is unavailable: {opt$micromamba}")
}
if (!dir.exists(opt$vep_prefix)) {
  die("VEP environment prefix does not exist: {opt$vep_prefix}")
}
vep_prefix <- normalizePath(opt$vep_prefix)
vep_command <- function(...) {
  blit::conda(
    "run",
    "-p",
    vep_prefix,
    "vep",
    ...,
    conda = micromamba
  )
}

vep_info <- tempfile(fileext = ".txt")
temporary_files <- c(temporary_files, vep_info)
info_rc <- vep_command("--help") |>
  blit::cmd_run(
    stdout = vep_info,
    stderr = "2>&1",
    stdin = NULL,
    verbose = FALSE
  )
if (info_rc != 0L || !file.exists(vep_info)) {
  die("cannot identify VEP in {vep_prefix}")
}
vep_info_lines <- readLines(vep_info, warn = FALSE)
component_version <- function(component) {
  line <- grep(
    glue("^\\s*{component}\\s*:"),
    vep_info_lines,
    value = TRUE
  )
  if (length(line) != 1L) {
    die("VEP did not report {component} exactly once")
  }
  trimws(sub("^[^:]+:", "", line))
}
oracle_version <- component_version("ensembl-vep")
if (!identical(oracle_version, "116.0")) {
  die("expected Ensembl VEP 116.0, found {oracle_version}")
}
oracle_details <- if (identical(oracle_mode, "cache")) {
  details <- c(
    "oracle=cache",
    glue("assembly={opt$assembly}"),
    glue("species={opt$species}")
  )
  if (nzchar(opt$cache_version)) {
    details <- c(details, glue("cache_version={opt$cache_version}"))
  }
  details
} else {
  "oracle=gff"
}
oracle_build <- paste(
  c(
    glue("core={component_version('ensembl')}"),
    glue("variation={component_version('ensembl-variation')}"),
    glue("vep={oracle_version}"),
    oracle_details,
    if (nmd_oracle_enabled) {
      c(
        glue("plugin=NMD.pm;plugin_sha256={nmd_plugin_sha256}"),
        glue("state_plugin_sha256={nmd_state_plugin_sha256}")
      )
    }
  ),
  collapse = ";"
)

vep_json <- tempfile(fileext = ".json")
temporary_files <- c(temporary_files, vep_json)
vep_args <- c(
  "-i",
  sample_vcf,
  "--fasta",
  opt$fasta,
  "--distance",
  as.character(opt$distance),
  "--json",
  "-o",
  vep_json,
  "--fork",
  opt$fork,
  "--force_overwrite",
  "--no_stats"
)
if (identical(oracle_mode, "cache")) {
  vep_args <- c(
    vep_args,
    "--cache",
    "--offline",
    "--dir_cache",
    normalizePath(opt$cache_dir),
    "--assembly",
    opt$assembly,
    "--species",
    opt$species
  )
  if (nzchar(opt$cache_version)) {
    vep_args <- c(vep_args, "--cache_version", opt$cache_version)
  }
} else {
  vep_args <- c(vep_args, "--gff", gff_for_vep)
}
if (nmd_oracle_enabled) {
  vep_args <- c(
    vep_args,
    "--plugin",
    "NMD",
    "--plugin",
    "DuckVEPNMDState",
    "--dir_plugins",
    nmd_plugin_run_dir
  )
}
vep_time <- system.time({
  rc <- do.call(vep_command, as.list(vep_args)) |>
    blit::cmd_run(stdout = "", stderr = "", stdin = NULL, verbose = FALSE)
})
if (rc != 0L || !file.exists(vep_json) || file.info(vep_json)$size == 0) {
  die("VEP failed with exit status {rc}")
}

vep_nmd_sql <- if (nmd_oracle_enabled) {
  "CASE
     WHEN NOT (
       list_contains(tc.consequence_terms, 'stop_gained') OR
       list_contains(tc.consequence_terms, 'frameshift_variant') OR
       list_contains(tc.consequence_terms, 'splice_donor_variant') OR
       list_contains(tc.consequence_terms, 'splice_acceptor_variant')
     ) THEN 'not_applicable'
     WHEN coalesce(
            json_extract_string(to_json(tc), '$.duckvep_nmd_cds'),
            'undefined'
          ) = 'undefined' THEN 'unresolved'
     WHEN json_extract_string(to_json(tc), '$.nmd') =
          'NMD_escaping_variant' THEN 'escaping'
     ELSE 'triggering'
   END"
} else {
  "'not_measured'::VARCHAR"
}
invisible(dbExecute(
  con,
  glue(
    "CREATE OR REPLACE TEMP TABLE vep_annotation AS
     SELECT
       j.id AS variant_id,
       tc.transcript_id AS tx,
       coalesce(
         list_aggregate(list_sort(list_distinct(tc.consequence_terms)), 'string_agg', '&'),
         ''
       ) AS consequence,
       coalesce(tc.impact, '') AS impact,
       'oracle'::VARCHAR AS status,
       NULL::VARCHAR AS reason,
       {vep_nmd_sql} AS nmd_prediction
     FROM read_json({sql_q(vep_json)}, format = 'newline_delimited', sample_size = -1) j,
     UNNEST(j.transcript_consequences) u(tc)"
  )
))

run_date <- as.character(Sys.Date())
invisible(dbExecute(
  con,
  glue(
    "CREATE OR REPLACE TEMP TABLE duckvep_annotation_dump AS
     SELECT
       {sql_q(run_date)} AS run_date,
       {sql_q(opt$corpus)} AS corpus,
       {sql_q(opt$model_name)} AS model,
       {sql_q(oracle_version)} AS oracle_version,
       {sql_q(oracle_build)} AS oracle_build,
       'vep'::VARCHAR AS source,
       a.variant_id,
       v.chrom,
       v.position AS pos,
       v.reference AS ref,
       v.alternate AS alt,
       v.var_type,
       v.length_bin,
       a.tx,
       a.consequence,
       a.impact,
       a.status,
       a.reason,
       a.nmd_prediction
     FROM vep_annotation a JOIN duckvep_sample v USING (variant_id)
     UNION ALL
     SELECT
       {sql_q(run_date)}, {sql_q(opt$corpus)}, {sql_q(opt$model_name)},
       {sql_q(oracle_version)}, {sql_q(oracle_build)}, 'duckvep',
       a.variant_id, v.chrom, v.position, v.reference, v.alternate,
       v.var_type, v.length_bin, a.tx, a.consequence, a.impact, a.status, a.reason,
       a.nmd_prediction
     FROM duckvep_annotation a JOIN duckvep_sample v USING (variant_id)"
  )
))
dir.create(dirname(opt$annotations_out), recursive = TRUE, showWarnings = FALSE)
invisible(dbExecute(
  con,
  glue(
    "COPY duckvep_annotation_dump TO {sql_q(opt$annotations_out)}
     (FORMAT parquet, COMPRESSION zstd)"
  )
))

counts <- dbGetQuery(
  con,
  "SELECT source, count(*) AS transcript_rows
   FROM duckvep_annotation_dump GROUP BY source ORDER BY source"
)
invisible(dbGetQuery(
  con,
  glue("SELECT duckvep_model_drop({sql_q(opt$model_name)})")
))

cat(glue("sampled variants: {sample_count}"), "\n", sep = "")
for (i in seq_len(nrow(counts))) {
  cat(
    glue("{counts$source[i]} transcript rows: {counts$transcript_rows[i]}"),
    "\n",
    sep = ""
  )
}
cat(
  glue("DuckVEP annotation: {sprintf('%.3f', engine_time[['elapsed']])} s"),
  "\n",
  sep = ""
)
cat(
  glue(
    "VEP {oracle_mode} annotation: ",
    "{sprintf('%.3f', vep_time[['elapsed']])} s"
  ),
  "\n",
  sep = ""
)
cat(glue("annotations: {opt$annotations_out}"), "\n", sep = "")

report <- file.path(
  root,
  "test",
  "duckvep",
  "conformance",
  "statistical_conformance.R"
)
rc <- system2("Rscript", c(report, "--annotations", opt$annotations_out))
if (rc != 0L) {
  die("statistical report failed with exit status {rc}")
}
