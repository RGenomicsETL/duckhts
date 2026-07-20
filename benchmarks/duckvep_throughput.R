#!/usr/bin/env Rscript

# Record the stable DuckDB C-API path over a sorted variant stream. With no
# database this uses the checked-in one-transcript fixture. --database expects
# pre-staged bench_regions, bench_transcripts, and bench_exons. Variants may
# reside in the model database or in a separately staged read-only database.
# Complete production models expose sequence_length plus full transcript flanks;
# --regulatory additionally loads a prepared regulation-feature relation. Older
# topology-only staging databases remain accepted as partial models.

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
op <- add_option(op, "--database", default = "")
op <- add_option(
  op,
  "--variants-database",
  dest = "variants_database",
  default = "",
  help = "optional database containing --variants-table [%default]"
)
op <- add_option(
  op,
  "--variants-table",
  dest = "variants_table",
  default = "bench_variants",
  help = "prepared variant relation in --database [%default]"
)
op <- add_option(
  op,
  "--event-mode",
  dest = "event_mode",
  default = "small",
  help = "prepared event shape: small or breakend [%default]"
)
op <- add_option(op, "--output", default = "rich")
op <- add_option(
  op,
  "--reference-fasta",
  dest = "reference_fasta",
  default = "",
  help = "indexed FASTA required by --output hgvs [%default]"
)
op <- add_option(
  op,
  "--regulatory",
  action = "store_true",
  default = FALSE,
  help = "load bench_regulation_features or duckvep_regulation_features"
)
op <- add_option(op, "--workload-name", dest = "workload_name", default = "")
op <- add_option(op, "--variants", type = "double", default = 10000000)
op <- add_option(op, "--passes", type = "integer", default = 3L)
op <- add_option(op, "--warmup", type = "double", default = 100000)
op <- add_option(op, "--threads", type = "integer", default = 1L)
op <- add_option(
  op,
  "--input-partitions",
  dest = "input_partitions",
  type = "integer",
  default = 1L,
  help = paste(
    "number of disjoint ordered input partitions exposed to DuckDB",
    "(1..1024) [%default]"
  )
)
op <- add_option(
  op,
  "--transcript-distance",
  dest = "transcript_distance",
  type = "double",
  default = 5000,
  help = "symmetric upstream/downstream transcript window in bases [%default]"
)
op <- add_option(
  op,
  "--history",
  default = file.path(root, "benchmarks", "data", "duckvep_throughput.csv")
)
op <- add_option(
  op,
  "--composition",
  default = "",
  help = "optional untimed compact SO/region/object composition CSV [%default]"
)
op <- add_option(
  op,
  "--fingerprint",
  default = "",
  help = "optional untimed full-public-row fingerprint CSV [%default]"
)
op <- add_option(
  op,
  "--expected-fingerprint",
  dest = "expected_fingerprint",
  default = "",
  help = "optional checked fingerprint baseline to verify [%default]"
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
whole_nonnegative <- function(x, name, maximum = Inf) {
  if (!is.finite(x) || x < 0 || x != floor(x) || x > maximum) {
    die("{name} must be an integer from zero through {maximum}")
  }
  format(x, scientific = FALSE, trim = TRUE)
}
variant_sql <- whole_count(opt$variants, "--variants")
distance_sql <- whole_nonnegative(
  opt$transcript_distance,
  "--transcript-distance",
  2^32 - 1
)
if (opt$passes < 1L) {
  die("--passes must be positive")
}
if (opt$threads < 1L) {
  die("--threads must be positive")
}
if (opt$input_partitions < 1L || opt$input_partitions > 1024L) {
  die("--input-partitions must be from 1 through 1024")
}
if (!opt$output %in% c("rich", "compact", "hgvs")) {
  die("--output must be rich, compact, or hgvs")
}
if (!opt$event_mode %in% c("small", "breakend")) {
  die("--event-mode must be small or breakend")
}
production <- nzchar(opt$database)
if (isTRUE(opt$regulatory) && !production) {
  die("--regulatory requires a production --database")
}
if (identical(opt$output, "hgvs") && !identical(opt$event_mode, "small")) {
  die("--output hgvs currently requires --event-mode small")
}
if (identical(opt$output, "hgvs") && isTRUE(opt$regulatory)) {
  die("--output hgvs cannot be combined with --regulatory")
}
if (identical(opt$output, "hgvs") && !nzchar(opt$reference_fasta)) {
  if (production) {
    die("--output hgvs with --database requires --reference-fasta")
  }
  opt$reference_fasta <- file.path(
    root,
    "test",
    "data",
    "duckvep",
    "minimal.fa"
  )
}
sha256_file <- function(path) {
  value <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(value, "status")
  if (!is.null(status) && status != 0L) die("sha256sum failed for {path}")
  strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
inputs <- c(
  opt$extension,
  if (production) opt$database else opt$model_sql,
  if (nzchar(opt$variants_database)) opt$variants_database,
  if (nzchar(opt$reference_fasta)) opt$reference_fasta
)
if (nzchar(opt$expected_fingerprint)) {
  inputs <- c(inputs, opt$expected_fingerprint)
}
missing <- inputs[!file.exists(inputs)]
if (length(missing) != 0L) {
  die("missing input(s):\n{paste(missing, collapse = '\n')}")
}
extension_sha256 <- sha256_file(normalizePath(opt$extension))
model_database_sha256 <- if (production) {
  sha256_file(normalizePath(opt$database))
} else {
  ""
}
variants_database_sha256 <- if (nzchar(opt$variants_database)) {
  sha256_file(normalizePath(opt$variants_database))
} else {
  ""
}

drv <- duckdb(config = list(allow_unsigned_extensions = "true"))
con <- if (production) {
  dbConnect(drv, dbdir = normalizePath(opt$database))
} else {
  dbConnect(drv)
}
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
sql_q <- function(x) as.character(dbQuoteString(con, x))
identifier_q <- function(x) as.character(dbQuoteIdentifier(con, x))
invisible(dbExecute(con, glue("SET threads = {opt$threads}")))
invisible(dbExecute(con, glue("LOAD {sql_q(normalizePath(opt$extension))}")))

variants_catalog <- "duckvep_benchmark_variants"
corpus_receipt <- NULL
model_logical_sha256 <- ""
model_assembly <- ""
model_source_version <- ""
model_region_ordinal_sha256 <- ""
if (nzchar(opt$variants_database)) {
  if (!production) {
    die("--variants-database requires --database")
  }
  invisible(dbExecute(
    con,
    glue(
      "ATTACH {sql_q(normalizePath(opt$variants_database))} AS ",
      "{identifier_q(variants_catalog)} (READ_ONLY)"
    )
  ))
}

if (production) {
  main_relations <- dbGetQuery(
    con,
    "SELECT table_name FROM information_schema.tables
     WHERE table_catalog = current_database() AND table_schema = 'main'"
  )$table_name
  if ("model_receipt" %in% main_relations) {
    model_receipt <- dbGetQuery(
      con,
      "SELECT source_version, assembly, model_sha256 FROM model_receipt"
    )
    if (nrow(model_receipt) != 1L) {
      die("model_receipt must contain exactly one row")
    }
    model_source_version <- as.character(model_receipt$source_version[[1L]])
    model_assembly <- as.character(model_receipt$assembly[[1L]])
    model_logical_sha256 <- as.character(model_receipt$model_sha256[[1L]])
  }
}

if (nzchar(opt$variants_database)) {
  sidecar_relations <- dbGetQuery(
    con,
    glue(
      "SELECT table_name FROM information_schema.tables
       WHERE table_catalog = {sql_q(variants_catalog)}
         AND table_schema = 'main'"
    )
  )$table_name
  has_dense_receipt <- "dense_corpus_receipt" %in% sidecar_relations
  if (startsWith(opt$variants_table, "bench_dense_") && !has_dense_receipt) {
    die("dense variant table requires dense_corpus_receipt")
  }
  if (has_dense_receipt) {
    corpus_receipt <- dbGetQuery(
      con,
      glue(
        "SELECT * FROM {identifier_q(variants_catalog)}.main.dense_corpus_receipt"
      )
    )
    if (nrow(corpus_receipt) != 1L ||
        corpus_receipt$schema_version[[1L]] != 2L) {
      die("dense_corpus_receipt must contain one schema-version-2 row")
    }
    if (!nzchar(model_logical_sha256)) {
      die("dense corpus requires a model_receipt")
    }
    model_region_ordinal_sha256 <- dbGetQuery(
      con,
      "SELECT sha256(string_agg(
         seq_region::VARCHAR || ':' || seq_region_name || ':' ||
         sequence_length::VARCHAR,
         '|' ORDER BY seq_region
       )) AS digest
       FROM model_regions"
    )$digest[[1L]]
    if (!identical(
          as.character(corpus_receipt$model_sha256[[1L]]),
          model_logical_sha256
        ) ||
        !identical(
          as.character(corpus_receipt$assembly[[1L]]),
          model_assembly
        ) ||
        !identical(
          as.character(corpus_receipt$model_source_version[[1L]]),
          model_source_version
        ) ||
        !identical(
          as.character(corpus_receipt$region_ordinal_sha256[[1L]]),
          as.character(model_region_ordinal_sha256)
        )) {
      die("dense corpus receipt does not match the benchmark model")
    }
  }
}

model_name <- "throughput"
if (production) {
  variants_table <- if (nzchar(opt$variants_database)) {
    paste(
      identifier_q(variants_catalog),
      identifier_q("main"),
      identifier_q(opt$variants_table),
      sep = "."
    )
  } else {
    identifier_q(opt$variants_table)
  }
  variants_columns <- dbGetQuery(
    con,
    glue("DESCRIBE SELECT * FROM {variants_table}")
  )$column_name
  stable_source_order <- all(
    c("source_record_index", "source_alt_index") %in% variants_columns
  )
  event_projection <- if (identical(opt$event_mode, "breakend")) {
    "seq_region, \"position\", mate_seq_region, mate_position"
  } else {
    "seq_region, \"position\", \"reference\", \"alternate\""
  }
  input_order <- if (identical(opt$event_mode, "breakend")) {
    "seq_region, \"position\", mate_seq_region, mate_position"
  } else {
    paste(
      "seq_region, \"position\", \"reference\", \"alternate\"",
      if (stable_source_order) {
        ", source_record_index, source_alt_index"
      } else {
        ""
      }
    )
  }
  invisible(dbExecute(
    con,
    glue(
      "CREATE OR REPLACE TEMP TABLE duckvep_throughput_variants AS
       SELECT
         (row_number() OVER (ORDER BY {input_order}) - 1)::UBIGINT
           AS input_variant_index,
         {event_projection}
       FROM {variants_table}
       ORDER BY {input_order}
       LIMIT {variant_sql}"
    )
  ))
  model_tables <- c("bench_regions", "bench_transcripts", "bench_exons")
  counts <- vapply(
    model_tables,
    function(table) dbGetQuery(con, glue("SELECT count(*) n FROM {table}"))$n[[1L]],
    numeric(1)
  )
  variant_count <- dbGetQuery(
    con,
    "SELECT count(*) n FROM duckvep_throughput_variants"
  )$n[[1L]]
  workload <- if (nzchar(opt$workload_name)) {
    opt$workload_name
  } else if (!is.null(corpus_receipt)) {
    paste0(
      "ensembl", corpus_receipt$model_source_version[[1L]], "_",
      tolower(corpus_receipt$assembly[[1L]]), "_",
      tolower(corpus_receipt$source_name[[1L]]),
      "_annotation_dense_v", corpus_receipt$schema_version[[1L]]
    )
  } else {
    paste0("production_", opt$variants_table)
  }
  region_columns <- dbGetQuery(
    con,
    "SELECT column_name FROM duckdb_columns() WHERE table_name = 'bench_regions'"
  )$column_name
  region_name_column <- if ("seq_region_name" %in% region_columns) {
    "seq_region_name"
  } else if ("name" %in% region_columns) {
    "name"
  } else if ("chrom" %in% region_columns) {
    "chrom"
  } else {
    ""
  }
  transcript_columns <- dbGetQuery(
    con,
    "SELECT column_name FROM duckdb_columns() WHERE table_name = 'bench_transcripts'"
  )$column_name
  complete_coverage <- "sequence_length" %in% region_columns
  if (
    identical(opt$output, "hgvs") &&
      (!complete_coverage || !nzchar(region_name_column))
  ) {
    die(
      "--output hgvs requires sequence_length and a seq_region_name, name, ",
      "or chrom column ",
      "in bench_regions"
    )
  }
  complete_flanks <- all(
    c("pre_cds_sequence", "post_cds_sequence") %in% transcript_columns
  )
  load_queries <- c(
    if (identical(opt$output, "hgvs")) {
      paste(
        "SELECT seq_region, sequence_length,",
        region_name_column,
        "AS seq_region_name FROM bench_regions ORDER BY seq_region"
      )
    } else {
      paste(
        "SELECT seq_region",
        if (complete_coverage) ", sequence_length" else "",
        "FROM bench_regions ORDER BY seq_region"
      )
    },
    paste(
      "SELECT transcript_index, seq_region, transcript_start, transcript_end,",
      "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table",
      if (complete_flanks) ", pre_cds_sequence, post_cds_sequence" else "",
      "FROM bench_transcripts ORDER BY transcript_index"
    ),
    paste(
      "SELECT transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end,",
      "phase, end_phase FROM bench_exons ORDER BY transcript_index, exon_cdna_start"
    )
  )
  present_relations <- dbGetQuery(
    con,
    paste(
      "SELECT table_name FROM information_schema.tables",
      "WHERE table_catalog = current_database() AND table_schema = 'main'"
    )
  )$table_name
  mature_mirna_query <- if ("bench_mature_mirna" %in% present_relations) {
    paste(
      "SELECT transcript_index, mature_mirna_start, mature_mirna_end",
      "FROM bench_mature_mirna ORDER BY transcript_index, mature_mirna_start"
    )
  } else {
    NULL
  }
  peptide_edit_query <- if ("bench_peptide_edits" %in% present_relations) {
    paste(
      "SELECT transcript_index, protein_position, alternate_amino_acid",
      "FROM bench_peptide_edits ORDER BY transcript_index, protein_position"
    )
  } else {
    NULL
  }
  regulation_relation <- if ("bench_regulation_features" %in% present_relations) {
    "bench_regulation_features"
  } else if ("duckvep_regulation_features" %in% present_relations) {
    "duckvep_regulation_features"
  } else {
    NULL
  }
  if (isTRUE(opt$regulatory) && is.null(regulation_relation)) {
    die(
      "--regulatory requires bench_regulation_features or ",
      "duckvep_regulation_features"
    )
  }
  interval_feature_query <- if (isTRUE(opt$regulatory)) {
    paste(
      "SELECT regulation_feature_index, seq_region, feature_start, feature_end,",
      "feature_kind FROM", regulation_relation,
      "ORDER BY seq_region, feature_start, regulation_feature_index"
    )
  } else {
    NULL
  }
  regulation_feature_count <- if (isTRUE(opt$regulatory)) {
    dbGetQuery(
      con,
      glue("SELECT count(*) AS n FROM {regulation_relation}")
    )$n[[1L]]
  } else {
    0
  }
  region_count <- counts[[1L]]
  transcript_count <- counts[[2L]]
  exon_count <- counts[[3L]]
} else {
  invisible(dbExecute(
    con,
    paste(readLines(opt$model_sql, warn = FALSE), collapse = "\n")
  ))
  variant_count <- opt$variants
  workload <- "fixture_one_transcript_sorted"
  region_count <- 1
  transcript_count <- 1
  exon_count <- 2
  complete_coverage <- identical(opt$output, "hgvs")
  complete_flanks <- FALSE
  load_queries <- c(
    if (identical(opt$output, "hgvs")) {
      paste(
        "SELECT seq_region,",
        "CASE seq_region WHEN 0 THEN 1 ELSE 260 END::UBIGINT AS sequence_length,",
        "name AS seq_region_name",
        "FROM duckvep_sequence_regions ORDER BY seq_region"
      )
    } else {
      "SELECT seq_region FROM duckvep_sequence_regions ORDER BY seq_region"
    },
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
  mature_mirna_query <- NULL
  peptide_edit_query <- NULL
  interval_feature_query <- NULL
  regulation_feature_count <- 0
}
if (variant_count < 1) {
  die("benchmark input contains no variants")
}
load_options <- character()
if (!is.null(mature_mirna_query)) {
  load_options <- c(
    load_options,
    glue("mature_mirna_query := {sql_q(mature_mirna_query)}")
  )
}
if (!is.null(peptide_edit_query)) {
  load_options <- c(
    load_options,
    glue("peptide_edit_query := {sql_q(peptide_edit_query)}")
  )
}
if (!is.null(interval_feature_query)) {
  load_options <- c(
    load_options,
    glue("interval_feature_query := {sql_q(interval_feature_query)}")
  )
}
if (nzchar(opt$reference_fasta)) {
  load_options <- c(
    load_options,
    glue("reference_fasta := {sql_q(normalizePath(opt$reference_fasta))}")
  )
}
if (complete_coverage) {
  load_options <- c(load_options, "transcript_coverage_complete := TRUE")
}
model_load_sql <- glue(
  "SELECT loaded FROM duckvep_model_load(
     {sql_q(model_name)}, {sql_q(load_queries[1])},
     {sql_q(load_queries[2])}, {sql_q(load_queries[3])}
     {if (length(load_options)) paste0(', ', paste(load_options, collapse = ', ')) else ''})"
)
model_load_timing <- system.time({
  loaded <- dbGetQuery(con, model_load_sql)$loaded
})
if (!identical(loaded, TRUE)) {
  die("model load failed")
}

annotation_input <- function(begin, count) {
  begin_sql <- whole_nonnegative(begin, "input partition start")
  end_sql <- whole_nonnegative(begin + count, "input partition end")
  if (production) {
    columns <- if (identical(opt$event_mode, "breakend")) {
      "seq_region, \"position\", mate_seq_region, mate_position"
    } else {
      "seq_region, \"position\", \"reference\", \"alternate\""
    }
    glue(
      "SELECT input_variant_index, {columns}
       FROM duckvep_throughput_variants
       WHERE input_variant_index >= {begin_sql}
         AND input_variant_index < {end_sql}
       ORDER BY input_variant_index"
    )
  } else if (identical(opt$event_mode, "breakend")) {
    glue(
      "SELECT i::UBIGINT AS input_variant_index,
              1::UINTEGER AS seq_region, 159::UBIGINT AS \"position\",
              1::UINTEGER AS mate_seq_region, 170::UBIGINT AS mate_position
       FROM range({begin_sql}, {end_sql}) r(i)"
    )
  } else {
    glue(
      "SELECT i::UBIGINT AS input_variant_index,
              1::UINTEGER AS seq_region, 124::UBIGINT AS \"position\",
              'T' AS \"reference\",
              CASE WHEN i % 2 = 0 THEN 'C' ELSE 'G' END AS \"alternate\"
       FROM range({begin_sql}, {end_sql}) r(i)"
    )
  }
}

annotation_function_name <- function(output_mode) {
  if (identical(opt$event_mode, "breakend")) {
    if (output_mode == "compact") {
      "duckvep_annotate_breakend_compact"
    } else {
      "duckvep_annotate_breakend"
    }
  } else if (output_mode == "compact") {
    "duckvep_annotate_compact"
  } else if (output_mode == "hgvs") {
    "duckvep_annotate_hgvs"
  } else {
    "duckvep_annotate"
  }
}

annotation_cte <- function(n, output_mode) {
  n <- as.numeric(n)
  function_name <- annotation_function_name(output_mode)
  function_arguments <- if (identical(opt$event_mode, "breakend")) {
    "seq_region, \"position\", mate_seq_region, mate_position"
  } else {
    "seq_region, \"position\", \"reference\", \"alternate\""
  }
  partition_count <- min(opt$input_partitions, n)
  partition_begin <- floor((0:(partition_count - 1)) * n / partition_count)
  partition_end <- floor((1:partition_count) * n / partition_count)
  branches <- vapply(
    seq_len(partition_count),
    function(index) {
      input <- annotation_input(
        partition_begin[[index]],
        partition_end[[index]] - partition_begin[[index]]
      )
      glue(
        "SELECT
           input_variant_index,
           seq_region,
           \"position\",
           {if (identical(opt$event_mode, 'breakend')) 'mate_seq_region, mate_position,' else '\"reference\", \"alternate\",'}
           {function_name}(
             {sql_q(model_name)}, {function_arguments},
             {distance_sql}::UBIGINT
           ) AS annotations
         FROM ({input}) AS variants"
      )
    },
    character(1L)
  )
  glue(
    "WITH annotation_lists AS (
       {paste(branches, collapse = '\nUNION ALL\n')}
     ), annotated AS (
       SELECT
         * EXCLUDE (annotations),
         generate_subscripts(annotations, 1)::UINTEGER AS annotation_index,
         unnest(annotations) AS annotation
       FROM annotation_lists
     )"
  )
}

annotation_query <- function(n) {
  checksum <- if (opt$output == "compact") {
    "CAST(sum(annotation.consequence_mask) AS VARCHAR)"
  } else if (opt$output == "hgvs") {
    paste(
      "CAST(sum(",
      "length(coalesce(annotation.transcript_hgvs, '')) +",
      "length(coalesce(annotation.protein_hgvs, '')) +",
      "length(annotation.transcript_hgvs_status) +",
      "length(annotation.protein_hgvs_status)",
      ") AS VARCHAR)"
    )
  } else {
    "CAST(sum(length(annotation.consequence)) AS VARCHAR)"
  }
  glue(
    "{annotation_cte(n, opt$output)}
     SELECT count(*)::DOUBLE AS annotated_rows, {checksum} AS checksum
     FROM annotated"
  )
}

composition_query <- function(n) {
  glue(
    "{annotation_cte(n, 'compact')}
     SELECT
       annotation.consequence_mask,
       annotation.region_mask,
       annotation.impact_code,
       annotation.status_code,
       annotation.reason_code,
       annotation.overlap_object_code,
       count(*)::DOUBLE AS output_rows
     FROM annotated
     GROUP BY ALL
     ORDER BY ALL"
  )
}

fingerprint_query <- function(n) {
  identity_fields <- if (identical(opt$event_mode, "breakend")) {
    c(
      "input_variant_index", "seq_region", "\"position\"",
      "mate_seq_region", "mate_position", "annotation_index"
    )
  } else {
    c(
      "input_variant_index", "seq_region", "\"position\"",
      "\"reference\"", "\"alternate\"", "annotation_index"
    )
  }
  annotation_fields <- if (opt$output == "rich") {
    c(
      "transcript_index", "gene_index", "consequence", "impact", "region",
      "status", "reason", "cdna_position", "cds_position", "protein_position",
      "reference_amino_acid", "alternate_amino_acid", "nmd_prediction",
      "nmd_escape_intronless", "nmd_escape_early_cds", "nmd_escape_last_exon",
      "nmd_escape_penultimate_exon_end", "regulation_feature_index",
      "overlap_object"
    )
  } else {
    compact_fields <- c(
      "transcript_index", "gene_index", "consequence_mask", "region_mask",
      "impact_code", "status_code", "reason_code", "cdna_position",
      "cds_position", "protein_position", "reference_amino_acid_code",
      "alternate_amino_acid_code", "nmd_prediction_code",
      "nmd_escape_reasons", "regulation_feature_index", "overlap_object_code"
    )
    if (opt$output == "hgvs") {
      c(
        compact_fields, "transcript_hgvs", "protein_hgvs", "hgvs_shift",
        "transcript_hgvs_status", "transcript_hgvs_reason",
        "protein_hgvs_status", "protein_hgvs_reason"
      )
    } else {
      compact_fields
    }
  }
  hash_arguments <- paste(
    c(identity_fields, paste0("annotation.", annotation_fields)),
    collapse = ", "
  )
  glue(
    "{annotation_cte(n, opt$output)}, row_fingerprints AS (
       SELECT hash({hash_arguments}) AS row_hash
       FROM annotated
     )
     SELECT
       count(*)::UBIGINT AS output_rows,
       bit_xor(row_hash)::VARCHAR AS xor_hash,
       sum(row_hash)::VARCHAR AS sum_hash
     FROM row_fingerprints"
  )
}

warmup_count <- min(opt$warmup, variant_count)
warmup_sql <- whole_count(warmup_count, "--warmup")
variant_sql <- whole_count(variant_count, "--variants")
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
if (length(unique(check$annotated_rows)) != 1L || check$annotated_rows[[1L]] < 1) {
  die("annotation cardinality changed across benchmark passes")
}
if (length(unique(check$checksum)) != 1L) {
  die("annotation checksum changed across benchmark passes")
}

if (nzchar(opt$composition)) {
  masks <- dbGetQuery(con, composition_query(variant_sql))
  composition_rows <- sum(masks$output_rows)
  if (!identical(as.numeric(composition_rows), as.numeric(check$annotated_rows[[1L]]))) {
    die("compact composition cardinality differs from timed output")
  }

  aggregate_codes <- function(codes, values, labels) {
    counts <- aggregate(
      masks$output_rows,
      by = list(code = codes),
      FUN = sum
    )
    names(counts)[[2L]] <- "output_rows"
    counts$category <- unname(labels[as.character(counts$code)])
    if (anyNA(counts$category)) {
      die("composition contains an unknown numeric code")
    }
    counts[, c("category", "output_rows"), drop = FALSE]
  }
  append_dimension <- function(parts, dimension, values) {
    if (!nrow(values)) return(parts)
    values$dimension <- dimension
    rbind(parts, values[, c("dimension", "category", "output_rows")])
  }

  composition <- data.frame(
    dimension = character(),
    category = character(),
    output_rows = numeric(),
    stringsAsFactors = FALSE
  )
  object_labels <- c(
    `0` = "transcript",
    `1` = "regulatory_region",
    `2` = "transcription_factor_binding_site"
  )
  impact_labels <- c(
    `0` = "modifier", `1` = "low", `2` = "moderate", `3` = "high"
  )
  status_labels <- c(`0` = "supported", `1` = "unresolved")
  reason_labels <- c(
    `0` = "none",
    `1` = "no_feature_in_loaded_model",
    `2` = "missing_sequence",
    `3` = "ambiguous_sequence",
    `4` = "reference_mismatch",
    `5` = "non_contiguous_cds_edit",
    `6` = "unsupported_compound_consequence",
    `7` = "invalid_model_projection",
    `8` = "internal_capacity_error",
    `9` = "missing_transcript_tail",
    `10` = "missing_transcript_flank"
  )
  composition <- append_dimension(
    composition,
    "overlap_object",
    aggregate_codes(masks$overlap_object_code, masks$output_rows, object_labels)
  )
  composition <- append_dimension(
    composition,
    "impact",
    aggregate_codes(masks$impact_code, masks$output_rows, impact_labels)
  )
  composition <- append_dimension(
    composition,
    "status",
    aggregate_codes(masks$status_code, masks$output_rows, status_labels)
  )
  composition <- append_dimension(
    composition,
    "reason",
    aggregate_codes(masks$reason_code, masks$output_rows, reason_labels)
  )

  region_labels <- c(
    "upstream", "downstream", "intron", "exon", "CDS", "UTR", "splice"
  )
  region_masks <- as.numeric(masks$region_mask)
  for (bit in seq_along(region_labels) - 1L) {
    present <- floor(region_masks / (2^bit)) %% 2 == 1
    if (any(present)) {
      composition <- append_dimension(
        composition,
        "region",
        data.frame(
          category = region_labels[[bit + 1L]],
          output_rows = sum(masks$output_rows[present]),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  if (any(region_masks == 0)) {
    composition <- append_dimension(
      composition,
      "region",
      data.frame(
        category = "none",
        output_rows = sum(masks$output_rows[region_masks == 0]),
        stringsAsFactors = FALSE
      )
    )
  }

  so_bindings <- utils::read.delim(
    file.path(
      root,
      "test",
      "duckvep",
      "conformance",
      "data",
      "so_bit_bindings.tsv"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  consequence_masks <- as.numeric(masks$consequence_mask)
  for (bit in seq_len(nrow(so_bindings)) - 1L) {
    present <- floor(consequence_masks / (2^bit)) %% 2 == 1
    if (any(present)) {
      composition <- append_dimension(
        composition,
        "SO_term",
        data.frame(
          category = so_bindings$SO_term[[bit + 1L]],
          output_rows = sum(masks$output_rows[present]),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  composition$workload <- workload
  composition$timed_output_mode <- opt$output
  composition$transcript_distance <- opt$transcript_distance
  composition$input_variants <- variant_count
  composition$denominator_output_rows <- composition_rows
  composition$share_of_output_rows <- composition$output_rows / composition_rows
  composition <- composition[
    order(composition$dimension, -composition$output_rows, composition$category),
    c(
      "workload", "timed_output_mode", "transcript_distance", "input_variants",
      "denominator_output_rows", "dimension", "category", "output_rows",
      "share_of_output_rows"
    ),
    drop = FALSE
  ]
  dir.create(dirname(opt$composition), recursive = TRUE, showWarnings = FALSE)
  tmp_composition <- tempfile(
    "duckvep-composition-",
    dirname(opt$composition),
    ".csv"
  )
  utils::write.csv(composition, tmp_composition, row.names = FALSE)
  if (!file.rename(tmp_composition, opt$composition)) {
    unlink(tmp_composition)
    die("cannot replace composition: {opt$composition}")
  }
}

if (nzchar(opt$fingerprint) || nzchar(opt$expected_fingerprint)) {
  fingerprint_values <- dbGetQuery(con, fingerprint_query(variant_sql))
  if (nrow(fingerprint_values) != 1L ||
      fingerprint_values$output_rows[[1L]] != check$annotated_rows[[1L]]) {
    die("full-row fingerprint cardinality differs from timed output")
  }
  fingerprint <- data.frame(
    schema_version = 2L,
    workload = workload,
    output_mode = opt$output,
    transcript_distance = as.character(opt$transcript_distance),
    input_variants = as.integer(variant_count),
    output_rows = as.character(fingerprint_values$output_rows[[1L]]),
    xor_hash = as.character(fingerprint_values$xor_hash[[1L]]),
    sum_hash = as.character(fingerprint_values$sum_hash[[1L]]),
    model_sha256 = model_logical_sha256,
    region_ordinal_sha256 = model_region_ordinal_sha256,
    variants_database_sha256 = variants_database_sha256,
    stringsAsFactors = FALSE
  )
  if (nzchar(opt$expected_fingerprint)) {
    expected <- utils::read.csv(
      opt$expected_fingerprint,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      colClasses = "character"
    )
    observed <- fingerprint
    observed[] <- lapply(observed, as.character)
    if (!identical(names(expected), names(observed)) ||
        nrow(expected) != 1L || !identical(expected, observed)) {
      die("full-row fingerprint differs from baseline: {opt$expected_fingerprint}")
    }
  }
  if (nzchar(opt$fingerprint)) {
    dir.create(dirname(opt$fingerprint), recursive = TRUE, showWarnings = FALSE)
    tmp_fingerprint <- tempfile(
      "duckvep-fingerprint-",
      dirname(opt$fingerprint),
      ".csv"
    )
    utils::write.csv(fingerprint, tmp_fingerprint, row.names = FALSE)
    if (!file.rename(tmp_fingerprint, opt$fingerprint)) {
      unlink(tmp_fingerprint)
      die("cannot replace fingerprint: {opt$fingerprint}")
    }
  }
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
cpu_affinity <- if (file.exists("/proc/self/status")) {
  line <- grep(
    "^Cpus_allowed_list[[:space:]]*:",
    readLines("/proc/self/status"),
    value = TRUE
  )
  if (length(line)) trimws(sub("^[^:]+:", "", line[[1L]])) else "unknown"
} else {
  "unknown"
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
  workload = workload,
  output_mode = opt$output,
  input_order = if (opt$input_partitions == 1L) {
    "nondecreasing_seq_region_position"
  } else {
    "disjoint_nondecreasing_seq_region_position"
  },
  threads = opt$threads,
  input_partitions = opt$input_partitions,
  transcript_distance = as.numeric(opt$transcript_distance),
  variants = as.integer(variant_count),
  regions = as.integer(region_count),
  transcripts = as.integer(transcript_count),
  exons = as.integer(exon_count),
  regulation_features = as.integer(regulation_feature_count),
  passes = opt$passes,
  warmup_variants = as.integer(warmup_count),
  min_seconds = min(elapsed),
  median_seconds = median_seconds,
  max_seconds = max(elapsed),
  variants_per_second = variant_count / median_seconds,
  ns_per_variant = median_seconds * 1e9 / variant_count,
  annotated_rows = as.integer(check$annotated_rows[[1L]]),
  checksum_kind = if (opt$output == "compact") {
    "consequence_mask_sum"
  } else if (opt$output == "hgvs") {
    "hgvs_text_status_bytes"
  } else {
    "consequence_text_bytes"
  },
  checksum_value = check$checksum[[1L]],
  cpu_affinity = cpu_affinity,
  stringsAsFactors = FALSE
)

if (nzchar(opt$history)) {
  dir.create(dirname(opt$history), recursive = TRUE, showWarnings = FALSE)
  rows <- row
  if (file.exists(opt$history)) {
    old <- utils::read.csv(
      opt$history,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      colClasses = c(
        source_revision = "character",
        checksum_value = "character"
      )
    )
    if (!identical(names(old), names(row))) {
      die("history schema does not match: {opt$history}")
    }
    same_affinity <- !is.na(old$cpu_affinity) &
      old$cpu_affinity == row$cpu_affinity
    same <- old$source_revision == row$source_revision &
      old$host == row$host &
      old$workload == row$workload &
      old$output_mode == row$output_mode &
      old$threads == row$threads &
      old$input_partitions == row$input_partitions &
      old$transcript_distance == row$transcript_distance &
      old$variants == row$variants &
      same_affinity
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
    "{format(variant_count, big.mark = ',', scientific = FALSE)} sorted variants; median ",
    "{sprintf('%.3f', median_seconds)} s; ",
    "{sprintf('%.1f', row$ns_per_variant)} ns/variant; ",
    "{format(check$annotated_rows[[1L]], big.mark = ',', scientific = FALSE)} rows; ",
    "{sprintf('%.2f', check$annotated_rows[[1L]] / variant_count)} rows/variant"
  ),
  "\n",
  sep = ""
)
cat(
  glue(
    "model load: {sprintf('%.3f', unname(model_load_timing[['elapsed']]))} s; ",
    "input partitions: {format(opt$input_partitions, big.mark = ',')}; ",
    "transcript distance: {format(opt$transcript_distance, big.mark = ',')} bases; ",
    "complete coverage: {tolower(as.character(complete_coverage))}; ",
    "complete transcript flanks: {tolower(as.character(complete_flanks))}; ",
    "regulation features: {format(regulation_feature_count, big.mark = ',')}"
  ),
  "\n",
  sep = ""
)
if (nzchar(opt$history)) {
  cat(glue("history -> {opt$history}"), "\n", sep = "")
}
if (nzchar(opt$composition)) {
  cat(glue("composition -> {opt$composition}"), "\n", sep = "")
}
