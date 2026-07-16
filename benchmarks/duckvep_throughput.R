#!/usr/bin/env Rscript

# Record the stable DuckDB C-API path over a sorted variant stream. With no
# database this uses the checked-in one-transcript fixture. --database expects
# pre-staged bench_regions, bench_transcripts, bench_exons, and bench_variants.
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
if (opt$passes < 1L) {
  die("--passes must be positive")
}
if (opt$threads < 1L) {
  die("--threads must be positive")
}
if (!opt$output %in% c("rich", "compact")) {
  die("--output must be rich or compact")
}
if (!opt$event_mode %in% c("small", "breakend")) {
  die("--event-mode must be small or breakend")
}
production <- nzchar(opt$database)
if (isTRUE(opt$regulatory) && !production) {
  die("--regulatory requires a production --database")
}
inputs <- c(opt$extension, if (production) opt$database else opt$model_sql)
missing <- inputs[!file.exists(inputs)]
if (length(missing) != 0L) {
  die("missing input(s):\n{paste(missing, collapse = '\n')}")
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

model_name <- "throughput"
if (production) {
  variants_table <- identifier_q(opt$variants_table)
  event_projection <- if (identical(opt$event_mode, "breakend")) {
    "seq_region, \"position\", mate_seq_region, mate_position"
  } else {
    "seq_region, \"position\", \"reference\", \"alternate\""
  }
  invisible(dbExecute(
    con,
    glue(
      "CREATE OR REPLACE TEMP TABLE duckvep_throughput_variants AS
       SELECT {event_projection}
       FROM {variants_table}
       ORDER BY seq_region, \"position\"
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
  } else {
    "ensembl116_grch38_giab_sites_hash40"
  }
  region_columns <- dbGetQuery(
    con,
    "SELECT column_name FROM duckdb_columns() WHERE table_name = 'bench_regions'"
  )$column_name
  transcript_columns <- dbGetQuery(
    con,
    "SELECT column_name FROM duckdb_columns() WHERE table_name = 'bench_transcripts'"
  )$column_name
  complete_coverage <- "sequence_length" %in% region_columns
  complete_flanks <- all(
    c("pre_cds_sequence", "post_cds_sequence") %in% transcript_columns
  )
  load_queries <- c(
    paste(
      "SELECT seq_region",
      if (complete_coverage) ", sequence_length" else "",
      "FROM bench_regions ORDER BY seq_region"
    ),
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
  complete_coverage <- FALSE
  complete_flanks <- TRUE
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

annotation_query <- function(n) {
  input <- if (production) {
    columns <- if (identical(opt$event_mode, "breakend")) {
      "seq_region, \"position\", mate_seq_region, mate_position"
    } else {
      "seq_region, \"position\", \"reference\", \"alternate\""
    }
    glue(
      "SELECT {columns} FROM duckvep_throughput_variants
       ORDER BY seq_region, \"position\" LIMIT {n}"
    )
  } else if (identical(opt$event_mode, "breakend")) {
    glue(
      "SELECT 1::UINTEGER AS seq_region, 159::UBIGINT AS \"position\",
              1::UINTEGER AS mate_seq_region, 170::UBIGINT AS mate_position
       FROM range({n})"
    )
  } else {
    glue(
      "SELECT 1::UINTEGER AS seq_region, 124::UBIGINT AS \"position\",
              'T' AS \"reference\",
              CASE WHEN i % 2 = 0 THEN 'C' ELSE 'G' END AS \"alternate\"
       FROM range({n}) r(i)"
    )
  }
  function_name <- if (identical(opt$event_mode, "breakend")) {
    if (opt$output == "compact") {
      "duckvep_annotate_breakend_compact"
    } else {
      "duckvep_annotate_breakend"
    }
  } else if (opt$output == "compact") {
    "duckvep_annotate_compact"
  } else {
    "duckvep_annotate"
  }
  checksum <- if (opt$output == "compact") {
    "CAST(sum(annotation.consequence_mask) AS VARCHAR)"
  } else {
    "CAST(sum(length(annotation.consequence)) AS VARCHAR)"
  }
  function_arguments <- if (identical(opt$event_mode, "breakend")) {
    "seq_region, \"position\", mate_seq_region, mate_position"
  } else {
    "seq_region, \"position\", \"reference\", \"alternate\""
  }
  glue(
    "WITH variants AS ({input}), annotated AS (
       SELECT unnest({function_name}(
         {sql_q(model_name)}, {function_arguments},
         5000::UBIGINT
       )) AS annotation
       FROM variants
     )
     SELECT count(*)::DOUBLE AS annotated_rows, {checksum} AS checksum
     FROM annotated"
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
  workload = workload,
  output_mode = opt$output,
  input_order = "nondecreasing_seq_region_position",
  threads = opt$threads,
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
  } else {
    "consequence_text_bytes"
  },
  checksum_value = check$checksum[[1L]],
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
    same <- old$source_revision == row$source_revision &
      old$host == row$host &
      old$workload == row$workload &
      old$output_mode == row$output_mode &
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
    "{format(variant_count, big.mark = ',', scientific = FALSE)} sorted variants; median ",
    "{sprintf('%.3f', median_seconds)} s; ",
    "{sprintf('%.1f', row$ns_per_variant)} ns/variant"
  ),
  "\n",
  sep = ""
)
cat(
  glue(
    "model load: {sprintf('%.3f', unname(model_load_timing[['elapsed']]))} s; ",
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
