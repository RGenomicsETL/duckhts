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
  help = "event family: small, structural, or breakend [default: %default]"
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
op <- add_option(
  op,
  "--max-allele-length",
  dest = "max_allele_length",
  type = "integer",
  default = 50L,
  help = paste(
    "maximum literal REF or ALT length in the small-event lane;",
    "the compact kernel limit is 65535 [%default]"
  )
)
op <- add_option(
  op,
  "--split-multiallelic",
  dest = "split_multiallelic",
  action = "store_true",
  default = FALSE,
  help = paste(
    "compare each ALT from a multiallelic record independently; genotype",
    "fields are not used or rewritten [default: %default]"
  )
)
op <- add_option(
  op,
  "--stratify-raw-allele-length",
  dest = "stratify_raw_allele_length",
  action = "store_true",
  default = FALSE,
  help = paste(
    "also partition deterministic small-event sampling by the larger raw",
    "REF/ALT length (1, 2..50, 51..1000, 1001..10000, >10000);",
    "this preserves long uploaded representations independently of the",
    "minimized edit-length bin [default: %default]"
  )
)
op <- add_option(op, "--seed", type = "integer", default = 1L)
op <- add_option(op, "--chrom", default = "")
op <- add_option(op, "--distance", type = "integer", default = 5000L)
op <- add_option(
  op,
  "--max-sv-size",
  dest = "max_sv_size",
  type = "double",
  default = 10000000,
  help = paste(
    "VEP --max_sv_size for structural-mode records; the runner passes 10 Mb",
    "by default so exact spans are not silently skipped, and records the value",
    "in the oracle receipt. This is not VEP's own 5 kb default [%default]"
  )
)
op <- add_option(
  op,
  "--regulatory",
  action = "store_true",
  default = FALSE,
  help = paste(
    "load duckvep_regulation_features and compare VEP --regulatory output",
    "[default: %default]"
  )
)
op <- add_option(
  op,
  "--hgvs",
  action = "store_true",
  default = FALSE,
  help = paste(
    "also compare independent-event HGVSc/HGVSn and HGVSp suffixes from",
    "duckvep_annotate_hgvs(...) with VEP --hgvs [default: %default]"
  )
)
op <- add_option(op, "--hgvs-out", dest = "hgvs_out", default = "")
op <- add_option(
  op,
  "--hgvs-pairs-out",
  dest = "hgvs_pairs_out",
  default = ""
)
op <- add_option(
  op,
  "--allow-hgvs-discordance",
  dest = "allow_hgvs_discordance",
  action = "store_true",
  default = FALSE,
  help = paste(
    "write investigative HGVS artifacts without failing on unresolved,",
    "mismatched, missing, or extra rows [default: fail closed]"
  )
)
op <- add_option(
  op,
  "--annotations-out",
  dest = "annotations_out",
  default = ""
)
op <- add_option(
  op,
  "--eligibility-out",
  dest = "eligibility_out",
  default = "",
  help = "optional CSV receipt for input allele classes and exclusions [%default]"
)
op <- add_option(
  op,
  "--source-audit-only",
  dest = "source_audit_only",
  action = "store_true",
  default = FALSE,
  help = paste(
    "scan the complete small-event source, write --eligibility-out, and stop",
    "before FASTA/VEP annotation; this records source support without",
    "claiming executable conformance [default: %default]"
  )
)
op <- add_option(
  op,
  "--source-url",
  dest = "source_url",
  default = "",
  help = "published URL from which --vcf was obtained [%default]"
)
op <- add_option(
  op,
  "--source-version",
  dest = "source_version",
  default = "",
  help = "release, object version, or accession identifying --vcf [%default]"
)
op <- add_option(
  op,
  "--source-checksum",
  dest = "source_checksum",
  default = "",
  help = "expected checksum as md5:HEX or sha256:HEX [%default]"
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
  "--vep-buffer-size",
  dest = "vep_buffer_size",
  type = "integer",
  default = 5000L,
  help = paste(
    "VEP input buffer size for small and structural modes; breakend mode",
    "remains isolated at one because neighboring BNDs change VEP semantics",
    "[%default]"
  )
)
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
if (!(opt$event_mode %in% c("small", "structural", "breakend"))) {
  die("--event-mode must be small, structural, or breakend")
}
if (isTRUE(opt$hgvs) && !identical(opt$event_mode, "small")) {
  die("--hgvs currently applies only to --event-mode small")
}
if (isTRUE(opt$hgvs) && isTRUE(opt$regulatory)) {
  die("--hgvs compares transcript rows and cannot be combined with --regulatory")
}
if (
  !isTRUE(opt$source_audit_only) &&
    (
      !requireNamespace("blit", quietly = TRUE) ||
        utils::packageVersion("blit") < "0.2.0.9000"
    )
) {
  die(
    "the current WangLabCSU/blit checkout is required ",
    "(blit >= 0.2.0.9000)"
  )
}
oracle_mode <- if (nzchar(opt$cache_dir)) "cache" else "gff"
fasta_index <- paste0(opt$fasta, ".fai")
required_files <- opt$extension
if (!isTRUE(opt$source_audit_only)) {
  required_files <- c(opt$fasta, fasta_index, required_files)
}
external_model_database <- !nzchar(opt$model_sql) &&
  !identical(opt$database, ":memory:")
if (external_model_database) {
  required_files <- c(required_files, opt$database)
}
if (!isTRUE(opt$source_audit_only) && identical(oracle_mode, "gff")) {
  required_files <- c(opt$gff, required_files)
} else if (
  !isTRUE(opt$source_audit_only) &&
    !dir.exists(opt$cache_dir)
) {
  die("VEP cache root does not exist: {opt$cache_dir}")
}
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) != 0L) {
  die("missing input(s):\n{paste(missing_files, collapse = '\n')}")
}

fasta_regions <- if (isTRUE(opt$source_audit_only)) {
  data.frame(chrom = character(), sequence_length = numeric())
} else {
  regions <- utils::read.delim(
    fasta_index,
    header = FALSE,
    colClasses = c("character", "numeric", rep("NULL", 3L)),
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  names(regions) <- c("chrom", "sequence_length")
  regions
}
fasta_chroms <- fasta_regions$chrom
if (
  !isTRUE(opt$source_audit_only) &&
    (
      length(fasta_chroms) == 0L ||
        any(!nzchar(fasta_chroms)) ||
        anyDuplicated(fasta_chroms) ||
        any(!is.finite(fasta_regions$sequence_length)) ||
        any(fasta_regions$sequence_length < 1) ||
        any(fasta_regions$sequence_length > 4294967295) ||
        any(fasta_regions$sequence_length != floor(fasta_regions$sequence_length))
    )
) {
  die("FASTA index has invalid names or lengths: {fasta_index}")
}
if (opt$sample_per_shape < 0L) {
  die("--sample-per-shape must be non-negative")
}
if (opt$max_allele_length < 1L || opt$max_allele_length > 65535L) {
  die("--max-allele-length must be between 1 and 65535")
}
if (opt$distance < 0L) {
  die("--distance must be non-negative")
}
if (opt$vep_buffer_size < 1L) {
  die("--vep-buffer-size must be positive")
}
if (
  !is.finite(opt$max_sv_size) ||
    opt$max_sv_size < 1 ||
    opt$max_sv_size > 4294967295 ||
    opt$max_sv_size != floor(opt$max_sv_size)
) {
  die("--max-sv-size must be an integer between 1 and 4294967295")
}
if (
  isTRUE(opt$source_audit_only) &&
    (
      !identical(opt$event_mode, "small") ||
        !nzchar(opt$eligibility_out)
    )
) {
  die("--source-audit-only requires small mode and --eligibility-out")
}
if (
  !identical(opt$event_mode, "small") &&
    (
      isTRUE(opt$split_multiallelic) ||
        isTRUE(opt$stratify_raw_allele_length) ||
        nzchar(opt$eligibility_out)
    )
) {
  die(
    "--split-multiallelic, --stratify-raw-allele-length, and ",
    "--eligibility-out apply only to --event-mode small"
  )
}
if (nzchar(opt$source_version) && !nzchar(opt$source_url)) {
  die("--source-version requires --source-url")
}
if (isTRUE(opt$regulatory) && identical(oracle_mode, "gff")) {
  die("--regulatory requires the indexed VEP cache, not the GFF oracle")
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
  if (isTRUE(opt$source_audit_only)) {
    die("--source-audit-only does not execute the NMD plugin")
  }
  if (!identical(opt$event_mode, "small")) {
    die("structural and breakend differentials do not compare NMD plugin output")
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
if (isTRUE(opt$hgvs) && !nzchar(opt$hgvs_out)) {
  opt$hgvs_out <- file.path(
    results_dir,
    glue("{opt$corpus}_hgvs_conformance.csv")
  )
}
if (isTRUE(opt$hgvs) && !nzchar(opt$hgvs_pairs_out)) {
  opt$hgvs_pairs_out <- file.path(
    results_dir,
    glue("{opt$corpus}_hgvs_pairs.parquet")
  )
}

canonical_output_path <- function(path) {
  normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
}
annotation_stem <- sub(
  "_annotations[.]parquet$",
  "",
  basename(opt$annotations_out)
)
statistical_outputs <- c(
  statistical_conformance = file.path(
    dirname(opt$annotations_out),
    glue("{annotation_stem}_statistical_conformance.csv")
  ),
  methodology_audit = file.path(
    dirname(opt$annotations_out),
    glue("{annotation_stem}_methodology_audit.csv")
  ),
  consequence_pairs = file.path(
    dirname(opt$annotations_out),
    glue("{annotation_stem}_pairs.parquet")
  ),
  nmd_conformance = file.path(
    dirname(opt$annotations_out),
    glue("{annotation_stem}_nmd_conformance.csv")
  )
)
declared_outputs <- c(annotations = opt$annotations_out, statistical_outputs)
hgvs_discordance_count <- 0L
if (isTRUE(opt$hgvs)) {
  declared_outputs <- c(
    declared_outputs,
    hgvs_summary = opt$hgvs_out,
    hgvs_pairs = opt$hgvs_pairs_out
  )
}
canonical_outputs <- vapply(
  declared_outputs,
  canonical_output_path,
  character(1L)
)
colliding_outputs <- unique(canonical_outputs[duplicated(canonical_outputs)])
if (length(colliding_outputs) != 0L) {
  labels <- vapply(
    colliding_outputs,
    function(path) paste(names(canonical_outputs)[canonical_outputs == path],
      collapse = ", "
    ),
    character(1L)
  )
  die(
    "output paths must be distinct; collision between ",
    paste(labels, collapse = "; ")
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
invisible(dbWriteTable(
  con,
  "duckvep_oracle_sequence_order",
  data.frame(
    chrom = fasta_chroms,
    sequence_length = fasta_regions$sequence_length,
    fasta_order = seq_along(fasta_chroms),
    stringsAsFactors = FALSE
  ),
  temporary = TRUE
))
if (isTRUE(opt$hgvs)) {
  so_bit_path <- file.path(
    root,
    "test",
    "duckvep",
    "conformance",
    "data",
    "so_bit_bindings.tsv"
  )
  if (!file.exists(so_bit_path)) {
    die("missing generated SO-bit binding table: {so_bit_path}")
  }
  so_bits <- utils::read.delim(
    so_bit_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (
    !identical(names(so_bits), c("SO_term", "so_enum")) ||
      nrow(so_bits) == 0L || nrow(so_bits) > 64L ||
      anyNA(so_bits$SO_term) || any(!nzchar(so_bits$SO_term)) ||
      anyDuplicated(so_bits$SO_term)
  ) {
    die("invalid generated SO-bit binding table: {so_bit_path}")
  }
  so_bits$bit_index <- seq_len(nrow(so_bits)) - 1L
  invisible(dbWriteTable(
    con,
    "duckvep_so_bits",
    so_bits[, c("SO_term", "bit_index")],
    temporary = TRUE
  ))
}

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
if (isTRUE(opt$regulatory)) {
  needed_relations <- c(needed_relations, "duckvep_regulation_features")
}
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

if (isTRUE(opt$regulatory)) {
  regulation_columns <- dbGetQuery(
    con,
    "SELECT column_name FROM duckdb_columns()
     WHERE table_name = 'duckvep_regulation_features'"
  )$column_name
  required_regulation_columns <- c(
    "regulation_feature_index",
    "seq_region",
    "feature_start",
    "feature_end",
    "feature_kind",
    "stable_id"
  )
  missing_regulation_columns <- setdiff(
    required_regulation_columns,
    regulation_columns
  )
  if (length(missing_regulation_columns) != 0L) {
    die(
      "duckvep_regulation_features is missing column(s): ",
      "{paste(missing_regulation_columns, collapse = ', ')}"
    )
  }
}

region_columns <- dbGetQuery(
  con,
  "SELECT column_name FROM duckdb_columns()
   WHERE table_name = 'duckvep_sequence_regions'"
)$column_name
complete_coverage <- "sequence_length" %in% region_columns
model_sequence_length_sql <- if (complete_coverage) {
  "r.sequence_length"
} else {
  "4294967294::UBIGINT"
}
model_query_relation <- function(relation) {
  relation_id <- as.character(dbQuoteIdentifier(con, relation))
  if (external_model_database) {
    return(glue("duckvep_model_source.main.{relation_id}"))
  }
  relation_id
}
hgvs_region_query <- ""
if (isTRUE(opt$hgvs)) {
  if (!("name" %in% region_columns)) {
    die("--hgvs requires name in duckvep_sequence_regions")
  }
  region_length_projection <- if (complete_coverage) {
    ", sequence_length"
  } else {
    ""
  }
  model_regions <- dbGetQuery(
    con,
    paste(
      "SELECT seq_region, name",
      region_length_projection,
      "FROM", model_query_relation("duckvep_sequence_regions"),
      "ORDER BY seq_region"
    )
  )
  if (
    nrow(model_regions) == 0L || anyNA(model_regions[, c("seq_region", "name")]) ||
      any(!nzchar(model_regions$name)) || anyDuplicated(model_regions$seq_region) ||
      anyDuplicated(model_regions$name)
  ) {
    die("--hgvs requires unique non-NULL sequence-region ordinals and names")
  }
  fasta_match <- match(model_regions$name, fasta_regions$chrom)
  if (anyNA(fasta_match)) {
    die(
      "model sequence region(s) are absent from the FASTA index: ",
      "{paste(model_regions$name[is.na(fasta_match)], collapse = ', ')}"
    )
  }
  fasta_lengths <- fasta_regions$sequence_length[fasta_match]
  if (
    complete_coverage &&
      any(as.numeric(model_regions$sequence_length) != fasta_lengths)
  ) {
    die("model sequence lengths disagree with the FASTA index")
  }
  value_rows <- vapply(
    seq_len(nrow(model_regions)),
    function(i) {
      glue(
        "({format(model_regions$seq_region[i], scientific = FALSE)}::UINTEGER, ",
        "{format(fasta_lengths[i], scientific = FALSE)}::UBIGINT, ",
        "{sql_q(model_regions$name[i])}::VARCHAR)"
      )
    },
    character(1)
  )
  hgvs_region_query <- paste0(
    "SELECT * FROM (VALUES ",
    paste(value_rows, collapse = ", "),
    ") t(seq_region, sequence_length, seq_region_name) ORDER BY seq_region"
  )
}
load_queries <- c(
  if (isTRUE(opt$hgvs)) {
    hgvs_region_query
  } else {
    paste(
      "SELECT seq_region",
      if (complete_coverage) ", sequence_length" else "",
      "FROM", model_query_relation("duckvep_sequence_regions"),
      "ORDER BY seq_region"
    )
  },
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
if (isTRUE(opt$regulatory)) {
  interval_feature_query <- paste(
    "SELECT regulation_feature_index, seq_region, feature_start, feature_end,",
    "feature_kind FROM", model_query_relation("duckvep_regulation_features"),
    "ORDER BY seq_region, feature_start, regulation_feature_index"
  )
  load_options <- c(
    load_options,
    glue("interval_feature_query := {sql_q(interval_feature_query)}")
  )
}
if (complete_coverage) {
  load_options <- c(load_options, "transcript_coverage_complete := TRUE")
}
if (isTRUE(opt$hgvs)) {
  load_options <- c(
    load_options,
    glue("reference_fasta := {sql_q(normalizePath(opt$fasta))}")
  )
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
generate_breakend <- identical(opt$event_mode, "breakend") &&
  !nzchar(source_vcf)
if (!nzchar(source_vcf) && !generate_structural && !generate_breakend) {
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

if (!nzchar(source_vcf) && any(nzchar(c(opt$source_url, opt$source_checksum)))) {
  die("source provenance options require --vcf")
}
source_sha256 <- ""
verified_source_checksum <- ""
source_bytes <- NA_real_
if (nzchar(source_vcf)) {
  source_vcf <- normalizePath(source_vcf)
  source_bytes <- unname(file.info(source_vcf)$size)
  sha_line <- system2("sha256sum", source_vcf, stdout = TRUE, stderr = FALSE)
  if (length(sha_line) != 1L) {
    die("cannot checksum {source_vcf}")
  }
  source_sha256 <- tolower(
    strsplit(trimws(sha_line), "[[:space:]]+")[[1L]][1L]
  )
  if (!grepl("^[0-9a-f]{64}$", source_sha256)) {
    die("sha256sum returned an invalid digest for {source_vcf}")
  }
  if (nzchar(opt$source_checksum)) {
    checksum_parts <- strsplit(opt$source_checksum, ":", fixed = TRUE)[[1L]]
    if (length(checksum_parts) != 2L) {
      die("--source-checksum must be md5:HEX or sha256:HEX")
    }
    checksum_algorithm <- tolower(checksum_parts[[1L]])
    expected_checksum <- tolower(checksum_parts[[2L]])
    actual_checksum <- switch(
      checksum_algorithm,
      md5 = tolower(unname(tools::md5sum(source_vcf))),
      sha256 = source_sha256,
      die("--source-checksum algorithm must be md5 or sha256")
    )
    expected_width <- if (identical(checksum_algorithm, "md5")) 32L else 64L
    if (
      nchar(expected_checksum) != expected_width ||
        !grepl("^[0-9a-f]+$", expected_checksum)
    ) {
      die("--source-checksum has an invalid {checksum_algorithm} digest")
    }
    if (!identical(actual_checksum, expected_checksum)) {
      die(
        "source checksum mismatch: expected {expected_checksum}, ",
        "found {actual_checksum}"
      )
    }
    verified_source_checksum <- paste(checksum_algorithm, actual_checksum, sep = ":")
  }
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
small_source_sql <- ""
if (identical(opt$event_mode, "small")) {
  small_source_sql <- if (isTRUE(opt$split_multiallelic)) {
    glue(
      "SELECT
         v.CHROM AS chrom,
         v.POS::UBIGINT AS position,
         v.REF AS reference,
         a.alternate AS alternate
       FROM read_bcf({sql_q(source_vcf)}, scan_mode := 'sequential') v
       CROSS JOIN UNNEST(v.ALT) WITH ORDINALITY AS a(alternate, alt_index)
       WHERE a.alternate IS NOT NULL
         {chrom_filter}"
    )
  } else {
    glue(
      "SELECT
         CHROM AS chrom,
         POS::UBIGINT AS position,
         REF AS reference,
         ALT[1] AS alternate
       FROM read_bcf({sql_q(source_vcf)}, scan_mode := 'sequential')
       WHERE len(ALT) = 1
         {chrom_filter}"
    )
  }
}
sample_filter <- if (opt$sample_per_shape == 0L) {
  ""
} else {
  glue("WHERE sample_rank <= {opt$sample_per_shape}")
}
small_sample_partition <- if (isTRUE(opt$stratify_raw_allele_length)) {
  "var_type, length_bin, raw_allele_length_bin"
} else {
  "var_type, length_bin"
}
generated_filter <- if (opt$sample_per_shape == 0L) {
  ""
} else {
  glue("WHERE geometry_rank <= {opt$sample_per_shape}")
}
regulation_span_geometries <- if (isTRUE(opt$regulatory)) {
  glue(
    "UNION ALL
     SELECT
       concat(
         CASE f.feature_kind WHEN 1 THEN 'regulatory' ELSE 'motif' END,
         '_exact'
       ),
       f.regulation_feature_index AS source_index, r.name AS chrom,
       f.feature_start AS event_start, f.feature_end AS event_end
     FROM duckvep_regulation_features f
     JOIN duckvep_sequence_regions r USING (seq_region)
     WHERE f.feature_start > 1
       AND f.feature_end <= {model_sequence_length_sql}
       {generated_chrom_filter}
     UNION ALL
     SELECT
       concat(
         CASE f.feature_kind WHEN 1 THEN 'regulatory' ELSE 'motif' END,
         '_containing'
       ),
       f.regulation_feature_index, r.name,
       CASE WHEN f.feature_start > 18 THEN f.feature_start - 17 ELSE 2 END,
       least({model_sequence_length_sql}, f.feature_end + 17)
     FROM duckvep_regulation_features f
     JOIN duckvep_sequence_regions r USING (seq_region)
     WHERE f.feature_start > 1
       AND f.feature_end <= {model_sequence_length_sql}
       {generated_chrom_filter}
     UNION ALL
     SELECT
       concat(
         CASE f.feature_kind WHEN 1 THEN 'regulatory' ELSE 'motif' END,
         '_left_partial'
       ),
       f.regulation_feature_index, r.name,
       CASE WHEN f.feature_start > 18 THEN f.feature_start - 17 ELSE 2 END,
       least(f.feature_end - 1, f.feature_start + 17)
     FROM duckvep_regulation_features f
     JOIN duckvep_sequence_regions r USING (seq_region)
     WHERE f.feature_start > 1 AND f.feature_end > f.feature_start
       AND f.feature_end <= {model_sequence_length_sql}
       {generated_chrom_filter}
     UNION ALL
     SELECT
       concat(
         CASE f.feature_kind WHEN 1 THEN 'regulatory' ELSE 'motif' END,
         '_right_partial'
       ),
       f.regulation_feature_index, r.name,
       greatest(f.feature_start + 1, f.feature_end - 17),
       least({model_sequence_length_sql}, f.feature_end + 17)
     FROM duckvep_regulation_features f
     JOIN duckvep_sequence_regions r USING (seq_region)
     WHERE f.feature_start > 1 AND f.feature_end > f.feature_start
       AND f.feature_end <= {model_sequence_length_sql}
       {generated_chrom_filter}"
  )
} else {
  ""
}

# BND regulation observes two different point transforms: VEP removes the local
# VCF anchor (POS + 1), while the mate coordinate is used verbatim. Generate
# both raw representations at feature starts, midpoints, and ends so a seeded
# campaign cannot appear exact merely because transcript-derived endpoints
# happened not to land in a RegulatoryFeature or MotifFeature.
regulation_breakend_geometries <- if (isTRUE(opt$regulatory)) {
  glue(
    "UNION ALL
     SELECT
       concat(
         CASE f.feature_kind WHEN 1 THEN 'regulatory' ELSE 'motif' END,
         '_', endpoint_role, '_', point_name
       ) AS event_state,
       f.regulation_feature_index AS transcript_index,
       f.seq_region, r.name AS chrom,
       CAST(feature_position + position_shift AS UINTEGER) AS position
     FROM duckvep_regulation_features f
     JOIN duckvep_sequence_regions r USING (seq_region)
     CROSS JOIN LATERAL (VALUES
       ('start', f.feature_start),
       ('mid', (f.feature_start + f.feature_end) // 2),
       ('end', f.feature_end)
     ) feature_points(point_name, feature_position)
     CROSS JOIN (VALUES
       ('local', -1),
       ('mate', 0)
     ) endpoint_roles(endpoint_role, position_shift)
     WHERE feature_position + position_shift > 0
       AND feature_position + position_shift < 4294967295
       {generated_chrom_filter}"
  )
} else {
  ""
}

# Force the two same-chromosome states that shuffled endpoints do not guarantee:
# both transformed endpoint points lie in the same object, or the mate lies in
# the object while the transformed local point is one base before it. The first
# state retains the local base term. The second exposes StructuralVariationOverlap's
# fixed 5000-base endpoint admission and yields the object-level union of
# intergenic_variant with feature_truncation.
regulation_breakend_both_pairs <- if (isTRUE(opt$regulatory)) {
  both_pair_limit <- if (opt$sample_per_shape == 0L) {
    ""
  } else {
    glue(
      "QUALIFY row_number() OVER (
         PARTITION BY f.feature_kind, f.seq_region
         ORDER BY hash(f.regulation_feature_index, {opt$seed}, 73)
       ) <= {opt$sample_per_shape}"
    )
  }
  glue(
    "UNION ALL
     SELECT
       'intra' AS pair_class,
       concat(
         CASE selected.feature_kind
           WHEN 1 THEN 'regulatory' ELSE 'motif'
         END,
         state.local_suffix
       ) AS local_state,
       concat(
         CASE selected.feature_kind
           WHEN 1 THEN 'regulatory' ELSE 'motif'
         END,
         state.mate_suffix
       ) AS mate_state,
       selected.regulation_feature_index AS local_transcript_index,
       selected.regulation_feature_index AS mate_transcript_index,
       selected.seq_region, selected.chrom,
       CAST(
         CAST(
           CASE state.local_anchor
             WHEN 'start' THEN selected.feature_start
             ELSE selected.feature_end
           END AS BIGINT
         ) + state.local_raw_shift AS UINTEGER
       ) AS position,
       selected.seq_region AS mate_seq_region,
       selected.chrom AS mate_chrom,
       CAST(
         (selected.feature_start + selected.feature_end) // 2 AS UINTEGER
       )
         AS mate_position
     FROM (
       SELECT f.*, r.name AS chrom
       FROM duckvep_regulation_features f
       JOIN duckvep_sequence_regions r USING (seq_region)
       WHERE f.feature_start > 2
         AND CAST(f.feature_end AS UBIGINT) + 5001 <=
             {model_sequence_length_sql}
         {generated_chrom_filter}
       {both_pair_limit}
     ) selected
     CROSS JOIN (VALUES
       ('_both_local_start', '_both_mate_mid', 'start', -1::BIGINT),
       ('_mate_only_close_local_before', '_mate_exact_mid', 'start', -2::BIGINT),
       ('_mate_only_exact_cap_after', '_mate_exact_mid', 'end', 4999::BIGINT),
       ('_mate_only_beyond_cap_after', '_mate_exact_mid', 'end', 5000::BIGINT)
     ) state(local_suffix, mate_suffix, local_anchor, local_raw_shift)"
  )
} else {
  ""
}

if (generate_structural) {
  invisible(dbExecute(
    con,
    glue(
      "CREATE OR REPLACE TEMP TABLE duckvep_generated_structural_source AS
       WITH transcripts AS MATERIALIZED (
         SELECT t.*, r.name AS chrom,
                {model_sequence_length_sql} AS sequence_length
         FROM duckvep_transcripts t
         JOIN duckvep_sequence_regions r USING (seq_region)
         WHERE t.transcript_start > 1
           AND t.transcript_end <= {model_sequence_length_sql}
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
           'transcript_exact' AS event_state,
           transcript_index AS source_index, chrom,
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
         {regulation_span_geometries}
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
           'insertion_coding' AS event_state,
           transcript_index AS source_index, chrom,
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
             chrom, event_start, event_end, source_index, {opt$seed}
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
             'duckvep-generated:', event_state, ':', source_index, ':',
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
             'duckvep-generated:', event_state, ':', source_index, ':',
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

if (generate_breakend) {
  invisible(dbExecute(
    con,
    glue(
      "CREATE OR REPLACE TEMP TABLE duckvep_generated_breakend_source AS
       WITH transcripts AS MATERIALIZED (
         SELECT t.*, r.name AS chrom,
                {model_sequence_length_sql} AS sequence_length
         FROM duckvep_transcripts t
         JOIN duckvep_sequence_regions r USING (seq_region)
         WHERE t.transcript_start > 1
           AND t.transcript_end < {model_sequence_length_sql}
           {generated_chrom_filter}
       ), exons AS MATERIALIZED (
         SELECT
           e.*, t.seq_region, t.chrom, t.sequence_length,
           t.transcript_start, t.transcript_end, t.cds_start, t.cds_end,
           lead(e.exon_start) OVER (
             PARTITION BY e.transcript_index ORDER BY e.exon_start, e.exon_end
           ) AS next_exon_start
         FROM duckvep_exons e
         JOIN transcripts t USING (transcript_index)
       ), endpoint_geometries AS (
         SELECT 'transcript_start' AS event_state, transcript_index,
                seq_region, chrom, transcript_start AS position
         FROM transcripts
         UNION ALL
         SELECT 'transcript_mid', transcript_index, seq_region, chrom,
                (transcript_start + transcript_end) // 2
         FROM transcripts
         UNION ALL
         SELECT 'transcript_end', transcript_index, seq_region, chrom,
                transcript_end
         FROM transcripts
         UNION ALL
         SELECT 'upstream_flank', transcript_index, seq_region, chrom,
                transcript_start - 1
         FROM transcripts
         UNION ALL
         SELECT 'downstream_flank', transcript_index, seq_region, chrom,
                transcript_end + 1
         FROM transcripts
         UNION ALL
         SELECT 'cds_start', transcript_index, seq_region, chrom, cds_start
         FROM transcripts WHERE cds_start IS NOT NULL
         UNION ALL
         SELECT 'cds_end', transcript_index, seq_region, chrom, cds_end
         FROM transcripts WHERE cds_end IS NOT NULL
         UNION ALL
         SELECT 'exon_mid', transcript_index, seq_region, chrom,
                (exon_start + exon_end) // 2
         FROM exons
         UNION ALL
         SELECT 'exon_start', transcript_index, seq_region, chrom, exon_start
         FROM exons
         UNION ALL
         SELECT 'before_exon_start', transcript_index, seq_region, chrom,
                exon_start - 1
         FROM exons
         UNION ALL
         SELECT 'exon_end', transcript_index, seq_region, chrom, exon_end
         FROM exons
         UNION ALL
         SELECT 'after_exon_end', transcript_index, seq_region, chrom,
                exon_end + 1
         FROM exons
         UNION ALL
         SELECT 'intron_mid', transcript_index, seq_region, chrom,
                (exon_end + next_exon_start) // 2
         FROM exons
         WHERE next_exon_start IS NOT NULL AND next_exon_start > exon_end + 1
         {regulation_breakend_geometries}
       ), distinct_geometries AS (
         SELECT DISTINCT event_state, transcript_index, seq_region, chrom, position
         FROM endpoint_geometries
         WHERE position > 0 AND position < 4294967295
       ), ranked_geometries AS (
         SELECT *, row_number() OVER (
           PARTITION BY event_state, seq_region
           ORDER BY hash(
             chrom, position, transcript_index, event_state, {opt$seed}
           )
         ) AS geometry_rank
         FROM distinct_geometries
       ), selected_geometries AS (
         SELECT * FROM ranked_geometries {generated_filter}
       ), numbered_endpoints AS (
         SELECT *,
           row_number() OVER (
             PARTITION BY seq_region
             ORDER BY hash(
               event_state, position, transcript_index, {opt$seed}, 41
             )
           ) AS endpoint_rank,
           count(*) OVER (PARTITION BY seq_region) AS endpoint_count
         FROM selected_geometries
       ), chrom_order AS (
         SELECT seq_region,
           coalesce(
             lead(seq_region) OVER ordered_chroms,
             first_value(seq_region) OVER ordered_chroms
           ) AS mate_seq_region
         FROM (SELECT DISTINCT seq_region FROM numbered_endpoints)
         WINDOW ordered_chroms AS (
           ORDER BY seq_region ROWS BETWEEN UNBOUNDED PRECEDING AND UNBOUNDED FOLLOWING
         )
       ), endpoint_pairs AS (
         SELECT
           'intra' AS pair_class,
           l.event_state AS local_state, m.event_state AS mate_state,
           l.transcript_index AS local_transcript_index,
           m.transcript_index AS mate_transcript_index,
           l.seq_region, l.chrom, l.position,
           m.seq_region AS mate_seq_region, m.chrom AS mate_chrom,
           m.position AS mate_position
         FROM numbered_endpoints l
         JOIN numbered_endpoints m
           ON m.seq_region = l.seq_region
          AND m.endpoint_rank = CASE
            WHEN l.endpoint_rank = l.endpoint_count THEN 1
            ELSE l.endpoint_rank + 1
          END
         WHERE l.endpoint_count > 1
         UNION ALL
         SELECT
           'inter', l.event_state, m.event_state,
           l.transcript_index, m.transcript_index,
           l.seq_region, l.chrom, l.position,
           m.seq_region, m.chrom, m.position
         FROM numbered_endpoints l
         JOIN chrom_order c USING (seq_region)
         JOIN numbered_endpoints m
           ON m.seq_region = c.mate_seq_region
          AND m.endpoint_rank = 1 + ((l.endpoint_rank - 1) % m.endpoint_count)
         WHERE c.mate_seq_region <> l.seq_region
         {regulation_breakend_both_pairs}
       ), oriented_pairs AS (
         SELECT p.*, orientation
         FROM endpoint_pairs p
         CROSS JOIN (VALUES
           ('N]M]'), ('N[M['), (']M]N'), ('[M[N')
         ) orientations(orientation)
       )
       SELECT
         concat(
           'duckvep-bnd:', pair_class, ':', local_state, ':', mate_state, ':',
           local_transcript_index, ':', mate_transcript_index, ':',
           chrom, ':', position, ':', mate_chrom, ':', mate_position, ':',
           orientation
         ) AS source_id,
         chrom, position AS raw_position, position AS raw_end,
         'N' AS reference,
         CASE orientation
           WHEN 'N]M]' THEN concat('N]', mate_chrom, ':', mate_position, ']')
           WHEN 'N[M[' THEN concat('N[', mate_chrom, ':', mate_position, '[')
           WHEN ']M]N' THEN concat(']', mate_chrom, ':', mate_position, ']N')
           ELSE concat('[', mate_chrom, ':', mate_position, '[N')
         END AS alternate,
         'BND' AS source_svtype, 'BND' AS structural_type,
         concat(pair_class, '/', local_state, '/', mate_state, '/', orientation)
           AS event_state,
         mate_chrom, mate_position, orientation
       FROM oriented_pairs"
    )
  ))
}

source_bcf_columns <- character()
if (
  identical(opt$event_mode, "structural") &&
    !generate_structural && nzchar(source_vcf)
) {
  source_bcf_columns <- dbGetQuery(
    con,
    glue(
      "DESCRIBE SELECT * FROM read_bcf(
         {sql_q(normalizePath(source_vcf))}, scan_mode := 'sequential'
       )"
    )
  )$column_name
}
source_optional_column <- function(name, missing_sql) {
  if (name %in% source_bcf_columns) name else missing_sql
}

structural_source_sql <- if (generate_structural) {
  paste(
    "SELECT *, NULL::INTEGER[] AS cipos, NULL::INTEGER[] AS ciend,",
    "false AS imprecise FROM duckvep_generated_structural_source"
  )
} else if (generate_breakend) {
  "SELECT * FROM duckvep_generated_breakend_source"
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
       {source_optional_column('INFO_CIPOS', 'NULL::INTEGER[]')} AS cipos,
       {source_optional_column('INFO_CIEND', 'NULL::INTEGER[]')} AS ciend,
       coalesce(
         {source_optional_column('INFO_IMPRECISE', 'false')}, false
       ) AS imprecise,
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
} else if (identical(opt$event_mode, "breakend")) {
  glue(
    "SELECT DISTINCT
       coalesce(nullif(ID, ''), concat('duckvep-bnd:', CHROM, ':', POS, ':', ALT[1]))
         AS source_id,
       CHROM AS chrom,
       POS::UBIGINT AS raw_position,
       POS::UBIGINT AS raw_end,
       REF AS reference,
       ALT[1] AS alternate,
       'BND' AS source_svtype,
       'BND' AS structural_type,
       'source' AS event_state,
       regexp_extract(ALT[1], '([A-Za-z0-9_.-]+):([0-9]+)', 1) AS mate_chrom,
       CAST(regexp_extract(ALT[1], '([A-Za-z0-9_.-]+):([0-9]+)', 2) AS UBIGINT)
         AS mate_position,
       CASE
         WHEN regexp_full_match(
           ALT[1], '[^\\[\\]]+\\][^\\[\\]]+:[0-9]+\\]'
         ) THEN 'N]M]'
         WHEN regexp_full_match(
           ALT[1], '[^\\[\\]]+\\[[^\\[\\]]+:[0-9]+\\['
         ) THEN 'N[M['
         WHEN regexp_full_match(
           ALT[1], '\\][^\\[\\]]+:[0-9]+\\][^\\[\\]]+'
         ) THEN ']M]N'
         WHEN regexp_full_match(
           ALT[1], '\\[[^\\[\\]]+:[0-9]+\\[[^\\[\\]]+'
         ) THEN '[M[N'
         ELSE 'source'
       END AS orientation
     FROM read_bcf({sql_q(normalizePath(source_vcf))}, scan_mode := 'sequential')
     WHERE len(ALT) = 1
       AND (upper(INFO_SVTYPE) = 'BND' OR regexp_matches(ALT[1], '[\\[\\]]'))
       AND regexp_matches(ALT[1], '[A-Za-z0-9_.-]+:[0-9]+')
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
         chrom,
         position,
         reference,
         alternate
       FROM ({small_source_sql}) source
       WHERE reference <> alternate
         AND length(reference) BETWEEN 1 AND {opt$max_allele_length}
         AND length(alternate) BETWEEN 1 AND {opt$max_allele_length}
         AND regexp_full_match(upper(reference), '[ACGTN]+')
         AND regexp_full_match(upper(alternate), '[ACGTN]+')
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
           WHEN abs(length_change) <= 50 THEN concat(if(length_change > 0, '+', '-'), '11..50')
           WHEN abs(length_change) <= 1000 THEN concat(if(length_change > 0, '+', '-'), '51..1000')
           WHEN abs(length_change) <= 10000 THEN concat(if(length_change > 0, '+', '-'), '1001..10000')
           ELSE concat(if(length_change > 0, '+', '-'), '>10000')
         END AS length_bin,
         CASE
           WHEN greatest(length(reference), length(alternate)) = 1 THEN '1'
           WHEN greatest(length(reference), length(alternate)) <= 50 THEN '2..50'
           WHEN greatest(length(reference), length(alternate)) <= 1000 THEN '51..1000'
           WHEN greatest(length(reference), length(alternate)) <= 10000 THEN '1001..10000'
           ELSE '>10000'
         END AS raw_allele_length_bin
       FROM shaped
     ), ranked AS (
       SELECT *, row_number() OVER (
         PARTITION BY {small_sample_partition}
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
} else if (identical(opt$event_mode, "breakend")) {
  glue(
    "CREATE OR REPLACE TEMP TABLE duckvep_breakend_events AS
     WITH source_events AS (
       {structural_source_sql}
     ), prepared AS (
       SELECT
         CASE
           WHEN source_id LIKE 'duckvep-bnd:%' THEN source_id
           ELSE concat(
             'duckvep-bnd:', source_id, ':', chrom, ':', raw_position, ':',
             mate_chrom, ':', mate_position, ':', orientation
           )
         END AS variant_id,
         lr.seq_region, chrom, raw_position AS position,
         mr.seq_region AS mate_seq_region, mate_chrom, mate_position,
         raw_position, raw_end, reference, alternate,
         structural_type, 'UNKNOWN' AS copy_change, source_svtype,
         'bnd' AS var_type, event_state AS length_bin, orientation,
         row_number() OVER (
           PARTITION BY event_state, orientation
           ORDER BY hash(
             chrom, raw_position, mate_chrom, mate_position,
             alternate, source_id, {opt$seed}
           )
         ) AS sample_rank
       FROM source_events s
       JOIN duckvep_sequence_regions lr ON lr.name = s.chrom
       JOIN duckvep_sequence_regions mr ON mr.name = s.mate_chrom
       WHERE raw_position > 0 AND raw_position < 4294967295
         AND mate_position > 0 AND mate_position < 4294967295
     )
     SELECT * EXCLUDE (sample_rank) FROM prepared {sample_filter};

     CREATE OR REPLACE TEMP TABLE duckvep_sample AS
     SELECT
       variant_id, seq_region, chrom, position,
       position AS event_start, position AS event_end,
       raw_position, raw_end, reference, alternate,
       structural_type, copy_change, source_svtype, var_type, length_bin,
       'paired' AS endpoint_role, orientation,
       mate_seq_region, mate_chrom, mate_position
     FROM duckvep_breakend_events
     ORDER BY seq_region, position, variant_id"
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
           'duckvep-sv:', source_id, ':', chrom, ':', raw_position, ':', raw_end, ':',
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
       cipos,
       ciend,
       imprecise,
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
if (identical(opt$event_mode, "breakend")) {
  sample_count <- dbGetQuery(
    con,
    "SELECT count(DISTINCT variant_id) AS n FROM duckvep_sample"
  )$n[[1L]]
}
if (identical(opt$event_mode, "small") && nzchar(opt$eligibility_out)) {
  split_predicate <- if (isTRUE(opt$split_multiallelic)) {
    "TRUE"
  } else {
    "alt_count = 1"
  }
  eligibility <- dbGetQuery(
    con,
    glue(
      "WITH records AS (
         SELECT CHROM AS chrom, POS::UBIGINT AS position,
                REF AS reference, ALT AS alternate_list
         FROM read_bcf({sql_q(source_vcf)}, scan_mode := 'sequential')
         WHERE TRUE {chrom_filter}
       ), record_counts AS (
         SELECT
           count(*)::UBIGINT AS source_records,
           count(*) FILTER (
             WHERE alternate_list IS NULL OR len(alternate_list) = 0
           )::UBIGINT AS records_without_alt,
           count(*) FILTER (
             WHERE len(alternate_list) > 1
           )::UBIGINT AS multiallelic_records
         FROM records
       ), alleles AS (
         SELECT r.chrom, r.position, r.reference, a.alternate,
                len(r.alternate_list) AS alt_count
         FROM records r
         CROSS JOIN UNNEST(r.alternate_list) AS a(alternate)
         WHERE a.alternate IS NOT NULL
       ), classified AS (
         SELECT *,
           reference = alternate AS ref_equal,
           regexp_full_match(alternate, '<[^>]+>') AS symbolic,
           regexp_matches(alternate, '[\\[\\]]') AS breakend,
           alternate = '*' AS spanning_deletion,
           regexp_full_match(upper(reference), '[ACGTN]+') AND
             regexp_full_match(upper(alternate), '[ACGTN]+') AS literal_acgtn
         FROM alleles
       ), allele_counts AS (
         SELECT
           count(*)::UBIGINT AS source_alt_alleles,
           count(*) FILTER (WHERE alt_count > 1)::UBIGINT
             AS multiallelic_alt_alleles,
           count(*) FILTER (WHERE ref_equal)::UBIGINT AS ref_equal_alleles,
           count(*) FILTER (WHERE symbolic)::UBIGINT AS symbolic_alleles,
           count(*) FILTER (WHERE breakend)::UBIGINT AS breakend_alleles,
           count(*) FILTER (WHERE spanning_deletion)::UBIGINT
             AS spanning_deletion_alleles,
           count(*) FILTER (WHERE literal_acgtn)::UBIGINT
             AS literal_acgtn_alleles,
           count(*) FILTER (
             WHERE NOT symbolic AND NOT breakend AND NOT spanning_deletion AND
                   NOT literal_acgtn
           )::UBIGINT AS other_alleles,
           count(*) FILTER (
             WHERE literal_acgtn AND NOT ref_equal AND
                   (length(reference) > {opt$max_allele_length} OR
                    length(alternate) > {opt$max_allele_length})
           )::UBIGINT AS literal_over_requested_limit,
           count(*) FILTER (
             WHERE literal_acgtn AND NOT ref_equal AND
                   (length(reference) > 65535 OR length(alternate) > 65535)
           )::UBIGINT AS literal_over_kernel_limit,
           max(length(reference)) FILTER (WHERE literal_acgtn)::UBIGINT
             AS max_literal_ref_length,
           max(length(alternate)) FILTER (WHERE literal_acgtn)::UBIGINT
             AS max_literal_alt_length
         FROM classified
       ), eligible_pre_model AS (
         SELECT * FROM classified
         WHERE {split_predicate}
           AND NOT ref_equal
           AND literal_acgtn
           AND length(reference) BETWEEN 1 AND {opt$max_allele_length}
           AND length(alternate) BETWEEN 1 AND {opt$max_allele_length}
       ), eligible_model AS (
         SELECT e.*
         FROM eligible_pre_model e
         JOIN duckvep_sequence_regions r ON r.name = e.chrom
       ), eligible_distinct AS (
         SELECT DISTINCT chrom, position, reference, alternate
         FROM eligible_model
       )
       SELECT
         rc.*,
         ac.*,
         {if (isTRUE(opt$split_multiallelic)) '0' else 'ac.multiallelic_alt_alleles'}::UBIGINT
           AS multiallelic_alt_alleles_excluded,
         (SELECT count(*) FROM eligible_pre_model)::UBIGINT
           AS eligible_rows_before_model,
         ((SELECT count(*) FROM eligible_pre_model) -
          (SELECT count(*) FROM eligible_model))::UBIGINT
           AS eligible_rows_outside_model_contigs,
         (SELECT count(*) FROM eligible_model)::UBIGINT AS eligible_model_rows,
         (SELECT count(*) FROM eligible_distinct)::UBIGINT
           AS eligible_distinct_alleles,
         ((SELECT count(*) FROM eligible_model) -
          (SELECT count(*) FROM eligible_distinct))::UBIGINT
           AS duplicate_eligible_rows_removed
       FROM record_counts rc CROSS JOIN allele_counts ac"
    )
  )
  eligibility <- cbind(
    data.frame(
      run_date = as.character(Sys.Date()),
      corpus = opt$corpus,
      event_mode = opt$event_mode,
      source_path = source_vcf,
      source_url = opt$source_url,
      source_version = opt$source_version,
      source_sha256 = source_sha256,
      verified_source_checksum = verified_source_checksum,
      source_bytes = source_bytes,
      max_allele_length = opt$max_allele_length,
      split_multiallelic = isTRUE(opt$split_multiallelic),
      stratify_raw_allele_length = isTRUE(opt$stratify_raw_allele_length),
      sample_per_shape = opt$sample_per_shape,
      selected_events = sample_count,
      stringsAsFactors = FALSE
    ),
    eligibility
  )
  dir.create(dirname(opt$eligibility_out), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(eligibility, opt$eligibility_out, row.names = FALSE, na = "")
}
if (isTRUE(opt$source_audit_only)) {
  invisible(dbGetQuery(
    con,
    glue("SELECT duckvep_model_drop({sql_q(opt$model_name)})")
  ))
  cat(glue("source records audited: {eligibility$source_records[[1L]]}"), "\n")
  cat(glue("eligible distinct alleles: {eligibility$eligible_distinct_alleles[[1L]]}"), "\n")
  cat(glue("eligibility receipt: {opt$eligibility_out}"), "\n")
  quit(save = "no", status = 0L)
}
if (sample_count == 0) {
  die("the sampler found no eligible {opt$event_mode} events in model regions")
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
if (!identical(opt$event_mode, "small")) {
  vcf_header <- c(
    vcf_header,
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">"
  )
}
if (identical(opt$event_mode, "structural")) {
  vcf_header <- c(
    vcf_header,
    "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">",
    "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">",
    "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">"
  )
}
writeLines(
  c(vcf_header, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"),
  vc
)
sample_chrom_sql <- if (identical(opt$event_mode, "breakend")) {
  "SELECT chrom FROM duckvep_sample UNION SELECT mate_chrom FROM duckvep_sample"
} else {
  "SELECT chrom FROM duckvep_sample"
}
missing_fasta_chroms <- dbGetQuery(
  con,
  glue(
    "SELECT DISTINCT chrom FROM ({sample_chrom_sql}) sampled
     EXCEPT SELECT chrom FROM duckvep_oracle_sequence_order
     ORDER BY chrom"
  )
)$chrom
if (length(missing_fasta_chroms) != 0L) {
  close(vc)
  die(
    "sampled chromosome(s) are absent from the FASTA index: ",
    "{paste(missing_fasta_chroms, collapse = ', ')}"
  )
}

# Keep chromosomes contiguous and positions increasing for VEP. The FASTA index
# supplies a deterministic chromosome order; model ordinals are independent.
sample_vcf_query <- if (identical(opt$event_mode, "small")) {
  paste(
    "SELECT v.chrom, v.position, v.variant_id, v.reference, v.alternate",
    "FROM duckvep_sample v",
    "JOIN duckvep_oracle_sequence_order o USING (chrom)",
    "ORDER BY o.fasta_order, v.position, v.reference, v.alternate"
  )
} else if (identical(opt$event_mode, "breakend")) {
  paste(
    "SELECT v.chrom, v.position, v.variant_id, v.reference, v.alternate",
    "FROM duckvep_breakend_events v",
    "JOIN duckvep_oracle_sequence_order o USING (chrom)",
    "ORDER BY o.fasta_order, v.position, v.variant_id"
  )
} else {
  paste(
    "SELECT v.chrom, v.raw_position AS position, v.raw_end, v.source_svtype,",
    "v.cipos, v.ciend, v.imprecise,",
    "v.variant_id, v.reference, v.alternate FROM duckvep_sample v",
    "JOIN duckvep_oracle_sequence_order o USING (chrom)",
    "ORDER BY o.fasta_order, v.event_start, v.event_end, v.structural_type"
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
  } else if (identical(opt$event_mode, "breakend")) {
    rep("SVTYPE=BND", nrow(chunk))
  } else {
    vapply(
      seq_len(nrow(chunk)),
      function(index) {
        fields <- c(
          paste0("END=", chunk$raw_end[[index]]),
          paste0("SVTYPE=", chunk$source_svtype[[index]])
        )
        if (isTRUE(chunk$imprecise[[index]])) {
          fields <- c(fields, "IMPRECISE")
        }
        cipos <- chunk$cipos[[index]]
        if (length(cipos) != 0L && !all(is.na(cipos))) {
          cipos_text <- ifelse(is.na(cipos), ".", as.character(cipos))
          fields <- c(
            fields,
            paste0("CIPOS=", paste(cipos_text, collapse = ","))
          )
        }
        ciend <- chunk$ciend[[index]]
        if (length(ciend) != 0L && !all(is.na(ciend))) {
          ciend_text <- ifelse(is.na(ciend), ".", as.character(ciend))
          fields <- c(
            fields,
            paste0("CIEND=", paste(ciend_text, collapse = ","))
          )
        }
        paste(fields, collapse = ";")
      },
      character(1L)
    )
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
  annotate_function <- if (isTRUE(opt$hgvs)) {
    "duckvep_annotate_hgvs"
  } else {
    "duckvep_annotate"
  }
  glue(
    "{annotate_function}(
       {sql_q(opt$model_name)}, v.seq_region, v.position,
       v.reference, v.alternate, {opt$distance}::UBIGINT
     )"
  )
} else if (identical(opt$event_mode, "breakend")) {
  glue(
    "duckvep_annotate_breakend(
       {sql_q(opt$model_name)}, v.seq_region, v.position,
       v.mate_seq_region, v.mate_position, {opt$distance}::UBIGINT
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
    if (isTRUE(opt$hgvs)) {
      "CASE v.annotation.nmd_prediction_code
         WHEN 1 THEN 'triggering'
         WHEN 2 THEN 'escaping'
         WHEN 3 THEN 'unresolved'
         ELSE 'not_applicable'
       END"
    } else {
      "coalesce(v.annotation.nmd_prediction, 'not_applicable')"
    }
  } else {
    "'not_measured'::VARCHAR"
  }
  engine_consequence_sql <- if (isTRUE(opt$hgvs)) {
    "coalesce((
       SELECT string_agg(b.SO_term, '&' ORDER BY b.SO_term)
       FROM duckvep_so_bits b
       WHERE (
         v.annotation.consequence_mask &
         (1::UBIGINT << b.bit_index)
       ) != 0
     ), '')"
  } else {
    "list_aggregate(
       list_sort(list_distinct(string_split(v.annotation.consequence, '&'))),
       'string_agg', '&'
     )"
  }
  engine_impact_sql <- if (isTRUE(opt$hgvs)) {
    "CASE v.annotation.impact_code
       WHEN 3 THEN 'HIGH'
       WHEN 2 THEN 'MODERATE'
       WHEN 1 THEN 'LOW'
       ELSE 'MODIFIER'
     END"
  } else {
    "v.annotation.impact"
  }
  engine_status_sql <- if (isTRUE(opt$hgvs)) {
    "CASE v.annotation.status_code
       WHEN 0 THEN 'supported'
       ELSE 'unresolved'
     END"
  } else {
    "v.annotation.status"
  }
  engine_reason_sql <- if (isTRUE(opt$hgvs)) {
    "CASE
       WHEN v.annotation.status_code = 0 THEN NULL
       WHEN v.annotation.transcript_index IS NULL AND
            v.annotation.reason_code = 1 THEN 'no_feature_in_loaded_model'
       WHEN v.annotation.reason_code = 2 THEN 'missing_sequence'
       WHEN v.annotation.reason_code = 3 THEN 'ambiguous_sequence'
       WHEN v.annotation.reason_code = 4 THEN 'reference_mismatch'
       WHEN v.annotation.reason_code = 5 THEN 'non_contiguous_cds_edit'
       WHEN v.annotation.reason_code = 6 THEN 'unsupported_compound_consequence'
       WHEN v.annotation.reason_code = 7 THEN 'invalid_model_projection'
       WHEN v.annotation.reason_code = 8 THEN 'internal_capacity_error'
       WHEN v.annotation.reason_code = 9 THEN 'missing_transcript_tail'
       WHEN v.annotation.reason_code = 10 THEN 'missing_transcript_flank'
       ELSE 'invalid_status'
     END"
  } else {
    "v.annotation.reason"
  }
  engine_hgvs_sql <- if (isTRUE(opt$hgvs)) {
    "v.annotation.transcript_hgvs AS hgvsc,
     v.annotation.protein_hgvs AS hgvsp,
     v.annotation.transcript_hgvs_status AS hgvsc_status,
     v.annotation.transcript_hgvs_reason AS hgvsc_reason,
     v.annotation.protein_hgvs_status AS hgvsp_status,
     v.annotation.protein_hgvs_reason AS hgvsp_reason,
     v.annotation.hgvs_shift"
  } else {
    "NULL::VARCHAR AS hgvsc,
     NULL::VARCHAR AS hgvsp,
     'not_measured'::VARCHAR AS hgvsc_status,
     NULL::VARCHAR AS hgvsc_reason,
     'not_measured'::VARCHAR AS hgvsp_status,
     NULL::VARCHAR AS hgvsp_reason,
     NULL::UINTEGER AS hgvs_shift"
  }
  regulation_annotation_sql <- if (isTRUE(opt$regulatory)) {
    glue(
      "UNION ALL
       SELECT
         v.variant_id, v.seq_region, v.position,
         1::UTINYINT AS object_order,
         f.regulation_feature_index AS object_index,
         f.stable_id AS tx,
         list_aggregate(
           list_sort(list_distinct(string_split(v.annotation.consequence, '&'))),
           'string_agg', '&'
         ) AS consequence,
         v.annotation.impact,
         v.annotation.status,
         v.annotation.reason,
         {engine_nmd_sql} AS nmd_prediction,
         NULL::VARCHAR AS hgvsc,
         NULL::VARCHAR AS hgvsp,
         'not_applicable'::VARCHAR AS hgvsc_status,
         NULL::VARCHAR AS hgvsc_reason,
         'not_applicable'::VARCHAR AS hgvsp_status,
         NULL::VARCHAR AS hgvsp_reason,
         NULL::UINTEGER AS hgvs_shift
       FROM annotated v
       JOIN duckvep_regulation_features f
         ON f.regulation_feature_index = v.annotation.regulation_feature_index"
    )
  } else {
    ""
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
         FROM (
           SELECT * FROM duckvep_sample
           ORDER BY seq_region, position, variant_id
         ) v
       ), named AS (
         SELECT
           v.variant_id, v.seq_region, v.position,
           0::UTINYINT AS object_order,
           n.transcript_index AS object_index,
           n.transcript_id AS tx,
           {engine_consequence_sql} AS consequence,
           {engine_impact_sql} AS impact,
           {engine_status_sql} AS status,
           {engine_reason_sql} AS reason,
           {engine_nmd_sql} AS nmd_prediction,
           {engine_hgvs_sql}
         FROM annotated v
         JOIN duckvep_transcript_names n
           ON n.transcript_index = v.annotation.transcript_index
         {regulation_annotation_sql}
       )
       SELECT
         variant_id, tx, consequence, impact, status, reason, nmd_prediction,
         hgvsc, hgvsp, hgvsc_status, hgvsc_reason,
         hgvsp_status, hgvsp_reason, hgvs_shift
       FROM named
       ORDER BY seq_region, position, object_order, object_index"
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
if (identical(opt$event_mode, "breakend")) {
  oracle_details <- c(
    oracle_details,
    "breakend_buffer_size=1",
    "input_order=fasta_index"
  )
} else {
  oracle_details <- c(
    oracle_details,
    glue("vep_buffer_size={opt$vep_buffer_size}")
  )
}
if (identical(opt$event_mode, "structural")) {
  oracle_details <- c(
    oracle_details,
    glue("max_sv_size={format(opt$max_sv_size, scientific = FALSE)}")
  )
}
if (isTRUE(opt$regulatory)) {
  oracle_details <- c(oracle_details, "regulatory=true")
}
if (nzchar(source_sha256)) {
  oracle_details <- c(
    oracle_details,
    glue("source_sha256={source_sha256}"),
    glue("source_bytes={format(source_bytes, scientific = FALSE, trim = TRUE)}")
  )
}
if (nzchar(opt$source_url)) {
  oracle_details <- c(oracle_details, glue("source_url={opt$source_url}"))
}
if (nzchar(opt$source_version)) {
  oracle_details <- c(
    oracle_details,
    glue("source_version={opt$source_version}")
  )
}
if (identical(opt$event_mode, "small")) {
  oracle_details <- c(
    oracle_details,
    glue("max_allele_length={opt$max_allele_length}"),
    glue("split_multiallelic={tolower(isTRUE(opt$split_multiallelic))}"),
    glue(
      "stratify_raw_allele_length=",
      "{tolower(isTRUE(opt$stratify_raw_allele_length))}"
    )
  )
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
  "--force_overwrite",
  "--no_stats"
)
if (identical(opt$event_mode, "breakend")) {
  # VEP's BND interval tree stores mate positions without a chromosome key.
  # Isolate each event so neighboring cross-chromosome records cannot change
  # the executable oracle's transcript set.
  vep_args <- c(vep_args, "--buffer_size", "1")
} else {
  vep_args <- c(
    vep_args,
    "--buffer_size",
    as.character(opt$vep_buffer_size),
    "--fork",
    opt$fork
  )
}
if (identical(opt$event_mode, "structural")) {
  vep_args <- c(
    vep_args,
    "--max_sv_size",
    format(opt$max_sv_size, scientific = FALSE)
  )
}
if (isTRUE(opt$regulatory)) {
  vep_args <- c(vep_args, "--regulatory")
}
if (isTRUE(opt$hgvs)) {
  vep_args <- c(vep_args, "--hgvs")
}
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
       list_contains(
         CAST(json_extract(tc.value, '$.consequence_terms') AS VARCHAR[]),
         'stop_gained'
       ) OR
       list_contains(
         CAST(json_extract(tc.value, '$.consequence_terms') AS VARCHAR[]),
         'frameshift_variant'
       ) OR
       list_contains(
         CAST(json_extract(tc.value, '$.consequence_terms') AS VARCHAR[]),
         'splice_donor_variant'
       ) OR
       list_contains(
         CAST(json_extract(tc.value, '$.consequence_terms') AS VARCHAR[]),
         'splice_acceptor_variant'
       )
     ) THEN 'not_applicable'
     WHEN coalesce(
            json_extract_string(tc.value, '$.duckvep_nmd_cds'),
            'undefined'
          ) = 'undefined' THEN 'unresolved'
     WHEN json_extract_string(tc.value, '$.nmd') =
          'NMD_escaping_variant' THEN 'escaping'
     ELSE 'triggering'
   END"
} else {
  "'not_measured'::VARCHAR"
}
sample_event_relation <- if (identical(opt$event_mode, "breakend")) {
  "duckvep_breakend_events"
} else {
  "duckvep_sample"
}
vep_hgvs_sql <- if (isTRUE(opt$hgvs)) {
  "regexp_replace(
     json_extract_string(tc.value, '$.hgvsc'), '^[^:]+:', ''
   ) AS hgvsc,
   regexp_replace(
     json_extract_string(tc.value, '$.hgvsp'), '^[^:]+:', ''
   ) AS hgvsp,
   CASE
     WHEN json_extract_string(tc.value, '$.hgvsc') IS NULL
       THEN 'not_applicable'
     ELSE 'supported'
   END AS hgvsc_status,
   NULL::VARCHAR AS hgvsc_reason,
   CASE
     WHEN json_extract_string(tc.value, '$.hgvsp') IS NULL
       THEN 'not_applicable'
     ELSE 'supported'
   END AS hgvsp_status,
   NULL::VARCHAR AS hgvsp_reason,
   NULL::UINTEGER AS hgvs_shift"
} else {
  "NULL::VARCHAR AS hgvsc,
   NULL::VARCHAR AS hgvsp,
   'not_measured'::VARCHAR AS hgvsc_status,
   NULL::VARCHAR AS hgvsc_reason,
   'not_measured'::VARCHAR AS hgvsp_status,
   NULL::VARCHAR AS hgvsp_reason,
   NULL::UINTEGER AS hgvs_shift"
}
invisible(dbExecute(
  con,
  glue(
    "CREATE OR REPLACE TEMP TABLE vep_annotation AS
     SELECT
       json_extract_string(to_json(j), '$.id') AS variant_id,
       json_extract_string(tc.value, '$.transcript_id') AS tx,
       coalesce(
         list_aggregate(
           list_sort(list_distinct(
             CAST(json_extract(tc.value, '$.consequence_terms') AS VARCHAR[])
           )),
           'string_agg', '&'
         ),
         ''
       ) AS consequence,
       coalesce(json_extract_string(tc.value, '$.impact'), '') AS impact,
       'oracle'::VARCHAR AS status,
       NULL::VARCHAR AS reason,
       {vep_nmd_sql} AS nmd_prediction,
       {vep_hgvs_sql}
     FROM read_json(
       {sql_q(vep_json)}, format = 'newline_delimited', sample_size = -1
     ) j,
     LATERAL json_each(to_json(j), '$.transcript_consequences') tc
     WHERE json_extract_string(tc.value, '$.transcript_id') IS NOT NULL"
  )
))
if (isTRUE(opt$regulatory)) {
  feature_nmd_prediction <- if (nmd_oracle_enabled) {
    "'not_applicable'::VARCHAR"
  } else {
    "'not_measured'::VARCHAR"
  }
  invisible(dbExecute(
    con,
    glue(
      "INSERT INTO vep_annotation
       WITH records AS (
         SELECT to_json(j) AS record
         FROM read_json(
           {sql_q(vep_json)}, format = 'newline_delimited', sample_size = -1
         ) j
       ), feature_rows AS (
         SELECT
           json_extract_string(record, '$.id') AS variant_id,
           json_extract_string(fc.value, '$.regulatory_feature_id') AS tx,
           CAST(json_extract(fc.value, '$.consequence_terms') AS VARCHAR[])
             AS consequence_terms,
           json_extract_string(fc.value, '$.impact') AS impact
         FROM records,
         LATERAL json_each(record, '$.regulatory_feature_consequences') fc
         WHERE json_extract_string(fc.value, '$.regulatory_feature_id') IS NOT NULL
         UNION ALL
         SELECT
           json_extract_string(record, '$.id'),
           json_extract_string(fc.value, '$.motif_feature_id'),
           CAST(json_extract(fc.value, '$.consequence_terms') AS VARCHAR[]),
           json_extract_string(fc.value, '$.impact')
         FROM records,
         LATERAL json_each(record, '$.motif_feature_consequences') fc
         WHERE json_extract_string(fc.value, '$.motif_feature_id') IS NOT NULL
       )
       SELECT
         variant_id,
         tx,
         coalesce(
           list_aggregate(
             list_sort(list_distinct(consequence_terms)), 'string_agg', '&'
           ),
           ''
         ),
         coalesce(impact, ''),
         'oracle'::VARCHAR,
         NULL::VARCHAR,
         {feature_nmd_prediction},
         NULL::VARCHAR,
         NULL::VARCHAR,
         'not_applicable'::VARCHAR,
         NULL::VARCHAR,
         'not_applicable'::VARCHAR,
         NULL::VARCHAR,
         NULL::UINTEGER
       FROM feature_rows"
    )
  ))
}

# InputBuffer::get_overlapping_vfs can construct the same structural overlap
# more than once. Each overlap also carries every endpoint no farther than 5000
# bases from the object. VEP therefore emits both byte-identical rows and
# distinct allele-level consequence rows for one (variant, object). The public
# conformance unit is the unioned consequence set for that pair, just as it is
# for transcript BNDs. Non-consequence state must still agree before unioning.
vep_state_conflicts <- dbGetQuery(
  con,
  "SELECT variant_id, tx
   FROM vep_annotation
   GROUP BY variant_id, tx
   HAVING count(DISTINCT (
     status, reason, nmd_prediction,
     hgvsc, hgvsp, hgvsc_status, hgvsc_reason, hgvsp_status, hgvsp_reason
   )) > 1
   LIMIT 1"
)
if (nrow(vep_state_conflicts) != 0L) {
  die(
    "VEP emitted conflicting non-consequence state for (",
    "{vep_state_conflicts$variant_id[1]}, {vep_state_conflicts$tx[1]})"
  )
}
vep_duplicate_rows <- dbGetQuery(
  con,
  "SELECT count(*) - count(DISTINCT (
     variant_id, tx, consequence, impact, status, reason, nmd_prediction,
     hgvsc, hgvsp, hgvsc_status, hgvsc_reason, hgvsp_status, hgvsp_reason
   )) AS duplicate_rows
   FROM vep_annotation"
)$duplicate_rows[[1L]]
vep_distinct_object_rows_merged <- dbGetQuery(
  con,
  "SELECT count(DISTINCT (
     variant_id, tx, consequence, impact, status, reason, nmd_prediction,
     hgvsc, hgvsp, hgvsc_status, hgvsc_reason, hgvsp_status, hgvsp_reason
   )) - count(DISTINCT (variant_id, tx)) AS merged_rows
   FROM vep_annotation"
)$merged_rows[[1L]]
invisible(dbExecute(
  con,
  "CREATE OR REPLACE TEMP TABLE vep_annotation AS
   WITH exploded AS (
     SELECT
       variant_id, tx,
       unnest(string_split(consequence, '&')) AS consequence_term,
       CASE impact
         WHEN 'HIGH' THEN 4
         WHEN 'MODERATE' THEN 3
         WHEN 'LOW' THEN 2
         WHEN 'MODIFIER' THEN 1
         ELSE 0
       END AS impact_rank,
       status, reason, nmd_prediction,
       hgvsc, hgvsp, hgvsc_status, hgvsc_reason,
       hgvsp_status, hgvsp_reason, hgvs_shift
     FROM vep_annotation
   )
   SELECT
     variant_id,
     tx,
     coalesce(
       list_aggregate(
         list_sort(list_distinct(list(consequence_term))),
         'string_agg', '&'
       ),
       ''
     ) AS consequence,
     CASE max(impact_rank)
       WHEN 4 THEN 'HIGH'
       WHEN 3 THEN 'MODERATE'
       WHEN 2 THEN 'LOW'
       WHEN 1 THEN 'MODIFIER'
       ELSE ''
     END AS impact,
     any_value(status) AS status,
     any_value(reason) AS reason,
     any_value(nmd_prediction) AS nmd_prediction,
     any_value(hgvsc) AS hgvsc,
     any_value(hgvsp) AS hgvsp,
     any_value(hgvsc_status) AS hgvsc_status,
     any_value(hgvsc_reason) AS hgvsc_reason,
     any_value(hgvsp_status) AS hgvsp_status,
     any_value(hgvsp_reason) AS hgvsp_reason,
     any_value(hgvs_shift) AS hgvs_shift
   FROM exploded
   GROUP BY variant_id, tx"
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
       a.nmd_prediction,
       a.hgvsc,
       a.hgvsp,
       a.hgvsc_status,
       a.hgvsc_reason,
       a.hgvsp_status,
       a.hgvsp_reason,
       a.hgvs_shift
     FROM vep_annotation a JOIN {sample_event_relation} v USING (variant_id)
     UNION ALL
     SELECT
       {sql_q(run_date)}, {sql_q(opt$corpus)}, {sql_q(opt$model_name)},
       {sql_q(oracle_version)}, {sql_q(oracle_build)}, 'duckvep',
       a.variant_id, v.chrom, v.position, v.reference, v.alternate,
       v.var_type, v.length_bin, a.tx, a.consequence, a.impact, a.status, a.reason,
       a.nmd_prediction, a.hgvsc, a.hgvsp, a.hgvsc_status, a.hgvsc_reason,
       a.hgvsp_status, a.hgvsp_reason, a.hgvs_shift
     FROM duckvep_annotation a JOIN {sample_event_relation} v USING (variant_id)"
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

if (isTRUE(opt$hgvs)) {
  invisible(dbExecute(
    con,
    glue(
      "CREATE OR REPLACE TEMP TABLE duckvep_hgvs_pairs AS
       WITH pair_keys AS (
         SELECT variant_id, tx FROM vep_annotation
         UNION
         SELECT variant_id, tx FROM duckvep_annotation
       )
       SELECT
         {sql_q(run_date)} AS run_date,
         {sql_q(opt$corpus)} AS corpus,
         {sql_q(opt$model_name)} AS model,
         {sql_q(oracle_version)} AS oracle_version,
         {sql_q(oracle_build)} AS oracle_build,
         k.variant_id,
         k.tx,
         s.chrom,
         s.position AS pos,
         s.reference AS ref,
         s.alternate AS alt,
         s.var_type,
         s.length_bin,
         v.consequence AS vep_consequence,
         d.consequence AS duckvep_consequence,
         v.hgvsc AS vep_hgvsc,
         d.hgvsc AS duckvep_hgvsc,
         coalesce(d.hgvsc_status, 'missing') AS duckvep_hgvsc_status,
         d.hgvsc_reason AS duckvep_hgvsc_reason,
         v.hgvsp AS vep_hgvsp,
         d.hgvsp AS duckvep_hgvsp,
         coalesce(d.hgvsp_status, 'missing') AS duckvep_hgvsp_status,
         d.hgvsp_reason AS duckvep_hgvsp_reason,
         d.hgvs_shift AS duckvep_hgvs_shift,
         CASE
           WHEN v.variant_id IS NULL THEN 'engine_extra_row'
           WHEN d.variant_id IS NULL THEN 'engine_missing_row'
           WHEN d.hgvsc_status = 'unresolved' THEN 'engine_unresolved'
           WHEN v.hgvsc IS NULL AND d.hgvsc IS NULL THEN 'both_absent'
           WHEN v.hgvsc IS NULL THEN 'engine_extra'
           WHEN d.hgvsc IS NULL THEN 'engine_missing'
           WHEN v.hgvsc = d.hgvsc THEN 'match'
           ELSE 'mismatch'
         END AS hgvsc_comparison,
         CASE
           WHEN v.variant_id IS NULL THEN 'engine_extra_row'
           WHEN d.variant_id IS NULL THEN 'engine_missing_row'
           WHEN d.hgvsp_status = 'unresolved' THEN 'engine_unresolved'
           WHEN v.hgvsp IS NULL AND d.hgvsp IS NULL THEN 'both_absent'
           WHEN v.hgvsp IS NULL THEN 'engine_extra'
           WHEN d.hgvsp IS NULL THEN 'engine_missing'
           WHEN v.hgvsp = d.hgvsp THEN 'match'
           ELSE 'mismatch'
         END AS hgvsp_comparison
       FROM pair_keys k
       LEFT JOIN vep_annotation v USING (variant_id, tx)
       LEFT JOIN duckvep_annotation d USING (variant_id, tx)
       JOIN {sample_event_relation} s USING (variant_id)"
    )
  ))
  hgvs_pair_count <- dbGetQuery(
    con,
    "SELECT count(*) AS n FROM duckvep_hgvs_pairs"
  )$n[[1L]]
  if (hgvs_pair_count == 0) {
    die("HGVS differential produced no transcript pairs")
  }
  dir.create(
    dirname(opt$hgvs_pairs_out),
    recursive = TRUE,
    showWarnings = FALSE
  )
  invisible(dbExecute(
    con,
    glue(
      "COPY duckvep_hgvs_pairs TO {sql_q(opt$hgvs_pairs_out)}
       (FORMAT parquet, COMPRESSION zstd)"
    )
  ))
  hgvs_summary <- dbGetQuery(
    con,
    "WITH comparisons AS (
       SELECT
         'hgvsc'::VARCHAR AS metric,
         hgvsc_comparison AS comparison,
         var_type,
         length_bin,
         coalesce(vep_consequence, '(no_vep_emission)') AS consequence_class,
         duckvep_hgvsc_status AS engine_status,
         coalesce(duckvep_hgvsc_reason, '') AS engine_reason
       FROM duckvep_hgvs_pairs
       UNION ALL
       SELECT
         'hgvsp', hgvsp_comparison, var_type, length_bin,
         coalesce(vep_consequence, '(no_vep_emission)'),
         duckvep_hgvsp_status,
         coalesce(duckvep_hgvsp_reason, '')
       FROM duckvep_hgvs_pairs
     )
     SELECT
       metric,
       comparison,
       var_type,
       length_bin,
       consequence_class,
       engine_status,
       engine_reason,
       count(*) AS n
     FROM comparisons
     GROUP BY ALL
     ORDER BY metric, comparison, n DESC, var_type, length_bin,
              consequence_class, engine_status, engine_reason"
  )
  dir.create(dirname(opt$hgvs_out), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(hgvs_summary, opt$hgvs_out, row.names = FALSE)
  hgvs_discordance_count <- dbGetQuery(
    con,
    "SELECT count(*) AS n
     FROM (
       SELECT hgvsc_comparison AS comparison FROM duckvep_hgvs_pairs
       UNION ALL
       SELECT hgvsp_comparison FROM duckvep_hgvs_pairs
     )
     WHERE comparison IN (
       'engine_unresolved', 'mismatch', 'engine_extra', 'engine_missing',
       'engine_extra_row', 'engine_missing_row'
     )"
  )$n[[1L]]
}

counts <- dbGetQuery(
  con,
  "SELECT source, count(*) AS annotation_rows
   FROM duckvep_annotation_dump GROUP BY source ORDER BY source"
)
invisible(dbGetQuery(
  con,
  glue("SELECT duckvep_model_drop({sql_q(opt$model_name)})")
))

cat(glue("sampled variants: {sample_count}"), "\n", sep = "")
for (i in seq_len(nrow(counts))) {
  cat(
    glue("{counts$source[i]} annotation rows: {counts$annotation_rows[i]}"),
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
cat(glue("VEP duplicate object rows collapsed: {vep_duplicate_rows}"), "\n", sep = "")
cat(
  glue(
    "VEP distinct allele-level object rows unioned: ",
    "{vep_distinct_object_rows_merged}"
  ),
  "\n",
  sep = ""
)
cat(glue("annotations: {opt$annotations_out}"), "\n", sep = "")
if (isTRUE(opt$hgvs)) {
  cat(glue("HGVS pairs: {opt$hgvs_pairs_out}"), "\n", sep = "")
  cat(glue("HGVS summary: {opt$hgvs_out}"), "\n", sep = "")
}
if (nzchar(opt$eligibility_out)) {
  cat(glue("eligibility receipt: {opt$eligibility_out}"), "\n", sep = "")
}

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
if (isTRUE(opt$hgvs) && hgvs_discordance_count != 0L &&
    !isTRUE(opt$allow_hgvs_discordance)) {
  die(glue(
    "HGVS differential found {hgvs_discordance_count} unresolved, ",
    "mismatch, missing, or extra comparisons; inspect {opt$hgvs_pairs_out} ",
    "and {opt$hgvs_out}. Use --allow-hgvs-discordance only for an ",
    "explicitly investigative run."
  ))
}
