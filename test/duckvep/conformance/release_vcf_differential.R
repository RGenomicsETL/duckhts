#!/usr/bin/env Rscript

# Compare DuckVEP consequences with the lossless VE rows in an official Ensembl
# Variation release VCF. This audits that published release product; it is not a
# substitute for comparison with the VEP executable/cache combination. VE retains
# Consequence|Index|Feature_type|Feature_id, and Index identifies the original GVF
# Variant_seq and its corresponding VCF ALT. CSQ is a further lossy presentation:
# gvf2vcf.pl overwrites repeated terms for one allele/feature while building it.

suppressMessages({
  library(DBI)
  library(duckdb)
  library(glue)
  library(optparse)
})

root <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)[[1L]]
root <- normalizePath(root, mustWork = TRUE)
source(file.path(root, "scripts", "duckvep_evidence.R"), local = TRUE)

op <- OptionParser()
op <- add_option(op, "--input", help = "official Ensembl release VCF.gz")
op <- add_option(op, "--database", help = "receipt-hashed DuckVEP model database")
op <- add_option(
  op,
  "--extension",
  default = file.path(root, "build", "release", "duckhts.duckdb_extension")
)
op <- add_option(op, "--release", default = "116")
op <- add_option(op, "--assembly", default = "GRCh38")
op <- add_option(op, "--chromosome", default = "22")
op <- add_option(op, "--source-url", dest = "source_url", default = "")
op <- add_option(
  op,
  "--source-checksum",
  dest = "source_checksum",
  default = "",
  help = "expected sha256:HEX checksum"
)
op <- add_option(
  op,
  "--variation-source-revision",
  dest = "variation_source_revision",
  default = "2fb834b987ede3824e200197a838ce11e91aeb4b"
)
op <- add_option(op, "--threads", type = "integer", default = 1L)
op <- add_option(
  op,
  "--limit",
  type = "double",
  default = 0,
  help = "maximum eligible SNV alleles; zero keeps the complete chromosome"
)
op <- add_option(op, "--output", help = "one-row CSV receipt")
op <- add_option(
  op,
  "--differences-out",
  dest = "differences_out",
  default = ""
)
op <- add_option(
  op,
  "--skip-extension-build",
  dest = "skip_extension_build",
  action = "store_true",
  default = FALSE
)
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
required <- c("input", "database", "output")
missing <- required[vapply(required, function(x) {
  is.null(opt[[x]]) || !nzchar(opt[[x]])
}, logical(1))]
if (length(missing)) die("missing required option(s): {paste(missing, collapse = ', ')}")
if (!file.exists(opt$input)) die("input does not exist: {opt$input}")
if (!file.exists(opt$database)) die("database does not exist: {opt$database}")
if (opt$threads < 1L) die("--threads must be positive")
if (!is.finite(opt$limit) || opt$limit < 0 || opt$limit != floor(opt$limit)) {
  die("--limit must be a whole non-negative number")
}

input <- normalizePath(opt$input, mustWork = TRUE)
database <- normalizePath(opt$database, mustWork = TRUE)
output <- normalizePath(opt$output, mustWork = FALSE)
differences_out <- if (nzchar(opt$differences_out)) {
  normalizePath(opt$differences_out, mustWork = FALSE)
} else {
  ""
}
allowed_outputs <- c(output, differences_out)
revision <- duckvep_evidence_revision(root)

if (isTRUE(opt$skip_extension_build)) {
  extension <- normalizePath(opt$extension, mustWork = TRUE)
  extension_binding <- "diagnostic_existing_extension"
  extension_sha256 <- duckvep_evidence_sha256(extension)
  checked_evidence <- FALSE
} else {
  duckvep_evidence_assert_checkout(
    root,
    revision,
    allowed_outputs,
    context = "official Ensembl release differential"
  )
  build <- duckvep_evidence_build_extension(
    root,
    opt$extension,
    revision,
    allowed_outputs
  )
  extension <- build$path
  extension_binding <- build$binding
  extension_sha256 <- build$sha256
  checked_evidence <- TRUE
}

input_sha256 <- duckvep_evidence_sha256(input)
if (nzchar(opt$source_checksum)) {
  expected <- tolower(sub("^sha256:", "", opt$source_checksum))
  if (!grepl("^[0-9a-f]{64}$", expected)) {
    die("--source-checksum must be sha256:HEX")
  }
  if (!identical(input_sha256, expected)) {
    die("input SHA-256 does not match --source-checksum")
  }
}
model_database_sha256 <- duckvep_evidence_sha256(database)

drv <- duckdb(config = list(allow_unsigned_extensions = "true"))
con <- dbConnect(drv, dbdir = database)
on.exit({
  try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
}, add = TRUE)

sql_q <- function(x) as.character(dbQuoteString(con, x))
dbExecute(con, glue("LOAD {sql_q(extension)}"))
dbExecute(con, glue("SET threads = {opt$threads}"))

required_relations <- c(
  "bench_regions", "bench_transcripts", "bench_exons", "model_regions",
  "model_transcripts", "model_receipt"
)
present_relations <- dbGetQuery(
  con,
  "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'"
)$table_name
missing_relations <- setdiff(required_relations, present_relations)
if (length(missing_relations)) {
  die("model database lacks relation(s): {paste(missing_relations, collapse = ', ')}")
}

region <- dbGetQuery(
  con,
  glue(
    "SELECT seq_region, sequence_length
       FROM model_regions
      WHERE seq_region_name = {sql_q(opt$chromosome)}"
  )
)
if (nrow(region) != 1L) {
  die("model must contain exactly one sequence region named {opt$chromosome}")
}
seq_region <- region$seq_region[[1L]]

model_receipt <- dbGetQuery(con, "SELECT * FROM model_receipt")
if (nrow(model_receipt) != 1L) die("model_receipt must contain exactly one row")
if (!identical(as.character(model_receipt$source_version[[1L]]), opt$release)) {
  die("model release does not match --release")
}
if (!identical(as.character(model_receipt$assembly[[1L]]), opt$assembly)) {
  die("model assembly does not match --assembly")
}

load_options <- c("transcript_coverage_complete := TRUE")
if ("bench_mature_mirna" %in% present_relations) {
  load_options <- c(
    load_options,
    paste0(
      "mature_mirna_query := ",
      sql_q(paste(
        "SELECT transcript_index, mature_mirna_start, mature_mirna_end",
        "FROM bench_mature_mirna",
        "ORDER BY transcript_index, mature_mirna_start"
      ))
    )
  )
}
if ("bench_peptide_edits" %in% present_relations) {
  load_options <- c(
    load_options,
    paste0(
      "peptide_edit_query := ",
      sql_q(paste(
        "SELECT transcript_index, protein_position, alternate_amino_acid",
        "FROM bench_peptide_edits",
        "ORDER BY transcript_index, protein_position"
      ))
    )
  )
}
model_name <- "official_release_differential"
load_sql <- glue(
  "SELECT loaded FROM duckvep_model_load(
     {sql_q(model_name)},
     {sql_q('SELECT seq_region, sequence_length FROM bench_regions ORDER BY seq_region')},
     {sql_q(paste(
       'SELECT transcript_index, seq_region, transcript_start, transcript_end,',
       'strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence,',
       'codon_table, pre_cds_sequence, post_cds_sequence',
       'FROM bench_transcripts ORDER BY transcript_index'
     ))},
     {sql_q(paste(
       'SELECT transcript_index, exon_start, exon_end, exon_cdna_start,',
       'exon_cdna_end, phase, end_phase FROM bench_exons',
       'ORDER BY transcript_index, exon_cdna_start'
     ))},
     {paste(load_options, collapse = ', ')})"
)
model_load_time <- system.time({
  loaded <- dbGetQuery(con, load_sql)$loaded
})[["elapsed"]]
if (length(loaded) != 1L || !isTRUE(loaded[[1L]])) {
  die("DuckVEP model load failed")
}

source_counts_time <- system.time({
  source_counts <- dbGetQuery(
    con,
    glue(
      "WITH records AS (
         SELECT ALT
           FROM read_bcf_v2({sql_q(input)}, include_format := FALSE)
       )
       SELECT count(*)::UBIGINT AS source_records,
              coalesce(sum(len(ALT)), 0)::UBIGINT AS source_alt_alleles
         FROM records"
    )
  )
})[["elapsed"]]

limit_sql <- if (opt$limit > 0) {
  glue("LIMIT {format(opt$limit, scientific = FALSE, trim = TRUE)}")
} else {
  ""
}
stage_time <- system.time({
  dbExecute(
    con,
    glue(
      "CREATE TEMP TABLE release_events AS
       WITH records AS (
         SELECT POS::UBIGINT AS position,
                REF AS reference,
                ALT,
                INFO_VE
           FROM read_bcf_v2({sql_q(input)}, include_format := FALSE)
          WHERE CHROM = {sql_q(opt$chromosome)}
       ),
       alleles AS (
         SELECT position, reference,
                alt.ordinality::UINTEGER - 1 AS alt_index,
                alt.alternate,
                INFO_VE
           FROM records,
                LATERAL unnest(ALT) WITH ORDINALITY AS alt(alternate, ordinality)
          WHERE length(reference) = 1
            AND length(alt.alternate) = 1
            AND regexp_full_match(reference, '[ACGTNacgtn]')
            AND regexp_full_match(alt.alternate, '[ACGTNacgtn]')
       ),
       selected AS (
         SELECT *
           FROM alleles
          ORDER BY position, alt_index, reference, alternate
          {limit_sql}
       )
       SELECT (row_number() OVER (
                  ORDER BY position, alt_index, reference, alternate
              ) - 1)::UBIGINT AS event_index,
              alt_index,
              {seq_region}::UINTEGER AS seq_region,
              position, upper(reference) AS reference,
              upper(alternate) AS alternate,
              INFO_VE
         FROM selected
        ORDER BY position, alt_index, reference, alternate, event_index
       "
    )
  )
  dbExecute(
    con,
    "CREATE TEMP TABLE release_event_input AS
     SELECT event_index, seq_region, position, reference, alternate,
            NULL::UBIGINT AS end_position,
            NULL::VARCHAR AS structural_type,
            NULL::VARCHAR AS copy_change,
            NULL::UINTEGER AS mate_seq_region,
            NULL::UBIGINT AS mate_position
       FROM release_events
      ORDER BY position, event_index"
  )
})[["elapsed"]]

eligible_events <- dbGetQuery(
  con,
  "SELECT count(*)::UBIGINT AS n FROM release_events"
)$n[[1L]]
if (eligible_events == 0) die("the selected release stratum contains no SNVs")

oracle_time <- system.time({
  dbExecute(
    con,
    "CREATE TEMP TABLE release_oracle_terms AS
     WITH entries AS (
       SELECT event_index, alt_index,
              string_split(raw_ve, '|') AS fields
         FROM release_events,
              LATERAL unnest(INFO_VE) AS u(raw_ve)
     ),
     selected AS (
       SELECT event_index,
              CASE WHEN len(fields) = 1 THEN '' ELSE fields[4] END AS feature,
              fields[1] AS term
         FROM entries
        WHERE (len(fields) = 1 AND fields[1] = 'intergenic_variant')
           OR (len(fields) = 4 AND try_cast(fields[2] AS UINTEGER) = alt_index)
     )
     SELECT DISTINCT event_index, feature, term FROM selected"
  )
  dbExecute(
    con,
    "CREATE TEMP TABLE release_oracle_pairs AS
     SELECT event_index, feature,
            string_agg(term, '&' ORDER BY term) AS consequence
       FROM release_oracle_terms
      GROUP BY event_index, feature"
  )
})[["elapsed"]]

annotation_time <- system.time({
  dbExecute(
    con,
    glue(
      "CREATE TEMP TABLE release_duck_rows AS
       SELECT a.event_index,
              coalesce(
                t.transcript_stable_id,
                CASE WHEN a.transcript_index IS NULL THEN ''
                     ELSE '<missing-transcript-id>' END
              ) AS feature,
              a.consequence,
              a.duckvep_status,
              a.duckvep_reason
         FROM duckvep_annotate(
                'release_event_input', {sql_q(model_name)},
                upstream_distance := 0, downstream_distance := 0,
                rich := TRUE
              ) AS a
         LEFT JOIN model_transcripts AS t USING (transcript_index)"
    )
  )
  dbExecute(
    con,
    "CREATE TEMP TABLE release_duck_terms AS
     SELECT DISTINCT event_index, feature, term
       FROM release_duck_rows,
            LATERAL unnest(
              CASE WHEN consequence IS NULL THEN ['<NULL>']
                   ELSE string_split(consequence, '&') END
            ) AS u(term)"
  )
  dbExecute(
    con,
    "CREATE TEMP TABLE release_duck_pairs AS
     SELECT event_index, feature,
            string_agg(term, '&' ORDER BY term) AS consequence
       FROM release_duck_terms
      GROUP BY event_index, feature"
  )
})[["elapsed"]]

comparison_time <- system.time({
  dbExecute(
    con,
    "CREATE TEMP TABLE release_comparison AS
     SELECT coalesce(o.event_index, d.event_index) AS event_index,
            coalesce(o.feature, d.feature) AS feature,
            o.consequence AS oracle_consequence,
            d.consequence AS duckvep_consequence,
            CASE
              WHEN o.event_index IS NULL THEN 'extra'
              WHEN d.event_index IS NULL THEN 'missing'
              WHEN o.consequence = d.consequence THEN 'exact'
              ELSE 'mismatch'
            END AS comparison
       FROM release_oracle_pairs AS o
       FULL OUTER JOIN release_duck_pairs AS d
         USING (event_index, feature)"
  )
})[["elapsed"]]

counts <- dbGetQuery(
  con,
  "SELECT count(*)::UBIGINT AS compared_pairs,
          count(*) FILTER (comparison = 'exact')::UBIGINT AS exact_pairs,
          count(*) FILTER (comparison = 'mismatch')::UBIGINT AS mismatched_pairs,
          count(*) FILTER (comparison = 'missing')::UBIGINT AS missing_pairs,
          count(*) FILTER (comparison = 'extra')::UBIGINT AS extra_pairs,
          count(DISTINCT event_index)::UBIGINT AS compared_events
     FROM release_comparison"
)
pair_counts <- dbGetQuery(
  con,
  "SELECT
     (SELECT count(*) FROM release_oracle_pairs)::UBIGINT AS oracle_pairs,
     (SELECT count(*) FROM release_duck_pairs)::UBIGINT AS duckvep_pairs,
     (SELECT count(*) FROM release_oracle_terms)::UBIGINT AS oracle_terms,
     (SELECT count(*) FROM release_duck_terms)::UBIGINT AS duckvep_terms,
     (SELECT count(*) FROM release_events e
       WHERE NOT EXISTS (
         SELECT 1 FROM release_oracle_pairs o WHERE o.event_index = e.event_index
       ))::UBIGINT AS events_without_oracle,
     (SELECT count(*) FROM release_duck_rows
       WHERE duckvep_status <> 'supported')::UBIGINT AS unsupported_duckvep_rows"
)
exact_events <- dbGetQuery(
  con,
  "SELECT count(*)::UBIGINT AS n
     FROM (
       SELECT event_index
         FROM release_comparison
        GROUP BY event_index
       HAVING bool_and(comparison = 'exact')
     )"
)$n[[1L]]

if (nzchar(differences_out)) {
  dir.create(dirname(differences_out), recursive = TRUE, showWarnings = FALSE)
  differences <- dbGetQuery(
    con,
    "SELECT * FROM release_comparison
      WHERE comparison <> 'exact'
      ORDER BY event_index, feature
      LIMIT 1000"
  )
  write.csv(differences, differences_out, row.names = FALSE, na = "")
}

loaded_version <- dbGetQuery(
  con,
  "SELECT extension_version FROM duckdb_extensions()
    WHERE extension_name = 'duckhts' AND loaded"
)$extension_version[[1L]]
duckdb_version <- dbGetQuery(con, "SELECT version() AS version")$version[[1L]]

receipt <- data.frame(
  run_date = as.character(Sys.Date()),
  source_revision = revision,
  checked_evidence = checked_evidence,
  extension_build_binding = extension_binding,
  extension_sha256 = extension_sha256,
  loaded_extension_version = loaded_version,
  duckdb_version = duckdb_version,
  release = opt$release,
  assembly = opt$assembly,
  chromosome = opt$chromosome,
  event_shape = "literal_snv",
  threads = opt$threads,
  requested_limit = format(opt$limit, scientific = FALSE, trim = TRUE),
  source_url = opt$source_url,
  source_sha256 = input_sha256,
  source_records = source_counts$source_records[[1L]],
  source_alt_alleles = source_counts$source_alt_alleles[[1L]],
  eligible_events = eligible_events,
  oracle_pairs = pair_counts$oracle_pairs[[1L]],
  duckvep_pairs = pair_counts$duckvep_pairs[[1L]],
  oracle_terms = pair_counts$oracle_terms[[1L]],
  duckvep_terms = pair_counts$duckvep_terms[[1L]],
  compared_pairs = counts$compared_pairs[[1L]],
  exact_pairs = counts$exact_pairs[[1L]],
  mismatched_pairs = counts$mismatched_pairs[[1L]],
  missing_pairs = counts$missing_pairs[[1L]],
  extra_pairs = counts$extra_pairs[[1L]],
  compared_events = counts$compared_events[[1L]],
  exact_events = exact_events,
  events_without_oracle = pair_counts$events_without_oracle[[1L]],
  unsupported_duckvep_rows = pair_counts$unsupported_duckvep_rows[[1L]],
  model_database_sha256 = model_database_sha256,
  model_logical_sha256 = as.character(model_receipt$model_sha256[[1L]]),
  model_source_manifest_sha256 =
    as.character(model_receipt$source_manifest_sha256[[1L]]),
  variation_source_revision = opt$variation_source_revision,
  oracle_field = "VE",
  comparison_authority =
    "Ensembl Variation release product; not the VEP executable oracle",
  allele_mapping = "VE zero-based Index identifies VCF ALT ordinal",
  upstream_distance = 0L,
  downstream_distance = 0L,
  model_load_seconds = model_load_time,
  source_count_seconds = source_counts_time,
  stage_seconds = stage_time,
  oracle_seconds = oracle_time,
  annotation_seconds = annotation_time,
  comparison_seconds = comparison_time,
  stringsAsFactors = FALSE
)

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(receipt, output, row.names = FALSE, na = "")
if (isTRUE(checked_evidence)) {
  duckvep_evidence_assert_checkout(
    root,
    revision,
    allowed_outputs,
    context = "official Ensembl release differential"
  )
}

cat(glue("eligible SNV alleles: {format(eligible_events, big.mark = ',')}"), "\n")
cat(glue(
  "exact transcript/object pairs: {format(counts$exact_pairs[[1L]], big.mark = ',')} / ",
  "{format(counts$compared_pairs[[1L]], big.mark = ',')}"
), "\n", sep = "")
cat(glue(
  "mismatch={counts$mismatched_pairs[[1L]]}; missing={counts$missing_pairs[[1L]]}; ",
  "extra={counts$extra_pairs[[1L]]}; unsupported={pair_counts$unsupported_duckvep_rows[[1L]]}"
), "\n", sep = "")
cat(glue("receipt: {output}"), "\n")

if (
  counts$mismatched_pairs[[1L]] != 0 ||
    counts$missing_pairs[[1L]] != 0 ||
    counts$extra_pairs[[1L]] != 0 ||
    pair_counts$events_without_oracle[[1L]] != 0 ||
    pair_counts$unsupported_duckvep_rows[[1L]] != 0
) {
  quit(status = 1L)
}
