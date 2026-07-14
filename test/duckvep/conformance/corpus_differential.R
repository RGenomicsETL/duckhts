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
  "--gff",
  default = file.path(root, "test", "data", "duckvep", "minimal.gff3")
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
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
if (
  !requireNamespace("blit", quietly = TRUE) ||
    utils::packageVersion("blit") < "0.2.0.9000"
) {
  die(
    "the current WangLabCSU/blit checkout is required ",
    "(blit >= 0.2.0.9000)"
  )
}
required_files <- c(opt$gff, opt$fasta, opt$extension)
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
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
if (!nzchar(opt$annotations_out)) {
  opt$annotations_out <- file.path(
    results_dir,
    glue("{opt$corpus}_annotations.parquet")
  )
}

temporary_files <- character()
cleanup <- function() unlink(temporary_files, recursive = FALSE, force = TRUE)
on.exit(cleanup(), add = TRUE)

drv <- duckdb(
  dbdir = opt$database,
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
  "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'"
)$table_name
missing_relations <- setdiff(needed_relations, present_relations)
if (length(missing_relations) != 0L) {
  die(
    "model database is missing relation(s): {paste(missing_relations, collapse = ', ')}"
  )
}

load_queries <- c(
  "SELECT seq_region FROM duckvep_sequence_regions ORDER BY seq_region",
  paste(
    "SELECT transcript_index, seq_region, transcript_start, transcript_end,",
    "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table,",
    "pre_cds_sequence, post_cds_sequence",
    "FROM duckvep_transcripts ORDER BY seq_region, transcript_start, transcript_index"
  ),
  paste(
    "SELECT transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end,",
    "phase, end_phase FROM duckvep_exons",
    "ORDER BY transcript_index, exon_cdna_start"
  )
)
loaded <- dbGetQuery(
  con,
  glue(
    "SELECT loaded FROM duckvep_model_load(
       {sql_q(opt$model_name)}, {sql_q(load_queries[1])},
       {sql_q(load_queries[2])}, {sql_q(load_queries[3])})"
  )
)$loaded
if (length(loaded) != 1L || !isTRUE(loaded[[1L]])) {
  die("DuckVEP model load failed")
}

source_vcf <- opt$vcf
if (!nzchar(source_vcf)) {
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
}
if (!file.exists(source_vcf)) {
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
sample_filter <- if (opt$sample_per_shape == 0L) {
  ""
} else {
  glue("WHERE sample_rank <= {opt$sample_per_shape}")
}

invisible(dbExecute(
  con,
  glue(
    "CREATE OR REPLACE TABLE duckvep_sample AS
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
))
sample_count <- dbGetQuery(con, "SELECT count(*) AS n FROM duckvep_sample")$n[[
  1L
]]
if (sample_count == 0) {
  die("the sampler found no biallelic A/C/G/T variants in model regions")
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
writeLines(
  c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"),
  vc
)
res <- dbSendQuery(
  con,
  paste(
    "SELECT chrom, position, variant_id, reference, alternate",
    "FROM duckvep_sample ORDER BY seq_region, position, reference, alternate"
  )
)
repeat {
  chunk <- dbFetch(res, n = 100000L)
  if (nrow(chunk) == 0L) {
    break
  }
  writeLines(
    paste(
      chunk$chrom,
      chunk$position,
      chunk$variant_id,
      chunk$reference,
      chunk$alternate,
      ".",
      "PASS",
      ".",
      sep = "\t"
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
gff_for_vep <- stage_gff(opt$gff)

engine_time <- system.time({
  dbExecute(
    con,
    glue(
      "CREATE OR REPLACE TABLE duckvep_annotation AS
       WITH annotated AS (
         SELECT
           v.variant_id,
           v.seq_region,
           v.position,
           unnest(duckvep_annotate(
             {sql_q(opt$model_name)}, v.seq_region, v.position,
             v.reference, v.alternate, {opt$distance}::UBIGINT
           )) AS annotation
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
         v.annotation.reason
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
oracle_build <- paste(
  glue("core={component_version('ensembl')}"),
  glue("variation={component_version('ensembl-variation')}"),
  glue("vep={oracle_version}"),
  sep = ";"
)

vep_json <- tempfile(fileext = ".json")
temporary_files <- c(temporary_files, vep_json)
vep_time <- system.time({
  rc <- vep_command(
    "-i",
    sample_vcf,
    "--gff",
    gff_for_vep,
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
  ) |>
    blit::cmd_run(stdout = "", stderr = "", stdin = NULL, verbose = FALSE)
})
if (rc != 0L || !file.exists(vep_json) || file.info(vep_json)$size == 0) {
  die("VEP failed with exit status {rc}")
}

invisible(dbExecute(
  con,
  glue(
    "CREATE OR REPLACE TABLE vep_annotation AS
     SELECT
       j.id AS variant_id,
       tc.transcript_id AS tx,
       coalesce(
         list_aggregate(list_sort(list_distinct(tc.consequence_terms)), 'string_agg', '&'),
         ''
       ) AS consequence,
       coalesce(tc.impact, '') AS impact,
       'oracle'::VARCHAR AS status,
       NULL::VARCHAR AS reason
     FROM read_json({sql_q(vep_json)}, format = 'newline_delimited', sample_size = -1) j,
     UNNEST(j.transcript_consequences) u(tc)"
  )
))

run_date <- as.character(Sys.Date())
invisible(dbExecute(
  con,
  glue(
    "CREATE OR REPLACE TABLE duckvep_annotation_dump AS
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
       a.reason
     FROM vep_annotation a JOIN duckvep_sample v USING (variant_id)
     UNION ALL
     SELECT
       {sql_q(run_date)}, {sql_q(opt$corpus)}, {sql_q(opt$model_name)},
       {sql_q(oracle_version)}, {sql_q(oracle_build)}, 'duckvep',
       a.variant_id, v.chrom, v.position, v.reference, v.alternate,
       v.var_type, v.length_bin, a.tx, a.consequence, a.impact, a.status, a.reason
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
  glue("VEP annotation: {sprintf('%.3f', vep_time[['elapsed']])} s"),
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
