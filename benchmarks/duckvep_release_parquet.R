#!/usr/bin/env Rscript

# Materialize an Ensembl variation consequence VCF as typed Parquet and record
# both the complete reader projection and the narrow release-product view.

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
op <- add_option(op, "--input", default = "", help = "local Ensembl VCF/BCF")
op <- add_option(op, "--source-url", dest = "source_url", default = "")
op <- add_option(op, "--output-dir", dest = "output_dir", default = "")
op <- add_option(op, "--assembly", default = "GRCh38")
op <- add_option(op, "--release", type = "integer", default = 116L)
op <- add_option(op, "--chromosome", default = "")
op <- add_option(op, "--threads", type = "integer", default = 4L)
op <- add_option(
  op,
  "--row-group-size",
  dest = "row_group_size",
  type = "integer",
  default = 100000L
)
op <- add_option(op, "--compression", default = "zstd")
op <- add_option(
  op,
  "--extension",
  default = file.path(root, "build", "release", "duckhts.duckdb_extension")
)
op <- add_option(
  op,
  "--history",
  default = file.path(
    root,
    "benchmarks",
    "data",
    "duckvep_release_parquet.csv"
  )
)
op <- add_option(
  op,
  "--source-revision",
  dest = "source_revision",
  default = ""
)
op <- add_option(
  op,
  "--overwrite",
  action = "store_true",
  default = FALSE
)
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
if (!nzchar(opt$input)) die("--input is required")
if (!nzchar(opt$output_dir)) die("--output-dir is required")
if (!file.exists(opt$input)) die("input does not exist: {opt$input}")
if (!file.exists(opt$extension)) die("extension does not exist: {opt$extension}")
if (opt$release < 1L) die("--release must be positive")
if (opt$threads < 1L) die("--threads must be positive")
if (opt$row_group_size < 2048L) die("--row-group-size must be at least 2048")
if (!opt$compression %in% c("zstd", "snappy", "gzip", "lz4_raw", "uncompressed")) {
  die("unsupported Parquet compression: {opt$compression}")
}

input <- normalizePath(opt$input, mustWork = TRUE)
extension <- normalizePath(opt$extension, mustWork = TRUE)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(opt$output_dir, mustWork = TRUE)

sha256 <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) die("sha256sum failed for {path}")
  sub("[[:space:]].*$", "", output[[1L]])
}

if (!nzchar(opt$source_revision)) {
  revision <- system2(
    "git",
    c("-C", root, "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  )
  status <- attr(revision, "status")
  if (!is.null(status) && status != 0L) die("cannot determine source revision")
  opt$source_revision <- trimws(revision[[1L]])
}

drv <- duckdb(config = list(allow_unsigned_extensions = "true"))
con <- dbConnect(drv)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
sql_q <- function(x) as.character(dbQuoteString(con, x))
invisible(dbExecute(con, glue("SET threads = {opt$threads}")))
invisible(dbExecute(con, glue("LOAD {sql_q(extension)}")))

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
    "({declared_version}); run `make release` first"
  )
}

input_sql <- sql_q(input)
projections <- list(
  full_typed = glue(
    "SELECT * FROM read_bcf_v2(
       {input_sql}, include_format := false, scan_mode := 'sequential'
     )"
  ),
  consequence = glue(
    "SELECT * FROM read_bcf_v2(
       {input_sql},
       info_fields := 'VE',
       include_format := false,
       vep_fields := 'Allele,Consequence,Feature_type,Feature,Amino_acids,SIFT',
       scan_mode := 'sequential'
     )"
  )
)

stem <- sub("[.]vcf([.]bgz|[.]gz)?$", "", basename(input), ignore.case = TRUE)
stem <- gsub("[^A-Za-z0-9._-]", "_", stem)
artifacts <- file.path(
  output_dir,
  paste0(stem, ".", names(projections), ".parquet")
)
names(artifacts) <- names(projections)
existing <- artifacts[file.exists(artifacts)]
if (length(existing) != 0L && !isTRUE(opt$overwrite)) {
  die("output exists; pass --overwrite:\n{paste(existing, collapse = '\n')}")
}
if (isTRUE(opt$overwrite)) unlink(existing)

source_bytes <- unname(file.info(input)$size)
source_hash <- sha256(input)
source_name <- if (nzchar(opt$source_url)) opt$source_url else input
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

rows <- vector("list", length(projections))
for (i in seq_along(projections)) {
  projection <- names(projections)[[i]]
  artifact <- artifacts[[projection]]
  timing <- system.time({
    invisible(dbExecute(
      con,
      glue(
        "COPY ({projections[[projection]]}) TO {sql_q(artifact)}
         (FORMAT PARQUET, COMPRESSION {opt$compression},
          ROW_GROUP_SIZE {opt$row_group_size})"
      )
    ))
  })
  counts <- dbGetQuery(
    con,
    glue(
      "SELECT
         count(*)::DOUBLE AS records,
         coalesce(sum(list_count(ALT)), 0)::DOUBLE AS alternate_alleles,
         coalesce(sum(list_count(VEP_Consequence)), 0)::DOUBLE AS csq_entries
       FROM read_parquet({sql_q(artifact)})"
    )
  )
  column_count <- nrow(dbGetQuery(
    con,
    glue("DESCRIBE SELECT * FROM read_parquet({sql_q(artifact)})")
  ))
  parquet_bytes <- unname(file.info(artifact)$size)
  elapsed <- unname(timing[["elapsed"]])
  rows[[i]] <- data.frame(
    run_date = as.character(Sys.Date()),
    source_revision = opt$source_revision,
    declared_extension_version = declared_version,
    loaded_extension_version = loaded_version,
    duckdb_version = duckdb_version,
    host = unname(Sys.info()[["nodename"]]),
    cpu = cpu,
    assembly = opt$assembly,
    ensembl_release = opt$release,
    chromosome = opt$chromosome,
    source = source_name,
    source_sha256 = source_hash,
    source_bytes = source_bytes,
    projection = projection,
    columns = column_count,
    records = counts$records[[1L]],
    alternate_alleles = counts$alternate_alleles[[1L]],
    csq_entries = counts$csq_entries[[1L]],
    threads = opt$threads,
    compression = opt$compression,
    row_group_size = opt$row_group_size,
    elapsed_seconds = elapsed,
    records_per_second = counts$records[[1L]] / elapsed,
    artifact = basename(artifact),
    artifact_sha256 = sha256(artifact),
    parquet_bytes = parquet_bytes,
    parquet_bytes_per_record = parquet_bytes / counts$records[[1L]],
    parquet_fraction_of_source = parquet_bytes / source_bytes,
    source_to_parquet_ratio = source_bytes / parquet_bytes,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
rows <- do.call(rbind, rows)

if (nzchar(opt$history)) {
  dir.create(dirname(opt$history), recursive = TRUE, showWarnings = FALSE)
  history <- rows
  if (file.exists(opt$history)) {
    old <- utils::read.csv(
      opt$history,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (!identical(names(old), names(rows))) {
      die("history schema does not match: {opt$history}")
    }
    same <- old$source_revision == opt$source_revision &
      old$source_sha256 == source_hash &
      old$projection %in% rows$projection &
      old$compression == opt$compression &
      old$row_group_size == opt$row_group_size
    history <- rbind(old[!same, , drop = FALSE], rows)
  }
  history <- history[
    order(history$run_date, history$source_revision, history$source, history$projection),
    ,
    drop = FALSE
  ]
  temporary <- tempfile("duckvep-release-parquet-", dirname(opt$history), ".csv")
  utils::write.csv(history, temporary, row.names = FALSE)
  if (!file.rename(temporary, opt$history)) {
    unlink(temporary)
    die("cannot replace history: {opt$history}")
  }
}

for (i in seq_len(nrow(rows))) {
  cat(
    glue(
      "{rows$projection[[i]]}: ",
      "{format(rows$records[[i]], big.mark = ',', scientific = FALSE)} records, ",
      "{format(rows$parquet_bytes[[i]], big.mark = ',', scientific = FALSE)} bytes, ",
      "{sprintf('%.3f', rows$parquet_fraction_of_source[[i]])}x source size, ",
      "{sprintf('%.3f', rows$elapsed_seconds[[i]])} s"
    ),
    "\n",
    sep = ""
  )
}
cat(glue("source sha256: {source_hash}"), "\n", sep = "")
if (nzchar(opt$history)) cat(glue("history: {opt$history}"), "\n", sep = "")
