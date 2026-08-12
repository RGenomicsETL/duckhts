#!/usr/bin/env Rscript
# Benchmark DuckHTS cgranges scalar probes against bedtk and bedtools.
# This is an R-owned orchestration layer; DuckDB executes the SQL workload.

if (!requireNamespace("optparse", quietly = TRUE)) stop("optparse is required", call. = FALSE)
if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("duckdb", quietly = TRUE)) {
  stop("DBI and duckdb are required", call. = FALSE)
}

option_list <- list(
  optparse::make_option("--extension", type = "character", default = "build/release/duckhts.duckdb_extension"),
  optparse::make_option("--bedtk", type = "character", default = ".sync/bedtk/bedtk"),
  optparse::make_option("--bedtools", type = "character", default = "bedtools"),
  optparse::make_option("--subjects", type = "integer", default = 50000L),
  optparse::make_option("--queries", type = "integer", default = 5000L),
  optparse::make_option("--seed", type = "integer", default = 1L),
  optparse::make_option("--passes", type = "integer", default = 3L),
  optparse::make_option("--subject-bed", type = "character"),
  optparse::make_option("--query-bed", type = "character"),
  optparse::make_option("--limit-subjects", type = "integer"),
  optparse::make_option("--limit-queries", type = "integer"),
  optparse::make_option("--label", type = "character"),
  optparse::make_option("--out-dir", type = "character", default = file.path(path.expand(Sys.getenv("DUCKHTS_CACHE_DIR", unset = "~/.cache/duckhts")), "benchmarks", "cgranges")),
  optparse::make_option("--timeout", type = "integer", default = 3600L)
)
options <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
if (options$subjects < 1L || options$queries < 1L || options$passes < 1L) {
  stop("subjects, queries, and passes must be positive", call. = FALSE)
}

sql_string <- function(value) gsub("'", "''", value, fixed = TRUE)
command_path <- function(value) {
  if (file.exists(value) && file.access(value, 1L) == 0L) return(normalizePath(value))
  found <- Sys.which(value)
  if (nzchar(found)) normalizePath(found) else ""
}
count_bed_rows <- function(path) {
  sum(nzchar(lines <- readLines(path, warn = FALSE)) & !startsWith(lines, "#"))
}
copy_bed_prefix <- function(source, destination, limit) {
  lines <- readLines(source, warn = FALSE)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (!is.null(limit) && !is.na(limit)) lines <- head(lines, limit)
  writeLines(lines, destination, useBytes = TRUE)
  length(lines)
}
write_synthetic_bed <- function(path, n, seed, query = FALSE) {
  set.seed(seed)
  chroms <- paste0("chr", 1:4)
  index <- seq_len(n) - 1L
  width <- if (query) 30L + index %% 120L else 50L + index %% 200L
  starts <- floor(stats::runif(n, 0, 5000000L - width))
  chrom <- chroms[if (query) (index * 3L) %% length(chroms) + 1L else index %% length(chroms) + 1L]
  write.table(data.frame(chrom, starts, starts + width, if (query) paste0("q", index) else paste0("s", index)),
    path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
prepare_inputs <- function() {
  if (nzchar(options$`subject-bed` %||% "") && nzchar(options$`query-bed` %||% "")) {
    subject <- file.path(options$`out-dir`, "subject.bed")
    query <- file.path(options$`out-dir`, "query.bed")
    subject_n <- copy_bed_prefix(options$`subject-bed`, subject, options$`limit-subjects`)
    query_n <- copy_bed_prefix(options$`query-bed`, query, options$`limit-queries`)
    return(list(subject = subject, query = query, subject_n = subject_n, query_n = query_n,
      label = options$label %||% "external_bed"))
  }
  subject <- file.path(options$`out-dir`, "subject.bed")
  query <- file.path(options$`out-dir`, "query.bed")
  write_synthetic_bed(subject, options$subjects, options$seed, FALSE)
  write_synthetic_bed(query, options$queries, options$seed + 1L, TRUE)
  list(subject = subject, query = query, subject_n = options$subjects, query_n = options$queries,
    label = options$label %||% "synthetic")
}
`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) y else x
measure_runs <- function(fun) {
  values <- numeric(options$passes)
  result <- NULL
  for (i in seq_len(options$passes)) {
    started <- proc.time()[["elapsed"]]
    result <- fun()
    values[[i]] <- proc.time()[["elapsed"]] - started
  }
  list(result = result, seconds = values)
}
run_command <- function(command, args) {
  output <- system2(command, args, stdout = TRUE, stderr = TRUE, timeout = options$timeout)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop(paste(output, collapse = "\n"), call. = FALSE)
  output
}

out_dir <- normalizePath(options$`out-dir`, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
inputs <- prepare_inputs()
extension <- normalizePath(options$extension, mustWork = TRUE)
con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")), dbdir = ":memory:")
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
DBI::dbExecute(con, sprintf("LOAD '%s'", sql_string(extension)))
subject_sql <- sql_string(normalizePath(inputs$subject))
query_sql <- sql_string(normalizePath(inputs$query))
index_name <- "bench_idx"
build_started <- proc.time()[["elapsed"]]
DBI::dbGetQuery(con, sprintf(
  "SELECT duckhts_cgranges_from_query('%s', 'SELECT chrom, start, \"end\" FROM read_bed(''%s'')', 'chrom', 'start', 'end')",
  index_name, subject_sql
))
DBI::dbGetQuery(con, sprintf("SELECT duckhts_cgranges_index('%s')", index_name))
build_seconds <- proc.time()[["elapsed"]] - build_started
on.exit(DBI::dbGetQuery(con, sprintf("SELECT duckhts_cgranges_destroy('%s')", index_name)), add = TRUE)

rows <- list()
add_row <- function(tool, variant, measured, hits = NA_real_) {
  rows[[length(rows) + 1L]] <<- data.frame(
    tool = tool, variant = variant, subject_intervals = inputs$subject_n,
    query_intervals = inputs$query_n, passes = options$passes,
    build_index_sec = if (tool == "duckhts") build_seconds else 0,
    query_total_sec = sum(measured$seconds), query_pass_1_sec = measured$seconds[[1L]],
    total_elapsed_sec = if (tool == "duckhts") build_seconds + sum(measured$seconds) else sum(measured$seconds),
    peak_rss_mb = NA_real_, matched_query_intervals = measured$result[[1L]], total_hits = hits,
    stringsAsFactors = FALSE
  )
}
filter_query <- sprintf("SELECT count(*)::BIGINT AS n FROM read_bed('%s') q WHERE duckhts_cgranges_has_overlap('%s', q.chrom, q.start::BIGINT, q.\"end\"::BIGINT)", query_sql, index_name)
count_query <- sprintf("WITH counts AS (SELECT duckhts_cgranges_count_overlaps('%s', q.chrom, q.start::BIGINT, q.\"end\"::BIGINT) AS n FROM read_bed('%s') q) SELECT count(*) FILTER (WHERE n > 0)::BIGINT AS matched, COALESCE(sum(n), 0)::BIGINT AS hits FROM counts", index_name, query_sql)
expand_query <- sprintf("WITH q AS (SELECT row_number() OVER () AS qid, chrom, start, \"end\" FROM read_bed('%s')), hits AS (SELECT qid FROM q CROSS JOIN UNNEST(duckhts_cgranges_overlaps_list('%s', chrom, start::BIGINT, \"end\"::BIGINT)) AS u(hit)) SELECT count(DISTINCT qid)::BIGINT AS matched, count(*)::BIGINT AS hits FROM hits", query_sql, index_name)
filter <- measure_runs(function() DBI::dbGetQuery(con, filter_query)[[1L]])
add_row("duckhts", "scalar_filter", filter)
count <- measure_runs(function() unlist(DBI::dbGetQuery(con, count_query)[1, ], use.names = FALSE))
add_row("duckhts", "scalar_count", count, count$result[[2L]])
expand <- measure_runs(function() unlist(DBI::dbGetQuery(con, expand_query)[1, ], use.names = FALSE))
add_row("duckhts", "scalar_expand", expand, expand$result[[2L]])

bedtk <- command_path(options$bedtk)
if (nzchar(bedtk)) {
  measured <- measure_runs(function() sum(nzchar(run_command(bedtk, c("flt", inputs$subject, inputs$query)))))
  add_row("bedtk", "flt", measured)
}
bedtools <- command_path(options$bedtools)
if (nzchar(bedtools)) {
  for (mode in c("u", "c", "wa_wb")) {
    measured <- measure_runs(function() {
      args <- c("intersect", "-a", inputs$query, "-b", inputs$subject, if (mode == "u") "-u" else if (mode == "c") "-c" else c("-wa", "-wb"))
      output <- run_command(bedtools, args)
      if (mode == "u") return(c(sum(nzchar(output)), NA_real_))
      if (mode == "c") {
        counts <- as.numeric(sub(".*\t", "", output))
        return(c(sum(counts > 0), sum(counts)))
      }
      fields <- strsplit(output[nzchar(output)], "\t", fixed = TRUE)
      qfields <- length(strsplit(readLines(inputs$query, n = 1L), "\t", fixed = TRUE)[[1L]])
      c(length(unique(vapply(fields, function(x) paste(x[seq_len(qfields)], collapse = "\t"), character(1L)))), length(fields))
    })
    add_row("bedtools", paste0("intersect_", mode), measured, measured$result[[2L]])
  }
}
results <- do.call(rbind, rows)
filter_rows <- results[results$variant %in% c("scalar_filter", "flt", "intersect_u"), ]
if (nrow(filter_rows) > 1L && length(unique(filter_rows$matched_query_intervals)) != 1L) stop("overlap-existence mismatch", call. = FALSE)
count_rows <- results[results$variant %in% c("scalar_count", "scalar_expand", "intersect_c", "intersect_wa_wb"), ]
if (nrow(count_rows) > 1L && length(unique(count_rows$total_hits)) != 1L) stop("overlap-count mismatch", call. = FALSE)
utils::write.csv(results, file.path(out_dir, "summary.csv"), row.names = FALSE, na = "")
metadata <- list(
  source_revision = trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
  extension = extension, extension_md5 = unname(tools::md5sum(extension)),
  dataset = inputs$label, subject_bed = normalizePath(inputs$subject), query_bed = normalizePath(inputs$query),
  subject_intervals = inputs$subject_n, query_intervals = inputs$query_n, passes = options$passes,
  bedtk_available = nzchar(bedtk), bedtools_available = nzchar(bedtools)
)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required", call. = FALSE)
jsonlite::write_json(metadata, file.path(out_dir, "metadata.json"), pretty = TRUE, auto_unbox = TRUE)
print(results, row.names = FALSE)
