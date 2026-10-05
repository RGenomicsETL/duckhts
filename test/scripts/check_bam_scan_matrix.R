# Maintainer check of read_bam over a matrix of configurations: every indexed
# BAM and CRAM fixture, each of its index variants, a set of region lists,
# three projections, two DuckDB thread counts, and htslib decompression with
# and without worker threads. It needs samtools and a built
# extension, and reads only committed fixtures. Run from the repository root:
#
#   Rscript test/scripts/check_bam_scan_matrix.R [extension]
#
# Each cell is compared with an oracle that never uses the build under test:
#   - no region, and ".":  samtools view of the file, streamed without an index;
#   - named regions:       samtools view with the file's standard sidecar index;
#   - "*":                 the streamed records whose RNAME is "*".
# The index variant of a cell (without the trailing count, without statistics,
# ...) is given only to read_bam, so a defective variant cannot make the oracle
# wrong. A cell has no oracle when samtools cannot answer it: an unknown contig,
# or a CRAM whose reference is not available.
#
# The check fails when a cell differs from its oracle, when a result depends on
# the thread count, or when a region query reports a planner estimate.

arguments <- commandArgs(trailingOnly = TRUE)
extension <- normalizePath(if (length(arguments)) arguments[[1L]] else "build/release/duckhts.duckdb_extension")
data_dir <- "test/data"
stopifnot(file.exists("src/bam_reader.c"), nzchar(Sys.which("samtools")),
          requireNamespace("DBI", quietly = TRUE), requireNamespace("duckdb", quietly = TRUE),
          requireNamespace("digest", quietly = TRUE))

samtools <- function(arguments) {
  value <- suppressWarnings(system2("samtools", shQuote(arguments), stdout = TRUE, stderr = FALSE))
  if (is.null(attr(value, "status"))) value else NULL
}

# QNAME, FLAG, RNAME and POS of each record, sorted: the identity a cell compares.
record_keys <- function(lines) {
  if (is.null(lines)) return(NULL)
  if (!length(lines)) return(character())
  fields <- strsplit(lines, "\t", fixed = TRUE)
  sort(vapply(fields, function(field) paste(field[1:4], collapse = "\t"), character(1L)))
}

# One row per file and index variant; index NA means sidecar discovery.
index_variants <- function() {
  files <- list.files(data_dir, full.names = TRUE)
  alignments <- files[grepl("\\.(bam|cram)$", files)]
  do.call(rbind, lapply(alignments, function(path) {
    indexes <- files[startsWith(files, paste0(path, ".")) & grepl("\\.(bai|csi|crai)$", files)]
    standard <- indexes[!grepl("legacy|nostats", indexes)][1L]
    if (is.na(standard)) return(NULL)
    data.frame(file = path, index = c(NA_character_, indexes), standard = standard)
  }))
}

# The region lists of one file, chosen from its own header and records.
region_lists <- function(path) {
  header <- samtools(c("view", "-H", path))
  contigs <- sub("^@SQ\t.*SN:([^\t]+).*$", "\\1", header[startsWith(header, "@SQ")])
  records <- samtools(c("view", path))
  fields <- strsplit(if (is.null(records)) character() else records, "\t", fixed = TRUE)
  rname <- vapply(fields, `[[`, character(1L), 3L)
  position <- as.numeric(vapply(fields, `[[`, character(1L), 4L))
  with_records <- names(sort(table(rname[rname != "*"]), decreasing = TRUE))
  lists <- list(none = character(), whole_file = ".", unplaced = "*", unknown = "no_such_contig")
  if (length(with_records)) {
    first <- with_records[[1L]]
    low <- max(1, floor(stats::quantile(position[rname == first], 0.25)))
    high <- max(low, ceiling(stats::quantile(position[rname == first], 0.75)))
    interval <- sprintf("%s:%d-%d", first, low, high)
    lists$contig <- first
    lists$interval <- interval
    lists$overlapping <- c(interval, sprintf("%s:%d-%d", first, max(1, low - 5), high + 5))
    lists$contig_and_unplaced <- c(first, "*")
    lists$unknown_and_contig <- c("no_such_contig", first)
    if (length(with_records) > 1L) lists$two_contigs <- with_records[1:2]
  }
  without_records <- setdiff(contigs, with_records)
  if (length(without_records)) lists$empty_contig <- without_records[[1L]]
  list(lists = lists, records = records)
}

# The expected record keys of one region list, or NULL when there is no oracle.
expected_keys <- function(path, standard, region, records) {
  if (is.null(records) || "no_such_contig" %in% region) return(NULL)
  if (!length(region) || "." %in% region) return(record_keys(records))
  keys <- character()
  named <- setdiff(region, "*")
  if (length(named)) {
    # The multi-region mode splits a contig name that contains a colon, so one
    # named item uses the single-region parser, which looks the whole name up first.
    found <- samtools(c("view", if (length(named) > 1L) "-M", "-X", path, standard, named))
    if (is.null(found)) return(NULL)
    keys <- record_keys(found)
  }
  if ("*" %in% region) {
    rname <- vapply(strsplit(records, "\t", fixed = TRUE), `[[`, character(1L), 3L)
    keys <- c(keys, record_keys(records[rname == "*"]))
  }
  sort(keys)
}

con <- DBI::dbConnect(duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE,
  config = list(allow_unsigned_extensions = "true")))
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
invisible(DBI::dbExecute(con, sprintf("LOAD '%s'", extension)))
quote_sql <- function(value) as.character(DBI::dbQuoteString(con, value))
attempt <- function(sql) {
  tryCatch(DBI::dbGetQuery(con, sql), error = function(error) {
    structure(list(), failure = substr(gsub("\\s+", " ", conditionMessage(error)), 1L, 100L))
  })
}
failure_of <- function(value) if (is.null(attr(value, "failure"))) "" else attr(value, "failure")

variants <- index_variants()
plans <- list()
cells <- list()
for (row in seq_len(nrow(variants))) {
  path <- variants$file[[row]]
  index <- variants$index[[row]]
  if (is.null(plans[[path]])) plans[[path]] <- region_lists(path)
  plan <- plans[[path]]
  for (name in names(plan$lists)) {
    region <- plan$lists[[name]]
    expected <- expected_keys(path, variants$standard[[row]], region, plan$records)
    for (configuration in list(c(threads = 1L, workers = 0L), c(threads = 4L, workers = 0L),
                               c(threads = 1L, workers = 2L), c(threads = 4L, workers = 2L))) {
      threads <- configuration[["threads"]]
      source <- sprintf("read_bam(%s%s%s, decompression_threads := %d)", quote_sql(path),
        if (is.na(index)) "" else sprintf(", index_path := %s", quote_sql(index)),
        if (length(region)) sprintf(", region := %s", quote_sql(paste(region, collapse = ","))) else "",
        configuration[["workers"]])
      invisible(DBI::dbExecute(con, sprintf("SET threads=%d", threads)))
      star <- attempt(sprintf("SELECT count(*) AS n FROM %s", source))
      column <- attempt(sprintf("SELECT count(QNAME) AS n FROM %s", source))
      rows <- attempt(sprintf(
        "SELECT QNAME || chr(9) || FLAG || chr(9) || RNAME || chr(9) || POS AS key FROM %s", source))
      plan_text <- attempt(sprintf("EXPLAIN (FORMAT json) SELECT * FROM %s", source))
      failed <- any(nzchar(c(failure_of(star), failure_of(column), failure_of(rows))))
      keys <- if (failed) character() else sort(rows$key)
      estimate <- if (nzchar(failure_of(plan_text))) character() else {
        regmatches(plan_text[[2L]], regexpr('"Estimated Cardinality": "[0-9]+"', plan_text[[2L]]))
      }
      cells[[length(cells) + 1L]] <- data.frame(
        file = basename(path),
        index = if (is.na(index)) "(sidecar)" else sub(basename(path), "", basename(index), fixed = TRUE),
        region = name, threads = threads, workers = configuration[["workers"]], failed = failed,
        star = if (failed) NA_real_ else as.numeric(star$n),
        column = if (failed) NA_real_ else as.numeric(column$n),
        rows = length(keys), digest = digest::digest(keys, algo = "md5"),
        estimate = if (length(estimate)) as.numeric(gsub("[^0-9]", "", estimate)) else NA_real_,
        has_oracle = !is.null(expected),
        expected_rows = if (is.null(expected)) NA_real_ else length(expected),
        expected_digest = if (is.null(expected)) NA_character_ else digest::digest(expected, algo = "md5"),
        stringsAsFactors = FALSE)
    }
  }
}
cells <- do.call(rbind, cells)

with_oracle <- cells[cells$has_oracle, ]
wrong <- with_oracle[with_oracle$failed | with_oracle$star != with_oracle$expected_rows |
                     with_oracle$column != with_oracle$expected_rows |
                     with_oracle$digest != with_oracle$expected_digest, ]
# The four configurations of one file, index and region list must agree.
first <- cells[cells$threads == 1L & cells$workers == 0L, ]
thread_dependent <- do.call(rbind, lapply(list(c(4L, 0L), c(1L, 2L), c(4L, 2L)), function(other) {
  compared <- cells[cells$threads == other[[1L]] & cells$workers == other[[2L]], ]
  compared[compared$digest != first$digest | compared$failed != first$failed, ]
}))
# DuckDB shows 1 for a scan that reports no estimate.
region_estimates <- cells[cells$region != "none" & !is.na(cells$estimate) & cells$estimate != 1, ]

cat(sprintf("read_bam matrix: %d cells, %d files, %d file and index variants, %d region lists\n",
            nrow(cells), length(unique(cells$file)), nrow(variants), length(unique(cells$region))))
cat(sprintf("  with an oracle: %d; equal to it: %d; different: %d\n",
            nrow(with_oracle), nrow(with_oracle) - nrow(wrong), nrow(wrong)))
cat(sprintf("  without an oracle: %d (%d of them are errors)\n",
            sum(!cells$has_oracle), sum(!cells$has_oracle & cells$failed)))
cat(sprintf("  results that depend on DuckDB threads or decompression workers: %d\n", nrow(thread_dependent)))
cat(sprintf("  region queries that report an estimate: %d\n", nrow(region_estimates)))
if (nrow(wrong)) {
  print(unique(wrong[c("file", "index", "region", "workers", "failed", "rows", "expected_rows")]), row.names = FALSE)
}
if (nrow(wrong) || nrow(thread_dependent) || nrow(region_estimates)) {
  stop("read_bam differs from the oracle in the configuration matrix", call. = FALSE)
}
cat("read_bam matrix: OK\n")
