#!/usr/bin/env Rscript

# Localhost Range-I/O parity for panel-site BAM/CRAM extraction.

fail <- function(...) stop(..., call. = FALSE)

main <- function(args) {
  if (length(args) != 1L) {
    fail("usage: Rscript somalier_extraction_http.R <duckhts-extension>")
  }
  for (package in c("DBI", "duckdb", "callr", "httpuv")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      fail("required R package is unavailable: ", package)
    }
  }
  extension <- normalizePath(args[[1L]], mustWork = TRUE)
  fixtures <- normalizePath("test/data", mustWork = TRUE)
  directory <- tempfile("duckhts-somalier-http-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  served <- file.path(directory, "served")
  cache <- file.path(directory, "cache")
  dir.create(served)
  dir.create(cache)
  files <- c(
    "range.bam", "range.bam.bai", "range.cram", "range.cram.crai",
    "ce.fa", "ce.fa.fai"
  )
  copied <- file.copy(file.path(fixtures, files), file.path(served, files))
  if (!all(copied)) fail("could not stage localhost extraction fixtures")

  port <- httpuv::randomPort()
  server <- callr::r_bg(function(directory, port) {
    handle <- httpuv::startServer("127.0.0.1", port, list(
      call = function(request) {
        list(status = 404L, headers = list(), body = "missing")
      },
      staticPaths = list("/" = httpuv::staticPath(
        directory, headers = list("Cache-Control" = "no-store")
      ))
    ))
    on.exit(httpuv::stopServer(handle))
    cat("ready\n")
    flush.console()
    repeat httpuv::service(100)
  }, args = list(directory = served, port = port), stdout = "|", stderr = "|")
  on.exit(server$kill(), add = TRUE, after = FALSE)
  ready <- FALSE
  for (attempt in seq_len(200L)) {
    if (!server$is_alive()) fail("localhost Range server stopped")
    if ("ready" %in% server$read_output_lines()) {
      ready <- TRUE
      break
    }
    Sys.sleep(0.025)
  }
  if (!ready) fail("localhost Range server did not start")

  previous <- setwd(cache)
  on.exit(setwd(previous), add = TRUE, after = FALSE)
  con <- DBI::dbConnect(duckdb::duckdb(
    config = list(allow_unsigned_extensions = "true")
  ))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE, after = FALSE)
  quote <- function(value) as.character(DBI::dbQuoteString(con, value))
  DBI::dbExecute(con, paste("LOAD", quote(extension)))
  DBI::dbExecute(con, "SET threads = 4")
  DBI::dbExecute(con, paste(
    "CREATE TABLE panel AS SELECT * FROM (VALUES",
    "('WBcel235', 0::UBIGINT, 'CHROMOSOME_I', 1::UBIGINT, 'A', 'G'),",
    "('WBcel235', 1::UBIGINT, 'CHROMOSOME_I', 914::UBIGINT, 'A', 'C'),",
    "('WBcel235', 2::UBIGINT, 'CHROMOSOME_I', 2::UBIGINT, 'A', 'G'),",
    "('WBcel235', 3::UBIGINT, 'NOT_IN_REFERENCE', 1::UBIGINT, 'A', 'C'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))
  panel_file <- tempfile(fileext = ".parquet")
  on.exit(unlink(panel_file), add = TRUE, after = FALSE)
  DBI::dbExecute(con, paste("COPY panel TO", quote(panel_file), "(FORMAT parquet)"))

  extract <- function(
    source, reference, index = NULL, reference_index = NULL, worker_count = 1L
  ) {
    arguments <- c(
      quote(source), "NULL", "'sample-1'", quote(reference),
      paste0("panel_parquet := ", quote(panel_file)),
      if (!is.null(index)) paste0("index_path := ", quote(index)),
      if (!is.null(reference_index)) {
        paste0("reference_index_path := ", quote(reference_index))
      },
      paste0("worker_count := ", worker_count),
      "remote_block_bytes := 65536", "remote_cache_bytes := 1048576",
      "reference_cache_bytes := 1048576"
    )
    DBI::dbGetQuery(con, paste0(
      "SELECT * EXCLUDE(source_path) FROM duckhts_somalier_bam_counts(",
      paste(arguments, collapse = ", "), ") ORDER BY site_index"
    ))
  }
  url <- function(name) paste0("http://127.0.0.1:", port, "/", name)
  reference <- file.path(fixtures, "ce.fa")
  for (format in c("bam", "cram")) {
    source_name <- paste0("range.", format)
    index_name <- paste0(source_name, if (format == "bam") ".bai" else ".crai")
    expected <- extract(
      file.path(fixtures, source_name), reference,
      file.path(fixtures, index_name), file.path(fixtures, "ce.fa.fai")
    )
    local_workers <- extract(
      file.path(fixtures, source_name), reference,
      file.path(fixtures, index_name), file.path(fixtures, "ce.fa.fai"),
      worker_count = 4L
    )
    inferred <- extract(url(source_name), url("ce.fa"), worker_count = 4L)
    explicit <- extract(
      url(source_name), url("ce.fa"), url(index_name), url("ce.fa.fai"),
      worker_count = 4L
    )
    if (!identical(expected, local_workers) ||
        !identical(expected, inferred) || !identical(expected, explicit)) {
      fail(format, " local/remote count rows differ")
    }
  }
  cat(
    "Somalier BAM/CRAM localhost Range extraction: OK ",
    "(1/4 workers; auto and explicit indexes; complete typed rows)\n",
    sep = ""
  )
}

main(commandArgs(trailingOnly = TRUE))
