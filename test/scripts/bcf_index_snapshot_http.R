# Localhost-only transport conformance; no public network or package installation.
# Rscript test/scripts/bcf_index_snapshot_http.R build/release/duckhts.duckdb_extension
# Uses the benchmark environment's DBI/duckdb/callr plus httpuv's real Range server.
test_bcf_index_snapshot_http <- function(extension) {
  sys.source("r/Rduckhts/inst/tinytest/bcf_snapshot_query.R", envir = environment())
  stopifnot(requireNamespace("httpuv", quietly = TRUE),
            requireNamespace("callr", quietly = TRUE))
  extension <- normalizePath(extension, mustWork = TRUE)
  fixtures <- normalizePath("test/data", mustWork = TRUE)
  directory <- tempfile("duckhts-bcf-http-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  served <- file.path(directory, "served")
  cache <- file.path(directory, "cache")
  dir.create(served)
  dir.create(cache)
  port <- httpuv::randomPort()
  server <- callr::r_bg(function(directory, port) {
    handle <- httpuv::startServer("127.0.0.1", port, list(
      call = function(request) list(status = 404L, headers = list(), body = "missing"),
      staticPaths = list("/" = httpuv::staticPath(directory,
                        headers = list("Cache-Control" = "no-store")))))
    on.exit(httpuv::stopServer(handle))
    cat("ready\n")
    flush.console()
    repeat httpuv::service(100)
  }, args = list(directory = served, port = port), stdout = "|", stderr = "|")
  on.exit(server$kill(), add = TRUE, after = FALSE)
  ready <- FALSE
  for (attempt in seq_len(200)) {
    stopifnot(server$is_alive())
    if ("ready" %in% server$read_output_lines()) { ready <- TRUE; break }
    Sys.sleep(0.025)
  }
  stopifnot(ready)
  previous <- setwd(cache) # HTSlib's remote-index cache never touches fixtures.
  on.exit(setwd(previous), add = TRUE, after = FALSE)
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE, after = FALSE)
  quote <- function(value) as.character(DBI::dbQuoteString(con, value))
  DBI::dbExecute(con, paste("LOAD", quote(extension)))
  DBI::dbExecute(con, "SET threads=4")
  fixture <- function(suffix) file.path(fixtures, paste0("bcf_scan_contigs.", suffix))
  bind_checks <- 0L
  control_checks <- 0L
  result_checks <- 0L
  for (kind in c("bcf", "tbi", "csi")) for (explicit in c(FALSE, TRUE)) {
    suffix <- if (kind == "bcf") "bcf" else "full.vcf.gz"
    source <- fixture(suffix)
    local_path <- file.path(served, paste0("input.", suffix))
    local_index <- if (explicit) file.path(served, "private.index") else
      paste0(local_path, if (kind == "tbi") ".tbi" else ".csi")
    original_index <- paste0(source, if (kind == "bcf") ".csi" else paste0(".index.", kind))
    stopifnot(file.copy(source, local_path, overwrite = TRUE))
    url <- paste0("http://127.0.0.1:", port, "/", basename(local_path))
    index_url <- paste0("http://127.0.0.1:", port, "/", basename(local_index))
    projection <- "SELECT POS, * EXCLUDE (POS) FROM "
    expected <- DBI::dbGetQuery(con, paste0(projection, "read_bcf(", quote(source),
      ",tidy_format:=true,scan_mode:='sequential') ORDER BY ALL"))
    stopifnot(nrow(expected) == 10L)
    for (mode in c("auto", "region", "sequential")) {
      stopifnot(file.copy(original_index, local_index, overwrite = TRUE))
      unlink(list.files(cache, full.names = TRUE))
      reader <- paste0("read_bcf(", quote(url), ",tidy_format:=true,decompression_threads:=0",
        if (explicit) paste0(",index_path:=", quote(index_url)) else "",
        switch(mode, auto = "", region = ",region:='chr3:20-40,chr3:20-20'",
               sequential = ",scan_mode:='sequential'"), ")")
      oracle <- if (mode == "region") expected[expected$CHROM == "chr3", ] else expected
      rownames(oracle) <- NULL
      # Populate the remote-index cache so both the served index and every
      # downloaded copy participate in the after-bind rewrite.
      stopifnot(identical(DBI::dbGetQuery(con, paste0(projection, reader, " ORDER BY ALL")), oracle))
      index_paths <- c(local_index, list.files(cache, full.names = TRUE))
      min_shift <- if (kind == "tbi") 0L else 14L
      for (replacement in c("empty", "shifted", "identical")) {
        replacement_path <- if (replacement == "identical") local_path else
          fixture(paste0(replacement, ".", if (kind == "bcf") "bcf" else "vcf.gz"))
        for (count in c(FALSE, TRUE)) {
          statements <- .bcf_snapshot_statements(con, reader, local_path, index_paths,
                                                 replacement_path, min_shift, count)
          DBI::dbGetQuery(con, statements$restore)
          observed <- DBI::dbGetQuery(con, statements$after_bind)
          if (count) stopifnot(observed$n[[1]] == nrow(oracle))
          else stopifnot(identical(observed, oracle))
        }
        bind_checks <- bind_checks + 1L
      }
      shifted_path <- fixture(paste0("shifted.", if (kind == "bcf") "bcf" else "vcf.gz"))
      statements <- .bcf_snapshot_statements(con, reader, local_path, index_paths,
                                             shifted_path, min_shift)
      # The window exists: a scan on the probe side reads the replaced index.
      DBI::dbGetQuery(con, statements$replace)
      shifted_digest <- unname(tools::md5sum(local_index))
      DBI::dbGetQuery(con, statements$restore)
      stopifnot(!identical(unname(tools::md5sum(local_index)), shifted_digest))
      stopifnot(identical(DBI::dbGetQuery(con, statements$window)$digest, shifted_digest))
      # Indexed reads must fail when the shifted offsets are in place before the
      # reader binds: invalid BCF block headers or unsatisfiable VCF HTTP ranges.
      DBI::dbGetQuery(con, statements$restore)
      DBI::dbGetQuery(con, statements$replace)
      sql <- paste0(projection, reader, " ORDER BY ALL")
      if (mode == "sequential") {
        stopifnot(identical(DBI::dbGetQuery(con, sql), oracle))
      } else {
        error <- tryCatch(DBI::dbGetQuery(con, sql), error = identity)
        stopifnot(inherits(error, "error"))
        stopifnot(grepl("read_bcf: failed to read or parse BCF/VCF record",
                        conditionMessage(error), fixed = TRUE))
      }
      control_checks <- control_checks + 1L
      # Preserve completed-result lifetime checks separately from scan checks.
      fetch_after_index_change <- function(query, damage) {
        stopifnot(file.copy(original_index, local_index, overwrite = TRUE))
        unlink(list.files(cache, full.names = TRUE))
        result <- DBI::dbSendQuery(con, query)
        on.exit(DBI::dbClearResult(result))
        if (damage == "missing") {
          unlink(local_index)
        } else if (damage == "corrupt") {
          writeLines("broken index", local_index)
        } else {
          replacement <- if (damage == "identical") original_index else
            fixture(paste0(damage, ".", if (kind == "bcf") "bcf.csi" else paste0("vcf.gz.index.", kind)))
          stopifnot(file.copy(replacement, local_index, overwrite = TRUE))
        }
        # A downloaded copy must not hide a reader that reloads the changed index.
        unlink(list.files(cache, full.names = TRUE))
        DBI::dbFetch(result)
      }
      for (damage in c("missing", "corrupt", "empty", "shifted", "identical")) {
        rows <- fetch_after_index_change(paste0(projection, reader, " ORDER BY ALL"), damage)
        count <- fetch_after_index_change(paste("SELECT count(*) AS n FROM", reader), damage)
        stopifnot(identical(rows, oracle), count$n[[1]] == nrow(oracle))
        result_checks <- result_checks + 1L
      }
    }
    unlink(c(local_path, local_index))
  }
  cat("BCF HTTP:", bind_checks, "after-bind rewrite cases,", control_checks,
      "before-bind controls,", result_checks,
      "result-lifetime cases; exact rows/counts, 4 requested workers: OK\n")
}

args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
test_bcf_index_snapshot_http(args[[1]])
