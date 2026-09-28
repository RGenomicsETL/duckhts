#' Stage the keyed bigsnpr ancestry reference as Parquet.
#'
#' The two checksum-validated CSV products are the upstream authority. The
#' derived product stores 16 DOUBLE loadings and 21 DOUBLE frequencies per
#' (chromosome, position, allele_a, allele_b), ordered by that key. A receipt
#' binds its checksum to both source checksums, the registry derivation and the
#' DuckDB writer version, and certifies one normalized locus per row with
#' complete finite numeric fields and group frequencies in [0, 1].
#'
#' @return The cached Parquet path.
#' @export
duckhts_bench_stage_ancestry_parquet <- function() {
  source_ids <- c("ancestry_ref_freqs", "ancestry_projection")
  registry <- duckhts_bench_registry()
  source_rows <- registry[match(source_ids, registry$id), , drop = FALSE]
  if (anyNA(source_rows$id)) stop("missing ancestry CSV registry entry", call. = FALSE)
  source_hashes <- sub(".*(?:^|;)sha256=([[:xdigit:]]{64}).*", "\\1",
                       source_rows$supplier_identity, perl = TRUE)
  if (any(nchar(source_hashes) != 64L)) {
    stop("ancestry CSV registry requires pinned SHA-256", call. = FALSE)
  }
  product <- registry[match("ancestry_reference_parquet", registry$id), , drop = FALSE]
  if (nrow(product) != 1L || is.na(product$transform) ||
      !nzchar(product$transform)) {
    stop("missing ancestry Parquet derivation", call. = FALSE)
  }
  derivation <- digest::digest(product$transform, algo = "sha256", serialize = FALSE)
  output <- duckhts_bench_artifact_path("ancestry_reference_parquet")
  receipt <- paste0(output, ".sources.tsv")
  version <- as.character(utils::packageVersion("duckdb"))
  expected <- data.frame(id = c(source_ids, "duckdb", "derivation", "validation", "parquet"),
                         sha256 = c(source_hashes, version, derivation,
                                    "frequency_unit_unique_v2", ""))
  if (file.exists(output) && file.exists(receipt)) {
    stored <- tryCatch(utils::read.delim(receipt, colClasses = "character"),
                       error = function(e) NULL)
    if (!is.null(stored) && identical(stored$id, expected$id) &&
        identical(stored$sha256[1:5], expected$sha256[1:5]) &&
        identical(unname(digest::digest(file = output, algo = "sha256")),
                  stored$sha256[[6L]])) {
      return(output)
    }
    stop("cached ancestry Parquet identity does not match its sources: ", output,
         call. = FALSE)
  }
  if (file.exists(output) || file.exists(receipt)) {
    stop("incomplete ancestry Parquet staging: ", output, call. = FALSE)
  }
  sources <- vapply(source_ids, duckhts_bench_fetch, character(1L))
  driver_args <- list(dbdir = ":memory:")
  if ("shared_home" %in% names(formals(duckdb::duckdb))) {
    driver_args$shared_home <- FALSE
  }
  con <- DBI::dbConnect(do.call(duckdb::duckdb, driver_args))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  headers <- lapply(sources, function(path) {
    stream <- gzfile(path)
    on.exit(close(stream))
    line <- readLines(stream, n = 1L, warn = FALSE)
    names(utils::read.csv(text = line, nrows = 0L, check.names = FALSE))
  })
  keys <- c("chr", "pos", "a0", "a1")
  groups <- setdiff(headers[[1L]], c(keys, "rsid"))
  pcs <- setdiff(headers[[2L]], c(keys, "rsid"))
  if (length(groups) != 21L || length(pcs) != 16L ||
      !identical(pcs, paste0("PC", seq_len(16L))) ||
      anyDuplicated(groups) || !all(keys %in% headers[[1L]]) ||
      !all(keys %in% headers[[2L]])) {
    stop("unexpected ancestry reference product schema", call. = FALSE)
  }
  quote_id <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  quote_str <- function(x) as.character(DBI::dbQuoteString(con, x))
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(output, ".partial-", Sys.getpid(), ".parquet")
  intermediate <- paste0(output, ".columns-", Sys.getpid(), c("-ref.parquet", "-pc.parquet"))
  spill <- paste0(output, ".spill-", Sys.getpid())
  temporary_receipt <- paste0(receipt, ".partial-", Sys.getpid())
  on.exit(unlink(c(temporary, temporary_receipt, intermediate, spill),
                 recursive = TRUE), add = TRUE)
  DBI::dbExecute(con, "SET threads=4")
  DBI::dbExecute(con, "SET memory_limit='1800MB'")
  DBI::dbExecute(con, "SET preserve_insertion_order=false")
  DBI::dbExecute(con, paste0("SET temp_directory=", quote_str(spill)))
  csv_sources <- paste0("read_csv_auto(", quote_str(sources),
                        ", types={'a0':'VARCHAR','a1':'VARCHAR'})")
  for (index in seq_along(sources)) {
    columns <- if (index == 1L) groups else pcs
    selected <- c("chr::UTINYINT AS chromosome", "pos::UINTEGER AS position",
                  "a0 AS allele_a", "a1 AS allele_b",
                  paste0(quote_id(columns), "::DOUBLE AS ", quote_id(columns)))
    DBI::dbExecute(con, paste0(
      "COPY (SELECT ", paste(selected, collapse = ", "), " FROM ", csv_sources[[index]],
      ") TO ", quote_str(intermediate[[index]]),
      " (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 32768)"))
  }
  row_counts <- vapply(intermediate, function(path) DBI::dbGetQuery(con, paste0(
    "SELECT count(*) AS n FROM read_parquet(", quote_str(path), ")"))$n, numeric(1L))
  selected <- c("p.chromosome", "p.position", "p.allele_a", "p.allele_b",
                paste0("p.", quote_id(pcs)), paste0("r.", quote_id(groups)))
  query <- paste0("COPY (SELECT ", paste(selected, collapse = ", "),
                  " FROM read_parquet(", quote_str(intermediate[[2L]]), ") p JOIN ",
                  "read_parquet(", quote_str(intermediate[[1L]]), ") r ",
                  "USING (chromosome, position, allele_a, allele_b) ",
                  "ORDER BY chromosome, position, allele_a, allele_b) TO ",
                  quote_str(temporary),
                  " (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 32768)")
  DBI::dbExecute(con, query)
  numeric_columns <- c(pcs, groups)
  complete <- paste0("(", quote_id(numeric_columns), " IS NULL OR NOT isfinite(",
                     quote_id(numeric_columns), "))", collapse = " OR ")
  frequency_range <- paste0("(", quote_id(groups), " < 0 OR ", quote_id(groups),
                            " > 1)", collapse = " OR ")
  counts <- DBI::dbGetQuery(con, paste0(
    "SELECT count(*) AS n, count(DISTINCT (chromosome, position)) ",
    "AS unique_loci, count(*) FILTER (WHERE chromosome IS NULL OR position IS NULL ",
    "OR allele_a IS NULL OR allele_b IS NULL OR ", complete, " OR ", frequency_range,
    ") AS invalid FROM read_parquet(", quote_str(temporary), ")"))
  if (counts$n != row_counts[[1L]] || counts$n != row_counts[[2L]] ||
      counts$unique_loci != counts$n || counts$invalid != 0) {
    stop("ancestry Parquet requires one complete keyed row per locus",
         call. = FALSE)
  }
  expected$sha256[[6L]] <- unname(digest::digest(file = temporary, algo = "sha256"))
  utils::write.table(expected, temporary_receipt, sep = "\t",
                     row.names = FALSE, quote = FALSE)
  if (!file.rename(temporary, output)) stop("cannot publish ancestry Parquet", call. = FALSE)
  if (!file.rename(temporary_receipt, receipt)) {
    unlink(output)
    stop("cannot publish ancestry Parquet receipt", call. = FALSE)
  }
  output
}
