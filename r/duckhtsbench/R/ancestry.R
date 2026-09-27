#' Stage the keyed bigsnpr ancestry reference as Parquet.
#'
#' The two checksum-validated CSV products are the upstream authority. The
#' derived product stores 16 DOUBLE loadings and 21 DOUBLE frequencies per
#' (chromosome, position, allele_a, allele_b), ordered by that key. A receipt
#' binds its checksum to both source checksums and the DuckDB writer version.
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
  output <- duckhts_bench_artifact_path("ancestry_reference_parquet")
  receipt <- paste0(output, ".sources.tsv")
  version <- as.character(utils::packageVersion("duckdb"))
  expected <- data.frame(id = c(source_ids, "duckdb", "parquet"),
                         sha256 = c(source_hashes, version, ""))
  if (file.exists(output) && file.exists(receipt)) {
    stored <- tryCatch(utils::read.delim(receipt, colClasses = "character"),
                       error = function(e) NULL)
    if (!is.null(stored) && identical(stored$id, expected$id) &&
        identical(stored$sha256[1:3], expected$sha256[1:3]) &&
        identical(unname(digest::digest(file = output, algo = "sha256")),
                  stored$sha256[[4L]])) {
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
    names(utils::read.csv(gzfile(path), nrows = 0L, check.names = FALSE))
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
  spill <- paste0(output, ".spill-", Sys.getpid())
  on.exit(unlink(c(temporary, spill), recursive = TRUE), add = TRUE)
  DBI::dbExecute(con, "SET memory_limit='3GB'")
  DBI::dbExecute(con, paste0("SET temp_directory=", quote_str(spill)))
  selected <- c("p.chr::INTEGER AS chromosome", "p.pos::BIGINT AS position",
                "p.a0::VARCHAR AS allele_a", "p.a1::VARCHAR AS allele_b",
                paste0("p.", quote_id(pcs), "::DOUBLE AS ", quote_id(pcs)),
                paste0("r.", quote_id(groups), "::DOUBLE AS ", quote_id(groups)))
  csv_sources <- paste0("read_csv_auto(", quote_str(sources),
                        ", types={'a0':'VARCHAR','a1':'VARCHAR'})")
  query <- paste0("COPY (SELECT ", paste(selected, collapse = ", "),
                  " FROM ", csv_sources[[2L]], " p JOIN ", csv_sources[[1L]], " r ",
                  "ON p.chr=r.chr AND p.pos=r.pos AND p.a0=r.a0 AND p.a1=r.a1 ",
                  "ORDER BY p.chr, p.pos, p.a0, p.a1) TO ", quote_str(temporary),
                  " (FORMAT PARQUET, COMPRESSION ZSTD)")
  DBI::dbExecute(con, query)
  row_counts <- vapply(csv_sources, function(source) DBI::dbGetQuery(con, paste0(
    "SELECT count(*) AS n FROM ", source))$n, numeric(1L))
  counts <- DBI::dbGetQuery(con, paste0(
    "SELECT count(*) AS n, count(DISTINCT (chromosome, position, allele_a, allele_b)) ",
    "AS unique_keys FROM read_parquet(", quote_str(temporary), ")"))
  if (counts$n != row_counts[[1L]] || counts$n != row_counts[[2L]] ||
      counts$unique_keys != counts$n) {
    stop("ancestry Parquet requires a unique, complete key match between CSVs",
         call. = FALSE)
  }
  expected$sha256[[4L]] <- unname(digest::digest(file = temporary, algo = "sha256"))
  if (!file.rename(temporary, output)) stop("cannot publish ancestry Parquet", call. = FALSE)
  utils::write.table(expected, receipt, sep = "\t", row.names = FALSE, quote = FALSE)
  output
}
