#' Stage the read-name bucket counts of the depth-fraction benchmark.
#'
#' For each registry row of workload `depth-fraction` whose transform is
#' `duckhts_bench_stage_depth_fraction`, reads the regions of one alignment
#' file and writes the number of reads for each bin and read-name bucket as a
#' Parquet file. A read is counted at its leftmost position when no bit of
#' `exclude_flags` is set and its MAPQ is at least `min_mapq`. Its bucket is
#' the MD5 of its read name, as a number, modulo `buckets`; buckets from 0 to
#' `kept_buckets - 1` are kept. Both mates of a pair share a bucket, and a set
#' of `k` buckets is a read-name hash subsample of `k / buckets` of the reads,
#' so sums of disjoint bucket sets give downsampled data and mixtures at known
#' shares without reading the alignments again.
#'
#' The row's `locator` names the alignment artifact and the reference FASTA
#' (`artifact:<alignments>;artifact:<reference>`). The alignment file is read
#' from its cached copy when there is one, otherwise from its registered URL
#' with range requests. The row's `supplier_identity` names `sample`,
#' `regions` (comma-separated), `bin_width`, `buckets`, `kept_buckets`,
#' `min_mapq` and `exclude_flags`; when it also names `rows` and `reads`, the
#' staged file must have them before it is published.
#'
#' Each output gets a receipt, `<output>.sources.tsv` with columns `field` and
#' `value`, binding the artifact identity, the derivation SQL, the DuckDB
#' version, the output SHA-256 and the counts. A cached output is reused when
#' its receipt matches; otherwise it is rebuilt. A failed derivation publishes
#' nothing.
#'
#' @param ids Artifact ids to stage; `NULL` stages every row of the workload.
#' @param extension_path DuckHTS extension; `NULL` uses the extension bundled
#'   with Rduckhts.
#' @param threads DuckDB threads.
#' @return A data frame with one row per artifact: `id`, `path`, `rows`,
#'   `reads` and `reused`.
#' @export
duckhts_bench_stage_depth_fraction <- function(ids = NULL, extension_path = NULL, threads = 4L) {
  for (package in c("DBI", "duckdb", "Rduckhts")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(package, " is required for depth-fraction staging", call. = FALSE)
    }
  }
  transform <- "duckhts_bench_stage_depth_fraction"
  registry <- duckhts_bench_registry()
  rows <- registry[registry$workload == "depth-fraction" & registry$transform == transform, ,
                   drop = FALSE]
  if (!is.null(ids)) {
    if (!is.character(ids) || anyNA(ids) || !all(ids %in% rows$id)) {
      stop("ids must name depth-fraction artifacts staged by ", transform, call. = FALSE)
    }
    rows <- rows[rows$id %in% ids, , drop = FALSE]
  }
  if (!nrow(rows)) stop("no depth-fraction artifact to stage", call. = FALSE)
  required <- c("sample", "regions", "bin_width", "buckets", "kept_buckets", "min_mapq",
                "exclude_flags")
  whole <- c("bin_width", "buckets", "kept_buckets", "min_mapq", "exclude_flags")
  # Check every row before any work.
  plan <- lapply(seq_len(nrow(rows)), function(i) {
    row <- rows[i, , drop = FALSE]
    fields <- duckhts_bench_identity_fields(row$supplier_identity)
    if (!all(required %in% names(fields))) {
      stop("the identity of ", row$id, " lacks: ",
           paste(setdiff(required, names(fields)), collapse = ", "), call. = FALSE)
    }
    numbers <- suppressWarnings(as.numeric(fields[whole]))
    if (anyNA(numbers) || any(numbers != floor(numbers)) || any(numbers < 0) ||
        numbers[[1L]] < 1 || numbers[[2L]] < 1 || numbers[[3L]] < 1 ||
        numbers[[3L]] > numbers[[2L]]) {
      stop("the identity of ", row$id, " needs whole numbers with 1 <= kept_buckets <= buckets",
           call. = FALSE)
    }
    if (!grepl("^[A-Za-z0-9_.]+(,[A-Za-z0-9_.]+)*$", fields[["regions"]])) {
      stop("the regions of ", row$id, " must be comma-separated contig names", call. = FALSE)
    }
    inputs <- strsplit(row$locator, ";", fixed = TRUE)[[1L]]
    if (length(inputs) != 2L || !all(startsWith(inputs, "artifact:"))) {
      stop(row$id, " must name its alignments and its reference as artifacts", call. = FALSE)
    }
    inputs <- sub("^artifact:", "", inputs)
    if (!all(inputs %in% registry$id)) {
      stop(row$id, " names an unknown input", call. = FALSE)
    }
    list(id = row$id, fields = fields, alignments = inputs[[1L]], reference = inputs[[2L]])
  })

  sha256 <- function(path) unname(digest::digest(file = path, algo = "sha256"))
  is_url <- function(path) grepl("^[A-Za-z][A-Za-z0-9+.-]*://", path)
  source_of <- function(id) {
    cached <- duckhts_bench_artifact_path(id)
    if (file.exists(cached)) return(normalizePath(cached))
    locator <- registry$locator[registry$id == id]
    if (!is_url(locator)) {
      stop("stage ", id, " first, or register a URL for it: ", cached, call. = FALSE)
    }
    locator
  }
  if (is.null(extension_path)) {
    extension_path <- system.file("duckhts_extension", "build", "duckhts.duckdb_extension",
                                  package = "Rduckhts")
  }
  if (length(extension_path) != 1L || !nzchar(extension_path) || !file.exists(extension_path)) {
    stop("DuckHTS extension not found: ", extension_path, call. = FALSE)
  }
  extension_path <- normalizePath(extension_path)

  con <- Rduckhts::rduckhts_connect(extension_path = extension_path)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  quote_str <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbExecute(con, paste0("SET threads = ", as.integer(threads)))
  duckdb_version <- DBI::dbGetQuery(con, "SELECT version() AS v")$v
  # htslib saves the index of a remote file in the working directory.
  scratch <- tempfile("depth-fraction-stage-")
  dir.create(scratch)
  previous <- setwd(scratch)
  on.exit({
    setwd(previous)
    unlink(scratch, recursive = TRUE)
  }, add = TRUE)

  # The derivation. Its SHA-256 is part of each receipt, so a change rebuilds
  # the outputs. The placeholders are the sample, the source, the reference,
  # the regions, the bin width (twice), the buckets (twice), the exclusion
  # flags, the minimum MAPQ, the kept buckets and the output.
  derivation_sql <- paste0(
    "COPY (SELECT %s AS sample_id, RNAME AS chrom, ",
    "(POS - 1) // %s * %s AS start, ",
    "CAST(md5_number(QNAME) %% %s AS INTEGER) AS bucket, count(*) AS reads ",
    "FROM read_bam(%s, reference := %s, region := %s) ",
    "WHERE (FLAG & %s) = 0 AND MAPQ >= %s AND md5_number(QNAME) %% %s < %s ",
    "GROUP BY ALL ORDER BY chrom, start, bucket) TO %s (FORMAT PARQUET, COMPRESSION ZSTD)")
  derivation_sha256 <- digest::digest(derivation_sql, algo = "sha256", serialize = FALSE)
  field_frame <- function(names, values) {
    data.frame(field = names, value = as.character(values), stringsAsFactors = FALSE)
  }
  receipt_current <- function(output, receipt, expected) {
    if (!file.exists(output) || !file.exists(receipt)) return(FALSE)
    stored <- tryCatch(utils::read.delim(receipt, colClasses = "character"),
                       error = function(error) NULL)
    !is.null(stored) && nrow(stored) > nrow(expected) &&
      identical(stored$field[seq_len(nrow(expected))], expected$field) &&
      identical(stored$value[seq_len(nrow(expected))], expected$value) &&
      identical(stored$value[stored$field == "output_sha256"], sha256(output))
  }
  check_counts <- function(id, fields, counts) {
    for (field in c("rows", "reads")) {
      if (field %in% names(fields) && as.numeric(fields[[field]]) != counts[[field]]) {
        stop(id, " differs from its registered identity: ", field, " ",
             format(counts[[field]], scientific = FALSE), " instead of ", fields[[field]],
             call. = FALSE)
      }
    }
  }
  receipt_value <- function(receipt, field) {
    stored <- utils::read.delim(receipt, colClasses = "character")
    as.numeric(stored$value[stored$field == field])
  }

  temporaries <- character()
  on.exit(unlink(temporaries), add = TRUE)
  staged <- list()
  for (item in plan) {
    fields <- item$fields
    reference <- duckhts_bench_artifact_path(item$reference)
    if (!file.exists(reference)) {
      stop("stage ", item$reference, " first: ", reference, call. = FALSE)
    }
    reference <- normalizePath(reference)
    source <- source_of(item$alignments)
    output <- duckhts_bench_artifact_path(item$id)
    receipt <- paste0(output, ".sources.tsv")
    # The receipt binds the derivation parameters, not the expected counts, so
    # that recording the counts in the registry does not rebuild the output.
    parameters <- paste0(required, "=", fields[required], collapse = ";")
    expected <- field_frame(
      c("artifact_id", "supplier_identity", "source", "derivation_sha256", "duckdb_version"),
      c(item$id, parameters, registry$locator[registry$id == item$alignments], derivation_sha256,
        duckdb_version))
    reused <- receipt_current(output, receipt, expected)
    if (reused) {
      counts <- data.frame(rows = receipt_value(receipt, "rows"),
                           reads = receipt_value(receipt, "reads"))
      check_counts(item$id, fields, counts)
    } else {
      dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
      output <- file.path(normalizePath(dirname(output)), basename(output))
      receipt <- paste0(output, ".sources.tsv")
      temporary <- paste0(output, ".partial-", Sys.getpid(), ".parquet")
      temporaries <- c(temporaries, temporary)
      DBI::dbExecute(con, sprintf(
        derivation_sql, quote_str(fields[["sample"]]), fields[["bin_width"]],
        fields[["bin_width"]], fields[["buckets"]], quote_str(source), quote_str(reference),
        quote_str(fields[["regions"]]), fields[["exclude_flags"]], fields[["min_mapq"]],
        fields[["buckets"]], fields[["kept_buckets"]], quote_str(temporary)))
      counts <- DBI::dbGetQuery(con, paste0(
        "SELECT count(*)::DOUBLE AS rows, coalesce(sum(reads), 0)::DOUBLE AS reads ",
        "FROM read_parquet(", quote_str(temporary), ")"))
      check_counts(item$id, fields, counts)
      temporary_receipt <- paste0(receipt, ".partial-", Sys.getpid())
      temporaries <- c(temporaries, temporary_receipt)
      utils::write.table(
        rbind(expected, field_frame(c("output_sha256", "rows", "reads"),
                                    c(sha256(temporary), format(counts$rows, scientific = FALSE),
                                      format(counts$reads, scientific = FALSE)))),
        temporary_receipt, sep = "\t", row.names = FALSE, quote = FALSE)
      duckhts_bench_publish_files(c(temporary, temporary_receipt), c(output, receipt))
    }
    staged[[length(staged) + 1L]] <- data.frame(
      id = item$id, path = output, rows = counts$rows, reads = counts$reads, reused = reused,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, staged)
}
