#' Stage the per-chromosome reference relations of the ROH ancestry evaluation.
#'
#' For each autosome, writes the cache files of registry artifacts
#' `roh_ancestry_chr<N>_reference_long` and `roh_ancestry_chr<N>_af_sites`, the
#' two frequency relations that the genome-wide ROH ancestry evaluation reads.
#'
#' * The long reference holds the rows of `ancestry_reference_grch38_parquet`
#'   on the chromosome, with the 21 group frequency columns unpivoted into
#'   `group_id` and `frequency` and the 16 PC columns left out, sorted by
#'   position, alleles and group.
#' * The AF sites hold the records of `roh_ancestry_chr<N>_children_bcf` on the
#'   chromosome with one ALT allele, whose position and alleles match a
#'   reference site that is not palindromic (A/T or C/G) in either orientation,
#'   and whose five INFO fields AF, AF_AFR, AF_AMR, AF_EAS and AF_EUR are all
#'   present: `chrom`, `pos`, `ref`, `alt` and the five INFO arrays, one row
#'   per site, sorted by position and alleles.
#'
#' This is the derivation that `benchmarks/roh_ancestry/prepare_chr20.R`
#' applies to chr20, with an explicit chromosome filter and a sort order; the
#' chr20 artifacts keep their own registry rows. Both registry rows of a
#' chromosome must name this function as their transform and name their
#' sources in `locator`.
#'
#' The children BCF must match the SHA-256 of its registry identity, and the
#' reference the `output_sha256` of its receipt (`<reference>.sources.tsv`).
#' Each output gets a receipt, `<output>.sources.tsv` with columns `field` and
#' `value`, that binds the SHA-256 of its sources (for the AF sites also of the
#' extension that reads the BCF), the SHA-256 of the derivation SQL, the DuckDB
#' version, the output SHA-256 and its row counts. A cached output is reused
#' only when its receipt matches; otherwise it is rebuilt and replaced. A
#' failed derivation publishes nothing.
#'
#' @param chromosomes Distinct autosome numbers.
#' @param extension_path DuckHTS extension that reads the BCF; `NULL` uses the
#'   extension bundled with Rduckhts.
#' @return A data frame with one row per artifact: `chromosome`, `id`, `path`,
#'   `rows` and `reused`.
#' @export
duckhts_bench_stage_roh_ancestry_references <- function(chromosomes, extension_path = NULL) {
  for (package in c("DBI", "duckdb", "Rduckhts")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(package, " is required for ROH ancestry reference staging", call. = FALSE)
    }
  }
  if (!is.numeric(chromosomes) || !length(chromosomes) || anyNA(chromosomes) ||
      any(chromosomes != round(chromosomes)) || any(chromosomes < 1 | chromosomes > 22) ||
      anyDuplicated(chromosomes)) {
    stop("chromosomes must be distinct autosome numbers from 1 to 22", call. = FALSE)
  }
  chromosomes <- as.integer(chromosomes)
  transform <- "duckhts_bench_stage_roh_ancestry_references"
  reference_id <- "ancestry_reference_grch38_parquet"
  registry <- duckhts_bench_registry()
  registry_row <- function(id) {
    row <- registry[registry$id == id, , drop = FALSE]
    if (nrow(row) != 1L) stop("unknown or non-unique benchmark artifact: ", id, call. = FALSE)
    row
  }
  # Check every requested chromosome's registry rows before any work.
  plan <- lapply(chromosomes, function(chromosome) {
    ids <- c(bcf = sprintf("roh_ancestry_chr%d_children_bcf", chromosome),
             long = sprintf("roh_ancestry_chr%d_reference_long", chromosome),
             sites = sprintf("roh_ancestry_chr%d_af_sites", chromosome))
    rows <- lapply(ids, registry_row)
    if (!identical(c(rows$long$transform, rows$sites$transform), c(transform, transform))) {
      stop("the registry rows of chromosome ", chromosome, " are not staged by ", transform,
           call. = FALSE)
    }
    if (!identical(rows$long$locator, paste0("artifact:", reference_id)) ||
        !identical(rows$sites$locator,
                   paste0("artifact:", ids[["bcf"]], ";artifact:", reference_id))) {
      stop("the registry rows of chromosome ", chromosome,
           " must name the reference and the children BCF as their sources", call. = FALSE)
    }
    bcf_sha256 <- duckhts_bench_identity_fields(rows$bcf$supplier_identity)["sha256"]
    if (is.na(bcf_sha256)) {
      stop("the registry identity of ", ids[["bcf"]], " lacks a SHA-256", call. = FALSE)
    }
    list(chromosome = chromosome, ids = ids, bcf_sha256 = unname(bcf_sha256))
  })

  sha256 <- function(path) unname(digest::digest(file = path, algo = "sha256"))
  reference <- duckhts_bench_artifact_path(reference_id)
  if (!file.exists(reference)) {
    stop("stage ", reference_id, " first: ", reference, call. = FALSE)
  }
  reference_sha256 <- sha256(reference)
  reference_receipt <- tryCatch(
    utils::read.delim(paste0(reference, ".sources.tsv"), colClasses = "character"),
    error = function(error) NULL)
  if (is.null(reference_receipt) ||
      !identical(reference_receipt$value[reference_receipt$field == "output_sha256"],
                 reference_sha256)) {
    stop("the ancestry reference does not match the output SHA-256 of its receipt: ",
         reference, call. = FALSE)
  }
  if (is.null(extension_path)) {
    extension_path <- system.file("duckhts_extension", "build", "duckhts.duckdb_extension",
                                  package = "Rduckhts")
  }
  if (length(extension_path) != 1L || !nzchar(extension_path) || !file.exists(extension_path)) {
    stop("DuckHTS extension not found: ", extension_path, call. = FALSE)
  }
  extension_path <- normalizePath(extension_path)
  extension_sha256 <- sha256(extension_path)

  con <- Rduckhts::rduckhts_connect(extension_path = extension_path)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  quote_str <- function(x) as.character(DBI::dbQuoteString(con, x))
  quote_id <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  spill <- file.path(dirname(duckhts_bench_artifact_path(plan[[1L]]$ids[["long"]])),
                     paste0("roh-ancestry-references.spill-", Sys.getpid()))
  # Partial outputs are renamed when published; any left at exit are removed.
  temporaries <- character()
  on.exit(unlink(c(spill, temporaries), recursive = TRUE), add = TRUE)
  DBI::dbExecute(con, "SET threads = 4")
  DBI::dbExecute(con, "SET memory_limit = '6GB'")
  DBI::dbExecute(con, paste0("SET temp_directory = ", quote_str(spill)))
  duckdb_version <- DBI::dbGetQuery(con, "SELECT version() AS v")$v
  columns <- DBI::dbGetQuery(con, paste0("DESCRIBE SELECT * FROM read_parquet(",
                                         quote_str(reference), ")"))$column_name
  keys <- c("chromosome", "position", "allele_a", "allele_b")
  groups <- setdiff(columns, c(keys, paste0("PC", seq_len(16L))))
  if (!all(c(keys, paste0("PC", seq_len(16L))) %in% columns) || length(groups) != 21L) {
    stop("the ancestry reference needs its keys, PC1-PC16 and 21 group columns", call. = FALSE)
  }

  # The derivations. Their SHA-256 is part of each receipt, so a change to
  # either SQL rebuilds the outputs it makes.
  long_sql <- paste0(
    "COPY (SELECT chromosome::VARCHAR AS chromosome, position::BIGINT AS position, ",
    "allele_a, allele_b, group_id, frequency::DOUBLE AS frequency FROM ",
    "(UNPIVOT reference_chromosome ON %s INTO NAME group_id VALUE frequency) ",
    "ORDER BY position, allele_a, allele_b, group_id) ",
    "TO %s (FORMAT PARQUET, COMPRESSION ZSTD)")
  palindromic <- paste0(
    "((r.allele_a IN ('A', 'T') AND r.allele_b IN ('A', 'T')) OR ",
    "(r.allele_a IN ('C', 'G') AND r.allele_b IN ('C', 'G')))")
  sites_sql <- paste0(
    "COPY (SELECT v.CHROM AS chrom, v.POS AS pos, v.REF AS ref, v.ALT[1] AS alt, ",
    "v.INFO_AF, v.INFO_AF_AFR, v.INFO_AF_AMR, v.INFO_AF_EAS, v.INFO_AF_EUR ",
    "FROM read_bcf(%s) v JOIN (SELECT DISTINCT position, allele_a, allele_b ",
    "FROM reference_chromosome) r ON v.POS = r.position AND ",
    "((v.REF = r.allele_a AND v.ALT[1] = r.allele_b) OR ",
    "(v.REF = r.allele_b AND v.ALT[1] = r.allele_a)) ",
    "WHERE try_cast(regexp_replace(v.CHROM, '^chr', '') AS INTEGER) = %d ",
    "AND len(v.ALT) = 1 AND NOT ", palindromic, " ",
    "AND v.INFO_AF IS NOT NULL AND v.INFO_AF_AFR IS NOT NULL ",
    "AND v.INFO_AF_AMR IS NOT NULL AND v.INFO_AF_EAS IS NOT NULL ",
    "AND v.INFO_AF_EUR IS NOT NULL ORDER BY pos, ref, alt) ",
    "TO %s (FORMAT PARQUET, COMPRESSION ZSTD)")
  derivation <- function(sql) digest::digest(sql, algo = "sha256", serialize = FALSE)

  # A receipt holds the expected fields, then output_sha256 and the counts.
  receipt_current <- function(output, receipt, expected) {
    if (!file.exists(output) || !file.exists(receipt)) return(FALSE)
    stored <- tryCatch(utils::read.delim(receipt, colClasses = "character"),
                       error = function(error) NULL)
    !is.null(stored) && nrow(stored) > nrow(expected) &&
      identical(stored$field[seq_len(nrow(expected))], expected$field) &&
      identical(stored$value[seq_len(nrow(expected))], expected$value) &&
      identical(stored$value[stored$field == "output_sha256"], sha256(output))
  }
  receipt_rows <- function(receipt) {
    stored <- utils::read.delim(receipt, colClasses = "character")
    as.numeric(stored$value[stored$field == "rows"])
  }
  # Output and receipt are published together, replacing stale copies; on a
  # failure the previous pair is restored (file.rename does not replace an
  # existing file on Windows).
  publish <- function(temporary, output, fields, receipt) {
    temporary_receipt <- paste0(receipt, ".partial-", Sys.getpid())
    on.exit(unlink(temporary_receipt), add = TRUE)
    utils::write.table(fields, temporary_receipt, sep = "\t", row.names = FALSE, quote = FALSE)
    duckhts_bench_publish_files(c(temporary, temporary_receipt), c(output, receipt))
  }
  count_fields <- function(names, values) {
    data.frame(field = names, value = as.character(values), stringsAsFactors = FALSE)
  }

  staged <- list()
  for (item in plan) {
    chromosome <- item$chromosome
    bcf <- duckhts_bench_artifact_path(item$ids[["bcf"]])
    if (!file.exists(bcf)) stop("stage ", item$ids[["bcf"]], " first: ", bcf, call. = FALSE)
    duckhts_bench_validate_identity(item$ids[["bcf"]], bcf)
    DBI::dbExecute(con, paste0(
      "CREATE OR REPLACE TEMP VIEW reference_chromosome AS SELECT * FROM read_parquet(",
      quote_str(reference), ") WHERE try_cast(regexp_replace(CAST(chromosome AS VARCHAR), ",
      "'^chr', '') AS INTEGER) = ", chromosome))
    reference_sites <- DBI::dbGetQuery(con, paste0(
      "SELECT count(*) AS sites, count(*) FILTER (WHERE NOT ", palindromic,
      ") AS non_palindromic FROM reference_chromosome r"))

    long_output <- duckhts_bench_artifact_path(item$ids[["long"]])
    long_receipt <- paste0(long_output, ".sources.tsv")
    long_expected <- count_fields(
      c("artifact_id", "chromosome", "reference_sha256", "derivation_sha256", "duckdb_version"),
      c(item$ids[["long"]], chromosome, reference_sha256, derivation(long_sql), duckdb_version))
    long_reused <- receipt_current(long_output, long_receipt, long_expected)
    if (long_reused) {
      long_rows <- receipt_rows(long_receipt)
    } else {
      dir.create(dirname(long_output), recursive = TRUE, showWarnings = FALSE)
      temporary <- paste0(long_output, ".partial-", Sys.getpid(), ".parquet")
      temporaries <- c(temporaries, temporary)
      DBI::dbExecute(con, sprintf(long_sql, paste(quote_id(groups), collapse = ", "),
                                  quote_str(temporary)))
      written <- DBI::dbGetQuery(con, paste0(
        "SELECT count(*) AS rows, count(DISTINCT (position, allele_a, allele_b)) AS sites, ",
        "count(DISTINCT (position, allele_a, allele_b, group_id)) AS keyed_rows, ",
        "count(*) FILTER (WHERE frequency IS NULL) AS missing FROM read_parquet(",
        quote_str(temporary), ")"))
      if (written$sites != reference_sites$sites || written$sites < 1 ||
          written$rows != 21 * written$sites || written$keyed_rows != written$rows ||
          written$missing != 0) {
        stop("the long reference of chromosome ", chromosome,
             " needs one row per reference site and group", call. = FALSE)
      }
      long_rows <- written$rows
      publish(temporary, long_output, rbind(long_expected, count_fields(
        c("output_sha256", "rows", "sites"),
        c(sha256(temporary), written$rows, written$sites))), long_receipt)
    }

    sites_output <- duckhts_bench_artifact_path(item$ids[["sites"]])
    sites_receipt <- paste0(sites_output, ".sources.tsv")
    sites_expected <- count_fields(
      c("artifact_id", "chromosome", "reference_sha256", "children_bcf_sha256",
        "extension_sha256", "derivation_sha256", "duckdb_version"),
      c(item$ids[["sites"]], chromosome, reference_sha256, item$bcf_sha256, extension_sha256,
        derivation(sites_sql), duckdb_version))
    sites_reused <- receipt_current(sites_output, sites_receipt, sites_expected)
    if (sites_reused) {
      sites_rows <- receipt_rows(sites_receipt)
    } else {
      dir.create(dirname(sites_output), recursive = TRUE, showWarnings = FALSE)
      temporary <- paste0(sites_output, ".partial-", Sys.getpid(), ".parquet")
      temporaries <- c(temporaries, temporary)
      DBI::dbExecute(con, sprintf(sites_sql, quote_str(bcf), chromosome, quote_str(temporary)))
      written <- DBI::dbGetQuery(con, paste0(
        "SELECT count(*) AS rows, count(DISTINCT (pos, ref, alt)) AS keyed_rows, ",
        "count(DISTINCT chrom) AS contigs FROM read_parquet(", quote_str(temporary), ")"))
      if (written$rows < 1 || written$keyed_rows != written$rows || written$contigs != 1) {
        stop("the AF sites of chromosome ", chromosome,
             " need one row per site on one contig", call. = FALSE)
      }
      sites_rows <- written$rows
      publish(temporary, sites_output, rbind(sites_expected, count_fields(
        c("output_sha256", "rows", "reference_sites"),
        c(sha256(temporary), written$rows, reference_sites$non_palindromic))), sites_receipt)
    }
    staged[[length(staged) + 1L]] <- data.frame(
      chromosome = chromosome, id = unname(item$ids[c("long", "sites")]),
      path = c(long_output, sites_output), rows = c(long_rows, sites_rows),
      reused = c(long_reused, sites_reused), stringsAsFactors = FALSE)
  }
  do.call(rbind, staged)
}
