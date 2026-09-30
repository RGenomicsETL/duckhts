#' Stage the GRCh38 derivation of the keyed ancestry reference as Parquet.
#'
#' Every locus of the keyed GRCh37 product (see
#' [duckhts_bench_stage_ancestry_parquet()]) is lifted with `duckdb_liftover`
#' using the registered GRCh37-to-GRCh38 chain and the registered source and
#' destination FASTAs. The coded allele of a locus is the source `allele_b`, the
#' allele that the group frequencies and PC loadings describe. Its destination
#' identity is `allele_b` on a forward chain block and its complement on a
#' reverse block, as `reverse_complemented` reports; the other allele follows the
#' same rule. The derived key is `allele_a` = the other allele and `allele_b` =
#' the coded allele in destination spelling, so every frequency and loading keeps
#' describing the same physical allele and its value is copied unchanged. This
#' is the reversal convention of `rduckhts_ancestry_proportions()`, which
#' matches a reversed input to an unchanged reference row. When the coded allele
#' is the destination reference, liftover has swapped the allele roles and the
#' locus is counted as swapped. Rewriting such a locus as `1 - f` with negated
#' loadings would not be equivalent under the per-PC correction coefficients, so
#' it is not done. Loci that do not map, are not single-base substitutions (the
#' ancestry contract selects biallelic SNVs, and an insertion or deletion has no
#' coded allele that survives left-alignment), land off chromosomes 1-22, keep
#' neither destination allele, or share a destination position with another
#' locus are dropped, and each is counted.
#'
#' The receipt (`<output>.sources.tsv`, columns `field` and `value`) binds the
#' SHA-256 of the source Parquet, chain, source FASTA and destination FASTA, the
#' registry derivation, the DuckDB, DuckHTS htslib and Rduckhts versions and the
#' output SHA-256 to the denominators: input loci, mapped, rejected by reason,
#' reverse-complemented, swapped, duplicate-destination dropped and output loci.
#' The companion `<output>` with suffix `.liftover.parquet` maps each output locus
#' to its GRCh37 source locus, source alleles and orientation flags; the receipt
#' binds its SHA-256 as `liftover_map_sha256`.
#' A cached product whose receipt does not match its sources is an error.
#'
#' @return The cached Parquet path.
#' @export
duckhts_bench_stage_ancestry_grch38_parquet <- function() {
  for (package in c("DBI", "duckdb", "Rduckhts")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(package, " is required for GRCh38 ancestry reference staging", call. = FALSE)
    }
  }
  registry <- duckhts_bench_registry()
  product <- registry[match("ancestry_reference_grch38_parquet", registry$id), , drop = FALSE]
  if (nrow(product) != 1L || is.na(product$transform) || !nzchar(product$transform)) {
    stop("missing GRCh38 ancestry Parquet derivation", call. = FALSE)
  }
  source_parquet <- duckhts_bench_stage_ancestry_parquet()
  chain <- duckhts_bench_fetch("liftover_grch37_grch38_chain")
  destination_fasta <- duckhts_bench_fetch("liftover_grch38_fasta")
  source_fasta <- duckhts_bench_artifact_path("liftover_grch37_fasta")
  if (!file.exists(source_fasta)) {
    stop("stage the liftover bundle with duckhts_bench_stage_liftover() first: ",
         source_fasta, call. = FALSE)
  }
  duckhts_bench_validate_identity("liftover_grch37_fasta", source_fasta)
  for (fasta in c(source_fasta, destination_fasta)) {
    if (!file.exists(paste0(fasta, ".fai"))) {
      stop("missing FASTA index: ", fasta, ".fai", call. = FALSE)
    }
  }
  sha256 <- function(path) unname(digest::digest(file = path, algo = "sha256"))
  con <- Rduckhts::rduckhts_connect()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  runtime <- DBI::dbGetQuery(con, paste(
    "SELECT version() AS duckdb, duckhts_htslib_version() AS htslib"))
  derivation <- digest::digest(product$transform, algo = "sha256", serialize = FALSE)
  output <- duckhts_bench_artifact_path("ancestry_reference_grch38_parquet")
  receipt <- paste0(output, ".sources.tsv")
  map_output <- sub("\\.parquet$", ".liftover.parquet", output)
  expected <- data.frame(
    field = c("source_parquet_sha256", "chain_sha256", "source_fasta_sha256",
              "destination_fasta_sha256", "derivation_sha256", "validation",
              "duckdb_version", "duckdb_r_version", "duckhts_htslib_version",
              "rduckhts_version"),
    value = c(sha256(source_parquet), sha256(chain), sha256(source_fasta),
              sha256(destination_fasta), derivation, "lifted_keyed_unique_v1",
              runtime$duckdb, as.character(utils::packageVersion("duckdb")),
              runtime$htslib, as.character(utils::packageVersion("Rduckhts"))),
    stringsAsFactors = FALSE)
  if (file.exists(output) && file.exists(receipt)) {
    stored <- tryCatch(utils::read.delim(receipt, colClasses = "character"),
                       error = function(e) NULL)
    if (!is.null(stored) &&
        identical(stored$field[seq_len(nrow(expected))], expected$field) &&
        identical(stored$value[seq_len(nrow(expected))], expected$value) &&
        "output_sha256" %in% stored$field &&
        identical(sha256(output), stored$value[stored$field == "output_sha256"]) &&
        file.exists(map_output) &&
        identical(sha256(map_output), stored$value[stored$field == "liftover_map_sha256"])) {
      return(output)
    }
    stop("cached GRCh38 ancestry Parquet identity does not match its sources: ",
         output, call. = FALSE)
  }
  if (file.exists(output) || file.exists(receipt) || file.exists(map_output)) {
    stop("incomplete GRCh38 ancestry Parquet staging: ", output, call. = FALSE)
  }
  quote_str <- function(x) as.character(DBI::dbQuoteString(con, x))
  quote_id <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(output, ".partial-", Sys.getpid(), ".parquet")
  spill <- paste0(output, ".spill-", Sys.getpid())
  temporary_map <- paste0(map_output, ".partial-", Sys.getpid(), ".parquet")
  temporary_receipt <- paste0(receipt, ".partial-", Sys.getpid())
  on.exit(unlink(c(temporary, temporary_map, temporary_receipt, spill), recursive = TRUE),
          add = TRUE)
  DBI::dbExecute(con, "SET threads=4")
  DBI::dbExecute(con, "SET memory_limit='6GB'")
  DBI::dbExecute(con, "SET preserve_insertion_order=false")
  DBI::dbExecute(con, paste0("SET temp_directory=", quote_str(spill)))
  source <- paste0("read_parquet(", quote_str(source_parquet), ")")
  columns <- DBI::dbGetQuery(con, paste0("DESCRIBE SELECT * FROM ", source))$column_name
  keys <- c("chromosome", "position", "allele_a", "allele_b")
  pcs <- paste0("PC", seq_len(16L))
  groups <- setdiff(columns, c(keys, pcs))
  if (!identical(columns, c(keys, pcs, groups)) || length(groups) != 21L) {
    stop("unexpected keyed ancestry reference schema", call. = FALSE)
  }
  DBI::dbExecute(con, paste0(
    "CREATE TEMP TABLE lift_input AS SELECT chromosome::VARCHAR AS chrom, ",
    "position::BIGINT AS pos, allele_a AS ref, allele_b AS alt FROM ", source))
  DBI::dbExecute(con, paste0(
    "CREATE TEMP TABLE lifted AS SELECT * FROM duckdb_liftover('lift_input', ",
    "'chrom', 'pos', ref_col := 'ref', alt_col := 'alt', chain_path := ",
    quote_str(chain), ", dst_fasta_ref := ", quote_str(destination_fasta),
    ", src_fasta_ref := ", quote_str(source_fasta), ")"))
  # The coded allele is allele_b; its destination identity follows the strand.
  complement <- function(x) paste0("translate(", x, ", 'ACGT', 'TGCA')")
  DBI::dbExecute(con, paste0(
    "CREATE TEMP TABLE classified AS SELECT src_chrom::UTINYINT AS chromosome, ",
    "src_pos::UINTEGER AS position, dest_chrom, dest_pos, dest_ref, dest_alt, ",
    "reverse_complemented, src_ref, src_alt, ",
    "regexp_extract(dest_chrom, '^(?:chr)?([0-9]+)$', 1) AS dest_number, ",
    "CASE WHEN reverse_complemented THEN ", complement("src_alt"), " ELSE src_alt END AS coded, ",
    "CASE WHEN reverse_complemented THEN ", complement("src_ref"), " ELSE src_ref END AS other, ",
    "CASE WHEN NOT mapped THEN 'liftover_' || coalesce(nullif(reject_reason, ''), 'unmapped') ",
    "WHEN length(src_ref) != 1 OR length(src_alt) != 1 THEN 'non_snv_source' ",
    "WHEN NOT regexp_matches(dest_chrom, '^(?:chr)?([1-9]|1[0-9]|2[0-2])$') ",
    "THEN 'non_autosomal_destination' ",
    "WHEN NOT ((dest_ref = CASE WHEN reverse_complemented THEN ", complement("src_ref"),
    " ELSE src_ref END AND dest_alt = CASE WHEN reverse_complemented THEN ",
    complement("src_alt"), " ELSE src_alt END) OR (dest_ref = CASE WHEN reverse_complemented ",
    "THEN ", complement("src_alt"), " ELSE src_alt END AND dest_alt = CASE WHEN ",
    "reverse_complemented THEN ", complement("src_ref"), " ELSE src_ref END)) ",
    "THEN 'destination_allele_mismatch' END AS rejection, ",
    "dest_ref = CASE WHEN reverse_complemented THEN ", complement("src_alt"),
    " ELSE src_alt END AS swapped FROM lifted"))
  DBI::dbExecute(con, paste0(
    "CREATE TEMP TABLE accepted AS SELECT c.*, ",
    "count(*) OVER (PARTITION BY dest_number, dest_pos) AS copies ",
    "FROM classified c WHERE rejection IS NULL"))
  counts <- DBI::dbGetQuery(con, paste0(
    "SELECT (SELECT count(*) FROM lift_input) AS input_loci, ",
    "(SELECT count(*) FROM lifted WHERE mapped) AS mapped, ",
    "(SELECT count(*) FROM classified WHERE rejection IS NOT NULL) AS rejected, ",
    "(SELECT count(*) FROM accepted WHERE reverse_complemented) AS reverse_complemented, ",
    "(SELECT count(*) FROM accepted WHERE swapped) AS swapped, ",
    "(SELECT count(*) FROM accepted WHERE copies > 1) AS duplicate_destination_dropped, ",
    "(SELECT count(*) FROM accepted WHERE copies = 1) AS output_loci"))
  reasons <- DBI::dbGetQuery(con,
    "SELECT rejection, count(*) AS n FROM classified WHERE rejection IS NOT NULL GROUP BY rejection ORDER BY rejection")
  lifted_rows <- DBI::dbGetQuery(con, "SELECT count(*) AS n FROM lifted")$n
  if (counts$input_loci != lifted_rows || sum(reasons$n) != counts$rejected ||
      counts$rejected + counts$duplicate_destination_dropped + counts$output_loci !=
        counts$input_loci) {
    stop("GRCh38 ancestry liftover denominators do not reconcile", call. = FALSE)
  }
  copied <- paste0("r.", quote_id(c(pcs, groups)))
  selected <- c("a.dest_number::UTINYINT AS chromosome", "a.dest_pos::UINTEGER AS position",
                "a.other AS allele_a", "a.coded AS allele_b", copied)
  DBI::dbExecute(con, paste0(
    "COPY (SELECT ", paste(selected, collapse = ", "), " FROM accepted a JOIN ",
    source, " r ON r.chromosome = a.chromosome AND r.position = a.position ",
    "WHERE a.copies = 1 ORDER BY 1, 2, 3, 4) TO ", quote_str(temporary),
    " (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 32768)"))
  DBI::dbExecute(con, paste0(
    "COPY (SELECT dest_number::UTINYINT AS chromosome, dest_pos::UINTEGER AS position, ",
    "chromosome AS source_chromosome, position AS source_position, ",
    "src_ref AS source_allele_a, src_alt AS source_allele_b, swapped, ",
    "reverse_complemented FROM accepted WHERE copies = 1 ORDER BY 1, 2) TO ",
    quote_str(temporary_map), " (FORMAT PARQUET, COMPRESSION ZSTD)"))
  numeric_columns <- c(pcs, groups)
  invalid <- paste0("(", quote_id(numeric_columns), " IS NULL OR NOT isfinite(",
                    quote_id(numeric_columns), "))", collapse = " OR ")
  frequency_range <- paste0("(", quote_id(groups), " < 0 OR ", quote_id(groups),
                            " > 1)", collapse = " OR ")
  checked <- DBI::dbGetQuery(con, paste0(
    "SELECT count(*) AS n, count(DISTINCT (chromosome, position)) AS unique_loci, ",
    "count(*) FILTER (WHERE chromosome IS NULL OR position IS NULL OR allele_a IS NULL ",
    "OR allele_b IS NULL OR allele_a = allele_b OR ", invalid, " OR ", frequency_range,
    ") AS invalid FROM read_parquet(", quote_str(temporary), ")"))
  if (checked$n != counts$output_loci || checked$unique_loci != checked$n ||
      checked$invalid != 0) {
    stop("GRCh38 ancestry Parquet requires one complete keyed row per locus",
         call. = FALSE)
  }
  denominators <- data.frame(
    field = c("input_loci", "mapped", "rejected",
              paste0("rejected_", reasons$rejection),
              "reverse_complemented", "swapped", "duplicate_destination_dropped",
              "output_loci"),
    value = as.character(c(counts$input_loci, counts$mapped, counts$rejected, reasons$n,
                           counts$reverse_complemented, counts$swapped,
                           counts$duplicate_destination_dropped, counts$output_loci)),
    stringsAsFactors = FALSE)
  full <- rbind(expected,
                data.frame(field = c("output_sha256", "liftover_map_sha256"),
                           value = c(sha256(temporary), sha256(temporary_map)),
                           stringsAsFactors = FALSE),
                denominators)
  utils::write.table(full, temporary_receipt, sep = "\t", row.names = FALSE, quote = FALSE)
  if (!file.rename(temporary, output)) {
    stop("cannot publish GRCh38 ancestry Parquet", call. = FALSE)
  }
  if (!file.rename(temporary_map, map_output)) {
    unlink(output)
    stop("cannot publish GRCh38 ancestry liftover map", call. = FALSE)
  }
  if (!file.rename(temporary_receipt, receipt)) {
    unlink(c(output, map_output))
    stop("cannot publish GRCh38 ancestry Parquet receipt", call. = FALSE)
  }
  output
}
