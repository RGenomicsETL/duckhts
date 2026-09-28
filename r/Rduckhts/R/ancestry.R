#' Estimate projected ancestry proportions from frequency or dosage relations
#'
#' Input columns are sample_id, chromosome, position, allele_a, allele_b and
#' frequency (or dosage), describing allele_b. Reference data may be a keyed
#' wide relation with consecutive PC1..PCn columns and one frequency column per group, or a long
#' relation with group_id and frequency plus a keyed pc/loading relation.
#' Corrections contain pc and coefficient. Groups (1..30), PCs (1..64),
#' frequencies and loadings must be complete and finite at each contributing
#' reference locus; group frequencies must lie in [0, 1]. Duplicate loci
#' matching the input are rejected. Long references require complete keyed
#' rows during pivoting. Inputs and
#' references must share an assembly, uppercase biallelic SNV alleles and
#' reference orientation. Numeric chromosomes match
#' after removing a `chr` prefix; other chromosome names match literally.
#' Input duplicates at a sample/locus, palindromic alleles and missing
#' frequencies are dropped; other alleles match directly, reversed, strand
#' complemented or both. Reversed alleles use 1-frequency. Missing genotypes
#' are not imputed. Audit counts include dropped physical input rows.
#' The solver and correlation gates use full precision; returned proportions
#' are rounded to seven decimals. Aligned scratch Parquet in `tempdir()` is
#' removed on return. Correlation failures return NULL proportions and a status.
#'
#' @param con DuckDB connection with DuckHTS loaded.
#' @param input_table,reference_table Caller-owned relation names.
#' @param loadings_table Long-format PC loading relation; for wide references,
#'   pass the correction relation here with `correction_table = NULL`.
#' @param correction_table PC correction relation for long references.
#' @param input_kind `frequency` or diploid `dosage` (divided by two).
#' @param sum_to_one Require coefficients to sum to one; otherwise at most one.
#' @param min_cor Minimum predicted-frequency correlation, default 0.4 as in
#'   bigsnpr 1.12.21. A stricter gate may suit a specific panel.
#' @param table_name Optional destination; NULL returns a data frame.
#' @param overwrite Replace an existing destination.
#' @return A row per sample and reference group, with status and matching audit.
#'   Data-frame results carry `aligned_bytes`, the size of temporary aligned
#'   Parquet written during the call.
#' @export
rduckhts_ancestry_proportions <- function(
  con, input_table, reference_table, loadings_table, correction_table = NULL,
  input_kind = c("frequency", "dosage"), sum_to_one = TRUE, min_cor = 0.4,
  table_name = NULL, overwrite = FALSE
) {
  input_kind <- match.arg(input_kind)
  if (!is.logical(sum_to_one) || length(sum_to_one) != 1L || is.na(sum_to_one)) {
    stop("sum_to_one must be one non-missing logical value", call. = FALSE)
  }
  if (!is.numeric(min_cor) || length(min_cor) != 1L || !is.finite(min_cor) ||
      min_cor < -1 || min_cor > 1) {
    stop("min_cor must be a finite correlation in [-1, 1]", call. = FALSE)
  }
  .somalier_validate_output(con, table_name, overwrite)
  names <- list(input_table, reference_table, loadings_table)
  if (!is.null(correction_table)) names <- c(names, list(correction_table))
  for (name in names) .somalier_validate_name(name, "relation")
  if (is.null(correction_table)) {
    return(.ancestry_bounded_proportions(
      con, input_table, reference_table, loadings_table, input_kind,
      sum_to_one, min_cor, table_name, overwrite))
  }
  wide <- basename(tempfile("ancestry_reference_"))
  on.exit(invisible(try(DBI::dbExecute(con, paste0("DROP VIEW IF EXISTS ",
    sql_quote_identifier(con, wide))), silent = TRUE)), add = TRUE)
  .ancestry_long_adapter(con, reference_table, loadings_table, wide)
  .ancestry_bounded_proportions(con, input_table, wide, correction_table,
                                input_kind, sum_to_one, min_cor, table_name, overwrite)
}

.ancestry_long_adapter <- function(con, reference_table, loadings_table, wide) {
  quote_id <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  r <- quote_id(reference_table)
  l <- quote_id(loadings_table)
  groups <- DBI::dbGetQuery(con, paste0("SELECT DISTINCT group_id FROM ", r,
                                        " ORDER BY group_id"))$group_id
  pcs <- DBI::dbGetQuery(con, paste0("SELECT DISTINCT pc FROM ", l,
                                      " ORDER BY pc"))$pc
  if (length(groups) < 1L || length(groups) > 30L || anyNA(groups) ||
      length(pcs) < 1L || length(pcs) > 64L ||
      !identical(as.integer(pcs), seq_along(pcs))) {
    stop("reference requires 1..30 groups and consecutive 1..64 PCs", call. = FALSE)
  }
  chromosome <- paste0("coalesce(try_cast(regexp_replace(chromosome::VARCHAR, '^chr', '') ",
                       "AS INTEGER)::VARCHAR, chromosome::VARCHAR)")
  ref_sites <- paste0("SELECT ", chromosome, " AS chromosome, position, allele_a, ",
                      "allele_b, group_id, frequency FROM ", r)
  pc_sites <- paste0("SELECT ", chromosome, " AS chromosome, position, allele_a, ",
                     "allele_b, pc, loading FROM ", l)
  keys <- "chromosome, position, allele_a, allele_b"
  invalid <- DBI::dbGetQuery(con, paste0(
    "WITH r AS (", ref_sites, "), l AS (", pc_sites, "), ",
    "rg AS (SELECT ", keys, ", count(*) AS n, count(DISTINCT group_id) AS distinct_n ",
    "FROM r GROUP BY ", keys, "), ",
    "lp AS (SELECT ", keys, ", count(*) AS n, count(DISTINCT pc) AS distinct_n ",
    "FROM l GROUP BY ", keys, ") ",
    "SELECT EXISTS (SELECT 1 FROM r WHERE chromosome IS NULL OR position IS NULL ",
    "OR allele_a IS NULL OR allele_b IS NULL OR group_id IS NULL ",
    "OR frequency IS NULL OR NOT isfinite(frequency)) OR EXISTS (SELECT 1 FROM rg ",
    "WHERE n != ", length(groups), " OR distinct_n != ", length(groups),
    ") AS bad_reference, EXISTS (SELECT 1 FROM l WHERE chromosome IS NULL ",
    "OR position IS NULL OR allele_a IS NULL OR allele_b IS NULL OR pc IS NULL ",
    "OR loading IS NULL OR NOT isfinite(loading)) OR EXISTS (SELECT 1 FROM lp WHERE n != ",
    length(pcs), " OR distinct_n != ", length(pcs), ") OR EXISTS (",
    "SELECT 1 FROM rg LEFT JOIN lp USING (", keys, ") WHERE lp.n IS NULL) ",
    "AS bad_loadings, EXISTS (SELECT 1 FROM r WHERE frequency < 0 OR frequency > 1) ",
    "AS out_of_range"))
  if (invalid$bad_reference) stop("reference sites require one finite frequency per group", call. = FALSE)
  if (invalid$bad_loadings) stop("reference sites require one finite loading per PC", call. = FALSE)
  if (invalid$out_of_range) stop("reference sites require frequencies in [0, 1]", call. = FALSE)
  quote_str <- function(x) as.character(DBI::dbQuoteString(con, x))
  frequency <- vapply(groups, function(g) paste0(
    "max(frequency) FILTER (WHERE group_id = ", quote_str(g), ") AS ", quote_id(g)),
    character(1L))
  loading <- vapply(seq_along(pcs), function(i) paste0(
    "max(loading) FILTER (WHERE pc = ", i, ") AS PC", i), character(1L))
  query <- paste0("CREATE TEMP VIEW ", quote_id(wide), " AS WITH r AS (", ref_sites,
    "), l AS (", pc_sites, ") SELECT r.*, l.* EXCLUDE (chromosome, position, allele_a, allele_b) ",
    "FROM (SELECT ", keys, ", ", paste(frequency, collapse = ", "),
    " FROM r GROUP BY ", keys, ") r JOIN (SELECT ", keys, ", ",
    paste(loading, collapse = ", "), " FROM l GROUP BY ", keys,
    ") l USING (", keys, ")")
  DBI::dbExecute(con, query)
}
