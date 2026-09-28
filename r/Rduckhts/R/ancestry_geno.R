#' Estimate ancestry from a materialized `read_geno` relation
#'
#' Converts each diploid, biallelic call to ALT dosage / 2. A haploid,
#' partially missing, or multiallelic call has NULL dosage and is counted in
#' the missing-site audit. Original sample indices identify the calls unless a
#' `read_bcf_samples` relation is provided for original-header sample names.
#' That relation must give every selected sample index one distinct non-null name.
#'
#' @param con Connection with DuckHTS loaded.
#' @param geno_table Materialized `read_geno` rows with CHROM, POS, REF, ALT
#'   and calls columns. Each site must contain a call for every selected sample,
#'   including homozygous-reference calls.
#' @param reference_table,loadings_table,correction_table Reference products.
#' @param samples_table Optional `read_bcf_samples` relation with sample_index
#'   and sample_name.
#' @param sum_to_one,min_cor Quality/constraint arguments for the solver.
#' @param non_reference_only Required declaration of the `read_geno` setting used
#'   to create `geno_table`. Must be `FALSE`: sparse call lists omit zero dosages
#'   and cannot be identified from the materialized relation's schema.
#' @return A row per sample and group with the matching audit.
#' @export
rduckhts_ancestry_geno <- function(
  con, geno_table, reference_table, loadings_table, correction_table,
  samples_table = NULL, sum_to_one = TRUE, min_cor = 0.4,
  non_reference_only
) {
  if (missing(non_reference_only) || !identical(non_reference_only, FALSE)) {
    stop("ancestry requires dense read_geno calls: specify non_reference_only = FALSE; TRUE is unsupported",
         call. = FALSE)
  }
  .somalier_validate_name(geno_table, "geno_table")
  if (!is.null(samples_table)) {
    .somalier_validate_name(samples_table, "samples_table")
    samples <- sql_quote_identifier(con, samples_table)
    geno <- sql_quote_identifier(con, geno_table)
    invalid <- DBI::dbGetQuery(con, paste0(
      "WITH selected AS (SELECT DISTINCT c.sample_index FROM ", geno,
      " g, unnest(g.calls) u(c)), mapped AS (SELECT sample_index, ",
      "count(*) AS n, count(sample_name) AS named FROM ", samples,
      " GROUP BY sample_index) SELECT EXISTS (SELECT 1 FROM ", samples,
      " WHERE sample_index IS NULL OR sample_name IS NULL) OR EXISTS ",
      "(SELECT 1 FROM mapped WHERE n != 1 OR named != 1) OR ",
      "(SELECT count(*) != count(DISTINCT sample_name) FROM ", samples,
      ") OR EXISTS (SELECT 1 FROM selected i LEFT JOIN mapped m ",
      "ON i.sample_index = m.sample_index WHERE m.n IS NULL) AS bad"))$bad
    if (invalid) stop("samples_table must map each selected sample index to a distinct non-null name",
                      call. = FALSE)
  }
  frequencies <- basename(tempfile("rduckhts_ancestry_geno_"))
  on.exit(invisible(try(DBI::dbExecute(con, paste("DROP VIEW IF EXISTS",
    sql_quote_identifier(con, frequencies))), silent = TRUE)), add = TRUE)
  sample_join <- if (is.null(samples_table)) "" else paste0(
    " LEFT JOIN ", sql_quote_identifier(con, samples_table),
    " s ON s.sample_index = c.sample_index"
  )
  sample_id <- if (is.null(samples_table)) {
    "c.sample_index::VARCHAR"
  } else {
    "s.sample_name"
  }
  query <- paste0(
    "CREATE TEMP VIEW ", sql_quote_identifier(con, frequencies),
    " AS SELECT ", sample_id, " AS sample_id, g.CHROM AS chromosome, ",
    "g.POS AS position, g.REF AS allele_a, g.ALT[1] AS allele_b, ",
    "CASE WHEN len(g.ALT) = 1 AND len(c.alleles) = 2 ",
    "AND c.alleles[1] IN (0,1) AND c.alleles[2] IN (0,1) ",
    "THEN (c.alleles[1] + c.alleles[2])::DOUBLE END AS dosage ",
    "FROM ", sql_quote_identifier(con, geno_table), " g, unnest(g.calls) u(c)",
    sample_join
  )
  DBI::dbExecute(con, query)
  rduckhts_ancestry_proportions(
    con, frequencies, reference_table, loadings_table, correction_table,
    input_kind = "dosage", sum_to_one = sum_to_one, min_cor = min_cor
  )
}
