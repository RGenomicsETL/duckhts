#' Estimate ancestry from a materialized `read_geno` relation
#'
#' Converts each diploid, biallelic call to ALT dosage / 2. A haploid,
#' partially missing, or multiallelic call has NULL dosage and is counted in
#' the missing-site audit. Original sample indices identify the calls unless a
#' `read_bcf_samples` relation is provided for original-header sample names.
#'
#' @param con Connection with DuckHTS loaded.
#' @param geno_table Materialized `read_geno` rows with CHROM, POS, REF, ALT
#'   and calls columns.
#' @param reference_table,loadings_table,correction_table Reference products.
#' @param samples_table Optional `read_bcf_samples` relation with sample_index
#'   and sample_name.
#' @param sum_to_one,min_cor Quality/constraint arguments for the solver.
#' @return A row per sample and group with the matching audit.
#' @export
rduckhts_ancestry_geno <- function(
  con, geno_table, reference_table, loadings_table, correction_table,
  samples_table = NULL, sum_to_one = TRUE, min_cor = 0.4
) {
  .somalier_validate_name(geno_table, "geno_table")
  if (!is.null(samples_table)) .somalier_validate_name(samples_table, "samples_table")
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
