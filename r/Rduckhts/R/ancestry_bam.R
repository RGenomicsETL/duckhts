#' Project ancestry directly from an indexed BAM or CRAM
#'
#' The panel must be a committed, versioned canonical panel intersected with
#' the reference products. It has Somalier panel columns assembly, site_index,
#' region, position, allele_a, and allele_b. Each count is extracted by the
#' DuckHTS panel reader on the caller's connection. The count relation and
#' frequency view are discarded on return; no private connection is created.
#'
#' @param con Connection with DuckHTS loaded.
#' @param source_path Indexed BAM/CRAM path.
#' @param sample_id Sample identifier in the output.
#' @param reference_path FASTA path for reference checks and CRAM decoding.
#' @param panel_table Committed panel relation with canonical site identity.
#' @param reference_table Keyed wide or long reference relation.
#' @param loadings_table PC loadings for long references; pass the correction
#'   relation here for wide references.
#' @param correction_table PC correction relation for long references; NULL for wide.
#' @param frequency_method `"allele_fraction"` for B/(A+B) or `"called_genotype"`
#'   for the Somalier balance-rule genotype divided by two.
#' @param min_depth Minimum measured A+B depth to retain a site.
#' @param min_cor Minimum predicted-frequency correlation, default 0.4.
#' @param ... Additional arguments to `rduckhts_somalier_bam_counts`, such as
#'   index_path, map/base quality thresholds and worker_count.
#' @return Proportions and matching audit with frequency_method and min_depth.
#' @export
rduckhts_ancestry_bam <- function(
  con, source_path, sample_id, reference_path, panel_table,
  reference_table, loadings_table, correction_table = NULL,
  frequency_method = c("allele_fraction", "called_genotype"),
  min_depth = 7, min_cor = 0.4, ...
) {
  frequency_method <- match.arg(frequency_method)
  min_depth <- .somalier_bounded_whole_number(min_depth, "min_depth", 1, 100000000)
  counts <- basename(tempfile("rduckhts_ancestry_counts_"))
  frequencies <- basename(tempfile("rduckhts_ancestry_frequency_"))
  on.exit({
    invisible(try(DBI::dbExecute(con, paste("DROP VIEW IF EXISTS",
      sql_quote_identifier(con, frequencies))), silent = TRUE))
    invisible(try(DBI::dbExecute(con, paste("DROP TABLE IF EXISTS",
      sql_quote_identifier(con, counts))), silent = TRUE))
  }, add = TRUE)
  rduckhts_somalier_bam_counts(
    con, source_path, sample_id, reference_path, panel_table = panel_table,
    table_name = counts, ...
  )
  frequency <- if (frequency_method == "allele_fraction") {
    sprintf("CASE WHEN a IS NOT NULL AND b IS NOT NULL AND a+b >= %d AND a+b > 0 THEN b::DOUBLE/(a+b) END", min_depth)
  } else {
    sprintf(paste0(
      "CASE WHEN a IS NOT NULL AND b IS NOT NULL AND a+b >= %d THEN ",
      "NULLIF(duckhts_somalier_classify(a,b,other,%d,0.3,0.01).genotype,-1)/2.0 END"
    ), min_depth, min_depth)
  }
  DBI::dbExecute(con, sprintf(
    "CREATE TEMP VIEW %s AS SELECT sample_id, region AS chromosome, position, allele_a, allele_b, %s AS frequency FROM %s",
    sql_quote_identifier(con, frequencies), frequency, sql_quote_identifier(con, counts)
  ))
  result <- rduckhts_ancestry_proportions(
    con, frequencies, reference_table, loadings_table, correction_table,
    min_cor = min_cor
  )
  result$frequency_method <- frequency_method
  result$min_depth <- min_depth
  result
}
