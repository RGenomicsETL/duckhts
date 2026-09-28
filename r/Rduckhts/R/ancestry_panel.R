#' Create a deterministic ancestry site panel for BAM/CRAM counts
#'
#' Sites are biallelic unambiguous SNVs present in every reference group and
#' every PC loading. With `candidate_table`, only matching candidate loci and
#' alleles are retained. Otherwise one site is selected per chromosome/spaced
#' genomic window. Results are sorted by chromosome and position and assigned
#' dense zero-based panel ordinals. Keep the returned panel SHA-256 and the
#' reference release with the materialized panel as its versioned identity.
#'
#' @param con A connection with DuckHTS loaded.
#' @param reference_table,loadings_table Reference products on the connection.
#' @param table_name Destination committed panel relation, visible to the
#'   BAM/CRAM panel preparation path.
#' @param assembly Reference genome assembly label.
#' @param candidate_table Optional site panel, e.g. a Somalier panel with
#'   region, position, allele_a and allele_b.
#' @param spacing_bp Genomic window width for deterministic spacing.
#' @param max_sites Maximum number of selected sites.
#' @param overwrite Replace an existing destination.
#' @return A one-row data frame with panel SHA-256 and selected site count.
#' @export
rduckhts_ancestry_panel <- function(
  con, reference_table, loadings_table, table_name, assembly,
  candidate_table = NULL, spacing_bp = 5000, max_sites = 17000,
  overwrite = FALSE
) {
  .somalier_validate_output(con, table_name, overwrite)
  .somalier_validate_name(table_name, "table_name")
  .somalier_validate_name(reference_table, "reference_table")
  .somalier_validate_name(loadings_table, "loadings_table")
  if (!is.null(candidate_table)) .somalier_validate_name(candidate_table, "candidate_table")
  .somalier_validate_name(assembly, "assembly")
  spacing_bp <- .somalier_bounded_whole_number(spacing_bp, "spacing_bp", 1, 100000000)
  max_sites <- .somalier_bounded_whole_number(max_sites, "max_sites", 1, 1000000)
  reference <- sql_quote_identifier(con, reference_table)
  loadings <- sql_quote_identifier(con, loadings_table)
  candidates <- if (is.null(candidate_table)) {
    sprintf(paste0(
      "SELECT chromosome AS region, position, least(ra,rb) AS allele_a, ",
      "greatest(ra,rb) AS allele_b FROM sites QUALIFY ",
      "row_number() OVER (PARTITION BY chromosome, floor((position-1)/%d) ",
      "ORDER BY position, ra, rb) = 1"
    ), spacing_bp)
  } else {
    paste0(
      "SELECT s.chromosome AS region, s.position, s.ra AS allele_a, s.rb AS allele_b FROM sites s ",
      "JOIN (SELECT region, position, min(allele_a) AS allele_a, min(allele_b) AS allele_b ",
      "FROM ", sql_quote_identifier(con, candidate_table),
      " GROUP BY region, position HAVING count(*) = 1) p ",
      "ON p.region = s.chromosome AND p.position = s.position ",
      "WHERE (p.allele_a = s.ra AND p.allele_b = s.rb) ",
      "OR (p.allele_a = s.rb AND p.allele_b = s.ra) ",
      "OR (translate(p.allele_a,'ACGT','TGCA') = s.ra ",
      "AND translate(p.allele_b,'ACGT','TGCA') = s.rb) ",
      "OR (translate(p.allele_a,'ACGT','TGCA') = s.rb ",
      "AND translate(p.allele_b,'ACGT','TGCA') = s.ra)"
    )
  }
  query <- paste0(
    "WITH site_groups AS (SELECT chromosome, position, min(allele_a) AS ra, ",
    "min(allele_b) AS rb FROM ", reference, " GROUP BY chromosome, position ",
    "HAVING count(DISTINCT allele_a || '>' || allele_b) = 1 ",
    "AND count(*) = count(DISTINCT group_id) ",
    "AND count(DISTINCT group_id) = (SELECT count(DISTINCT group_id) FROM ", reference, ") ",
    "AND length(min(allele_a)) = 1 AND length(min(allele_b)) = 1 ",
    "AND min(allele_a) IN ('A','C','G','T') AND min(allele_b) IN ('A','C','G','T') ",
    "AND min(allele_a) != min(allele_b) ",
    "AND min(allele_a) || min(allele_b) NOT IN ('AT','TA','CG','GC')), ",
    "sites AS (SELECT s.* FROM site_groups s WHERE EXISTS (SELECT 1 FROM ", loadings, " l ",
    "WHERE l.chromosome = s.chromosome AND l.position = s.position ",
    "AND l.allele_a = s.ra AND l.allele_b = s.rb ",
    "GROUP BY l.chromosome, l.position HAVING count(*) = count(DISTINCT l.pc) ",
    "AND count(DISTINCT l.pc) = (SELECT count(DISTINCT pc) FROM ", loadings, "))), ",
    "chosen AS (", candidates, "), limited AS (SELECT region, position, ",
    "least(allele_a,allele_b) AS allele_a, greatest(allele_a,allele_b) AS allele_b ",
    "FROM chosen ORDER BY region, position, allele_a, allele_b LIMIT ", max_sites, ") ",
    "SELECT ", sql_quote_string(con, assembly), " AS assembly, ",
    "(row_number() OVER (ORDER BY region, position, allele_a, allele_b)-1)::UBIGINT AS site_index, ",
    "region, position::UBIGINT AS position, allele_a, allele_b FROM limited"
  )
  # Validate in a scratch TEMP table before publishing, so a failed panel leaves no
  # destination table; this works inside or outside a caller's transaction.
  .duckhts_check_table_target(con, table_name, overwrite)
  scratch <- gsub("[^A-Za-z0-9_]", "_", basename(tempfile("__duckhts_ancestry_panel_")))
  on.exit(try(DBI::dbExecute(con, paste("DROP TABLE IF EXISTS",
                                        sql_quote_identifier(con, scratch))),
              silent = TRUE), add = TRUE)
  DBI::dbExecute(con, paste("CREATE TEMP TABLE", sql_quote_identifier(con, scratch), "AS", query))
  summary <- DBI::dbGetQuery(con, sprintf(
    "SELECT count(*) AS sites, duckhts_somalier_panel_sha256(%s) AS panel_sha256 FROM %s",
    sql_quote_string(con, scratch), sql_quote_identifier(con, scratch)
  ))
  .duckhts_create_table(con, table_name,
                        paste("SELECT * FROM", sql_quote_identifier(con, scratch)), overwrite)
  summary
}
