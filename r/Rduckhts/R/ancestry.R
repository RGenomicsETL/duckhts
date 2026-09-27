#' Estimate projected reference-group proportions from frequency or dosage relations
#'
#' `input_table` has sample_id, chromosome, position, allele_a, allele_b,
#' and frequency (or dosage for `input_kind = "dosage"`). Frequencies describe
#' allele_b. `reference_table` has the same site columns, group_id and frequency;
#' `loadings_table` has the site columns, pc and loading; `correction_table`
#' has pc and coefficient. Reference groups and PC loadings must be complete
#' and unique at every participating site. All tables must use the same genome
#' assembly, uppercase single-base alleles, and a common PC numbering scheme.
#' Reference orientation must agree with the loadings. Input sites with missing
#' frequencies/dosages, duplicates at a sample/locus, or palindromic allele pairs
#' are dropped; other biallelic SNVs can be matched directly, reversed, strand
#' complemented or both. Missing genotypes are dropped, not imputed as in
#' bigsnpr. Audit columns count all input rows including dropped rows.
#' Correlation gates return NULL proportions and a named status, never a forced
#' assignment. The native solver uses at most 30 groups and 64 PCs.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param input_table,reference_table,loadings_table,correction_table Caller-owned
#'   relation names on the current connection.
#' @param input_kind `"frequency"` (cohort or individual) or `"dosage"` (diploid
#'   genotype with dosage 0, 1, or 2; converted to dosage / 2).
#' @param sum_to_one Require coefficients to sum to one; otherwise at most one.
#' @param min_cor Minimum predicted-frequency correlation; defaults to 0.4,
#'   as in bigsnpr 1.12.21. Higher gates may be chosen for a particular panel.
#' @param table_name Optional destination table. NULL returns a data frame.
#' @param overwrite Whether an existing destination may be replaced.
#' @return A row per sample and reference group, with quality status and audit.
#' @export
rduckhts_ancestry_proportions <- function(
  con, input_table, reference_table, loadings_table, correction_table,
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
  names <- list(input_table, reference_table, loadings_table, correction_table)
  for (i in seq_along(names)) .somalier_validate_name(names[[i]], "relation")
  relations <- vapply(names, function(x) sql_quote_identifier(con, x), character(1L))
  dimensions <- DBI::dbGetQuery(con, paste0(
    "SELECT (SELECT count(DISTINCT group_id) FROM ", relations[[2L]], ") AS groups, ",
    "(SELECT count(DISTINCT pc) FROM ", relations[[3L]], ") AS pcs, ",
    "(SELECT count(*) FROM ", relations[[4L]], ") AS corrections, ",
    "(SELECT count(DISTINCT pc) FROM ", relations[[4L]], ") AS corrected_pcs"
  ))
  if (dimensions$groups < 1 || dimensions$groups > 30 ||
      dimensions$pcs < 1 || dimensions$pcs > 64 ||
      dimensions$corrections != dimensions$pcs ||
      dimensions$corrected_pcs != dimensions$pcs) {
    stop("reference requires 1..30 groups and 1..64 PCs with one correction per PC",
         call. = FALSE)
  }
  frequency <- if (input_kind == "dosage") "dosage / 2.0" else "frequency"
  query <- .ancestry_query(relations, frequency,
                           if (sum_to_one) "true" else "false",
                           .somalier_quote_number(con, min_cor))
  .somalier_publish_query(con, query, table_name, overwrite)
}

.ancestry_query <- function(relations, frequency, equality, gate) {
  input <- relations[[1L]]
  reference <- relations[[2L]]
  loadings <- relations[[3L]]
  correction <- relations[[4L]]
  paste0(
    "WITH raw AS (SELECT sample_id, chromosome, position, allele_a, allele_b, ",
    frequency, " AS f, count(*) OVER (PARTITION BY sample_id, chromosome, position) AS copies FROM ", input, "), ",
    "sites AS (SELECT chromosome, position, min(allele_a) AS ra, min(allele_b) AS rb ",
    "FROM ", reference, " GROUP BY chromosome, position ",
    "HAVING count(DISTINCT allele_a || '>' || allele_b) = 1 ",
    "AND count(*) = count(DISTINCT group_id) ",
    "AND count(DISTINCT group_id) = (SELECT count(DISTINCT group_id) FROM ", reference, ")), ",
    "classified AS (SELECT i.*, s.ra, s.rb, CASE ",
    "WHEN i.copies > 1 THEN 'duplicate' ",
    "WHEN i.allele_a IS NULL OR i.allele_b IS NULL OR length(i.allele_a) != 1 ",
    "OR length(i.allele_b) != 1 OR i.allele_a NOT IN ('A','C','G','T') ",
    "OR i.allele_b NOT IN ('A','C','G','T') OR i.allele_a = i.allele_b ",
    "THEN 'invalid_alleles' ",
    "WHEN i.f IS NULL OR NOT isfinite(i.f) OR i.f < 0 OR i.f > 1 THEN 'missing_frequency' ",
    "WHEN s.ra IS NULL THEN 'missing_reference' ",
    "WHEN s.ra || s.rb IN ('AT','TA','CG','GC') THEN 'ambiguous' ",
    "WHEN i.allele_a = s.ra AND i.allele_b = s.rb THEN 'direct' ",
    "WHEN i.allele_a = s.rb AND i.allele_b = s.ra THEN 'reversed' ",
    "WHEN translate(i.allele_a, 'ACGT', 'TGCA') = s.ra ",
    "AND translate(i.allele_b, 'ACGT', 'TGCA') = s.rb THEN 'flipped' ",
    "WHEN translate(i.allele_a, 'ACGT', 'TGCA') = s.rb ",
    "AND translate(i.allele_b, 'ACGT', 'TGCA') = s.ra THEN 'flipped_reversed' ",
    "ELSE 'allele_mismatch' END AS orientation ",
    "FROM raw i LEFT JOIN sites s USING (chromosome, position)), ",
    "audit_rows AS (SELECT *, CASE WHEN orientation IN ",
    "('direct','reversed','flipped','flipped_reversed') ",
    "AND NOT EXISTS (SELECT 1 FROM ", loadings, " l ",
    "WHERE l.chromosome = classified.chromosome ",
    "AND l.position = classified.position AND l.allele_a = classified.ra ",
    "AND l.allele_b = classified.rb GROUP BY l.chromosome, l.position ",
    "HAVING count(*) = count(DISTINCT l.pc) AND count(DISTINCT l.pc) = ",
    "(SELECT count(DISTINCT pc) FROM ", loadings, ")) ",
    "THEN 'missing_loading' ELSE orientation END AS match_status FROM classified), ",
    "aligned AS (SELECT *, CASE WHEN match_status IN ('reversed','flipped_reversed') ",
    "THEN 1 - f ELSE f END AS aligned_f FROM audit_rows ",
    "WHERE match_status IN ('direct','reversed','flipped','flipped_reversed')), ",
    "audit AS (SELECT sample_id, count(*) AS input_variants, ",
    "count(*) FILTER (WHERE match_status IN ('direct','reversed','flipped','flipped_reversed')) AS used_variants, ",
    "count(*) FILTER (WHERE match_status = 'reversed') AS reversed_variants, ",
    "count(*) FILTER (WHERE match_status IN ('flipped','flipped_reversed')) AS flipped_variants, ",
    "count(*) FILTER (WHERE match_status = 'duplicate') AS duplicate_variants, ",
    "count(*) FILTER (WHERE match_status = 'ambiguous') AS ambiguous_variants, ",
    "count(*) FILTER (WHERE match_status IN ('missing_reference','missing_loading','allele_mismatch')) ",
    "AS unmatched_variants, ",
    "count(*) FILTER (WHERE match_status = 'invalid_alleles') AS invalid_variants, ",
    "count(*) FILTER (WHERE match_status = 'missing_frequency') AS missing_variants ",
    "FROM audit_rows GROUP BY sample_id), ",
    "basis AS (SELECT a.sample_id, a.chromosome, a.position, a.aligned_f, ",
    "l.pc, l.loading FROM aligned a JOIN ", loadings, " l ",
    "ON l.chromosome = a.chromosome AND l.position = a.position ",
    "AND l.allele_a = a.ra AND l.allele_b = a.rb), ",
    "x AS (SELECT b.sample_id, b.pc, r.group_id, sum(b.loading * r.frequency) AS value ",
    "FROM basis b JOIN ", reference, " r ON r.chromosome = b.chromosome ",
    "AND r.position = b.position GROUP BY b.sample_id, b.pc, r.group_id), ",
    "y AS (SELECT b.sample_id, b.pc, sum(b.loading * b.aligned_f) * max(c.coefficient) AS value ",
    "FROM basis b JOIN ", correction, " c ON c.pc = b.pc GROUP BY b.sample_id, b.pc), ",
    "groups AS (SELECT sample_id, group_id, row_number() OVER ",
    "(PARTITION BY sample_id ORDER BY group_id) AS group_index FROM ",
    "(SELECT DISTINCT sample_id, group_id FROM x)), ",
    "xvec AS (SELECT g.sample_id, count(DISTINCT g.group_id)::INTEGER AS k, ",
    "list(coalesce(x.value, 0) ORDER BY y.pc, g.group_id) AS values_x ",
    "FROM groups g JOIN y ON y.sample_id = g.sample_id ",
    "LEFT JOIN x ON x.sample_id = g.sample_id AND x.group_id = g.group_id AND x.pc = y.pc ",
    "GROUP BY g.sample_id), ",
    "yvec AS (SELECT sample_id, list(value ORDER BY pc) AS values_y FROM y GROUP BY sample_id), ",
    "solved AS (SELECT xvec.sample_id, duckhts_ancestry_proportions(",
    "xvec.values_x, yvec.values_y, xvec.k, ", equality, ") AS q ",
    "FROM xvec JOIN yvec USING (sample_id)), ",
    "each_cor AS (SELECT a.sample_id, r.group_id, corr(a.aligned_f, r.frequency) AS cor_each ",
    "FROM aligned a JOIN ", reference, " r USING (chromosome, position) ",
    "GROUP BY a.sample_id, r.group_id), ",
    "pred_sites AS (SELECT a.sample_id, a.chromosome, a.position, ",
    "any_value(a.aligned_f) AS observed, sum(r.frequency * s.q[g.group_index]) AS predicted ",
    "FROM aligned a JOIN ", reference, " r USING (chromosome, position) ",
    "JOIN groups g ON g.sample_id = a.sample_id AND g.group_id = r.group_id ",
    "JOIN solved s ON s.sample_id = a.sample_id ",
    "GROUP BY a.sample_id, a.chromosome, a.position), ",
    "pred_cor AS (SELECT sample_id, corr(observed, predicted) AS cor_pred ",
    "FROM pred_sites GROUP BY sample_id), ",
    "quality AS (SELECT a.sample_id, p.cor_pred, CASE ",
    "WHEN a.used_variants = 0 THEN 'no_matched_variants' ",
    "WHEN (SELECT avg(e.cor_each) FROM each_cor e WHERE e.sample_id = a.sample_id) < -0.2 ",
    "THEN 'reversed_alleles' WHEN p.cor_pred IS NULL OR NOT isfinite(p.cor_pred) ",
    "OR p.cor_pred < ", gate,
    " THEN 'low_correlation' ELSE 'ok' END AS status FROM audit a ",
    "LEFT JOIN pred_cor p USING (sample_id)) ",
    "SELECT a.sample_id, g.group_id, CASE WHEN q.status = 'ok' ",
    "THEN s.q[g.group_index] ELSE NULL END AS proportion, ",
    "q.cor_pred, e.cor_each, q.status, a.input_variants, a.used_variants, ",
    "a.input_variants - a.used_variants AS dropped_variants, ",
    "a.reversed_variants, a.flipped_variants, a.duplicate_variants, ",
    "a.ambiguous_variants, a.unmatched_variants, a.invalid_variants, a.missing_variants, ",
    equality, " AS sum_to_one, ", gate, "::DOUBLE AS min_cor ",
    "FROM audit a JOIN quality q USING (sample_id) ",
    "LEFT JOIN groups g ON g.sample_id = a.sample_id ",
    "LEFT JOIN solved s ON s.sample_id = a.sample_id ",
    "LEFT JOIN each_cor e ON e.sample_id = a.sample_id AND e.group_id = g.group_id ",
    "ORDER BY a.sample_id, g.group_id"
  )
}
