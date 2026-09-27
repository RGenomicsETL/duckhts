#' Infer ancestry proportions from a keyed wide reference relation.
#'
#' The reference relation has `chromosome` (integer), `position` (one-based),
#' `allele_a`, `allele_b`, `PC1` through `PC16`, and one DOUBLE frequency column
#' per reference group. It can be a view over the Parquet product staged by
#' `duckhtsbench::duckhts_bench_stage_ancestry_parquet()`. Each site must have
#' exactly one reference row. The input and correction schemas, audit statuses,
#' and QP solver contract are those of [rduckhts_ancestry_proportions()].
#'
#' @param con DuckDB DBI connection with DuckHTS loaded.
#' @param input_table Relation with per-sample allele frequencies or dosages.
#' @param reference_table Wide reference relation (or view over staged Parquet).
#' @param correction_table Relation of PC number and correction coefficient.
#' @param input_kind Either `frequency` or `dosage`.
#' @param sum_to_one Constrain the nonnegative solution to sum to one.
#' @param min_cor Minimum predicted-frequency correlation for an `ok` status.
#' @param table_name Optional destination table; NULL returns a data frame.
#' @param overwrite Whether an existing destination may be replaced.
#' @return A row per sample and reference group, with quality status and audit.
#' @export
rduckhts_ancestry_proportions_wide <- function(
  con, input_table, reference_table, correction_table,
  input_kind = c("frequency", "dosage"), sum_to_one = TRUE, min_cor = 0.4,
  table_name = NULL, overwrite = FALSE
) {
  input_kind <- match.arg(input_kind)
  .ancestry_validate_options(sum_to_one, min_cor)
  .somalier_validate_output(con, table_name, overwrite)
  names <- list(input_table, reference_table, correction_table)
  for (name in names) .somalier_validate_name(name, "relation")
  relations <- vapply(names, function(x) sql_quote_identifier(con, x), character(1L))
  fields <- DBI::dbListFields(con, reference_table)
  keys <- c("chromosome", "position", "allele_a", "allele_b")
  pcs <- paste0("PC", seq_len(16L))
  groups <- sort(setdiff(fields, c(keys, pcs)))
  if (!all(c(keys, pcs) %in% fields) || length(groups) < 1L ||
      length(groups) > 30L || anyDuplicated(fields)) {
    stop("wide reference requires keyed sites, PC1..PC16 and 1..30 groups", call. = FALSE)
  }
  corrections <- DBI::dbGetQuery(con, paste0("SELECT pc, coefficient FROM ", relations[[3L]],
                                             " ORDER BY pc"))
  if (nrow(corrections) != length(pcs) ||
      !identical(as.integer(corrections$pc), seq_along(pcs)) ||
      any(!is.finite(corrections$coefficient))) {
    stop("correction requires one finite coefficient for every PC", call. = FALSE)
  }
  frequency <- if (input_kind == "dosage") "dosage / 2.0" else "frequency"
  query <- .ancestry_wide_query(con, relations, groups, corrections$coefficient,
                                frequency, sum_to_one, min_cor)
  .somalier_publish_query(con, query, table_name, overwrite)
}

.ancestry_wide_query <- function(con, relations, groups, coefficients,
                                 frequency, sum_to_one, min_cor) {
  input <- relations[[1L]]
  reference <- relations[[2L]]
  quote_id <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  quote_str <- function(x) as.character(DBI::dbQuoteString(con, x))
  pc <- paste0("r.", quote_id(paste0("PC", seq_along(coefficients))))
  group <- paste0("r.", quote_id(groups))
  # The solver reads X in PC-major, group-minor order.
  x_names <- unlist(lapply(seq_along(pc), function(k) {
    paste0("x", k, "_", seq_along(group))
  }), use.names = FALSE)
  x_terms <- unlist(lapply(seq_along(pc), function(k) {
    paste0("sum(", pc[[k]], " * ", group, ") AS ",
           paste0("x", k, "_", seq_along(group)))
  }), use.names = FALSE)
  y_names <- paste0("y", seq_along(pc))
  y_terms <- paste0("sum(", pc, " * a.aligned_f) * ",
                    vapply(coefficients, function(x) .somalier_quote_number(con, x),
                           character(1L)), " AS ", y_names)
  cor_names <- paste0("cor", seq_along(group))
  cor_terms <- paste0("corr(a.aligned_f, ", group, ") AS ", cor_names)
  values_x <- paste0("m.", x_names, collapse = ", ")
  values_y <- paste0("m.", y_names, collapse = ", ")
  values_cor <- paste0("m.", cor_names, collapse = ", ")
  predicted <- paste0("(", group, " * s.q[", seq_along(group), "])",
                      collapse = " + ")
  group_rows <- paste0("(", seq_along(groups), ", ",
                       vapply(groups, quote_str, character(1L)), ")", collapse = ", ")
  equality <- if (sum_to_one) "true" else "false"
  gate <- .somalier_quote_number(con, min_cor)
  paste0(
    "WITH raw AS (SELECT sample_id, chromosome, position, allele_a, allele_b, ",
    "try_cast(regexp_replace(chromosome, '^chr', '') AS INTEGER) AS ref_chr, ",
    frequency, " AS f, count(*) OVER (PARTITION BY sample_id, chromosome, position) AS copies ",
    "FROM ", input, "), ",
    "sites AS (SELECT chromosome, position, min(allele_a) AS ra, min(allele_b) AS rb ",
    "FROM ", reference, " GROUP BY chromosome, position HAVING count(*) = 1), ",
    "classified AS (SELECT i.*, s.ra, s.rb, CASE ",
    "WHEN i.copies > 1 THEN 'duplicate' ",
    "WHEN i.allele_a IS NULL OR i.allele_b IS NULL OR length(i.allele_a) != 1 ",
    "OR length(i.allele_b) != 1 OR i.allele_a NOT IN ('A','C','G','T') ",
    "OR i.allele_b NOT IN ('A','C','G','T') OR i.allele_a = i.allele_b ",
    "THEN 'invalid_alleles' ",
    "WHEN i.allele_a || i.allele_b IN ('AT','TA','CG','GC') ",
    "OR s.ra || s.rb IN ('AT','TA','CG','GC') THEN 'ambiguous' ",
    "WHEN i.f IS NULL OR NOT isfinite(i.f) OR i.f < 0 OR i.f > 1 THEN 'missing_frequency' ",
    "WHEN s.ra IS NULL THEN 'missing_reference' ",
    "WHEN i.allele_a = s.ra AND i.allele_b = s.rb THEN 'direct' ",
    "WHEN i.allele_a = s.rb AND i.allele_b = s.ra THEN 'reversed' ",
    "WHEN translate(i.allele_a, 'ACGT', 'TGCA') = s.ra ",
    "AND translate(i.allele_b, 'ACGT', 'TGCA') = s.rb THEN 'flipped' ",
    "WHEN translate(i.allele_a, 'ACGT', 'TGCA') = s.rb ",
    "AND translate(i.allele_b, 'ACGT', 'TGCA') = s.ra THEN 'flipped_reversed' ",
    "ELSE 'allele_mismatch' END AS match_status FROM raw i LEFT JOIN sites s ",
    "ON i.ref_chr = s.chromosome AND i.position = s.position), ",
    "aligned AS (SELECT sample_id, ref_chr, position, ra, rb, ",
    "CASE WHEN match_status IN ('reversed','flipped_reversed') THEN 1-f ELSE f END AS aligned_f ",
    "FROM classified WHERE match_status IN ('direct','reversed','flipped','flipped_reversed')), ",
    "audit AS (SELECT sample_id, count(*) AS input_variants, ",
    "count(*) FILTER (WHERE match_status IN ('direct','reversed','flipped','flipped_reversed')) AS used_variants, ",
    "count(*) FILTER (WHERE match_status = 'reversed') AS reversed_variants, ",
    "count(*) FILTER (WHERE match_status IN ('flipped','flipped_reversed')) AS flipped_variants, ",
    "count(*) FILTER (WHERE match_status = 'duplicate') AS duplicate_variants, ",
    "count(*) FILTER (WHERE match_status = 'ambiguous') AS ambiguous_variants, ",
    "count(*) FILTER (WHERE match_status IN ('missing_reference','allele_mismatch')) AS unmatched_variants, ",
    "count(*) FILTER (WHERE match_status = 'invalid_alleles') AS invalid_variants, ",
    "count(*) FILTER (WHERE match_status = 'missing_frequency') AS missing_variants ",
    "FROM classified GROUP BY sample_id), ",
    "moments AS (SELECT a.sample_id, ", paste(c(x_terms, y_terms, cor_terms), collapse = ", "),
    " FROM aligned a JOIN ", reference, " r ON a.ref_chr = r.chromosome ",
    "AND a.position = r.position AND a.ra = r.allele_a AND a.rb = r.allele_b ",
    "GROUP BY a.sample_id), ",
    "solved AS (SELECT m.sample_id, duckhts_ancestry_proportions(",
    "list_value(", values_x, "), list_value(", values_y, "), ",
    length(groups), ", ", equality, ") AS q FROM moments m), ",
    "each_cor AS (SELECT m.sample_id, c.group_index, c.cor_each FROM moments m, ",
    "UNNEST(list_value(", values_cor, ")) WITH ORDINALITY AS c(cor_each, group_index)), ",
    "pred_cor AS (SELECT a.sample_id, corr(a.aligned_f, ", predicted, ") AS cor_pred ",
    "FROM aligned a JOIN ", reference, " r ON a.ref_chr = r.chromosome ",
    "AND a.position = r.position AND a.ra = r.allele_a AND a.rb = r.allele_b ",
    "JOIN solved s ON s.sample_id = a.sample_id GROUP BY a.sample_id), ",
    "quality AS (SELECT a.sample_id, p.cor_pred, CASE ",
    "WHEN a.used_variants = 0 THEN 'no_matched_variants' ",
    "WHEN (SELECT avg(e.cor_each) FROM each_cor e WHERE e.sample_id = a.sample_id) < -0.2 ",
    "THEN 'reversed_alleles' WHEN p.cor_pred IS NULL OR NOT isfinite(p.cor_pred) ",
    "OR p.cor_pred < ", gate, " THEN 'low_correlation' ELSE 'ok' END AS status ",
    "FROM audit a LEFT JOIN pred_cor p USING (sample_id)) ",
    "SELECT a.sample_id, g.group_id, CASE WHEN q.status = 'ok' ",
    "THEN s.q[g.group_index] ELSE NULL END AS proportion, ",
    "q.cor_pred, e.cor_each, q.status, a.input_variants, a.used_variants, ",
    "a.input_variants - a.used_variants AS dropped_variants, ",
    "a.reversed_variants, a.flipped_variants, a.duplicate_variants, ",
    "a.ambiguous_variants, a.unmatched_variants, a.invalid_variants, a.missing_variants, ",
    equality, " AS sum_to_one, ", gate, "::DOUBLE AS min_cor ",
    "FROM audit a JOIN quality q USING (sample_id) ",
    "LEFT JOIN solved s ON s.sample_id = a.sample_id ",
    "LEFT JOIN (VALUES ", group_rows, ") g(group_index, group_id) ON s.sample_id IS NOT NULL ",
    "LEFT JOIN each_cor e ON e.sample_id = a.sample_id AND e.group_index = g.group_index ",
    "ORDER BY a.sample_id, g.group_id"
  )
}
