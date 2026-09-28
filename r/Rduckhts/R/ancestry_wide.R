#' Infer ancestry proportions from a keyed wide reference relation.
#'
#' The reference relation has `chromosome` (integer), `position` (one-based),
#' `allele_a`, `allele_b`, `PC1` through `PC16`, and one DOUBLE frequency column
#' per reference group. It can be a view over the Parquet product staged by
#' `duckhtsbench::duckhts_bench_stage_ancestry_parquet()`. Each site must have
#' exactly one reference row. The input and correction schemas, audit statuses,
#' and QP solver contract are those of [rduckhts_ancestry_proportions()].
#' Allele-aligned input frequencies are written to a temporary Parquet file in
#' `tempdir()` and removed before the function returns.
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
  quote_id <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  quote_str <- function(x) as.character(DBI::dbQuoteString(con, x))
  audit <- basename(tempfile("ancestry_audit_"))
  moments <- basename(tempfile("ancestry_moments_"))
  solved <- basename(tempfile("ancestry_solved_"))
  prediction <- basename(tempfile("ancestry_prediction_"))
  aligned <- tempfile("ancestry_aligned_", fileext = ".parquet")
  on.exit({
    unlink(aligned)
    DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ", quote_id(prediction)))
    DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ", quote_id(solved)))
    DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ", quote_id(moments)))
    DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ", quote_id(audit)))
  }, add = TRUE)
  classification <- .ancestry_classification_query(relations, frequency)
  DBI::dbExecute(con, paste0("CREATE TEMP TABLE ", quote_id(audit),
    " AS ", .ancestry_indexed_audit_query(classification)))
  DBI::dbExecute(con, paste0("COPY (", .ancestry_aligned_query(classification),
    ") TO ", quote_str(aligned),
    " (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 32768)"))
  aligned_source <- paste0("read_parquet(", quote_str(aligned), ")")
  one_sample <- DBI::dbGetQuery(con, paste0("SELECT count(*) AS n FROM ",
                                         quote_id(audit)))$n == 1L
  indexed_source <- if (one_sample) aligned_source else paste0(
    "(SELECT a.*, i.sample_index FROM ", aligned_source, " a JOIN ",
    quote_id(audit), " i USING (sample_id))")
  DBI::dbExecute(con, paste0("CREATE TEMP TABLE ", quote_id(moments), " AS ",
    .ancestry_moments_query(con, relations[[2L]], groups, corrections$coefficient,
                            indexed_source, one_sample)))
  DBI::dbExecute(con, paste0("CREATE TEMP TABLE ", quote_id(solved), " AS SELECT ",
    "m.sample_id, ", .ancestry_solver_expression(groups, sum_to_one),
    " AS q FROM ", quote_id(moments), " m"))
  single_proportions <- NULL
  if (one_sample) {
    proportions <- DBI::dbGetQuery(con, paste0("SELECT q FROM ", quote_id(solved)))
    if (nrow(proportions)) {
      single_proportions <- proportions$q[[1L]]
    }
  }
  DBI::dbExecute(con, paste0("CREATE TEMP TABLE ", quote_id(prediction), " AS ",
    .ancestry_prediction_query(con, relations[[2L]], groups, aligned_source,
                               quote_id(solved), one_sample, single_proportions)))
  query <- .ancestry_wide_query(con, groups, quote_id(audit), quote_id(moments),
                                quote_id(solved), quote_id(prediction), sum_to_one, min_cor)
  .somalier_publish_query(con, query, table_name, overwrite)
}

.ancestry_classification_query <- function(relations, frequency) {
  input <- relations[[1L]]
  reference <- relations[[2L]]
  paste0(
    "raw AS (SELECT sample_id, chromosome, position, allele_a, allele_b, ",
    "try_cast(regexp_replace(chromosome, '^chr', '') AS INTEGER) AS ref_chr, ",
    frequency, " AS f, count(*) OVER (PARTITION BY sample_id, chromosome, position) AS copies ",
    "FROM ", input, "), ",
    "sites AS (SELECT chromosome, position, allele_a AS ra, allele_b AS rb ",
    "FROM ", reference, "), ",
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
    "ON i.ref_chr = s.chromosome AND i.position = s.position)")
}

.ancestry_aligned_query <- function(classification) {
  paste0("WITH ", classification,
    " SELECT sample_id, ref_chr AS chromosome, position, ra AS allele_a, rb AS allele_b, ",
    "CASE WHEN match_status IN ('reversed','flipped_reversed') THEN 1-f ELSE f END ",
    "AS aligned_f FROM classified WHERE match_status IN ",
    "('direct','reversed','flipped','flipped_reversed')")
}

.ancestry_audit_query <- function() {
  paste0(
    "SELECT sample_id, count(*) AS input_variants, ",
    "count(*) FILTER (WHERE match_status IN ('direct','reversed','flipped','flipped_reversed')) AS used_variants, ",
    "count(*) FILTER (WHERE match_status = 'reversed') AS reversed_variants, ",
    "count(*) FILTER (WHERE match_status IN ('flipped','flipped_reversed')) AS flipped_variants, ",
    "count(*) FILTER (WHERE match_status = 'duplicate') AS duplicate_variants, ",
    "count(*) FILTER (WHERE match_status = 'ambiguous') AS ambiguous_variants, ",
    "count(*) FILTER (WHERE match_status IN ('missing_reference','allele_mismatch')) AS unmatched_variants, ",
    "count(*) FILTER (WHERE match_status = 'invalid_alleles') AS invalid_variants, ",
    "count(*) FILTER (WHERE match_status = 'missing_frequency') AS missing_variants ",
    "FROM classified GROUP BY sample_id")
}

.ancestry_indexed_audit_query <- function(classification) {
  paste0("WITH ", classification, ", totals AS (", .ancestry_audit_query(),
         ") SELECT totals.*, row_number() OVER (ORDER BY sample_id)::INTEGER ",
         "AS sample_index FROM totals")
}

.ancestry_moments_query <- function(con, reference, groups, coefficients, aligned,
                                    one_sample) {
  quote_id <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  pc <- paste0("r.", quote_id(paste0("PC", seq_along(coefficients))))
  group <- paste0("r.", quote_id(groups))
  # The solver reads X in PC-major, group-minor order.
  x_terms <- unlist(lapply(seq_along(pc), function(k) {
    paste0("sum(", pc[[k]], " * ", group, ") AS ",
           paste0("x", k, "_", seq_along(group)))
  }), use.names = FALSE)
  y_terms <- paste0("sum(", pc, " * a.aligned_f) * ",
                    vapply(coefficients, function(x) .somalier_quote_number(con, x),
                           character(1L)), " AS y", seq_along(pc))
  cor_terms <- paste0("corr(a.aligned_f, ", group, ") AS cor", seq_along(group))
  sample <- "any_value(a.sample_id) AS sample_id"
  aggregate <- if (one_sample) "HAVING count(*) > 0" else "GROUP BY a.sample_index"
  paste0("SELECT ", sample, ", ", paste(c(x_terms, y_terms, cor_terms), collapse = ", "),
         " FROM ", aligned, " a JOIN ", reference, " r ON a.chromosome = r.chromosome ",
         "AND a.position = r.position AND a.allele_a = r.allele_a ",
         "AND a.allele_b = r.allele_b ", aggregate)
}

.ancestry_solver_expression <- function(groups, sum_to_one) {
  x_names <- unlist(lapply(seq_len(16L), function(k) {
    paste0("m.x", k, "_", seq_along(groups))
  }), use.names = FALSE)
  paste0("duckhts_ancestry_proportions(list_value(", paste(x_names, collapse = ", "),
         "), list_value(", paste0("m.y", seq_len(16L), collapse = ", "), "), ",
         length(groups), ", ", if (sum_to_one) "true" else "false", ")")
}

.ancestry_prediction_query <- function(con, reference, groups, aligned,
                                       solved, one_sample, single_proportions) {
  quote_id <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  group <- paste0("r.", quote_id(groups))
  if (one_sample && length(single_proportions)) {
    weights <- vapply(single_proportions, function(x) .somalier_quote_number(con, x),
                      character(1L))
  } else if (one_sample) {
    weights <- rep("NULL::DOUBLE", length(groups))
  } else {
    weights <- paste0("s.q[", seq_along(groups), "]")
  }
  predicted <- paste0("(", group, " * ", weights, ")", collapse = " + ")
  sample <- if (one_sample) "any_value(a.sample_id)" else "a.sample_id"
  aggregate <- if (one_sample) "HAVING count(*) > 0" else "GROUP BY a.sample_id"
  solved_join <- if (one_sample) "" else paste0("JOIN ", solved,
                                               " s ON s.sample_id = a.sample_id ")
  paste0("SELECT ", sample, " AS sample_id, corr(a.aligned_f, ", predicted,
    ") AS cor_pred FROM ", aligned, " a JOIN ", reference,
    " r USING (chromosome, position, allele_a, allele_b) ", solved_join, aggregate)
}

.ancestry_wide_query <- function(con, groups, audit, moments, solved, prediction,
                                 sum_to_one, min_cor) {
  quote_str <- function(x) as.character(DBI::dbQuoteString(con, x))
  values_cor <- paste0("m.cor", seq_along(groups), collapse = ", ")
  group_rows <- paste0("(", seq_along(groups), ", ",
                       vapply(groups, quote_str, character(1L)), ")", collapse = ", ")
  equality <- if (sum_to_one) "true" else "false"
  gate <- .somalier_quote_number(con, min_cor)
  paste0(
    "WITH each_cor AS (SELECT m.sample_id, c.group_index, c.cor_each FROM ", moments,
    " m, UNNEST(list_value(", values_cor,
    ")) WITH ORDINALITY AS c(cor_each, group_index)), ",
    "quality AS (SELECT a.sample_id, p.cor_pred, CASE ",
    "WHEN a.used_variants = 0 THEN 'no_matched_variants' ",
    "WHEN (SELECT avg(e.cor_each) FROM each_cor e WHERE e.sample_id = a.sample_id) < -0.2 ",
    "THEN 'reversed_alleles' WHEN p.cor_pred IS NULL OR NOT isfinite(p.cor_pred) ",
    "OR p.cor_pred < ", gate, " THEN 'low_correlation' ELSE 'ok' END AS status ",
    "FROM ", audit, " a LEFT JOIN ", prediction, " p USING (sample_id)) ",
    "SELECT a.sample_id, g.group_id, CASE WHEN q.status = 'ok' ",
    "THEN s.q[g.group_index] ELSE NULL END AS proportion, ",
    "q.cor_pred, e.cor_each, q.status, a.input_variants, a.used_variants, ",
    "a.input_variants - a.used_variants AS dropped_variants, ",
    "a.reversed_variants, a.flipped_variants, a.duplicate_variants, ",
    "a.ambiguous_variants, a.unmatched_variants, a.invalid_variants, a.missing_variants, ",
    equality, " AS sum_to_one, ", gate, "::DOUBLE AS min_cor ",
    "FROM ", audit, " a JOIN quality q USING (sample_id) ",
    "LEFT JOIN ", solved, " s ON s.sample_id = a.sample_id ",
    "LEFT JOIN (VALUES ", group_rows, ") g(group_index, group_id) ON s.sample_id IS NOT NULL ",
    "LEFT JOIN each_cor e ON e.sample_id = a.sample_id AND e.group_index = g.group_index ",
    "ORDER BY a.sample_id, g.group_id"
  )
}
