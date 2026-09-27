#' Select a Somalier Panel from Population Variants
#'
#' Select typed population variants on the caller's DuckDB connection. A table
#' or view may be populated by `read_bcf()`; Parquet and VCF/BCF paths are also
#' supported. Input columns follow `read_bcf()` (`CHROM`, `POS`, `REF`, `ALT`,
#' `FILTER`, and declared `INFO_*` fields). Interval relations have `chrom`,
#' `start`, `end` in zero-based half-open BED coordinates, as returned by
#' `read_bed()`. A gnotate exclusion
#' relation has `chrom`, `pos`, `ref`, `alt` in one-based VCF coordinates.
#'
#' `somalier_v0.3.4` names the pinned human-contig rules. Upstream skips
#' autosomal REF=C and restricts X to the inclusive zero-based range
#' 2781479..154931044. By default X/Y selections do not enter the spacing
#' state, as in upstream v0.3.4. `sex_spacing = "enforced"` applies the X/Y
#' distances to every selection. Upstream's stable AF-score sort retains input
#' order on ties; this selector uses scan order by default, or lexical position
#' and allele order with `tie_order = "lexical"`. For a table or view, scan order
#' is only reliable when its query explicitly orders the records. Only
#' uppercase A/C/G/T biallelic SNVs enter the canonical panel. `generic` does
#' not interpret contig names as human sex chromosomes. The result can include
#' X/Y, whereas `duckhts_somalier_import_sites()` admits autosomes only;
#' filter `region` for autosome-only consumers. Gnotate zip input is not
#' supported: use an allele-specific exclusion relation. `assembly` labels
#' the input coordinates; this function does not run the liftover steps in
#' Somalier's `scripts/find_sites.sh`.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param source_table A typed source table or view, including a view over read_bcf().
#' @param source_parquet Parquet source path instead of source_table.
#' @param source_vcf VCF/BCF path instead of source_table.
#' @param assembly Nonempty assembly identifier.
#' @param min_af,min_an,af_field,an_field,snp_dist,target_af Somalier find-sites
#'   thresholds and INFO field names; defaults match Somalier v0.3.4.
#' @param include_table,exclude_table BED interval relations; NULL omits the gate.
#' @param gnotate_exclude_table Optional allele-specific exclusion relation in
#'   place of Somalier's gnotate zip.
#' @param compatibility_mode Named contig and REF=C policy.
#' @param sex_spacing Whether to reproduce upstream's unspaced X/Y selection
#'   or enforce X/Y distances.
#' @param tie_order Retain input scan order or use lexical position and allele
#'   order for equal AF scores.
#' @param max_autosomal,max_x,max_y Site caps after greedy selection.
#' @param table_name Optional output table, otherwise return a data frame.
#' @param overwrite Whether to replace an existing output table.
#' @return A canonical panel relation with selection score, source AN and
#'   spacing-distance diagnostics.
#' @export
rduckhts_somalier_find_sites <- function(
  con, source_table = NULL, source_parquet = NULL, source_vcf = NULL,
  assembly, min_af = 0.15, min_an = 115000, af_field = "AF", an_field = "AN",
  snp_dist = 10000, target_af = 0.48, include_table = NULL,
  exclude_table = NULL, gnotate_exclude_table = NULL,
  compatibility_mode = c("somalier_v0.3.4", "generic"),
  max_autosomal = 65535, max_x = 10001, max_y = 5001,
  table_name = NULL, overwrite = FALSE,
  sex_spacing = c("upstream", "enforced"), tie_order = c("input", "lexical")
) {
  .somalier_validate_output(con, table_name, overwrite)
  .somalier_scalar_text(assembly, "assembly")
  if (nchar(assembly, type = "bytes") > 1024L) {
    stop("assembly must be at most 1024 bytes", call. = FALSE)
  }
  compatibility_mode <- match.arg(compatibility_mode)
  sex_spacing <- match.arg(sex_spacing)
  tie_order <- match.arg(tie_order)
  min_af <- .somalier_fraction(min_af, "min_af")
  target_af <- .somalier_fraction(target_af, "target_af")
  min_an <- .somalier_bounded_whole_number(min_an, "min_an", 0, 2147483647)
  snp_dist <- .somalier_bounded_whole_number(snp_dist, "snp_dist", 1, 2147483647)
  max_autosomal <- .somalier_bounded_whole_number(
    max_autosomal, "max_autosomal", 1, 1000000
  )
  max_x <- .somalier_bounded_whole_number(max_x, "max_x", 1, 1000000)
  max_y <- .somalier_bounded_whole_number(max_y, "max_y", 1, 1000000)
  if (length(af_field) != 1L || length(an_field) != 1L) {
    stop("INFO fields must be scalar", call. = FALSE)
  }
  for (field in c(af_field, an_field)) .somalier_scalar_text(field, "INFO field")
  if (any(nchar(c(af_field, an_field), type = "bytes") > 256L)) {
    stop("INFO field names must be at most 256 bytes", call. = FALSE)
  }
  sources <- sum(!vapply(list(source_table, source_parquet, source_vcf),
                        is.null, logical(1)))
  if (sources != 1L) {
    stop("supply exactly one population source", call. = FALSE)
  }
  if (!is.null(source_vcf)) {
    .somalier_scalar_text(source_vcf, "source_vcf")
    view <- .somalier_temp_view_name()
    DBI::dbExecute(con, sprintf(
      "CREATE TEMP VIEW %s AS SELECT * FROM read_bcf(%s, samples := '', scan_mode := 'sequential')",
      sql_quote_identifier(con, view), sql_quote_string(con, source_vcf)
    ))
    source <- list(name = view, temporary = TRUE)
  } else {
    source <- .somalier_source(con, source_table, source_parquet, "population")
  }
  if (source$temporary) {
    on.exit(.somalier_drop_view(con, source$name), add = TRUE, after = FALSE)
  }
  interval_sources <- list(include = include_table, exclude = exclude_table,
                           gnotate = gnotate_exclude_table)
  for (name in names(interval_sources)) {
    value <- interval_sources[[name]]
    if (!is.null(value)) {
      .somalier_validate_name(value, paste0(name, "_table"))
    }
  }
  query <- .somalier_find_sites_query(
    con, source$name, assembly, min_af, min_an, af_field, an_field,
    snp_dist, target_af, interval_sources, compatibility_mode,
    max_autosomal, max_x, max_y,
    if (!is.null(source_table)) source_table else if (!is.null(source_parquet))
      source_parquet else source_vcf, sex_spacing = sex_spacing,
    tie_order = tie_order
  )
  .somalier_publish_query(con, query, table_name, overwrite)
}

.somalier_find_sites_field <- function(con, columns, name, result_type) {
  column <- paste0("INFO_", name)
  if (!(column %in% columns$column_name)) return(paste0("NULL::", result_type))
  type <- columns$column_type[match(column, columns$column_name)]
  expression <- sql_quote_identifier(con, column)
  if (grepl("\\[\\]$", type)) expression <- paste0(expression, "[1]")
  paste0("TRY_CAST(", expression, " AS ", result_type, ")")
}

.somalier_find_sites_query <- function(
  con, source, assembly, min_af, min_an, af_field, an_field,
  snp_dist, target_af, intervals, mode, max_autosomal, max_x, max_y,
  source_label = source, diagnostics = FALSE,
  sex_spacing = "upstream", tie_order = "input"
) {
  quote <- function(value) sql_quote_string(con, value)
  qsource <- sql_quote_identifier(con, source)
  columns <- DBI::dbGetQuery(con, paste0("DESCRIBE SELECT * FROM ", qsource))
  required <- c("CHROM", "POS", "REF", "ALT", "FILTER")
  if (!all(required %in% columns$column_name)) {
    stop("population relation needs CHROM, POS, REF, ALT, FILTER", call. = FALSE)
  }
  field <- function(name, type) .somalier_find_sites_field(
    con, columns, name, type
  )
  af <- field(af_field, "FLOAT")
  an <- field(an_field, "BIGINT")
  old_annotation <- function(name) {
    column <- paste0("INFO_", name)
    type <- columns$column_type[match(column, columns$column_name)]
    if (length(type) && !is.na(type) && type == "BOOLEAN") {
      return(paste0("CASE WHEN ", sql_quote_identifier(con, column),
                    " THEN 'present' END"))
    }
    field(name, "VARCHAR")
  }
  annotation <- c(
    paste0(field("AS_FilterStatus", "VARCHAR"), " AS as_status"),
    paste0(old_annotation("OLD_MULTIALLELIC"), " AS old_multiallelic"),
    paste0(old_annotation("OLD_VARIANT"), " AS old_variant"),
    paste0(field("segdup", "BOOLEAN"), " AS segdup"),
    paste0(field("lcr", "BOOLEAN"), " AS lcr"),
    vapply(c("BaseQRankSum", "MQRankSum", "ClippingRankSum",
             "ReadPosRankSum", "FS", "QD", "MQ"), function(name) {
      paste0(field(name, "DOUBLE"), " AS ", sql_quote_identifier(con, name))
    }, character(1))
  )
  interval_sql <- function(name, names, types) {
    table <- intervals[[name]]
    quoted <- vapply(names, function(x) sql_quote_identifier(con, x), character(1))
    if (is.null(table)) {
      fields <- paste0("NULL::", types, " AS ", quoted)
      return(paste0("SELECT ", paste(fields, collapse = ", "), " WHERE false"))
    }
    paste0("SELECT ", paste(sprintf("CAST(%s AS %s) AS %s", quoted, types, quoted),
                            collapse = ", "), " FROM ",
           sql_quote_identifier(con, table))
  }
  include <- interval_sql("include", c("chrom", "start", "end"),
                          c("VARCHAR", "BIGINT", "BIGINT"))
  exclude <- interval_sql("exclude", c("chrom", "start", "end"),
                          c("VARCHAR", "BIGINT", "BIGINT"))
  gnotate <- interval_sql("gnotate", c("chrom", "pos", "ref", "alt"),
                          c("VARCHAR", "BIGINT", "VARCHAR", "VARCHAR"))
  human <- identical(mode, "somalier_v0.3.4")
  x <- "('X','chrX','NC_000023.10','NC_000023.11')"
  y <- "('Y','chrY','NC_000024.9','NC_000024.10')"
  sex <- if (human) paste0(
    "CASE WHEN chrom IN ", x, " THEN 'X' WHEN chrom IN ", y,
    " THEN 'Y' ELSE 'autosome' END"
  ) else "'autosome'"
  af_literal <- format(min_af, digits = 17, scientific = FALSE)
  target <- format(target_af, digits = 17, scientific = FALSE)
  rank_tie <- if (tie_order == "input") "source_ordinal" else
    "pos, ref, alts[1], af, an, source_ordinal"
  sex_mask <- if (sex_spacing == "upstream" && human) paste0(
    "CASE WHEN max(sex) IN ('X','Y') THEN ",
    "list_transform(list(CAST(pos AS UBIGINT) ORDER BY chrom_rank), ",
    "lambda p: true) ELSE "
  ) else ""
  snv_pred <- paste0(
    "len(ref) = 1 AND len(alts) = 1 AND len(alts[1]) = 1 ",
    "AND ref IN ('A','C','G','T') AND alts[1] IN ('A','C','G','T') ",
    "AND ref != alts[1] AND filters = ['PASS']"
  )
  af_an_pred <- paste0(
    "(an IS NULL OR sex != 'autosome' OR an >= ", min_an, ") ",
    "AND isfinite(af) AND ",
    "CAST(af AS DOUBLE) BETWEEN CASE WHEN sex = 'autosome' THEN ", af_literal,
    " ELSE 0.04 END AND CASE WHEN sex = 'autosome' THEN 1 - ",
    af_literal, " ELSE 0.96 END"
  )
  indel_pred <- paste0(
    "(len(ref) != 1 OR len(alts) != 1 OR len(alts[1]) != 1) AND af > 0.02"
  )
  snp_pred <- paste0(
    "len(ref) = 1 AND len(alts) = 1 AND len(alts[1]) = 1 ",
    "AND filters = ['PASS'] AND af > 0.01"
  )
  event_filter <- if (diagnostics) "" else paste0(
    " WHERE (", indel_pred, ") OR (", snp_pred, ") OR ((",
    snv_pred, ") AND (", af_an_pred, "))"
  )
  inventory_source <- if (diagnostics) "eligible" else "events"
  materialization <- if (diagnostics) "MATERIALIZED" else "NOT MATERIALIZED"
  query <- paste0(
    "WITH src AS ", materialization, " (SELECT ",
    "row_number() OVER () AS source_ordinal, ",
    "CAST(CHROM AS VARCHAR) AS chrom, ",
    "CAST(POS AS BIGINT) AS pos, CAST(REF AS VARCHAR) AS ref, ",
    "CAST(ALT AS VARCHAR[]) AS alts, CAST(FILTER AS VARCHAR[]) AS filters, ",
    af, " AS af_source, ", an, " AS an, ", paste(annotation, collapse = ", "),
    " FROM ", qsource, "), ",
    "incl AS (", include, "), excl AS (", exclude, "), ",
    "gno AS (", gnotate, "), ",
    "classified AS NOT MATERIALIZED (SELECT *, coalesce(af_source, 0::FLOAT) ",
    "AS af, ", sex, " AS sex FROM src), ",
    "eligible AS ", materialization,
    " (SELECT * FROM classified s WHERE ",
    "pos > 0 AND ref IS NOT NULL AND alts IS NOT NULL AND len(alts) > 0 ",
    if (human) "AND (sex != 'autosome' OR ref != 'C') " else "",
    if (human) "AND (sex != 'X' OR pos - 1 BETWEEN 2781479 AND 154931044) " else "",
    "), ",
    "events AS MATERIALIZED (SELECT * FROM eligible", event_filter, "), ",
    "indels AS MATERIALIZED (SELECT chrom, greatest(0, pos - 8) AS lo, ",
    "pos - 1 + len(ref) + 7 AS hi FROM ", inventory_source, " WHERE ",
    indel_pred, "), ",
    "snps AS MATERIALIZED (SELECT chrom, pos FROM ", inventory_source,
    " WHERE ", snp_pred, "), ",
    "snp_positions AS (SELECT chrom, pos, count(*) AS n FROM snps ",
    "GROUP BY chrom, pos), ",
    "snp_neighbors AS (SELECT chrom, pos, sum(n) OVER (PARTITION BY chrom ",
    "ORDER BY pos RANGE BETWEEN 2 PRECEDING AND 2 FOLLOWING) AS n ",
    "FROM snp_positions), ",
    "pass_snv AS (SELECT * FROM ", inventory_source, " WHERE ",
    snv_pred, "), ",
    "af_an AS (SELECT * FROM pass_snv WHERE ", af_an_pred, "), ",
    "annotation_gate AS (SELECT * FROM af_an WHERE ",
    "(as_status IS NULL OR as_status = 'PASS') ",
    "AND old_multiallelic IS NULL AND old_variant IS NULL ",
    "AND (sex = 'Y' OR NOT coalesce(segdup, false)) ",
    "AND NOT coalesce(lcr, false)), ",
    # Somalier v0.3.4 findsites.nim:214-224 takes abs() of every QC statistic,
    # QD included (rejects abs(QD) < 12), so QD is symmetric here too.
    "qc_gate AS (SELECT * FROM annotation_gate WHERE ",
    "sex != 'autosome' OR (",
    paste(sprintf("(%s IS NULL OR abs(%s) <= 2.4)",
                  c("BaseQRankSum", "MQRankSum", "ClippingRankSum",
                    "ReadPosRankSum", "FS"),
                  c("BaseQRankSum", "MQRankSum", "ClippingRankSum",
                    "ReadPosRankSum", "FS")), collapse = " AND "),
    " AND (QD IS NULL OR abs(QD) >= 12) AND (MQ IS NULL OR MQ >= 50))), ",
    "interval_gate AS (SELECT * FROM qc_gate s WHERE ",
    "NOT EXISTS (SELECT 1 FROM excl e WHERE e.chrom = s.chrom ",
    "AND e.start < s.pos + 5 AND e.\"end\" > greatest(0, s.pos - 6)) ",
    if (!is.null(intervals$include)) paste0(
      "AND EXISTS (SELECT 1 FROM incl i WHERE i.chrom = s.chrom ",
      "AND i.start < s.pos AND i.\"end\" > s.pos - 1) "
    ) else "",
    "), gated AS MATERIALIZED (SELECT * FROM interval_gate s WHERE ",
    "s.sex = 'X' OR NOT EXISTS (SELECT 1 FROM gno g WHERE ",
    "g.chrom = s.chrom AND g.pos = s.pos AND g.ref = s.ref ",
    "AND g.alt = s.alts[1])), ",
    "indel_events AS (SELECT chrom, lo AS pos, hi, 0 AS kind, ",
    "NULL::BIGINT AS source_ordinal FROM indels UNION ALL ",
    "SELECT chrom, pos, NULL::BIGINT, 1, source_ordinal FROM gated), ",
    "indel_coverage AS (SELECT source_ordinal, reach_hi FROM (SELECT ",
    "source_ordinal, kind, max(hi) OVER (PARTITION BY chrom ORDER BY ",
    "pos, kind ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW) AS reach_hi ",
    "FROM indel_events) WHERE kind = 1), ",
    "indel_clear AS MATERIALIZED (SELECT s.* FROM gated s ",
    "JOIN indel_coverage i USING (source_ordinal) WHERE i.reach_hi IS NULL ",
    "OR i.reach_hi <= greatest(0, s.pos - 2)), ",
    "neighbors AS MATERIALIZED (SELECT s.* FROM indel_clear s ",
    "LEFT JOIN snp_neighbors n ON n.chrom = s.chrom AND n.pos = s.pos ",
    "WHERE coalesce(n.n, 0) <= 1), ",
    "ranked AS MATERIALIZED (SELECT *, ",
    "abs(CAST(af AS DOUBLE) - CAST(", target, " AS DOUBLE)) AS selection_score, ",
    "row_number() OVER (PARTITION BY chrom ORDER BY selection_score, ",
    rank_tie, ") AS chrom_rank FROM neighbors), ",
    "spaced AS (SELECT chrom, generate_subscripts(mask, 1) AS chrom_rank, ",
    "unnest(mask) AS selected FROM (SELECT chrom, ", sex_mask,
    "duckhts_somalier_spacing(list(CAST(pos AS UBIGINT) ORDER BY chrom_rank), ",
    "CAST(max(CASE WHEN sex = 'X' THEN 1000 WHEN sex = 'Y' THEN 200 ",
    "ELSE ", snp_dist, " END) AS UBIGINT)) ",
    if (nzchar(sex_mask)) "END " else "", "AS mask FROM ranked GROUP BY chrom)), ",
    "kept AS (SELECT r.* FROM ranked r JOIN spaced p USING(chrom, chrom_rank) ",
    "WHERE p.selected), ",
    "capped AS (SELECT *, row_number() OVER (PARTITION BY sex ",
    "ORDER BY selection_score, ",
    if (tie_order == "input") "source_ordinal" else
      "chrom, pos, ref, alts[1], af, an, source_ordinal",
    ") AS cap_rank ",
    "FROM kept) "
  )
  if (diagnostics) {
    stages <- c("src", "eligible", "pass_snv", "af_an",
                "annotation_gate", "qc_gate", "interval_gate", "gated",
                "indels", "snps", "indel_clear", "neighbors", "kept")
    counts <- paste0("SELECT '", stages, "' AS gate, count(*) AS records FROM ",
                     stages)
    counts <- c(counts, paste0(
      "SELECT 'selected' AS gate, count(*) AS records FROM capped WHERE ",
      "cap_rank <= CASE WHEN sex = 'X' THEN ", max_x,
      " WHEN sex = 'Y' THEN ", max_y, " ELSE ", max_autosomal, " END"
    ))
    return(paste0(query, paste(counts, collapse = " UNION ALL ")))
  }
  paste0(query, "SELECT ", quote(assembly), " AS assembly, ",
    "CAST(row_number() OVER (ORDER BY chrom, pos, ref, alts[1]) - 1 ",
    "AS UBIGINT) AS site_index, chrom AS region, CAST(pos AS UBIGINT) AS position, ",
    "least(ref, alts[1]) AS allele_a, greatest(ref, alts[1]) AS allele_b, ",
    "CAST(CASE WHEN ref < alts[1] THEN af ELSE 1 - af END AS DOUBLE) ",
    "AS population_b_af, ", quote(source_label), " AS source_path, ",
    "ref AS source_ref, alts[1] AS source_alt, ",
    "CAST(af_source AS DOUBLE) AS source_alt_af, filters AS source_filter, ",
    "an AS source_an, CAST(selection_score AS DOUBLE) AS selection_score, ",
    "CASE WHEN sex = 'X' THEN 1000 WHEN sex = 'Y' THEN 200 ELSE ",
    snp_dist, " END AS spacing_distance FROM capped WHERE cap_rank <= ",
    "CASE WHEN sex = 'X' THEN ", max_x, " WHEN sex = 'Y' THEN ",
    max_y, " ELSE ", max_autosomal, " END ORDER BY region, position"
  )
}
