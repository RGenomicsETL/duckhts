#' Select a Somalier Panel from Population Variants
#'
#' Select typed population variants on the caller's DuckDB connection. A table
#' or view may be populated by `read_bcf()`; Parquet and VCF/BCF paths are also
#' supported. Input columns follow `read_bcf()` (`CHROM`, `POS`, `REF`, `ALT`,
#' `FILTER`, and declared `INFO_*` fields). Interval relations have `chrom`,
#' `start`, `stop` in zero-based half-open BED coordinates. A gnotate exclusion
#' relation has `chrom`, `pos`, `ref`, `alt` in one-based VCF coordinates.
#'
#' `somalier_v0.3.4` names the pinned human-contig rules. Upstream skips
#' autosomal REF=C and restricts X to the inclusive zero-based range
#' 2781479..154931044; this selector applies the documented X/Y spacing,
#' whereas the upstream implementation does not insert X/Y sites into its
#' spacing state. Equal AF scores use lexical chromosome, position, and allele
#' ordering, rather than the upstream sort's unspecified tie order. Only
#' uppercase A/C/G/T biallelic SNVs enter the canonical panel. `generic` does
#' not interpret contig names as human sex chromosomes. The result can include
#' X/Y, whereas `duckhts_somalier_import_sites()` admits autosomes only;
#' filter `region` for autosome-only consumers. Gnotate zip input is not
#' supported: use an allele-specific exclusion relation.
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
  table_name = NULL, overwrite = FALSE
) {
  .somalier_validate_output(con, table_name, overwrite)
  .somalier_scalar_text(assembly, "assembly")
  if (nchar(assembly, type = "bytes") > 1024L) {
    stop("assembly must be at most 1024 bytes", call. = FALSE)
  }
  compatibility_mode <- match.arg(compatibility_mode)
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
      source_parquet else source_vcf
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
  source_label = source, diagnostics = FALSE
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
    if (is.null(table)) {
      fields <- paste0("NULL::", types, " AS ", names)
      return(paste0("SELECT ", paste(fields, collapse = ", "), " WHERE false"))
    }
    paste0("SELECT ", paste(sprintf("CAST(%s AS %s) AS %s", names, types, names),
                            collapse = ", "), " FROM ",
           sql_quote_identifier(con, table))
  }
  include <- interval_sql("include", c("chrom", "start", "stop"),
                          c("VARCHAR", "BIGINT", "BIGINT"))
  exclude <- interval_sql("exclude", c("chrom", "start", "stop"),
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
  query <- paste0(
    "WITH src AS MATERIALIZED (SELECT CAST(CHROM AS VARCHAR) AS chrom, ",
    "CAST(POS AS BIGINT) AS pos, CAST(REF AS VARCHAR) AS ref, ",
    "CAST(ALT AS VARCHAR[]) AS alts, CAST(FILTER AS VARCHAR[]) AS filters, ",
    af, " AS af_source, ", an, " AS an, ", paste(annotation, collapse = ", "),
    " FROM ", qsource, "), ",
    "incl AS (", include, "), excl AS (", exclude, "), ",
    "gno AS (", gnotate, "), ",
    "classified AS MATERIALIZED (SELECT *, coalesce(af_source, 0::FLOAT) ",
    "AS af, ", sex, " AS sex FROM src), ",
    "eligible AS MATERIALIZED (SELECT * FROM classified s WHERE ",
    "pos > 0 AND ref IS NOT NULL AND alts IS NOT NULL AND len(alts) > 0 ",
    if (human) "AND (sex != 'autosome' OR ref != 'C') " else "",
    if (human) "AND (sex != 'X' OR pos - 1 BETWEEN 2781479 AND 154931044) " else "",
    "), ",
    "indels AS MATERIALIZED (SELECT chrom, greatest(0, pos - 8) AS lo, ",
    "pos - 1 + len(ref) + 7 AS hi FROM eligible WHERE ",
    "(len(ref) != 1 OR len(alts) != 1 OR len(alts[1]) != 1) AND af > 0.02), ",
    "snps AS MATERIALIZED (SELECT chrom, pos FROM eligible WHERE ",
    "len(ref) = 1 AND len(alts) = 1 AND len(alts[1]) = 1 AND ",
    "filters = ['PASS'] AND af > 0.01), ",
    "gated AS MATERIALIZED (SELECT * FROM eligible s WHERE ",
    "len(ref) = 1 AND len(alts) = 1 AND len(alts[1]) = 1 ",
    "AND ref IN ('A','C','G','T') AND alts[1] IN ('A','C','G','T') ",
    "AND ref != alts[1] AND filters = ['PASS'] ",
    "AND (an IS NULL OR sex != 'autosome' OR an >= ", min_an, ") ",
    "AND isfinite(af) AND ",
    "CAST(af AS DOUBLE) BETWEEN CASE WHEN sex = 'autosome' THEN ", af_literal,
    " ELSE 0.04 END AND CASE WHEN sex = 'autosome' THEN 1 - ",
    af_literal, " ELSE 0.96 END ",
    "AND (as_status IS NULL OR as_status = 'PASS') ",
    "AND old_multiallelic IS NULL AND old_variant IS NULL ",
    "AND (sex = 'Y' OR NOT coalesce(segdup, false)) ",
    "AND NOT coalesce(lcr, false) ",
    "AND (sex != 'autosome' OR (",
    paste(sprintf("(%s IS NULL OR abs(%s) <= 2.4)",
                  c("BaseQRankSum", "MQRankSum", "ClippingRankSum",
                    "ReadPosRankSum", "FS"),
                  c("BaseQRankSum", "MQRankSum", "ClippingRankSum",
                    "ReadPosRankSum", "FS")), collapse = " AND "),
    " AND (QD IS NULL OR abs(QD) >= 12) AND (MQ IS NULL OR MQ >= 50))) ",
    "AND NOT EXISTS (SELECT 1 FROM excl e WHERE e.chrom = s.chrom ",
    "AND e.start < s.pos + 5 AND e.stop > greatest(0, s.pos - 6)) ",
    if (!is.null(intervals$include)) paste0(
      "AND EXISTS (SELECT 1 FROM incl i WHERE i.chrom = s.chrom ",
      "AND i.start < s.pos AND i.stop > s.pos - 1) "
    ) else "",
    "AND (s.sex = 'X' OR NOT EXISTS (SELECT 1 FROM gno g WHERE ",
    "g.chrom = s.chrom AND g.pos = s.pos AND g.ref = s.ref ",
    "AND g.alt = s.alts[1]))), ",
    "neighbors AS MATERIALIZED (SELECT * FROM gated s WHERE ",
    "NOT EXISTS (SELECT 1 FROM indels i WHERE i.chrom = s.chrom ",
    "AND i.lo < s.pos + 1 AND i.hi > greatest(0, s.pos - 2)) ",
    "AND (SELECT count(*) FROM snps n WHERE n.chrom = s.chrom ",
    "AND n.pos BETWEEN s.pos - 2 AND s.pos + 2) <= 1), ",
    "ranked AS MATERIALIZED (SELECT *, ",
    "abs(CAST(af AS DOUBLE) - CAST(", target, " AS DOUBLE)) AS selection_score, ",
    "row_number() OVER (PARTITION BY chrom ORDER BY selection_score, ",
    "pos, ref, alts[1], af, an) AS chrom_rank FROM neighbors), ",
    "spaced AS (SELECT chrom, generate_subscripts(mask, 1) AS chrom_rank, ",
    "unnest(mask) AS selected FROM (SELECT chrom, ",
    "duckhts_somalier_spacing(list(CAST(pos AS UBIGINT) ORDER BY chrom_rank), ",
    "CAST(max(CASE WHEN sex = 'X' THEN 1000 WHEN sex = 'Y' THEN 200 ",
    "ELSE ", snp_dist, " END) AS UBIGINT)) AS mask FROM ranked GROUP BY chrom)), ",
    "kept AS (SELECT r.* FROM ranked r JOIN spaced p USING(chrom, chrom_rank) ",
    "WHERE p.selected), ",
    "capped AS (SELECT *, row_number() OVER (PARTITION BY sex ",
    "ORDER BY selection_score, chrom, pos, ref, alts[1], af, an) AS cap_rank ",
    "FROM kept) "
  )
  if (diagnostics) {
    stages <- c("src", "eligible", "indels", "snps", "gated",
                "neighbors", "kept")
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
