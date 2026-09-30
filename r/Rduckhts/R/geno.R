#' Read Record-Major Genotypes
#'
#' Read typed GT/PS calls without repeating variant text per sample. Each
#' record has a zero-based scan-local `record_index` and a list of calls with
#' original-header sample indices, nullable allele indices, per-slot phase bits,
#' and nullable scalar phase sets. An absent GT has NULL allele/phase lists;
#' a missing allele still occupies a slot. Phase bits follow HTSlib decoding,
#' including its leading-slot convention for VCF versions before 4.4.
#' PS cardinality excludes vector-end padding retained after sample selection.
#'
#' Full scans preserve the input stream when assigning ordinals. Indexed regions
#' use HTSlib's union order and start a new ordinal at zero. Use SQL `ORDER BY
#' record_index` when order matters; the ordinal is not a persistent file locator.
#' Empty sample selection and sparse calls preserve every selected variant row.
#'
#' @inheritParams rduckhts_bcf
#' @param table_name Optional table to create; `NULL` returns a data frame.
#' @param non_reference_only Omit calls without any called alternate allele.
#'   This does not infer phase or remove variant records.
#' @param format_fields Character vector of extra FORMAT tags, for example
#'   `c("AD", "DP", "GQ")`. Selected fields are typed members of each call's
#'   `format` struct using the header's Type and Number. Missing elements retain
#'   their positions; no allele normalization or depth inference is performed.
#'   NULL or an empty vector keeps the default GT/PS schema. Unknown, empty,
#'   missing and case-insensitively duplicate names error; GT and PS are already
#'   exposed by the typed call fields and cannot be selected again.
#'   Tag lookup uses exact header spelling: declared lowercase `gt` and `ps`
#'   are distinct extra fields. Selected names must not collide under DuckDB's
#'   case-insensitive struct-member lookup.
#' @param include_filter Add the physical record's `FILTER` as a `VARCHAR[]`
#'   column. `PASS` is `c("PASS")`; an unapplied `.` filter is `NULL`; named
#'   failing filters retain header order. The default keeps the existing schema.
#' @param raw_gt Retain exact original VCF genotype text in each call's `raw_gt`
#'   member. The default `FALSE` keeps the typed-call schema unchanged. `TRUE`
#'   preserves leading phase markers, mixed separators and allele spelling;
#'   an absent GT is `NULL`, while a literal missing `.` remains text. BCF input
#'   errors because its encoded genotypes do not retain the original text.
#' @param regions_var Optional name of a session variable on `con` holding typed
#'   intervals, `STRUCT(chrom VARCHAR, start BIGINT, "end" BIGINT)[]`, with
#'   0-based half-open coordinates. See [rduckhts_bcf()] for the contract; the
#'   same rules apply. [rduckhts_geno_sites()] prepares such a plan from
#'   requested sites and joins exact alleles.
#' @return A data frame when `table_name` is `NULL`, otherwise invisible `TRUE`.
#' @seealso [rduckhts_bcf_samples()], [rduckhts_geno_sites()]
#' @examples
#' con <- rduckhts_connect()
#' path <- system.file("extdata", "geno_calls.bcf", package = "Rduckhts")
#' rduckhts_geno(con, "calls", path, non_reference_only = TRUE)
#' DBI::dbGetQuery(con, "SELECT record_index, len(calls) AS n FROM calls ORDER BY record_index")
#' DBI::dbDisconnect(con, shutdown = TRUE)
#' @export
rduckhts_geno <- function(con, table_name = NULL, path, region = NULL,
                          index_path = NULL, samples = NULL,
                          non_reference_only = FALSE, scan_mode = "auto",
                          decompression_threads = 0, decode_error_policy = "null",
                          overwrite = FALSE, format_fields = NULL, raw_gt = FALSE,
                          include_filter = FALSE, regions_var = NULL) {
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE", call. = FALSE)
  }
  .duckhts_check_table_target(con, table_name, overwrite)
  query <- .duckhts_geno_query(
    con, path, region = region, index_path = index_path, samples = samples,
    non_reference_only = non_reference_only, scan_mode = scan_mode,
    decompression_threads = decompression_threads,
    decode_error_policy = decode_error_policy, format_fields = format_fields,
    raw_gt = raw_gt, include_filter = include_filter, regions_var = regions_var)
  if (is.null(table_name)) return(DBI::dbGetQuery(con, query))
  .duckhts_create_table(con, table_name, query, overwrite)
  invisible(TRUE)
}

.duckhts_geno_query <- function(con, path, region = NULL, index_path = NULL,
                                samples = NULL, non_reference_only = FALSE,
                                scan_mode = "auto", decompression_threads = 0,
                                decode_error_policy = "null", format_fields = NULL,
                                raw_gt = FALSE, include_filter = FALSE,
                                regions_var = NULL) {
  params <- list()
  if (!is.null(region)) params$region <- sql_quote_string(con, region)
  if (!is.null(regions_var)) params$regions <- .duckhts_regions_expression(con, regions_var)
  if (!is.null(index_path)) params$index_path <- sql_quote_string(con, index_path)
  if (!is.null(samples)) params$samples <- sql_quote_string(con, samples)
  if (!is.logical(non_reference_only) || length(non_reference_only) != 1L || is.na(non_reference_only)) {
    stop("non_reference_only must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(raw_gt) || length(raw_gt) != 1L || is.na(raw_gt)) {
    stop("raw_gt must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(include_filter) || length(include_filter) != 1L || is.na(include_filter)) {
    stop("include_filter must be TRUE or FALSE", call. = FALSE)
  }
  params$non_reference_only <- if (non_reference_only) "true" else "false"
  params$scan_mode <- sql_quote_string(con, .validate_scan_mode_param(scan_mode))
  params$decompression_threads <- .validate_nonnegative_integer_param(
    decompression_threads, "decompression_threads")
  params$decode_error_policy <- sql_quote_string(con, decode_error_policy)
  params$format_fields <- sql_varchar_list_literal(con, format_fields, "format_fields")
  params$include_filter <- if (include_filter) "true" else "false"
  params$raw_gt <- if (raw_gt) "true" else "false"
  paste0("SELECT * FROM read_geno(", sql_quote_string(con, path), build_param_str(params), ")")
}

#' Read the VCF/BCF Sample Catalog
#'
#' Return one row per selected sample with its original-header zero-based
#' `sample_index` and `sample_name`. Join this relation to `read_geno()` calls
#' from the same unchanged file; names are not duplicated into every call.
#'
#' @inheritParams rduckhts_bcf
#' @return A data frame with `sample_index` and `sample_name` columns.
#' @export
rduckhts_bcf_samples <- function(con, path, samples = NULL) {
  params <- if (is.null(samples)) list() else list(samples = sql_quote_string(con, samples))
  DBI::dbGetQuery(con, paste0("SELECT * FROM read_bcf_samples(", sql_quote_string(con, path),
                             build_param_str(params), ") ORDER BY sample_index"))
}

.geno_sites_columns <- c("request_id", "build", "chrom", "pos", "ref", "alt")
.geno_sites_output_columns <- c("record_index", "full_alt", "alt_index", "calls", "match_status")

# A unique, syntactically plain name for session variables and temp relations.
.geno_sites_scratch_name <- function(kind) {
  gsub("[^A-Za-z0-9_]", "_", basename(tempfile(paste0("duckhts_geno_sites_", kind, "_"))))
}

.geno_sites_source_sql <- function(con, sites) {
  if (inherits(sites, "SQL")) {
    if (length(sites) != 1L || is.na(sites) || !nzchar(sites)) {
      stop("sites must be a relation name, a DBI::Id, or a single DBI::SQL query", call. = FALSE)
    }
    return(paste0("(", as.character(sites), ")"))
  }
  if (inherits(sites, "Id") ||
      (is.character(sites) && length(sites) == 1L && !is.na(sites) && nzchar(sites))) {
    return(as.character(DBI::dbQuoteIdentifier(con, sites)))
  }
  stop("sites must be a relation name, a DBI::Id, or a single DBI::SQL query", call. = FALSE)
}

#' Fetch Exact Alleles at Requested Sites
#'
#' Fetch full-site genotypes for requested (chromosome, position, REF, ALT)
#' alleles from an indexed VCF/BCF. Positions become 0-based half-open
#' intervals for one native indexed scan through `read_geno(regions := ...)`,
#' then exact alleles are restored with an equality join. Everything runs on
#' `con`: the request snapshot and the plan variable are session-local, are
#' removed on exit even after an error, and no private connection is opened.
#'
#' Coordinates in `sites` are 1-based VCF positions. Requested alleles are
#' compared verbatim: nothing is trimmed, left-aligned, strand-flipped, split or
#' lifted over, and `source_build` is only a caller-asserted label that every
#' request must carry (contig names cannot prove assembly identity). The
#' `match_status` labels are position-level diagnostics, not claims about
#' normalization equivalence.
#'
#' The result holds every request column followed by:
#' \describe{
#'   \item{`record_index`}{Scan-local record ordinal (not a persistent locator).}
#'   \item{`full_alt`}{The record's complete `ALT` list in file order.}
#'   \item{`alt_index`}{1-based ordinal of the matching ALT slot.}
#'   \item{`calls`}{The full-site `calls` list of [rduckhts_geno()], including
#'     `format_fields` members such as `GP`, `DS` or `HS`.}
#'   \item{`match_status`}{`"matched"`: chrom, pos, REF and ALT equal a record's
#'     ALT slot (one row per matching record and ALT slot, so physical
#'     duplicate records and duplicated ALT slots each return a row).
#'     `"allele_not_at_site"`: a record at chrom and pos has an equal REF but
#'     none of its ALTs equals `alt`. `"ref_mismatch"`: records exist at chrom
#'     and pos but none has an equal REF. `"absent"`: no record starts at chrom
#'     and pos (an unknown contig is absent).}
#' }
#' Unmatched requests have one row with `NULL` record and call fields; no
#' reference genotype is fabricated. Reference and missing calls at a matched
#' site are returned because `non_reference_only` is fixed to `FALSE`. Rows carry
#' no implicit order: use `ORDER BY request_id, record_index, alt_index`.
#'
#' @inheritParams rduckhts_geno
#' @param sites A relation name (character), a [DBI::Id()], or a [DBI::SQL()]
#'   query evaluated on `con`, with columns `request_id` (unique, non-NULL),
#'   `build`, `chrom`, `pos` (integer type, at least 1), `ref` and `alt` (one
#'   allele, non-empty, no comma). Extra columns are preserved; the names
#'   `record_index`, `full_alt`, `alt_index`, `calls` and `match_status` are
#'   reserved. Duplicate (chrom, pos, ref, alt) rows are separate requests.
#' @param source_build Single string that every `sites$build` must equal.
#' @param table_name Optional table to create with the result (a regular table,
#'   replaced only when `overwrite = TRUE`); its name is returned invisibly.
#'   `NULL` returns a data frame.
#' @param samples,index_path,format_fields,decode_error_policy See
#'   [rduckhts_geno()]. `format_fields` such as `c("GP", "DS", "HS")` are kept
#'   for the whole site without any allele remapping. The decode policy defaults
#'   to `"error"` here.
#' @return A data frame, or invisibly `table_name` when it is given.
#' @seealso [rduckhts_geno()]
#' @examples
#' con <- rduckhts_connect()
#' path <- system.file("extdata", "geno_sites.bcf", package = "Rduckhts")
#' DBI::dbWriteTable(con, "requests", data.frame(
#'   request_id = 1:2, build = "GRCh38", chrom = "chr1", pos = c(100L, 100L),
#'   ref = "A", alt = c("G", "T")))
#' rduckhts_geno_sites(con, path, "requests", "GRCh38", format_fields = c("GP", "DS"))
#' DBI::dbDisconnect(con, shutdown = TRUE)
#' @export
rduckhts_geno_sites <- function(con, path, sites, source_build, table_name = NULL,
                                samples = NULL, format_fields = NULL, index_path = NULL,
                                decode_error_policy = "error", overwrite = FALSE) {
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.character(source_build) || length(source_build) != 1L || is.na(source_build) ||
      !nzchar(source_build)) {
    stop("source_build must be a single non-empty string", call. = FALSE)
  }
  if (!is.null(table_name) &&
      (!is.character(table_name) || length(table_name) != 1L || is.na(table_name) || !nzchar(table_name))) {
    stop("table_name must be NULL or a single non-empty string", call. = FALSE)
  }
  source_sql <- .geno_sites_source_sql(con, sites)
  .duckhts_check_table_target(con, table_name, overwrite)

  snapshot <- .geno_sites_scratch_name("requests")
  variable <- .geno_sites_scratch_name("plan")
  quoted_snapshot <- sql_quote_identifier(con, snapshot)
  on.exit({
    try(DBI::dbExecute(con, paste0("RESET VARIABLE ", sql_quote_identifier(con, variable))), silent = TRUE)
    try(DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS temp.", quoted_snapshot)), silent = TRUE)
  }, add = TRUE)

  # One consistent snapshot: the requests are read here once, on this connection.
  DBI::dbExecute(con, paste0("CREATE TEMP TABLE ", quoted_snapshot,
                             " AS SELECT * FROM ", source_sql))
  described <- DBI::dbGetQuery(con, paste0("SELECT column_name, column_type FROM (DESCRIBE ",
                                           quoted_snapshot, ")"))
  missing <- setdiff(.geno_sites_columns, described$column_name)
  if (length(missing)) {
    stop("sites is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  reserved <- intersect(.geno_sites_output_columns, described$column_name)
  if (length(reserved)) {
    stop("sites uses reserved output column name(s): ", paste(reserved, collapse = ", "), call. = FALSE)
  }
  type_of <- function(name) described$column_type[match(name, described$column_name)]
  integer_types <- c("TINYINT", "SMALLINT", "INTEGER", "BIGINT", "HUGEINT", "UTINYINT",
                     "USMALLINT", "UINTEGER", "UBIGINT", "UHUGEINT")
  if (!type_of("pos") %in% integer_types) {
    stop("sites$pos must have an integer type (1-based VCF position), not ", type_of("pos"),
         call. = FALSE)
  }
  for (name in c("build", "chrom", "ref", "alt")) {
    if (!identical(type_of(name), "VARCHAR")) {
      stop("sites$", name, " must be VARCHAR, not ", type_of(name), call. = FALSE)
    }
  }
  q <- function(name) sql_quote_identifier(con, name)
  checks <- DBI::dbGetQuery(con, paste0(
    "SELECT count(*) AS n, ",
    "count(*) - count(DISTINCT ", q("request_id"), ") AS duplicate_ids, ",
    "count(*) FILTER (WHERE ", q("request_id"), " IS NULL) AS null_ids, ",
    "count(*) FILTER (WHERE ", q("build"), " IS NULL OR ", q("build"), " <> ",
      sql_quote_string(con, source_build), ") AS wrong_build, ",
    "count(*) FILTER (WHERE ", q("chrom"), " IS NULL OR ", q("chrom"), " = '') AS bad_chrom, ",
    "count(*) FILTER (WHERE ", q("pos"), " IS NULL OR ", q("pos"), " < 1) AS bad_pos, ",
    "count(*) FILTER (WHERE ", q("ref"), " IS NULL OR ", q("ref"), " = '') AS bad_ref, ",
    "count(*) FILTER (WHERE ", q("alt"), " IS NULL OR ", q("alt"), " = '' OR position(',' IN ",
      q("alt"), ") > 0) AS bad_alt FROM ", quoted_snapshot))
  problems <- c(
    if (checks$null_ids > 0 || checks$duplicate_ids > 0)
      "request_id must be non-NULL and unique",
    if (checks$wrong_build > 0)
      paste0("build must equal source_build ('", source_build, "') on every request; ",
             checks$wrong_build, " request(s) differ or are NULL"),
    if (checks$bad_chrom > 0) "chrom must be non-NULL and non-empty",
    if (checks$bad_pos > 0) "pos must be non-NULL and at least 1 (1-based VCF position)",
    if (checks$bad_ref > 0) "ref must be non-NULL and non-empty",
    if (checks$bad_alt > 0) "alt must be one non-empty allele (non-NULL, no comma)")
  if (length(problems)) stop("invalid sites: ", paste(problems, collapse = "; "), call. = FALSE)

  # Distinct positions are only I/O planning; the requests themselves are kept.
  # list() over zero rows is NULL, which would mean an unrestricted scan.
  DBI::dbExecute(con, paste0(
    "SET VARIABLE ", sql_quote_identifier(con, variable), " = (SELECT coalesce(",
    "list({'chrom': chrom, 'start': pos - 1, 'end': pos} ORDER BY chrom, pos), ",
    "[]::STRUCT(chrom VARCHAR, start BIGINT, \"end\" BIGINT)[]) FROM (SELECT DISTINCT ",
    q("chrom"), " AS chrom, ", q("pos"), "::BIGINT AS pos FROM ", quoted_snapshot, "))"))
  genotypes <- .duckhts_geno_query(
    con, path, index_path = index_path, samples = samples, non_reference_only = FALSE,
    decode_error_policy = decode_error_policy, format_fields = format_fields,
    regions_var = variable)

  helper <- c("__at_pos", "__ref_ok", "__is_match", "__any_match", "__any_ref", "__rn")
  query <- paste0(
    "WITH alleles AS (SELECT record_index, CHROM, POS, REF, ALT AS full_alt, calls, ",
    "unnest(slots) AS match_alt, generate_subscripts(slots, 1) AS alt_index ",
    "FROM (SELECT *, CASE WHEN len(ALT) = 0 THEN [NULL::VARCHAR] ELSE ALT END AS slots FROM (",
    genotypes, "))), ",
    "joined AS (SELECT s.*, a.record_index AS __rec, a.full_alt AS __full_alt, ",
    "a.alt_index AS __alt_index, a.calls AS __calls, a.CHROM IS NOT NULL AS __at_pos, ",
    "a.REF = s.", q("ref"), " AS __ref_ok, ",
    "a.REF = s.", q("ref"), " AND a.match_alt = s.", q("alt"), " AS __is_match ",
    "FROM ", quoted_snapshot, " s LEFT JOIN alleles a ON s.", q("chrom"), " = a.CHROM AND s.",
    q("pos"), " = a.POS), ",
    "flagged AS (SELECT *, coalesce(bool_or(__is_match) OVER w, false) AS __any_match, ",
    "coalesce(bool_or(__ref_ok) OVER w, false) AS __any_ref, ",
    "row_number() OVER (w ORDER BY __rec, __alt_index) AS __rn ",
    "FROM joined WINDOW w AS (PARTITION BY ", q("request_id"), ")) ",
    "SELECT * EXCLUDE (", paste(c(helper, "__rec", "__full_alt", "__alt_index", "__calls"), collapse = ", "), "), ",
    "CASE WHEN __is_match THEN __rec END AS record_index, ",
    "CASE WHEN __is_match THEN __full_alt END AS full_alt, ",
    "CASE WHEN __is_match THEN __alt_index END AS alt_index, ",
    "CASE WHEN __is_match THEN __calls END AS calls, ",
    "CASE WHEN __is_match THEN 'matched' WHEN __any_ref THEN 'allele_not_at_site' ",
    "WHEN __at_pos THEN 'ref_mismatch' ELSE 'absent' END AS match_status ",
    "FROM flagged WHERE coalesce(__is_match, false) OR (NOT __any_match AND __rn = 1)")
  if (is.null(table_name)) return(DBI::dbGetQuery(con, query))
  .duckhts_create_table(con, table_name, query, overwrite)
  invisible(table_name)
}
