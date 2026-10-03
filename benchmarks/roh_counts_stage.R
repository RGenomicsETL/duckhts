# Stage read counts for duckhts_roh_counts from one indexed BAM or CRAM.
#
# duckhts_somalier_bam_counts counts the bases at the sites of an AF-site
# Parquet (chrom, pos, ref, alt, INFO_AF). The counts are oriented so that
# alt_count counts the allele whose frequency is af, and af is clamped to
# [1e-3, 0.999] as in the ROH ancestry evaluation. When an expected identity
# is given, the counts are compared with it before anything is published.
#
# The identity is the row count and the SHA-256 of a canonical text of the
# counts in (chrom, pos) order; it does not depend on Parquet encoding.

roh_counts_connect <- function(extension) {
  driver <- duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE,
    allow_extensions = TRUE, config = list(threads = "4",
      allow_unsigned_extensions = "true", autoinstall_known_extensions = "false",
      autoload_known_extensions = "false"))
  con <- DBI::dbConnect(driver)
  DBI::dbExecute(con, sprintf("LOAD %s", DBI::dbQuoteString(con, normalizePath(extension))))
  con
}

roh_counts_identity <- function(con, path) {
  identity <- DBI::dbGetQuery(con, sprintf(paste0(
    "SELECT count(*)::BIGINT AS rows, sha256(coalesce(string_agg(concat_ws('|', ",
    "sample_id, chrom, pos, ref, alt, coalesce(ref_count::VARCHAR, 'NA'), ",
    "coalesce(alt_count::VARCHAR, 'NA'), coalesce(other_count::VARCHAR, 'NA'), ",
    "printf('%%.9e', af), status), chr(10) ORDER BY chrom, pos), '')) AS counts_sha256 ",
    "FROM read_parquet(%s)"), DBI::dbQuoteString(con, path)))
  list(rows = as.numeric(identity$rows[[1L]]), counts_sha256 = identity$counts_sha256[[1L]])
}

stage_roh_counts <- function(source, sample_id, af_sites, reference, output, extension,
                             assembly, min_mapq, min_baseq, exclude_flags,
                             overlap_policy, worker_count = 4L, expected = NULL) {
  integers <- c(min_mapq = min_mapq, min_baseq = min_baseq,
                exclude_flags = exclude_flags, worker_count = worker_count)
  if (anyNA(integers) || any(integers < 0) || any(integers != round(integers))) {
    stop("count settings must be non-negative integers", call. = FALSE)
  }
  for (path in c(af_sites, reference)) {
    if (!file.exists(path)) stop("missing staging input: ", path, call. = FALSE)
  }
  con <- roh_counts_connect(extension)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  q <- function(value) as.character(DBI::dbQuoteString(con, value))
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  panel <- paste0(output, ".panel-", Sys.getpid(), ".parquet")
  temporary <- paste0(output, ".partial-", Sys.getpid(), ".parquet")
  on.exit(unlink(c(panel, temporary), force = TRUE), add = TRUE)

  DBI::dbExecute(con, sprintf(paste0(
    "COPY (SELECT %s AS assembly, ",
    "(row_number() OVER (ORDER BY chrom, pos) - 1)::UBIGINT AS site_index, ",
    "chrom AS region, pos::UBIGINT AS position, least(ref, alt) AS allele_a, ",
    "greatest(ref, alt) AS allele_b FROM read_parquet(%s) ORDER BY chrom, pos) ",
    "TO %s (FORMAT PARQUET)"), q(assembly), q(af_sites), q(panel)))
  DBI::dbExecute(con, sprintf(paste0(
    "COPY (WITH counts AS (SELECT * FROM duckhts_somalier_bam_counts(%s, NULL, %s, %s, ",
    "panel_parquet := %s, min_mapq := %d, min_baseq := %d, exclude_flags := %d, ",
    "overlap_policy := %s, worker_count := %d)), ",
    "sites AS (SELECT chrom, pos, ref, alt, ",
    "greatest(1e-3, least(0.999, INFO_AF[1]))::DOUBLE AS af FROM read_parquet(%s)) ",
    "SELECT c.sample_id, s.chrom, s.pos, s.ref, s.alt, ",
    "(CASE WHEN s.alt = c.allele_b THEN c.a ELSE c.b END)::INTEGER AS ref_count, ",
    "(CASE WHEN s.alt = c.allele_b THEN c.b ELSE c.a END)::INTEGER AS alt_count, ",
    "c.other::INTEGER AS other_count, s.af, c.status ",
    "FROM counts AS c JOIN sites AS s ON s.chrom = c.region AND s.pos = c.position ",
    "ORDER BY s.chrom, s.pos) TO %s (FORMAT PARQUET)"),
    q(source), q(sample_id), q(reference), q(panel), as.integer(min_mapq),
    as.integer(min_baseq), as.integer(exclude_flags), q(overlap_policy),
    as.integer(worker_count), q(af_sites), q(temporary)))

  identity <- roh_counts_identity(con, temporary)
  sites <- DBI::dbGetQuery(con, sprintf("SELECT count(*) AS n FROM read_parquet(%s)", q(af_sites)))$n
  if (identity$rows != sites) {
    stop("read counts cover ", identity$rows, " of ", sites, " sites", call. = FALSE)
  }
  if (!is.null(expected) &&
      (!identical(identity$rows, as.numeric(expected$rows)) ||
       !identical(identity$counts_sha256, expected$counts_sha256))) {
    stop("read counts differ from their registered identity: ", identity$rows,
         " rows, counts SHA-256 ", identity$counts_sha256, call. = FALSE)
  }
  if (!file.rename(temporary, output)) stop("could not publish the read counts", call. = FALSE)
  c(identity, list(output = output))
}

# A counts artifact's locator names three artifacts in order: the alignment
# source, the AF-site Parquet and the reference FASTA. Its supplier identity
# holds the sample, the count settings and the expected rows and counts_sha256.
stage_roh_counts_from_registry <- function(id, extension, worker_count = 4L) {
  registry <- duckhtsbench::duckhts_bench_registry()
  row <- registry[registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("unknown or non-unique benchmark artifact: ", id, call. = FALSE)
  inputs <- sub("^artifact:", "", strsplit(row$locator[[1L]], ";", fixed = TRUE)[[1L]])
  if (length(inputs) != 3L || !all(inputs %in% registry$id)) {
    stop("counts artifact must name a source, an AF-site Parquet and a reference", call. = FALSE)
  }
  fields <- duckhtsbench:::duckhts_bench_identity_fields(row$supplier_identity[[1L]])
  required <- c("sample", "assembly", "min_mapq", "min_baseq", "exclude_flags",
                "overlap_policy", "rows", "counts_sha256")
  if (!all(required %in% names(fields))) {
    stop("counts artifact identity lacks: ",
         paste(setdiff(required, names(fields)), collapse = ", "), call. = FALSE)
  }
  source_row <- registry[registry$id == inputs[[1L]], , drop = FALSE]
  # A registered public CRAM is read where it is, by range requests; a local
  # source is its cached file.
  source <- if (grepl("^https?://", source_row$locator[[1L]])) {
    source_row$locator[[1L]]
  } else {
    duckhtsbench::duckhts_bench_artifact_path(inputs[[1L]])
  }
  stage_roh_counts(
    source = source, sample_id = fields[["sample"]],
    af_sites = duckhtsbench::duckhts_bench_artifact_path(inputs[[2L]]),
    reference = duckhtsbench::duckhts_bench_artifact_path(inputs[[3L]]),
    output = duckhtsbench::duckhts_bench_artifact_path(id), extension = extension,
    assembly = fields[["assembly"]], min_mapq = as.integer(fields[["min_mapq"]]),
    min_baseq = as.integer(fields[["min_baseq"]]),
    exclude_flags = as.integer(fields[["exclude_flags"]]),
    overlap_policy = fields[["overlap_policy"]], worker_count = worker_count,
    expected = list(rows = fields[["rows"]], counts_sha256 = fields[["counts_sha256"]]))
}
