#' Stage a synthetic Somalier-shaped Parquet workload without network access.
#'
#' The fixed panel has 17,000 distinct biallelic GRCh38-shaped sites. Its
#' coordinates and counts are arithmetic fixtures, not Somalier upstream sites
#' or biological observations. Evidence has one row per sample and panel site.
#' The default 250-sample paths belong to the registry. Other sample counts
#' and selected-pair topologies are parameterized instances of those artifacts.
#'
#' @param samples Number of synthetic samples, from 2 through 1000.
#' @param selected_samples Number of distinct samples in requested pairs.
#' @param selected_pairs Number of distinct ordered requested pairs.
#' @return Named paths for panel, frequency, evidence, and selected pairs.
#' @export
duckhts_bench_stage_somalier <- function(samples = 250L,
                                        selected_samples = min(samples, 3L),
                                        selected_pairs = selected_samples - 1L) {
  if (length(samples) != 1L || is.na(samples) || !is.numeric(samples) ||
      samples != floor(samples) || samples < 2 || samples > 1000) {
    stop("samples must be one whole number from 2 through 1000", call. = FALSE)
  }
  samples <- as.integer(samples)
  if (length(selected_samples) != 1L || is.na(selected_samples) ||
      !is.numeric(selected_samples) || selected_samples != floor(selected_samples) ||
      selected_samples < 2 || selected_samples > samples) {
    stop("selected_samples must be one whole number from 2 through samples", call. = FALSE)
  }
  selected_samples <- as.integer(selected_samples)
  pair_limit <- min(selected_samples * (selected_samples - 1), 100000L)
  if (length(selected_pairs) != 1L || is.na(selected_pairs) ||
      !is.numeric(selected_pairs) || selected_pairs != floor(selected_pairs) ||
      selected_pairs < selected_samples - 1L || selected_pairs > pair_limit) {
    stop("selected_pairs must cover every selected sample and be at most ", pair_limit,
         call. = FALSE)
  }
  selected_pairs <- as.integer(selected_pairs)
  plan <- duckhts_bench_stage_plan("somalier-synthetic")
  ids <- c(panel = "somalier_synthetic_panel", frequency = "somalier_synthetic_frequency",
           evidence = "somalier_synthetic_evidence",
           pairs = "somalier_synthetic_pairs")
  if (!identical(plan$id, unname(ids)) ||
      !identical(plan$transform,
        c("generate_synthetic_panel_v2", "generate_synthetic_frequency_v2",
          "generate_synthetic_evidence_v2", "generate_synthetic_pairs_v2"))) {
    stop("Somalier registry plan is incomplete or reordered", call. = FALSE)
  }
  relative <- stats::setNames(plan$cache_relpath, names(ids))
  registered <- "samples-250/"
  if (!grepl(registered, relative[["evidence"]], fixed = TRUE) ||
      !grepl("samples-250/selected-3/pairs-2/",
             relative[["pairs"]], fixed = TRUE)) {
    stop("Somalier evidence and pair paths must name the registered default instance",
         call. = FALSE)
  }
  for (name in c("evidence", "pairs")) {
    relative[[name]] <- sub("samples-250", paste0("samples-", samples), relative[[name]],
                            fixed = TRUE)
  }
  relative[["pairs"]] <- sub("selected-3", paste0("selected-", selected_samples),
                              relative[["pairs"]], fixed = TRUE)
  relative[["pairs"]] <- sub("pairs-2", paste0("pairs-", selected_pairs),
                              relative[["pairs"]], fixed = TRUE)
  paths <- stats::setNames(vapply(relative, duckhts_bench_cache_path, character(1)), names(ids))
  if (!requireNamespace("DBI", quietly = TRUE) ||
      !requireNamespace("duckdb", quietly = TRUE)) {
    stop("DBI and duckdb are required for Somalier Parquet staging", call. = FALSE)
  }
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, "SET threads=1")
  DBI::dbExecute(con, "CREATE TEMP TABLE panel AS SELECT 'GRCh38'::VARCHAR AS assembly,
    site_index::UBIGINT AS site_index,
    ('chr' || ((site_index % 22) + 1)::VARCHAR)::VARCHAR AS region,
    (1000000 + (site_index // 22) * 101)::UBIGINT AS position,
    CASE site_index % 4 WHEN 0 THEN 'A' WHEN 1 THEN 'A'
      WHEN 2 THEN 'C' ELSE 'G' END::VARCHAR AS allele_a,
    CASE site_index % 4 WHEN 0 THEN 'C' WHEN 1 THEN 'G'
      WHEN 2 THEN 'T' ELSE 'T' END::VARCHAR AS allele_b
    FROM range(17000) sites(site_index)")
  DBI::dbExecute(con, "CREATE TEMP TABLE frequency AS SELECT p.*,
    ((10 + (site_index * 17) % 39)::DOUBLE / 100)::DOUBLE AS population_b_af
    FROM panel p")
  DBI::dbExecute(con, sprintf("CREATE TEMP TABLE evidence AS
    WITH sample_numbers AS (SELECT sample_number FROM range(1, %d) s(sample_number)),
    states AS (SELECT p.*, sample_number,
      (p.site_index * 7 + sample_number * 13) %% 10 AS genotype_code,
      (p.site_index + sample_number * 17) %% 97 = 0 AS unavailable,
      (p.site_index + sample_number * 11) %% 23 = 0 AS minor_count
      FROM panel p CROSS JOIN sample_numbers)
    SELECT printf('sample%%03d', sample_number)::VARCHAR AS sample_id,
      assembly, site_index, region, position, allele_a, allele_b,
      CASE WHEN unavailable THEN NULL::UBIGINT
        WHEN genotype_code <= 3 THEN (30 - minor_count::INTEGER)::UBIGINT
        WHEN genotype_code <= 6 THEN 15::UBIGINT
        ELSE minor_count::UBIGINT END AS a,
      CASE WHEN unavailable THEN NULL::UBIGINT
        WHEN genotype_code <= 3 THEN minor_count::UBIGINT
        WHEN genotype_code <= 6 THEN 15::UBIGINT
        ELSE (30 - minor_count::INTEGER)::UBIGINT END AS b,
      CASE WHEN unavailable THEN NULL::UBIGINT ELSE 0::UBIGINT END AS other
    FROM states", samples + 1L))
  DBI::dbExecute(con, sprintf("CREATE TEMP TABLE pairs AS
    WITH star AS (SELECT receiver_number, 1 AS anchor_number
      FROM range(2, %d) s(receiver_number)),
    extras AS (SELECT receiver_number, anchor_number
      FROM range(1, %d) r(receiver_number)
      CROSS JOIN range(1, %d) a(anchor_number)
      WHERE receiver_number != anchor_number AND anchor_number != 1
      ORDER BY receiver_number, anchor_number LIMIT %d),
    selected AS (SELECT * FROM star UNION ALL SELECT * FROM extras)
    SELECT printf('sample%%03d', receiver_number)::VARCHAR AS sample_a,
      printf('sample%%03d', anchor_number)::VARCHAR AS sample_b,
      printf('sample%%03d', receiver_number)::VARCHAR AS receiver_id,
      printf('sample%%03d', anchor_number)::VARCHAR AS anchor_id
    FROM selected ORDER BY receiver_number, anchor_number",
    selected_samples + 1L, selected_samples + 1L, selected_samples + 1L,
    selected_pairs - selected_samples + 1L))
  tables <- c(panel = "panel", frequency = "frequency", evidence = "evidence", pairs = "pairs")
  expected_rows <- c(panel = 17000, frequency = 17000,
                     evidence = 17000 * samples, pairs = selected_pairs)
  generator_sha256 <- digest::digest(paste(deparse(body(duckhts_bench_stage_somalier)),
                                         collapse = "\n"), algo = "sha256")
  duckdb_version <- as.character(utils::packageVersion("duckdb"))
  writer_options <- "FORMAT=PARQUET;COMPRESSION=ZSTD;ROW_GROUP_SIZE=100000;THREADS=1"
  r_version <- as.character(getRversion())
  instance_keys <- c(panel = "fixed", frequency = "fixed",
                     evidence = paste0("samples=", samples),
                     pairs = paste0("samples=", samples, ";selected_samples=",
                                    selected_samples, ";selected_pairs=", selected_pairs))
  for (name in names(tables)) {
    output <- paths[[name]]
    receipt <- paste0(output, ".provenance.tsv")
    if (file.exists(output) || file.exists(receipt)) {
      if (!file.exists(output) || !file.exists(receipt)) {
        stop("incomplete existing Somalier artifact: ", output, call. = FALSE)
      }
      retained <- utils::read.delim(receipt, colClasses = "character", check.names = FALSE)
      expected <- duckhts_bench_provenance_fields(ids[[name]], output)
      if (!identical(names(retained), c("field", "value")) || anyDuplicated(retained$field) ||
          !all(expected$field %in% retained$field) ||
          !identical(retained$value[match(expected$field, retained$field)], expected$value)) {
        stop("existing Somalier artifact provenance does not match the registry: ",
             output, call. = FALSE)
      }
      fields <- stats::setNames(retained$value, retained$field)
      if (name %in% c("evidence", "pairs") &&
          !identical(fields[["samples"]], as.character(samples))) {
        stop("existing Somalier artifact has another sample count: ", output,
             call. = FALSE)
      }
      if (!identical(fields[["sites"]], "17000") ||
          !identical(fields[["rows"]], as.character(expected_rows[[name]])) ||
          !identical(fields[["bytes"]], as.character(file.info(output)$size)) ||
          !identical(fields[["sha256"]], digest::digest(file = output, algo = "sha256")) ||
          !identical(fields[["generator"]], "somalier-synthetic-v2") ||
          !identical(fields[["instance_key"]], instance_keys[[name]]) ||
          !identical(fields[["generator_sha256"]], generator_sha256) ||
          !identical(fields[["duckdb_version"]], duckdb_version) ||
          !identical(fields[["r_version"]], r_version) ||
          !identical(fields[["writer_options"]], writer_options)) {
        stop("existing Somalier artifact receipt does not match bytes: ", output,
             call. = FALSE)
      }
    } else {
      dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
      temporary <- paste0(output, ".partial-", Sys.getpid())
      if (file.exists(temporary)) stop("Somalier staging temporary exists: ", temporary,
                                       call. = FALSE)
      on.exit(unlink(temporary), add = TRUE)
      DBI::dbExecute(con, sprintf("COPY %s TO %s
        (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)",
        tables[[name]], as.character(DBI::dbQuoteString(con, temporary))))
      rows <- DBI::dbGetQuery(con, sprintf("SELECT count(*) AS rows FROM read_parquet(%s)",
        as.character(DBI::dbQuoteString(con, temporary))))$rows[[1L]]
      if (rows != expected_rows[[name]]) stop("Somalier Parquet row count changed", call. = FALSE)
      if (!file.rename(temporary, output)) stop("could not publish Somalier artifact", call. = FALSE)
      provenance <- duckhts_bench_provenance_fields(ids[[name]], output)
      provenance <- rbind(provenance, data.frame(
        field = c("samples", "sites", "rows", "bytes", "sha256", "generator",
                  "instance_key", "generator_sha256", "duckdb_version", "r_version",
                  "writer_options"),
        value = c(if (name %in% c("evidence", "pairs")) samples else "not_applicable",
                  17000, rows, file.info(output)$size,
                  digest::digest(file = output, algo = "sha256"), "somalier-synthetic-v2",
                  instance_keys[[name]], generator_sha256, duckdb_version, r_version,
                  writer_options),
        stringsAsFactors = FALSE))
      utils::write.table(provenance, receipt, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
  observed <- DBI::dbGetQuery(con, "SELECT count(*) AS rows,
    count(DISTINCT sample_id) AS samples,
    count(*) FILTER (WHERE a IS NULL AND b IS NULL AND other IS NULL) AS unavailable,
    count(*) FILTER (WHERE a IS NOT NULL AND a + b = 30 AND other = 0) AS measured
    FROM evidence")
  if (observed$rows[[1L]] != expected_rows[["evidence"]] ||
      observed$samples[[1L]] != samples ||
      observed$unavailable[[1L]] + observed$measured[[1L]] != expected_rows[["evidence"]]) {
    stop("Somalier source count evidence violates its fixture contract", call. = FALSE)
  }
  paths
}
