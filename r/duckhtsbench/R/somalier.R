duckhts_bench_stage_somalier_v034 <- function() {
  system <- Sys.info()
  if (!identical(unname(system[["sysname"]]), "Linux") ||
      !unname(system[["machine"]]) %in% c("x86_64", "amd64")) {
    stop("the pinned Somalier v0.3.4 executable requires Linux x86-64",
         call. = FALSE)
  }
  executable <- duckhts_bench_fetch("somalier_v034_linux_x86_64")
  Sys.chmod(executable, mode = "0755")
  help <- suppressWarnings(system2(executable, "--help", stdout = TRUE,
                                   stderr = TRUE))
  status <- attr(help, "status")
  if ((!is.null(status) && status != 0L) ||
      !any(grepl("somalier version: 0.3.4", help, fixed = TRUE))) {
    stop("the staged executable did not identify itself as Somalier v0.3.4",
         call. = FALSE)
  }
  executable
}

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

duckhts_bench_somalier_site_panel <- function() {
  site_index <- 0:16999
  allele_index <- site_index %% 4L + 1L
  data.frame(
    assembly = rep("Synthetic_GRCh38_shaped_v1", length(site_index)),
    site_index = site_index,
    region = paste0("chr", site_index %% 22L + 1L),
    position = 1001L + (site_index %/% 22L) * 150L,
    allele_a = c("A", "A", "C", "G")[allele_index],
    allele_b = c("C", "G", "T", "T")[allele_index],
    stringsAsFactors = FALSE
  )
}

duckhts_bench_somalier_samtools <- function(samtools, arguments, error,
                                            capture = FALSE) {
  if (capture) {
    # Parse stdout only: samtools warnings (for example an index older than
    # its data file) go to stderr and must not become parsed rows.
    stderr_file <- tempfile("duckhtsbench-samtools-")
    on.exit(unlink(stderr_file), add = TRUE)
    output <- suppressWarnings(system2(
      samtools, arguments, stdout = TRUE, stderr = stderr_file
    ))
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    if (status != 0L) {
      messages <- c(output, readLines(stderr_file, warn = FALSE))
      stop(error, if (length(messages)) paste0(": ", utils::tail(messages, 1L)) else "",
           call. = FALSE)
    }
    return(output)
  }
  status <- system2(samtools, arguments)
  if (status != 0L) stop(error, call. = FALSE)
  invisible(character())
}

duckhts_bench_write_somalier_reference <- function(path, panel) {
  regions <- paste0("chr", seq_len(22L))
  sequences <- stats::setNames(vector("list", length(regions)), regions)
  connection <- file(path, open = "wb")
  on.exit(close(connection))
  for (region in regions) {
    sites <- panel[panel$region == region, , drop = FALSE]
    sequence_length <- max(sites$position) + 50L
    sequence <- charToRaw(strrep("A", sequence_length))
    sequence[sites$position] <- charToRaw(paste0(sites$allele_a, collapse = ""))
    sequences[[region]] <- rawToChar(sequence)
    starts <- seq.int(1L, sequence_length, by = 80L)
    lines <- substring(sequences[[region]], starts, pmin(starts + 79L, sequence_length))
    writeLines(c(paste0(">", region), lines), connection, sep = "\n", useBytes = TRUE)
  }
  sequences
}

duckhts_bench_write_somalier_sam <- function(path, panel, sequences) {
  regions <- names(sequences)
  lengths <- nchar(unname(sequences), type = "bytes")
  header <- c(
    "@HD\tVN:1.6\tSO:coordinate",
    paste0("@SQ\tSN:", regions, "\tLN:", lengths),
    "@RG\tID:synthetic\tSM:synthetic"
  )
  connection <- file(path, open = "wb")
  on.exit(close(connection))
  writeLines(header, connection, sep = "\n", useBytes = TRUE)

  ordered <- panel[order(match(panel$region, regions), panel$position), , drop = FALSE]
  quality <- strrep("I", 101L)
  for (offset in seq.int(1L, nrow(ordered), by = 1000L)) {
    rows <- ordered[offset:min(offset + 999L, nrow(ordered)), , drop = FALSE]
    lines <- vapply(seq_len(nrow(rows)), function(i) {
      site <- rows[i, , drop = FALSE]
      start <- site$position - 50L
      sequence <- charToRaw(substring(
        sequences[[site$region]], start, start + 100L
      ))
      observed <- if (site$site_index %% 2L == 0L) site$allele_a else site$allele_b
      sequence[[51L]] <- charToRaw(observed)
      paste(
        sprintf("site%05d", site$site_index), "0", site$region, start,
        "60", "101M", "*", "0", "0", rawToChar(sequence), quality,
        "RG:Z:synthetic", sep = "\t"
      )
    }, character(1L))
    writeLines(lines, connection, sep = "\n", useBytes = TRUE)
  }
}

duckhts_bench_validate_somalier_site_extraction <- function(paths, samtools) {
  if (any(!file.exists(paths)) || any(file.info(paths)$size <= 0)) {
    stop("Somalier site-extraction staging is incomplete", call. = FALSE)
  }
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  panel_path <- as.character(DBI::dbQuoteString(con, paths[["panel"]]))
  observed <- DBI::dbGetQuery(con, paste(
    "SELECT count(*) AS sites, count(DISTINCT (region, position)) AS positions,",
    "min(site_index) AS first_site, max(site_index) AS last_site",
    "FROM read_parquet(", panel_path, ")"
  ))
  if (observed$sites[[1L]] != 17000 || observed$positions[[1L]] != 17000 ||
      observed$first_site[[1L]] != 0 || observed$last_site[[1L]] != 16999) {
    stop("Somalier site-extraction panel contract changed", call. = FALSE)
  }

  fai <- utils::read.delim(paths[["reference_fai"]], header = FALSE,
                           colClasses = "character", check.names = FALSE)
  if (nrow(fai) != 22L || !identical(fai[[1L]], paste0("chr", seq_len(22L)))) {
    stop("Somalier site-extraction reference index contract changed", call. = FALSE)
  }
  duckhts_bench_somalier_samtools(
    samtools,
    c("quickcheck", "-v", shQuote(paths[["bam"]]), shQuote(paths[["cram"]])),
    "Somalier BAM/CRAM failed samtools quickcheck"
  )
  for (format in c("bam", "cram")) {
    arguments <- c("view", "-c")
    if (format == "cram") {
      arguments <- c(arguments, "-T", shQuote(paths[["reference"]]))
    }
    count <- duckhts_bench_somalier_samtools(
      samtools, c(arguments, shQuote(paths[[format]])),
      paste("could not count staged Somalier", toupper(format)), capture = TRUE
    )
    if (length(count) != 1L || suppressWarnings(as.numeric(count)) != 17000) {
      stop("Somalier site-extraction source-read denominator changed", call. = FALSE)
    }
    index_rows <- duckhts_bench_somalier_samtools(
      samtools, c("idxstats", shQuote(paths[[format]])),
      paste("could not inspect staged Somalier", toupper(format), "index"),
      capture = TRUE
    )
    fields <- strsplit(index_rows, "\t", fixed = TRUE)
    mapped <- sum(vapply(fields, function(row) as.numeric(row[[3L]]), numeric(1L)))
    if (mapped != 17000) {
      stop("Somalier site-extraction alignment index denominator changed", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Internal staging entry point for the reproducible alignment-extraction report.
duckhts_bench_stage_somalier_site_extraction <- function(
    samtools = Sys.which("samtools")) {
  if (length(samtools) != 1L || is.na(samtools) || !nzchar(samtools)) {
    stop("samtools is required for Somalier site-extraction staging", call. = FALSE)
  }
  if (!requireNamespace("DBI", quietly = TRUE) ||
      !requireNamespace("duckdb", quietly = TRUE)) {
    stop("DBI and duckdb are required for Somalier site-extraction staging",
         call. = FALSE)
  }
  plan <- duckhts_bench_stage_plan("somalier-site-extraction")
  ids <- c(
    panel = "somalier_site_panel",
    reference = "somalier_site_reference",
    reference_fai = "somalier_site_reference_fai",
    bam = "somalier_site_bam",
    bam_bai = "somalier_site_bam_bai",
    cram = "somalier_site_cram",
    cram_crai = "somalier_site_cram_crai"
  )
  transforms <- c(
    "generate_somalier_site_panel_v1", "generate_somalier_site_reference_v1",
    "samtools_faidx", "generate_somalier_site_bam_v1", "samtools_bam_index",
    "samtools_bam_to_cram", "samtools_cram_index"
  )
  if (!identical(plan$id, unname(ids)) || !identical(plan$transform, transforms) ||
      !all(plan$access %in% c("local_generated", "local_derived"))) {
    stop("Somalier site-extraction registry plan is incomplete or reordered",
         call. = FALSE)
  }
  paths <- stats::setNames(
    vapply(plan$id, duckhts_bench_artifact_path, character(1L)), names(ids)
  )
  receipts <- stats::setNames(paste0(paths, ".provenance.tsv"), names(paths))
  present <- file.exists(c(paths, receipts))
  if (any(present) && !all(present)) {
    stop("incomplete existing Somalier site-extraction workload", call. = FALSE)
  }
  if (all(present)) {
    for (name in names(paths)) {
      retained <- utils::read.delim(
        receipts[[name]], colClasses = "character", check.names = FALSE,
        stringsAsFactors = FALSE
      )
      expected <- duckhts_bench_provenance_fields(ids[[name]], paths[[name]])
      fields <- stats::setNames(retained$value, retained$field)
      if (!identical(names(retained), c("field", "value")) ||
          anyDuplicated(retained$field) || !all(expected$field %in% retained$field) ||
          !identical(retained$value[match(expected$field, retained$field)], expected$value) ||
          !identical(fields[["generator"]], "somalier-site-extraction-v1") ||
          !identical(fields[["panel_sites"]], "17000") ||
          !identical(fields[["source_reads"]], "17000") ||
          !identical(fields[["artifact_bytes"]],
                     as.character(file.info(paths[[name]])$size)) ||
          !identical(fields[["artifact_sha256"]],
                     digest::digest(file = paths[[name]], algo = "sha256"))) {
        stop("existing Somalier site-extraction provenance is invalid: ",
             paths[[name]], call. = FALSE)
      }
    }
    duckhts_bench_validate_somalier_site_extraction(paths, samtools)
    return(paths)
  }

  destination <- dirname(paths[["panel"]])
  if (!all(dirname(paths) == destination)) {
    stop("Somalier site-extraction artifacts must share one cache directory",
         call. = FALSE)
  }
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  work <- tempfile("somalier-site-extraction-", tmpdir = destination)
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
  temporary <- stats::setNames(file.path(work, basename(paths)), names(paths))
  panel <- duckhts_bench_somalier_site_panel()
  sequences <- duckhts_bench_write_somalier_reference(temporary[["reference"]], panel)
  sam <- file.path(work, "source.sam")
  duckhts_bench_write_somalier_sam(sam, panel, sequences)

  con <- DBI::dbConnect(duckdb::duckdb())
  DBI::dbWriteTable(con, "somalier_site_panel_source", panel, temporary = TRUE)
  panel_output <- as.character(DBI::dbQuoteString(con, temporary[["panel"]]))
  DBI::dbExecute(con, paste0(
    "COPY (SELECT CAST(assembly AS VARCHAR) AS assembly, ",
    "CAST(site_index AS UBIGINT) AS site_index, CAST(region AS VARCHAR) AS region, ",
    "CAST(position AS UBIGINT) AS position, CAST(allele_a AS VARCHAR) AS allele_a, ",
    "CAST(allele_b AS VARCHAR) AS allele_b FROM somalier_site_panel_source ",
    "ORDER BY site_index) TO ", panel_output,
    " (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 20000)"
  ))
  DBI::dbDisconnect(con, shutdown = TRUE)

  duckhts_bench_somalier_samtools(
    samtools, c("faidx", shQuote(temporary[["reference"]])),
    "could not index synthetic Somalier reference"
  )
  duckhts_bench_somalier_samtools(
    samtools, c("view", "--no-PG", "-@", "1", "-b", "-o",
                shQuote(temporary[["bam"]]), shQuote(sam)),
    "could not create synthetic Somalier BAM"
  )
  duckhts_bench_somalier_samtools(
    samtools, c("index", "-@", "1", "-b", shQuote(temporary[["bam"]]),
                shQuote(temporary[["bam_bai"]])),
    "could not index synthetic Somalier BAM"
  )
  duckhts_bench_somalier_samtools(
    samtools, c("view", "--no-PG", "-@", "1", "-C", "-T",
                shQuote(temporary[["reference"]]), "-o", shQuote(temporary[["cram"]]),
                shQuote(temporary[["bam"]])),
    "could not create synthetic Somalier CRAM"
  )
  duckhts_bench_somalier_samtools(
    samtools, c("index", "-@", "1", shQuote(temporary[["cram"]]),
                shQuote(temporary[["cram_crai"]])),
    "could not index synthetic Somalier CRAM"
  )
  duckhts_bench_validate_somalier_site_extraction(temporary, samtools)

  published <- character()
  complete <- FALSE
  on.exit(if (!complete) {
    unlink(c(published, paste0(published, ".provenance.tsv")), force = TRUE)
  }, add = TRUE)
  for (name in names(paths)) {
    if (!file.rename(temporary[[name]], paths[[name]])) {
      stop("could not publish Somalier site-extraction artifact: ", name,
           call. = FALSE)
    }
    published <- c(published, paths[[name]])
  }
  samtools_version <- duckhts_bench_somalier_samtools(
    samtools, "--version", "could not identify samtools", capture = TRUE
  )[[1L]]
  for (name in names(paths)) {
    provenance <- duckhts_bench_provenance_fields(ids[[name]], paths[[name]])
    provenance <- rbind(provenance, data.frame(
      field = c("generator", "panel_sites", "source_reads", "contigs",
                "read_length", "artifact_bytes", "artifact_sha256",
                "samtools_version", "r_version"),
      value = c("somalier-site-extraction-v1", "17000", "17000", "22", "101",
                as.character(file.info(paths[[name]])$size),
                digest::digest(file = paths[[name]], algo = "sha256"), samtools_version,
                as.character(getRversion())),
      stringsAsFactors = FALSE
    ))
    utils::write.table(provenance, receipts[[name]], sep = "\t", row.names = FALSE,
                       quote = FALSE)
  }
  duckhts_bench_validate_somalier_site_extraction(paths, samtools)
  complete <- TRUE
  paths
}
