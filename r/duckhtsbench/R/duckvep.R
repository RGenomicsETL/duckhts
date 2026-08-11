duckhts_bench_duckvep_sql_literal <- function(value) {
  paste0("'", gsub("'", "''", value, fixed = TRUE), "'")
}

duckhts_bench_duckvep_sha256_file <- function(path) {
  command <- Sys.which("sha256sum")
  arguments <- shQuote(path)
  if (!nzchar(command)) {
    command <- Sys.which("shasum")
    arguments <- c("-a", "256", shQuote(path))
  }
  if (!nzchar(command)) {
    stop("sha256sum or shasum is required to stage the DuckVEP model", call. = FALSE)
  }
  output <- trimws(system2(command, arguments, stdout = TRUE))
  if (!length(output)) stop("could not hash DuckVEP model input: ", path, call. = FALSE)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

duckhts_bench_duckvep_sha256_text <- function(value) {
  path <- tempfile("duckhtsbench-sha256-")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  connection <- file(path, open = "wb")
  writeBin(charToRaw(enc2utf8(value)), connection)
  close(connection)
  duckhts_bench_duckvep_sha256_file(path)
}

duckhts_bench_duckvep_identity <- function(value) {
  if (!length(value) || is.na(value) || !nzchar(value)) return(character())
  fields <- strsplit(value, ";", fixed = TRUE)[[1L]]
  stats::setNames(sub("^[^=]+=", "", fields), sub("=.*$", "", fields))
}

duckhts_bench_duckvep_sources <- function() {
  registry <- duckhts_bench_registry()
  rows <- registry[
    grepl("duckvep_ensembl116_model", registry$consumer, fixed = TRUE) &
      grepl("^source_(manifest|schema|table):(core|funcgen)(:|$)", registry$role),
    , drop = FALSE
  ]
  if (!nrow(rows)) stop("registry has no Ensembl 116 DuckVEP model sources", call. = FALSE)

  parts <- strsplit(rows$role, ":", fixed = TRUE)
  valid <- lengths(parts) %in% c(2L, 3L)
  if (!all(valid)) stop("invalid Ensembl 116 model source role in registry", call. = FALSE)
  rows$source_kind <- vapply(parts, `[[`, character(1L), 1L)
  rows$source_group <- vapply(parts, `[[`, character(1L), 2L)
  rows$source_table <- vapply(parts, function(part) if (length(part) == 3L) part[[3L]] else NA_character_, character(1L))

  required <- list(
    core = c(
      "meta", "coord_system", "seq_region", "gene", "transcript", "exon",
      "exon_transcript", "translation", "attrib_type", "seq_region_attrib",
      "transcript_attrib", "translation_attrib"
    ),
    funcgen = c("meta", "regulatory_feature", "feature_type", "motif_feature")
  )
  releases <- c(
    core = "Ensembl_116_homo_sapiens_core_116_38",
    funcgen = "Ensembl_116_homo_sapiens_funcgen_116_38"
  )
  databases <- c(
    core = "homo_sapiens_core_116_38",
    funcgen = "homo_sapiens_funcgen_116_38"
  )

  for (group in names(required)) {
    selected <- rows[rows$source_group == group, , drop = FALSE]
    if (nrow(selected[selected$source_kind == "source_manifest", , drop = FALSE]) != 1L ||
        nrow(selected[selected$source_kind == "source_schema", , drop = FALSE]) != 1L ||
        !identical(sort(selected$source_table[selected$source_kind == "source_table"]), sort(required[[group]])) ||
        any(selected$release != releases[[group]])) {
      stop("registry does not declare the exact Ensembl 116 ", group, " model inputs", call. = FALSE)
    }
  }
  if (anyDuplicated(rows$id) || anyDuplicated(rows$role)) {
    stop("registry Ensembl 116 model source IDs and roles must be unique", call. = FALSE)
  }
  rows$source_database <- unname(databases[rows$source_group])
  rows
}

duckhts_bench_duckvep_checksums <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fields <- strsplit(trimws(lines), "[[:space:]]+")
  valid <- lengths(fields) == 3L & vapply(fields, function(x) all(grepl("^[0-9]+$", x[1:2])), logical(1L))
  if (!length(lines) || !all(valid)) stop("invalid Ensembl CHECKSUMS manifest: ", path, call. = FALSE)
  data.frame(
    sum = vapply(fields, `[[`, character(1L), 1L),
    blocks = vapply(fields, `[[`, character(1L), 2L),
    filename = vapply(fields, `[[`, character(1L), 3L),
    stringsAsFactors = FALSE
  )
}

duckhts_bench_duckvep_validate_manifests <- function(rows) {
  manifest_hashes <- character()
  for (group in c("core", "funcgen")) {
    manifest_row <- rows[rows$source_group == group & rows$source_kind == "source_manifest", , drop = FALSE]
    manifest_path <- duckhts_bench_artifact_path(manifest_row$id[[1L]])
    duckhts_bench_validate_identity(manifest_row$id[[1L]], manifest_path)
    manifest <- duckhts_bench_duckvep_checksums(manifest_path)
    inputs <- rows[rows$source_group == group & rows$source_kind != "source_manifest", , drop = FALSE]
    for (index in seq_len(nrow(inputs))) {
      input <- inputs[index, , drop = FALSE]
      path <- duckhts_bench_artifact_path(input$id[[1L]])
      duckhts_bench_validate_identity(input$id[[1L]], path)
      filename <- basename(path)
      expected <- manifest[manifest$filename == filename, , drop = FALSE]
      identity <- duckhts_bench_duckvep_identity(input$supplier_identity[[1L]])
      if (nrow(expected) != 1L || !all(c("sum", "blocks") %in% names(identity)) ||
          identity[["sum"]] != expected$sum[[1L]] || identity[["blocks"]] != expected$blocks[[1L]]) {
        stop("registry identity does not match Ensembl CHECKSUMS for ", input$id[[1L]], call. = FALSE)
      }
    }
    manifest_hashes[[manifest_row$source_database[[1L]]]] <- duckhts_bench_duckvep_sha256_file(manifest_path)
  }
  manifest_order <- order(names(manifest_hashes))
  content <- paste0(
    names(manifest_hashes)[manifest_order], "\t", unname(manifest_hashes[manifest_order]), "\n",
    collapse = ""
  )
  list(
    source_manifest_sha256 = duckhts_bench_duckvep_sha256_text(content),
    database_manifest_sha256 = manifest_hashes
  )
}

duckhts_bench_duckvep_schema <- function(path, database, tables) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  lines <- readLines(connection, warn = FALSE)
  if (!any(grepl(paste0("^-- Host: .* Database: ", database, "$"), lines))) {
    stop("Ensembl schema does not declare database ", database, call. = FALSE)
  }

  result <- vector("list", length(tables))
  names(result) <- tables
  for (table in tables) {
    start <- which(lines == paste0("CREATE TABLE `", table, "` ("))
    if (length(start) != 1L) stop("Ensembl schema lacks one table definition for ", table, call. = FALSE)
    finish <- which(seq_along(lines) > start & grepl("^\\) ENGINE=", lines))
    if (!length(finish)) stop("unterminated Ensembl table definition for ", table, call. = FALSE)
    definition <- lines[(start + 1L):(finish[[1L]] - 1L)]
    definition <- definition[grepl("^  `", definition)]
    column <- sub("^  `([^`]+)`.*$", "\\1", definition)
    mysql_type <- tolower(sub("^  `[^`]+`[[:space:]]+([^ (,]+).*$", "\\1", definition))
    duckdb_type <- ifelse(
      mysql_type %in% c("tinyint", "smallint", "mediumint", "int", "integer", "bigint", "year"),
      "BIGINT",
      ifelse(mysql_type %in% c("float", "double", "real", "decimal", "numeric"), "DOUBLE", "VARCHAR")
    )
    if (!length(column) || anyDuplicated(column)) stop("invalid Ensembl columns for ", table, call. = FALSE)
    result[[table]] <- stats::setNames(duckdb_type, column)
  }
  result
}

duckhts_bench_duckvep_import_table <- function(connection, schema, table, columns, path) {
  input <- gzfile(path, open = "rb")
  first_byte <- readBin(input, what = "raw", n = 1L)
  close(input)
  if (!length(first_byte)) {
    definitions <- paste(sprintf("\"%s\" %s", names(columns), unname(columns)), collapse = ", ")
    DBI::dbExecute(connection, sprintf(
      "CREATE TABLE \"%s\".\"%s\" (%s)", schema, table, definitions
    ))
    return(invisible(NULL))
  }

  column_spec <- paste(
    sprintf("%s: %s", duckhts_bench_duckvep_sql_literal(names(columns)), duckhts_bench_duckvep_sql_literal(unname(columns))),
    collapse = ", "
  )
  sql <- sprintf(
    paste0(
      "CREATE TABLE \"%s\".\"%s\" AS SELECT * FROM read_csv(%s, ",
      "delim = '\\t', header = false, columns = {%s}, nullstr = '\\N', ",
      "quote = '', compression = 'gzip', strict_mode = true, parallel = true)"
    ),
    schema, table, duckhts_bench_duckvep_sql_literal(path), column_spec
  )
  DBI::dbExecute(connection, sql)
}

duckhts_bench_duckvep_require_metadata <- function(connection, schema, key, value, species_id = NULL) {
  query <- sprintf(
    "SELECT DISTINCT CAST(meta_value AS VARCHAR) AS value FROM \"%s\".meta WHERE meta_key = %s%s",
    schema,
    duckhts_bench_duckvep_sql_literal(key),
    if (is.null(species_id)) "" else paste0(" AND species_id = ", as.integer(species_id))
  )
  observed <- DBI::dbGetQuery(connection, query)$value
  if (!identical(observed, value)) {
    stop("Ensembl ", schema, " metadata does not declare ", key, "=", value, call. = FALSE)
  }
}

duckhts_bench_duckvep_compare_receipts <- function(stored, rebuilt) {
  if (nrow(stored) != 1L || nrow(rebuilt) != 1L || !identical(names(stored), names(rebuilt))) return(FALSE)
  all(vapply(names(stored), function(name) {
    left <- ifelse(is.na(stored[[name]]), "<NA>", as.character(stored[[name]]))
    right <- ifelse(is.na(rebuilt[[name]]), "<NA>", as.character(rebuilt[[name]]))
    identical(left, right)
  }, logical(1L)))
}

duckhts_bench_validate_duckvep_ensembl116_model <- function(
  model,
  extension,
  source_manifest_sha256
) {
  driver <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  connection <- DBI::dbConnect(driver, dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(connection, paste("LOAD", duckhts_bench_duckvep_sql_literal(extension)))
  DBI::dbExecute(connection, sprintf(
    "ATTACH %s AS staged (READ_ONLY)", duckhts_bench_duckvep_sql_literal(model)
  ))

  metadata <- DBI::dbGetQuery(connection, "SELECT * FROM staged.main.model_metadata")
  expected_metadata <- c(
    source_name = "Ensembl", source_release = "116", vep_release = "116",
    species = "homo_sapiens", species_id = "1", assembly = "GRCh38"
  )
  if (nrow(metadata) != 1L || any(vapply(names(expected_metadata), function(name) {
    !name %in% names(metadata) || as.character(metadata[[name]][[1L]]) != expected_metadata[[name]]
  }, logical(1L))) || metadata$source_manifest_sha256[[1L]] != source_manifest_sha256) {
    stop("DuckVEP model metadata does not match Ensembl 116 homo_sapiens GRCh38", call. = FALSE)
  }

  required <- data.frame(
    table_schema = c(rep("ensembl_core", 12L), rep("ensembl_funcgen", 4L)),
    table_name = c(
      "meta", "coord_system", "seq_region", "gene", "transcript", "exon",
      "exon_transcript", "translation", "attrib_type", "seq_region_attrib",
      "transcript_attrib", "translation_attrib", "meta", "regulatory_feature",
      "feature_type", "motif_feature"
    ),
    stringsAsFactors = FALSE
  )
  present <- DBI::dbGetQuery(
    connection,
    paste(
      "SELECT table_schema, table_name FROM information_schema.tables",
      "WHERE table_catalog = 'staged' AND table_type = 'BASE TABLE'"
    )
  )
  if (nrow(merge(required, present, by = c("table_schema", "table_name"))) != nrow(required)) {
    stop("DuckVEP model does not preserve every required Ensembl source table", call. = FALSE)
  }

  registered_sources <- duckhts_bench_duckvep_sources()
  expected_manifest <- registered_sources[order(registered_sources$id), c(
    "id", "role", "release", "locator", "access", "cache_relpath", "transform", "supplier_identity"
  )]
  names(expected_manifest) <- c(
    "artifact_id", "role", "release", "source_locator", "access", "cache_relpath", "transform", "supplier_identity"
  )
  row.names(expected_manifest) <- NULL
  stored_manifest <- DBI::dbGetQuery(connection, paste(
    "SELECT artifact_id, role, release, source_locator, access, cache_relpath, transform, supplier_identity",
    "FROM staged.main.model_source_manifest ORDER BY artifact_id"
  ))
  if (!identical(expected_manifest, stored_manifest)) {
    stop("DuckVEP model source manifest does not match the registered Ensembl artifacts", call. = FALSE)
  }

  DBI::dbExecute(connection, "CREATE TEMP VIEW model_regions AS SELECT * FROM staged.main.model_regions")
  DBI::dbExecute(connection, "CREATE TEMP VIEW model_transcripts AS SELECT * FROM staged.main.model_transcripts")
  DBI::dbExecute(connection, paste(
    "CREATE TEMP VIEW duckvep_regulation_features AS",
    "SELECT * FROM staged.main.duckvep_regulation_features"
  ))
  stored <- DBI::dbGetQuery(connection, "SELECT * FROM staged.main.model_receipt")
  rebuilt <- DBI::dbGetQuery(connection, sprintf(
    paste0(
      "SELECT * FROM duckvep_model_receipt(",
      "'model_regions', 'model_transcripts', 'Ensembl', '116', 'GRCh38', %s, %s, %s, ",
      "regulation_features_table := 'duckvep_regulation_features')"
    ),
    duckhts_bench_duckvep_sql_literal(metadata$source_manifest_sha256[[1L]]),
    duckhts_bench_duckvep_sql_literal(metadata$reference_sha256[[1L]]),
    duckhts_bench_duckvep_sql_literal(metadata$transcript_filter[[1L]])
  ))
  if (!duckhts_bench_duckvep_compare_receipts(stored, rebuilt)) {
    stop("DuckVEP model relations do not reproduce the stored deterministic receipt", call. = FALSE)
  }

  registry <- duckhts_bench_registry()
  model_row <- registry[registry$id == "duckvep_ensembl116_model", , drop = FALSE]
  identity <- duckhts_bench_duckvep_identity(model_row$supplier_identity[[1L]])
  observed <- c(
    model_sha256 = as.character(stored$model_sha256[[1L]]),
    source_manifest_sha256 = as.character(stored$source_manifest_sha256[[1L]]),
    reference_sha256 = as.character(stored$reference_sha256[[1L]])
  )
  declared <- intersect(names(identity), names(observed))
  if (nrow(model_row) != 1L || any(observed[declared] != identity[declared])) {
    stop("DuckVEP model receipt does not match the registered logical identity", call. = FALSE)
  }
  invisible(stored)
}

#' Stage the public Ensembl 116 GRCh38 DuckVEP model.
#'
#' The producer imports the checksum-pinned Ensembl MySQL dump relations named by
#' the benchmark registry, compiles them with DuckHTS's DuckVEP model macros and
#' the matching primary-assembly FASTA, and verifies the deterministic receipt.
#'
#' @param extension Path to a locally built DuckHTS extension.
#' @param output Registered model path.
#' @return The staged model path, invisibly.
#' @export
duckhts_bench_stage_duckvep_ensembl116_model <- function(
  extension,
  output = duckhts_bench_artifact_path("duckvep_ensembl116_model")
) {
  if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("duckdb", quietly = TRUE)) {
    stop("staging the DuckVEP model requires the installed R packages DBI and duckdb", call. = FALSE)
  }
  extension <- normalizePath(extension, mustWork = TRUE)
  rows <- duckhts_bench_duckvep_sources()
  for (id in rows$id) duckhts_bench_fetch(id)
  manifests <- duckhts_bench_duckvep_validate_manifests(rows)

  if (file.exists(output) && file.info(output)$size > 0L) {
    duckhts_bench_validate_duckvep_ensembl116_model(
      output, extension, manifests$source_manifest_sha256
    )
    duckhts_bench_write_provenance("duckvep_ensembl116_model", output)
    return(invisible(output))
  }
  unlink(c(output, paste0(output, ".wal")), force = TRUE)

  fasta <- duckhts_bench_artifact_path("ensembl116_grch38_fasta_fa")
  if (!file.exists(fasta) || !file.info(fasta)$size || !file.exists(paste0(fasta, ".fai"))) {
    stop("registered Ensembl 116 GRCh38 FASTA and index must be staged before the model", call. = FALSE)
  }

  schema_definitions <- list()
  for (group in c("core", "funcgen")) {
    selected <- rows[rows$source_group == group, , drop = FALSE]
    schema_row <- selected[selected$source_kind == "source_schema", , drop = FALSE]
    tables <- selected$source_table[selected$source_kind == "source_table"]
    schema_definitions[[group]] <- duckhts_bench_duckvep_schema(
      duckhts_bench_artifact_path(schema_row$id[[1L]]),
      schema_row$source_database[[1L]],
      tables
    )
  }

  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(output, ".partial-", Sys.getpid())
  unlink(c(temporary, paste0(temporary, ".wal")), force = TRUE)
  success <- FALSE
  on.exit(if (!success) unlink(c(temporary, paste0(temporary, ".wal")), force = TRUE), add = TRUE)

  driver <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  connection <- DBI::dbConnect(driver, dbdir = temporary)
  connected <- TRUE
  on.exit(if (connected) DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(connection, "SET threads = 1")
  DBI::dbExecute(connection, paste("LOAD", duckhts_bench_duckvep_sql_literal(extension)))
  DBI::dbExecute(connection, "CREATE SCHEMA ensembl_core")
  DBI::dbExecute(connection, "CREATE SCHEMA ensembl_funcgen")

  for (group in c("core", "funcgen")) {
    schema <- paste0("ensembl_", group)
    selected <- rows[rows$source_group == group & rows$source_kind == "source_table", , drop = FALSE]
    selected <- selected[order(selected$source_table), , drop = FALSE]
    for (index in seq_len(nrow(selected))) {
      table <- selected$source_table[[index]]
      duckhts_bench_duckvep_import_table(
        connection, schema, table, schema_definitions[[group]][[table]],
        duckhts_bench_artifact_path(selected$id[[index]])
      )
    }
  }

  duckhts_bench_duckvep_require_metadata(connection, "ensembl_core", "schema_version", "116")
  duckhts_bench_duckvep_require_metadata(connection, "ensembl_core", "assembly.default", "GRCh38", 1L)
  duckhts_bench_duckvep_require_metadata(connection, "ensembl_core", "species.production_name", "homo_sapiens", 1L)
  duckhts_bench_duckvep_require_metadata(connection, "ensembl_core", "species.scientific_name", "Homo sapiens", 1L)
  duckhts_bench_duckvep_require_metadata(connection, "ensembl_funcgen", "schema_version", "116")
  duckhts_bench_duckvep_require_metadata(connection, "ensembl_funcgen", "species.production_name", "homo_sapiens", 1L)
  assembly_count <- DBI::dbGetQuery(
    connection,
    "SELECT count(*) AS n FROM ensembl_core.coord_system WHERE species_id = 1 AND version = 'GRCh38'"
  )$n[[1L]]
  if (assembly_count < 1) stop("Ensembl core dump lacks species 1 GRCh38 coordinate systems", call. = FALSE)

  DBI::dbExecute(connection, sprintf(
    paste0(
      "CREATE TABLE reference_chunks AS SELECT chrom, start, \"end\", seq ",
      "FROM fasta_nuc(%s, bin_width := 1048576, include_seq := true)"
    ),
    duckhts_bench_duckvep_sql_literal(fasta)
  ))
  DBI::dbExecute(
    connection,
    "CREATE TABLE model_regions AS SELECT * FROM duckvep_ensembl_regions('ensembl_core', 'reference_chunks', 'GRCh38', species_id := 1)"
  )
  DBI::dbExecute(connection, paste(
    "CREATE TABLE model_reference_hashes AS",
    "SELECT chrom, sha256(string_agg(seq, '' ORDER BY start)) AS sequence_sha256",
    "FROM reference_chunks GROUP BY chrom ORDER BY chrom"
  ))
  reference_sha256 <- DBI::dbGetQuery(
    connection,
    "SELECT sha256(string_agg(sequence_sha256, '' ORDER BY chrom)) AS reference_sha256 FROM model_reference_hashes"
  )$reference_sha256[[1L]]

  DBI::dbExecute(
    connection,
    "CREATE TABLE model_transcripts AS SELECT * FROM duckvep_ensembl_transcripts('ensembl_core', 'reference_chunks', 'GRCh38', species_id := 1)"
  )
  DBI::dbExecute(
    connection,
    "CREATE TABLE duckvep_regulation_features AS SELECT * FROM duckvep_ensembl_regulation_features('ensembl_funcgen', 'model_regions')"
  )

  transcript_filter <- "VEP 116 core selection on Ensembl 116 homo_sapiens GRCh38 primary-assembly FASTA"
  metadata <- data.frame(
    source_name = "Ensembl", source_release = "116", vep_release = "116",
    species = "homo_sapiens", species_id = 1L, assembly = "GRCh38",
    core_database = "homo_sapiens_core_116_38",
    funcgen_database = "homo_sapiens_funcgen_116_38",
    transcript_filter = transcript_filter,
    source_manifest_sha256 = manifests$source_manifest_sha256,
    reference_sha256 = reference_sha256,
    stringsAsFactors = FALSE
  )
  DBI::dbWriteTable(connection, "model_metadata", metadata)

  manifest_rows <- rows[order(rows$id), c(
    "id", "role", "release", "locator", "access", "cache_relpath", "transform", "supplier_identity"
  )]
  names(manifest_rows) <- c(
    "artifact_id", "role", "release", "source_locator", "access", "cache_relpath", "transform", "supplier_identity"
  )
  manifest_rows$source_database <- rows$source_database[match(manifest_rows$artifact_id, rows$id)]
  manifest_rows$source_manifest_sha256 <- manifests$database_manifest_sha256[manifest_rows$source_database]
  DBI::dbWriteTable(connection, "model_source_manifest", manifest_rows)

  DBI::dbExecute(connection, sprintf(
    paste0(
      "CREATE TABLE model_receipt AS SELECT * FROM duckvep_model_receipt(",
      "'model_regions', 'model_transcripts', 'Ensembl', '116', 'GRCh38', %s, %s, %s, ",
      "regulation_features_table := 'duckvep_regulation_features')"
    ),
    duckhts_bench_duckvep_sql_literal(manifests$source_manifest_sha256),
    duckhts_bench_duckvep_sql_literal(reference_sha256),
    duckhts_bench_duckvep_sql_literal(transcript_filter)
  ))
  DBI::dbExecute(connection, "CHECKPOINT")
  DBI::dbDisconnect(connection, shutdown = TRUE)
  connected <- FALSE

  duckhts_bench_validate_duckvep_ensembl116_model(
    temporary, extension, manifests$source_manifest_sha256
  )
  if (file.exists(output) && unlink(output, force = TRUE) != 0L) {
    stop("could not replace staged DuckVEP model", call. = FALSE)
  }
  if (!file.rename(temporary, output)) stop("could not publish staged DuckVEP model", call. = FALSE)
  success <- TRUE
  duckhts_bench_write_provenance("duckvep_ensembl116_model", output)
  invisible(output)
}
