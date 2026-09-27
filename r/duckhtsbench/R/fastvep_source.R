# Resolve retained evidence only. Rendering never derives or repairs an object.
duckhts_bench_read_fastvep_source_map <- function(expected_sha256 = NULL,
    expected_input_sha256 = NULL, expected_extension_sha256 = NULL) {
  id <- "fastvep_giab_hg002_v421_source_map"
  row <- duckhts_bench_registry()
  row <- row[row$id == id, , drop = FALSE]
  if (nrow(row) != 1L || row$locator != "artifact:variantkey_giab_hg002_v421" ||
      row$transform != "map_physical_alt_ordinals") {
    stop("registered GIAB source-map derivation is missing or different", call. = FALSE)
  }
  declared <- duckhts_bench_identity_fields(row$supplier_identity)
  output <- file.path(duckhts_bench_artifact_path(id), "source_alleles.parquet")
  result <- c(map = output, receipt = paste0(output, ".provenance.tsv"))
  if (!all(file.exists(result))) stop("registered GIAB source-map bundle is incomplete", call. = FALSE)
  receipt <- utils::read.delim(result[["receipt"]], colClasses = "character",
    quote = "", comment.char = "", check.names = FALSE)
  required <- c("artifact_id", "transform", "supplier_identity", "source_id", "input_sha256",
    "generator_sha256", "generator_source_sha256", "extension_sha256", "duckdb_version",
    "source_map_sha256", "records", "alleles", "eligible_alleles")
  if (!identical(names(receipt), c("field", "value")) || anyDuplicated(receipt$field) ||
      !all(required %in% receipt$field) || anyNA(receipt$value) || any(!nzchar(receipt$value))) {
    stop("invalid GIAB source-map receipt", call. = FALSE)
  }
  values <- stats::setNames(receipt$value, receipt$field)
  digest_fields <- c("input_sha256", "generator_sha256", "generator_source_sha256",
    "extension_sha256", "source_map_sha256")
  count_fields <- c("records", "alleles", "eligible_alleles")
  if (values[["artifact_id"]] != id || values[["source_id"]] != "variantkey_giab_hg002_v421" ||
      values[["transform"]] != row$transform || values[["supplier_identity"]] != row$supplier_identity ||
      !all(c("schema", "input_sha256", count_fields) %in% names(declared)) ||
      declared[["schema"]] != "physical_alt_source_v1" ||
      any(!grepl("^[0-9a-f]{64}$", values[digest_fields])) ||
      any(!grepl("^[1-9][0-9]*$", values[count_fields])) ||
      !identical(unname(values[c("input_sha256", count_fields)]),
        unname(declared[c("input_sha256", count_fields)]))) {
    stop("GIAB source-map receipt differs from its registered identity", call. = FALSE)
  }
  for (expected in list(expected_sha256, expected_input_sha256, expected_extension_sha256)) {
    if (!is.null(expected) && (length(expected) != 1L || is.na(expected) ||
        !grepl("^[0-9a-f]{64}$", expected))) stop("expected source-map digests must be SHA256 values", call. = FALSE)
  }
  if (!is.null(expected_extension_sha256) &&
      expected_extension_sha256 != values[["extension_sha256"]]) {
    stop("retained GIAB source-map extension differs from the expected digest", call. = FALSE)
  }
  if ((!is.null(expected_sha256) && expected_sha256 != values[["source_map_sha256"]]) ||
      (!is.null(expected_input_sha256) && expected_input_sha256 != values[["input_sha256"]]) ||
      !identical(duckhts_bench_duckvep_sha256_file(output), values[["source_map_sha256"]])) {
    stop("retained GIAB source map differs from the expected digest", call. = FALSE)
  }
  attr(result, "identity") <- values
  result
}

# The generator identity covers the two source-map functions, not unrelated
# annotation projections that share their R file. The complete file is also receipted.
duckhts_bench_stage_fastvep_source_map <- function(repo, extension,
    memory_limit = "4GB", max_spill = "8GB") {
  id <- "fastvep_giab_hg002_v421_source_map"
  source_id <- "variantkey_giab_hg002_v421"
  repo <- normalizePath(repo, mustWork = TRUE)
  generator_path <- file.path(repo, "benchmarks/benchmark_duckvep_fastvep_fields.R")
  generator <- new.env(parent = baseenv())
  source(generator_path, local = generator)
  functions <- c("duckvep_fastvep_prepare_source", "duckvep_fastvep_write_source_map")
  code <- vapply(functions, function(name) {
    if (!is.function(generator[[name]])) stop("source-map generator is missing: ", name, call. = FALSE)
    paste(name, paste(deparse(generator[[name]], width.cutoff = 500L), collapse = "\n"), sep = "\n")
  }, character(1L))
  generator_hash <- duckhts_bench_duckvep_sha256_text(paste(code, collapse = "\n"))
  destination <- duckhts_bench_artifact_path(id)
  registry <- duckhts_bench_registry()
  row <- registry[registry$id == id, , drop = FALSE]
  declared <- duckhts_bench_identity_fields(row$supplier_identity)
  if (nrow(row) != 1L || row$locator != paste0("artifact:", source_id) ||
      row$transform != "map_physical_alt_ordinals" || declared[["schema"]] != "physical_alt_source_v1") {
    stop("registered GIAB source-map derivation is missing or different", call. = FALSE)
  }
  input <- duckhts_bench_artifact_path(source_id)
  if (!file.exists(input)) stop("registered GIAB source must be staged before its map", call. = FALSE)
  duckhts_bench_validate_identity(source_id, input)
  extension <- normalizePath(extension, mustWork = TRUE)
  hash <- duckhts_bench_duckvep_sha256_file
  inputs <- c(input_sha256 = input, generator_source_sha256 = generator_path, extension_sha256 = extension)
  hashes <- vapply(inputs, hash, character(1L))
  if (hashes[["input_sha256"]] != declared[["input_sha256"]]) {
    stop("GIAB source bytes differ from the registered source-map input", call. = FALSE)
  }
  if (file.exists(destination)) {
    result <- duckhts_bench_read_fastvep_source_map(
      expected_input_sha256 = hashes[["input_sha256"]],
      expected_extension_sha256 = hashes[["extension_sha256"]])
    if (attr(result, "identity")[["generator_sha256"]] != generator_hash) {
      stop("retained GIAB source-map generator differs; existing bundle preserved", call. = FALSE)
    }
    return(result)
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  staging <- tempfile(".source-map-", tmpdir = dirname(destination))
  if (!dir.create(staging)) stop("could not create source-map staging directory", call. = FALSE)
  on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE), add = TRUE)
  output <- file.path(staging, "source_alleles.parquet")
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
  connected <- TRUE
  on.exit(if (connected) DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE, after = FALSE)
  q <- function(value) as.character(DBI::dbQuoteString(con, value))
  DBI::dbExecute(con, paste("LOAD", q(extension)))
  DBI::dbExecute(con, "SET threads = 1")
  DBI::dbExecute(con, paste("SET memory_limit =", q(memory_limit)))
  DBI::dbExecute(con, paste("SET max_temp_directory_size =", q(max_spill)))
  DBI::dbExecute(con, paste("SET temp_directory =", q(file.path(staging, "spill"))))
  counts <- DBI::dbGetQuery(con, paste0("SELECT count(*)::VARCHAR AS records,
    sum(len(ALT))::VARCHAR AS alleles,
    sum(len(list_filter(ALT, lambda a: regexp_full_match(REF, '[ACGTNacgtn]+')
      AND regexp_full_match(a, '[ACGTNacgtn]+') AND upper(REF) <> upper(a))))::VARCHAR
      AS eligible_alleles FROM read_bcf(", q(input),
    ", scan_mode := 'sequential', decompression_threads := 0)"))
  count_values <- unlist(counts, use.names = TRUE)
  if (!identical(unname(count_values), unname(declared[names(count_values)]))) {
    stop("raw GIAB counts differ from the registered source-map denominator", call. = FALSE)
  }
  generator$duckvep_fastvep_write_source_map(con, input, output)
  mapped <- DBI::dbGetQuery(con, paste0("SELECT count(DISTINCT record_index)::VARCHAR AS records,
    count(*)::VARCHAR AS alleles,
    count(*) FILTER (WHERE eligible)::VARCHAR AS eligible_alleles FROM read_parquet(", q(output), ")"))
  if (!identical(as.character(mapped[1, ]), as.character(counts[1, names(mapped)]))) {
    stop("GIAB source map does not reconcile to the raw input", call. = FALSE)
  }
  invalid <- DBI::dbGetQuery(con, paste0("SELECT count(*) AS n FROM (
    SELECT record_index, alt_index FROM read_parquet(", q(output), ") GROUP BY ALL
    HAVING count(*) != 1 OR record_index IS NULL OR alt_index IS NULL
      OR record_index = 0 OR alt_index = 0)"))$n
  if (invalid != 0) stop("GIAB source map has invalid physical ALT keys", call. = FALSE)
  DBI::dbDisconnect(con, shutdown = TRUE)
  connected <- FALSE
  if (!identical(hashes, vapply(inputs, hash, character(1L)))) {
    stop("GIAB source-map input or generator changed during staging", call. = FALSE)
  }
  receipt <- duckhts_bench_write_provenance(id, output)
  retained <- utils::read.delim(receipt, colClasses = "character", quote = "", comment.char = "")
  retained$value[retained$field == "cached_output"] <- file.path(destination, basename(output))
  evidence <- c(source_id = source_id, hashes, generator_sha256 = generator_hash,
    duckdb_version = as.character(utils::packageVersion("duckdb")),
    source_map_sha256 = hash(output), count_values)
  retained <- rbind(retained, data.frame(field = names(evidence), value = unname(evidence)))
  utils::write.table(retained, receipt, sep = "\t", row.names = FALSE, quote = FALSE)
  if (file.exists(destination) || !file.rename(staging, destination)) {
    stop("could not publish source-map bundle without replacing existing data", call. = FALSE)
  }
  duckhts_bench_read_fastvep_source_map(evidence[["source_map_sha256"]],
    evidence[["input_sha256"]], evidence[["extension_sha256"]])
}
