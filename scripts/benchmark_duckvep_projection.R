# Full SQL presentation materialization, one fresh process per model/worker count.
# Input staging and annotation are outside the timer; no fused comparator exists.
benchmark_duckvep_projection <- function(root, extension, output, extension_receipt = NULL,
                                        copies = 8L, iterations = 3L, threads = c(1L, 4L),
                                        cases = c("forward", "reverse", "three_exon_phase2")) {
  root <- normalizePath(root, mustWork = TRUE)
  extension <- normalizePath(extension, mustWork = TRUE)
  source(file.path(root, "scripts/duckvep_evidence.R"), local = TRUE)
  stopifnot(copies > 0L, copies == floor(copies), iterations > 0L,
    iterations == floor(iterations), all(threads > 0L), all(threads == floor(threads)))
  revision <- duckvep_evidence_revision(root)
  extension_hash <- duckvep_evidence_sha256(extension)
  binding <- "diagnostic_unbound"
  if (!is.null(extension_receipt)) {
    duckvep_evidence_assert_checkout(root, revision, allowed_outputs = output)
    binding <- duckvep_evidence_read_extension_receipt(extension_receipt, root,
      extension, revision)$binding
  }
  results <- list()
  for (case in cases) for (workers in threads) {
    result <- callr::r(function(root, extension, case, workers, copies, iterations) {
      source(file.path(root, "scripts/duckvep_evidence.R"), local = TRUE)
      source(file.path(root, "test/duckvep/conformance/projection_fixtures.R"), local = TRUE)
      inputs <- duckhtsbench::duckhts_bench_stage_repository_fixtures(root, "duckvep-projection")
      con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
      on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
      DBI::dbExecute(con, paste("LOAD", DBI::dbQuoteString(con, extension)))
      DBI::dbExecute(con, sprintf("SET threads=%d", workers))
      directory <- tempfile("projection-benchmark-")
      dir.create(directory)
      on.exit(unlink(directory, recursive = TRUE), add = TRUE)
      gff <- duckvep_projection_fixture(con, root, inputs, case, directory)
      generated <- file.path(directory, "generated.vcf")
      vcf <- file.path(directory, "events.vcf")
      duckvep_evidence_command("Rscript", c(file.path(root,
        "test/duckvep/conformance/generate_witnesses.R"), "--gff", gff,
        "--fasta", inputs[["projection_reference"]], "--ext", extension,
        "--out", generated, "--random-cases", "1000", "--seed", "173"), "projection workload generation")
      duckvep_projection_label_records(generated, vcf)
      duckvep_projection_input(con, vcf)
      original <- DBI::dbGetQuery(con, "SELECT count(*) n FROM projection_events")$n
      DBI::dbExecute(con, sprintf("CREATE OR REPLACE TABLE projection_events AS
        SELECT e.* REPLACE((event_index + %d*i)::UBIGINT AS event_index)
        FROM projection_events e, range(%d) r(i) ORDER BY position, event_index", original, copies))
      DBI::dbExecute(con, "CREATE OR REPLACE TABLE projection_annotations AS
        SELECT * FROM duckvep_annotate('projection_events', 'projection')")
      expanded <- DBI::dbGetQuery(con, "SELECT count(*) n FROM projection_annotations")$n
      materialize <- "CREATE OR REPLACE TEMP TABLE measured AS SELECT *
        FROM duckvep_transcript_projection('projection_events',
          'projection_annotations', 'projection_transcripts')"
      times <- numeric(iterations)
      fingerprint <- NULL
      for (pass in 0:iterations) {
        start <- proc.time()[["elapsed"]]
        DBI::dbExecute(con, materialize)
        elapsed <- proc.time()[["elapsed"]] - start
        # The full schema and typed row multiset, including duplicates and NULLs,
        # are checked after every pass outside the materialization timer.
        observed <- digest::digest(list(schema = DBI::dbGetQuery(con, "DESCRIBE measured"),
          rows = DBI::dbGetQuery(con, "SELECT * FROM measured ORDER BY ALL")),
          algo = "sha256", serializeVersion = 3)
        if (pass == 0L) fingerprint <- observed else {
          stopifnot(identical(observed, fingerprint))
          times[[pass]] <- elapsed
        }
      }
      denominators <- DBI::dbGetQuery(con, "SELECT count(*) output_rows,
        sum(octet_length(encode(to_json(m)))) json_utf8_bytes FROM measured m")
      stopifnot(denominators$output_rows == expanded)
      # Linux kernel process high-water mark includes setup and verification.
      # This is not a claim about DuckDB-only memory or materialization-only RSS.
      status <- readLines("/proc/self/status")
      peak_rss <- as.numeric(sub("^VmHWM:[[:space:]]*([0-9]+).*", "\\1",
        grep("^VmHWM:", status, value = TRUE))) * 1024
      stopifnot(length(peak_rss) == 1L, is.finite(peak_rss))
      cpu <- sub("^[^:]*:[[:space:]]*", "",
        grep("^model name[[:space:]]*:", readLines("/proc/cpuinfo"), value = TRUE)[[1L]])
      data.frame(case, threads = workers, seed = 173L, random_cases = 1000L,
        unique_input_records = original, copies, input_alleles = original * copies,
        expanded_transcript_rows = expanded, denominators, iterations,
        min_seconds = min(times), median_seconds = median(times), max_seconds = max(times),
        timings_seconds = paste(times, collapse = ";"), peak_rss_bytes = peak_rss,
        rss_scope = "Linux VmHWM: fresh R/DuckDB process including setup and full-output verification",
        output_sha256 = fingerprint, input_vcf_sha256 = duckvep_evidence_sha256(vcf),
        model_gff_sha256 = duckvep_evidence_sha256(gff),
        reference_sha256 = duckvep_evidence_sha256(inputs[["projection_reference"]]),
        duckdb_version = as.character(utils::packageVersion("duckdb")),
        machine = paste(Sys.info()[c("sysname", "release", "machine", "nodename")], collapse = " "), cpu)
    }, args = list(root = root, extension = extension, case = case, workers = workers,
      copies = copies, iterations = iterations))
    print(result[c("case", "threads", "input_alleles", "output_rows", "median_seconds", "peak_rss_bytes")])
    results[[length(results) + 1L]] <- result
  }
  results <- do.call(rbind, results)
  for (case in unique(results$case)) {
    rows <- results[results$case == case, ]
    stopifnot(length(unique(rows$output_sha256)) == 1L,
      length(unique(rows$input_vcf_sha256)) == 1L, length(unique(rows$output_rows)) == 1L)
  }
  stopifnot(identical(extension_hash, duckvep_evidence_sha256(extension)))
  if (!is.null(extension_receipt)) duckvep_evidence_assert_checkout(root, revision, allowed_outputs = output)
  results <- cbind(run_date = as.character(Sys.Date()), source_revision = revision,
    build_binding = binding, extension_sha256 = extension_hash, results)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  history <- if (file.exists(output)) rbind(
    utils::read.csv(output, stringsAsFactors = FALSE), results) else results
  utils::write.csv(history, output, row.names = FALSE)
  invisible(results)
}
