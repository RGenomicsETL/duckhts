#!/usr/bin/env Rscript
# One invocation loads exactly one extension. See the companion Rmd for commands.
suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

materialize_input <- function(con, sql, oracles, input) {
  dbExecute(con, paste("CREATE OR REPLACE TABLE records AS", sql))
  dbExecute(con, "ALTER TABLE records ADD COLUMN text VARCHAR")
  dbExecute(con, paste0(
    "UPDATE records SET text = CASE WHEN packed IS NULL THEN NULL ",
    "WHEN len(packed) = 0 THEN '*' ELSE array_to_string(list_transform(packed, ",
    "lambda c: (c >> 4)::VARCHAR || substr('MIDNSHP=X', (c & 15) + 1, 1)), '') END"
  ))
  denominator <- dbGetQuery(con, paste0(
    "SELECT count(*)::VARCHAR AS records, ",
    "COALESCE(sum(len(packed)), 0)::VARCHAR AS operations, ",
    "COALESCE(sum(len(list_filter(packed, lambda c: (c & 15) IN (0, 7, 8)))), 0)::VARCHAR AS blocks, ",
    "count(*) FILTER (WHERE packed IS NULL)::VARCHAR AS null_cigars, ",
    "count(*) FILTER (WHERE len(packed) = 0)::VARCHAR AS empty_cigars, ",
    "max(len(packed))::VARCHAR AS max_ops_per_record, ",
    "hex(bit_xor(hash(rid, QNAME, FLAG, RNAME, POS, packed, text))) AS input_fingerprint ",
    "FROM records"
  ))
  if (as.numeric(denominator$records) == 0) stop("Empty workload: ", input)
  invalid <- dbGetQuery(con, paste0(
    "SELECT count(*) AS n FROM records WHERE len(list_filter(packed, ",
    "lambda c: c IS NULL OR (c & 15) > 8 OR (c >> 4) = 0)) > 0"
  ))$n
  if (invalid != 0) stop("Unsupported packed input; no records were removed")
  oracle_columns <- paste(sprintf("%s AS %s", oracles, names(oracles)), collapse = ", ")
  dbExecute(con, paste("CREATE OR REPLACE TABLE expected AS SELECT rid,", oracle_columns,
                       "FROM records"))
  denominator
}

run_case <- function(con, input, workload, representation, mode, repetitions, mismatch_path) {
  run <- function(sql) dbGetQuery(con, sql)
  strict <- switch(mode, default = "", false = ", false", true = ", true")
  expression <- if (workload %in% c("has_op_M", "has_op_P")) {
    sprintf("cigar_has_op(%s, '%s'%s)", representation,
            if (workload == "has_op_M") "M" else "P", strict)
  } else if (workload == "cigar_aligned_blocks") {
    sprintf("cigar_aligned_blocks(%s, POS%s)", representation, strict)
  } else {
    sprintf("%s(%s%s)", workload, representation, strict)
  }
  # Exact row comparisons include NULL and all block lists. On failure retain
  # the first counterexamples in the declared output namespace, then stop.
  mismatches <- run(sprintf(paste0(
    "SELECT rid, QNAME, FLAG, RNAME, POS, packed::VARCHAR AS packed, text, ",
    "(%s)::VARCHAR AS observed, e.%s::VARCHAR AS expected ",
    "FROM records JOIN expected e USING (rid) ",
    "WHERE (%s) IS DISTINCT FROM e.%s LIMIT 10"
  ), expression, workload, expression, workload))
  if (nrow(mismatches)) {
    dir.create(dirname(mismatch_path), recursive = TRUE, showWarnings = FALSE)
    write.csv(cbind(input, workload, representation, mode, mismatches),
              mismatch_path, row.names = FALSE)
    stop("Oracle disagreement: ", input, " / ", workload, " / ", representation, " / ", mode)
  }
  # Both summaries consume the full output and bind its hash to physical row
  # identity. Two reductions avoid cancellation of duplicate biological names.
  is_blocks <- workload == "cigar_aligned_blocks"
  total <- if (is_blocks) "COALESCE(len(value.width), 0)" else "value::BIGINT"
  summary_sql <- function(value, relation) sprintf(paste0(
    "SELECT count(*)::VARCHAR AS output_records, ",
    "count(value)::VARCHAR AS nonnull_outputs, ",
    "COALESCE(sum(%s), 0)::VARCHAR AS output_total, ",
    "hex(bit_xor(hash(rid, QNAME, FLAG, RNAME, POS, value))) AS fingerprint_xor, ",
    "sum(hash(rid, QNAME, FLAG, RNAME, POS, value)::HUGEINT)::VARCHAR AS fingerprint_sum ",
    "FROM (SELECT rid, QNAME, FLAG, RNAME, POS, %s AS value FROM %s)"
  ), total, value, relation)
  expected <- run(summary_sql(paste0("e.", workload), "records JOIN expected e USING (rid)"))
  query <- summary_sql(expression, "records")
  # The full query warms the helper and materialized data; this pass is untimed.
  if (!identical(run(query), expected)) stop("Warm-up summary disagreement")
  seconds <- numeric(repetitions)
  for (repeat_id in seq_len(repetitions)) {
    started <- proc.time()[["elapsed"]]
    observed <- run(query)
    seconds[[repeat_id]] <- proc.time()[["elapsed"]] - started
    if (!identical(observed, expected)) stop("Timed output disagreement")
  }
  key <- data.frame(input = input, workload = workload,
                    representation = representation, mode = mode)
  list(
    summary = cbind(key, expected, repetitions = repetitions, min_seconds = min(seconds),
                    median_seconds = median(seconds), max_seconds = max(seconds)),
    timings = cbind(key, repeat_id = seq_len(repetitions), seconds)
  )
}

main <- function(args) {
  if (length(args) != 9L) {
    stop(paste(
      "Usage: Rscript scripts/benchmark_cigar_validation.R",
      "baseline|candidate EXTENSION SOURCE_DIR REVISION SRC_TREE clean|uncommitted",
      "OUTPUT_PREFIX full|smoke REPEATS"
    ), call. = FALSE)
  }
  role <- match.arg(args[[1L]], c("baseline", "candidate"))
  extension <- normalizePath(args[[2L]], mustWork = TRUE)
  source_dir <- normalizePath(args[[3L]], mustWork = TRUE)
  revision <- args[[4L]]
  src_tree <- args[[5L]]
  source_state <- match.arg(args[[6L]], c("clean", "uncommitted"))
  prefix <- args[[7L]]
  size <- match.arg(args[[8L]], c("full", "smoke"))
  repetitions <- as.numeric(args[[9L]])
  if (!repetitions %in% 1:20) stop("REPEATS must be an integer in 1..20")
  if (!all(grepl("^[0-9a-f]{40}$", c(revision, src_tree)))) {
    stop("REVISION and SRC_TREE must be full Git object IDs supplied by the caller")
  }

  # Preserve the caller's affinity; use the same taskset mask for both invocations.
  affinity <- sub("^Cpus_allowed_list:[[:space:]]*", "",
                  grep("^Cpus_allowed_list:", readLines("/proc/self/status"), value = TRUE))
  if (length(affinity) != 1L) stop("Cannot record Linux CPU affinity")
  cpu <- sub("^model name[[:space:]]*:[[:space:]]*", "",
             grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[[1L]])
  sha256 <- function(path) digest::digest(file = path, algo = "sha256")
  extension_sha256 <- sha256(extension)
  con <- dbConnect(duckdb(shared_home = FALSE,
                          config = list(allow_unsigned_extensions = "true", threads = "1")))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  invisible(dbExecute(con, paste("LOAD", dbQuoteString(con, extension))))
  run <- function(sql) dbGetQuery(con, sql)

  # Resolve only cached artifacts. Missing inputs fail before any timing or output write.
  bam <- duckhtsbench::duckhts_bench_artifact_path("ont_ecoli_k12_bam")
  receipt_path <- paste0(bam, ".provenance.tsv")
  if (!all(file.exists(c(bam, receipt_path)))) stop("ONT BAM and receipt must be staged")
  receipt <- read.delim(receipt_path, colClasses = "character", quote = "", comment.char = "")
  receipt_value <- function(field) {
    value <- receipt$value[receipt$field == field]
    if (length(value) != 1L) stop("Missing ONT receipt field: ", field)
    value
  }
  bam_sha256 <- sha256(bam)
  limit <- c(full = "", smoke = " LIMIT 256")[[size]]
  ont_sql <- sprintf(paste0(
    "SELECT row_number() OVER ()::BIGINT AS rid, QNAME, FLAG, RNAME, POS, CIGAR AS packed ",
    "FROM read_bam(%s, scan_mode := 'sequential', cigar_representation := 'binary')%s"
  ), dbQuoteString(con, bam), limit)

  # BAM op codes follow third_party/htslib/htslib/sam.h: MIDNSHP=X = 0..8.
  # Vary every record's lengths and endpoints; store columns before invoking helpers.
  # M is the first aligned op, P is absent, and H/S endpoints alternate by record.
  rows <- c(full = 8192L, smoke = 4096L)[[size]]
  cycles <- c(full = 32L, smoke = 8L)[[size]]
  synthetic_sql <- sprintf(paste0(
    "SELECT i + 1 AS rid, 'synthetic-' || i AS QNAME, 0::INTEGER AS FLAG, ",
    "'literal' AS RNAME, (i * 17)::BIGINT AS POS, ",
    "list_concat([(((i %% 7 + 1) << 4) | CASE WHEN i %% 2 = 0 THEN 4 ELSE 5 END)::UINTEGER], ",
    "flatten(list_transform(range(0, %d), lambda j: ",
    "[((i %% 31 + j + 1) << 4)::UINTEGER, 17::UINTEGER, 34::UINTEGER, ",
    "51::UINTEGER, 71::UINTEGER, 24::UINTEGER])), ",
    "[(((i %% 5 + 1) << 4) | CASE WHEN i %% 2 = 0 THEN 5 ELSE 4 END)::UINTEGER]) AS packed ",
    "FROM range(%d) t(i)"
  ), cycles, rows)

  # Independent packed-op sums and literal first/last-op soft clipping rules.
  op_sum <- function(codes) sprintf(paste0(
    "COALESCE(list_sum(list_transform(packed, lambda c: ",
    "CASE WHEN (c & 15) IN (%s) THEN (c >> 4)::BIGINT ELSE 0 END)), 0)::BIGINT"
  ), codes)
  metric_oracles <- c(
    cigar_has_soft_clip = "list_contains(list_transform(packed, lambda c: c & 15), 4)",
    cigar_has_hard_clip = "list_contains(list_transform(packed, lambda c: c & 15), 5)",
    cigar_left_soft_clip = "CASE WHEN (packed[1] & 15) = 4 THEN (packed[1] >> 4)::BIGINT ELSE 0 END",
    cigar_right_soft_clip = "CASE WHEN (packed[-1] & 15) = 4 THEN (packed[-1] >> 4)::BIGINT ELSE 0 END",
    cigar_query_length = op_sum("0, 1, 4, 7, 8"),
    cigar_aligned_query_length = op_sum("0, 7, 8"),
    cigar_reference_length = op_sum("0, 2, 3, 7, 8")
  )
  metric_oracles <- paste0("CASE WHEN len(packed) > 0 THEN ", metric_oracles, " END") |>
    setNames(names(metric_oracles))
  # Prefix-slice geometry from benchmark_cigar_aligned_blocks.Rmd. It is evaluated
  # once outside timers and does not call any native CIGAR helper.
  ref_starts <- paste0(
    "list_transform(list_filter(list_transform(packed, lambda c, i: {'kind': c & 15, ",
    "'v': POS + COALESCE(list_sum(list_transform(list_slice(packed, 1, i - 1), ",
    "lambda p: CASE WHEN (p & 15) IN (0, 2, 3, 7, 8) ",
    "THEN (p >> 4)::BIGINT ELSE 0 END)), 0)::BIGINT}), ",
    "lambda x: x.kind IN (0, 7, 8)), lambda x: x.v)"
  )
  query_starts <- paste0(
    "list_transform(list_filter(list_transform(packed, lambda c, i: {'kind': c & 15, ",
    "'v': COALESCE(list_sum(list_transform(list_slice(packed, 1, i - 1), ",
    "lambda p: CASE WHEN (p & 15) IN (0, 1, 4, 7, 8) ",
    "THEN (p >> 4)::BIGINT ELSE 0 END)), 0)::BIGINT}), ",
    "lambda x: x.kind IN (0, 7, 8)), lambda x: x.v)"
  )
  widths <- paste0("list_transform(list_filter(packed, lambda c: (c & 15) IN (0, 7, 8)), ",
                   "lambda c: (c >> 4)::BIGINT)")
  block_oracle <- sprintf(paste0(
    "CASE WHEN len(packed) > 0 THEN ",
    "{'ref_start': %s, 'query_start': %s, 'width': %s} END"
  ), ref_starts, query_starts, widths)
  oracles <- c(metric_oracles,
               has_op_M = "list_contains(list_transform(packed, lambda c: c & 15), 0)",
               has_op_P = "list_contains(list_transform(packed, lambda c: c & 15), 6)",
               cigar_aligned_blocks = block_oracle)
  modes <- list(baseline = "default", candidate = c("default", "false", "true"))[[role]]
  cases <- expand.grid(workload = names(oracles), representation = c("text", "packed"),
                       mode = modes, stringsAsFactors = FALSE)
  summaries <- timings <- inputs <- list()

  input_queries <- c(ont = ont_sql, long_cigar = synthetic_sql)
  for (input in names(input_queries)) {
    message(role, ": materialize and validate ", input)
    denominator <- materialize_input(con, input_queries[[input]], oracles, input)
    inputs[[input]] <- cbind(data.frame(input = input, size = size), denominator)

    for (case in seq_len(nrow(cases))) {
      result <- run_case(
        con, input, cases$workload[[case]], cases$representation[[case]], cases$mode[[case]],
        repetitions, paste0(prefix, "_", role, "_mismatches.csv")
      )
      summaries[[length(summaries) + 1L]] <- result$summary
      timings[[length(timings) + 1L]] <- result$timings
    }
  }
  if (!identical(c(sha256(extension), sha256(bam)), c(extension_sha256, bam_sha256))) {
    stop("Extension or ONT input changed during the run")
  }
  metadata <- c(
    role = role, source_dir = source_dir, revision = revision, src_tree = src_tree,
    source_state = source_state, extension_sha256 = extension_sha256,
    run_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), host = Sys.info()[["nodename"]],
    system = paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "), cpu = cpu,
    affinity = affinity, threads = "1", r_version = R.version.string,
    duckdb_version = run("SELECT version() AS v")$v,
    htslib_version = run("SELECT duckhts_htslib_version() AS v")$v,
    driver_sha256 = sha256("scripts/benchmark_cigar_validation.R"),
    ont_artifact = "ont_ecoli_k12_bam", ont_sha256 = bam_sha256,
    ont_release = receipt_value("release"), ont_sources = receipt_value("source_locator"),
    ont_transform = receipt_value("transform"), ont_aligner = receipt_value("aligner"),
    ont_sorter = receipt_value("sorter"), size = size,
    synthetic_rows = as.character(rows), synthetic_cycles = as.character(cycles),
    riker = if (file.exists(duckhtsbench::duckhts_bench_artifact_path("riker_hg00188_bam"))) {
      "Not measured by this focused workload"
    } else "Not staged; not measured (approximately 17 GB acquisition)"
  )
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  outputs <- list(metadata = data.frame(property = names(metadata), value = unname(metadata)),
                  inputs = do.call(rbind, inputs), results = do.call(rbind, summaries),
                  timings = do.call(rbind, timings))
  for (name in names(outputs)) {
    write.csv(outputs[[name]], paste0(prefix, "_", role, "_", name, ".csv"), row.names = FALSE)
  }
  message(role, ": all oracle and timed-output checks passed; wrote ", prefix, "_", role, "_*.csv")
}

main(commandArgs(trailingOnly = TRUE))
