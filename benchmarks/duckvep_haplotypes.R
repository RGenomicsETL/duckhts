#!/usr/bin/env Rscript
# Two explicitly different scopes over one registered, factorized input:
# initialized native replay/count sink; full public SQL sort/replay/materialization.

fixture <- function(paths) {
  fasta <- readLines(paths[["haplotype_benchmark_reference"]])
  stopifnot(length(fasta) == 2L, nchar(fasta[2L]) == 180L)
  events <- read.delim(paths[["haplotype_benchmark_events"]])
  stopifnot(identical(events$event_index, 1:4), all(diff(events$position) > 0))
  masks <- c(3L, 12L, 15L)
  cds <- vapply(masks, function(mask) {
    sequence <- fasta[2L]
    for (i in rev(which(bitwAnd(mask, bitwShiftL(1L, 0:3)) != 0L))) {
      p <- events$position[i]
      r <- nchar(events$reference[i])
      stopifnot(substr(sequence, p, p + r - 1L) == events$reference[i])
      sequence <- paste0(substr(sequence, 1L, p - 1L), events$alternate[i],
        substring(sequence, p + r))
    }
    sequence
  }, "")
  protein <- as.character(Biostrings::translate(Biostrings::DNAStringSet(cds)))
  protein <- sub("(\\*).*", "\\1", protein)
  list(reference = fasta[2L], events = events,
    expected = data.frame(mask = masks, cds, protein))
}

metric_names <- c("seconds", "workspace_bytes", "model_bytes", "input_records",
  "projected_events", "input_calls", "peak_transcripts", "peak_carriers", "peak_prefixes",
  "prefixes_created", "output_leaves", "translated_bases", "output_carriers", "cds_bytes",
  "protein_bytes", "physical_edits", "contributors", "blocks", "peak_events",
  "peak_projections", "peak_allele_bytes")

native <- function(job, input) {
  dll <- dyn.load(job$shared)
  on.exit(dyn.unload(job$shared))
  symbol <- getNativeSymbolInfo("duckhts_bench_haplotype_stream", PACKAGE = dll)
  invoke <- function(expected = input$expected$cds) .C(symbol,
    input$reference, as.integer(job$transcripts), as.integer(job$samples), as.integer(job$overlap),
    input$events$position, input$events$reference, input$events$alternate,
    as.integer(t(as.matrix(input$events[5:12]))), expected, input$expected$protein,
    metrics = double(length(metric_names)), status = integer(1L))
  answer <- invoke()
  stopifnot(answer$status == 0L)
  result <- as.list(setNames(answer$metrics, metric_names))
  if (job$verify) {
    broken <- input$expected$cds
    broken[1L] <- paste0("C", substring(broken[1L], 2L))
    stopifnot(invoke(broken)$status != 0L)
  }
  stopifnot(result$input_records == job$transcripts / job$overlap * 4,
    result$projected_events == job$transcripts * 4,
    result$input_calls == job$transcripts * job$samples * 4,
    result$output_leaves == job$transcripts * 3,
    result$output_carriers == job$transcripts * job$samples * 7 / 4,
    result$translated_bases == job$transcripts * sum(nchar(input$expected$cds)),
    result$cds_bytes == result$translated_bases,
    result$protein_bytes == job$transcripts * sum(nchar(input$expected$protein)),
    result$physical_edits == job$transcripts * 8,
    result$contributors == job$transcripts * 8)
  result
}

sql <- function(job, input) {
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbExecute(con, paste("LOAD", q(job$extension)))
  DBI::dbExecute(con, "SET threads=1")
  DBI::dbExecute(con, "SET memory_limit='4GB'")
  DBI::dbWriteTable(con, "fixture_events", input$events)
  DBI::dbWriteTable(con, "expected_sequences", input$expected)
  DBI::dbExecute(con, sprintf("CREATE TABLE tx AS SELECT i::UINTEGER transcript_index,
    0::UINTEGER seq_region,(i//%d*1000+100)::UBIGINT transcript_start,
    (transcript_start+179)::UBIGINT transcript_end,1::TINYINT strand,i::UINTEGER gene_index,
    3::UBIGINT transcript_flags,transcript_start cds_start,transcript_end cds_end,
    %s::BLOB cds_sequence,1::UTINYINT codon_table FROM range(%d) t(i)",
    job$overlap, q(input$reference), job$transcripts))
  exons <- "SELECT transcript_index,transcript_start exon_start,transcript_end exon_end,
    1::UBIGINT exon_cdna_start,180::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase FROM tx"
  stopifnot(DBI::dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('bench',",
    q("SELECT 0::UINTEGER seq_region"), ",", q("SELECT * FROM tx"), ",", q(exons), ")"))$loaded)
  DBI::dbExecute(con, sprintf("CREATE TABLE calls AS SELECT
    (transcript_index//%d*4+e.event_index)::UBIGINT event_index,0::UINTEGER seq_region,
    (transcript_start+e.position-1)::UBIGINT AS position,reference,alternate,1::UINTEGER alt_index,
    transcript_index,s.i::UINTEGER sample_index,
    CASE s.i%%4 WHEN 0 THEN [g0_h1,g0_h2] WHEN 1 THEN [g1_h1,g1_h2]
      WHEN 2 THEN [g2_h1,g2_h2] ELSE [g3_h1,g3_h2] END::INTEGER[] alleles,
    [true,true] phase_before,10::BIGINT phase_set
    FROM tx CROSS JOIN fixture_events e CROSS JOIN range(%d) s(i)
    ORDER BY hash(event_index,transcript_index,sample_index)", job$overlap, job$samples))
  n <- DBI::dbGetQuery(con, "SELECT count(*) n FROM calls")$n
  stopifnot(n == job$transcripts * job$samples * 4)
  # Identical native pool capacities, plus separately bounded SQL alignment scratch.
  limits <- c(max_active_events = 5, max_active_transcripts = job$overlap,
    max_active_carriers = 2 * job$samples * job$overlap, max_active_prefixes = 16 * job$overlap,
    max_active_projections = 4 * job$overlap + 1, max_allele_bytes = 64,
    max_leaf_events = 4, max_leaf_edits = 4, max_sequence_bases = 184, max_ploidy = 2,
    max_phase_sets = 1, max_alignment_cells = 65536, max_leaf_differences = 64,
    workspace_limit = 64 * 1024^2)
  query <- paste0("CREATE OR REPLACE TABLE measured AS SELECT * FROM duckvep_haplotypes(",
    q("SELECT * FROM calls"), ",'bench',", paste(names(limits), ":=", limits, collapse = ","), ")")
  DBI::dbExecute(con, query) # Full warm-up, including the mandatory internal sort.
  start <- proc.time()[["elapsed"]]
  DBI::dbExecute(con, query)
  elapsed <- proc.time()[["elapsed"]] - start
  fingerprint <- DBI::dbGetQuery(con, "SELECT count(*) output_leaves,
    sum(carrier_count) output_carriers,sum(length(cds)) cds_bytes,sum(length(protein)) protein_bytes,
    sum(edit_count) physical_edits,sum(length(contributors)) contributors,sum(length(coding_blocks)) blocks,
    sum(octet_length(encode(to_json(h)))) json_bytes,
    bit_xor(hash(to_json(h)))::VARCHAR xor_hash,sum(hash(to_json(h))::HUGEINT)::VARCHAR sum_hash
    FROM (SELECT * REPLACE(list_sort(carriers) AS carriers,list_sort(contributors) AS contributors)
      FROM measured) h")
  stopifnot(fingerprint$output_leaves == job$transcripts * 3,
    fingerprint$output_carriers == job$transcripts * job$samples * 7 / 4,
    fingerprint$cds_bytes == job$transcripts * sum(nchar(input$expected$cds)),
    fingerprint$protein_bytes == job$transcripts * sum(nchar(input$expected$protein)),
    fingerprint$physical_edits == job$transcripts * 8, fingerprint$contributors == job$transcripts * 8)
  if (job$verify) {
    # Compare the complete expected carrier relation, not just row counts/sums.
    DBI::dbExecute(con, "CREATE TABLE expected_carriers AS SELECT transcript_index,sample_index,lane,
      list_sort(list(event_index)) events FROM calls CROSS JOIN range(1,3) l(lane)
      WHERE alleles[lane]=1 GROUP BY ALL")
    DBI::dbExecute(con, "CREATE VIEW actual_carriers AS SELECT transcript_index,c.sample_index,
      c.haplotype_lane lane,list_sort(list_transform(contributors,x->x.event_index)) events
      FROM measured,unnest(carriers) u(c)")
    exact <- function() DBI::dbGetQuery(con, "SELECT count(*)=0 ok FROM (
      (SELECT * FROM expected_carriers EXCEPT ALL SELECT * FROM actual_carriers)
      UNION ALL (SELECT * FROM actual_carriers EXCEPT ALL SELECT * FROM expected_carriers))")$ok
    stopifnot(exact(), DBI::dbGetQuery(con, sprintf("SELECT bool_and(
      m.cds=e.cds AND m.protein=e.protein AND m.sequence_status='ok' AND m.projection_status='ok'
      AND m.evidence_flags=1 AND list_unique(list_transform(m.contributors,x->x.event_index))=length(m.contributors)
      AND list_sum(list_transform(m.coding_blocks,x->length(x.event_indices)))=m.edit_count
      AND list_bool_and(list_transform(m.carriers,x->x.ploidy=2 AND x.phase_set=10))) ok
      FROM measured m JOIN expected_sequences e ON e.mask=list_sum(list_transform(m.contributors,
        x->(1::BIGINT << ((x.event_index-1)%%4))))"))$ok)
    # A duplicate output cannot be hidden by a set-only comparison.
    DBI::dbExecute(con, "INSERT INTO measured SELECT * FROM measured LIMIT 1")
    stopifnot(!exact())
  }
  c(list(seconds = elapsed, duckdb_version = as.character(utils::packageVersion("duckdb"))),
    as.list(fingerprint[1L, ]))
}

main <- function() {
  options(rlang_backtrace_on_error = "none")
  args <- commandArgs(TRUE)
  if (length(args) == 2L && args[1L] == "--job") {
    job <- readRDS(args[2L])
    input <- job$input
    result <- if (job$mode == "native") native(job, input) else sql(job, input)
    saveRDS(result, job$result)
    return(invisible(NULL))
  }
  options <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--transcripts", type = "integer", default = 1024L),
    optparse::make_option("--samples", type = "integer", default = 64L),
    optparse::make_option("--overlap", type = "integer", default = 16L),
    optparse::make_option("--passes", type = "integer", default = 3L),
    optparse::make_option("--cpu", type = "integer", default = 2L),
    optparse::make_option("--diagnostic", action = "store_true", default = FALSE),
    optparse::make_option("--extension-receipt", dest = "extension_receipt", default = NULL)
  )))
  stopifnot(options$transcripts > 0L, options$transcripts <= 100000L,
    options$samples >= 4L, options$samples <= 4096L, options$samples %% 4L == 0L,
    options$overlap > 0L, options$overlap <= 512L, options$transcripts %% options$overlap == 0L,
    options$passes >= 1L, options$passes <= 20L, options$cpu >= 0L)
  root <- normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
  source(file.path(root, "scripts/duckvep_evidence.R"), local = TRUE)
  revision <- duckvep_evidence_revision(root)
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(root, "r/duckhtsbench/inst/benchmark_registry.tsv"))
  paths <- duckhtsbench::duckhts_bench_stage_repository_fixtures(root, "duckvep-haplotypes")
  input <- fixture(paths)
  extension <- file.path(root, "build/release/duckhts.duckdb_extension")
  binding <- "diagnostic_unbound"
  if (!options$diagnostic) {
    duckvep_evidence_assert_checkout(root, revision)
    if (is.null(options$extension_receipt)) stop("supply a clean-build --extension-receipt")
    binding <- duckvep_evidence_read_extension_receipt(options$extension_receipt,
      root, extension, revision)$binding
  }
  out <- tempfile("haplotype_benchmark_", tmpdir = file.path(root, "test/duckvep/conformance/results"))
  dir.create(out)
  message("Benchmark artifacts: ", out)
  sources <- c("benchmarks/duckvep_haplotype_stream.c", paste0("src/duckvep/kernel/src/duckvep_",
    c("haplotype", "carriers", "phase", "haplotype_stream", "classify", "codon", "coding", "projection", "delta"), ".c"))
  shared <- file.path(out, paste0("stream", .Platform$dynlib.ext))
  compiler <- Sys.getenv("CC", "cc")
  flags <- c("-std=c99", "-O3", "-DNDEBUG", "-Wall", "-Wextra", "-fPIC", "-shared",
    "-I", "src/duckvep/kernel/src", "-I", "src/duckvep/kernel/include")
  stopifnot(system2(compiler, shQuote(c(flags, sources, "-o", shared)),
    stdout = file.path(out, "compiler.log"), stderr = file.path(out, "compiler.log")) == 0L)
  identities <- c(paths, extension, shared, sources,
    list.files("src/duckvep/kernel/src", "\\.(h|inc)$", full.names = TRUE),
    "src/duckvep/kernel/include/duckvep_kernel.h", "benchmarks/duckvep_haplotypes.R",
    "r/duckhtsbench/inst/benchmark_registry.tsv")
  hashes <- vapply(identities, duckvep_evidence_sha256, "")
  results <- list()
  for (mode in c("native", "sql")) {
    fingerprint <- NULL
    for (pass in 0:options$passes) {
      stem <- file.path(out, paste0(mode, "_", pass))
      job <- c(options[c("transcripts", "samples", "overlap")], list(mode = mode, verify = pass == 0L,
        input = input, extension = extension, shared = shared, result = paste0(stem, ".rds")))
      saveRDS(job, paste0(stem, "_job.rds"))
      command <- c("-v", "-o", paste0(stem, "_time.txt"), "taskset", "-c", as.character(options$cpu),
        file.path(R.home("bin"), "Rscript"), file.path(root, "benchmarks/duckvep_haplotypes.R"),
        "--job", paste0(stem, "_job.rds"))
      status <- system2("/usr/bin/time", shQuote(command), stdout = paste0(stem, ".log"),
        stderr = paste0(stem, ".log"))
      if (status != 0L) stop("worker failed: ", stem, ".log")
      value <- readRDS(job$result)
      observed <- value[setdiff(names(value), "seconds")]
      if (pass == 0L) fingerprint <- observed else stopifnot(identical(observed, fingerprint))
      if (pass == 0L) next
      time <- readLines(paste0(stem, "_time.txt"))
      rss <- as.numeric(sub(".*: ", "", time[grepl("Maximum resident set size", time, fixed = TRUE)]))
      stopifnot(length(rss) == 1L, is.finite(rss), rss > 0)
      results[[length(results) + 1L]] <- data.frame(source_revision = revision,
        extension_build_binding = binding, mode, pass, transcripts = options$transcripts,
        samples = options$samples, overlap = options$overlap, threads = 1L, cpu_affinity = options$cpu,
        peak_process_rss_kib = rss, as.data.frame(value), check.names = FALSE)
    }
  }
  columns <- unique(unlist(lapply(results, names)))
  results <- do.call(rbind, lapply(results, function(x) {
    x[setdiff(columns, names(x))] <- NA
    x[columns]
  }))
  shared_counts <- c("output_leaves", "output_carriers", "cds_bytes", "protein_bytes",
    "physical_edits", "contributors", "blocks")
  stopifnot(all(vapply(results[shared_counts], function(x) length(unique(x)) == 1L, TRUE)))
  stopifnot(identical(hashes, vapply(identities, duckvep_evidence_sha256, "")),
    identical(revision, duckvep_evidence_revision(root)))
  if (!options$diagnostic) duckvep_evidence_assert_checkout(root, revision)
  write.csv(results, file.path(out, "results.csv"), row.names = FALSE)
  jsonlite::write_json(list(source_revision = revision, extension_build_binding = binding,
    scope = "phased_literal_sequence_not_compound_so_hgvs", options = options,
    compiler = system2(compiler, "--version", stdout = TRUE), compiler_flags = flags,
    worker_command = "/usr/bin/time -v -o TIME taskset -c CPU Rscript benchmarks/duckvep_haplotypes.R --job JOB",
    cpu = system2("lscpu", stdout = TRUE), session = capture.output(sessionInfo()),
    sha256 = as.list(vapply(unique(c(identities, list.files(out, full.names = TRUE))),
      duckvep_evidence_sha256, ""))),
    file.path(out, "receipt.json"), pretty = TRUE, auto_unbox = TRUE)
  print(results[c("mode", "pass", "seconds", "peak_process_rss_kib", "output_leaves", "output_carriers")],
    row.names = FALSE)
}
main()
