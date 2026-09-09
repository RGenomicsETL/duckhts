#!/usr/bin/env Rscript
# Initialized native replay/count sink and full public SQL over registered inputs.

fixture <- function(paths) {
  fasta <- readLines(paths[["haplotype_benchmark_reference"]])
  stopifnot(length(fasta) == 2L, nchar(fasta[2L]) == 180L)
  events <- read.delim(paths[["haplotype_benchmark_events"]])
  stopifnot(identical(events$event_index, 1:4), all(diff(events$position) > 0))
  masks <- c(3L, 12L, 15L, 1L, 2L, 4L, 8L)
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
  protein <- as.character(Biostrings::translate(Biostrings::DNAStringSet(
    substring(cds, 1L, nchar(cds) %/% 3L * 3L))))
  protein <- sub("(\\*).*", "\\1", protein)
  list(reference = fasta[2L], events = events,
    expected = data.frame(mask = masks, cds, protein)[1:3, ],
    singletons = data.frame(mask = masks, cds, protein)[4:7, ],
    expected_reference = data.frame(mask = 0L, cds = fasta[2L], protein = sub("(\\*).*", "\\1",
      as.character(Biostrings::translate(Biostrings::DNAString(fasta[2L]))))))
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
  records <- identical(job$mode, "sql_records")
  singletons <- job$mode %in% c("sql_singletons", "sql_singletons_hgvs")
  hgvs <- identical(job$mode, "sql_singletons_hgvs")
  expected <- if (singletons) input$singletons else if (records)
    rbind(input$expected_reference, input$expected) else input$expected
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbExecute(con, paste("LOAD", q(job$extension)))
  DBI::dbExecute(con, "SET threads=1")
  DBI::dbExecute(con, "SET memory_limit='4GB'")
  DBI::dbWriteTable(con, "fixture_events", input$events)
  DBI::dbWriteTable(con, "expected_sequences", expected)
  DBI::dbExecute(con, sprintf("CREATE TABLE tx AS SELECT i::UINTEGER transcript_index,
    0::UINTEGER seq_region,(i//%d*1000+100)::UBIGINT transcript_start,
    (transcript_start+179)::UBIGINT transcript_end,1::TINYINT strand,i::UINTEGER gene_index,
    3::UBIGINT transcript_flags,transcript_start cds_start,transcript_end cds_end,
    %s::BLOB cds_sequence,1::UTINYINT codon_table%s FROM range(%d) t(i)",
    job$overlap, q(input$reference), if (singletons)
      ",''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence" else "", job$transcripts))
  exons <- "SELECT transcript_index,transcript_start exon_start,transcript_end exon_end,
    1::UBIGINT exon_cdna_start,180::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase FROM tx"
  regions <- if (singletons) sprintf("SELECT 0::UINTEGER seq_region,
    %d::UBIGINT sequence_length,'bench'::VARCHAR seq_region_name", job$reference_length) else
    "SELECT 0::UINTEGER seq_region"
  stopifnot(DBI::dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('bench',",
    q(regions), ",", q("SELECT * FROM tx"), ",", q(exons),
    if (singletons) paste0(",reference_fasta:=", q(job$reference_fasta)), ")"))$loaded)
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
  if (singletons) DBI::dbExecute(con, "UPDATE calls SET alleles=
    CASE WHEN sample_index%4=(event_index-1)%4 THEN [0,1] ELSE [0,0] END::INTEGER[]")
  denominators <- list()
  if (records) {
    DBI::dbExecute(con, "CREATE TABLE source_records AS SELECT
      event_index,seq_region,position,reference,[alternate] alternates,
      transcript_index,sample_index,alleles[1]::VARCHAR || '|' || alleles[2]::VARCHAR gt
      FROM calls")
  }
  if (records || singletons) {
    denominators <- as.list(DBI::dbGetQuery(con, paste0("SELECT
      count(DISTINCT event_index) input_physical_records,
      count(DISTINCT (event_index,sample_index)) input_record_sample_calls,
      count(*) input_candidate_sample_rows FROM ", if (records) "source_records" else "calls"))[1L, ])
    stopifnot(denominators$input_physical_records == job$transcripts / job$overlap * 4,
      denominators$input_record_sample_calls == job$transcripts / job$overlap * job$samples * 4,
      denominators$input_candidate_sample_rows == n)
  }
  # Decoded pools match the native lane. Raw preparation visits REF, ALT and
  # undefined-slot descriptors for each biallelic record; these are not records.
  limits <- c(max_active_events = 5, max_active_transcripts = job$overlap,
    max_active_carriers = 2 * job$samples * job$overlap, max_active_prefixes = 16 * job$overlap,
    max_active_projections = 4 * job$overlap + 1, max_allele_bytes = 64,
    max_leaf_events = 4, max_leaf_edits = 4, max_sequence_bases = 184, max_ploidy = 2,
    max_phase_sets = 1, max_alignment_cells = 65536, max_leaf_differences = 64,
    workspace_limit = 64 * 1024^2)
  if (records) limits[c("max_active_events", "max_active_projections", "max_allele_bytes")] <-
    c(13, 12 * job$overlap + 1, 192)
  if (singletons) limits <- c(limits, max_hgvs_operations = 4, max_hgvs_bytes = 128,
    max_hgvs_reference_bytes = 16384)
  query <- paste0("CREATE OR REPLACE TABLE measured AS SELECT * FROM duckvep_haplotypes(",
    q(if (records) "SELECT * FROM source_records" else "SELECT * FROM calls"), ",'bench',",
    if (records) "input_mode:='source_records',phase_policy:='vep116_compat'," else "",
    if (hgvs) "hgvs:=true," else "",
    paste(names(limits), ":=", limits, collapse = ","), ")")
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
  # Full output, local coding-block output, and literal replay are separate
  # contracts. Narrow projections never replace the full-output check.
  local_so <- DBI::dbGetQuery(con, "SELECT sum(octet_length(encode(to_json(h)))) local_so_json_bytes,
    bit_xor(hash(to_json(h)))::VARCHAR local_so_xor_hash,
    sum(hash(to_json(h))::HUGEINT)::VARCHAR local_so_sum_hash FROM (
      SELECT * EXCLUDE(hgvsp,hgvsp_status)
        REPLACE(list_sort(carriers) AS carriers,list_sort(contributors) AS contributors)
      FROM measured) h")
  replay <- DBI::dbGetQuery(con, "SELECT sum(octet_length(encode(to_json(h)))) replay_json_bytes,
    bit_xor(hash(to_json(h)))::VARCHAR replay_xor_hash,
    sum(hash(to_json(h))::HUGEINT)::VARCHAR replay_sum_hash FROM (
      SELECT * EXCLUDE(hgvsp,hgvsp_status)
        REPLACE(list_sort(carriers) AS carriers,list_sort(contributors) AS contributors,
        list_transform(coding_blocks,b->struct_pack(cds_start:=b.cds_start,reference:=b.reference,
          alternate:=b.alternate,alt_start0:=b.alt_start0,length_change:=b.length_change,
          sequence_flags:=b.sequence_flags,event_indices:=b.event_indices)) AS coding_blocks)
      FROM measured) h")
  stopifnot(fingerprint$output_leaves == job$transcripts * nrow(expected),
    fingerprint$output_carriers == job$transcripts * job$samples *
      (if (singletons) 1 else if (records) 2 else 7 / 4),
    fingerprint$cds_bytes == job$transcripts * sum(nchar(expected$cds)),
    fingerprint$protein_bytes == job$transcripts * sum(nchar(expected$protein)),
    fingerprint$physical_edits == job$transcripts * (if (singletons) 4 else 8),
    fingerprint$contributors == job$transcripts * (if (singletons) 4 else 8))
  if (job$verify) {
    # Compare the complete expected carrier relation, not just row counts/sums.
    DBI::dbExecute(con, paste0("CREATE TABLE expected_carriers AS SELECT transcript_index,sample_index,lane,
      ", if (records) "NULL::BIGINT" else if (singletons) "10::BIGINT" else
        "CASE WHEN sample_index%4=2 THEN NULL::BIGINT ELSE 10::BIGINT END", " phase_set,
      2::UINTEGER ploidy,
      coalesce(list_sort(list(event_index) FILTER(WHERE alleles[lane]=1)),[]::UBIGINT[]) events
      FROM calls CROSS JOIN range(1,3) l(lane) ",
      if (records) "" else "WHERE alleles[lane]=1 ", "GROUP BY ALL"))
    DBI::dbExecute(con, "CREATE VIEW actual_carriers AS SELECT transcript_index,c.sample_index,
      c.haplotype_lane lane,c.phase_set,c.ploidy,
      list_sort(list_transform(contributors,x->x.event_index)) events
      FROM measured,unnest(carriers) u(c)")
    exact <- function() DBI::dbGetQuery(con, "SELECT count(*)=0 ok FROM (
      (SELECT * FROM expected_carriers EXCEPT ALL SELECT * FROM actual_carriers)
      UNION ALL (SELECT * FROM actual_carriers EXCEPT ALL SELECT * FROM expected_carriers))")$ok
    sequences <- function() DBI::dbGetQuery(con, sprintf("SELECT bool_and(coalesce(
      m.cds=e.cds AND m.protein=e.protein AND m.sequence_status='ok' AND m.projection_status='ok'
      AND %s
      AND m.evidence_flags=CASE WHEN e.mask=0 THEN 0 ELSE 1 END
      AND list_unique(list_transform(m.contributors,x->x.event_index))=length(m.contributors)
      AND coalesce(list_sum(list_transform(m.coding_blocks,x->length(x.event_indices))),0)=m.edit_count
      AND %s,false)) ok FROM measured m JOIN expected_sequences e
      ON e.mask=coalesce(list_sum(list_transform(m.contributors,
        x->(1::BIGINT << ((x.event_index-1)%%4)))),0)",
      if (hgvs) "m.hgvsp IS NOT NULL AND m.hgvsp_status='ok'" else
        "m.hgvsp IS NULL AND m.hgvsp_status='not_requested'",
      if (records) "len(list_filter(m.contributors,x->x.alt_index IS DISTINCT FROM 1))=0"
        else "true"))$ok
    hgvs_equal <- function() TRUE
    if (hgvs) {
      DBI::dbExecute(con, "CREATE TABLE independent_events AS SELECT DISTINCT event_index,
        seq_region,position,reference,alternate,NULL::UBIGINT end_position,
        NULL::VARCHAR structural_type,NULL::VARCHAR copy_change,NULL::UINTEGER mate_seq_region,
        NULL::UBIGINT mate_position FROM calls")
      DBI::dbExecute(con, "CREATE TABLE independent AS SELECT * FROM duckvep_annotate(
        'independent_events','bench',hgvs:=true,upstream_distance:=0,downstream_distance:=0)")
      stopifnot(DBI::dbGetQuery(con, "SELECT count(*) n FROM independent")$n == job$transcripts * 4)
      # The independent renderer omits prediction parentheses. Compare every
      # event/transcript key and exact suffix; this is an internal path check,
      # not a substitute for the executable VEP differential.
      DBI::dbExecute(con, "CREATE VIEW actual_hgvs AS SELECT transcript_index,
        contributors[1].event_index event_index,
        CASE WHEN starts_with(hgvsp,'p.(') AND ends_with(hgvsp,')')
          THEN 'p.' || substr(hgvsp,4,length(hgvsp)-4) ELSE hgvsp END protein_hgvs
        FROM measured")
      hgvs_equal <- function() DBI::dbGetQuery(con, "SELECT count(*)=0 ok FROM (
        (SELECT * FROM actual_hgvs EXCEPT ALL
          SELECT transcript_index,event_index,protein_hgvs FROM independent)
        UNION ALL (SELECT transcript_index,event_index,protein_hgvs FROM independent
          EXCEPT ALL SELECT * FROM actual_hgvs))")$ok
      write.csv(DBI::dbGetQuery(con, "SELECT m.transcript_index,
        contributors[1].event_index event_index,hgvsp,hgvsp_status,
        i.protein_hgvs,i.protein_hgvs_status,i.protein_hgvs_reason
        FROM measured m FULL JOIN independent i ON m.transcript_index=i.transcript_index
          AND m.contributors[1].event_index=i.event_index ORDER BY event_index,m.transcript_index"),
        sub("\\.rds$", "_hgvs.csv", job$result), row.names = FALSE)
    }
    stopifnot(exact(), isTRUE(sequences()), isTRUE(hgvs_equal()))
    if (records || singletons) {
      # Every mutation is restored before the next independent control.
      DBI::dbExecute(con, "CREATE TABLE verified AS SELECT * FROM measured")
      mutations <- c(cds = "UPDATE measured SET cds='C' || substr(cds,2)",
        protein = "UPDATE measured SET protein='L' || substr(protein,2)",
        missing_cds = "UPDATE measured SET cds=NULL WHERE transcript_index=0",
        invented_hgvs = "UPDATE measured SET hgvsp='p.(Met1Leu)'",
        hgvs_status = paste0("UPDATE measured SET hgvsp_status='",
          if (hgvs) "not_requested" else "ok", "'"),
        carrier = "UPDATE measured SET carriers=list_transform(carriers,c->struct_update(c,sample_index:=4294967295::UINTEGER))",
        missing_ploidy = "UPDATE measured SET carriers=list_transform(carriers,c->struct_update(c,ploidy:=NULL::USMALLINT))",
        contributor = "UPDATE measured SET contributors=list_transform(contributors,c->struct_update(c,event_index:=0::UBIGINT))",
        dropped_row = "DELETE FROM measured WHERE transcript_index=0")
      if (records) mutations <- c(mutations,
        missing_alt_index = "UPDATE measured SET contributors=list_transform(contributors,c->struct_update(c,alt_index:=NULL::UINTEGER))")
      if (hgvs) mutations <- c(mutations,
        missing_hgvs = "UPDATE measured SET hgvsp=NULL",
        malformed_hgvs = "UPDATE measured SET hgvsp=replace(hgvsp,'Ala','(Ala)')")
      for (mutation in mutations) {
        DBI::dbExecute(con, mutation)
        stopifnot(!exact() || !isTRUE(sequences()) || !isTRUE(hgvs_equal()))
        DBI::dbExecute(con, "CREATE OR REPLACE TABLE measured AS SELECT * FROM verified")
      }
      if (records) {
        DBI::dbExecute(con, "UPDATE source_records SET gt='2|0'")
        invalid <- tryCatch({DBI::dbExecute(con, query); NULL}, error = identity)
        stopifnot(inherits(invalid, "error"), grepl("raw GT status", conditionMessage(invalid)))
        DBI::dbExecute(con, "CREATE OR REPLACE TABLE measured AS SELECT * FROM verified")
      }
    }
    # A duplicate output cannot be hidden by a set-only comparison.
    DBI::dbExecute(con, "INSERT INTO measured SELECT * FROM measured LIMIT 1")
    stopifnot(!exact())
  }
  c(list(seconds = elapsed, duckdb_version = as.character(utils::packageVersion("duckdb")),
    output_contract = if (singletons) paste0("singletons_local_coding_block_so_",
      if (hgvs) "hgvs" else "hgvs_status") else if (records)
      "source_records_local_coding_block_so_hgvs_status" else
      "local_coding_block_so_hgvs_status"),
    denominators,
    as.list(fingerprint[1L, ]), as.list(local_so[1L, ]), as.list(replay[1L, ]))
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
    optparse::make_option("--modes", default = "native,sql"),
    optparse::make_option("--diagnostic", action = "store_true", default = FALSE),
    optparse::make_option("--extension-receipt", dest = "extension_receipt", default = NULL)
  )))
  stopifnot(options$transcripts > 0L, options$transcripts <= 100000L,
    options$samples >= 4L, options$samples <= 4096L, options$samples %% 4L == 0L,
    options$overlap > 0L, options$overlap <= 512L, options$transcripts %% options$overlap == 0L,
    options$passes >= 1L, options$passes <= 20L, options$cpu >= 0L)
  modes <- strsplit(options$modes, ",", fixed = TRUE)[[1L]]
  stopifnot(length(modes) > 0L, !anyDuplicated(modes),
    all(modes %in% c("native", "sql", "sql_records", "sql_singletons", "sql_singletons_hgvs")))
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
  reference_fasta <- NULL
  reference_length <- NULL
  if (any(modes %in% c("sql_singletons", "sql_singletons_hgvs"))) {
    # Repeat the registered CDS at exactly the SQL model's genomic coordinates.
    # This per-run derivation and its faidx are retained with the worker jobs.
    sequence <- paste0(strrep("A", 99L), strrep(paste0(input$reference, strrep("A", 820L)),
      options$transcripts / options$overlap))
    reference_length <- nchar(sequence)
    reference_fasta <- file.path(out, "reference.fa")
    writeLines(c(">bench", sequence), reference_fasta)
    stopifnot(system2("samtools", c("faidx", shQuote(reference_fasta))) == 0L)
  }
  sources <- c("benchmarks/duckvep_haplotype_stream.c", paste0("src/duckvep/kernel/src/duckvep_",
    c("haplotype", "carriers", "phase", "haplotype_stream", "classify", "codon", "coding", "projection", "delta"), ".c"))
  shared <- file.path(out, paste0("stream", .Platform$dynlib.ext))
  compiler <- Sys.getenv("CC", "cc")
  flags <- c("-std=c99", "-O3", "-DNDEBUG", "-Wall", "-Wextra", "-fPIC", "-shared",
    "-I", "src/duckvep/kernel/src", "-I", "src/duckvep/kernel/include")
  stopifnot(system2(compiler, shQuote(c(flags, sources, "-o", shared)),
    stdout = file.path(out, "compiler.log"), stderr = file.path(out, "compiler.log")) == 0L)
  identities <- c(paths, if (!is.null(reference_fasta)) c(reference_fasta, paste0(reference_fasta, ".fai")),
    extension, shared, sources,
    list.files("src/duckvep/kernel/src", "\\.(h|inc)$", full.names = TRUE),
    "src/duckvep/kernel/include/duckvep_kernel.h", "benchmarks/duckvep_haplotypes.R",
    "r/duckhtsbench/inst/benchmark_registry.tsv")
  hashes <- vapply(identities, duckvep_evidence_sha256, "")
  results <- list()
  for (mode in modes) {
    fingerprint <- NULL
    for (pass in 0:options$passes) {
      stem <- file.path(out, paste0(mode, "_", pass))
      job <- c(options[c("transcripts", "samples", "overlap")], list(mode = mode, verify = pass == 0L,
        input = input, extension = extension, shared = shared, reference_fasta = reference_fasta,
        reference_length = reference_length, result = paste0(stem, ".rds")))
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
  decoded <- results[results$mode %in% c("native", "sql"), shared_counts, drop = FALSE]
  stopifnot(all(vapply(decoded, function(x) length(unique(x)) <= 1L, TRUE)))
  if (any(results$mode %in% c("sql_singletons", "sql_singletons_hgvs"))) {
    singleton_rows <- results[results$mode %in% c("sql_singletons", "sql_singletons_hgvs"),
      c(shared_counts, "local_so_json_bytes", "local_so_xor_hash", "local_so_sum_hash"), drop = FALSE]
    stopifnot(all(vapply(singleton_rows, function(x) length(unique(x)) <= 1L, TRUE)))
  }
  stopifnot(identical(hashes, vapply(identities, duckvep_evidence_sha256, "")),
    identical(revision, duckvep_evidence_revision(root)))
  if (!options$diagnostic) duckvep_evidence_assert_checkout(root, revision)
  write.csv(results, file.path(out, "results.csv"), row.names = FALSE)
  jsonlite::write_json(list(source_revision = revision, extension_build_binding = binding,
    scope = "initialized_native_decoded_replay_and_selected_sql_input_modes_not_whole_haplotype_so_hgvs",
    options = options,
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
