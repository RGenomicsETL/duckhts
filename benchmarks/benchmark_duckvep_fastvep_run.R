#!/usr/bin/env Rscript

# Paired fresh-process observations. All engines include startup, model loading
# and real-file output; CSQ includes flattening the upstream VCF transport.
duckvep_fastvep_run_job <- function(job, root) {
  run <- function(executable, args, env = character()) {
    status <- system2(executable, shQuote(args), env = env)
    if (status != 0L) stop("benchmark command failed: ", executable, " (", status, ")")
  }
  if (job$engine == "duckvep") {
    args <- c(file.path(root, "benchmarks/benchmark_duckvep_fastvep_worker.R"),
      "--extension", job$extension, "--model", job$model, "--input", job$input,
      "--output", job$output, "--output-contract", job$contract,
      "--threads", job$threads, "--distance", job$distance, "--memory-limit", job$memory_limit,
      "--max-spill", job$max_spill)
    if (job$contract != "operational17") args <- c(args, "--gff3", job$gff3)
    if (job$contract == "vep_csq") args <- c(args, "--fasta", job$fasta, "--include-identity")
    run("Rscript", args)
  } else {
    output <- if (job$contract == "vep_csq") paste0(job$output, ".vcf") else job$output
    args <- c("annotate", "--input", job$input, "--output", output,
      "--output-format", if (job$contract == "vep_csq") "vcf" else "tab",
      "--transcript-cache", job$cache, "--distance", job$distance, "--no-progress")
    if (job$contract == "vep_csq") args <- c(args, "--hgvs", "--fasta", job$fasta)
    run(job$fastvep, args, paste0("RAYON_NUM_THREADS=", job$threads))
    if (job$contract == "vep_csq") {
      run("Rscript", c(file.path(root, "benchmarks/benchmark_duckvep_fastvep_extract.R"),
        "--input", output, "--output", job$output, "--threads", job$threads,
        "--memory-limit", job$memory_limit, "--max-spill", job$max_spill,
        "--source-map", job$source_map))
      unlink(output)
    }
  }
}

main <- function() {
  root <- normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
  parser <- optparse::OptionParser(option_list = list(
    optparse::make_option("--job", default = ""),
    optparse::make_option("--output", default = ""),
    optparse::make_option("--work-dir", dest = "work_dir", default = dirname(tempdir())),
    optparse::make_option("--input-id", dest = "input_id", default = "variantkey_giab_hg002_v421"),
    optparse::make_option("--extension", default = "build/release/duckhts.duckdb_extension"),
    optparse::make_option("--extension-receipt", dest = "extension_receipt", default = ""),
    optparse::make_option("--fastvep", default = ".sync/fastVEP/target/release/fastvep"),
    optparse::make_option("--fastvep-build-receipt", dest = "fastvep_build_receipt", default = NULL),
    optparse::make_option("--checkout", default = ".sync/fastVEP"),
    optparse::make_option("--affinity-one", dest = "affinity_one", default = "2"),
    optparse::make_option("--affinity-four", dest = "affinity_four", default = "2,4,6,8"),
    optparse::make_option("--memory-limit", dest = "memory_limit", default = "4GB"),
    optparse::make_option("--max-spill", dest = "max_spill", default = "8GB"),
    optparse::make_option("--repetitions", type = "integer", default = 3L),
    optparse::make_option("--diagnostic", action = "store_true", default = FALSE)))
  opt <- optparse::parse_args(parser)
  if (nzchar(opt$job)) return(duckvep_fastvep_run_job(readRDS(opt$job), root))
  if (!nzchar(opt$output) || file.exists(opt$output) || opt$repetitions < 1L) {
    stop("a new --output directory and positive repetitions are required")
  }
  if (!opt$diagnostic && is.null(opt$fastvep_build_receipt)) {
    stop("published observations require --fastvep-build-receipt from a fresh build of the pinned source")
  }
  for (count in c(1L, 4L)) {
    affinity <- if (count == 1L) opt$affinity_one else opt$affinity_four
    cpus <- strsplit(affinity, ",", fixed = TRUE)[[1L]]
    if (length(cpus) != count || anyDuplicated(cpus) || any(!grepl("^[0-9]+$", cpus))) {
      stop("affinity must name exactly ", count, " distinct CPU IDs")
    }
  }
  source(file.path(root, "scripts/duckvep_evidence.R"))
  source(file.path(root, "benchmarks/benchmark_duckvep_fastvep_fields.R"))
  for (name in c("registry", "stage", "duckvep", "fastvep", "fastvep_source")) {
    source(file.path(root, "r/duckhtsbench/R", paste0(name, ".R")))
  }
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(root, "r/duckhtsbench/inst/benchmark_registry.tsv"))
  registry <- duckhts_bench_registry()
  cache_id <- "fastvep_ensembl116_duckvep_matched_cache"
  fastvep_identity <- duckhts_bench_identity_fields(
    registry$supplier_identity[registry$id == cache_id])
  fastvep_build <- if (!is.null(opt$fastvep_build_receipt)) {
    duckhts_bench_read_fastvep_build(opt$fastvep_build_receipt,
      fastvep_identity[["source_commit"]], opt$fastvep)
  } else NULL
  revision <- duckvep_evidence_revision(root)
  extension <- normalizePath(opt$extension, mustWork = TRUE)
  if (!opt$diagnostic) {
    if (!nzchar(opt$extension_receipt)) stop("published observations require --extension-receipt")
    duckvep_evidence_assert_checkout(root, revision)
    duckvep_evidence_read_extension_receipt(opt$extension_receipt, root, extension, revision)
    if (opt$repetitions < 3L) stop("published paired observations require at least three repetitions")
  }
  staged <- duckhts_bench_stage_fastvep(root, opt$checkout, opt$fastvep,
    cache_id = cache_id, extension = extension)
  inputs <- c(input = opt$input_id, model = "duckvep_ensembl116_model",
    fasta = "ensembl116_grch38_fasta_fa", gff3 = "fastvep_ensembl116_duckvep_gff3")
  paths <- vapply(inputs, duckhts_bench_artifact_path, character(1L))
  if (!all(file.exists(paths))) stop("all registered benchmark inputs must be staged first")
  for (id in inputs) duckhts_bench_validate_identity(id, duckhts_bench_artifact_path(id))
  identity <- duckhts_bench_identity_fields(
    registry$supplier_identity[registry$id == inputs[["model"]]])
  model_receipt <- duckhts_bench_validate_duckvep_ensembl116_model(
    paths[["model"]], extension, identity[["source_manifest_sha256"]]
  )
  matched_receipt <- duckhts_bench_fastvep_model_gff_receipt(paths[["gff3"]])
  if (matched_receipt[["model_sha256"]] != as.character(model_receipt$model_sha256[[1L]]) ||
      matched_receipt[["model_sha256"]] != fastvep_identity[["model_sha256"]] ||
      matched_receipt[["transcript_count"]] != fastvep_identity[["transcripts"]]) {
    stop("FastVEP GFF3 is not bound to the exact DuckVEP model inventory")
  }
  if (!dir.create(opt$output, recursive = TRUE)) stop("could not create a new output directory")
  work <- tempfile("duckvep-fastvep-", tmpdir = normalizePath(opt$work_dir, mustWork = TRUE))
  dir.create(work)
  message("Working files: ", work)
  source_bundle <- if (opt$input_id == "variantkey_giab_hg002_v421") {
    duckhts_bench_stage_fastvep_source_map(root, extension, opt$memory_limit, opt$max_spill)
  } else NULL
  source_map <- if (is.null(source_bundle)) file.path(opt$output, "source_alleles.parquet")
    else source_bundle[["map"]]
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
  DBI::dbExecute(con, paste("LOAD", DBI::dbQuoteString(con, extension)))
  DBI::dbExecute(con, paste("SET memory_limit =", DBI::dbQuoteString(con, opt$memory_limit)))
  DBI::dbExecute(con, paste("SET max_temp_directory_size =", DBI::dbQuoteString(con, opt$max_spill)))
  DBI::dbExecute(con, paste("SET temp_directory =", DBI::dbQuoteString(con, file.path(work, "source_spill"))))
  if (is.null(source_bundle)) duckvep_fastvep_write_source_map(con, paths[["input"]], source_map)
  mapped <- DBI::dbGetQuery(con, paste0("SELECT count(DISTINCT record_index)::VARCHAR records,
    count(*)::VARCHAR alt_alleles, count(*) FILTER (WHERE eligible)::VARCHAR
      eligible_literal_alleles FROM read_parquet(", DBI::dbQuoteString(con, source_map), ")"))
  counts <- DBI::dbGetQuery(con, paste0("SELECT count(*)::VARCHAR records,
    count(*) FILTER (WHERE len(ALT) > 0)::VARCHAR records_with_alt,
    sum(len(ALT))::VARCHAR alt_alleles,
    sum(len(list_filter(ALT, lambda a: regexp_full_match(REF, '[ACGTNacgtn]+')
      AND regexp_full_match(a, '[ACGTNacgtn]+') AND upper(REF) <> upper(a))))::VARCHAR
      eligible_literal_alleles FROM read_bcf(", DBI::dbQuoteString(con, paths[["input"]]),
    ", scan_mode := 'sequential', decompression_threads := 0)"))
  DBI::dbDisconnect(con, shutdown = TRUE)
  stopifnot(identical(mapped$records, counts$records_with_alt),
    identical(mapped$alt_alleles, counts$alt_alleles),
    identical(mapped$eligible_literal_alleles, counts$eligible_literal_alleles))
  # Keep a failed job and its log for diagnosis; remove successful output files
  # only after a full-row receipt has been written.
  source_files <- file.path(root, "benchmarks", paste0("benchmark_duckvep_fastvep_",
    c("run", "worker", "fields", "extract", "receipt"), ".R"))
  bound_files <- c(paths, fasta_index = paste0(paths[["fasta"]], ".fai"),
    extension = extension, fastvep = normalizePath(opt$fastvep),
    cache = staged[["cache"]], cache_receipt = staged[["receipt"]], source_map = source_map,
    stats::setNames(source_files, basename(source_files)))
  matched_gff_receipt <- file.path(opt$output, "matched_gff_receipt.tsv")
  if (!file.copy(staged[["gff3_receipt"]], matched_gff_receipt)) {
    stop("could not retain model-matched GFF3 provenance")
  }
  bound_files <- c(bound_files, matched_gff_receipt = matched_gff_receipt)
  if (!is.null(source_bundle)) {
    if (!file.copy(source_bundle[["receipt"]], file.path(opt$output, "source_map_receipt.tsv"))) {
      stop("could not retain registered source-map provenance")
    }
    bound_files <- c(bound_files, source_map_receipt = source_bundle[["receipt"]])
  }
  if (!is.null(fastvep_build)) {
    build_files <- c(fastvep_build = opt$fastvep_build_receipt,
      fastvep_build_log = file.path(dirname(opt$fastvep_build_receipt), fastvep_build[["log"]]),
      fastvep_source_tree = file.path(dirname(opt$fastvep_build_receipt), fastvep_build[["source_tree"]]),
      fastvep_source_commit_object = file.path(dirname(opt$fastvep_build_receipt),
        fastvep_build[["source_commit_object"]]))
    retained <- file.path(opt$output, c("fastvep_build.tsv", fastvep_build[["log"]],
      fastvep_build[["source_tree"]], fastvep_build[["source_commit_object"]]))
    if (!all(file.copy(build_files, retained))) stop("could not retain FastVEP build provenance")
    bound_files <- c(bound_files, build_files)
  }
  hashes <- vapply(bound_files, duckvep_evidence_sha256, character(1L))
  if (!is.null(fastvep_build) && !identical(hashes[["fastvep"]], fastvep_build[["executable_sha256"]])) {
    stop("FastVEP executable changed after build-receipt validation")
  }
  utils::write.csv(data.frame(artifact = names(bound_files), path = unname(bound_files), sha256 = hashes),
    file.path(opt$output, "inputs.csv"), row.names = FALSE)
  if (!file.copy(staged[["receipt"]], file.path(opt$output, "fastvep_cache_receipt.tsv"))) {
    stop("could not retain cache provenance")
  }
  if (is.na(counts$records) || as.numeric(counts$records) == 0 ||
      is.na(counts$eligible_literal_alleles) || as.numeric(counts$eligible_literal_alleles) == 0) {
    stop("benchmark input must contain records and eligible literal ALT alleles")
  }
  metadata <- c(source_revision = revision,
    binding = if (opt$diagnostic) "diagnostic_unbound" else "source_bound",
    fastvep_binding = if (is.null(fastvep_build)) "diagnostic_binary_unbound" else fastvep_build[["binding"]],
    fastvep_source_revision = duckvep_evidence_command("git", c("-C", opt$checkout,
      "rev-parse", "HEAD"), "cannot identify FastVEP source"),
    fastvep_version = duckvep_evidence_command(opt$fastvep, "--version", "cannot identify FastVEP executable"),
    fastvep_cache_id = cache_id,
    matched_model_sha256 = duckhts_bench_fastvep_model_gff_receipt(paths[["gff3"]])[["model_sha256"]],
    matched_transcripts = fastvep_identity[["transcripts"]],
    input_id = opt$input_id, input_records = counts$records, input_alt_alleles = counts$alt_alleles,
    eligible_literal_alleles = counts$eligible_literal_alleles,
    r_version = R.version.string, duckdb_version = as.character(utils::packageVersion("duckdb")),
    output_filesystem = duckvep_evidence_command("stat", c("-f", "-c", "%T", work),
      "cannot identify output filesystem"),
    distance = "5000", memory_limit = opt$memory_limit, max_spill = opt$max_spill,
    supplementary_providers = "none",
    affinity_one = opt$affinity_one, affinity_four = opt$affinity_four)
  utils::write.csv(data.frame(field = names(metadata), value = unname(metadata)),
    file.path(opt$output, "metadata.csv"), row.names = FALSE)
  configurations <- data.frame(engine = c(rep("duckvep", 3), rep("fastvep", 2)),
    contract = c("operational17", "native_tab17", "vep_csq", "native_tab17", "vep_csq"))
  completed <- character()
  coverage_files <- character()
  for (threads in c(1L, 4L)) for (run in seq_len(opt$repetitions)) {
    # Alternate tool order across repetitions to expose drift in a shared host.
    order <- if (run %% 2L) seq_len(nrow(configurations)) else rev(seq_len(nrow(configurations)))
    for (i in order) {
      job <- c(as.list(configurations[i, ]), as.list(paths), list(extension = extension,
        fastvep = normalizePath(opt$fastvep), cache = staged[["cache"]], source_map = source_map,
        threads = threads,
        distance = 5000L, memory_limit = opt$memory_limit, max_spill = opt$max_spill,
        output = file.path(work, "output.tsv")))
      label <- paste(job$engine, job$contract, threads, run, sep = "_")
      message("Running ", label)
      job_path <- file.path(work, "job.rds")
      saveRDS(job, job_path)
      timing <- file.path(opt$output, paste0(label, ".time"))
      log <- file.path(opt$output, paste0(label, ".log"))
      affinity <- if (threads == 1L) opt$affinity_one else opt$affinity_four
      args <- c("-v", "-o", timing, "taskset", "-c", affinity, "Rscript",
        file.path(root, "benchmarks/benchmark_duckvep_fastvep_run.R"), "--job", job_path)
      status <- system2("/usr/bin/time", shQuote(args), stdout = log, stderr = log)
      if (status != 0L) stop("benchmark failed; files retained at ", work, "; log: ", log)
      skip <- duckvep_fastvep_tab_header(job$output, duckvep_fastvep_transport_fields(job$contract))
      receipt <- file.path(opt$output, paste0(label, ".csv"))
      coverage <- file.path(opt$output, paste0(label, ".coverage.csv"))
      args <- c(file.path(root, "benchmarks/benchmark_duckvep_fastvep_receipt.R"),
        "--input", job$output, "--tool", job$engine, "--output-contract", job$contract,
        "--threads", threads, "--run", run, "--timing-file", timing, "--output", receipt)
      args <- c(args, "--skip-lines", skip, "--source-map", source_map,
        "--source-sha256", hashes[["input"]], "--coverage-output", coverage,
        "--memory-limit", opt$memory_limit, "--max-spill", opt$max_spill)
      if (system2("Rscript", shQuote(args)) != 0L) stop("receipt failed; output retained at ", work)
      result <- utils::read.csv(receipt, colClasses = "character")
      if (nrow(result) != 1L || is.na(result$row_count) || as.numeric(result$row_count) == 0) {
        stop("empty or invalid output receipt; output retained at ", work)
      }
      completed <- c(completed, receipt)
      coverage_files <- c(coverage_files, coverage)
      unlink(job$output)
    }
  }
  if (!identical(hashes, vapply(bound_files, duckvep_evidence_sha256, character(1L)))) {
    stop("benchmark source, executable or input changed during measurement")
  }
  if (!opt$diagnostic) duckvep_evidence_assert_checkout(root, revision)
  duckhts_bench_fastvep_source(opt$checkout, metadata[["fastvep_source_revision"]])
  stopifnot(length(completed) == nrow(configurations) * 2L * opt$repetitions)
  observations <- do.call(rbind, lapply(completed, utils::read.csv, colClasses = "character"))
  coverage <- do.call(rbind, lapply(coverage_files, utils::read.csv, colClasses = "character"))
  stopifnot(nrow(coverage) == nrow(observations), all(coverage$passed == "TRUE"),
    all(coverage$source_map_sha256 == hashes[["source_map"]]),
    all(coverage$source_alleles == metadata[["eligible_literal_alleles"]]))
  utils::write.csv(coverage, file.path(opt$output, "allele_coverage.csv"), row.names = FALSE)
  groups <- split(observations, paste(observations$tool, observations$output_contract))
  fingerprints <- c("row_count", "xor_hash", "low32_sum", "high32_sum")
  for (name in names(groups)) {
    if (nrow(unique(groups[[name]][fingerprints])) != 1L) {
      stop("output multiset changes across repeats or thread counts: ", name)
    }
  }
  artifacts <- list.files(opt$output, full.names = TRUE)
  utils::write.csv(data.frame(file = basename(artifacts),
    sha256 = vapply(artifacts, duckvep_evidence_sha256, character(1L))),
    file.path(opt$output, "artifacts.csv"), row.names = FALSE)
  utils::write.csv(data.frame(source_revision = revision, completed = TRUE,
    binding = metadata[["binding"]], observations = length(completed), repetitions = opt$repetitions,
    allele_coverage = "verified",
    artifacts_sha256 = duckvep_evidence_sha256(file.path(opt$output, "artifacts.csv"))),
    file.path(opt$output, "completion.csv"), row.names = FALSE)
  unlink(work, recursive = TRUE)
  message("Paired observations and final source-allele coverage retained: ", opt$output)
}

if (sys.nframe() == 0L) main()
