# Run from the repository root. Each timed query uses a fresh DuckDB process.
args <- commandArgs(trailingOnly = TRUE)
source("r/duckhtsbench/R/registry.R")
source("r/duckhtsbench/R/stage.R")
Sys.setenv(DUCKHTSBENCH_REGISTRY = normalizePath(
  "r/duckhtsbench/inst/benchmark_registry.tsv", mustWork = TRUE))

if (length(args) == 6L && args[[1L]] == "--child") {
  library(DBI)
  case <- args[[2L]]
  threads <- as.integer(args[[3L]])
  multiplier <- as.integer(args[[4L]])
  extension <- args[[5L]]
  input <- args[[6L]]
  con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true"),
                                  shared_home = FALSE))
  quote_string <- function(value) as.character(dbQuoteString(con, value))
  dbExecute(con, paste("LOAD", quote_string(extension)))
  dbExecute(con, sprintf("SET threads=%d", threads))
  if (case == "vcf_counts") {
    dbExecute(con, paste0(
      "CREATE TEMP TABLE panel AS WITH sites AS (",
      "SELECT CHROM AS region, POS AS position, REF AS ref, ALT[1] AS alt ",
      "FROM read_bcf(", quote_string(input), ", samples := '', ",
      "scan_mode := 'sequential') WHERE len(REF) = 1 AND len(ALT) = 1 ",
      "AND len(ALT[1]) = 1 AND REF IN ('A','C','G','T') ",
      "AND ALT[1] IN ('A','C','G','T') AND REF != ALT[1] ",
      "ORDER BY POS LIMIT 2000) ",
      "SELECT 'GRCh38' AS assembly, ",
      "CAST(row_number() OVER (ORDER BY position) - 1 AS UBIGINT) AS site_index, ",
      "region, CAST(position AS UBIGINT) AS position, ",
      "least(ref,alt) AS allele_a, greatest(ref,alt) AS allele_b FROM sites"))
    query <- paste0("CREATE TEMP TABLE result AS SELECT * FROM ",
                    "duckhts_somalier_vcf_counts(", quote_string(input), ", 'panel')")
  } else if (case == "import_sites") {
    query <- paste0("CREATE TEMP TABLE result AS SELECT * FROM ",
                    "duckhts_somalier_import_sites(", quote_string(input),
                    ", 'GRCh37', max_sites := 2000000)")
  } else if (case %in% c("munge", "munge_metal", "r_munge", "r_munge_metal")) {
    limit <- 200000L * multiplier
    dbExecute(con, paste0(
      "CREATE TEMP VIEW source AS SELECT CAST(CHR AS VARCHAR) AS CHR, ",
      "* EXCLUDE CHR FROM read_csv(", quote_string(input),
      ", delim := '\t', header := true) LIMIT ", limit))
    column_map <- paste0(
      "map(['CHR','BP','SNP','A1','A2','P','Z','FRQ','NEFF','HET_I2',",
      "'HET_P','DIRE'], ['CHR','BP','MarkerName','Allele1','Allele2',",
      "'P-value','Zscore','Freq1','Weight','HetISq','HetPVal','Direction'])")
    macro <- if (case %in% c("munge", "r_munge")) "duckdb_munge" else "duckdb_munge_metal"
    fasta <- duckhts_bench_artifact_path("liftover_grch37_fasta")
    stopifnot(file.exists(fasta))
    query <- paste0("CREATE TEMP TABLE result AS SELECT * FROM ", macro,
                    "('source', column_map := ", column_map,
                    ", fasta_ref := ", quote_string(fasta), ")")
  } else if (case == "liftover") {
    limit <- 250000L * multiplier
    chain <- duckhts_bench_artifact_path("liftover_grch37_grch38_chain")
    dst <- duckhts_bench_artifact_path("liftover_grch38_fasta")
    src <- duckhts_bench_artifact_path("liftover_grch37_fasta")
    stopifnot(file.exists(chain), file.exists(dst), file.exists(src))
    dbExecute(con, paste0(
      "CREATE TEMP VIEW source AS SELECT CHROM AS chrom, POS AS pos, ",
      "REF AS ref, ALT[1] AS alt FROM read_bcf(", quote_string(input),
      ", samples := '', scan_mode := 'sequential') LIMIT ", limit))
    query <- paste0(
      "CREATE TEMP TABLE result AS SELECT * FROM duckdb_liftover(",
      "'source', 'chrom', 'pos', ref_col := 'ref', alt_col := 'alt', ",
      "chain_path := ", quote_string(chain), ", dst_fasta_ref := ",
      quote_string(dst), ", src_fasta_ref := ", quote_string(src), ")")
  } else {
    stop("unknown case: ", case)
  }
  profile_file <- tempfile(fileext = ".json")
  dbExecute(con, "PRAGMA enable_profiling='json'")
  dbExecute(con, paste("PRAGMA profiling_output=", quote_string(profile_file)))
  started <- proc.time()[["elapsed"]]
  if (startsWith(case, "r_")) {
    mapping <- c(CHR = "CHR", BP = "BP", SNP = "MarkerName",
                 A1 = "Allele1", A2 = "Allele2", P = "P-value",
                 Z = "Zscore", FRQ = "Freq1", NEFF = "Weight")
    if (case == "r_munge_metal") {
      mapping <- c(mapping, HET_I2 = "HetISq", HET_P = "HetPVal",
                   DIRE = "Direction")
    }
    result <- Rduckhts::rduckhts_munge(con, "source", column_map = mapping,
                                      fasta_ref = fasta)
    seconds <- proc.time()[["elapsed"]] - started
    rows <- nrow(result)
    rm(result)
  } else {
    dbExecute(con, query)
    seconds <- proc.time()[["elapsed"]] - started
    rows <- dbGetQuery(con, "SELECT count(*) AS n FROM result")$n[[1L]]
  }
  dbExecute(con, "PRAGMA disable_profiling")
  profile <- jsonlite::fromJSON(profile_file)
  saved_result <- Sys.getenv("DUCKHTS_SCALING_RESULT", "")
  if (nzchar(saved_result)) {
    dbExecute(con, paste0("COPY (SELECT * FROM result ORDER BY sample_id, site_index) TO ",
                          quote_string(saved_result), " (FORMAT PARQUET)"))
  }
  rss <- readLines("/proc/self/status", warn = FALSE)
  peak_rss_kib <- as.numeric(gsub("[^0-9]", "", rss[startsWith(rss, "VmHWM:")]))
  input_records <- if (case %in% c("munge", "munge_metal", "r_munge",
                                  "r_munge_metal", "liftover")) {
    limit
  } else {
    as.integer(system2("bcftools", c("index", "-n", shQuote(input)),
                       stdout = TRUE))
  }
  write.table(data.frame(case, multiplier, threads, input_records, seconds,
                         peak_rss_kib,
                         peak_buffer_mib = profile$system_peak_buffer_memory / 1048576,
                         result_rows = rows), stdout(), sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  dbDisconnect(con, shutdown = TRUE)
  saved_profile <- Sys.getenv("DUCKHTS_SCALING_PROFILE", "")
  if (nzchar(saved_profile)) file.copy(profile_file, saved_profile, overwrite = TRUE)
  unlink(profile_file)
} else {
  stopifnot(length(args) %in% c(1L, 2L), file.exists(args[[1L]]))
  wrapper_run <- length(args) == 2L && identical(args[[2L]], "--wrappers")
  stopifnot(length(args) == 1L || wrapper_run)
  extension <- normalizePath(args[[1L]], mustWork = TRUE)
  output <- if (wrapper_run) "benchmarks/sql_scaling_audit_wrappers.tsv" else {
    "benchmarks/sql_scaling_audit_real.tsv"
  }
  cache <- file.path(duckhts_bench_cache_dir(), "benchmarks/sql-scaling-audit")
  inputs <- if (wrapper_run) {
    c(r_munge = "epilepsy", r_munge_metal = "epilepsy")
  } else {
    c(vcf_counts = "giab", import_sites = "phase3", munge = "epilepsy",
      munge_metal = "epilepsy", liftover = "phase3_source")
  }
  header <- paste(c("case", "multiplier", "threads", "input_records", "seconds",
                    "peak_rss_kib", "peak_buffer_mib", "result_rows"), collapse = "\t")
  writeLines(header, output)
  for (case in names(inputs)) {
    for (multiplier in c(1L, 2L, 4L)) {
      if (case %in% c("munge", "munge_metal", "r_munge", "r_munge_metal")) {
        input <- duckhts_bench_artifact_path("sql_scaling_epilepsy")
      } else if (case == "liftover") {
        input <- duckhts_bench_artifact_path("sql_scaling_1000g_phase3_chr22")
      } else {
        input <- file.path(cache, sprintf("%s-%dx.vcf.gz", inputs[[case]], multiplier))
        stopifnot(file.exists(paste0(input, ".tbi")))
      }
      stopifnot(file.exists(input))
      for (threads in c(1L, 4L)) {
        command <- c("benchmarks/benchmark_sql_scaling_real_run.R", "--child",
                     case, threads, multiplier, extension, input)
        stderr_file <- tempfile("sql-scaling-stderr-")
        result <- suppressWarnings(system2("Rscript", shQuote(command),
                                           stdout = TRUE, stderr = stderr_file))
        if (length(result) != 1L || !is.null(attr(result, "status"))) {
          stop(paste("Failed:", case, multiplier, threads,
                     paste(readLines(stderr_file, warn = FALSE), collapse = "\n")))
        }
        write(result, file = output, append = TRUE)
        message(result)
        unlink(stderr_file)
      }
    }
  }
}
