# Scaling of duckhts_somalier_sex() over a sample-by-panel count relation.
#
# Data:    Rscript benchmarks/benchmark_somalier_sex_run.R --data <dir> <lib>
# Driver:  Rscript benchmarks/benchmark_somalier_sex_run.R <out.tsv> <dir> <lib>
# Worker:  Rscript benchmarks/benchmark_somalier_sex_run.R --worker <lib> <dir> <samples> <threads>
#
# <lib> is a library holding the installed Rduckhts build under test. The count
# relation is synthetic (see the report): 20,000 panel sites, 17,000 autosomal,
# 2,500 X and 500 Y, for 500, 1,000 or 2,000 samples, stored in a DuckDB file
# that each worker attaches read-only. Every worker is a fresh process.
# samples = 0 measures the fixed overhead (packages, connection and extension
# load, attach, no query). DuckDB peak buffer memory and temporary-directory
# size come from the query profile; peak RSS is VmHWM.
args <- commandArgs(TRUE)
memory_limit <- "128MB"
panel_sites <- 20000L
scales <- c(500L, 1000L, 2000L)
threads_grid <- c(1L, 4L)
repetitions <- 3L

if (identical(args[1], "--data")) {
  dir <- args[2]
  .libPaths(c(args[3], .libPaths()))
  suppressPackageStartupMessages({ library(DBI); library(Rduckhts) })
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (samples in scales) {
    path <- file.path(dir, sprintf("counts_%d.duckdb", samples))
    unlink(path)
    con <- rduckhts_connect(dbdir = path)
    dbExecute(con, sprintf(paste0(
      "CREATE TABLE panel AS SELECT 'GRCh38' AS assembly, i::UBIGINT AS site_index, ",
      "CASE WHEN i < 17000 THEN 'chr' || (1 + i %% 22)::VARCHAR WHEN i < 19500 THEN 'chrX' ",
      "ELSE 'chrY' END AS region, (1000 + i)::UBIGINT AS position, 'A' AS allele_a, ",
      "'G' AS allele_b FROM range(%d) t(i)"), panel_sites))
    dbExecute(con, sprintf(paste0(
      "CREATE TABLE counts AS WITH g AS (SELECT s, p.*, hash(s, p.site_index) %% 100 AS r, ",
      "8 + hash(s * 7, p.site_index) %% 50 AS depth FROM range(%d) sa(s), panel p) ",
      "SELECT 'S' || lpad(s::VARCHAR, 6, '0') AS sample_id, assembly, site_index, region, ",
      "position, allele_a, allele_b, ",
      "(CASE WHEN r < 50 THEN depth WHEN r < 80 THEN depth // 2 ELSE 0 END)::UBIGINT AS a, ",
      "(CASE WHEN r < 50 THEN 0 WHEN r < 80 THEN depth - depth // 2 ELSE depth END)::UBIGINT AS b, ",
      "(r %% 3)::UBIGINT AS other FROM g"), samples))
    dbDisconnect(con, shutdown = TRUE)
  }
  quit(save = "no")
}

if (identical(args[1], "--worker")) {
  .libPaths(c(args[2], .libPaths()))
  dir <- args[3]; samples <- as.integer(args[4]); threads <- as.integer(args[5])
  suppressPackageStartupMessages({ library(DBI); library(Rduckhts) })
  loaded <- find.package("Rduckhts")
  extension <- list.files(file.path(loaded, "duckhts_extension"), pattern = "\\.duckdb_extension$",
                          recursive = TRUE, full.names = TRUE)
  stopifnot(length(extension) == 1L)
  cat(sprintf("BUILD\t%s\n", digest::digest(file = extension, algo = "sha256")))
  con <- rduckhts_connect()
  spill <- tempfile("sex-spill-")
  dir.create(spill)
  dbExecute(con, sprintf("SET threads = %d", threads))
  dbExecute(con, sprintf("SET memory_limit = '%s'", memory_limit))
  dbExecute(con, sprintf("SET temp_directory = %s", dbQuoteString(con, spill)))
  rss_kib <- function() as.numeric(strsplit(trimws(grep("^VmHWM:", readLines("/proc/self/status"),
                                                        value = TRUE)), "[[:space:]]+")[[1]][2])
  if (samples == 0L) {
    cat(sprintf("RESULT\t0\t0\t0\t0\t0\t%.0f\n", rss_kib()))
    quit(save = "no")
  }
  dbExecute(con, sprintf("ATTACH %s AS c (READ_ONLY)",
                         dbQuoteString(con, file.path(dir, sprintf("counts_%d.duckdb", samples)))))
  input_rows <- dbGetQuery(con, "SELECT count(*) AS n FROM c.counts")$n
  profile <- tempfile("sex-profile-", fileext = ".json")
  dbExecute(con, "PRAGMA enable_profiling = 'json'")
  dbExecute(con, sprintf("PRAGMA profiling_output = %s", dbQuoteString(con, profile)))
  started <- proc.time()[["elapsed"]]
  result <- dbGetQuery(con, "SELECT * FROM duckhts_somalier_sex('c.counts', 'c.panel')")
  seconds <- proc.time()[["elapsed"]] - started
  dbExecute(con, "PRAGMA disable_profiling")
  m <- jsonlite::fromJSON(profile, simplifyVector = FALSE)
  # The result is consumed: every sample has a status and a call.
  stopifnot(nrow(result) == samples, !anyNA(result$inferred_sex), !anyNA(result$status))
  cat(sprintf("RESULT\t%.0f\t%d\t%.2f\t%.0f\t%.0f\t%.0f\n", input_rows, nrow(result), seconds,
              as.numeric(m$system_peak_buffer_memory), as.numeric(m$system_peak_temp_dir_size),
              rss_kib()))
  dbDisconnect(con, shutdown = TRUE)
  unlink(c(profile, spill), recursive = TRUE)
  quit(save = "no")
}

stopifnot(length(args) == 3L)
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
out_path <- args[1]; dir <- args[2]; lib <- args[3]
git_sha <- system2("git", c("-C", dirname(dirname(script)), "rev-parse", "HEAD"), stdout = TRUE)
fields <- c("input_rows", "output_rows", "seconds", "peak_buffer_bytes", "peak_temp_bytes",
            "peak_rss_kib")
columns <- c("mode", "samples", "threads", "repetition", "memory_limit", "git_sha",
             "extension_sha256", fields)
cat(paste(columns, collapse = "\t"), "\n", sep = "", file = out_path)
run <- function(mode, samples, threads, rep) {
  out <- suppressWarnings(system2("Rscript", shQuote(c(script, "--worker", lib, dir, samples,
                                                       threads)), stdout = TRUE, stderr = TRUE))
  build <- sub("^BUILD\t", "", grep("^BUILD\t", out, value = TRUE))
  result <- strsplit(sub("^RESULT\t", "", grep("^RESULT\t", out, value = TRUE)), "\t")[[1]]
  stopifnot(length(build) == 1L, length(result) == length(fields))
  cat(paste(c(mode, samples, threads, rep, memory_limit, git_sha, build, result), collapse = "\t"),
      "\n", sep = "", file = out_path, append = TRUE)
}
# Repetitions are interleaved across cells so drift does not align with one cell.
for (rep in seq_len(repetitions)) {
  for (threads in threads_grid) {
    run("overhead", 0L, threads, rep)
    for (samples in scales) run("run", samples, threads, rep)
  }
}
