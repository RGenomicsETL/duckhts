# Panel-building scaling for rduckhts_ancestry_panel() on bigsnpr's staged reference.
#
# Driver:  Rscript benchmarks/benchmark_ancestry_panel_run.R <out.tsv> <name>=<lib> ...
# Worker:  Rscript benchmarks/benchmark_ancestry_panel_run.R --worker <lib> <quarters> <mode> <threads>
#
# Each worker is a fresh process. The reference is subsampled to 1, 2 or 4 quarters
# of its loci by hash, expanded to the long frequency and loading relations the
# builder reads, and consumed in full into a committed 17,000-site panel. Modes:
# "spaced" uses the builder's window spacing; "candidates" offers every locus as a
# candidate. quarters = 0 measures the fixed overhead (packages, connection and
# extension load, no panel). DuckDB peak buffer memory and temporary-directory size
# are the maxima over every statement the builder runs; peak RSS is VmHWM.
args <- commandArgs(TRUE)
memory_limit <- "16GB"
max_sites <- 17000L

if (identical(args[1], "--worker")) {
  .libPaths(c(args[2], .libPaths()))
  quarters <- as.integer(args[3]); mode <- args[4]; threads <- as.integer(args[5])
  suppressPackageStartupMessages({
    library(DBI); library(Rduckhts); library(duckhtsbench)
  })
  # Identify the Rduckhts build this process loaded, before any work that could fail.
  loaded <- find.package("Rduckhts")
  extension <- list.files(file.path(loaded, "duckhts_extension"), pattern = "\\.duckdb_extension$",
                          recursive = TRUE, full.names = TRUE)
  stopifnot(length(extension) == 1L)
  cat(sprintf("BUILD\t%s\t%s\n",
              digest::digest(file = file.path(loaded, "R", "Rduckhts.rdb"), algo = "sha256"),
              digest::digest(file = extension, algo = "sha256")))
  parquet <- duckhts_bench_stage_ancestry_parquet()
  con <- rduckhts_connect()
  spill <- tempfile("panel-spill-")
  dir.create(spill)
  dbExecute(con, sprintf("SET threads = %d", threads))
  dbExecute(con, sprintf("SET memory_limit = '%s'", memory_limit))
  dbExecute(con, sprintf("SET temp_directory = %s", dbQuoteString(con, spill)))
  rss_kib <- function() as.numeric(strsplit(trimws(grep("^VmHWM:", readLines("/proc/self/status"),
                                                        value = TRUE)), "[[:space:]]+")[[1]][2])
  if (quarters == 0L) {
    cat(sprintf("RESULT\t0\t0\t0\t0\t0\t\t0\t0\t0\t%.0f\n", rss_kib()))
    quit(save = "no")
  }
  columns <- dbGetQuery(con, sprintf("DESCRIBE SELECT * FROM read_parquet(%s)",
                                     dbQuoteString(con, parquet)))$column_name
  pcs <- grep("^PC[0-9]+$", columns, value = TRUE)
  groups <- setdiff(columns, c("chromosome", "position", "allele_a", "allele_b", pcs))
  quote_ids <- function(x) paste(vapply(x, function(v) as.character(dbQuoteIdentifier(con, v)), ""),
                                 collapse = ", ")
  loci_sql <- sprintf("SELECT * FROM read_parquet(%s) WHERE hash(chromosome, position) %% 4 < %d",
                      dbQuoteString(con, parquet), quarters)
  loci <- dbGetQuery(con, sprintf("SELECT count(*) AS n FROM (%s)", loci_sql))$n
  dbExecute(con, sprintf(paste0(
    "CREATE VIEW panel_reference AS SELECT chromosome::VARCHAR AS chromosome, position, ",
    "allele_a, allele_b, group_id, frequency FROM ",
    "(UNPIVOT (%s) ON %s INTO NAME group_id VALUE frequency)"), loci_sql, quote_ids(groups)))
  dbExecute(con, sprintf(paste0(
    "CREATE VIEW panel_loadings AS SELECT chromosome::VARCHAR AS chromosome, position, ",
    "allele_a, allele_b, CAST(substr(pc, 3) AS INTEGER) AS pc, loading FROM ",
    "(UNPIVOT (%s) ON %s INTO NAME pc VALUE loading)"), loci_sql, quote_ids(pcs)))
  candidates <- NULL
  if (identical(mode, "candidates")) {
    dbExecute(con, sprintf(paste0(
      "CREATE VIEW panel_candidates AS SELECT chromosome::VARCHAR AS region, position, ",
      "allele_a, allele_b FROM (%s)"), loci_sql))
    candidates <- "panel_candidates"
  }
  # Record DuckDB's per-statement peaks across every statement the builder runs.
  profile <- tempfile("panel-profile-", fileext = ".json")
  dbExecute(con, "PRAGMA enable_profiling = 'json'")
  dbExecute(con, sprintf("PRAGMA profiling_output = %s", dbQuoteString(con, profile)))
  peaks <- c(buffer = 0, temp = 0)
  record <- function() {
    if (!file.exists(profile)) return(invisible())
    m <- jsonlite::fromJSON(profile, simplifyVector = FALSE)
    peaks[["buffer"]] <<- max(peaks[["buffer"]], as.numeric(m$system_peak_buffer_memory))
    peaks[["temp"]] <<- max(peaks[["temp"]], as.numeric(m$system_peak_temp_dir_size))
  }
  suppressMessages({
    trace(DBI::dbExecute, exit = record, print = FALSE, where = asNamespace("DBI"))
    trace(DBI::dbGetQuery, exit = record, print = FALSE, where = asNamespace("DBI"))
  })
  started <- proc.time()[["elapsed"]]
  panel <- rduckhts_ancestry_panel(con, "panel_reference", "panel_loadings", "panel", "GRCh37",
                                   candidate_table = candidates, max_sites = max_sites)
  seconds <- proc.time()[["elapsed"]] - started
  suppressMessages({
    untrace(DBI::dbExecute, where = asNamespace("DBI"))
    untrace(DBI::dbGetQuery, where = asNamespace("DBI"))
  })
  dbExecute(con, "PRAGMA disable_profiling")
  shape <- dbGetQuery(con, "SELECT count(*) AS sites, count(DISTINCT region) AS contigs FROM panel")
  cat(sprintf("RESULT\t%d\t%d\t%d\t%d\t%d\t%s\t%.2f\t%.0f\t%.0f\t%.0f\n",
              loci, length(groups), length(pcs), shape$sites, shape$contigs, panel$panel_sha256,
              seconds, peaks[["buffer"]], peaks[["temp"]], rss_kib()))
  dbDisconnect(con, shutdown = TRUE)
  unlink(c(profile, spill), recursive = TRUE)
  quit(save = "no")
}

stopifnot(length(args) >= 2L)
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
out_path <- args[1]
libs <- vapply(strsplit(args[-1], "=", fixed = TRUE), `[`, "", 2)
names(libs) <- vapply(strsplit(args[-1], "=", fixed = TRUE), `[`, "", 1)
fields <- c("loci", "groups", "pcs", "sites", "contigs", "panel_sha256", "seconds",
            "peak_buffer_bytes", "peak_temp_bytes", "peak_rss_kib")
columns <- c("implementation", "mode", "quarters", "threads", "repetition", "memory_limit",
             "status", "package_code_sha256", "extension_sha256", fields)
# Rows are appended as they finish, so an interrupted matrix resumes where it stopped.
done <- if (file.exists(out_path)) utils::read.delim(out_path, colClasses = "character") else NULL
if (is.null(done)) cat(paste(columns, collapse = "\t"), "\n", sep = "", file = out_path)
run <- function(impl, quarters, mode, threads, rep) {
  key <- c(impl, mode, quarters, threads, rep)
  if (!is.null(done) && any(do.call(paste, done[columns[1:5]]) == paste(key, collapse = " "))) {
    return(invisible())
  }
  out <- suppressWarnings(system2("Rscript", shQuote(c(script, "--worker", libs[[impl]], quarters,
                                                       mode, threads)),
                                  stdout = TRUE, stderr = TRUE))
  build <- strsplit(grep("^BUILD\t", out, value = TRUE), "\t")[[1]][-1]
  stopifnot(length(build) == 2L)
  line <- grep("^RESULT\t", out, value = TRUE)
  if (length(line) == 1L) {
    values <- strsplit(line, "\t")[[1]][-1]
    status <- "ok"
  } else if (any(grepl("Out of Memory Error", out, fixed = TRUE))) {
    # A builder that exceeds the memory limit is a result, not a harness failure.
    values <- rep("NA", length(fields))
    status <- "out_of_memory"
  } else {
    stop(paste(out, collapse = "\n"))
  }
  row <- c(key, memory_limit, status, build, values)
  message(paste(row, collapse = " "))
  cat(paste(row, collapse = "\t"), "\n", sep = "", file = out_path, append = TRUE)
}
for (rep in 1:3) {
  for (impl in names(libs)) run(impl, 0L, "overhead", 1L, rep)
  for (threads in c(1L, 4L)) for (mode in c("spaced", "candidates")) for (quarters in c(1L, 2L, 4L)) {
    for (impl in names(libs)) run(impl, quarters, mode, threads, rep)
  }
}
