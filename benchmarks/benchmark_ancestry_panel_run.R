# Panel-building scaling for rduckhts_ancestry_panel() on bigsnpr's staged reference.
#
# Driver:  Rscript benchmarks/benchmark_ancestry_panel_run.R <old_lib> <new_lib> <out.tsv>
# Worker:  Rscript benchmarks/benchmark_ancestry_panel_run.R --worker <lib> <quarters> <mode> <max_sites>
#
# Each worker is a fresh single-thread process. The reference is subsampled to 1, 2
# or 4 quarters of its loci by hash, expanded to the long frequency and loading
# relations the panel builder reads, and consumed in full into a committed panel.
# Modes: "spaced" uses the builder's window spacing; "candidates" offers every locus
# as a candidate. Peak RSS comes from GNU time.
args <- commandArgs(TRUE)

if (identical(args[1], "--worker")) {
  .libPaths(c(args[2], .libPaths()))
  quarters <- as.integer(args[3]); mode <- args[4]; max_sites <- as.integer(args[5])
  suppressPackageStartupMessages({
    library(DBI); library(Rduckhts); library(duckhtsbench)
  })
  parquet <- duckhts_bench_stage_ancestry_parquet()
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, "SET threads = 1")
  columns <- dbGetQuery(con, sprintf("DESCRIBE SELECT * FROM read_parquet(%s)",
                                     dbQuoteString(con, parquet)))$column_name
  pcs <- grep("^PC[0-9]+$", columns, value = TRUE)
  groups <- setdiff(columns, c("chromosome", "position", "allele_a", "allele_b", pcs))
  quote_ids <- function(x) paste(vapply(x, function(v) as.character(dbQuoteIdentifier(con, v)), ""),
                                 collapse = ", ")
  loci <- sprintf("SELECT * FROM read_parquet(%s) WHERE hash(chromosome, position) %% 4 < %d",
                  dbQuoteString(con, parquet), quarters)
  dbExecute(con, sprintf(paste0(
    "CREATE VIEW panel_reference AS SELECT chromosome::VARCHAR AS chromosome, position, ",
    "allele_a, allele_b, group_id, frequency FROM ",
    "(UNPIVOT (%s) ON %s INTO NAME group_id VALUE frequency)"), loci, quote_ids(groups)))
  dbExecute(con, sprintf(paste0(
    "CREATE VIEW panel_loadings AS SELECT chromosome::VARCHAR AS chromosome, position, ",
    "allele_a, allele_b, CAST(substr(pc, 3) AS INTEGER) AS pc, loading FROM ",
    "(UNPIVOT (%s) ON %s INTO NAME pc VALUE loading)"), loci, quote_ids(pcs)))
  candidates <- NULL
  if (identical(mode, "candidates")) {
    dbExecute(con, sprintf(paste0(
      "CREATE VIEW panel_candidates AS SELECT chromosome::VARCHAR AS region, position, ",
      "allele_a, allele_b FROM (%s)"), loci))
    candidates <- "panel_candidates"
  }
  started <- proc.time()[["elapsed"]]
  rduckhts_ancestry_panel(con, "panel_reference", "panel_loadings", "panel", "GRCh37",
                          candidate_table = candidates, max_sites = max_sites)
  seconds <- proc.time()[["elapsed"]] - started
  shape <- dbGetQuery(con, "SELECT count(*) AS sites, count(DISTINCT region) AS contigs FROM panel")
  cat(sprintf("RESULT\t%d\t%d\t%.2f\n", shape$sites, shape$contigs, seconds))
  quit(save = "no")
}

stopifnot(length(args) == 3L, file.exists("/usr/bin/time"))
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
libs <- c(old = args[1], new = args[2])
rows <- list()
for (rep in 1:3) for (mode in c("spaced", "candidates")) for (quarters in c(1L, 2L, 4L)) {
  for (impl in names(libs)) {
    out <- system2("/usr/bin/time", shQuote(c("-f", "RSS\t%M", "Rscript", script, "--worker",
                                             libs[[impl]], quarters, mode, 17000L)),
                   stdout = TRUE, stderr = TRUE)
    result <- strsplit(grep("^RESULT\t", out, value = TRUE), "\t")[[1]]
    rss <- strsplit(grep("^RSS\t", out, value = TRUE), "\t")[[1]]
    stopifnot(length(result) == 4L, length(rss) == 2L)
    rows[[length(rows) + 1L]] <- data.frame(
      implementation = impl, mode = mode, quarters = quarters, repetition = rep,
      sites = as.integer(result[2]), contigs = as.integer(result[3]),
      seconds = as.numeric(result[4]), peak_rss_kb = as.numeric(rss[2]))
    message(paste(unlist(rows[[length(rows)]]), collapse = " "))
  }
}
utils::write.table(do.call(rbind, rows), args[3], sep = "\t", quote = FALSE, row.names = FALSE)
