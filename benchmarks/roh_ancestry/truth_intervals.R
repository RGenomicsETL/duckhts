# chr20 truth intervals from the truth windows written by truth_chr20.R.
#
# A truth run is a merged run of 100 kb windows with at most one heterozygous call,
# covering at least 1 Mb of actual bases. The amended rule also requires each window
# to hold at least `min_sites` called genotypes: a window without them is missing
# evidence, not homozygous evidence, so it is ineligible and breaks a run. The
# preregistered rule, which counted such windows as homozygous, is kept as retained
# evidence. A sensitivity rule fixed before it was computed requires
# `sensitivity_min_sites` called genotypes per window, after PLINK's 50-SNP
# --homozyg-window-snp default; it is not PLINK's algorithm.
library(DBI)
# Base-pair totals print as integers, not in scientific notation.
options(scipen = 999)
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
results <- "benchmarks/results/roh-ancestry/chr20"
windows_path <- duckhtsbench::duckhts_bench_artifact_path("roh_ancestry_chr20_truth_windows")
length_bp <- 64444167L
min_sites <- 1L
sensitivity_min_sites <- 50L
max_heterozygotes <- 1L
min_run_bp <- 1000000L

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:", config = list(threads = 2L))
q <- function(x) as.character(dbQuoteString(con, x))
dbExecute(con, sprintf(paste0(
  "CREATE TEMP TABLE windows AS SELECT sample_id, window_id, sites, heterozygotes, ",
  "least(100000, %d - window_id * 100000) AS window_bp FROM read_parquet(%s)"),
  length_bp, q(windows_path)))

truth_runs <- function(eligible) {
  dbGetQuery(con, sprintf(paste0(
    "WITH eligible AS (SELECT * FROM windows WHERE %s), ",
    "numbered AS (SELECT *, window_id - row_number() OVER ",
    "(PARTITION BY sample_id ORDER BY window_id) AS island FROM eligible), ",
    "runs AS (SELECT sample_id, min(window_id) AS start_window, max(window_id) AS end_window, ",
    "count(*) AS windows, sum(window_bp) AS run_bp FROM numbered ",
    "GROUP BY sample_id, island HAVING sum(window_bp) >= %d) ",
    "SELECT sample_id, start_window, end_window, windows, run_bp, ",
    "start_window * 100000 + 1 AS start, least((end_window + 1) * 100000, %d) AS \"end\" ",
    "FROM runs ORDER BY sample_id, start_window"), eligible, min_run_bp, length_bp))
}
preregistered <- truth_runs(sprintf("heterozygotes <= %d", max_heterozygotes))
amended <- truth_runs(sprintf("heterozygotes <= %d AND sites >= %d",
                              max_heterozygotes, min_sites))
write.csv(amended, file.path(cache, "chr20.truth_intervals.csv"), row.names = FALSE)
write.csv(preregistered, file.path(results, "truth_intervals_preregistered.csv"), row.names = FALSE)
write.csv(amended, file.path(results, "truth_intervals.csv"), row.names = FALSE)
sensitivity <- truth_runs(sprintf("heterozygotes <= %d AND sites >= %d",
                                  max_heterozygotes, sensitivity_min_sites))
interval_key <- c("sample_id", "start_window", "end_window")
same_as_amended <- function(runs) {
  isTRUE(all.equal(runs[interval_key], amended[interval_key], check.attributes = FALSE))
}
rule_summary <- function(rule, minimum, runs) {
  data.frame(rule = rule, min_called_genotypes = minimum, intervals = nrow(runs),
             children = length(unique(runs$sample_id)), truth_bp = sum(runs$run_bp),
             same_intervals_as_amended = same_as_amended(runs))
}
write.csv(rbind(rule_summary("preregistered", 0L, preregistered),
                rule_summary("amended", min_sites, amended),
                rule_summary("sensitivity", sensitivity_min_sites, sensitivity)),
          file.path(results, "truth_rule_comparison.csv"), row.names = FALSE)

# The window distribution over windows with enough records, and the windows without.
distribution <- dbGetQuery(con, sprintf(paste0(
  "SELECT heterozygotes, count(*) AS child_windows FROM windows WHERE sites >= %d ",
  "GROUP BY heterozygotes ORDER BY heterozygotes"), min_sites))
write.csv(distribution, file.path(results, "truth_window_distribution.csv"), row.names = FALSE)
uncalled_bin <- if (min_sites == 1L) {
  "no called genotype"
} else {
  sprintf("fewer than %d called genotypes", min_sites)
}
bins <- dbGetQuery(con, sprintf(paste0(
  "SELECT bin, count(*) AS child_windows FROM (SELECT CASE ",
  "WHEN sites < %d THEN '%s' ",
  "WHEN heterozygotes = 0 THEN '0' WHEN heterozygotes = 1 THEN '1' ",
  "WHEN heterozygotes <= 4 THEN '2-4' WHEN heterozygotes <= 9 THEN '5-9' ",
  "WHEN heterozygotes <= 19 THEN '10-19' WHEN heterozygotes <= 49 THEN '20-49' ",
  "ELSE '50+' END AS bin, CASE WHEN sites < %d THEN -1 ELSE heterozygotes END AS ord ",
  "FROM windows) GROUP BY bin ORDER BY min(ord)"), min_sites, uncalled_bin, min_sites))
write.csv(bins, file.path(results, "truth_window_bins.csv"), row.names = FALSE)

pedigree <- read.table(file.path(cache, "pedigree.txt"), header = TRUE, stringsAsFactors = FALSE)
children <- readLines(file.path(cache, "chr20.children.present.txt"))
population <- pedigree$Population[match(children, pedigree$SampleID)]
if (anyNA(population)) stop("every chr20 child needs a pedigree population", call. = FALSE)
by_population <- function(runs) {
  runs$Population <- population[match(runs$sample_id, children)]
  do.call(rbind, lapply(sort(unique(population)), function(p) {
    selected <- runs[runs$Population == p, , drop = FALSE]
    data.frame(Population = p, children = sum(population == p),
               truth_children = length(unique(selected$sample_id)),
               truth_runs = nrow(selected), truth_bp = sum(selected$run_bp))
  }))
}
write.csv(by_population(amended), file.path(results, "truth_runs_by_population.csv"),
          row.names = FALSE)
write.csv(by_population(preregistered),
          file.path(results, "truth_runs_by_population_preregistered.csv"), row.names = FALSE)
cat("truth intervals: preregistered", nrow(preregistered), "amended", nrow(amended), "\n")
dbDisconnect(con, shutdown = TRUE)
