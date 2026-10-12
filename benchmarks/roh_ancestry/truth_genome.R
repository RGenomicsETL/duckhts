# Genome-wide truth of the ROH ancestry evaluation: the rule of truth_chr20.R
# and truth_intervals.R on every autosome of the run. For each child and
# 100 kb window, `sites` counts the fully called diploid genotypes and
# `heterozygotes` the heterozygous ones. A window is eligible when it has at
# most one heterozygous call and at least one called genotype: a window
# without one is missing evidence, not homozygous evidence, so it breaks a
# run. A truth run is a merged run of eligible windows that covers at least
# 1 Mb of actual bases. The rule uses no allele frequencies, so it marks runs
# of low heterozygosity, not autozygosity.
#
# A chromosome's windows are reused only when their marker records the same
# SHA-256 for the extension, this script, the shared helpers, the pedigree,
# the children list, the FAI and the chromosome's children BCF. With chr20 in
# the run, its intervals must equal the amended chr20 truth committed in
# benchmarks/results/roh-ancestry/chr20/truth_intervals.csv.
#
# froh_assembly is truth length over the summed FAI length of the run's
# autosomes, gaps and centromeres included. No callable-base mask is
# registered, so it is a whole-assembly fraction and smaller than a FROH over
# callable bases.
source("benchmarks/roh_ancestry/genome_common.R")
# Base-pair totals print as integers, not in scientific notation.
options(scipen = 999)
window_bp <- 100000L
max_heterozygotes <- 1L
min_sites <- 1L
min_run_bp <- 1000000L
chromosomes <- genome_chromosomes()
lengths <- genome_lengths(chromosomes)
pedigree_path <- duckhtsbench::duckhts_bench_artifact_path("roh_ancestry_pedigree")
pedigree <- utils::read.table(pedigree_path, header = TRUE, stringsAsFactors = FALSE,
                              check.names = FALSE)
pedigree_children <- pedigree[pedigree$FatherID != "0" & pedigree$MotherID != "0", ]
children <- pedigree_children[pedigree_children$SampleID %in% readLines(genome_children_list),
                              c("SampleID", "Population", "Superpopulation")]
if (nrow(children) != 377L || anyDuplicated(children$SampleID)) {
  stop("pedigree child set does not contain 377 unique samples", call. = FALSE)
}
extension <- normalizePath("build/release/duckhts.duckdb_extension", mustWork = TRUE)
connection <- Rduckhts::rduckhts_connect(
  extension_path = extension,
  config = list(memory_limit = "12GB", threads = "4",
                temp_directory = file.path(genome_cache, "duckdb-tmp")))
quote_string <- function(value) as.character(DBI::dbQuoteString(connection, value))
DBI::dbWriteTable(connection, "truth_children", data.frame(
  sample_id = children$SampleID, population = children$Population,
  superpopulation = children$Superpopulation), temporary = TRUE)
DBI::dbWriteTable(connection, "truth_lengths", lengths, temporary = TRUE)
common_inputs <- c(
  extension = extension,
  driver = normalizePath("benchmarks/roh_ancestry/truth_genome.R", mustWork = TRUE),
  helpers = normalizePath("benchmarks/roh_ancestry/genome_common.R", mustWork = TRUE),
  pedigree = pedigree_path, children = genome_children_list,
  fai = normalizePath(genome_fai, mustWork = TRUE))
settings <- c(window_bp = window_bp, duckdb_r = as.character(utils::packageVersion("duckdb")))
# sites counts the child's fully called diploid genotypes in the window: a
# missing or partly missing genotype is no evidence either way. The input is
# biallelic.
called <- "('0/0', '0/1', '1/0', '1/1', '0|0', '0|1', '1|0', '1|1')"
heterozygous <- "('0/1', '1/0', '0|1', '1|0')"
window_paths <- character(length(chromosomes))
for (index in seq_along(chromosomes)) {
  chromosome <- chromosomes[[index]]
  length_bp <- lengths$length_bp[[index]]
  n_windows <- (length_bp - 1L) %/% window_bp + 1L
  bcf_path <- genome_children_bcf(chromosome)
  output <- file.path(genome_cache, sprintf("chr%d.truth_windows_genome.parquet", chromosome))
  window_paths[[index]] <- output
  identity <- genome_identity(c(common_inputs, bcf = bcf_path),
                              c(chromosome = chromosome, settings))
  if (genome_current(output, identity)) {
    cat("truth chromosome", chromosome, "reused windows with a matching identity\n")
    next
  }
  unlink(c(output, genome_marker(output)))
  DBI::dbExecute(connection, "DROP TABLE IF EXISTS truth_observed")
  DBI::dbExecute(connection, "DROP TABLE IF EXISTS truth_windows")
  DBI::dbExecute(connection, sprintf(paste0(
    "CREATE TEMP TABLE truth_observed AS ",
    "SELECT SAMPLE_ID AS sample_id, floor((POS - 1) / %d)::INTEGER AS window_id, ",
    "count(*) FILTER (WHERE FORMAT_GT IN %s)::INTEGER AS sites, ",
    "count(*) FILTER (WHERE FORMAT_GT IN %s)::INTEGER AS heterozygotes ",
    "FROM read_bcf(%s, tidy_format := true) GROUP BY SAMPLE_ID, window_id"),
    window_bp, called, heterozygous, quote_string(bcf_path)))
  outside <- DBI::dbGetQuery(connection, sprintf(
    "SELECT count(*) AS n FROM truth_observed WHERE window_id < 0 OR window_id >= %d",
    n_windows))$n
  if (outside != 0) {
    stop("chr", chromosome, " has records outside its FAI length", call. = FALSE)
  }
  DBI::dbExecute(connection, sprintf(
    "CREATE TEMP TABLE truth_windows AS SELECT range::INTEGER AS window_id FROM range(%d)",
    n_windows))
  temporary <- paste0(output, ".partial-", Sys.getpid(), ".parquet")
  DBI::dbExecute(connection, sprintf(paste0(
    "COPY (SELECT %d::INTEGER AS chromosome, c.sample_id, w.window_id, ",
    "coalesce(o.sites, 0)::INTEGER AS sites, ",
    "coalesce(o.heterozygotes, 0)::INTEGER AS heterozygotes ",
    "FROM truth_children c CROSS JOIN truth_windows w ",
    "LEFT JOIN truth_observed o USING (sample_id, window_id) ",
    "ORDER BY c.sample_id, w.window_id) TO %s (FORMAT PARQUET, COMPRESSION ZSTD)"),
    chromosome, quote_string(temporary)))
  written <- DBI::dbGetQuery(connection, sprintf(
    "SELECT count(*) AS n FROM read_parquet(%s)", quote_string(temporary)))$n
  if (written != nrow(children) * n_windows || !file.rename(temporary, output)) {
    unlink(temporary)
    stop("truth windows of chr", chromosome, " do not cover every child and window",
         call. = FALSE)
  }
  genome_record(output, identity)
  cat("truth chromosome", chromosome, "windows per child", n_windows, "parquet", output, "\n")
  flush.console()
}
path_list <- paste(quote_string(window_paths), collapse = ",")
windows_relation <- sprintf("read_parquet([%s])", path_list)
all_windows <- file.path(genome_cache, "autosome.truth_windows.parquet")
DBI::dbExecute(connection, sprintf(paste0(
  "COPY (SELECT * FROM %s ORDER BY chromosome, sample_id, window_id) TO %s ",
  "(FORMAT PARQUET, COMPRESSION ZSTD)"), windows_relation, quote_string(all_windows)))

# The heterozygote distribution over windows with a called genotype, and the
# windows without one, which no truth run can cross.
distribution <- DBI::dbGetQuery(connection, sprintf(paste0(
  "SELECT heterozygotes, count(*) AS child_windows, count(DISTINCT sample_id) AS children ",
  "FROM %s WHERE sites >= %d GROUP BY heterozygotes ORDER BY heterozygotes"),
  windows_relation, min_sites))
utils::write.csv(distribution, file.path(genome_results, "autosome_truth_window_distribution.csv"),
                 row.names = FALSE)
evidence <- DBI::dbGetQuery(connection, sprintf(paste0(
  "SELECT chromosome, count(*) AS child_windows, ",
  "count(*) FILTER (WHERE sites < %d) AS child_windows_without_called_genotype ",
  "FROM %s GROUP BY chromosome ORDER BY chromosome"), min_sites, windows_relation))
utils::write.csv(evidence, file.path(genome_results, "autosome_truth_window_evidence.csv"),
                 row.names = FALSE)
print(evidence)

DBI::dbExecute(connection, sprintf(paste0(
  "CREATE TEMP TABLE truth_runs AS ",
  "WITH eligible AS (SELECT w.chromosome, w.sample_id, w.window_id, l.length_bp, ",
  "least(%d, l.length_bp - w.window_id * %d) AS window_bp ",
  "FROM %s w JOIN truth_lengths l USING (chromosome) ",
  "WHERE w.heterozygotes <= %d AND w.sites >= %d), ",
  "numbered AS (SELECT *, window_id - row_number() OVER ",
  "(PARTITION BY chromosome, sample_id ORDER BY window_id) AS island FROM eligible), ",
  "runs AS (SELECT chromosome, sample_id, min(window_id) AS start_window, ",
  "max(window_id) AS end_window, count(*) AS windows, sum(window_bp) AS run_bp, ",
  "min(length_bp) AS length_bp FROM numbered GROUP BY chromosome, sample_id, island ",
  "HAVING sum(window_bp) >= %d) ",
  "SELECT chromosome, sample_id, start_window, end_window, windows, run_bp, ",
  "start_window * %d + 1 AS start, least((end_window + 1) * %d, length_bp) AS \"end\" ",
  "FROM runs"), window_bp, window_bp, windows_relation, max_heterozygotes, min_sites,
  min_run_bp, window_bp, window_bp))
intervals <- DBI::dbGetQuery(connection,
  "SELECT * FROM truth_runs ORDER BY chromosome, sample_id, start_window")
if (20L %in% chromosomes) {
  committed <- utils::read.csv(file.path(genome_results, "chr20", "truth_intervals.csv"),
                               stringsAsFactors = FALSE)
  keys <- c("sample_id", "start_window", "end_window", "windows", "run_bp", "start", "end")
  if (!isTRUE(all.equal(committed[keys], intervals[intervals$chromosome == 20L, keys],
                        check.attributes = FALSE))) {
    stop("chr20 truth intervals differ from the amended chr20 truth", call. = FALSE)
  }
}
intervals_path <- file.path(genome_cache, "autosome.truth_intervals.csv")
utils::write.csv(intervals, intervals_path, row.names = FALSE)
genome_record(intervals_path, c(genome_set_line(chromosomes), paste0(
  "windows_", chromosomes, "=", vapply(window_paths, genome_sha256, character(1L)))))

assembly_bp <- sum(lengths$length_bp)
summary <- DBI::dbGetQuery(connection, sprintf(paste0(
  "WITH run_summary AS (SELECT sample_id, count(*) AS truth_runs, ",
  "sum(run_bp) AS truth_bp, ",
  "count(*) FILTER (WHERE run_bp >= 2000000) AS runs_ge_2mb, ",
  "count(*) FILTER (WHERE run_bp >= 5000000) AS runs_ge_5mb ",
  "FROM truth_runs GROUP BY sample_id) ",
  "SELECT c.sample_id, c.population, coalesce(r.truth_runs, 0) AS truth_runs, ",
  "coalesce(r.truth_bp, 0) AS truth_bp, ",
  "coalesce(r.truth_bp, 0)::DOUBLE / %s AS froh_assembly, ",
  "coalesce(r.runs_ge_2mb, 0) AS runs_ge_2mb, ",
  "coalesce(r.runs_ge_5mb, 0) AS runs_ge_5mb ",
  "FROM truth_children c LEFT JOIN run_summary r USING (sample_id) ",
  "ORDER BY c.population, c.sample_id"),
  format(assembly_bp, scientific = FALSE, trim = TRUE)))
utils::write.csv(summary, file.path(genome_results, "autosome_truth_by_child.csv"),
                 row.names = FALSE)
population_summary <- do.call(rbind, lapply(split(summary, summary$population), function(x) {
  data.frame(population = x$population[[1L]], children = nrow(x),
             truth_runs = sum(x$truth_runs), truth_bp = sum(x$truth_bp),
             mean_froh_assembly = mean(x$froh_assembly), runs_ge_2mb = sum(x$runs_ge_2mb),
             runs_ge_5mb = sum(x$runs_ge_5mb))
}))
utils::write.csv(population_summary, file.path(genome_results, "autosome_truth_by_population.csv"),
                 row.names = FALSE)
cat("truth children", nrow(summary), "chromosomes", length(chromosomes),
    "intervals", nrow(intervals), "assembly bp", assembly_bp, "\n")
DBI::dbDisconnect(connection, shutdown = TRUE)
