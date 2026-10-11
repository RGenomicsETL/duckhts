# GRCh38 autosome lengths: lines 1-22 of the Ensembl 116 primary-assembly FASTA index
# (Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai, SHA-256
# 0998f61682f4041b11f0d156e1db6dae3e4c743e26643a3f45ea7faea70cb604), committed so the
# truth windows and FROH denominators do not depend on an unregistered local file.
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
reference_fai <- "benchmarks/roh_ancestry/grch38_autosomes.fai"
fai <- utils::read.delim(reference_fai, header = FALSE, stringsAsFactors = FALSE,
                         colClasses = c("character", "integer", "NULL", "NULL", "NULL"))
chromosomes <- 1:22
lengths <- data.frame(chromosome = chromosomes,
  length_bp = fai$V2[match(as.character(chromosomes), fai$V1)])
if (anyNA(lengths$length_bp)) stop("GRCh38 FAI lacks an autosome", call. = FALSE)
pedigree <- utils::read.table(file.path(cache, "pedigree.txt"), header = TRUE,
  stringsAsFactors = FALSE, check.names = FALSE)
pedigree_children <- pedigree[pedigree$FatherID != "0" & pedigree$MotherID != "0", ]
children <- pedigree_children[pedigree_children$SampleID %in%
  readLines(file.path(cache, "chr20.children.present.txt")),
  c("SampleID", "Population", "Superpopulation")]
if (nrow(children) != 377L || anyDuplicated(children$SampleID)) {
  stop("pedigree child set does not contain 377 unique samples", call. = FALSE)
}
connection <- Rduckhts::rduckhts_connect(
  extension_path = normalizePath("build/release/duckhts.duckdb_extension"),
  config = list(memory_limit = "12GB", threads = "4",
                temp_directory = file.path(cache, "duckdb-tmp")))
quote_string <- function(value) as.character(DBI::dbQuoteString(connection, value))
DBI::dbWriteTable(connection, "truth_children", data.frame(
  sample_id = children$SampleID, population = children$Population,
  superpopulation = children$Superpopulation), temporary = TRUE)
DBI::dbWriteTable(connection, "truth_lengths", lengths, temporary = TRUE)
window_paths <- character(length(chromosomes))
for (index in seq_along(chromosomes)) {
  chromosome <- chromosomes[[index]]
  chromosome_length <- lengths$length_bp[[index]]
  path <- file.path(cache, sprintf("chr%d.children.bcf", chromosome))
  output <- if (chromosome == 20L) {
    file.path(cache, "chr20.truth_windows_genome.parquet")
  } else {
    file.path(cache, sprintf("chr%d.truth_windows.parquet", chromosome))
  }
  if (chromosome == 20L && file.exists(file.path(cache, "chr20.truth_windows.parquet"))) {
    reused <- file.path(cache, "chr20.truth_windows.parquet")
    expected_windows <- 377L * as.integer(ceiling(chromosome_length / 100000))
    observed_windows <- DBI::dbGetQuery(connection, sprintf(
      "SELECT count(*) AS n FROM read_parquet(%s)", quote_string(reused)))$n
    if (observed_windows != expected_windows) {
      stop("cached chr20 truth windows have an unexpected row count", call. = FALSE)
    }
    unlink(output)
    DBI::dbExecute(connection, sprintf(paste0(
      "COPY (SELECT 20::INTEGER AS chromosome, sample_id, window_id, heterozygotes ",
      "FROM read_parquet(%s)) TO %s (FORMAT PARQUET, COMPRESSION ZSTD)"),
      quote_string(reused), quote_string(output)))
    window_paths[[index]] <- output
    cat("truth chromosome 20 reused completed windows", observed_windows, "rows\n")
    next
  }
  expected_windows <- 377L * as.integer(ceiling(chromosome_length / 100000))
  if (file.exists(output)) {
    observed_windows <- DBI::dbGetQuery(connection, sprintf(
      "SELECT count(*) AS n FROM read_parquet(%s)", quote_string(output)))$n
    if (observed_windows == expected_windows) {
      window_paths[[index]] <- output
      cat("truth chromosome", chromosome, "reused completed windows",
          observed_windows, "rows\n")
      next
    }
    unlink(output)
  }
  DBI::dbExecute(connection, "DROP TABLE IF EXISTS truth_observed")
  DBI::dbExecute(connection, "DROP TABLE IF EXISTS truth_windows")
  DBI::dbExecute(connection, sprintf(paste0(
    "CREATE TEMP TABLE truth_observed AS ",
    "SELECT SAMPLE_ID AS sample_id, floor((POS - 1) / 100000)::INTEGER AS window_id, ",
    "count(*) FILTER (WHERE FORMAT_GT IN ('0/1','1/0','0|1','1|0'))::INTEGER AS heterozygotes ",
    "FROM read_bcf(%s, tidy_format := true) GROUP BY SAMPLE_ID, window_id"),
    quote_string(path)))
  n_windows <- as.integer(ceiling(chromosome_length / 100000))
  DBI::dbExecute(connection, sprintf(
    "CREATE TEMP TABLE truth_windows AS SELECT range::INTEGER AS window_id FROM range(%d)",
    n_windows))
  sql <- sprintf(paste0(
    "COPY (SELECT %d::INTEGER AS chromosome, c.sample_id, w.window_id, ",
    "coalesce(o.heterozygotes, 0)::INTEGER AS heterozygotes ",
    "FROM truth_children c CROSS JOIN truth_windows w ",
    "LEFT JOIN truth_observed o USING(sample_id, window_id) ",
    "ORDER BY c.sample_id, w.window_id) TO %s ",
    "(FORMAT PARQUET, COMPRESSION ZSTD)"), chromosome, quote_string(output))
  DBI::dbExecute(connection, sql)
  window_paths[[index]] <- output
  cat("truth chromosome", chromosome, "windows per child", n_windows,
      "parquet", output, "\n")
  flush.console()
}
path_list <- paste(quote_string(window_paths), collapse = ",")
windows_relation <- sprintf("read_parquet([%s])", path_list)
all_windows <- file.path(cache, "autosome.truth_windows.parquet")
DBI::dbExecute(connection, sprintf(paste0(
  "COPY (SELECT * FROM %s ORDER BY chromosome, sample_id, window_id) TO %s ",
  "(FORMAT PARQUET, COMPRESSION ZSTD)"), windows_relation, quote_string(all_windows)))
distribution <- DBI::dbGetQuery(connection, sprintf(paste0(
  "SELECT heterozygotes, count(*) AS child_windows, count(DISTINCT sample_id) AS children ",
  "FROM %s GROUP BY heterozygotes ORDER BY heterozygotes"), windows_relation))
utils::write.csv(distribution,
  "benchmarks/results/roh-ancestry/autosome_truth_window_distribution.csv",
  row.names = FALSE)
print(distribution)
run_sql <- sprintf(paste0(
  "CREATE TEMP TABLE truth_runs AS ",
  "WITH windows AS (SELECT w.chromosome, w.sample_id, w.window_id, ",
  "w.heterozygotes, l.length_bp, least(100000, l.length_bp - w.window_id * 100000) AS window_bp ",
  "FROM %s w JOIN truth_lengths l USING(chromosome) ",
  "WHERE w.heterozygotes <= 1), ",
  "numbered AS (SELECT *, window_id - row_number() OVER ",
  "(PARTITION BY chromosome, sample_id ORDER BY window_id) AS island FROM windows), ",
  "runs AS (SELECT chromosome, sample_id, min(window_id) AS start_window, ",
  "max(window_id) AS end_window, count(*) AS windows, sum(window_bp) AS run_bp, ",
  "min(length_bp) AS length_bp FROM numbered GROUP BY chromosome, sample_id, island ",
  "HAVING sum(window_bp) >= 1000000) ",
  "SELECT chromosome, sample_id, start_window, end_window, windows, run_bp, ",
  "start_window * 100000 + 1 AS start, least((end_window + 1) * 100000, length_bp) AS end ",
  "FROM runs"), windows_relation)
DBI::dbExecute(connection, run_sql)
intervals_path <- file.path(cache, "autosome.truth_intervals.csv")
intervals <- DBI::dbGetQuery(connection,
  "SELECT * FROM truth_runs ORDER BY chromosome, sample_id, start_window")
previous_chr20 <- utils::read.csv(file.path(cache, "chr20.truth_intervals.csv"),
                                  stringsAsFactors = FALSE)
observed_chr20 <- intervals[intervals$chromosome == 20L, ]
chr20_keys <- c("sample_id", "start_window", "end_window", "windows", "start", "end")
if (!isTRUE(all.equal(previous_chr20[chr20_keys], observed_chr20[chr20_keys],
                      check.attributes = FALSE))) {
  stop("reused chr20 truth intervals differ from the prior evaluation", call. = FALSE)
}
utils::write.csv(intervals, intervals_path, row.names = FALSE)
autosome_length <- sum(lengths$length_bp)
summary <- DBI::dbGetQuery(connection, sprintf(paste0(
  "WITH run_summary AS (SELECT sample_id, count(*) AS truth_runs, ",
  "sum(run_bp) AS truth_bp, ",
  "count(*) FILTER (WHERE run_bp >= 2000000) AS runs_ge_2mb, ",
  "count(*) FILTER (WHERE run_bp >= 5000000) AS runs_ge_5mb ",
  "FROM truth_runs GROUP BY sample_id) ",
  "SELECT c.sample_id, c.population, coalesce(r.truth_runs,0) AS truth_runs, ",
  "coalesce(r.truth_bp,0) AS truth_bp, ",
  "coalesce(r.truth_bp,0)::DOUBLE / %s AS froh, ",
  "coalesce(r.runs_ge_2mb,0) AS runs_ge_2mb, ",
  "coalesce(r.runs_ge_5mb,0) AS runs_ge_5mb ",
  "FROM truth_children c LEFT JOIN run_summary r USING(sample_id) ",
  "ORDER BY c.population,c.sample_id"),
  format(autosome_length, scientific = FALSE, trim = TRUE)))
utils::write.csv(summary,
  "benchmarks/results/roh-ancestry/autosome_truth_by_child.csv", row.names = FALSE)
population_summary <- DBI::dbGetQuery(connection, sprintf(paste0(
  "SELECT population, count(*) AS children, sum(truth_runs) AS truth_runs, ",
  "sum(truth_bp) AS truth_bp, avg(froh) AS mean_froh, ",
  "sum(runs_ge_2mb) AS runs_ge_2mb, sum(runs_ge_5mb) AS runs_ge_5mb ",
  "FROM read_csv_auto(%s) GROUP BY population ORDER BY population"),
  quote_string("benchmarks/results/roh-ancestry/autosome_truth_by_child.csv")))
utils::write.csv(population_summary,
  "benchmarks/results/roh-ancestry/autosome_truth_by_population.csv", row.names = FALSE)
cat("truth children", nrow(summary), "intervals", nrow(intervals),
    "total windows", sum(distribution$child_windows),
    "1Mb+ runs", sum(summary$truth_runs), "\n")
DBI::dbDisconnect(connection, shutdown = TRUE)
