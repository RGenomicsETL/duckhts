# Diagnosis of (a), written after its result: why the read-count decode has
# fewer bases in runs than the VCF genotype decode at gt_error 30.
#
# It answers three questions with the staged inputs of run.R:
#   1. Does the genotype decode approach the read-count decode as gt_error rises?
#   2. Does seq_error move the read-count decode?
#   3. Do the genotypes and the read counts disagree at the sites, and what do
#      the parts of genotype runs that the read-count decode loses contain?
#
# Run from the repository root, after run.R's inputs are staged.

source("benchmarks/roh_counts_validation/run.R")

roh_diagnosis_gt_errors <- c(30, 60, 100, 150, 250)
roh_diagnosis_seq_errors <- c(1e-4, 1e-3, 1e-2, 1e-1)

# Read-count classes of one site, by the reads of the rarer allele.
roh_diagnosis_read_class_sql <- paste(
  "CASE WHEN ref_count + alt_count = 0 THEN 'no reads'",
  "WHEN least(ref_count, alt_count) = 0 THEN 'one allele only'",
  "WHEN least(ref_count, alt_count) < 0.1 * (ref_count + alt_count) THEN 'minor allele under 10%'",
  "WHEN least(ref_count, alt_count) < 0.25 * (ref_count + alt_count) THEN 'minor allele 10% to 25%'",
  "ELSE 'minor allele 25% or more' END")

roh_diagnosis_run <- function(extension, output_dir, threads = 2L,
                              declaration = roh_validation_declaration) {
  samples <- declaration$samples
  con <- roh_validation_connect(extension, threads)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  q <- function(value) as.character(DBI::dbQuoteString(con, value))
  count_paths <- vapply(declaration$count_artifacts[samples],
                        duckhtsbench::duckhts_bench_artifact_path, character(1L))
  genotype_path <- duckhtsbench::duckhts_bench_artifact_path(declaration$genotype_artifact)

  DBI::dbExecute(con, sprintf(paste0(
    "CREATE TABLE counts AS SELECT sample_id, chrom, pos, ref, alt, ref_count, alt_count, af ",
    "FROM read_parquet([%s])"), paste(q(count_paths), collapse = ", ")))
  DBI::dbExecute(con, sprintf(
    "CREATE VIEW af_sites AS SELECT chrom, pos, ref, alt, af FROM counts WHERE sample_id = %s",
    q(samples[[1L]])))
  span_bases <- DBI::dbGetQuery(con, "SELECT max(pos) - min(pos) + 1 AS span FROM counts")$span

  genotype_runs <- function(gt_error) DBI::dbGetQuery(con, sprintf(
    "SELECT * FROM duckhts_roh_af_table(%s, 'af_sites', gt_error := %s)",
    q(genotype_path), format(gt_error)))
  count_runs <- function(seq_error) DBI::dbGetQuery(con, sprintf(
    "SELECT * FROM duckhts_roh_counts('counts', seq_error := %s, contamination := 0.0)",
    format(seq_error, scientific = FALSE)))
  compare <- function(test, reference, label) {
    do.call(rbind, lapply(samples, function(sample) {
      metrics <- roh_run_comparison_by_class(
        roh_validation_runs_of(test, sample), roh_validation_runs_of(reference, sample),
        span_bases, declaration$long_run_bases)
      data.frame(label, sample = sample, metrics[metrics$run_class == "all", ])
    }))
  }

  # 1. The read-count decode (default seq_error) against the genotype decode,
  #    by gt_error.
  default_count_runs <- count_runs(0.001)
  by_gt_error <- do.call(rbind, lapply(roh_diagnosis_gt_errors, function(gt_error) {
    compare(default_count_runs, genotype_runs(gt_error), data.frame(gt_error = gt_error))
  }))

  # 2. The read-count decode by seq_error, against the genotype decode at the
  #    declared gt_error.
  declared_genotype_runs <- genotype_runs(declaration$gt_error)
  by_seq_error <- do.call(rbind, lapply(roh_diagnosis_seq_errors, function(seq_error) {
    compare(count_runs(seq_error), declared_genotype_runs, data.frame(seq_error = seq_error))
  }))

  # 3. Genotype against read-count class at every site, and inside the
  #    genotype runs split by whether the read-count decode keeps the site.
  DBI::dbExecute(con, sprintf(paste0(
    "CREATE TABLE paired AS SELECT c.sample_id, c.pos, ",
    "CASE g.FORMAT_GT WHEN '0|0' THEN 'homozygous' WHEN '1|1' THEN 'homozygous' ",
    "WHEN '0|1' THEN 'heterozygous' WHEN '1|0' THEN 'heterozygous' ELSE 'other' END AS genotype, ",
    "%s AS read_class FROM counts AS c ",
    "JOIN read_bcf(%s, tidy_format := true) AS g ON g.SAMPLE_ID = c.sample_id AND g.POS = c.pos"),
    roh_diagnosis_read_class_sql, q(genotype_path)))
  site_classes <- DBI::dbGetQuery(con, paste(
    "SELECT genotype, read_class, count(*) AS sites FROM paired",
    "GROUP BY genotype, read_class ORDER BY genotype, read_class"))
  DBI::dbWriteTable(con, "genotype_runs", declared_genotype_runs)
  DBI::dbWriteTable(con, "count_runs", default_count_runs)
  inside_runs <- DBI::dbGetQuery(con, paste(
    "WITH inside AS (SELECT p.*, EXISTS (SELECT 1 FROM count_runs AS r",
    "WHERE r.sample = p.sample_id AND p.pos BETWEEN r.start AND r.\"end\") AS kept",
    "FROM paired AS p JOIN genotype_runs AS v",
    "ON v.sample = p.sample_id AND p.pos BETWEEN v.start AND v.\"end\")",
    "SELECT sample_id AS sample, kept, count(*) AS sites,",
    "count(*) FILTER (WHERE genotype = 'heterozygous') AS heterozygous_genotypes,",
    "count(*) FILTER (WHERE genotype = 'homozygous' AND read_class = 'minor allele 25% or more')",
    "AS homozygous_with_balanced_reads",
    "FROM inside GROUP BY sample_id, kept ORDER BY sample_id, kept"))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(by_gt_error, file.path(output_dir, "diagnosis_by_gt_error.csv"), row.names = FALSE)
  utils::write.csv(by_seq_error, file.path(output_dir, "diagnosis_by_seq_error.csv"), row.names = FALSE)
  utils::write.csv(site_classes, file.path(output_dir, "diagnosis_site_classes.csv"), row.names = FALSE)
  utils::write.csv(inside_runs, file.path(output_dir, "diagnosis_inside_genotype_runs.csv"), row.names = FALSE)
  list(by_gt_error = by_gt_error, by_seq_error = by_seq_error,
       site_classes = site_classes, inside_runs = inside_runs)
}
