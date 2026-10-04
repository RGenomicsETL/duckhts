# Run one ROH arm on the synthetic chr20 children. The macro calls, the
# frequency clamp, gt_error, the batch sizes and the DuckDB settings are those
# of benchmarks/roh_ancestry/run_chr20_arm.R; only the BCF differs.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: run_arm.R ARM", call. = FALSE)
arm <- args[[1L]]
if (!arm %in% c("pooled", "single_AFR", "single_AMR", "single_EUR", "single_EAS", "ancestry")) {
  stop("unknown arm", call. = FALSE)
}
artifact <- duckhtsbench::duckhts_bench_artifact_path
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-synthetic")
dir.create(file.path(cache, "duckdb-tmp"), recursive = TRUE, showWarnings = FALSE)
driver <- duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE,
  allow_extensions = TRUE, config = list(memory_limit = "12GB", threads = "4",
    temp_directory = file.path(cache, "duckdb-tmp"),
    allow_unsigned_extensions = "true", autoinstall_known_extensions = "false",
    autoload_known_extensions = "false"))
connection <- DBI::dbConnect(driver)
DBI::dbExecute(connection, sprintf("LOAD %s",
  as.character(DBI::dbQuoteString(connection,
    normalizePath("build/release/duckhts.duckdb_extension")))))
q <- function(x) as.character(DBI::dbQuoteString(connection, x))
bcf <- artifact("roh_synthetic_chr20_children_bcf")
DBI::dbExecute(connection, sprintf(
  "CREATE VIEW af_source AS SELECT * FROM read_parquet(%s)",
  q(artifact("roh_ancestry_chr20_af_sites"))))
DBI::dbExecute(connection, sprintf(
  "CREATE VIEW ref_long AS SELECT * FROM read_parquet(%s)",
  q(artifact("roh_ancestry_chr20_reference_long"))))
q_data <- utils::read.csv(artifact("roh_ancestry_chr20_q"), stringsAsFactors = FALSE)
DBI::dbWriteTable(connection, "q_table",
  q_data[c("sample_id", "group_id", "proportion")], overwrite = TRUE, temporary = TRUE)
all_ids <- DBI::dbGetQuery(connection, sprintf(
  "SELECT DISTINCT sample_id FROM read_csv(%s) ORDER BY sample_id",
  q(artifact("roh_synthetic_chr20_truth"))))$sample_id
if (arm %in% c("ancestry", "pooled")) {
  ids <- all_ids
} else {
  population_data <- read.table(artifact("roh_ancestry_pedigree"), header = TRUE,
                                stringsAsFactors = FALSE)
  superpopulation <- sub("^single_", "", arm)
  populations <- switch(superpopulation, AFR = c("ACB", "ASW", "ESN", "YRI"),
                        AMR = c("CLM", "MXL", "PEL", "PUR"),
                        EUR = "CEU", EAS = "CHS")
  ids <- population_data$SampleID[
    population_data$Population %in% populations & population_data$SampleID %in% all_ids]
}
batch_size <- if (arm == "ancestry") 16L else 64L
batches <- split(ids, ceiling(seq_along(ids) / batch_size))
if (arm == "ancestry") {
  run_batch <- function(batch) {
    sql <- sprintf(paste0(
      "SELECT * FROM duckhts_roh_ancestry(%s, 'ref_long', 'q_table', ",
      "af_clamp := 1e-3, gt_error := 30, samples := %s)"),
      q(bcf), q(paste(batch, collapse = ",")))
    DBI::dbGetQuery(connection, sql)
  }
} else {
  field <- switch(arm, pooled = "INFO_AF", single_AFR = "INFO_AF_AFR",
                  single_AMR = "INFO_AF_AMR", single_EUR = "INFO_AF_EUR",
                  single_EAS = "INFO_AF_EAS")
  view <- paste0("af_", tolower(arm))
  DBI::dbExecute(connection, sprintf(paste0(
    "CREATE VIEW %s AS SELECT chrom, pos, ref, alt, ",
    "greatest(1e-3, least(0.999, list_extract(%s, 1)::DOUBLE)) AS af FROM af_source"),
    view, field))
  run_batch <- function(batch) {
    sql <- sprintf(paste0(
      "SELECT * FROM duckhts_roh_af_table(%s, '%s', gt_error := 30, samples := %s)"),
      q(bcf), view, q(paste(batch, collapse = ",")))
    DBI::dbGetQuery(connection, sql)
  }
}
timing <- system.time({
  results <- lapply(seq_along(batches), function(index) {
    cat("arm", arm, "batch", index, "of", length(batches), "samples",
        length(batches[[index]]), "\n")
    flush.console()
    run_batch(batches[[index]])
  })
})
result <- do.call(rbind, results)
out <- data.frame(arm = arm, result, stringsAsFactors = FALSE)
utils::write.csv(out, file.path(cache, sprintf("chr20.arm_%s.csv", arm)), row.names = FALSE)
cat("arm", arm, "samples", length(ids), "batches", length(batches),
    "segments", nrow(result), "elapsed", unname(timing[["elapsed"]]), "\n")
DBI::dbDisconnect(connection, shutdown = TRUE)
