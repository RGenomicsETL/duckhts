args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: run_chr20_arm.R ARM REP", call. = FALSE)
arm <- args[[1L]]
rep <- args[[2L]]
if (!arm %in% c("pooled", "single_AFR", "single_AMR", "single_EUR", "single_EAS", "ancestry")) {
  stop("unknown arm", call. = FALSE)
}
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
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
DBI::dbExecute(connection, sprintf(
  "CREATE VIEW af_source AS SELECT * FROM read_parquet(%s)",
  q(file.path(cache, "chr20.af_reference_sites.parquet"))))
DBI::dbExecute(connection, sprintf(
  "CREATE VIEW ref_long AS SELECT * FROM read_parquet(%s)",
  q(file.path(cache, "chr20.reference_long.parquet"))))
q_data <- utils::read.csv(file.path(cache, "chr20.q.csv"), stringsAsFactors = FALSE)
DBI::dbWriteTable(connection, "q_table",
  q_data[c("sample_id", "group_id", "proportion")], overwrite = TRUE, temporary = TRUE)
all_ids <- readLines(file.path(cache, "chr20.children.present.txt"))
if (arm %in% c("ancestry", "pooled")) {
  ids <- all_ids
} else {
  population_data <- read.table(file.path(cache, "pedigree.txt"), header = TRUE,
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
      q(file.path(cache, "chr20.children.bcf")),
      q(paste(batch, collapse = ",")))
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
      q(file.path(cache, "chr20.children.bcf")), view,
      q(paste(batch, collapse = ",")))
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
out <- data.frame(arm = arm, repetition = as.integer(rep), result,
                  stringsAsFactors = FALSE)
utils::write.csv(out, file.path(cache,
  sprintf("chr20.arm_%s_%s.csv", arm, rep)), row.names = FALSE)
cat("arm", arm, "rep", rep, "samples", length(ids), "batches", length(batches),
    "segments", nrow(result), "elapsed", unname(timing[["elapsed"]]), "\n")
