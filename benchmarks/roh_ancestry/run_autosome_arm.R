args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: run_autosome_arm.R ARM REP", call. = FALSE)
arm <- args[[1L]]
repetition <- as.integer(args[[2L]])
if (length(repetition) != 1L || is.na(repetition) || repetition < 1L || repetition > 3L) {
  stop("repetition must be 1, 2 or 3", call. = FALSE)
}
if (!arm %in% c("pooled", "single_AFR", "single_AMR", "single_EUR",
                "single_EAS", "ancestry")) {
  stop("unknown arm", call. = FALSE)
}
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
chromosomes <- if (arm == "ancestry") 1:22 else setdiff(1:22, 20L)
driver <- duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE,
  allow_extensions = TRUE, config = list(memory_limit = "12GB", threads = "4",
    temp_directory = file.path(cache, "duckdb-tmp"),
    allow_unsigned_extensions = "true", autoinstall_known_extensions = "false",
    autoload_known_extensions = "false"))
connection <- DBI::dbConnect(driver)
DBI::dbExecute(connection, sprintf("LOAD %s",
  as.character(DBI::dbQuoteString(connection,
    normalizePath("build/release/duckhts.duckdb_extension")))))
quote_string <- function(value) as.character(DBI::dbQuoteString(connection, value))
q_data <- utils::read.csv(file.path(cache, "autosome.q.csv"), stringsAsFactors = FALSE)
if (arm == "ancestry") {
  DBI::dbWriteTable(connection, "q_table",
    q_data[c("sample_id", "group_id", "proportion")], overwrite = TRUE, temporary = TRUE)
}
for (chromosome in chromosomes) {
  bcf_path <- file.path(cache, sprintf("chr%d.children.bcf", chromosome))
  af_path <- file.path(cache, sprintf("chr%d.af_reference_sites.parquet", chromosome))
  ref_path <- file.path(cache, sprintf("chr%d.reference_long.parquet", chromosome))
  DBI::dbExecute(connection, sprintf("CREATE VIEW af_source_%d AS SELECT * FROM read_parquet(%s)",
    chromosome, quote_string(af_path)))
  DBI::dbExecute(connection, sprintf("CREATE VIEW ref_long_%d AS SELECT * FROM read_parquet(%s)",
    chromosome, quote_string(ref_path)))
  if (arm != "ancestry") {
    field <- switch(arm, pooled = "INFO_AF", single_AFR = "INFO_AF_AFR",
      single_AMR = "INFO_AF_AMR", single_EUR = "INFO_AF_EUR", single_EAS = "INFO_AF_EAS")
    DBI::dbExecute(connection, sprintf(paste0(
      "CREATE VIEW af_table_%d AS SELECT chrom, pos, ref, alt, ",
      "greatest(1e-3, least(0.999, list_extract(%s, 1)::DOUBLE)) AS af ",
      "FROM af_source_%d"), chromosome, field, chromosome))
  }
}
all_ids <- readLines(file.path(cache, "chr20.children.present.txt"))
if (arm %in% c("ancestry", "pooled")) {
  ids <- all_ids
} else {
  population_data <- read.table(file.path(cache, "pedigree.txt"), header = TRUE,
                                stringsAsFactors = FALSE)
  superpopulation <- sub("^single_", "", arm)
  populations <- switch(superpopulation,
    AFR = c("ACB", "ASW", "ESN", "YRI"),
    AMR = c("CLM", "MXL", "PEL", "PUR"), EUR = "CEU", EAS = "CHS")
  ids <- population_data$SampleID[
    population_data$Population %in% populations & population_data$SampleID %in% all_ids]
}
batch_size <- 8L
batches <- split(ids, ceiling(seq_along(ids) / batch_size))
reads_for_batch <- function(batch) {
  sample_selector <- quote_string(paste(batch, collapse = ","))
  vapply(chromosomes, function(chromosome) {
    bcf_path <- file.path(cache, sprintf("chr%d.children.bcf", chromosome))
    if (arm == "ancestry") {
      sprintf(paste0(
        "SELECT * FROM duckhts_roh_ancestry(%s, 'ref_long_%d', 'q_table', ",
        "af_clamp := 1e-3, gt_error := 30, samples := %s)"),
        quote_string(bcf_path), chromosome, sample_selector)
    } else {
      sprintf(paste0(
        "SELECT * FROM duckhts_roh_af_table(%s, 'af_table_%d', ",
        "gt_error := 30, samples := %s)"),
        quote_string(bcf_path), chromosome, sample_selector)
    }
  }, character(1L))
}
timing <- system.time({
  results <- lapply(seq_along(batches), function(index) {
    cat("arm", arm, "rep", repetition, "batch", index, "of", length(batches),
        "chromosomes", length(chromosomes), "samples", length(batches[[index]]), "\n")
    flush.console()
    sql <- paste(reads_for_batch(batches[[index]]), collapse = " UNION ALL BY NAME ")
    DBI::dbGetQuery(connection, sql)
  })
})
result <- do.call(rbind, results)
output <- file.path(cache, sprintf("autosome.arm_%s_%d.csv", arm, repetition))
utils::write.csv(data.frame(arm = arm, repetition = repetition, result,
                            stringsAsFactors = FALSE), output, row.names = FALSE)
cat("arm", arm, "rep", repetition, "samples", length(ids), "batches",
    length(batches), "chromosomes", length(chromosomes), "segments", nrow(result),
    "elapsed", unname(timing[["elapsed"]]), "output", output, "\n")
DBI::dbDisconnect(connection, shutdown = TRUE)
