# One repetition of one genome-wide ROH arm over the autosomes of the run.
# Usage, from the repository root: Rscript benchmarks/roh_ancestry/run_autosome_arm.R ARM REP
# Each query decodes one chromosome for one batch of children with the shipped
# macro, duckhts_roh_ancestry or duckhts_roh_af_table, and the chromosomes run
# one after another, so the sites held at once follow one chromosome and one
# batch. The batch sizes are those of run_chr20_arm.R, which decodes one
# chromosome per query too: 16 children for the ancestry arm, 64 otherwise.
source("benchmarks/roh_ancestry/genome_common.R")
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
chromosomes <- genome_chromosomes()
driver <- duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE,
  allow_extensions = TRUE, config = list(memory_limit = "12GB", threads = "4",
    temp_directory = file.path(genome_cache, "duckdb-tmp"),
    allow_unsigned_extensions = "true", autoinstall_known_extensions = "false",
    autoload_known_extensions = "false"))
connection <- DBI::dbConnect(driver)
quote_string <- function(value) as.character(DBI::dbQuoteString(connection, value))
DBI::dbExecute(connection, sprintf("LOAD %s", quote_string(
  normalizePath("build/release/duckhts.duckdb_extension", mustWork = TRUE))))
if (arm == "ancestry") {
  q_data <- utils::read.csv(file.path(genome_cache, "autosome.q.csv"), stringsAsFactors = FALSE)
  DBI::dbWriteTable(connection, "q_table",
    q_data[c("sample_id", "group_id", "proportion")], overwrite = TRUE, temporary = TRUE)
}
field <- switch(arm, pooled = "INFO_AF", single_AFR = "INFO_AF_AFR",
  single_AMR = "INFO_AF_AMR", single_EUR = "INFO_AF_EUR", single_EAS = "INFO_AF_EAS",
  ancestry = NA_character_)
for (chromosome in chromosomes) {
  if (arm == "ancestry") {
    DBI::dbExecute(connection, sprintf(
      "CREATE VIEW ref_long_%d AS SELECT * FROM read_parquet(%s)",
      chromosome, quote_string(genome_reference_long(chromosome))))
  } else {
    DBI::dbExecute(connection, sprintf(paste0(
      "CREATE VIEW af_table_%d AS SELECT chrom, pos, ref, alt, ",
      "greatest(1e-3, least(0.999, list_extract(%s, 1)::DOUBLE)) AS af ",
      "FROM read_parquet(%s)"), chromosome, field,
      quote_string(genome_af_sites(chromosome))))
  }
}
all_ids <- readLines(genome_children_list)
if (arm %in% c("ancestry", "pooled")) {
  ids <- all_ids
} else {
  population_data <- read.table(duckhtsbench::duckhts_bench_artifact_path("roh_ancestry_pedigree"),
                                header = TRUE, stringsAsFactors = FALSE)
  superpopulation <- sub("^single_", "", arm)
  populations <- switch(superpopulation,
    AFR = c("ACB", "ASW", "ESN", "YRI"),
    AMR = c("CLM", "MXL", "PEL", "PUR"), EUR = "CEU", EAS = "CHS")
  ids <- population_data$SampleID[
    population_data$Population %in% populations & population_data$SampleID %in% all_ids]
}
batch_size <- if (arm == "ancestry") 16L else 64L
batches <- split(ids, ceiling(seq_along(ids) / batch_size))
decode_sql <- function(chromosome, batch) {
  bcf_path <- quote_string(genome_children_bcf(chromosome))
  sample_selector <- quote_string(paste(batch, collapse = ","))
  if (arm == "ancestry") {
    sprintf(paste0(
      "SELECT * FROM duckhts_roh_ancestry(%s, 'ref_long_%d', 'q_table', ",
      "af_clamp := 1e-3, gt_error := 30, samples := %s)"),
      bcf_path, chromosome, sample_selector)
  } else {
    sprintf(paste0(
      "SELECT * FROM duckhts_roh_af_table(%s, 'af_table_%d', ",
      "gt_error := 30, samples := %s)"),
      bcf_path, chromosome, sample_selector)
  }
}
results <- vector("list", length(chromosomes) * length(batches))
timing <- system.time({
  for (chromosome_index in seq_along(chromosomes)) {
    chromosome <- chromosomes[[chromosome_index]]
    for (batch_index in seq_along(batches)) {
      cat("arm", arm, "rep", repetition, "chromosome", chromosome, "batch", batch_index,
          "of", length(batches), "samples", length(batches[[batch_index]]), "\n")
      flush.console()
      slot <- (chromosome_index - 1L) * length(batches) + batch_index
      results[[slot]] <- DBI::dbGetQuery(connection,
                                         decode_sql(chromosome, batches[[batch_index]]))
    }
  }
})
result <- do.call(rbind, results)
output <- file.path(genome_cache, sprintf("autosome.arm_%s_%d.csv", arm, repetition))
utils::write.csv(data.frame(arm = arm, repetition = repetition, result,
                            stringsAsFactors = FALSE), output, row.names = FALSE)
cat("arm", arm, "rep", repetition, "samples", length(ids), "batches", length(batches),
    "chromosomes", length(chromosomes), "queries", length(results), "segments", nrow(result),
    "elapsed", unname(timing[["elapsed"]]), "output", output, "\n")
DBI::dbDisconnect(connection, shutdown = TRUE)
