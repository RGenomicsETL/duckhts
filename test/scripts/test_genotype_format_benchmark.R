#!/usr/bin/env Rscript
# Network-free checks for benchmark denominators and the full multiset comparator.
source("benchmarks/genotype_format_run.R")
lines <- c(
  "chrG\t10\ta\tA\tC,G\t.\tPASS\t.\tGT:PS:AD:DP:GQ\t./.:10:10,.,5:15:42",
  "chrG\t10\ta\tA\tC,G\t.\tPASS\t.\tGT:PS:AD:DP:GQ\t./.:10:10,.,5:15:42",
  "chrG\t20\tb\tA\tC\t.\tPASS\t.\tAD:DP\t.:.",
  "chrG\t30\tc\tA\tC\t.\tPASS\t.\tGT\t1")
expected <- c(records=4, calls=4, gt_slots=5, ps_values=2,
              ad_slots=7, ad_values=4, dp_values=2, gq_values=2)
stopifnot(identical(genotype_format_counts(lines), expected),
          identical(genotype_format_counts(character()), expected * 0))
for (cut in 1:3) {
  stopifnot(identical(genotype_format_counts(lines[seq_len(cut)]) +
    genotype_format_counts(lines[seq.int(cut + 1L, length(lines))]), expected))
}
fails <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error=identity), "error"))
fails(genotype_format_counts(paste0(lines, "\textra_sample")))
fails(genotype_format_counts(sub("./.", "|0/1", lines[1], fixed=TRUE)))
con <- DBI::dbConnect(duckdb::duckdb())
stopifnot(genotype_format_difference(con,
  "SELECT 1 AS id, [10,NULL,5] AS AD", "SELECT 1 AS id, [10,NULL,5] AS AD") == 0)
stopifnot(genotype_format_difference(con,
  "SELECT 1 AS id, [10,NULL,5] AS AD", "SELECT 1 AS id, [10,5] AS AD") == 2)
stopifnot(genotype_format_difference(con,
  "SELECT 1 AS id UNION ALL SELECT 1", "SELECT 1 AS id") == 1)
stopifnot(genotype_format_difference(con, "SELECT 15 AS DP", "SELECT 14 AS DP") == 2)
stopifnot(genotype_format_difference(con,
  "SELECT NULL::INTEGER AS GQ", "SELECT 0 AS GQ") == 2)
DBI::dbDisconnect(con, shutdown=TRUE)
comparison <- genotype_hprc_comparison("benchmarks/benchmark_genotypes.md",
                                      "benchmarks/benchmark_genotypes.md")
stopifnot(nrow(comparison) == 16L, all(comparison$elapsed_change_percent == 0))
cat("Genotype FORMAT benchmark: independent slot counts, batch invariance and corruption controls: OK\n")
