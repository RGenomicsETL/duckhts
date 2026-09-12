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
raw_lines <- c(lines, sub("./.", "|0/1", lines[1], fixed=TRUE),
  "chrG\t40\td\tA\tC\t.\tPASS\t.\tPS:GT\t10:/0|1/2",
  "chrG\t50\te\tA\tC\t.\tPASS\t.\tGT\t.",
  "chrG\t60\tf\tA\tC\t.\tPASS\t.\tGT\t",
  "chrG\t70\tg\tA\tC\t.\tPASS\t.\tPS:GT\t10")
raw_expected <- data.frame(record_index=0:8, sample_index=0L,
  raw_gt=c("./.","./.",NA,"1","|0/1","/0|1/2",".","",NA))
raw <- genotype_format_raw_gt(raw_lines)
stopifnot(genotype_format_raw_difference(raw_expected,raw) == 0,
          nrow(genotype_format_raw_gt(character())) == 0L)
for (cut in seq_len(length(raw_lines) - 1L)) {
  batched <- rbind(genotype_format_raw_gt(raw_lines[seq_len(cut)]),
    genotype_format_raw_gt(raw_lines[seq.int(cut + 1L,length(raw_lines))],cut))
  stopifnot(genotype_format_raw_difference(raw_expected,batched) == 0)
}
fails(genotype_format_raw_gt(paste0(lines,"\textra_sample")))
fails(genotype_format_raw_gt(sub("GT:PS", "GT:GT",lines[1],fixed=TRUE)))
raw_mutated <- raw_missing <- raw_dot <- raw_wrong_index <- raw_wrong_sample <- raw
raw_mutated$raw_gt[5] <- "0/1"
raw_missing$raw_gt[5] <- NA_character_
raw_dot$raw_gt[3] <- "."
raw_wrong_index$record_index[2] <- 0L
raw_wrong_sample$sample_index[1] <- 1L
raw_controls <- list(raw_mutated,raw_missing,raw_dot,raw_wrong_index,raw_wrong_sample,
                     raw[-1L,],rbind(raw,raw[1L,]))
stopifnot(all(vapply(raw_controls,function(x)
  genotype_format_raw_difference(raw_expected,x) > 0,logical(1))))
con <- DBI::dbConnect(duckdb::duckdb())
stopifnot(DBI::dbGetQuery(con,"SELECT octet_length(encode('|0/1')) AS n")$n == 4L)
stopifnot(genotype_format_difference(con,
  "SELECT 1 AS id, [10,NULL,5] AS AD", "SELECT 1 AS id, [10,NULL,5] AS AD") == 0)
stopifnot(genotype_format_difference(con,
  "SELECT 1 AS id, [10,NULL,5] AS AD", "SELECT 1 AS id, [10,5] AS AD") == 2)
stopifnot(genotype_format_difference(con,
  "SELECT 1 AS id UNION ALL SELECT 1", "SELECT 1 AS id") == 1)
stopifnot(genotype_format_difference(con, "SELECT 15 AS DP", "SELECT 14 AS DP") == 2)
stopifnot(genotype_format_difference(con,
  "SELECT NULL::INTEGER AS GQ", "SELECT 0 AS GQ") == 2)
input <- tempfile("genotype-format-source-",fileext=".vcf.gz")
stream <- gzfile(input,"wt")
writeLines(c("##fileformat=VCFv4.4",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",raw_lines),stream)
close(stream)
DBI::dbWriteTable(con,"raw_calls",raw)
compare_raw <- function(batch_size) {
  result <- DBI::dbSendQuery(con,"SELECT * FROM raw_calls ORDER BY record_index,sample_index")
  on.exit(DBI::dbClearResult(result))
  genotype_format_raw_compare(input,result,batch_size)
}
raw_counts <- c(records=9,raw_gt_values=7,raw_gt_bytes=18,differences=0)
for (batch_size in c(1L,2L,3L,65536L))
  stopifnot(identical(compare_raw(batch_size),raw_counts))
for (control in raw_controls) {
  DBI::dbWriteTable(con,"raw_calls",control,overwrite=TRUE)
  stopifnot(compare_raw(2L)[["differences"]] > 0)
}
baseline <- tempfile("genotype-format-baseline-",fileext=".md")
default <- expand.grid(selection=c("GT_PS","GT_PS_AD","GT_PS_AD_DP_GQ"),
                       calls_projected=c(TRUE,FALSE),stringsAsFactors=FALSE)
default$elapsed <- seq_len(6L)
default$peak_rss_kib <- 100L
default$raw_gt <- FALSE
current <- rbind(default,transform(default[1L,],raw_gt=TRUE))
writeLines(c(sprintf("Input artifact: test ; bytes: %s ; observed MD5: %s",
  file.info(input)$size,unname(tools::md5sum(input))),
  "| denominator | count |","|---|---|",
  sprintf("| %s | %s |",names(expected),expected),
  "","| selection | calls_projected | elapsed | peak_rss_kib |","|---|---|---|---|",
  sprintf("| %s | %s | %s | %s |",default$selection,default$calls_projected,
    default$elapsed,default$peak_rss_kib)),baseline)
comparison <- genotype_format_comparison(baseline,current,expected,"test",input)
stopifnot(nrow(comparison) == 6L,all(comparison$elapsed_change_percent == 0))
fails(genotype_format_comparison(baseline,current,expected+1,"test",input))
fails(genotype_format_comparison(baseline,current[-1L,],expected,"test",input))
fails(genotype_format_comparison(baseline,current,expected,"other",input))
unlink(baseline)
unlink(input)
DBI::dbDisconnect(con, shutdown=TRUE)
comparison <- genotype_hprc_comparison("benchmarks/benchmark_genotypes.md",
                                      "benchmarks/benchmark_genotypes.md")
stopifnot(nrow(comparison) == 16L, all(comparison$elapsed_change_percent == 0))
cat("Genotype FORMAT benchmark: slot counts, raw source bytes, batch invariance and corruption controls: OK\n")
