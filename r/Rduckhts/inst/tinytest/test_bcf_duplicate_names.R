# VCF keys and sample IDs may differ only in case; DuckDB column names ignore case.
library(tinytest)
library(DBI)

con <- rduckhts_connect()
vcf <- tempfile(fileext = ".vcf")
writeLines(c(
  "##fileformat=VCFv4.3",
  "##contig=<ID=chr1,length=100>",
  '##INFO=<ID=DP,Number=1,Type=Integer,Description="upper">',
  '##INFO=<ID=dp,Number=1,Type=Integer,Description="lower">',
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "s1", "S1", sep = "\t"),
  paste("chr1", "5", ".", "A", "G", ".", "PASS", "DP=1;dp=2", "GT", "0/1", "1/1", sep = "\t")
), vcf)

expect_error(rduckhts_bcf(con, "dup_default", vcf), "INFO_DP and INFO_dp name the same column")
expect_error(rduckhts_bcf(con, "dup_bad", vcf, duplicate_names = "rename"), "should be one of")

rduckhts_bcf(con, "dup_suffix", vcf, duplicate_names = "suffix")
row <- dbGetQuery(con, "SELECT INFO_DP, INFO_dp__2, FORMAT_GT_s1, FORMAT_GT_S1__2 FROM dup_suffix")
expect_equal(row$INFO_DP, 1L)
expect_equal(row$INFO_dp__2, 2L)
expect_equal(row$FORMAT_GT_s1, "0/1")
expect_equal(row$FORMAT_GT_S1__2, "1/1")

unlink(vcf)
dbDisconnect(con, shutdown = TRUE)
