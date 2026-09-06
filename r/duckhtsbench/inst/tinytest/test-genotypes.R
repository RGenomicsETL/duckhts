library(tinytest)

test_genotype_staging <- function() {
  bcftools <- Sys.which("bcftools")
  if (!nzchar(bcftools)) stop("bcftools is required for network-free genotype staging tests")
  directory <- tempfile("genotype-stage-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  source <- file.path(directory, "source.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2", "##contig=<ID=chr22,length=100>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB",
    "chr22\t10\ta\tA\tC,G\t.\tPASS\t.\tGT:PS\t0|2:10\t./1:.",
    "chr22\t20\tb\tA\tT\t.\tPASS\t.\tGT\t0/0\t1",
    "chr22\t90\tout\tT\tC\t.\tPASS\t.\tGT\t0/1\t1/1"), source)
  compressed <- paste0(source, ".gz")
  stopifnot(system2(bcftools, shQuote(c("view", "-Oz", "-o", compressed, source))) == 0L,
            system2(bcftools, shQuote(c("index", "-t", compressed))) == 0L)
  outputs <- c(vcf = file.path(directory, "cohort.vcf.gz"), bcf = file.path(directory, "cohort.bcf"))
  counts <- duckhtsbench:::duckhts_bench_genotype_pair(compressed, "chr22:1-20", outputs, bcftools)
  expect_equal(counts, list(records = 2L, samples = 2L))
  expected <- c("chr22\t10\ta\tA\tC,G\t.\tPASS\t.\tGT:PS\t0|2:10\t./1:.",
                "chr22\t20\tb\tA\tT\t.\tPASS\t.\tGT\t0/0\t1")
  for (path in outputs) {
    expect_equal(system2(bcftools, shQuote(c("view", "-H", path)), stdout = TRUE), expected)
    expect_equal(system2(bcftools, shQuote(c("query", "-l", path)), stdout = TRUE), c("A", "B"))
  }
  plan <- duckhts_bench_stage_plan("genotype-reader")
  expect_equal(plan$id, c("geno_hprc_vcfgz", "geno_hprc_bcf"))
  expect_true(all(grepl("genotypes_removed=false", plan$supplier_identity, fixed = TRUE)))
}
test_genotype_staging()
