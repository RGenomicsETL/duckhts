library(tinytest)
library(DBI)

test_parallel_empty_contigs <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  expect_silent(dbExecute(con, "PRAGMA threads=4"))

  bam_path <- system.file(
    "extdata",
    "parallel_empty_contigs.bam",
    package = "Rduckhts"
  )
  bam_index_path <- system.file(
    "extdata",
    "parallel_empty_contigs.bam.bai",
    package = "Rduckhts"
  )
  vcf_path <- system.file(
    "extdata",
    "parallel_empty_contigs.vcf.gz",
    package = "Rduckhts"
  )
  vcf_index_path <- system.file(
    "extdata",
    "parallel_empty_contigs.vcf.gz.tbi",
    package = "Rduckhts"
  )
  bcf_path <- system.file(
    "extdata",
    "parallel_empty_contigs.bcf",
    package = "Rduckhts"
  )
  bcf_index_path <- system.file(
    "extdata",
    "parallel_empty_contigs.bcf.csi",
    package = "Rduckhts"
  )

  expect_true(file.exists(bam_path))
  expect_true(file.exists(bam_index_path))
  expect_true(file.exists(vcf_path))
  expect_true(file.exists(vcf_index_path))
  expect_true(file.exists(bcf_path))
  expect_true(file.exists(bcf_index_path))

  bam_counts <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT count(*) AS n, min(RNAME) AS chrom",
        "FROM read_bam('%s', index_path := '%s')"
      ),
      bam_path,
      bam_index_path
    )
  )
  expect_equal(bam_counts$n[[1]], 2)
  expect_equal(bam_counts$chrom[[1]], "chr1")

  bcf_counts <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT count(*) AS n, min(CHROM) AS chrom",
        "FROM read_bcf('%s', index_path := '%s')"
      ),
      vcf_path,
      vcf_index_path
    )
  )
  expect_equal(bcf_counts$n[[1]], 2)
  expect_equal(bcf_counts$chrom[[1]], "chr1")

  regions <- paste(c(paste0("empty", 1:8), "chr1"), collapse = ",")
  for (path in c(bcf_path, vcf_path)) {
    quoted <- dbQuoteString(con, path)
    selected <- dbGetQuery(con, sprintf(paste(
      "SELECT * FROM read_bcf(%s, region := %s, tidy_format := true)",
      "ORDER BY CHROM,POS,SAMPLE_ID"), quoted, dbQuoteString(con, regions)))
    sequential <- dbGetQuery(con, sprintf(paste(
      "SELECT * FROM read_bcf(%s, scan_mode := 'sequential', tidy_format := true)",
      "ORDER BY CHROM,POS,SAMPLE_ID"), quoted))
    expect_equal(nrow(selected), 2L)
    expect_equal(selected, sequential)
  }
}

test_parallel_empty_contigs()
