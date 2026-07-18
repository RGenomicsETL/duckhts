library(tinytest)
library(DBI)

test_parallel_empty_contigs <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  expect_silent(rduckhts_load(con))
  expect_silent(dbExecute(con, "PRAGMA threads=4"))

  bam_path <- system.file("extdata", "parallel_empty_contigs.bam", package = "Rduckhts")
  bam_index_path <- system.file("extdata", "parallel_empty_contigs.bam.bai", package = "Rduckhts")
  vcf_path <- system.file("extdata", "parallel_empty_contigs.vcf.gz", package = "Rduckhts")
  vcf_index_path <- system.file("extdata", "parallel_empty_contigs.vcf.gz.tbi", package = "Rduckhts")
  bcf_path <- system.file("extdata", "parallel_empty_contigs.bcf", package = "Rduckhts")
  bcf_index_path <- system.file("extdata", "parallel_empty_contigs.bcf.csi", package = "Rduckhts")

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
  quoted_bcf <- DBI::dbQuoteString(con, bcf_path)
  for (threads in c(1L, 4L)) {
    target <- sprintf("parallel_empty_contigs_appender_%d", threads)
    result <- DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT rows_written FROM read_bcf_appender(",
          "%s, '%s', region := '%s', tidy_format := true,",
          "overwrite := true, include_file_offset := true, region_threads := %d)"
        ),
        quoted_bcf, target, regions, threads
      )
    )
    expect_equal(result$rows_written[[1]], 2)
  }

  appender_counts <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT count(*) AS n, count(DISTINCT FILE_OFFSET) AS offsets",
      "FROM parallel_empty_contigs_appender_4"
    )
  )
  expect_equal(appender_counts$n[[1]], 2)
  expect_equal(appender_counts$offsets[[1]], 2)

  appender_mismatch <- DBI::dbGetQuery(
    con,
    paste0(
      "SELECT count(*) AS n FROM (",
      "(SELECT * FROM parallel_empty_contigs_appender_1 ",
      " EXCEPT ALL SELECT * FROM parallel_empty_contigs_appender_4) UNION ALL ",
      "(SELECT * FROM parallel_empty_contigs_appender_4 ",
      " EXCEPT ALL SELECT * FROM parallel_empty_contigs_appender_1))"
    )
  )
  expect_equal(appender_mismatch$n[[1]], 0)
}

test_parallel_empty_contigs()
