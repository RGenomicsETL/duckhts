# Multi-file reading wrapper tests
library(tinytest)
library(DBI)

test_multi_read <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  expect_silent(rduckhts_load(con))

  # Locate bundled fixture files
  bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
  bam2_path <- system.file("extdata", "parallel_empty_contigs.bam", package = "Rduckhts")
  bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
  bcf2_path <- system.file("extdata", "formatcols.vcf.gz", package = "Rduckhts")
  fq_r1 <- system.file("extdata", "r1.fq", package = "Rduckhts")
  fq_r2 <- system.file("extdata", "r2.fq", package = "Rduckhts")
  fasta_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
  bed_path <- system.file("extdata", "targets.bed", package = "Rduckhts")
  gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")

  expect_true(file.exists(bam_path))
  expect_true(file.exists(bam2_path))
  expect_true(file.exists(bcf_path))
  expect_true(file.exists(bcf2_path))
  expect_true(file.exists(fq_r1))
  expect_true(file.exists(fq_r2))
  expect_true(file.exists(fasta_path))
  expect_true(file.exists(bed_path))
  expect_true(file.exists(gff_path))

  # -------------------------------------------------------
  # read_bam_multi: combine two BAM files
  # -------------------------------------------------------
  res <- read_bam_multi(con, c(bam_path, bam2_path))
  expect_true(is.data.frame(res))
  expect_true("filename" %in% names(res))
  expect_true(nrow(res) > 0)
  # Should have rows from both files
  expect_equal(length(unique(res$filename)), 2L)
  # Each filename should be one of the inputs
  expect_true(all(unique(res$filename) %in% c(bam_path, bam2_path)))

  # read_bam_multi: single file (verify filename column still present)
  res_single <- read_bam_multi(con, bam_path)
  expect_true(is.data.frame(res_single))
  expect_equal(unique(res_single$filename), bam_path)

  # -------------------------------------------------------
  # read_bcf_multi: combine two VCF/BCF files
  # -------------------------------------------------------
  res <- read_bcf_multi(con, c(bcf_path, bcf2_path))
  expect_true(is.data.frame(res))
  expect_true("filename" %in% names(res))
  expect_true(nrow(res) > 0)
  expect_equal(length(unique(res$filename)), 2L)

  # -------------------------------------------------------
  # read_fastq_multi: combine two FASTQ files
  # -------------------------------------------------------
  res <- read_fastq_multi(con, c(fq_r1, fq_r2))
  expect_true(is.data.frame(res))
  expect_true("filename" %in% names(res))
  expect_true(nrow(res) > 0)
  expect_equal(length(unique(res$filename)), 2L)

  # -------------------------------------------------------
  # read_fasta_multi: single file (only one FASTA in extdata)
  # -------------------------------------------------------
  res <- read_fasta_multi(con, fasta_path)
  expect_true(is.data.frame(res))
  expect_true("filename" %in% names(res))
  expect_true(nrow(res) > 0)
  expect_equal(unique(res$filename), fasta_path)

  # -------------------------------------------------------
  # read_bed_multi: single file
  # -------------------------------------------------------
  res <- read_bed_multi(con, bed_path)
  expect_true(is.data.frame(res))
  expect_true("filename" %in% names(res))
  expect_true(nrow(res) > 0)
  expect_equal(unique(res$filename), bed_path)

  # -------------------------------------------------------
  # read_gff_multi: single file
  # -------------------------------------------------------
  res <- read_gff_multi(con, gff_path)
  expect_true(is.data.frame(res))
  expect_true("filename" %in% names(res))
  expect_true(nrow(res) > 0)

  # -------------------------------------------------------
  # .params data.frame: per-file parameter overrides
  # -------------------------------------------------------
  # Use .params to read two FASTQ files separately (no uniform mate_path)
  params_df <- data.frame(file = c(fq_r1, fq_r2), stringsAsFactors = FALSE)
  res <- read_fastq_multi(con, files = character(0), .params = params_df)
  expect_true(is.data.frame(res))
  expect_true(nrow(res) > 0)
  expect_equal(length(unique(res$filename)), 2L)

  # -------------------------------------------------------
  # .params validation: missing file column
  # -------------------------------------------------------
  expect_error(read_bam_multi(con, bam_path, .params = data.frame(x = 1)),
               "file")

  # .params validation: not a data.frame
  expect_error(read_bam_multi(con, bam_path, .params = "bad"),
               "data.frame")

  # -------------------------------------------------------
  # hts_union_query SQL macro: verify it generates valid SQL
  # -------------------------------------------------------
  query_sql <- DBI::dbGetQuery(con, sprintf(
    "SELECT hts_union_query('read_bam', '%s') AS q",
    gsub("'", "''", bam_path)
  ))
  expect_true(is.data.frame(query_sql))
  expect_true(nrow(query_sql) == 1L)
  expect_true(grepl("read_bam", query_sql$q[1]))
  expect_true(grepl("filename", query_sql$q[1]))

  # Execute the generated query
  generated <- DBI::dbGetQuery(con, query_sql$q[1])
  expect_true(is.data.frame(generated))
  expect_true("filename" %in% names(generated))
  expect_true(nrow(generated) > 0)

  message("Multi-file reading tests passed!")
}

test_multi_read()
