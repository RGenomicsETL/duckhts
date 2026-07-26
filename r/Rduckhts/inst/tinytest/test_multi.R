# Multi-file reading wrapper tests
library(tinytest)
library(DBI)

test_multi_read <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

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
  # rduckhts_bam_multi: combine two BAM files
  # -------------------------------------------------------
  rduckhts_bam_multi(con, "bam_multi", c(bam_path, bam2_path))
  res <- dbGetQuery(con, "SELECT * FROM bam_multi")
  expect_true(is.data.frame(res))
  expect_true("filename" %in% names(res))
  expect_true(nrow(res) > 0)
  expect_equal(length(unique(res$filename)), 2L)
  expect_true(all(unique(res$filename) %in% c(bam_path, bam2_path)))

  # overwrite guard
  expect_error(rduckhts_bam_multi(con, "bam_multi", bam_path),
               "already exists")

  # overwrite = TRUE replaces
  rduckhts_bam_multi(con, "bam_multi", bam_path, overwrite = TRUE)
  res_single <- dbGetQuery(con, "SELECT DISTINCT filename FROM bam_multi")
  expect_equal(nrow(res_single), 1L)

  # -------------------------------------------------------
  # rduckhts_bcf_multi: combine two VCF/BCF files
  # -------------------------------------------------------
  rduckhts_bcf_multi(con, "bcf_multi", c(bcf_path, bcf2_path))
  res <- dbGetQuery(con, "SELECT * FROM bcf_multi")
  expect_true("filename" %in% names(res))
  expect_true(nrow(res) > 0)
  expect_equal(length(unique(res$filename)), 2L)

  # -------------------------------------------------------
  # rduckhts_fastq_multi: combine two FASTQ files
  # -------------------------------------------------------
  rduckhts_fastq_multi(con, "fq_multi", c(fq_r1, fq_r2))
  res <- dbGetQuery(con, "SELECT * FROM fq_multi")
  expect_true("filename" %in% names(res))
  expect_true(nrow(res) > 0)
  expect_equal(length(unique(res$filename)), 2L)

  # -------------------------------------------------------
  # rduckhts_fasta_multi: single file (only one FASTA in extdata)
  # -------------------------------------------------------
  rduckhts_fasta_multi(con, "fa_multi", fasta_path)
  res <- dbGetQuery(con, "SELECT * FROM fa_multi")
  expect_true("filename" %in% names(res))
  expect_true(nrow(res) > 0)
  expect_equal(unique(res$filename), fasta_path)

  # -------------------------------------------------------
  # rduckhts_bed_multi: single file
  # -------------------------------------------------------
  rduckhts_bed_multi(con, "bed_multi", bed_path)
  res <- dbGetQuery(con, "SELECT * FROM bed_multi")
  expect_true("filename" %in% names(res))
  expect_true(nrow(res) > 0)
  expect_equal(unique(res$filename), bed_path)

  # -------------------------------------------------------
  # rduckhts_gff_multi: single file
  # -------------------------------------------------------
  rduckhts_gff_multi(con, "gff_multi", gff_path)
  res <- dbGetQuery(con, "SELECT * FROM gff_multi")
  expect_true("filename" %in% names(res))
  expect_true(nrow(res) > 0)

  # -------------------------------------------------------
  # .params data.frame: per-file parameter overrides
  # -------------------------------------------------------
  params_df <- data.frame(file = c(fq_r1, fq_r2), stringsAsFactors = FALSE)
  rduckhts_fastq_multi(con, "fq_params", files = character(0),
                       .params = params_df)
  res <- dbGetQuery(con, "SELECT * FROM fq_params")
  expect_true(nrow(res) > 0)
  expect_equal(length(unique(res$filename)), 2L)

  # -------------------------------------------------------
  # .params validation: missing file column
  # -------------------------------------------------------
  expect_error(
    rduckhts_bam_multi(con, "t_bad", bam_path,
                       .params = data.frame(x = 1)),
    "file"
  )

  # .params validation: not a data.frame
  expect_error(
    rduckhts_bam_multi(con, "t_bad2", bam_path, .params = "bad"),
    "data.frame"
  )

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

  # Execute the generated query via SET VARIABLE pattern
  DBI::dbExecute(con, sprintf(
    "SET VARIABLE q = hts_union_query('read_bam', '%s')",
    gsub("'", "''", bam_path)
  ))
  generated <- DBI::dbGetQuery(con, "SELECT * FROM query(getvariable('q'))")
  expect_true(is.data.frame(generated))
  expect_true("filename" %in% names(generated))
  expect_true(nrow(generated) > 0)

  message("Multi-file reading tests passed!")
}

test_multi_read()
