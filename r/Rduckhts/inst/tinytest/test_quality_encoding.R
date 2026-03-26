# Quality encoding tests for Rduckhts
library(tinytest)
library(DBI)

test_quality_encoding <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  expect_silent(rduckhts_load(con))

  bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
  fastq_r1 <- system.file("extdata", "r1.fq", package = "Rduckhts")
  legacy_fastq <- system.file("extdata", "legacy_phred64.fq", package = "Rduckhts")

  expect_true(file.exists(bam_path))
  expect_true(file.exists(fastq_r1))
  expect_true(file.exists(legacy_fastq))

  # BAM QUAL can be materialized as raw Phred integers
  rduckhts_bam(con, "bam_qphred", bam_path,
    quality_representation = "phred",
    overwrite = TRUE)
  bam_q_type <- DBI::dbGetQuery(con,
    "SELECT typeof(QUAL) AS t FROM bam_qphred LIMIT 1")
  expect_equal(bam_q_type$t[1], "UTINYINT[]")

  # FASTQ QUALITY can be materialized as raw Phred integers
  rduckhts_fastq(con, "fq_qphred", fastq_r1,
    quality_representation = "phred",
    overwrite = TRUE)
  fq_qphred_type <- DBI::dbGetQuery(con,
    "SELECT typeof(QUALITY) AS t FROM fq_qphred LIMIT 1")
  expect_equal(fq_qphred_type$t[1], "UTINYINT[]")

  # Detector exposes the legacy fixture as phred64-compatible
  legacy_detect <- rduckhts_detect_quality_encoding(con, legacy_fastq)
  expect_equal(legacy_detect$format[1], "fastq")
  expect_equal(legacy_detect$guessed_encoding[1], "phred64")
  expect_equal(legacy_detect$observed_ascii_min[1], 104)
  expect_equal(legacy_detect$observed_ascii_max[1], 104)
  expect_equal(legacy_detect$is_ambiguous[1], TRUE)
  expect_true(grepl("phred64", legacy_detect$compatible_encodings[1], fixed = TRUE))

  # Default FASTQ decoding assumes modern Phred+33 input
  rduckhts_fastq(con, "legacy_fq", legacy_fastq, overwrite = TRUE)
  legacy_q <- DBI::dbGetQuery(con, "SELECT QUALITY FROM legacy_fq LIMIT 1")
  expect_equal(legacy_q$QUALITY[1], "hhhh")

  # Explicit auto-detection normalizes legacy text qualities to canonical Phred+33 text
  rduckhts_fastq(con, "legacy_fq_auto", legacy_fastq,
    input_quality_encoding = "auto",
    overwrite = TRUE)
  legacy_q_auto <- DBI::dbGetQuery(con, "SELECT QUALITY FROM legacy_fq_auto LIMIT 1")
  expect_equal(legacy_q_auto$QUALITY[1], "IIII")

  # Explicit legacy decoding + numeric representation returns raw Phred scores
  rduckhts_fastq(con, "legacy_fq_phred", legacy_fastq,
    quality_representation = "phred",
    input_quality_encoding = "phred64",
    overwrite = TRUE)
  legacy_q_phred <- DBI::dbGetQuery(con,
    "SELECT QUALITY::VARCHAR AS q FROM legacy_fq_phred LIMIT 1")
  expect_equal(legacy_q_phred$q[1], "[40, 40, 40, 40]")
}

test_quality_encoding()

message("Quality encoding tests passed!")
