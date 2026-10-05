library(tinytest)
library(DBI)

# The fixtures are written by test/scripts/prepare_bam_mismatch_fixtures.R; the
# expected counts are derived by hand in test/sql/bam_mismatch_counts.test.
test_bam_mismatch_counts <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  extdata <- function(name) system.file("extdata", name, package = "Rduckhts")
  bam <- extdata("bam_mismatch.bam")
  reference <- extdata("bam_mismatch.fa")
  mask <- extdata("bam_mismatch.mask.bcf")

  counts <- rduckhts_bam_mismatch_counts(con, bam, reference, mask = mask)
  expect_identical(names(counts),
                   c("mate", "cycle", "base_quality", "reference_base", "read_base", "bases"))
  expect_equal(sum(counts$bases), 99)
  expect_equal(sum(counts$bases[counts$reference_base != counts$read_base]), 4)
  # The read without base qualities gives 12 bases with a missing quality.
  expect_equal(sum(counts$bases[is.na(counts$base_quality)]), 12)
  # The wrapper returns what the SQL function returns.
  direct <- dbGetQuery(con, sprintf(
    "SELECT * FROM duckhts_bam_mismatch_counts(%s, %s, mask := %s)",
    dbQuoteString(con, bam), dbQuoteString(con, reference), dbQuoteString(con, mask)))
  expect_equal(counts, direct)

  # No flank and no mask count every aligned base that is A, C, G or T.
  every <- rduckhts_bam_mismatch_counts(con, bam, reference, indel_flank = 0)
  expect_equal(sum(every$bases), 141)
  expect_equal(sum(every$bases[every$reference_base != every$read_base]), 7)

  region <- rduckhts_bam_mismatch_counts(con, bam, reference, region = "ref1:1-45")
  expect_equal(sum(region$bases), 40)
  expect_equal(sort(unique(region$mate)), c(1, 2))

  unfiltered <- rduckhts_bam_mismatch_counts(con, bam, reference, mask = mask,
                                             exclude_flags = 0, min_mapq = 0)
  expect_equal(sum(unfiltered$bases), 139)

  expect_error(rduckhts_bam_mismatch_counts(con, bam, reference, indel_flank = -1),
               "indel_flank must be one whole number from 0 to 1000")
  expect_error(rduckhts_bam_mismatch_counts(con, bam, reference, min_mapq = 2.5),
               "min_mapq must be one whole number from 0 to 255")
  expect_error(rduckhts_bam_mismatch_counts(con, bam, reference, region = NA_character_),
               "one non-missing character string")
  expect_error(rduckhts_bam_mismatch_counts(con, bam, extdata("fixture_ref.fa")),
               "reference has no contig named ref1")
}

test_bam_mismatch_counts()
