library(tinytest)
library(DBI)

# A small coordinate-sorted SAM on the contigs of bam_mismatch.fa, where ref1
# repeats TGGCAA from position 1. test/sql/bam_mate_overlap_counts.test derives
# the counts of a larger fixture by hand.
#   ovl     40/20  10M at 11 and reverse 10M at 15: 15-20 shared, read 2 shows C for A at 18
#   lowmate 15/15  read 2 has MAPQ 5: counted only without the MAPQ floor
#   noqual  -/40   read 1 stores no qualities: 145-150 shared
#   eq      25/25  read 1 has '=' at 166: 165-170 shared, 166 only with the reference
mate_overlap_sam <- function(path) {
  lines <- c(
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:ref1\tLN:260",
    "@SQ\tSN:ref2\tLN:80",
    "ovl\t99\tref1\t11\t60\t10M\t=\t15\t14\tAATGGCAATG\tIIIIIIIIII",
    "ovl\t147\tref1\t15\t60\t10M\t=\t11\t-14\tGCACTGGCAA\t5555555555",
    "lowmate\t99\tref1\t41\t60\t10M\t=\t45\t14\tAATGGCAATG\t0000000000",
    "lowmate\t147\tref1\t45\t5\t10M\t=\t41\t-14\tGCAATGGCAA\t0000000000",
    "noqual\t99\tref1\t141\t60\t10M\t=\t145\t14\tGCAATGGCAA\t*",
    "noqual\t147\tref1\t145\t60\t10M\t=\t141\t-14\tTGGCAATGGC\tIIIIIIIIII",
    "eq\t99\tref1\t161\t60\t10M\t=\t165\t14\tAATGG=AATG\t::::::::::",
    "eq\t147\tref1\t165\t60\t10M\t=\t161\t-14\tGCAATGGCAA\t::::::::::"
  )
  writeLines(lines, path)
  path
}

test_bam_mate_overlap_counts <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  extdata <- function(name) system.file("extdata", name, package = "Rduckhts")
  sam <- mate_overlap_sam(tempfile(fileext = ".sam"))
  on.exit(unlink(sam), add = TRUE)

  counts <- rduckhts_bam_mate_overlap_counts(con, sam)
  expect_identical(names(counts),
                   c("read1_quality", "read2_quality", "read1_base", "read2_base", "positions"))
  # ovl 6, noqual 6 and eq 5: eq leaves out its '=' base without a reference.
  expect_equal(sum(counts$positions), 17)
  disagree <- counts[counts$read1_base != counts$read2_base, ]
  expect_equal(nrow(disagree), 1L)
  expect_equal(disagree$read1_quality, 40)
  expect_equal(disagree$read2_quality, 20)
  expect_identical(c(disagree$read1_base, disagree$read2_base), c("A", "C"))
  expect_equal(sum(counts$positions[is.na(counts$read1_quality)]), 6)
  # The wrapper returns what the SQL function returns.
  direct <- dbGetQuery(con, sprintf(
    "SELECT * FROM duckhts_bam_mate_overlap_counts(%s)", dbQuoteString(con, sam)))
  expect_equal(counts, direct)

  # The reference resolves the '=' base; the MAPQ floor drops lowmate.
  with_reference <- rduckhts_bam_mate_overlap_counts(con, sam, reference = extdata("bam_mismatch.fa"))
  expect_equal(sum(with_reference$positions), 18)
  unfiltered <- rduckhts_bam_mate_overlap_counts(con, sam, min_mapq = 0)
  expect_equal(sum(unfiltered$positions), 23)
  expect_equal(sum(rduckhts_bam_mate_overlap_counts(con, sam, max_pending = 1)$positions), 17)

  # A bundled paired-end BAM; the SQL test checks it against an independent count.
  bam <- extdata("empty-tids.bam")
  whole <- rduckhts_bam_mate_overlap_counts(con, bam, min_mapq = 0)
  expect_equal(sum(whole$positions), 2316)

  expect_error(rduckhts_bam_mate_overlap_counts(con, sam, max_pending = 0),
               "max_pending must be one whole number from 1 to 100000000")
  expect_error(rduckhts_bam_mate_overlap_counts(con, sam, min_mapq = 2.5),
               "min_mapq must be one whole number from 0 to 255")
  expect_error(rduckhts_bam_mate_overlap_counts(con, sam, exclude_flags = 70000),
               "exclude_flags must be one whole number from 0 to 65535")
  expect_error(rduckhts_bam_mate_overlap_counts(con, sam, region = NA_character_),
               "one non-missing character string")
  expect_error(rduckhts_bam_mate_overlap_counts(con, sam, region = "ref1:1-100"),
               "region needs an index")
  expect_error(rduckhts_bam_mate_overlap_counts(con, extdata("bam_mismatch.cram")),
               "needs reference")
}

test_bam_mate_overlap_counts()
