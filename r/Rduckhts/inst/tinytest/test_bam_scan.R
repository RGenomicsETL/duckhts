library(tinytest)
library(DBI)

test_bam_full_scan <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  expected_counts <- c(mixed = 5L, single = 3L, all_unplaced = 2053L, empty = 0L)
  for (threads in c(1L, 4L)) {
    dbExecute(con, sprintf("SET threads=%d", threads))
    for (name in names(expected_counts)) {
      for (index in c("bai", "csi", "crai")) {
        format <- if (index == "crai") "cram" else "bam"
        path <- system.file("extdata", paste0("bam_scan_", name, ".", format), package = "Rduckhts")
        expect_true(nzchar(path))
        quoted <- dbQuoteString(con, path)
        source <- sprintf("read_bam(%s, index_path := %s, decompression_threads := 0)",
                          quoted, dbQuoteString(con, paste0(path, ".", index)))
        counts <- dbGetQuery(con, sprintf("SELECT count(*) AS n, count(QNAME) AS projected FROM %s", source))
        expect_equal(as.integer(counts$n), expected_counts[[name]])
        expect_equal(as.integer(counts$projected), expected_counts[[name]])
        # FILE_OFFSET is a transport diagnostic, not a SAM field or CRAM locator.
        automatic <- dbGetQuery(con, sprintf("SELECT * EXCLUDE(FILE_OFFSET) FROM %s ORDER BY QNAME, POS", source))
        sequential <- dbGetQuery(con, sprintf(paste(
          "SELECT * EXCLUDE(FILE_OFFSET) FROM read_bam(%s,",
          "scan_mode := 'sequential', decompression_threads := 0) ORDER BY QNAME, POS"), quoted))
        expect_equal(automatic, sequential)
        if (name == "mixed") {
          expect_equal(automatic$QNAME, c("mapped1", "mapped2", "placed_unmapped", "unplaced", "unplaced"))
          region <- dbGetQuery(con, sprintf(paste(
            "SELECT QNAME FROM read_bam(%s, region := 'chr1',",
            "decompression_threads := 0) ORDER BY QNAME"), quoted))
          expect_equal(region$QNAME, c("mapped1", "placed_unmapped"))
        }
      }
    }
  }
  malformed <- system.file("extdata", "bam_scan_malformed.sam", package = "Rduckhts")
  expect_error(dbGetQuery(con, sprintf(paste(
    "SELECT QNAME FROM read_bam(%s, scan_mode := 'sequential',",
    "decompression_threads := 0)"), dbQuoteString(con, malformed))),
    pattern = "read_bam: failed to read SAM/BAM/CRAM record")
  expect_equal(dbGetQuery(con, "SELECT 4242 AS n")$n, 4242L)
}

test_bam_full_scan()
