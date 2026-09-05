library(tinytest)
library(DBI)

test_bam_full_scan <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  read_groups <- system.file("extdata", "bam_read_groups.sam", package = "Rduckhts")
  expect_true(nzchar(read_groups))
  groups <- dbGetQuery(con, sprintf(paste(
    "SELECT QNAME, READ_GROUP_ID, SAMPLE_ID FROM read_bam(%s,",
    "scan_mode := 'sequential', decompression_threads := 0) ORDER BY QNAME"),
    dbQuoteString(con, read_groups)))
  expect_equal(groups$QNAME, paste0("r", 1:8))
  expect_equal(groups$READ_GROUP_ID, c("one", "one", "missing", "missing", "three", NA, "unknown", "one"))
  expect_equal(groups$SAMPLE_ID, c("sample_one", "sample_one", NA, NA, "sample_three", NA, NA, "sample_one"))
  expected_counts <- c(mixed = 5L, single = 3L, all_unplaced = 2053L, empty = 0L)
  for (threads in c(1L, 4L)) {
    dbExecute(con, sprintf("SET threads=%d", threads))
    for (name in names(expected_counts)) {
      for (index in c("bai", "csi", "legacy.bai", "legacy.csi", "crai")) {
        format <- if (index == "crai") "cram" else "bam"
        path <- system.file("extdata", paste0("bam_scan_", name, ".", format), package = "Rduckhts")
        expect_true(nzchar(path))
        quoted <- dbQuoteString(con, path)
        source <- sprintf("read_bam(%s, index_path := %s, decompression_threads := 0)",
                          quoted, dbQuoteString(con, paste0(path, ".", index)))
        counts <- dbGetQuery(con, sprintf("SELECT count(*) AS n, count(QNAME) AS projected FROM %s", source))
        expect_equal(as.integer(counts$n), expected_counts[[name]])
        expect_equal(as.integer(counts$projected), expected_counts[[name]])
        # Same-file scans agree on every column, including BAM offsets / CRAM NULLs.
        automatic <- dbGetQuery(con, sprintf("SELECT * FROM %s ORDER BY QNAME, POS, FILE_OFFSET", source))
        sequential <- dbGetQuery(con, sprintf(paste(
          "SELECT * FROM read_bam(%s,",
          "scan_mode := 'sequential', decompression_threads := 0) ORDER BY QNAME, POS, FILE_OFFSET"), quoted))
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
