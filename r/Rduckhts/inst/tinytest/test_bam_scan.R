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

# The planner's row estimate of a whole-file indexed scan, read from the plan.
test_bam_row_estimate <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  estimate <- function(path) {
    plan <- dbGetQuery(con, sprintf("EXPLAIN (FORMAT json) SELECT * FROM read_bam(%s)",
                                    dbQuoteString(con, path)))[[2L]]
    found <- regmatches(plan, regexpr('"Estimated Cardinality": "[0-9]+"', plan))
    as.numeric(gsub("[^0-9]", "", found))
  }
  bundled <- function(name) {
    path <- system.file("extdata", name, package = "Rduckhts")
    expect_true(nzchar(path))
    path
  }
  # Every reference of this BAM is empty; the estimate is its unplaced reads.
  expect_equal(estimate(bundled("bam_scan_all_unplaced.bam")), 2053)
  # 186 of 199 references have no reads and no index statistics.
  expect_equal(estimate(bundled("empty-tids.bam")), 12495)
  # A CRAM index has no statistics: the count is right and no estimate is 5.
  cram <- bundled("bam_scan_mixed.cram")
  expect_equal(dbGetQuery(con, sprintf("SELECT count(*) AS n FROM read_bam(%s)",
                                       dbQuoteString(con, cram)))$n, 5)
  expect_false(identical(estimate(cram), 5))
  expect_equal(nrow(rduckhts_hts_index(con, cram)), 0L)
}

# An index that has alignments of a reference and no statistics for it cannot
# locate the reads without coordinates. A full scan with it must still return
# each record once, at any thread count.
test_bam_index_without_statistics <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  bam <- system.file("extdata", "bam_scan_mixed.bam", package = "Rduckhts")
  expect_true(nzchar(bam))
  expected <- dbGetQuery(con, sprintf(paste(
    "SELECT QNAME, FLAG, RNAME, POS FROM read_bam(%s, scan_mode := 'sequential',",
    "decompression_threads := 0) ORDER BY ALL"), dbQuoteString(con, bam)))
  expect_equal(nrow(expected), 5L)
  for (index in c("nostats.bai", "nostats.legacy.bai")) {
    index_path <- paste0(bam, ".", index)
    expect_true(file.exists(index_path))
    for (threads in c(1L, 2L, 4L, 8L)) {
      dbExecute(con, sprintf("SET threads=%d", threads))
      scanned <- dbGetQuery(con, sprintf(paste(
        "SELECT QNAME, FLAG, RNAME, POS FROM read_bam(%s, index_path := %s,",
        "decompression_threads := 0) ORDER BY ALL"),
        dbQuoteString(con, bam), dbQuoteString(con, index_path)))
      expect_equal(scanned, expected)
    }
    in_region <- dbGetQuery(con, sprintf(
      "SELECT count(*) AS n FROM read_bam(%s, index_path := %s, region := 'chr2')",
      dbQuoteString(con, bam), dbQuoteString(con, index_path)))
    expect_equal(as.integer(in_region$n), 1L)
  }
}

test_bam_full_scan()
test_bam_row_estimate()
test_bam_index_without_statistics()
