# Tests for FILE_OFFSET column in read_bam()
#
# FILE_OFFSET is a UBIGINT column that exposes the BGZF virtual file offset
# after each BGZF BAM record, not its start. Other transports return NA.
# It enables ORDER BY FILE_OFFSET within one unchanged BAM in SQL window functions
# to reproduce exact BAM file order, required for streaming dedup algorithms
# such as WisecondorX's larp/larp2 state machine.
library(tinytest)
library(DBI)

test_file_offset <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
  expect_true(file.exists(bam_path))

  # ---- Schema: FILE_OFFSET present as UBIGINT, immediately after SAMPLE_ID ----
  schema <- dbGetQuery(
    con,
    sprintf("DESCRIBE SELECT * FROM read_bam('%s')", bam_path)
  )
  col_names <- schema$column_name

  expect_true(
    "FILE_OFFSET" %in% col_names,
    info = "FILE_OFFSET column present in read_bam() schema"
  )
  fo_type <- schema$column_type[col_names == "FILE_OFFSET"]
  expect_identical(fo_type, "UBIGINT", info = "FILE_OFFSET is UBIGINT")

  sample_id_pos <- which(col_names == "SAMPLE_ID")
  file_offset_pos <- which(col_names == "FILE_OFFSET")
  expect_true(
    file_offset_pos == sample_id_pos + 1L,
    info = "FILE_OFFSET immediately follows SAMPLE_ID"
  )

  # ---- Values: positive, unique, monotonically non-decreasing ----
  offsets <- dbGetQuery(
    con,
    sprintf("SELECT FILE_OFFSET FROM read_bam('%s')", bam_path)
  )$FILE_OFFSET

  expect_true(all(!is.na(offsets)), info = "FILE_OFFSET has no NAs")
  expect_true(all(offsets > 0L), info = "FILE_OFFSET values are positive")
  expect_identical(
    length(unique(offsets)),
    length(offsets),
    info = "FILE_OFFSET is unique per record (no two records share an offset)"
  )
  expect_true(
    all(diff(sort(offsets)) > 0L),
    info = "FILE_OFFSET values are strictly increasing when explicitly sorted"
  )

  # ---- Window function: LAG(FILE_OFFSET) usable without errors ----
  # Verify the column works as an ORDER BY key in a window function; this is
  # the foundation of the larp/larp2 SQL dedup used by WisecondorX.
  #
  # NOTE: window functions cannot be nested inside aggregate FILTER clauses in
  # DuckDB ("Binder Error: aggregate function calls cannot contain window
  # function calls"). Compute the window in a CTE first, then aggregate.
  lag_q <- dbGetQuery(
    con,
    sprintf(
      "
    WITH windowed AS (
      SELECT
        FILE_OFFSET,
        LAG(FILE_OFFSET) OVER (ORDER BY FILE_OFFSET) AS prev_offset
      FROM read_bam('%s')
    )
    SELECT
      COUNT(*)::INTEGER                                             AS total,
      COUNT(DISTINCT FILE_OFFSET)::INTEGER                         AS distinct_offsets,
      COUNT(*) FILTER (WHERE prev_offset >= FILE_OFFSET)::INTEGER  AS order_violations
    FROM windowed
  ",
      bam_path
    )
  )

  expect_identical(
    lag_q$total[1],
    lag_q$distinct_offsets[1],
    info = "no duplicate FILE_OFFSET values"
  )
  expect_identical(
    lag_q$order_violations[1],
    0L,
    info = "no order violations: FILE_OFFSET strictly increases"
  )

  extdata <- function(name) {
    path <- system.file("extdata", name, package = "Rduckhts")
    expect_true(nzchar(path))
    path
  }
  oracle <- read.delim(extdata("bam_record_ends.tsv"), stringsAsFactors = FALSE)
  unsupported <- c("bam_scan_mixed.cram", "bam_offset.sam", "bam_offset.sam.gz",
                   "bam_offset.sam.bgz", "bam_offset_uncompressed.bam", "bam_offset_gzip.bam")
  biological <- dbGetQuery(con, sprintf(paste(
    "SELECT * EXCLUDE(FILE_OFFSET) FROM read_bam(%s)",
    "ORDER BY QNAME, POS"), dbQuoteString(con, extdata("bam_scan_mixed.bam"))))
  for (workers in c(1L, 4L)) {
    dbExecute(con, sprintf("SET threads=%d", workers))
    for (decompression in c(0L, 2L)) {
      for (mode in c("auto", "sequential")) {
        for (name in c("bam_scan_mixed.bam", "nanopore.bam", unsupported)) {
          sql <- sprintf(paste(
            "SELECT * FROM read_bam(%s, scan_mode := '%s', decompression_threads := %d)",
            "ORDER BY QNAME, POS, FILE_OFFSET"), dbQuoteString(con, extdata(name)), mode, decompression)
          rows <- dbGetQuery(con, sql)
          expect_equal(dbGetQuery(con, sql), rows, info = paste(name, mode, "repeat query"))
          if (name %in% unsupported) {
            expect_true(all(is.na(rows$FILE_OFFSET)), info = name)
            expect_equal(rows[names(rows) != "FILE_OFFSET"], biological, info = name)
          } else {
            expected <- oracle[oracle$file == name, c("QNAME", "FILE_OFFSET")]
            expected <- expected[order(expected$QNAME, expected$FILE_OFFSET), ]
            expect_equal(rows$QNAME, expected$QNAME, info = name)
            expect_equal(as.numeric(rows$FILE_OFFSET), expected$FILE_OFFSET, info = name)
          }
        }
      }
    }
    for (name in c("bam_scan_mixed.bam", "bam_scan_mixed.cram")) {
      rows <- dbGetQuery(con, sprintf(paste(
        "SELECT QNAME, FILE_OFFSET FROM read_bam(%s,",
        "region := 'chr1:20-20,chr1:1-30,chr1:20-20') ORDER BY QNAME"),
        dbQuoteString(con, extdata(name))))
      expect_equal(rows$QNAME, c("mapped1", "placed_unmapped"))
      if (endsWith(name, ".cram")) {
        expect_true(all(is.na(rows$FILE_OFFSET)))
      } else {
        expect_equal(as.numeric(rows$FILE_OFFSET), c(7536690, 7536744))
      }
    }
  }
}

test_file_offset()
