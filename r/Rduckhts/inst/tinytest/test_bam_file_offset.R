# Tests for FILE_OFFSET column in read_bam()
#
# FILE_OFFSET is a UBIGINT column that exposes the BGZF virtual file offset
# after each record.  It enables ORDER BY FILE_OFFSET in SQL window functions
# to reproduce exact BAM file order, required for streaming dedup algorithms
# such as WisecondorX's larp/larp2 state machine.
library(tinytest)
library(DBI)

test_file_offset <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)

  rduckhts_load(con)

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

  sample_id_pos   <- which(col_names == "SAMPLE_ID")
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
  expect_true(all(offsets > 0L),    info = "FILE_OFFSET values are positive")
  expect_identical(
    length(unique(offsets)), length(offsets),
    info = "FILE_OFFSET is unique per record (no two records share an offset)"
  )
  expect_true(
    all(diff(sort(offsets)) > 0L),
    info = "FILE_OFFSET values are strictly increasing in scan order"
  )

  # ---- Window function: LAG(FILE_OFFSET) usable without errors ----
  # Verify the column works as an ORDER BY key in a window function; this is
  # the foundation of the larp/larp2 SQL dedup used by WisecondorX.
  #
  # NOTE: window functions cannot be nested inside aggregate FILTER clauses in
  # DuckDB ("Binder Error: aggregate function calls cannot contain window
  # function calls"). Compute the window in a CTE first, then aggregate.
  lag_q <- dbGetQuery(con, sprintf("
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
  ", bam_path))

  expect_identical(lag_q$total[1], lag_q$distinct_offsets[1],
                   info = "no duplicate FILE_OFFSET values")
  expect_identical(lag_q$order_violations[1], 0L,
                   info = "no order violations: FILE_OFFSET strictly increases")

  dbDisconnect(con, shutdown = TRUE)
}

test_file_offset()
