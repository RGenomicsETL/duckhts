library(tinytest)
library(DBI)

test_cgranges_api <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  expect_silent(rduckhts_load(con))

  expect_equal(DBI::dbGetQuery(con, "SELECT duckhts_cgranges_create('r_idx') AS ok")$ok[1], TRUE)
  expect_equal(DBI::dbGetQuery(con, "SELECT duckhts_cgranges_add('r_idx', 'chr1', 10, 20) AS ok")$ok[1], TRUE)
  expect_equal(DBI::dbGetQuery(con, "SELECT duckhts_cgranges_add('r_idx', 'chr1', 30, 40) AS ok")$ok[1], TRUE)
  expect_equal(DBI::dbGetQuery(con, "SELECT duckhts_cgranges_index('r_idx') AS ok")$ok[1], TRUE)

  expect_equal(DBI::dbGetQuery(con, "SELECT duckhts_cgranges_create('r_str_idx') AS ok")$ok[1], TRUE)
  expect_equal(DBI::dbGetQuery(con, "SELECT duckhts_cgranges_add('r_str_idx', 'chr1', 50, 60, 'alpha') AS ok")$ok[1], TRUE)
  expect_equal(DBI::dbGetQuery(con, "SELECT duckhts_cgranges_add('r_str_idx', 'chr1', 70, 80, 'beta') AS ok")$ok[1], TRUE)
  expect_equal(DBI::dbGetQuery(con, "SELECT duckhts_cgranges_index('r_str_idx') AS ok")$ok[1], TRUE)

  hits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT interval_ordinal, label, interval_chrom, interval_start, interval_end",
      "FROM duckhts_cgranges_overlaps('r_idx', 'chr1', 35, 36, query_row_id := 9)"
    )
  )
  expect_equal(nrow(hits), 1)
  expect_equal(hits$interval_ordinal[1], 1)
  expect_equal(hits$label[1], 1)
  expect_equal(hits$interval_chrom[1], "chr1")
  expect_equal(hits$interval_start[1], 30)
  expect_equal(hits$interval_end[1], 40)

  str_hits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT interval_ordinal, label, interval_chrom, interval_start, interval_end",
      "FROM duckhts_cgranges_overlaps('r_str_idx', 'chr1', 75, 76)"
    )
  )
  expect_equal(nrow(str_hits), 1)
  expect_equal(str_hits$interval_ordinal[1], 1)
  expect_equal(str_hits$label[1], "beta")
  expect_equal(str_hits$interval_chrom[1], "chr1")
  expect_equal(str_hits$interval_start[1], 70)
  expect_equal(str_hits$interval_end[1], 80)

  DBI::dbExecute(
    con,
    paste(
      "CREATE TABLE cgr_src AS SELECT * FROM (VALUES",
      "('chr2', 100, 120, 11::BIGINT),",
      "('chr2', 140, 170, NULL::BIGINT)",
      ") AS t(chrom, start, \"end\", label)"
    )
  )

  expect_equal(
    DBI::dbGetQuery(
      con,
      paste(
        "SELECT duckhts_cgranges_from_query(",
        "  'qry_idx',",
        "  'SELECT chrom, start, \"end\", label FROM cgr_src',",
        "  'chrom', 'start', 'end', 'label'",
        ") AS ok"
      )
    )$ok[1],
    TRUE
  )
  expect_equal(DBI::dbGetQuery(con, "SELECT duckhts_cgranges_index('qry_idx') AS ok")$ok[1], TRUE)

  num_hits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT interval_ordinal, label, interval_start, interval_end",
      "FROM duckhts_cgranges_overlaps('qry_idx', 'chr2', 110, 111)"
    )
  )
  expect_equal(nrow(num_hits), 1)
  expect_equal(num_hits$interval_ordinal[1], 0)
  expect_equal(num_hits$label[1], 11)
  expect_equal(num_hits$interval_start[1], 100)
  expect_equal(num_hits$interval_end[1], 120)

  null_hits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT interval_ordinal, label, interval_start, interval_end",
      "FROM duckhts_cgranges_overlaps('qry_idx', 'chr2', 140, 170, mode := 'contain')"
    )
  )
  expect_equal(nrow(null_hits), 1)
  expect_equal(null_hits$interval_ordinal[1], 1)
  expect_true(is.na(null_hits$label[1]))
  expect_equal(null_hits$interval_start[1], 140)
  expect_equal(null_hits$interval_end[1], 170)

  expect_error(
    DBI::dbGetQuery(con, "SELECT * FROM duckhts_cgranges_overlaps('missing', 'chr1', 1, 2)"),
    pattern = "unknown index name"
  )
}

test_cgranges_api()
