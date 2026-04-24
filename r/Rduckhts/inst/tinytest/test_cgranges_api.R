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

  DBI::dbExecute(
    con,
    paste(
      "CREATE TABLE cgr_probe_src AS SELECT * FROM (VALUES",
      "(10::BIGINT, 'chr2', 100, 105),",
      "(20::BIGINT, 'chr2', 160, 161),",
      "(30::BIGINT, 'chr2', 500, 510),",
      "(40::BIGINT, 'chr2', 140, 170)",
      ") AS t(qid, chrom, start, \"end\")"
    )
  )

  scalar_probe <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT qid,",
      "  duckhts_cgranges_has_overlap('qry_idx', chrom, start, \"end\") AS has_overlap,",
      "  duckhts_cgranges_count_overlaps('qry_idx', chrom, start, \"end\") AS n_overlap",
      "FROM cgr_probe_src",
      "ORDER BY qid"
    )
  )
  expect_equal(scalar_probe$qid, c(10, 20, 30, 40))
  expect_equal(scalar_probe$has_overlap, c(TRUE, TRUE, FALSE, TRUE))
  expect_equal(scalar_probe$n_overlap, c(1, 1, 0, 1))

  scalar_contain <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT qid,",
      "  duckhts_cgranges_has_overlap('qry_idx', chrom, start, \"end\", 'contain') AS has_contained,",
      "  duckhts_cgranges_count_overlaps('qry_idx', chrom, start, \"end\", 'contain') AS n_contained",
      "FROM cgr_probe_src",
      "ORDER BY qid"
    )
  )
  expect_equal(scalar_contain$has_contained, c(FALSE, FALSE, FALSE, TRUE))
  expect_equal(scalar_contain$n_contained, c(0, 0, 0, 1))

  null_probe <- DBI::dbGetQuery(
    con,
    "SELECT duckhts_cgranges_has_overlap('qry_idx', NULL::VARCHAR, 1, 2) AS hit"
  )
  expect_true(is.na(null_probe$hit[1]))

  list_hits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT p.qid, hit.interval_ordinal, hit.label, hit.label_type,",
      "  hit.interval_chrom, hit.interval_start, hit.interval_end",
      "FROM cgr_probe_src AS p",
      "CROSS JOIN UNNEST(duckhts_cgranges_overlaps_list('qry_idx', p.chrom, p.start, p.\"end\")) AS u(hit)",
      "ORDER BY p.qid, hit.interval_ordinal"
    )
  )
  expect_equal(list_hits$qid, c(10, 20, 40))
  expect_equal(list_hits$interval_ordinal, c(0, 1, 1))
  expect_equal(list_hits$label, c("11", NA, NA))
  expect_equal(list_hits$label_type, c("BIGINT", "BIGINT", "BIGINT"))
  expect_equal(list_hits$interval_start, c(100, 140, 140))
  expect_equal(list_hits$interval_end, c(120, 170, 170))

  list_contain <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT p.qid, hit.interval_ordinal, hit.label, hit.label_type, hit.interval_start, hit.interval_end",
      "FROM cgr_probe_src AS p",
      "CROSS JOIN UNNEST(duckhts_cgranges_overlaps_list('qry_idx', p.chrom, p.start, p.\"end\", 'contain')) AS u(hit)",
      "ORDER BY p.qid, hit.interval_ordinal"
    )
  )
  expect_equal(list_contain$qid, 40)
  expect_equal(list_contain$interval_ordinal, 1)
  expect_true(is.na(list_contain$label[1]))
  expect_equal(list_contain$label_type, "BIGINT")

  null_list_probe <- DBI::dbGetQuery(
    con,
    "SELECT duckhts_cgranges_overlaps_list('qry_idx', NULL::VARCHAR, 1, 2) IS NULL AS is_null"
  )
  expect_equal(null_list_probe$is_null[1], TRUE)

  bed_path <- system.file("extdata", "targets.bed", package = "Rduckhts")
  bed_query <- paste0(
    "SELECT chrom, start, \"end\" FROM read_bed(",
    DBI::dbQuoteString(con, bed_path),
    ")"
  )
  expect_equal(
    DBI::dbGetQuery(
      con,
      paste0(
        "SELECT duckhts_cgranges_from_query('bed_from_query_idx', ",
        DBI::dbQuoteString(con, bed_query),
        ", 'chrom', 'start', 'end') AS ok"
      )
    )$ok[1],
    TRUE
  )
  bed_from_query_hits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT interval_ordinal, label, interval_chrom, interval_start, interval_end",
      "FROM duckhts_cgranges_overlaps('bed_from_query_idx', 'CHROMOSOME_I', 5, 15)",
      "ORDER BY interval_ordinal"
    )
  )
  expect_equal(bed_from_query_hits$interval_ordinal, c(0, 1))
  expect_equal(bed_from_query_hits$label, c(0, 1))
  expect_equal(bed_from_query_hits$interval_chrom, c("CHROMOSOME_I", "CHROMOSOME_I"))

  expect_equal(DBI::dbGetQuery(con, "SELECT duckhts_cgranges_create('bed_provider_idx') AS ok")$ok[1], TRUE)
  expect_equal(
    DBI::dbGetQuery(con, "SELECT duckhts_cgranges_add('bed_provider_idx', 'CHROMOSOME_I', 5, 15) AS ok")$ok[1],
    TRUE
  )
  expect_equal(DBI::dbGetQuery(con, "SELECT duckhts_cgranges_index('bed_provider_idx') AS ok")$ok[1], TRUE)
  bed_provider_hits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT count(*) AS overlapping_rows,",
      "  sum(duckhts_cgranges_count_overlaps('bed_provider_idx', chrom, start, \"end\"))::BIGINT AS total_hits",
      "FROM read_bed(", DBI::dbQuoteString(con, bed_path), ")",
      "WHERE duckhts_cgranges_has_overlap('bed_provider_idx', chrom, start, \"end\")"
    )
  )
  expect_equal(bed_provider_hits$overlapping_rows[1], 2)
  expect_equal(bed_provider_hits$total_hits[1], 2)

  bed_provider_list_hits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT q.chrom, q.start, q.\"end\", hit.interval_ordinal, hit.label, hit.label_type,",
      "  hit.interval_start, hit.interval_end",
      "FROM read_bed(", DBI::dbQuoteString(con, bed_path), ") AS q",
      "CROSS JOIN UNNEST(duckhts_cgranges_overlaps_list('bed_provider_idx', q.chrom, q.start, q.\"end\")) AS u(hit)",
      "ORDER BY q.chrom, q.start, hit.interval_ordinal"
    )
  )
  expect_equal(nrow(bed_provider_list_hits), 2)
  expect_equal(bed_provider_list_hits$chrom, c("CHROMOSOME_I", "CHROMOSOME_I"))
  expect_equal(bed_provider_list_hits$start, c(0, 10))
  expect_equal(bed_provider_list_hits$label, c("0", "0"))
  expect_equal(bed_provider_list_hits$label_type, c("ORDINAL", "ORDINAL"))
  expect_equal(bed_provider_list_hits$interval_start, c(5, 5))
  expect_equal(bed_provider_list_hits$interval_end, c(15, 15))

  bulk_hits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT query_row_id, interval_ordinal, label, interval_start, interval_end",
      "FROM duckhts_cgranges_overlaps_bulk(",
      "  'qry_idx',",
      "  'SELECT qid, chrom, start, \"end\" FROM cgr_probe_src',",
      "  'chrom', 'start', 'end',",
      "  query_row_id_col := 'qid'",
      ")",
      "ORDER BY query_row_id, interval_ordinal"
    )
  )
  expect_equal(nrow(bulk_hits), 3)
  expect_equal(bulk_hits$query_row_id, c(10, 20, 40))
  expect_equal(bulk_hits$interval_ordinal, c(0, 1, 1))
  expect_equal(as.character(bulk_hits$label), c("11", NA_character_, NA_character_))
  expect_equal(bulk_hits$interval_start, c(100, 140, 140))
  expect_equal(bulk_hits$interval_end, c(120, 170, 170))

  auto_qid_hits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT query_row_id, interval_ordinal, label, interval_start, interval_end",
      "FROM duckhts_cgranges_overlaps_bulk(",
      "  'qry_idx',",
      "  'SELECT chrom, start, \"end\" FROM cgr_probe_src',",
      "  'chrom', 'start', 'end'",
      ")",
      "ORDER BY query_row_id, interval_ordinal"
    )
  )
  expect_equal(auto_qid_hits$query_row_id, c(1, 2, 4))
  expect_equal(auto_qid_hits$interval_ordinal, c(0, 1, 1))

  contain_hits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT query_row_id, interval_ordinal, label, interval_start, interval_end",
      "FROM duckhts_cgranges_overlaps_bulk(",
      "  'qry_idx',",
      "  'SELECT qid, chrom, start, \"end\" FROM cgr_probe_src',",
      "  'chrom', 'start', 'end',",
      "  mode := 'contain',",
      "  query_row_id_col := 'qid'",
      ")",
      "ORDER BY query_row_id, interval_ordinal"
    )
  )
  expect_equal(contain_hits$query_row_id, 40)
  expect_equal(contain_hits$interval_ordinal, 1)
  expect_true(is.na(contain_hits$label[1]))

  expect_error(
    DBI::dbGetQuery(con, "SELECT * FROM duckhts_cgranges_overlaps('missing', 'chr1', 1, 2)"),
    pattern = "unknown index name"
  )
}

test_cgranges_api()
