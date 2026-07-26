library(tinytest)
library(DBI)

test_fastq_qc <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  fastq <- system.file("extdata", "r1.fq", package = "Rduckhts")
  long_qname_fastq <- system.file(
    "extdata", "fastq_long_qname.fq", package = "Rduckhts"
  )
  expect_true(file.exists(fastq))
  expect_true(file.exists(long_qname_fastq))
  quoted_fastq <- as.character(DBI::dbQuoteString(con, fastq))
  quoted_long_qname_fastq <- as.character(
    DBI::dbQuoteString(con, long_qname_fastq)
  )

  summary <- DBI::dbGetQuery(con, sprintf(paste(
    "WITH q AS (",
    "SELECT duckhts_fastq_qc(SEQUENCE, QUALITY) AS qc",
    "FROM read_fastq(%s)",
    ")",
    "SELECT qc.reads::DOUBLE AS reads, qc.bases::DOUBLE AS bases,",
    "qc.q20_bases::DOUBLE AS q20, qc.q30_bases::DOUBLE AS q30,",
    "qc.q40_bases::DOUBLE AS q40, qc.quality_sum::DOUBLE AS quality_sum,",
    "qc.a_bases::DOUBLE AS a, qc.c_bases::DOUBLE AS c,",
    "qc.g_bases::DOUBLE AS g, qc.t_bases::DOUBLE AS t,",
    "qc.n_bases::DOUBLE AS n, qc.other_bases::DOUBLE AS other,",
    "qc.max_read_length::INTEGER AS max_read_length,",
    "len(qc.cycles)::INTEGER AS cycle_count",
    "FROM q"
  ), quoted_fastq))
  expect_equal(
    unname(as.numeric(summary[1, c("reads", "bases", "q20", "q30", "q40")])),
    c(5, 500, 499, 475, 67)
  )
  expect_equal(summary$quality_sum[1], 18425)
  expect_equal(
    unname(as.numeric(summary[1, c("a", "c", "g", "t", "n", "other")])),
    c(152, 98, 94, 156, 0, 0)
  )
  expect_equal(summary$max_read_length[1], 100L)
  expect_equal(summary$cycle_count[1], 100L)

  cycles <- DBI::dbGetQuery(con, sprintf(paste(
    "WITH q AS (",
    "SELECT duckhts_fastq_qc(SEQUENCE, QUALITY) AS qc",
    "FROM read_fastq(%s)",
    ")",
    "SELECT cycle.cycle::INTEGER AS cycle, cycle.bases::DOUBLE AS bases,",
    "cycle.quality_sum::DOUBLE AS quality_sum,",
    "cycle.a_bases::DOUBLE AS a, cycle.c_bases::DOUBLE AS c,",
    "cycle.g_bases::DOUBLE AS g, cycle.t_bases::DOUBLE AS t",
    "FROM q, UNNEST(qc.cycles) AS u(cycle)",
    "WHERE cycle.cycle IN (1, 32, 33, 100)",
    "ORDER BY cycle.cycle"
  ), quoted_fastq))
  expect_equal(cycles$cycle, c(1L, 32L, 33L, 100L))
  expect_equal(cycles$bases, rep(5, 4))
  expect_equal(cycles$quality_sum, c(169, 194, 193, 184))

  long_qname_summary <- DBI::dbGetQuery(con, sprintf(paste(
    "WITH q AS (",
    "SELECT duckhts_fastq_qc(SEQUENCE, QUALITY) AS qc",
    "FROM read_fastq(%s)",
    ") SELECT qc.reads::INTEGER AS reads, qc.bases::INTEGER AS bases FROM q"
  ), quoted_long_qname_fastq))
  expect_equal(long_qname_summary$reads[1], 1L)
  expect_equal(long_qname_summary$bases[1], 4L)
  expect_error(
    DBI::dbGetQuery(
      con, sprintf("SELECT NAME FROM read_fastq(%s)", quoted_long_qname_fastq)
    ),
    "query name exceeds 254 bytes"
  )

  expect_identical(rduckhts_simd_set_backend(con, "scalar"), "scalar")
  DBI::dbExecute(con, paste(
    "CREATE TEMP TABLE fastq_qc_scalar AS",
    "WITH reads AS (",
    "SELECT i % 5 AS grp,",
    "substr(repeat('ACGTNMRWSYKVHDBXacgtn', 7), 1, 1 + i % 129) AS sequence,",
    "substr(repeat('!5?I~', 26), 1, 1 + i % 129) AS quality",
    "FROM range(0, 10000) AS r(i)",
    ")",
    "SELECT grp, duckhts_fastq_qc(sequence, quality) AS qc",
    "FROM reads GROUP BY grp"
  ))
  expect_true(nzchar(rduckhts_simd_set_backend(con, "auto")))
  differences <- DBI::dbGetQuery(con, paste(
    "WITH reads AS (",
    "SELECT i % 5 AS grp,",
    "substr(repeat('ACGTNMRWSYKVHDBXacgtn', 7), 1, 1 + i % 129) AS sequence,",
    "substr(repeat('!5?I~', 26), 1, 1 + i % 129) AS quality",
    "FROM range(0, 10000) AS r(i)",
    "), automatic AS (",
    "SELECT grp, duckhts_fastq_qc(sequence, quality) AS qc",
    "FROM reads GROUP BY grp",
    ")",
    "SELECT count(*)::INTEGER AS n",
    "FROM fastq_qc_scalar AS scalar",
    "FULL OUTER JOIN automatic USING (grp)",
    "WHERE scalar.qc IS DISTINCT FROM automatic.qc"
  ))
  expect_equal(differences$n[1], 0L)

  explicit_limit <- DBI::dbGetQuery(con, paste(
    "WITH q AS (",
    "SELECT duckhts_fastq_qc(repeat('A', 129), repeat('I', 129), 129) AS qc",
    ") SELECT qc.max_read_length::INTEGER AS n FROM q"
  ))
  expect_equal(explicit_limit$n[1], 129L)
  single_cycle_limit <- DBI::dbGetQuery(con, paste(
    "WITH q AS (SELECT duckhts_fastq_qc('A', 'I', 1) AS qc)",
    "SELECT qc.max_read_length::INTEGER AS n,",
    "len(qc.cycles)::INTEGER AS cycle_count FROM q"
  ))
  expect_equal(single_cycle_limit$n[1], 1L)
  expect_equal(single_cycle_limit$cycle_count[1], 1L)
  expect_error(
    DBI::dbGetQuery(con, "SELECT duckhts_fastq_qc('AC', 'I')"),
    "equal length"
  )
  expect_error(
    DBI::dbGetQuery(con, paste(
      "SELECT duckhts_fastq_qc(repeat('A', 129), repeat('I', 129), 128)"
    )),
    "exceeds max_cycles"
  )
  expect_true(nzchar(rduckhts_simd_set_backend(con, "auto")))
}

test_fastq_qc()
