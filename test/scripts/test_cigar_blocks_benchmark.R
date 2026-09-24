library(DBI)
library(duckdb)
source("scripts/cigar_aligned_blocks_benchmark.R")

check_cigar_benchmark <- function(extension) {
  con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = "true"), shared_home = FALSE))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  dbExecute(con, paste("LOAD", dbQuoteString(con, normalizePath(extension, mustWork = TRUE))))
  dbExecute(con, "CREATE TABLE records AS SELECT * FROM (VALUES
    ('duplicate', 0, 10::BIGINT, [80]::UINTEGER[]),
    ('duplicate', 0, 10::BIGINT, [80]::UINTEGER[]),
    ('empty', 4, 0::BIGINT, []::UINTEGER[]),
    ('null', 4, 0::BIGINT, NULL::UINTEGER[]),
    ('soft', 0, 10::BIGINT, [36]::UINTEGER[]),
    ('deletion', 0, 10::BIGINT, [34]::UINTEGER[]),
    ('position', 0, NULL::BIGINT, [80]::UINTEGER[])
    ) v(QNAME, FLAG, POS, CIGAR)")
  stopifnot(validate_blocks(con, "records")$records == "7")
  expected <- dbGetQuery(con, workload_sql("records", "blocks in SQL"))
  observed <- dbGetQuery(con, workload_sql("records", "cigar_aligned_blocks"))
  stopifnot(identical(observed, expected), expected$reads == "7", expected$blocks == "2")
  sentinel <- dbGetQuery(con, sprintf(
    "SELECT QNAME, ref_start IS NULL AS missing, len(width) AS blocks
     FROM (SELECT QNAME, %s FROM records) ORDER BY QNAME", sql_blocks))
  stopifnot(identical(sentinel$missing, c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE)),
            identical(sentinel$blocks, c(0, 1, 1, NA, NA, NA, 0)))

  dbExecute(con, "CREATE TABLE duplicates AS SELECT * FROM records WHERE QNAME = 'duplicate'")
  expected <- dbGetQuery(con, workload_sql("duplicates", "blocks in SQL"))
  dbExecute(con, "CREATE SCHEMA faulty")
  dbExecute(con, "CREATE MACRO faulty.cigar_aligned_blocks(cigar, pos) AS (
    SELECT {'ref_start':list_transform(b.ref_start, lambda p: p + 1),
            'query_start':b.query_start, 'width':b.width}
    FROM (SELECT main.cigar_aligned_blocks(cigar, pos) AS b))")
  dbExecute(con, "SET search_path = 'faulty,main'")
  observed <- dbGetQuery(con, workload_sql("duplicates", "cigar_aligned_blocks"))
  stopifnot(!identical(observed, expected), identical(observed$blocks, expected$blocks))
  failure <- tryCatch(validate_blocks(con, "duplicates"), error = conditionMessage)
  stopifnot(identical(failure, "2 physical-record geometry mismatches out of 2"))

  # Identical physical rows demonstrate the XOR cancellation independently.
  cancellation <- dbGetQuery(con, "SELECT
    bit_xor(hash(QNAME, FLAG, [POS], [0::BIGINT], [5::BIGINT])) =
    bit_xor(hash(QNAME, FLAG, [POS + 1], [0::BIGINT], [5::BIGINT])) AS cancels
    FROM duplicates")
  stopifnot(cancellation$cancels)

  dbExecute(con, "CREATE OR REPLACE MACRO faulty.cigar_aligned_blocks(cigar, pos) AS (
    CASE WHEN cigar IS NULL OR len(cigar) = 0 OR pos IS NULL
    THEN {'ref_start':[]::BIGINT[], 'query_start':[]::BIGINT[], 'width':[]::BIGINT[]}
    ELSE main.cigar_aligned_blocks(cigar, pos) END)")
  observed <- dbGetQuery(con, workload_sql("records", "cigar_aligned_blocks"))
  expected <- dbGetQuery(con, workload_sql("records", "blocks in SQL"))
  stopifnot(!identical(observed, expected))
  failure <- tryCatch(validate_blocks(con, "records"), error = conditionMessage)
  stopifnot(identical(failure, "3 physical-record geometry mismatches out of 7"))

  # A valid struct containing NULL children is not a NULL struct.
  dbExecute(con, "CREATE OR REPLACE MACRO faulty.cigar_aligned_blocks(cigar, pos) AS (
    CASE WHEN cigar IS NULL OR len(cigar) = 0 OR pos IS NULL
    THEN {'ref_start':NULL::BIGINT[], 'query_start':NULL::BIGINT[], 'width':NULL::BIGINT[]}
    ELSE main.cigar_aligned_blocks(cigar, pos) END)")
  stopifnot(!identical(dbGetQuery(con, workload_sql("records", "cigar_aligned_blocks")), expected))
  failure <- tryCatch(validate_blocks(con, "records"), error = conditionMessage)
  stopifnot(identical(failure, "3 physical-record geometry mismatches out of 7"))
  cat("CIGAR benchmark NULL/duplicate/failure controls passed\n")
}

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
check_cigar_benchmark(arguments[[1L]])
