library(tinytest)
library(DBI)

test_bam_materialize <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  path <- system.file("extdata", "bam_materialize.sam", package = "Rduckhts")
  expect_true(nzchar(path))
  for (workers in c(1L, 4L)) {
    dbExecute(con, sprintf("SET threads=%d", workers))
    rduckhts_bam(con, "bam_text", path, standard_tags = TRUE, auxiliary_tags = TRUE, overwrite = TRUE)
    text <- dbGetQuery(con, "SELECT QNAME, CIGAR, SEQ, QUAL, ML, FZ, CG FROM bam_text ORDER BY QNAME")
    expect_equal(text$QNAME, c("r1", "r2", "r3"))
    expect_equal(text$CIGAR, c("5M", "1S2M1I1M", "*"))
    expect_equal(text$SEQ, c("ACGTN", "TGCAN", "*"))
    expect_equal(text$QUAL, c("IIIII", "!!!!!", "*"))
    expect_equal(as.numeric(text$ML[[1]]), c(0, 255))
    expect_equal(as.numeric(text$ML[[2]]), c(0:8, 255))
    expect_equal(as.numeric(text$FZ[[1]]), c(0, 65535))
    expect_equal(length(text$FZ[[2]]), 0L)
    tags <- dbGetQuery(con, paste(
      "SELECT QNAME, map_extract_value(AUXILIARY_TAGS, 'XZ') AS text,",
      "map_extract_value(AUXILIARY_TAGS, 'X7') AS floats FROM bam_text ORDER BY QNAME"))
    expect_equal(tags$text, c("previous", "", NA_character_))
    expect_equal(tags$floats, c("f,0,1.25,1e+30", "f,-1.5,2.25", NA_character_))
    rduckhts_bam(con, "bam_numeric", path, cigar_representation = "binary",
                  sequence_encoding = "nt16", quality_representation = "phred", overwrite = TRUE)
    numeric <- dbGetQuery(con, "SELECT CIGAR, SEQ, QUAL FROM bam_numeric ORDER BY QNAME")
    expect_equal(numeric$CIGAR[[1]], 80)
    expect_equal(numeric$CIGAR[[2]], c(20, 32, 17, 16))
    expect_equal(length(numeric$CIGAR[[3]]), 0L)
    expect_equal(numeric$SEQ[[1]], c(1L, 2L, 4L, 8L, 15L))
    expect_equal(numeric$SEQ[[2]], c(8L, 4L, 2L, 1L, 15L))
    expect_equal(numeric$QUAL[[1]], rep(40L, 5L))
    expect_equal(numeric$QUAL[[2]], rep(0L, 5L))
    large <- system.file("extdata", "bam_scan_all_unplaced.bam", package = "Rduckhts")
    counts <- dbGetQuery(con, sprintf(paste(
      "SELECT count(*) AS n, count(*) FILTER (WHERE SEQ=[4] AND QUAL=[40] AND CIGAR=[]) AS valid",
      "FROM read_bam(%s, sequence_encoding:='nt16', quality_representation:='phred',",
      "cigar_representation:='binary')"), dbQuoteString(con, large)))
    expect_equal(as.numeric(counts$n), 2053)
    expect_equal(as.numeric(counts$valid), 2053)
  }
}

test_bam_materialize()
