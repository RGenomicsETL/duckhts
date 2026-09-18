# GenBank reader and wrapper tests
library(tinytest)
library(DBI)

test_genbank <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  genbank_path <- system.file(
    "extdata",
    "phix174.gb",
    package = "Rduckhts",
    mustWork = TRUE
  )

  expect_silent(rduckhts_genbank(
    con,
    "phix_features",
    genbank_path,
    overwrite = TRUE
  ))

  counts <- dbGetQuery(
    con,
    paste(
      "SELECT feature, count(*)::INTEGER AS n FROM phix_features",
      "GROUP BY feature ORDER BY feature"
    )
  )
  expect_equal(
    counts$feature,
    c("CDS", "gene", "mRNA", "misc_feature", "rep_origin", "variation")
  )
  expect_equal(counts$n, c(14L, 14L, 2L, 2L, 1L, 4L))

  # The record-level `source` feature describes the whole record, not an
  # interval within it, so it is dropped rather than emitted.
  expect_equal(
    dbGetQuery(
      con,
      paste(
        "SELECT count(*)::INTEGER AS n FROM phix_features",
        "WHERE feature = 'source'"
      )
    )$n,
    0L
  )

  # phiX174 is circular: join(3981..5386,1..136) wraps the origin and each
  # segment becomes its own row, tied together by the synthesized ID.
  wrapped <- dbGetQuery(
    con,
    paste(
      "SELECT start::INTEGER AS start, \"end\"::INTEGER AS end",
      "FROM phix_features",
      "WHERE regexp_extract(attributes, 'ID=([^;]*)', 1) = 'CDS-1'",
      "ORDER BY start"
    )
  )
  expect_equal(nrow(wrapped), 2L)
  expect_equal(wrapped$start, c(1L, 3981L))
  expect_equal(wrapped$end, c(136L, 5386L))

  expect_silent(rduckhts_genbank(
    con,
    "phix_mapped",
    genbank_path,
    attributes_map = TRUE,
    overwrite = TRUE
  ))
  expect_equal(
    dbGetQuery(
      con,
      paste(
        "SELECT attributes_map['locus_tag'] AS locus_tag FROM phix_mapped",
        "WHERE feature = 'CDS' AND start = 1001"
      )
    )$locus_tag,
    "phiX174p06"
  )

  fasta_path <- tempfile("rduckhts_genbank_", fileext = ".fa")
  on.exit(unlink(fasta_path), add = TRUE)
  written <- rduckhts_genbank_to_fasta(
    con,
    genbank_path,
    output_path = fasta_path,
    overwrite = TRUE
  )
  expect_true(written$success)
  expect_equal(as.integer(written$records_written), 1L)
  expect_true(file.exists(fasta_path))

  # The FASTA carries the same name read_genbank reports as seqname, so
  # feature coordinates land on the contig of that name.
  roundtrip <- dbGetQuery(
    con,
    paste0(
      "SELECT NAME, length(SEQUENCE)::INTEGER AS bp FROM read_fasta(",
      as.character(dbQuoteString(con, fasta_path)),
      ")"
    )
  )
  expect_equal(roundtrip$NAME, "NC_001422.1")
  expect_equal(roundtrip$bp, 5386L)
  expect_equal(
    dbGetQuery(
      con,
      "SELECT DISTINCT seqname FROM phix_features"
    )$seqname,
    "NC_001422.1"
  )

  expect_error(
    rduckhts_genbank(con, "phix_features", genbank_path),
    "already exists"
  )
  expect_error(
    rduckhts_genbank(con, "bad_path", character()),
    "path must be one non-empty character string"
  )
  expect_error(
    rduckhts_genbank(con, "bad_map", genbank_path, attributes_map = NA),
    "attributes_map must be TRUE or FALSE"
  )
  expect_error(
    rduckhts_genbank_to_fasta(con, genbank_path, line_width = 0),
    "line_width must be a positive whole number"
  )
}

test_genbank()
