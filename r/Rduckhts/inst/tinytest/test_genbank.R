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
      "WHERE regexp_extract(attributes, 'ID=([^;]*)', 1) = 'CDS-2'",
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
  expect_error(
    rduckhts_genbank_to_fasta(con, genbank_path, line_width = Inf),
    "line_width must be a positive whole number"
  )
  expect_error(
    rduckhts_genbank_to_fasta(con, genbank_path, line_width = 3e9),
    "line_width must be a positive whole number"
  )
  expect_error(
    rduckhts_genbank_to_fasta(con, genbank_path, output_path = ""),
    "output_path must be NULL or one non-empty character string"
  )
  expect_error(
    rduckhts_genbank_to_fasta(con, genbank_path, overwrite = NA),
    "overwrite must be TRUE or FALSE"
  )

  # A FEATURES table that reaches // without a sequence section is an error,
  # not a short result, for the reader and the converter alike.
  malformed_path <- tempfile("rduckhts_genbank_malformed_", fileext = ".gb")
  on.exit(unlink(malformed_path), add = TRUE)
  writeLines(
    c(
      "LOCUS       NOSECTION                100 bp    DNA     linear   PHG 01-JAN-2000",
      "FEATURES             Location/Qualifiers",
      "     gene            1..30",
      "COMMENT     no sequence section follows the features",
      "//"
    ),
    malformed_path
  )
  expect_error(
    rduckhts_genbank(con, "malformed_features", malformed_path),
    "FEATURES table is not followed by a sequence section"
  )
  malformed_fasta <- tempfile("rduckhts_genbank_malformed_", fileext = ".fa")
  expect_error(
    rduckhts_genbank_to_fasta(con, malformed_path, output_path = malformed_fasta),
    "FEATURES table is not followed by a sequence section"
  )
  expect_false(file.exists(malformed_fasta))

  # Refusing to overwrite leaves the existing output untouched.
  expect_error(
    rduckhts_genbank_to_fasta(con, genbank_path, output_path = fasta_path),
    "already exists"
  )
  expect_equal(
    dbGetQuery(
      con,
      paste0(
        "SELECT length(SEQUENCE)::INTEGER AS bp FROM read_fasta(",
        as.character(dbQuoteString(con, fasta_path)),
        ")"
      )
    )$bp,
    5386L
  )
}

test_genbank_semantics <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # Segments come out in biological order with the CDS phase carried across
  # them: complement(join(1..4,10..17)) reads 10..17 first (phase 0), then
  # 1..4 (phase 1); join(complement(10..17),complement(1..4)) with
  # /codon_start=2 starts at phase 1 and continues at 2.
  phase_path <- system.file(
    "extdata",
    "genbank_phase.gb",
    package = "Rduckhts",
    mustWork = TRUE
  )
  expect_silent(rduckhts_genbank(con, "phase_features", phase_path))
  rows <- dbGetQuery(
    con,
    paste(
      "SELECT regexp_extract(attributes, 'locus_tag=([^;]*)', 1) AS tag,",
      "start::INTEGER AS start, \"end\"::INTEGER AS end, strand, frame",
      "FROM phase_features"
    )
  )
  expect_equal(rows$tag, c("p1", "p1", "p2", "p2", "p3", "p3", "p4", "p4", "p5", "p5"))
  expect_equal(rows$start, c(1L, 10L, 10L, 1L, 10L, 1L, 1L, 10L, 10L, 1L))
  expect_equal(rows$strand, c("+", "+", "-", "-", "-", "-", "+", "+", "-", "-"))
  expect_equal(rows$frame, c("0", "2", "0", "1", "1", "2", "0", "2", ".", "."))

  # Parent links by /locus_tag whichever of the gene and its child comes first.
  order_path <- system.file(
    "extdata",
    "genbank_gene_order.gb",
    package = "Rduckhts",
    mustWork = TRUE
  )
  expect_silent(rduckhts_genbank(con, "order_features", order_path))
  parents <- dbGetQuery(
    con,
    paste(
      "SELECT feature, start::INTEGER AS start,",
      "regexp_extract(attributes, 'Parent=([^;]*)', 1) AS parent",
      "FROM order_features ORDER BY start, feature"
    )
  )
  expect_equal(parents$feature, c("CDS", "gene", "CDS", "gene", "tRNA"))
  expect_equal(parents$parent, c("gene-g1", "", "gene-g2", "", ""))
}

test_genbank_named_attributes <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  genbank_path <- system.file(
    "extdata",
    "phix174.gb",
    package = "Rduckhts",
    mustWork = TRUE
  )

  # Named qualifier columns equal attributes_map lookups, including absent keys.
  keys <- c("locus_tag", "product", "note", "db_xref", "ID", "Parent", "missing")
  expect_silent(rduckhts_genbank(
    con,
    "phix_named",
    genbank_path,
    attributes = keys,
    attributes_map = TRUE,
    overwrite = TRUE
  ))
  named <- dbGetQuery(con, "SELECT * FROM phix_named")
  expect_equal(names(named)[10:17], c("attributes_map", keys))
  for (key in keys) {
    same <- dbGetQuery(
      con,
      sprintf(
        "SELECT bool_and(%s IS NOT DISTINCT FROM attributes_map[%s]) AS same FROM phix_named",
        dbQuoteIdentifier(con, key),
        dbQuoteString(con, key)
      )
    )$same
    expect_true(same, info = key)
  }
  expect_true(any(!is.na(named$product)))
  expect_true(all(is.na(named$missing)))

  expect_silent(rduckhts_genbank(
    con,
    "phix_unmapped",
    genbank_path,
    attributes = "locus_tag",
    overwrite = TRUE
  ))
  expect_equal(names(dbGetQuery(con, "SELECT * FROM phix_unmapped LIMIT 1"))[10L], "locus_tag")
  expect_silent(rduckhts_genbank(
    con,
    "phix_no_keys",
    genbank_path,
    attributes = character(),
    overwrite = TRUE
  ))
  expect_equal(ncol(dbGetQuery(con, "SELECT * FROM phix_no_keys LIMIT 1")), 9L)

  expect_error(rduckhts_genbank(con, "bad", genbank_path, attributes = 1), "character vector")
  expect_error(rduckhts_genbank(con, "bad", genbank_path, attributes = NA_character_), "NA")
  expect_error(rduckhts_genbank(con, "bad", genbank_path, attributes = ""), "empty")
  expect_error(rduckhts_genbank(con, "bad", genbank_path, attributes = "START"), "collides")
  expect_error(rduckhts_genbank(con, "bad", genbank_path, attributes = c("ID", "id")), "distinct")
}

test_genbank()
test_genbank_semantics()
test_genbank_named_attributes()
