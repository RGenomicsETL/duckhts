library(tinytest)
library(DBI)

con <- dbConnect(
  duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
)
expect_silent(rduckhts_load(con))

dbExecute(
  con,
  "CREATE TABLE duckvep_r_regions AS SELECT 1::UINTEGER AS seq_region"
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_transcripts AS SELECT",
    "0::UINTEGER transcript_index, 1::UINTEGER seq_region,",
    "100::UBIGINT transcript_start, 250::UBIGINT transcript_end,",
    "1::TINYINT strand, 0::UINTEGER gene_index, 3::UBIGINT transcript_flags,",
    "120::UBIGINT cds_start, 240::UBIGINT cds_end,",
    "'ATGGTACGTACGTACGTACGTACGTACGTACTACGTACGTACGTACGTACGTACGTACGTACGTACTGGTAA'::BLOB cds_sequence,",
    "1::UTINYINT codon_table"
  )
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_exons AS SELECT * FROM (VALUES",
    "(0::UINTEGER, 100::UBIGINT, 150::UBIGINT, 1::UBIGINT, 51::UBIGINT, 0::TINYINT, 0::TINYINT),",
    "(0::UINTEGER, 200::UBIGINT, 250::UBIGINT, 52::UBIGINT, 102::UBIGINT, 0::TINYINT, 0::TINYINT)",
    ") t(transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end, phase, end_phase)"
  )
)

queries <- c(
  "SELECT seq_region FROM duckvep_r_regions ORDER BY seq_region",
  paste(
    "SELECT transcript_index, seq_region, transcript_start, transcript_end,",
    "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table",
    "FROM duckvep_r_transcripts ORDER BY seq_region, transcript_start, transcript_index"
  ),
  paste(
    "SELECT transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end,",
    "phase, end_phase FROM duckvep_r_exons ORDER BY transcript_index, exon_cdna_start"
  )
)
loaded <- dbGetQuery(
  con,
  paste0(
    "SELECT loaded FROM duckvep_model_load(",
    paste(
      vapply(
        c("r-test", queries),
        function(x) as.character(dbQuoteString(con, x)),
        character(1)
      ),
      collapse = ", "
    ),
    ")"
  )
)
expect_true(loaded$loaded)

annotation <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence, a.impact, a.status, a.cds_position, a.protein_position",
    "FROM unnest(duckvep_annotate(",
    "'r-test', 1::UINTEGER, 124::UBIGINT, 'T', 'C', 0::UBIGINT",
    ")) u(a)"
  )
)
expect_identical(annotation$consequence, "missense_variant")
expect_identical(annotation$impact, "MODERATE")
expect_identical(annotation$status, "supported")
expect_equal(annotation$cds_position, 5)
expect_equal(annotation$protein_position, 2)
expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('r-test') AS dropped")$dropped)

dbDisconnect(con, shutdown = TRUE)
