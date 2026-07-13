library(tinytest)
library(DBI)

con <- dbConnect(
  duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
)
expect_silent(rduckhts_load(con))

ensembl_fixture_sql <- c(
  "CREATE SCHEMA duckvep_r_core",
  paste(
    "CREATE TABLE duckvep_r_reference AS SELECT '1'::VARCHAR AS chrom,",
    "0::BIGINT AS \"start\", 12::BIGINT AS \"end\",",
    "'ATGAAATAACCC'::VARCHAR seq"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.coord_system AS SELECT",
    "1::BIGINT AS coord_system_id, 1::BIGINT AS species_id,",
    "'chromosome'::VARCHAR AS name, 'GRCh38'::VARCHAR AS version,",
    "1::BIGINT AS rank"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.seq_region AS SELECT",
    "1::BIGINT AS seq_region_id, '1'::VARCHAR AS name,",
    "1::BIGINT AS coord_system_id, 12::BIGINT AS length"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.gene AS SELECT 1::BIGINT AS gene_id,",
    "'protein_coding'::VARCHAR AS biotype, 1::BIGINT AS is_current,",
    "'ENSG_R_TEST'::VARCHAR AS stable_id, 1::BIGINT AS version"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.transcript AS SELECT",
    "1::BIGINT AS transcript_id, 1::BIGINT AS gene_id,",
    "1::BIGINT AS seq_region_id, 1::BIGINT AS seq_region_start,",
    "12::BIGINT AS seq_region_end, 1::BIGINT AS seq_region_strand,",
    "'protein_coding'::VARCHAR AS biotype, 1::BIGINT AS is_current,",
    "'ENST_R_TEST'::VARCHAR AS stable_id, 1::BIGINT AS version"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.exon AS SELECT 1::BIGINT AS exon_id,",
    "1::BIGINT AS seq_region_id, 1::BIGINT AS seq_region_start,",
    "12::BIGINT AS seq_region_end, 1::BIGINT AS seq_region_strand,",
    "-1::BIGINT AS phase, 0::BIGINT AS end_phase,",
    "1::BIGINT AS is_current, 'ENSE_R_TEST'::VARCHAR AS stable_id"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.exon_transcript AS SELECT",
    "1::BIGINT AS exon_id, 1::BIGINT AS transcript_id, 1::BIGINT AS rank"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.translation AS SELECT",
    "1::BIGINT AS translation_id, 1::BIGINT AS transcript_id,",
    "1::BIGINT AS seq_start, 1::BIGINT AS start_exon_id,",
    "9::BIGINT AS seq_end, 1::BIGINT AS end_exon_id,",
    "'ENSP_R_TEST'::VARCHAR AS stable_id, 1::BIGINT AS version"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.attrib_type AS SELECT",
    "NULL::BIGINT AS attrib_type_id, NULL::VARCHAR AS code WHERE false"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.transcript_attrib AS SELECT",
    "NULL::BIGINT AS transcript_id, NULL::BIGINT AS attrib_type_id,",
    "NULL::VARCHAR AS value WHERE false"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.translation_attrib AS SELECT",
    "NULL::BIGINT AS translation_id, NULL::BIGINT AS attrib_type_id,",
    "NULL::VARCHAR AS value WHERE false"
  )
)
for (sql in ensembl_fixture_sql) {
  dbExecute(con, sql)
}
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_ensembl_regions AS SELECT * FROM",
    "duckvep_ensembl_regions(",
    "'duckvep_r_core', 'duckvep_r_reference', 'GRCh38')"
  )
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_ensembl_transcripts AS SELECT * FROM",
    "duckvep_ensembl_transcripts(",
    "'duckvep_r_core', 'duckvep_r_reference', 'GRCh38')"
  )
)
ensembl_model <- dbGetQuery(
  con,
  paste(
    "SELECT CAST(cds_sequence AS VARCHAR) cds_sequence,",
    "CAST(post_cds_bases AS VARCHAR) post_cds_bases,",
    "transcript_flags, exons[1].phase phase",
    "FROM duckvep_r_ensembl_transcripts"
  )
)
expect_identical(ensembl_model$cds_sequence, "ATGAAATAA")
expect_identical(ensembl_model$post_cds_bases, "CCC")
expect_equal(ensembl_model$transcript_flags, 3)
expect_equal(ensembl_model$phase, -1)
ensembl_receipt <- dbGetQuery(
  con,
  paste(
    "SELECT * FROM duckvep_model_receipt(",
    "'duckvep_r_ensembl_regions', 'duckvep_r_ensembl_transcripts',",
    "'Ensembl', '116', 'GRCh38', repeat('a', 64), repeat('b', 64),",
    "'all current transcripts on FASTA-covered assembly regions')"
  )
)
expect_equal(ensembl_receipt$region_count, 1)
expect_equal(ensembl_receipt$reference_base_count, 12)
expect_equal(ensembl_receipt$transcript_count, 1)
expect_equal(nchar(ensembl_receipt$model_sha256), 64)

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
    "1::UTINYINT codon_table, 'ACG'::BLOB post_cds_bases"
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
    "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table,",
    "post_cds_bases",
    "FROM duckvep_r_transcripts ORDER BY seq_region, transcript_start, transcript_index"
  ),
  paste(
    "SELECT transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end,",
    "phase, end_phase FROM duckvep_r_exons ORDER BY transcript_index, exon_cdna_start"
  )
)
load_model <- function(name, model_queries, transcript_coverage_complete = FALSE) {
  arguments <- paste(
    vapply(
      c(name, model_queries),
      function(x) as.character(dbQuoteString(con, x)),
      character(1)
    ),
    collapse = ", "
  )
  if (transcript_coverage_complete) {
    arguments <- paste(
      arguments,
      "transcript_coverage_complete := TRUE",
      sep = ", "
    )
  }
  dbGetQuery(
    con,
    paste0(
      "SELECT loaded FROM duckvep_model_load(",
      arguments,
      ")"
    )
  )
}

loaded <- load_model("r-test", queries)
expect_true(loaded$loaded)

invalid_queries <- queries
invalid_queries[2] <- paste0(
  "SELECT * REPLACE (99::UBIGINT AS transcript_start) FROM (",
  invalid_queries[2],
  ")"
)
expect_error(
  load_model("r-invalid-envelope", invalid_queries),
  pattern = "transcript span is not the outer exon envelope"
)

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

reference_mismatch <- dbGetQuery(
  con,
  paste(
    "SELECT a.status, a.reason FROM unnest(duckvep_annotate(",
    "'r-test', 1::UINTEGER, 124::UBIGINT, 'A', 'C', 0::UBIGINT",
    ")) u(a)"
  )
)
expect_identical(reference_mismatch$status, "unresolved")
expect_identical(reference_mismatch$reason, "reference_mismatch")

frameshift <- dbGetQuery(
  con,
  paste(
    "SELECT a.protein_position, a.reference_amino_acid,",
    "a.alternate_amino_acid FROM unnest(duckvep_annotate(",
    "'r-test', 1::UINTEGER, 124::UBIGINT, 'T', 'TC', 0::UBIGINT",
    ")) u(a)"
  )
)
expect_equal(frameshift$protein_position, 2)
expect_true(is.na(frameshift$reference_amino_acid))
expect_true(is.na(frameshift$alternate_amino_acid))

boundary_insertions <- dbGetQuery(
  con,
  paste(
    "WITH variants(ord, position, reference, alternate) AS (VALUES",
    "(1, 199::UBIGINT, 'G', 'GT'),",
    "(2, 199::UBIGINT, 'G', 'GATG'),",
    "(3, 240::UBIGINT, 'A', 'AATG'),",
    "(4, 250::UBIGINT, 'C', 'CT'))",
    "SELECT ord, a.consequence, a.status, a.reason FROM variants,",
    "LATERAL unnest(duckvep_annotate(",
    "'r-test', 1::UINTEGER, position, reference, alternate, 1::UBIGINT",
    ")) u(a) ORDER BY ord"
  )
)
expect_equal(boundary_insertions$ord, 1:4)
expect_identical(
  boundary_insertions$consequence,
  c(
    "frameshift_variant&splice_region_variant",
    "protein_altering_variant&splice_region_variant",
    "3_prime_UTR_variant",
    "downstream_gene_variant"
  )
)
expect_identical(boundary_insertions$status, rep("supported", 4))
expect_true(all(is.na(boundary_insertions$reason)))

terminal_stop_edits <- dbGetQuery(
  con,
  paste(
    "WITH variants(ord, position, reference, alternate) AS (VALUES",
    "(1, 238::UBIGINT, 'T', 'AC'),",
    "(2, 238::UBIGINT, 'TA', 'T'),",
    "(3, 239::UBIGINT, 'AA', 'A'),",
    "(4, 240::UBIGINT, 'A', 'CG'))",
    "SELECT ord, a.consequence, a.status, a.reason FROM variants,",
    "LATERAL unnest(duckvep_annotate(",
    "'r-test', 1::UINTEGER, position, reference, alternate, 0::UBIGINT",
    ")) u(a) ORDER BY ord"
  )
)
expect_equal(terminal_stop_edits$ord, 1:4)
expect_identical(
  terminal_stop_edits$consequence,
  c(
    "stop_lost",
    "stop_retained_variant",
    "stop_retained_variant",
    "stop_lost"
  )
)
expect_identical(terminal_stop_edits$status, rep("supported", 4))
expect_true(all(is.na(terminal_stop_edits$reason)))

dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_transcripts_no_tail AS",
    "SELECT * EXCLUDE (post_cds_bases) FROM duckvep_r_transcripts"
  )
)
no_tail_queries <- queries
no_tail_queries[2] <- sub(
  "duckvep_r_transcripts",
  "duckvep_r_transcripts_no_tail",
  no_tail_queries[2],
  fixed = TRUE
)
no_tail_queries[2] <- sub(", post_cds_bases", "", no_tail_queries[2], fixed = TRUE)
expect_true(load_model("r-no-tail", no_tail_queries)$loaded)
missing_tail <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence, a.status, a.reason FROM unnest(duckvep_annotate(",
    "'r-no-tail', 1::UINTEGER, 238::UBIGINT, 'TA', 'T', 0::UBIGINT",
    ")) u(a)"
  )
)
expect_identical(missing_tail$consequence, "coding_sequence_variant")
expect_identical(missing_tail$status, "unresolved")
expect_identical(missing_tail$reason, "missing_transcript_tail")

prepared_events <- dbGetQuery(
  con,
  paste(
    "WITH variants(ord, position, reference, alternate) AS (VALUES",
    "(1, 123::UBIGINT, 'GTA', 'GCA'),",
    "(2, 1::UBIGINT, 'AC', 'C'),",
    "(3, 1::UBIGINT, 'C', 'AC'))",
    "SELECT ord, a.consequence, a.status, a.reason FROM variants,",
    "LATERAL unnest(duckvep_annotate(",
    "'r-test', 1::UINTEGER, position, reference, alternate, 0::UBIGINT",
    ")) u(a) ORDER BY ord"
  )
)
expect_equal(prepared_events$ord, 1:3)
expect_identical(prepared_events$consequence[1], "missense_variant")
expect_identical(prepared_events$consequence[2:3], rep("sequence_variant", 2))
expect_identical(prepared_events$status[2:3], rep("unresolved", 2))
expect_identical(
  prepared_events$reason[2:3],
  rep("no_feature_in_loaded_model", 2)
)

complete_queries <- queries
complete_queries[1] <- paste(
  "SELECT seq_region, 1000::UBIGINT AS sequence_length",
  "FROM duckvep_r_regions ORDER BY seq_region"
)
expect_true(load_model("r-complete", complete_queries, TRUE)$loaded)
complete_intergenic <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence, a.status FROM unnest(duckvep_annotate(",
    "'r-complete', 1::UINTEGER, 1::UBIGINT, 'A', 'C', 0::UBIGINT",
    ")) u(a)"
  )
)
expect_identical(complete_intergenic$consequence, "intergenic_variant")
expect_identical(complete_intergenic$status, "supported")
expect_error(
  dbGetQuery(
    con,
    paste(
      "SELECT duckvep_annotate(",
      "'r-complete', 1::UINTEGER, 1001::UBIGINT, 'A', 'C', 0::UBIGINT",
      ")"
    )
  ),
  pattern = "variant span exceeds sequence-region length"
)

dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_transcripts_nmd AS",
    "SELECT * REPLACE (7::UBIGINT AS transcript_flags)",
    "FROM duckvep_r_transcripts"
  )
)
nmd_queries <- queries
nmd_queries[2] <- sub(
  "duckvep_r_transcripts",
  "duckvep_r_transcripts_nmd",
  nmd_queries[2],
  fixed = TRUE
)
expect_true(load_model("r-nmd", nmd_queries)$loaded)
nmd <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence FROM unnest(duckvep_annotate(",
    "'r-nmd', 1::UINTEGER, 160::UBIGINT, 'T', 'C', 0::UBIGINT",
    ")) u(a)"
  )
)
expect_identical(
  nmd$consequence,
  "intron_variant&NMD_transcript_variant"
)
expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('r-test') AS dropped")$dropped)
expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('r-nmd') AS dropped")$dropped)
expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('r-complete') AS dropped")$dropped)

dbDisconnect(con, shutdown = TRUE)
