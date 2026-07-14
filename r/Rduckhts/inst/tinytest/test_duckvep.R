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
    "1::BIGINT AS attrib_type_id, 'MANE_Select'::VARCHAR AS code"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.transcript_attrib AS SELECT",
    "1::BIGINT AS transcript_id, 1::BIGINT AS attrib_type_id,",
    "'NM_R_TEST.1'::VARCHAR AS value"
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
    "CAST(pre_cds_sequence AS VARCHAR) pre_cds_sequence,",
    "CAST(post_cds_sequence AS VARCHAR) post_cds_sequence,",
    "transcript_flags, mane_select_refseq, exons[1].phase phase",
    "FROM duckvep_r_ensembl_transcripts"
  )
)
expect_identical(ensembl_model$cds_sequence, "ATGAAATAA")
expect_identical(ensembl_model$post_cds_bases, "CCC")
expect_identical(ensembl_model$pre_cds_sequence, "")
expect_identical(ensembl_model$post_cds_sequence, "CCC")
expect_equal(ensembl_model$transcript_flags, 1027)
expect_identical(ensembl_model$mane_select_refseq, "NM_R_TEST.1")
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
expect_equal(ensembl_receipt$transcript_flank_base_count, 3)
expect_equal(nchar(ensembl_receipt$model_sha256), 64)
expect_identical(
  names(ensembl_receipt)[1:6],
  c(
    "source_name",
    "source_version",
    "assembly",
    "source_manifest_sha256",
    "reference_sha256",
    "transcript_filter"
  )
)
expect_identical(ensembl_receipt$source_name, "Ensembl")
expect_identical(ensembl_receipt$source_version, "116")
expect_identical(ensembl_receipt$assembly, "GRCh38")
expect_identical(
  ensembl_receipt$transcript_filter,
  "all current transcripts on FASTA-covered assembly regions"
)

dbExecute(
  con,
  "INSERT INTO duckvep_r_core.attrib_type VALUES (2, '_rna_edit')"
)
dbExecute(
  con,
  "INSERT INTO duckvep_r_core.transcript_attrib VALUES (1, 2, '4 5 A')"
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_ensembl_rna_edit AS SELECT * FROM",
    "duckvep_ensembl_transcripts(",
    "'duckvep_r_core', 'duckvep_r_reference', 'GRCh38')"
  )
)
ensembl_rna_edit <- dbGetQuery(
  con,
  paste(
    "SELECT CAST(cds_sequence AS VARCHAR) cds_sequence,",
    "sequence_withheld_reason FROM duckvep_r_ensembl_rna_edit"
  )
)
expect_true(is.na(ensembl_rna_edit$cds_sequence))
expect_identical(ensembl_rna_edit$sequence_withheld_reason, "rna_edit")

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
    "1::UTINYINT codon_table,",
    "'TACGTACGTACGTACGTACG'::BLOB pre_cds_sequence,",
    "'ACGTACGTAC'::BLOB post_cds_sequence"
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
    "pre_cds_sequence, post_cds_sequence",
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

dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_partial_cds_transcripts AS SELECT",
    "0::UINTEGER transcript_index, 1::UINTEGER seq_region,",
    "100::UBIGINT transcript_start, 113::UBIGINT transcript_end,",
    "1::TINYINT strand, 0::UINTEGER gene_index,",
    "35::UBIGINT transcript_flags, 100::UBIGINT cds_start,",
    "113::UBIGINT cds_end, 'ATGAAACCCGGGTT'::BLOB cds_sequence,",
    "1::UTINYINT codon_table, ''::BLOB pre_cds_sequence,",
    "''::BLOB post_cds_sequence"
  )
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_partial_cds_exons AS SELECT",
    "0::UINTEGER transcript_index, 100::UBIGINT exon_start,",
    "113::UBIGINT exon_end, 1::UBIGINT exon_cdna_start,",
    "14::UBIGINT exon_cdna_end, 0::TINYINT phase,",
    "2::TINYINT end_phase"
  )
)
partial_cds_queries <- c(
  "SELECT seq_region FROM duckvep_r_regions ORDER BY seq_region",
  paste(
    "SELECT * FROM duckvep_r_partial_cds_transcripts",
    "ORDER BY seq_region, transcript_start, transcript_index"
  ),
  paste(
    "SELECT * FROM duckvep_r_partial_cds_exons",
    "ORDER BY transcript_index, exon_cdna_start"
  )
)
expect_true(load_model("r-partial-cds-end", partial_cds_queries)$loaded)
partial_cds_annotation <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence, a.status, a.reason",
    "FROM unnest(duckvep_annotate(",
    "'r-partial-cds-end', 1::UINTEGER, 106::UBIGINT,",
    "'C', 'CA', 0::UBIGINT)) u(a)"
  )
)
expect_identical(partial_cds_annotation$consequence, "frameshift_variant")
expect_identical(partial_cds_annotation$status, "supported")
expect_true(is.na(partial_cds_annotation$reason))

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

compact_annotation <- dbGetQuery(
  con,
  paste(
    "SELECT a.* FROM unnest(duckvep_annotate_compact(",
    "'r-test', 1::UINTEGER, 124::UBIGINT, 'T', 'C', 0::UBIGINT",
    ")) u(a)"
  )
)
expect_equal(compact_annotation$transcript_index, 0)
expect_equal(compact_annotation$gene_index, 0)
expect_equal(compact_annotation$consequence_mask, 8192)
expect_equal(compact_annotation$region_mask, 16)
expect_equal(compact_annotation$impact_code, 2)
expect_equal(compact_annotation$status_code, 0)
expect_equal(compact_annotation$reason_code, 0)
expect_equal(compact_annotation$cdna_position, 25)
expect_equal(compact_annotation$cds_position, 5)
expect_equal(compact_annotation$protein_position, 2)
expect_equal(compact_annotation$reference_amino_acid_code, utf8ToInt("V"))
expect_equal(compact_annotation$alternate_amino_acid_code, utf8ToInt("A"))
expect_equal(compact_annotation$nmd_prediction_code, 0)
expect_equal(compact_annotation$nmd_escape_reasons, 0)

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
    "(4, 240::UBIGINT, 'A', 'CG'),",
    "(5, 239::UBIGINT, 'AA', 'CA'),",
    "(6, 240::UBIGINT, 'AA', 'TA'))",
    "SELECT ord, a.consequence, a.status, a.reason FROM variants,",
    "LATERAL unnest(duckvep_annotate(",
    "'r-test', 1::UINTEGER, position, reference, alternate, 0::UBIGINT",
    ")) u(a) ORDER BY ord"
  )
)
expect_equal(terminal_stop_edits$ord, 1:6)
expect_identical(
  terminal_stop_edits$consequence,
  c(
    "stop_lost",
    "stop_retained_variant",
    "stop_retained_variant",
    "stop_lost",
    "stop_lost",
    "coding_sequence_variant&3_prime_UTR_variant"
  )
)
expect_identical(terminal_stop_edits$status, rep("supported", 6))
expect_true(all(is.na(terminal_stop_edits$reason)))

feature_window_substitutions <- dbGetQuery(
  con,
  paste(
    "WITH variants(ord, position, reference, alternate) AS (VALUES",
    "(1, 122::UBIGINT, 'GG', 'GA'),",
    "(2, 123::UBIGINT, 'G', 'A'),",
    "(3, 237::UBIGINT, 'GT', 'TT'),",
    "(4, 237::UBIGINT, 'G', 'T'),",
    "(5, 237::UBIGINT, 'GT', 'AT'),",
    "(6, 237::UBIGINT, 'G', 'A'),",
    "(7, 122::UBIGINT, 'AG', 'AA'),",
    "(8, 122::UBIGINT, 'NG', 'NA'))",
    "SELECT ord, a.consequence, a.status, a.reason FROM variants,",
    "LATERAL unnest(duckvep_annotate(",
    "'r-test', 1::UINTEGER, position, reference, alternate, 0::UBIGINT",
    ")) u(a) ORDER BY ord"
  )
)
expect_equal(feature_window_substitutions$ord, 1:8)
expect_identical(
  feature_window_substitutions$consequence,
  c(
    "start_lost&start_retained_variant",
    "missense_variant",
    "stop_retained_variant",
    "missense_variant",
    "missense_variant",
    "stop_gained",
    "coding_sequence_variant",
    "coding_sequence_variant"
  )
)
expect_identical(
  feature_window_substitutions$status,
  c(rep("supported", 6), "unresolved", "unresolved")
)
expect_identical(
  feature_window_substitutions$reason[7:8],
  c("reference_mismatch", "ambiguous_sequence")
)

utr5_start_boundary <- dbGetQuery(
  con,
  paste(
    "WITH variants(ord, position, reference, alternate) AS (VALUES",
    "(1, 119::UBIGINT, 'GA', 'AC'),",
    "(2, 118::UBIGINT, 'CGA', 'GAA'),",
    "(3, 119::UBIGINT, 'GATGGT', 'TATGAT'),",
    "(4, 119::UBIGINT, 'GATG', 'GTAA'))",
    "SELECT ord, a.consequence, a.status FROM variants,",
    "LATERAL unnest(duckvep_annotate(",
    "'r-test', 1::UINTEGER, position, reference, alternate, 0::UBIGINT",
    ")) u(a) ORDER BY ord"
  )
)
expect_equal(utr5_start_boundary$ord, 1:4)
expect_identical(
  utr5_start_boundary$consequence,
  c(
    "start_lost&5_prime_UTR_variant",
    "start_retained_variant&5_prime_UTR_variant",
    "start_retained_variant&5_prime_UTR_variant",
    "start_lost&5_prime_UTR_variant"
  )
)
expect_identical(utr5_start_boundary$status, rep("supported", 4))

length_changing_boundaries <- dbGetQuery(
  con,
  paste(
    "WITH variants(ord, position, reference, alternate) AS (VALUES",
    "(1, 117::UBIGINT, 'ACGA', 'A'),",
    "(2, 118::UBIGINT, 'CGA', 'C'),",
    "(3, 118::UBIGINT, 'CGA', 'CACAC'),",
    "(4, 239::UBIGINT, 'AAA', 'A'),",
    "(5, 239::UBIGINT, 'AAAC', 'A'),",
    "(6, 240::UBIGINT, 'AA', 'C'))",
    "SELECT ord, a.consequence, a.status, a.reason FROM variants,",
    "LATERAL unnest(duckvep_annotate(",
    "'r-test', 1::UINTEGER, position, reference, alternate, 0::UBIGINT",
    ")) u(a) ORDER BY ord"
  )
)
expect_equal(length_changing_boundaries$ord, 1:6)
expect_identical(
  length_changing_boundaries$consequence,
  c(
    "start_lost&start_retained_variant&5_prime_UTR_variant",
    "start_lost&5_prime_UTR_variant",
    "start_lost&5_prime_UTR_variant",
    "stop_lost&3_prime_UTR_variant",
    "stop_retained_variant&3_prime_UTR_variant",
    "stop_lost&3_prime_UTR_variant"
  )
)
expect_identical(length_changing_boundaries$status, rep("supported", 6))
expect_true(all(is.na(length_changing_boundaries$reason)))

dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_transcripts_no_tail AS",
    "SELECT * EXCLUDE (pre_cds_sequence, post_cds_sequence)",
    "FROM duckvep_r_transcripts"
  )
)
no_tail_queries <- queries
no_tail_queries[2] <- sub(
  "duckvep_r_transcripts",
  "duckvep_r_transcripts_no_tail",
  no_tail_queries[2],
  fixed = TRUE
)
no_tail_queries[2] <- sub(
  ", pre_cds_sequence, post_cds_sequence",
  "",
  no_tail_queries[2],
  fixed = TRUE
)
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

dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_transcripts_legacy_tail AS",
    "SELECT * EXCLUDE (pre_cds_sequence, post_cds_sequence),",
    "'ACG'::BLOB AS post_cds_bases FROM duckvep_r_transcripts"
  )
)
legacy_queries <- queries
legacy_queries[2] <- sub(
  "duckvep_r_transcripts",
  "duckvep_r_transcripts_legacy_tail",
  legacy_queries[2],
  fixed = TRUE
)
legacy_queries[2] <- sub(
  "pre_cds_sequence, post_cds_sequence",
  "post_cds_bases",
  legacy_queries[2],
  fixed = TRUE
)
expect_true(load_model("r-legacy-tail", legacy_queries)$loaded)
legacy_boundary <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence, a.status, a.reason FROM unnest(duckvep_annotate(",
    "'r-legacy-tail', 1::UINTEGER, 117::UBIGINT, 'ACGA', 'A', 0::UBIGINT",
    ")) u(a)"
  )
)
expect_identical(
  legacy_boundary$consequence,
  "coding_sequence_variant&5_prime_UTR_variant"
)
expect_identical(legacy_boundary$status, "unresolved")
expect_identical(legacy_boundary$reason, "missing_transcript_flank")

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

dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_nmd_prediction_transcripts AS SELECT",
    "0::UINTEGER transcript_index, 1::UINTEGER seq_region,",
    "100::UBIGINT transcript_start, 498::UBIGINT transcript_end,",
    "1::TINYINT strand, 0::UINTEGER gene_index, 3::UBIGINT transcript_flags,",
    "100::UBIGINT cds_start, 498::UBIGINT cds_end,",
    "('ATG' || repeat('CAA', 81) || 'TAA')::BLOB cds_sequence,",
    "1::UTINYINT codon_table"
  )
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_nmd_prediction_exons AS SELECT * FROM (VALUES",
    "(0::UINTEGER, 100::UBIGINT, 199::UBIGINT, 1::UBIGINT, 100::UBIGINT, 0::TINYINT, 1::TINYINT),",
    "(0::UINTEGER, 300::UBIGINT, 399::UBIGINT, 101::UBIGINT, 200::UBIGINT, 1::TINYINT, 2::TINYINT),",
    "(0::UINTEGER, 450::UBIGINT, 498::UBIGINT, 201::UBIGINT, 249::UBIGINT, 2::TINYINT, 0::TINYINT)",
    ") t(transcript_index, exon_start, exon_end, exon_cdna_start,",
    "exon_cdna_end, phase, end_phase)"
  )
)
nmd_prediction_queries <- c(
  queries[1],
  paste(
    "SELECT * FROM duckvep_r_nmd_prediction_transcripts",
    "ORDER BY seq_region, transcript_start, transcript_index"
  ),
  paste(
    "SELECT * FROM duckvep_r_nmd_prediction_exons",
    "ORDER BY transcript_index, exon_cdna_start"
  )
)
expect_true(load_model("r-nmd-prediction", nmd_prediction_queries)$loaded)
nmd_prediction <- dbGetQuery(
  con,
  paste(
    "WITH variants(ord, position) AS (VALUES",
    "(1, 305::UBIGINT), (2, 350::UBIGINT))",
    "SELECT ord, a.nmd_prediction, a.nmd_escape_intronless,",
    "a.nmd_escape_early_cds, a.nmd_escape_last_exon,",
    "a.nmd_escape_penultimate_exon_end FROM variants,",
    "LATERAL unnest(duckvep_annotate(",
    "'r-nmd-prediction', 1::UINTEGER, position, 'C', 'T', 0::UBIGINT",
    ")) u(a) ORDER BY ord"
  )
)
expect_identical(nmd_prediction$nmd_prediction, c("triggering", "escaping"))
expect_identical(
  nmd_prediction$nmd_escape_penultimate_exon_end,
  c(FALSE, TRUE)
)
expect_true(
  all(!nmd_prediction$nmd_escape_intronless) &&
    all(!nmd_prediction$nmd_escape_early_cds) &&
    all(!nmd_prediction$nmd_escape_last_exon)
)

fixture_root <- system.file(
  "extdata",
  "duckvep",
  "ensembl_core",
  package = "Rduckhts"
)
expect_true(nzchar(fixture_root))
fixture_tables <- c(
  "attrib_type",
  "coord_system",
  "exon",
  "exon_transcript",
  "gene",
  "seq_region",
  "transcript",
  "transcript_attrib",
  "translation",
  "translation_attrib"
)

prepare_ensembl_fixture <- function(directory, assembly) {
  fixture_dir <- file.path(fixture_root, directory)
  schema <- paste0("duckvep_r_", directory, "_core")
  reference <- paste0("duckvep_r_", directory, "_reference")
  regions <- paste0("duckvep_r_", directory, "_regions")
  transcripts <- paste0("duckvep_r_", directory, "_transcripts")
  model <- paste0("r-ensembl-", directory)

  dbExecute(con, paste("CREATE SCHEMA", schema))
  for (table in fixture_tables) {
    path <- normalizePath(
      file.path(fixture_dir, paste0(table, ".parquet")),
      mustWork = TRUE
    )
    dbExecute(
      con,
      paste0(
        "CREATE VIEW ", schema, ".",
        as.character(dbQuoteIdentifier(con, table)),
        " AS SELECT * FROM read_parquet(",
        as.character(dbQuoteString(con, path)),
        ")"
      )
    )
  }
  reference_path <- normalizePath(
    file.path(fixture_dir, "reference_chunks.parquet"),
    mustWork = TRUE
  )
  dbExecute(
    con,
    paste0(
      "CREATE VIEW ", reference, " AS SELECT * FROM read_parquet(",
      as.character(dbQuoteString(con, reference_path)),
      ")"
    )
  )
  dbExecute(
    con,
    paste0(
      "CREATE TABLE ", regions,
      " AS SELECT * FROM duckvep_ensembl_regions(",
      as.character(dbQuoteString(con, schema)), ", ",
      as.character(dbQuoteString(con, reference)), ", ",
      as.character(dbQuoteString(con, assembly)), ")"
    )
  )
  dbExecute(
    con,
    paste0(
      "CREATE TABLE ", transcripts,
      " AS SELECT * FROM duckvep_ensembl_transcripts(",
      as.character(dbQuoteString(con, schema)), ", ",
      as.character(dbQuoteString(con, reference)), ", ",
      as.character(dbQuoteString(con, assembly)), ")"
    )
  )
  summary <- dbGetQuery(
    con,
    paste(
      "SELECT count(*) transcript_count,",
      "count(*) FILTER (WHERE cds_sequence IS NOT NULL) sequence_backed,",
      "count(*) FILTER (WHERE sequence_withheld_reason IS NOT NULL) withheld,",
      "sum(length(exons)) exon_memberships,",
      "count(*) FILTER (WHERE mane_select_refseq IS NOT NULL) mane_count",
      "FROM", transcripts
    )
  )
  model_queries <- c(
    paste(
      "SELECT seq_region, sequence_length FROM", regions,
      "ORDER BY seq_region"
    ),
    paste(
      "SELECT transcript_index, seq_region, transcript_start, transcript_end,",
      "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence,",
      "codon_table, pre_cds_sequence, post_cds_sequence FROM", transcripts,
      "ORDER BY seq_region, transcript_start, transcript_index"
    ),
    paste(
      "SELECT transcript_index, exon.exon_start, exon.exon_end,",
      "exon.exon_cdna_start, exon.exon_cdna_end, exon.phase, exon.end_phase",
      "FROM", transcripts, ", LATERAL unnest(exons) AS u(exon)",
      "ORDER BY transcript_index, exon.exon_cdna_start"
    )
  )
  expect_true(load_model(model, model_queries, TRUE)$loaded)
  list(
    model = model,
    transcripts = transcripts,
    summary = summary
  )
}

grch38_fixture <- prepare_ensembl_fixture("grch38", "GRCh38")
grch37_fixture <- prepare_ensembl_fixture("grch37", "GRCh37")
expect_equal(
  unlist(grch38_fixture$summary[1, ], use.names = FALSE),
  c(39, 1, 13, 40, 1)
)
expect_equal(
  unlist(grch37_fixture$summary[1, ], use.names = FALSE),
  c(39, 2, 13, 48, 0)
)

grch38_mane <- dbGetQuery(
  con,
  paste(
    "SELECT transcript_stable_id, mane_select_refseq FROM",
    grch38_fixture$transcripts,
    "WHERE mane_select_refseq IS NOT NULL"
  )
)
expect_identical(grch38_mane$transcript_stable_id, "ENST00000715685")
expect_identical(grch38_mane$mane_select_refseq, "NM_032790.4")

grch38_coding <- dbGetQuery(
  con,
  paste(
    "SELECT a.transcript_index, a.consequence, a.status, a.cds_position,",
    "a.protein_position, a.reference_amino_acid, a.alternate_amino_acid",
    "FROM unnest(duckvep_annotate(",
    as.character(dbQuoteString(con, grch38_fixture$model)), ",",
    "1::UINTEGER, 32522::UBIGINT, 'C', 'T', 0::UBIGINT)) u(a)"
  )
)
expect_equal(grch38_coding$transcript_index, 38)
expect_identical(grch38_coding$consequence, "missense_variant")
expect_identical(grch38_coding$status, "supported")
expect_equal(grch38_coding$cds_position, 4)
expect_equal(grch38_coding$protein_position, 2)
expect_identical(grch38_coding$reference_amino_acid, "H")
expect_identical(grch38_coding$alternate_amino_acid, "Y")

grch37_coding <- dbGetQuery(
  con,
  paste(
    "SELECT a.transcript_index, a.consequence, a.status, a.cds_position,",
    "a.protein_position, a.reference_amino_acid, a.alternate_amino_acid",
    "FROM unnest(duckvep_annotate(",
    as.character(dbQuoteString(con, grch37_fixture$model)), ",",
    "1::UINTEGER, 9452::UBIGINT, 'T', 'C', 0::UBIGINT)) u(a)"
  )
)
expect_equal(grch37_coding$transcript_index, 37)
expect_identical(grch37_coding$consequence, "missense_variant")
expect_identical(grch37_coding$status, "supported")
expect_equal(grch37_coding$cds_position, 4)
expect_equal(grch37_coding$protein_position, 2)
expect_identical(grch37_coding$reference_amino_acid, "Y")
expect_identical(grch37_coding$alternate_amino_acid, "H")

fixture_mitochondrial <- dbGetQuery(
  con,
  paste(
    "SELECT a.transcript_index, a.consequence, a.status, a.reason",
    "FROM unnest(duckvep_annotate(",
    as.character(dbQuoteString(con, grch38_fixture$model)), ",",
    "0::UINTEGER, 3307::UBIGINT,",
    "'A', 'C', 0::UBIGINT)) u(a)"
  )
)
expect_equal(fixture_mitochondrial$transcript_index, 5)
expect_identical(
  fixture_mitochondrial$consequence,
  "coding_sequence_variant"
)
expect_identical(fixture_mitochondrial$status, "unresolved")
expect_identical(fixture_mitochondrial$reason, "missing_sequence")

expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('r-test') AS dropped")$dropped)
expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('r-no-tail') AS dropped")$dropped)
expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('r-legacy-tail') AS dropped")$dropped)
expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('r-nmd') AS dropped")$dropped)
expect_true(
  dbGetQuery(
    con,
    "SELECT duckvep_model_drop('r-nmd-prediction') AS dropped"
  )$dropped
)
expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('r-complete') AS dropped")$dropped)
expect_true(
  dbGetQuery(
    con,
    "SELECT duckvep_model_drop('r-partial-cds-end') AS dropped"
  )$dropped
)
expect_true(
  dbGetQuery(
    con,
    "SELECT duckvep_model_drop('r-ensembl-grch38') AS dropped"
  )$dropped
)
expect_true(
  dbGetQuery(
    con,
    "SELECT duckvep_model_drop('r-ensembl-grch37') AS dropped"
  )$dropped
)

dbDisconnect(con, shutdown = TRUE)
