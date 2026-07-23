library(tinytest)
library(DBI)

con <- dbConnect(
  duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
)
expect_silent(rduckhts_load(con))

geometry <- dbGetQuery(
  con,
  paste(
    "SELECT g.* FROM (SELECT",
    "duckvep_allele_geometry(100, 'a', 'atg') AS g)"
  )
)
expect_identical(geometry$kind_code, 1L)
expect_identical(geometry$interbase, TRUE)
expect_equal(geometry$raw_start0, 99)
expect_equal(geometry$raw_end0, 100)
expect_equal(geometry$feature_start0, 100)
expect_equal(geometry$feature_end0, 100)
expect_equal(geometry$edit_start0, 100)
expect_equal(geometry$edit_end0, 100)
expect_equal(geometry$insertion_boundary0, 100)
expect_identical(geometry$reference_difference_length, 0L)
expect_identical(geometry$alternate_difference_length, 2L)

right_anchor <- dbGetQuery(
  con,
  paste(
    "SELECT g.* FROM (SELECT",
    "duckvep_allele_geometry(1, 'A', 'CA') AS g)"
  )
)
expect_identical(right_anchor$kind_code, 1L)
expect_identical(right_anchor$anchor_side_code, 2L)
expect_equal(right_anchor$insertion_boundary0, 0)

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
    "CREATE TABLE duckvep_r_core.seq_region_attrib AS SELECT",
    "NULL::BIGINT AS seq_region_id, NULL::BIGINT AS attrib_type_id,",
    "NULL::VARCHAR AS value WHERE false"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.gene AS SELECT 1::BIGINT AS gene_id,",
    "'protein_coding'::VARCHAR AS biotype, 1::BIGINT AS is_current,",
    "'ENSG_R_TEST'::VARCHAR AS stable_id, 1::BIGINT AS version",
    "UNION ALL SELECT 2, 'miRNA', 1, 'ENSG_R_MIRNA', 1"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.transcript AS SELECT",
    "1::BIGINT AS transcript_id, 1::BIGINT AS gene_id,",
    "1::BIGINT AS seq_region_id, 1::BIGINT AS seq_region_start,",
    "12::BIGINT AS seq_region_end, 1::BIGINT AS seq_region_strand,",
    "'protein_coding'::VARCHAR AS biotype, 1::BIGINT AS is_current,",
    "'ENST_R_TEST'::VARCHAR AS stable_id, 1::BIGINT AS version",
    "UNION ALL SELECT 2, 2, 1, 1, 12, 1, 'miRNA', 1,",
    "'ENST_R_MIRNA', 1"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.exon AS SELECT 1::BIGINT AS exon_id,",
    "1::BIGINT AS seq_region_id, 1::BIGINT AS seq_region_start,",
    "12::BIGINT AS seq_region_end, 1::BIGINT AS seq_region_strand,",
    "-1::BIGINT AS phase, 0::BIGINT AS end_phase,",
    "1::BIGINT AS is_current, 'ENSE_R_TEST'::VARCHAR AS stable_id",
    "UNION ALL SELECT 2, 1, 1, 12, 1, -1, -1, 1, 'ENSE_R_MIRNA'"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.exon_transcript AS SELECT",
    "1::BIGINT AS exon_id, 1::BIGINT AS transcript_id, 1::BIGINT AS rank",
    "UNION ALL SELECT 2, 2, 1"
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
    "1::BIGINT AS attrib_type_id, 'MANE_Select'::VARCHAR AS code",
    "UNION ALL SELECT 2, 'miRNA'"
  ),
  paste(
    "CREATE TABLE duckvep_r_core.transcript_attrib AS SELECT",
    "1::BIGINT AS transcript_id, 1::BIGINT AS attrib_type_id,",
    "'NM_R_TEST.1'::VARCHAR AS value",
    "UNION ALL SELECT 2, 2, '2-5'"
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
funcgen_fixture_sql <- c(
  "CREATE SCHEMA duckvep_r_funcgen",
  paste(
    "CREATE TABLE duckvep_r_funcgen.feature_type AS SELECT",
    "70::BIGINT feature_type_id, 'Promoter'::VARCHAR \"name\",",
    "'SO:0000167'::VARCHAR so_accession, 'promoter'::VARCHAR so_term"
  ),
  paste(
    "CREATE TABLE duckvep_r_funcgen.regulatory_feature AS SELECT",
    "80::BIGINT regulatory_feature_id, 70::BIGINT feature_type_id,",
    "1::BIGINT seq_region_id, 0::BIGINT seq_region_strand,",
    "2::BIGINT seq_region_start, 6::BIGINT seq_region_end,",
    "'ENSR_R_TEST'::VARCHAR stable_id, 1::BIGINT regulatory_build_id"
  ),
  paste(
    "CREATE TABLE duckvep_r_funcgen.motif_feature AS SELECT",
    "90::BIGINT motif_feature_id, 100::BIGINT binding_matrix_id,",
    "1::BIGINT seq_region_id, 7::BIGINT seq_region_start,",
    "10::BIGINT seq_region_end, 1::BIGINT seq_region_strand,",
    "8.5::DOUBLE score, 'ENSM_R_TEST'::VARCHAR stable_id"
  )
)
for (sql in funcgen_fixture_sql) {
  dbExecute(con, sql)
}
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_ensembl_regulation AS SELECT * FROM",
    "duckvep_ensembl_regulation_features(",
    "'duckvep_r_funcgen', 'duckvep_r_ensembl_regions')"
  )
)
regulation_features <- dbGetQuery(
  con,
  paste(
    "SELECT feature_class, source_feature_table, stable_id",
    "FROM duckvep_r_ensembl_regulation",
    "ORDER BY regulation_feature_index"
  )
)
expect_identical(
  regulation_features$feature_class,
  c("regulatory_region", "transcription_factor_binding_site")
)
expect_identical(
  regulation_features$source_feature_table,
  c("regulatory_feature", "motif_feature")
)
expect_identical(
  regulation_features$stable_id,
  c("ENSR_R_TEST", "ENSM_R_TEST")
)
ensembl_model <- dbGetQuery(
  con,
  paste(
    "SELECT CAST(cds_sequence AS VARCHAR) cds_sequence,",
    "CAST(post_cds_bases AS VARCHAR) post_cds_bases,",
    "CAST(pre_cds_sequence AS VARCHAR) pre_cds_sequence,",
    "CAST(post_cds_sequence AS VARCHAR) post_cds_sequence,",
    "transcript_flags, mane_select_refseq, exons[1].phase phase",
    "FROM duckvep_r_ensembl_transcripts",
    "WHERE transcript_stable_id = 'ENST_R_TEST'"
  )
)
expect_identical(ensembl_model$cds_sequence, "ATGAAATAA")
expect_identical(ensembl_model$post_cds_bases, "CCC")
expect_identical(ensembl_model$pre_cds_sequence, "")
expect_identical(ensembl_model$post_cds_sequence, "CCC")
expect_equal(ensembl_model$transcript_flags, 1027)
expect_identical(ensembl_model$mane_select_refseq, "NM_R_TEST.1")
expect_equal(ensembl_model$phase, -1)
mature_mirna_region <- dbGetQuery(
  con,
  paste(
    "SELECT transcript_index, region.cdna_start, region.cdna_end,",
    "region.mature_mirna_start, region.mature_mirna_end",
    "FROM duckvep_r_ensembl_transcripts,",
    "LATERAL unnest(mature_mirna_regions) AS u(region)"
  )
)
expect_equal(mature_mirna_region$transcript_index, 1)
expect_equal(
  unlist(mature_mirna_region[1, -1], use.names = FALSE),
  c(2, 5, 2, 5)
)
ensembl_receipt <- dbGetQuery(
  con,
  paste(
    "SELECT * FROM duckvep_model_receipt(",
    "'duckvep_r_ensembl_regions', 'duckvep_r_ensembl_transcripts',",
    "'Ensembl', '116', 'GRCh38', repeat('a', 64), repeat('b', 64),",
    "'all current transcripts on FASTA-covered assembly regions',",
    "regulation_features_table := 'duckvep_r_ensembl_regulation')"
  )
)
expect_equal(ensembl_receipt$region_count, 1)
expect_equal(ensembl_receipt$reference_base_count, 12)
expect_equal(ensembl_receipt$transcript_count, 2)
expect_equal(ensembl_receipt$mature_mirna_transcript_count, 1)
expect_equal(ensembl_receipt$mature_mirna_segment_count, 1)
expect_equal(ensembl_receipt$transcript_flank_base_count, 3)
expect_equal(ensembl_receipt$regulation_feature_count, 2)
expect_equal(ensembl_receipt$regulatory_region_count, 1)
expect_equal(ensembl_receipt$motif_feature_count, 1)
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
  paste(
    "INSERT INTO duckvep_r_core.transcript VALUES",
    "(3, 1, 1, 1, 12, 1, 'protein_coding', 1, NULL, 1),",
    "(4, 1, 1, 1, 12, 1, 'artifact', 1, 'ENST_R_ARTIFACT', 1),",
    "(5, 1, 1, 1, 12, 1, 'protein_coding', 1, 'ENST_R_READTHROUGH', 1)"
  )
)
dbExecute(
  con,
  "INSERT INTO duckvep_r_core.attrib_type VALUES (3, 'readthrough_tra')"
)
dbExecute(
  con,
  "INSERT INTO duckvep_r_core.transcript_attrib VALUES (5, 3, '1')"
)
eligible_transcripts <- dbGetQuery(
  con,
  paste(
    "SELECT transcript_stable_id FROM duckvep_ensembl_transcripts(",
    "'duckvep_r_core', 'duckvep_r_reference', 'GRCh38')",
    "ORDER BY transcript_index"
  )
)
expect_identical(
  eligible_transcripts$transcript_stable_id,
  c("ENST_R_TEST", "ENST_R_MIRNA")
)

dbExecute(
  con,
  "INSERT INTO duckvep_r_core.attrib_type VALUES (4, '_rna_edit')"
)
dbExecute(
  con,
  "INSERT INTO duckvep_r_core.transcript_attrib VALUES (1, 4, '4 5 A')"
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
    "sequence_withheld_reason FROM duckvep_r_ensembl_rna_edit",
    "WHERE transcript_stable_id = 'ENST_R_TEST'"
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
load_model <- function(
  name,
  model_queries,
  transcript_coverage_complete = FALSE,
  mature_mirna_query = NULL,
  peptide_edit_query = NULL,
  interval_feature_query = NULL,
  reference_fasta = NULL
) {
  arguments <- paste(
    vapply(
      c(name, model_queries),
      function(x) as.character(dbQuoteString(con, x)),
      character(1)
    ),
    collapse = ", "
  )
  if (!is.null(mature_mirna_query)) {
    arguments <- paste(
      arguments,
      paste0(
        "mature_mirna_query := ",
        as.character(dbQuoteString(con, mature_mirna_query))
      ),
      sep = ", "
    )
  }
  if (!is.null(peptide_edit_query)) {
    arguments <- paste(
      arguments,
      paste0(
        "peptide_edit_query := ",
        as.character(dbQuoteString(con, peptide_edit_query))
      ),
      sep = ", "
    )
  }
  if (!is.null(interval_feature_query)) {
    arguments <- paste(
      arguments,
      paste0(
        "interval_feature_query := ",
        as.character(dbQuoteString(con, interval_feature_query))
      ),
      sep = ", "
    )
  }
  if (!is.null(reference_fasta)) {
    arguments <- paste(
      arguments,
      paste0(
        "reference_fasta := ",
        as.character(dbQuoteString(con, reference_fasta))
      ),
      sep = ", "
    )
  }
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

hgvs_reference <- system.file(
  "extdata",
  "fixture_ref.fa",
  package = "Rduckhts",
  mustWork = TRUE
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_hgvs_regions AS SELECT",
    "1::UINTEGER seq_region, 50000::UBIGINT sequence_length,",
    "'11'::VARCHAR seq_region_name"
  )
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_hgvs_transcripts AS SELECT",
    "0::UINTEGER transcript_index, 1::UINTEGER seq_region,",
    "100::UBIGINT transcript_start, 150::UBIGINT transcript_end,",
    "1::TINYINT strand, 0::UINTEGER gene_index,",
    "0::UBIGINT transcript_flags, NULL::UBIGINT cds_start,",
    "NULL::UBIGINT cds_end, NULL::BLOB cds_sequence,",
    "NULL::UTINYINT codon_table, NULL::BLOB pre_cds_sequence,",
    "NULL::BLOB post_cds_sequence"
  )
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_hgvs_exons AS SELECT",
    "0::UINTEGER transcript_index, 100::UBIGINT exon_start,",
    "150::UBIGINT exon_end, 1::UBIGINT exon_cdna_start,",
    "51::UBIGINT exon_cdna_end, -1::TINYINT phase,",
    "-1::TINYINT end_phase"
  )
)
expect_true(
  load_model(
    "r-hgvs",
    c(
      "SELECT * FROM duckvep_r_hgvs_regions ORDER BY seq_region",
      paste(
        "SELECT * FROM duckvep_r_hgvs_transcripts",
        "ORDER BY seq_region, transcript_start, transcript_index"
      ),
      paste(
        "SELECT * FROM duckvep_r_hgvs_exons",
        "ORDER BY transcript_index, exon_cdna_start"
      )
    ),
    reference_fasta = hgvs_reference
  )$loaded
)
hgvs_annotation <- dbGetQuery(
  con,
  paste(
    "WITH variants(ord, position, reference, alternate) AS (VALUES",
    "(1, 124::UBIGINT, 'A', 'G'),",
    "(2, 124::UBIGINT, 'A', 'AC'))",
    "SELECT ord, a.transcript_hgvs, a.hgvs_shift,",
    "a.transcript_hgvs_status, a.protein_hgvs_status",
    "FROM variants, LATERAL unnest(_duckvep_annotate_small_hgvs(",
    "'r-hgvs', 1::UINTEGER, position, reference, alternate, 0::UBIGINT",
    ")) AS u(a) ORDER BY ord"
  )
)
expect_identical(hgvs_annotation$transcript_hgvs, c("n.25A>G", "n.25_26insC"))
expect_equal(hgvs_annotation$hgvs_shift, c(0, 0))
expect_identical(
  hgvs_annotation$transcript_hgvs_status,
  c("supported", "supported")
)
expect_identical(
  hgvs_annotation$protein_hgvs_status,
  c("not_applicable", "not_applicable")
)

# The public relation macro owns event-family dispatch and optional HGVS. Keep
# one DBI regression on that surface rather than proving only its private lanes.
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_public_events AS SELECT * FROM (VALUES",
    "(1::UBIGINT, 1::UINTEGER, 124::UBIGINT, 'A'::VARCHAR, 'G'::VARCHAR,",
    "NULL::UBIGINT, NULL::VARCHAR, NULL::VARCHAR,",
    "NULL::UINTEGER, NULL::UBIGINT))",
    "e(event_index, seq_region, position, reference, alternate, end_position,",
    "structural_type, copy_change, mate_seq_region, mate_position)"
  )
)
public_annotation <- dbGetQuery(
  con,
  paste(
    "SELECT a.event_index, a.transcript_index,",
    "string_agg(t.consequence, '&' ORDER BY t.severity_rank) consequence,",
    "a.transcript_hgvs, a.protein_hgvs",
    "FROM duckvep_annotate('duckvep_r_public_events', 'r-hgvs',",
    "hgvs := true, upstream_distance := 0, downstream_distance := 0) a",
    "JOIN duckvep_so_terms() t",
    "ON (a.consequence_mask & t.consequence_mask) <> 0",
    "GROUP BY ALL ORDER BY a.event_index, a.transcript_index"
  )
)
expect_equal(public_annotation$event_index, 1)
expect_equal(public_annotation$transcript_index, 0)
expect_identical(
  public_annotation$consequence,
  "non_coding_transcript_exon_variant"
)
expect_identical(public_annotation$transcript_hgvs, "n.25A>G")
expect_true(is.na(public_annotation$protein_hgvs))

public_rich_annotation <- dbGetQuery(
  con,
  paste(
    "SELECT a.event_index, a.transcript_index, a.consequence, a.impact,",
    "a.protein_position, a.nmd_prediction, a.duckvep_status,",
    "a.transcript_hgvs, a.protein_hgvs",
    "FROM duckvep_annotate('duckvep_r_public_events', 'r-hgvs',",
    "hgvs := TRUE, upstream_distance := 0, downstream_distance := 0,",
    "rich := TRUE) a ORDER BY a.event_index, a.transcript_index"
  )
)
expect_equal(public_rich_annotation$event_index, 1)
expect_equal(public_rich_annotation$transcript_index, 0)
expect_identical(
  public_rich_annotation$consequence,
  "non_coding_transcript_exon_variant"
)
expect_identical(public_rich_annotation$impact, "MODIFIER")
expect_true(is.na(public_rich_annotation$protein_position))
expect_true(is.na(public_rich_annotation$nmd_prediction))
expect_identical(public_rich_annotation$duckvep_status, "supported")
expect_identical(public_rich_annotation$transcript_hgvs, "n.25A>G")
expect_true(is.na(public_rich_annotation$protein_hgvs))

# Exercise the rare terminal-CDS states through the public event relation, not
# only through the kernel and private scalar fixtures. The start-deletion model
# is synthetic; the terminal-partial model below uses exact remapped Ensembl-116
# ENST00000650713 exon/CDS/reference data.
partial_reference <- system.file(
  "extdata",
  "duckvep_terminal_partial.fa",
  package = "Rduckhts",
  mustWork = TRUE
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_public_partial_regions AS SELECT",
    "1::UINTEGER seq_region, 54::UBIGINT sequence_length,",
    "'partial'::VARCHAR seq_region_name"
  )
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_public_partial_transcripts AS SELECT",
    "0::UINTEGER transcript_index, 1::UINTEGER seq_region,",
    "21::UBIGINT transcript_start, 34::UBIGINT transcript_end,",
    "1::TINYINT strand, 0::UINTEGER gene_index,",
    "35::UBIGINT transcript_flags, 21::UBIGINT cds_start,",
    "34::UBIGINT cds_end, 'ATGCCCAAAGGGTC'::BLOB cds_sequence,",
    "1::UTINYINT codon_table, ''::BLOB pre_cds_sequence,",
    "''::BLOB post_cds_sequence"
  )
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_public_partial_exons AS SELECT",
    "0::UINTEGER transcript_index, 21::UBIGINT exon_start,",
    "34::UBIGINT exon_end, 1::UBIGINT exon_cdna_start,",
    "14::UBIGINT exon_cdna_end, 0::TINYINT phase,",
    "2::TINYINT end_phase"
  )
)
expect_true(
  load_model(
    "r-public-terminal-partial",
    c(
      "SELECT * FROM duckvep_r_public_partial_regions ORDER BY seq_region",
      paste(
        "SELECT * FROM duckvep_r_public_partial_transcripts",
        "ORDER BY seq_region, transcript_start, transcript_index"
      ),
      paste(
        "SELECT * FROM duckvep_r_public_partial_exons",
        "ORDER BY transcript_index, exon_cdna_start"
      )
    ),
    reference_fasta = partial_reference
  )$loaded
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_public_partial_events AS SELECT * FROM (VALUES",
    "(1::UBIGINT, 1::UINTEGER, 21::UBIGINT, 'ATGC'::VARCHAR, 'C'::VARCHAR,",
    "NULL::UBIGINT, NULL::VARCHAR, NULL::VARCHAR,",
    "NULL::UINTEGER, NULL::UBIGINT))",
    "e(event_index, seq_region, position, reference, alternate, end_position,",
    "structural_type, copy_change, mate_seq_region, mate_position)"
  )
)

clinvar_partial_reference <- system.file(
  "extdata",
  "duckvep_clinvar_terminal_partial.fa",
  package = "Rduckhts",
  mustWork = TRUE
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_clinvar_partial_regions AS SELECT",
    "1::UINTEGER seq_region, 2788::UBIGINT sequence_length,",
    "'clinvar_terminal_partial'::VARCHAR seq_region_name"
  )
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_clinvar_partial_transcripts AS SELECT",
    "0::UINTEGER transcript_index, 1::UINTEGER seq_region,",
    "1::UBIGINT transcript_start, 1788::UBIGINT transcript_end,",
    "1::TINYINT strand, 0::UINTEGER gene_index,",
    "35::UBIGINT transcript_flags, 1077::UBIGINT cds_start,",
    "1788::UBIGINT cds_end,",
    paste0(
      "'ATGCGCCAGGCAGGCCGTTACTTACCAGAGTTTAGGGAAACCCGGGCTGCCCAGGACTTTTTC",
      "AGCACGTGTCGCTCTCCTGAGGCCTGCTGTGAACTGACTCTGCAGCCACTGCGTCGCTTCCCT",
      "CTGGATGCTGCCATCATTTTCTCCGACATCCTTGTTGTACCCCAGGCACTGGGCATGGAGGTG",
      "ACCATGGTACCTGGCAAAGGACCCAGCTTCCCAGAGCCATTAAGAGAAGAGCAGGACCTAGAAC",
      "GCCTACGGGATCCAGAAGTGGTAGCCTC'::BLOB cds_sequence,"
    ),
    "1::UTINYINT codon_table,",
    paste0(
      "'GGCTGGATAAGACTGTTGGTATCATGAGTGGGACTTGCGCCAAGCCTCCGGATACCCAGACTGT",
      "CAGATGAGAACAAATTCCTCATGTCACCGTAAGATACATTTACAGCGGAGTTTTCTTTTGGGCC",
      "TTTGTTGTTGCGTCGCTACAGCAAACTTTACGGTGAAAAAAGACCTCAGGGTTTTCCGGAGCTG",
      "AAGAATGACACATTCCTGCGAGCAGCCTGGGGAGAGGAAACAGACTACACTCCCGTTTGGTGC'",
      "::BLOB pre_cds_sequence,"
    ),
    "''::BLOB post_cds_sequence"
  )
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_clinvar_partial_exons AS SELECT * FROM (VALUES",
    "(0::UINTEGER,1::UBIGINT,170::UBIGINT,1::UBIGINT,170::UBIGINT,-1::TINYINT,-1::TINYINT),",
    "(0::UINTEGER,992::UBIGINT,1104::UBIGINT,171::UBIGINT,283::UBIGINT,-1::TINYINT,1::TINYINT),",
    "(0::UINTEGER,1221::UBIGINT,1300::UBIGINT,284::UBIGINT,363::UBIGINT,1::TINYINT,0::TINYINT),",
    "(0::UINTEGER,1377::UBIGINT,1439::UBIGINT,364::UBIGINT,426::UBIGINT,0::TINYINT,0::TINYINT),",
    "(0::UINTEGER,1679::UBIGINT,1788::UBIGINT,427::UBIGINT,536::UBIGINT,0::TINYINT,2::TINYINT))",
    "e(transcript_index,exon_start,exon_end,exon_cdna_start,exon_cdna_end,phase,end_phase)"
  )
)
expect_true(
  load_model(
    "r-clinvar-terminal-partial",
    c(
      "SELECT * FROM duckvep_r_clinvar_partial_regions ORDER BY seq_region",
      paste(
        "SELECT * FROM duckvep_r_clinvar_partial_transcripts",
        "ORDER BY seq_region, transcript_start, transcript_index"
      ),
      paste(
        "SELECT * FROM duckvep_r_clinvar_partial_exons",
        "ORDER BY transcript_index, exon_cdna_start"
      )
    ),
    reference_fasta = clinvar_partial_reference
  )$loaded
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_clinvar_partial_events AS SELECT * FROM (VALUES",
    "(2::UBIGINT,1::UINTEGER,1786::UBIGINT,'C'::VARCHAR,'CTAG'::VARCHAR,",
    "NULL::UBIGINT,NULL::VARCHAR,NULL::VARCHAR,NULL::UINTEGER,NULL::UBIGINT))",
    "e(event_index,seq_region,position,reference,alternate,end_position,",
    "structural_type,copy_change,mate_seq_region,mate_position)"
  )
)
public_partial <- dbGetQuery(
  con,
  paste(
    "WITH annotations AS (",
    "SELECT * FROM duckvep_annotate('duckvep_r_public_partial_events',",
    "'r-public-terminal-partial', hgvs := TRUE, rich := TRUE,",
    "upstream_distance := 0, downstream_distance := 0)",
    "UNION ALL BY NAME",
    "SELECT * FROM duckvep_annotate('duckvep_r_clinvar_partial_events',",
    "'r-clinvar-terminal-partial', hgvs := TRUE, rich := TRUE,",
    "upstream_distance := 0, downstream_distance := 0))",
    "SELECT event_index, consequence, transcript_hgvs, protein_hgvs,",
    "transcript_hgvs_status, protein_hgvs_status, hgvs_shift",
    "FROM annotations",
    "ORDER BY event_index"
  )
)
expect_equal(public_partial$event_index, c(1, 2))
expect_identical(
  public_partial$consequence,
  c(
    "inframe_deletion",
    "stop_gained&inframe_insertion&incomplete_terminal_codon_variant"
  )
)
expect_false(grepl("start_lost", public_partial$consequence[[1]], fixed = TRUE))
expect_identical(public_partial$transcript_hgvs_status, c("supported", "supported"))
expect_identical(public_partial$protein_hgvs_status, c("supported", "supported"))
expect_identical(
  public_partial$transcript_hgvs,
  c("c.2_4del", "c.280_281insAGT")
)
expect_identical(
  public_partial$protein_hgvs,
  c("p.Met1_Pro2delinsThr", "p.Ter94=")
)
expect_equal(public_partial$hgvs_shift, c(1, 1))

public_annotation_without_hgvs <- dbGetQuery(
  con,
  paste(
    "SELECT a.event_index, a.transcript_hgvs IS NULL AS hgvs_is_null",
    "FROM duckvep_annotate('duckvep_r_public_events', 'r-hgvs',",
    "hgvs := NULL, upstream_distance := 0, downstream_distance := 0) a"
  )
)
expect_equal(public_annotation_without_hgvs$event_index, 1)
expect_true(public_annotation_without_hgvs$hgvs_is_null)

dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_contradictory_event AS SELECT * FROM (VALUES",
    "(2::UBIGINT, 1::UINTEGER, 90::UBIGINT, NULL::VARCHAR, '<DEL>'::VARCHAR,",
    "260::UBIGINT, 'DUP'::VARCHAR, NULL::VARCHAR,",
    "NULL::UINTEGER, NULL::UBIGINT))",
    "e(event_index, seq_region, position, reference, alternate, end_position,",
    "structural_type, copy_change, mate_seq_region, mate_position)"
  )
)
expect_error(
  dbGetQuery(
    con,
    "SELECT * FROM duckvep_annotate('duckvep_r_contradictory_event', 'r-hgvs')"
  ),
  pattern = "symbolic alternate and structural_type disagree"
)

# Output beyond the initial native HGVS scratch capacity must be retried at
# the reported size rather than exposing a truncated or invalid DuckDB string.
hgvs_long <- dbGetQuery(
  con,
  paste(
    "SELECT a.transcript_hgvs, a.transcript_hgvs_status",
    "FROM unnest(_duckvep_annotate_small_hgvs(",
    "'r-hgvs', 1::UINTEGER, 124::UBIGINT,",
    "'A', 'C' || repeat('A', 1405), 0::UBIGINT)) AS u(a)"
  )
)
expect_identical(hgvs_long$transcript_hgvs_status, "supported")
expect_true(nchar(hgvs_long$transcript_hgvs, type = "bytes") > 1400L)

# A BGZF reference exercises the retained .gzi descriptor and worker-local
# compressed-reference handle through the public DBI surface.
compressed_reference <- tempfile("duckvep-reference-", fileext = ".fa.gz")
compressed_index <- paste0(compressed_reference, ".fai")
compressed_gzi <- paste0(compressed_reference, ".gzi")
expect_true(
  rduckhts_bgzip(
    con,
    hgvs_reference,
    output_path = compressed_reference,
    threads = 1L,
    keep = TRUE,
    overwrite = TRUE
  )$success
)
expect_true(
  rduckhts_fasta_index(
    con,
    compressed_reference,
    index_path = compressed_index
  )$success
)
expect_true(file.exists(compressed_gzi))
expect_true(
  load_model(
    "r-hgvs-bgzf",
    c(
      "SELECT * FROM duckvep_r_hgvs_regions ORDER BY seq_region",
      paste(
        "SELECT * FROM duckvep_r_hgvs_transcripts",
        "ORDER BY seq_region, transcript_start, transcript_index"
      ),
      paste(
        "SELECT * FROM duckvep_r_hgvs_exons",
        "ORDER BY transcript_index, exon_cdna_start"
      )
    ),
    reference_fasta = compressed_reference
  )$loaded
)
compressed_annotation <- dbGetQuery(
  con,
  paste(
    "SELECT a.transcript_hgvs, a.transcript_hgvs_status",
    "FROM unnest(_duckvep_annotate_small_hgvs(",
    "'r-hgvs-bgzf', 1::UINTEGER, 124::UBIGINT,",
    "'A', 'G', 0::UBIGINT)) u(a)"
  )
)
expect_identical(compressed_annotation$transcript_hgvs, "n.25A>G")
expect_identical(compressed_annotation$transcript_hgvs_status, "supported")

# A named model pins the reference files validated at load. Worker-local faidx
# handles must not silently reopen replacement content.
identity_reference <- tempfile("duckvep-reference-", fileext = ".fa")
identity_index <- paste0(identity_reference, ".fai")
expect_true(file.copy(hgvs_reference, identity_reference))
expect_true(file.copy(paste0(hgvs_reference, ".fai"), identity_index))
expect_true(
  load_model(
    "r-hgvs-reference-identity",
    c(
      "SELECT * FROM duckvep_r_hgvs_regions ORDER BY seq_region",
      paste(
        "SELECT * FROM duckvep_r_hgvs_transcripts",
        "ORDER BY seq_region, transcript_start, transcript_index"
      ),
      paste(
        "SELECT * FROM duckvep_r_hgvs_exons",
        "ORDER BY transcript_index, exon_cdna_start"
      )
    ),
    reference_fasta = identity_reference
  )$loaded
)
identity_connection <- try(
  file(identity_reference, open = "r+b"),
  silent = TRUE
)
identity_query <- paste(
  "SELECT a.transcript_hgvs_status FROM unnest(_duckvep_annotate_small_hgvs(",
  "'r-hgvs-reference-identity', 1::UINTEGER, 124::UBIGINT,",
  "'A', 'G', 0::UBIGINT)) u(a)"
)
if (inherits(identity_connection, "try-error")) {
  # Windows retains deny-write handles and the still-pinned source remains usable.
  expect_identical(.Platform$OS.type, "windows")
  expect_identical(
    dbGetQuery(con, identity_query)$transcript_hgvs_status,
    "supported"
  )
} else {
  seek(identity_connection, where = 4L, origin = "start")
  writeBin(charToRaw("C"), identity_connection)
  close(identity_connection)
  expect_error(
    dbGetQuery(con, identity_query),
    pattern = "reference FASTA or index changed after model load"
  )
}

# Linux path replacement cannot redirect a later worker open. Windows keeps a
# deny-write/delete handle on the final resolved path. Other POSIX systems use
# an independently opened resolved path and reject an observed replacement.
replacement_reference <- tempfile("duckvep-reference-replacement-", fileext = ".fa")
replacement_index <- paste0(replacement_reference, ".fai")
replacement_reference_moved <- paste0(replacement_reference, ".loaded")
replacement_index_moved <- paste0(replacement_index, ".loaded")
expect_true(file.copy(hgvs_reference, replacement_reference))
expect_true(file.copy(paste0(hgvs_reference, ".fai"), replacement_index))
expect_true(
  load_model(
    "r-hgvs-reference-replacement",
    c(
      "SELECT * FROM duckvep_r_hgvs_regions ORDER BY seq_region",
      paste(
        "SELECT * FROM duckvep_r_hgvs_transcripts",
        "ORDER BY seq_region, transcript_start, transcript_index"
      ),
      paste(
        "SELECT * FROM duckvep_r_hgvs_exons",
        "ORDER BY transcript_index, exon_cdna_start"
      )
    ),
    reference_fasta = replacement_reference
  )$loaded
)
replacement_moved <- file.rename(
  replacement_reference,
  replacement_reference_moved
)
if (.Platform$OS.type == "windows") {
  expect_false(replacement_moved)
} else {
  expect_true(replacement_moved)
  expect_true(file.rename(replacement_index, replacement_index_moved))
  expect_true(file.copy(replacement_reference_moved, replacement_reference))
  expect_true(file.copy(replacement_index_moved, replacement_index))
  replacement_connection <- file(replacement_reference, open = "r+b")
  # fixture_ref.fa position 124 is byte offset 128 after its header/newline.
  seek(replacement_connection, where = 128L, origin = "start")
  writeBin(charToRaw("C"), replacement_connection)
  close(replacement_connection)
  replacement_query <- paste(
    "SELECT a.transcript_hgvs, a.transcript_hgvs_status",
    "FROM unnest(_duckvep_annotate_small_hgvs(",
    "'r-hgvs-reference-replacement', 1::UINTEGER, 124::UBIGINT,",
    "'A', 'G', 0::UBIGINT)) u(a)"
  )
  if (identical(Sys.info()[["sysname"]], "Linux")) {
    replacement_annotation <- dbGetQuery(con, replacement_query)
    expect_identical(replacement_annotation$transcript_hgvs, "n.25A>G")
    expect_identical(
      replacement_annotation$transcript_hgvs_status,
      "supported"
    )
  } else {
    expect_error(
      dbGetQuery(con, replacement_query),
      pattern = "reference FASTA or index changed after model load"
    )
  }
}

dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_hgvs_coding_transcripts AS SELECT",
    "0::UINTEGER transcript_index, 1::UINTEGER seq_region,",
    "100::UBIGINT transcript_start, 150::UBIGINT transcript_end,",
    "1::TINYINT strand, 0::UINTEGER gene_index,",
    "3::UBIGINT transcript_flags, 100::UBIGINT cds_start,",
    "150::UBIGINT cds_end, repeat('A', 51)::BLOB cds_sequence,",
    "1::UTINYINT codon_table, ''::BLOB pre_cds_sequence,",
    "''::BLOB post_cds_sequence"
  )
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_hgvs_coding_exons AS SELECT",
    "0::UINTEGER transcript_index, 100::UBIGINT exon_start,",
    "150::UBIGINT exon_end, 1::UBIGINT exon_cdna_start,",
    "51::UBIGINT exon_cdna_end, 0::TINYINT phase,",
    "0::TINYINT end_phase"
  )
)
expect_true(
  load_model(
    "r-hgvs-coding",
    c(
      "SELECT * FROM duckvep_r_hgvs_regions ORDER BY seq_region",
      paste(
        "SELECT * FROM duckvep_r_hgvs_coding_transcripts",
        "ORDER BY seq_region, transcript_start, transcript_index"
      ),
      paste(
        "SELECT * FROM duckvep_r_hgvs_coding_exons",
        "ORDER BY transcript_index, exon_cdna_start"
      )
    ),
    reference_fasta = hgvs_reference
  )$loaded
)
hgvs_coding <- dbGetQuery(
  con,
  paste(
    "SELECT a.transcript_hgvs, a.protein_hgvs,",
    "a.transcript_hgvs_status, a.protein_hgvs_status",
    "FROM unnest(_duckvep_annotate_small_hgvs(",
    "'r-hgvs-coding', 1::UINTEGER, 124::UBIGINT,",
    "'A', 'G', 0::UBIGINT",
    ")) AS u(a)"
  )
)
expect_identical(hgvs_coding$transcript_hgvs, "c.25A>G")
expect_identical(hgvs_coding$protein_hgvs, "p.Lys9Glu")
expect_identical(hgvs_coding$transcript_hgvs_status, "supported")
expect_identical(hgvs_coding$protein_hgvs_status, "supported")

# Protein HGVS uses its own exact-size retry after the initial native scratch
# fills. Exercise that adapter path through DBI with an in-frame insertion.
hgvs_long_protein <- dbGetQuery(
  con,
  paste(
    "SELECT a.transcript_hgvs, a.protein_hgvs,",
    "a.transcript_hgvs_status, a.protein_hgvs_status",
    "FROM unnest(_duckvep_annotate_small_hgvs(",
    "'r-hgvs-coding', 1::UINTEGER, 124::UBIGINT,",
    "'A', 'A' || repeat('GCT', 100), 0::UBIGINT",
    ")) AS u(a)"
  )
)
expect_true(nchar(hgvs_long_protein$transcript_hgvs, type = "bytes") > 256L)
expect_true(nchar(hgvs_long_protein$protein_hgvs, type = "bytes") > 256L)
expect_identical(hgvs_long_protein$transcript_hgvs_status, "supported")
expect_identical(hgvs_long_protein$protein_hgvs_status, "supported")

hgvs_directional <- dbGetQuery(
  con,
  paste(
    "WITH variants(ord, position) AS (VALUES",
    "(1, 90::UBIGINT), (2, 151::UBIGINT))",
    "SELECT ord, a.transcript_hgvs, a.transcript_hgvs_status,",
    "a.transcript_hgvs_reason",
    "FROM variants, LATERAL unnest(_duckvep_annotate_small_hgvs(",
    "'r-hgvs', 1::UINTEGER, position, 'A', 'C', 5000::UBIGINT",
    ")) AS u(a) ORDER BY ord"
  )
)
expect_true(all(is.na(hgvs_directional$transcript_hgvs)))
expect_identical(
  hgvs_directional$transcript_hgvs_status,
  c("not_applicable", "not_applicable")
)
expect_true(all(is.na(hgvs_directional$transcript_hgvs_reason)))

hgvs_terminal_insertion <- dbGetQuery(
  con,
  paste(
    "SELECT a.transcript_hgvs, a.transcript_hgvs_status,",
    "a.transcript_hgvs_reason, a.protein_hgvs_status",
    "FROM unnest(_duckvep_annotate_small_hgvs(",
    "'r-hgvs', 1::UINTEGER, 50000::UBIGINT, 'A', 'AC', 50000::UBIGINT",
    ")) AS u(a) WHERE a.transcript_index = 0"
  )
)
expect_true(is.na(hgvs_terminal_insertion$transcript_hgvs))
expect_identical(
  hgvs_terminal_insertion$transcript_hgvs_status,
  "not_applicable"
)
expect_true(is.na(hgvs_terminal_insertion$transcript_hgvs_reason))
expect_identical(hgvs_terminal_insertion$protein_hgvs_status, "not_applicable")

dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_hgvs_phase_transcripts AS SELECT",
    "0::UINTEGER transcript_index, 1::UINTEGER seq_region,",
    "100::UBIGINT transcript_start, 209::UBIGINT transcript_end,",
    "1::TINYINT strand, 0::UINTEGER gene_index,",
    "16::UBIGINT transcript_flags, 100::UBIGINT cds_start,",
    "209::UBIGINT cds_end, NULL::BLOB cds_sequence,",
    "NULL::UTINYINT codon_table, NULL::BLOB pre_cds_sequence,",
    "NULL::BLOB post_cds_sequence"
  )
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_hgvs_phase_exons AS SELECT * FROM (VALUES",
    "(0::UINTEGER, 100::UBIGINT, 109::UBIGINT, 1::UBIGINT,",
    "10::UBIGINT, 2::TINYINT, 0::TINYINT),",
    "(0::UINTEGER, 200::UBIGINT, 209::UBIGINT, 11::UBIGINT,",
    "20::UBIGINT, 0::TINYINT, 0::TINYINT)",
    ") t(transcript_index, exon_start, exon_end, exon_cdna_start,",
    "exon_cdna_end, phase, end_phase)"
  )
)
expect_true(
  load_model(
    "r-hgvs-phase",
    c(
      "SELECT * FROM duckvep_r_hgvs_regions ORDER BY seq_region",
      paste(
        "SELECT * FROM duckvep_r_hgvs_phase_transcripts",
        "ORDER BY seq_region, transcript_start, transcript_index"
      ),
      paste(
        "SELECT * FROM duckvep_r_hgvs_phase_exons",
        "ORDER BY transcript_index, exon_cdna_start"
      )
    ),
    reference_fasta = hgvs_reference
  )$loaded
)
hgvs_phase <- dbGetQuery(
  con,
  paste(
    "WITH variants(ord, position, reference, alternate) AS (VALUES",
    "(1, 100::UBIGINT, 'A', 'G'),",
    "(2, 100::UBIGINT, 'AA', 'AG'))",
    "SELECT ord, a.transcript_hgvs, a.transcript_hgvs_status",
    "FROM variants, LATERAL unnest(_duckvep_annotate_small_hgvs(",
    "'r-hgvs-phase', 1::UINTEGER, position, reference, alternate, 0::UBIGINT",
    ")) AS u(a) ORDER BY ord"
  )
)
expect_identical(hgvs_phase$transcript_hgvs, c("c.3A>G", "c.2A>G"))
expect_identical(
  hgvs_phase$transcript_hgvs_status,
  c("supported", "supported")
)
hgvs_protein_gap <- dbGetQuery(
  con,
  paste(
    "SELECT a.transcript_hgvs_status, a.protein_hgvs_status,",
    "a.protein_hgvs_reason FROM unnest(_duckvep_annotate_small_hgvs(",
    "'r-hgvs-phase', 1::UINTEGER, 109::UBIGINT, 'AAA', 'CCC', 0::UBIGINT",
    ")) AS u(a)"
  )
)
expect_identical(hgvs_protein_gap$transcript_hgvs_status, "supported")
expect_identical(hgvs_protein_gap$protein_hgvs_status, "not_applicable")
expect_true(is.na(hgvs_protein_gap$protein_hgvs_reason))

ensembl_mirna_queries <- c(
  paste(
    "SELECT seq_region, sequence_length",
    "FROM duckvep_r_ensembl_regions ORDER BY seq_region"
  ),
  paste(
    "SELECT transcript_index, seq_region, transcript_start, transcript_end,",
    "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence,",
    "codon_table, pre_cds_sequence, post_cds_sequence",
    "FROM duckvep_r_ensembl_transcripts",
    "ORDER BY seq_region, transcript_start, transcript_index"
  ),
  paste(
    "SELECT transcript_index, exon.exon_start, exon.exon_end,",
    "exon.exon_cdna_start, exon.exon_cdna_end, exon.phase, exon.end_phase",
    "FROM duckvep_r_ensembl_transcripts,",
    "LATERAL unnest(exons) AS u(exon)",
    "ORDER BY transcript_index, exon.exon_cdna_start"
  )
)
ensembl_mirna_query <- paste(
  "SELECT transcript_index, region.mature_mirna_start,",
  "region.mature_mirna_end FROM duckvep_r_ensembl_transcripts,",
  "LATERAL unnest(mature_mirna_regions) AS u(region)",
  "ORDER BY transcript_index, region.mature_mirna_start"
)
ensembl_regulation_query <- paste(
  "SELECT regulation_feature_index, seq_region, feature_start, feature_end,",
  "feature_kind FROM duckvep_r_ensembl_regulation",
  "ORDER BY seq_region, feature_start, regulation_feature_index"
)
expect_true(
  load_model(
    "r-ensembl-mirna",
    ensembl_mirna_queries,
    transcript_coverage_complete = TRUE,
    mature_mirna_query = ensembl_mirna_query,
    interval_feature_query = ensembl_regulation_query
  )$loaded
)
hgvs_retained_ref_without_fasta <- dbGetQuery(
  con,
  paste(
    "SELECT a.transcript_hgvs_status, a.transcript_hgvs_reason,",
    "a.protein_hgvs_status, a.protein_hgvs_reason",
    "FROM unnest(_duckvep_annotate_small_hgvs(",
    "'r-ensembl-mirna', 0::UINTEGER, 1::UBIGINT,",
    "'CT', 'CG', 0::UBIGINT)) u(a)",
    "WHERE a.transcript_index = 0"
  )
)
expect_identical(
  hgvs_retained_ref_without_fasta$transcript_hgvs_status,
  "unresolved"
)
expect_identical(
  hgvs_retained_ref_without_fasta$transcript_hgvs_reason,
  "missing_reference"
)
expect_identical(
  hgvs_retained_ref_without_fasta$protein_hgvs_status,
  "unresolved"
)
expect_identical(
  hgvs_retained_ref_without_fasta$protein_hgvs_reason,
  "missing_reference"
)
ensembl_mirna_annotation <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence, a.status FROM unnest(_duckvep_annotate_small_rich(",
    "'r-ensembl-mirna', 0::UINTEGER, 2::UBIGINT,",
    "'T', 'C', 0::UBIGINT)) u(a) WHERE a.transcript_index = 1"
  )
)
expect_identical(
  ensembl_mirna_annotation$consequence,
  "mature_miRNA_variant"
)
expect_identical(ensembl_mirna_annotation$status, "supported")

regulation_annotations <- dbGetQuery(
  con,
  paste(
    "WITH small AS (SELECT a.* FROM unnest(_duckvep_annotate_small_rich(",
    "'r-ensembl-mirna', 0::UINTEGER, 3::UBIGINT, 'G', 'A', 0::UBIGINT",
    ")) u(a) WHERE a.regulation_feature_index IS NOT NULL),",
    "structural AS (SELECT a.* FROM unnest(_duckvep_annotate_structural_rich(",
    "'r-ensembl-mirna', 0::UINTEGER, 1::UBIGINT, 12::UBIGINT,",
    "'DUP', 'GAIN', 0::UBIGINT)) u(a)",
    "WHERE a.regulation_feature_index IS NOT NULL)",
    "SELECT 1 ord, regulation_feature_index, consequence, overlap_object",
    "FROM small UNION ALL SELECT 2, regulation_feature_index, consequence,",
    "overlap_object FROM structural ORDER BY ord, regulation_feature_index"
  )
)
expect_equal(regulation_annotations$regulation_feature_index, c(0, 0, 1))
expect_identical(
  regulation_annotations$consequence,
  c(
    "regulatory_region_variant",
    "regulatory_region_amplification&regulatory_region_variant",
    "TFBS_amplification&TF_binding_site_variant"
  )
)
expect_identical(
  regulation_annotations$overlap_object,
  c(
    "regulatory_region",
    "regulatory_region",
    "transcription_factor_binding_site"
  )
)

hgvs_regulation <- dbGetQuery(
  con,
  paste(
    "SELECT a.regulation_feature_index, a.overlap_object_code,",
    "a.transcript_hgvs, a.protein_hgvs, a.transcript_hgvs_status,",
    "a.protein_hgvs_status FROM unnest(_duckvep_annotate_small_hgvs(",
    "'r-ensembl-mirna', 0::UINTEGER, 3::UBIGINT,",
    "'G', 'A', 0::UBIGINT)) u(a)",
    "WHERE a.regulation_feature_index IS NOT NULL"
  )
)
expect_equal(hgvs_regulation$regulation_feature_index, 0)
expect_equal(hgvs_regulation$overlap_object_code, 1)
expect_true(all(is.na(hgvs_regulation[, 3:6])))

breakend_regulation <- dbGetQuery(
  con,
  paste(
    "SELECT a.regulation_feature_index, a.consequence, a.overlap_object",
    "FROM unnest(_duckvep_annotate_breakend_rich(",
    "'r-ensembl-mirna', 0::UINTEGER, 1::UBIGINT,",
    "0::UINTEGER, 8::UBIGINT, 0::UBIGINT)) u(a)",
    "WHERE a.regulation_feature_index IS NOT NULL",
    "ORDER BY a.regulation_feature_index"
  )
)
expect_equal(breakend_regulation$regulation_feature_index, c(0, 1))
expect_identical(
  breakend_regulation$consequence,
  c(
    "regulatory_region_variant",
    "feature_truncation&intergenic_variant"
  )
)
expect_identical(
  breakend_regulation$overlap_object,
  c("regulatory_region", "transcription_factor_binding_site")
)

breakend_mixed_vectors <- dbGetQuery(
  con,
  paste(
    "WITH events(ord, local_position, mate_position) AS (VALUES",
    "(1, 1::UBIGINT, 8::UBIGINT), (2, 2::UBIGINT, 9::UBIGINT))",
    "SELECT ord, count(*)::INTEGER result_count,",
    "count(*) FILTER (WHERE a.transcript_index IS NOT NULL)::INTEGER",
    "transcript_count,",
    "count(*) FILTER (WHERE a.regulation_feature_index IS NOT NULL)::INTEGER",
    "regulation_count FROM events,",
    "LATERAL unnest(_duckvep_annotate_breakend_compact(",
    "'r-ensembl-mirna', 0::UINTEGER, local_position,",
    "0::UINTEGER, mate_position, 0::UBIGINT)) u(a)",
    "GROUP BY ord ORDER BY ord"
  )
)
expect_equal(breakend_mixed_vectors$ord, c(1, 2))
expect_equal(breakend_mixed_vectors$result_count, c(4, 4))
expect_equal(breakend_mixed_vectors$transcript_count, c(2, 2))
expect_equal(breakend_mixed_vectors$regulation_count, c(2, 2))

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
    "FROM unnest(_duckvep_annotate_small_rich(",
    "'r-partial-cds-end', 1::UINTEGER, 106::UBIGINT,",
    "'C', 'CA', 0::UBIGINT)) u(a)"
  )
)
expect_identical(partial_cds_annotation$consequence, "frameshift_variant")
expect_identical(partial_cds_annotation$status, "supported")
expect_true(is.na(partial_cds_annotation$reason))

partial_tail_deletion <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence, a.status, a.reason",
    "FROM unnest(_duckvep_annotate_small_rich(",
    "'r-partial-cds-end', 1::UINTEGER, 109::UBIGINT,",
    "'GGGT', 'G', 0::UBIGINT)) u(a)"
  )
)
expect_identical(partial_tail_deletion$consequence, "inframe_deletion")
expect_identical(partial_tail_deletion$status, "supported")
expect_true(is.na(partial_tail_deletion$reason))

dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_complete_span_transcripts AS SELECT",
    "0::UINTEGER transcript_index, 1::UINTEGER seq_region,",
    "100::UBIGINT transcript_start, 108::UBIGINT transcript_end,",
    "1::TINYINT strand, 0::UINTEGER gene_index,",
    "3::UBIGINT transcript_flags, 100::UBIGINT cds_start,",
    "108::UBIGINT cds_end, 'ATGAAATAA'::BLOB cds_sequence,",
    "1::UTINYINT codon_table, ''::BLOB pre_cds_sequence,",
    "''::BLOB post_cds_sequence"
  )
)
dbExecute(
  con,
  paste(
    "CREATE TABLE duckvep_r_complete_span_exons AS SELECT",
    "0::UINTEGER transcript_index, 100::UBIGINT exon_start,",
    "108::UBIGINT exon_end, 1::UBIGINT exon_cdna_start,",
    "9::UBIGINT exon_cdna_end, 0::TINYINT phase,",
    "0::TINYINT end_phase"
  )
)
complete_span_queries <- c(
  "SELECT seq_region FROM duckvep_r_regions ORDER BY seq_region",
  "SELECT * FROM duckvep_r_complete_span_transcripts",
  "SELECT * FROM duckvep_r_complete_span_exons"
)
expect_true(load_model("r-complete-span", complete_span_queries)$loaded)
complete_span_annotations <- dbGetQuery(
  con,
  paste(
    "WITH variants(ord, reference, alternate) AS (VALUES",
    "(1, repeat('A', 11), repeat('A', 5) || 'C' || repeat('A', 5)),",
    "(2, repeat('A', 11), 'A'))",
    "SELECT ord, a.consequence, a.status FROM variants,",
    "LATERAL unnest(_duckvep_annotate_small_rich(",
    "'r-complete-span', 1::UINTEGER, 99::UBIGINT,",
    "reference, alternate, 0::UBIGINT)) u(a) ORDER BY ord"
  )
)
expect_identical(
  complete_span_annotations$consequence,
  c(
    "5_prime_UTR_variant&3_prime_UTR_variant&coding_transcript_variant",
    "transcript_ablation"
  )
)
expect_true(all(complete_span_annotations$status == "supported"))

noncoding_boundary_queries <- queries
noncoding_boundary_queries[2] <- paste(
  "SELECT transcript_index, seq_region, transcript_start, transcript_end,",
  "strand, gene_index, 0::UBIGINT transcript_flags,",
  "NULL::UBIGINT cds_start, NULL::UBIGINT cds_end,",
  "NULL::BLOB cds_sequence, NULL::UTINYINT codon_table,",
  "NULL::BLOB pre_cds_sequence, NULL::BLOB post_cds_sequence",
  "FROM duckvep_r_transcripts ORDER BY seq_region, transcript_start, transcript_index"
)
expect_true(
  load_model("r-noncoding-boundary", noncoding_boundary_queries)$loaded
)
noncoding_boundary <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence, a.status FROM unnest(_duckvep_annotate_small_rich(",
    "'r-noncoding-boundary', 1::UINTEGER, 150::UBIGINT,",
    "'A', 'AT', 0::UBIGINT)) u(a)"
  )
)
expect_identical(
  sort(strsplit(noncoding_boundary$consequence, "&", fixed = TRUE)[[1]]),
  c("non_coding_transcript_variant", "splice_region_variant")
)
expect_identical(noncoding_boundary$status, "supported")

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
    "FROM unnest(_duckvep_annotate_small_rich(",
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
    "SELECT a.* FROM unnest(_duckvep_annotate_small_compact(",
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

structural_deletion <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence, a.status FROM unnest(_duckvep_annotate_structural_rich(",
    "'r-test', 1::UINTEGER, 90::UBIGINT, 260::UBIGINT,",
    "'DEL', 'LOSS', 0::UBIGINT)) u(a)"
  )
)
expect_identical(structural_deletion$consequence, "transcript_ablation")
expect_identical(structural_deletion$status, "supported")

structural_duplication <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence_mask, a.status_code",
    "FROM unnest(_duckvep_annotate_structural_compact(",
    "'r-test', 1::UINTEGER, 90::UBIGINT, 260::UBIGINT,",
    "'DUP', 'GAIN', 0::UBIGINT)) u(a)"
  )
)
expect_equal(structural_duplication$consequence_mask, 67108864)
expect_equal(structural_duplication$status_code, 0)

structural_tandem_repeat <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence, a.status",
    "FROM unnest(_duckvep_annotate_structural_rich(",
    "'r-test', 1::UINTEGER, 90::UBIGINT, 260::UBIGINT,",
    "'STR', 'UNKNOWN', 0::UBIGINT)) u(a)"
  )
)
expect_identical(
  structural_tandem_repeat$consequence,
  "transcript_amplification"
)
expect_identical(structural_tandem_repeat$status, "supported")

breakend <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence, a.region, a.status",
    "FROM unnest(_duckvep_annotate_breakend_rich(",
    "'r-test', 1::UINTEGER, 159::UBIGINT,",
    "1::UINTEGER, 170::UBIGINT, 0::UBIGINT)) u(a)"
  )
)
expect_identical(
  breakend$consequence,
  "feature_truncation&intron_variant"
)
expect_identical(breakend$region, "intron")
expect_identical(breakend$status, "supported")

breakend_compact <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence_mask, a.region_mask, a.status_code",
    "FROM unnest(_duckvep_annotate_breakend_compact(",
    "'r-test', 1::UINTEGER, 159::UBIGINT,",
    "1::UINTEGER, 170::UBIGINT, 0::UBIGINT)) u(a)"
  )
)
expect_equal(breakend_compact$consequence_mask, 268435460)
expect_equal(breakend_compact$region_mask, 4)
expect_equal(breakend_compact$status_code, 0)

mate_only_breakend <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence, a.region FROM unnest(_duckvep_annotate_breakend_rich(",
    "'r-test', 1::UINTEGER, 5250::UBIGINT,",
    "1::UINTEGER, 170::UBIGINT, 0::UBIGINT)) u(a)"
  )
)
expect_identical(mate_only_breakend$consequence, "feature_truncation")
expect_true(is.na(mate_only_breakend$region))

# VEP's configurable transcript window and its fixed structural-breakend
# allele-admission cap both default to 5000, but are not the same setting. A
# wider 10000-base caller window does not create a local endpoint allele at
# 5001, but ordinary predicates on the mate allele still inspect that local
# feature and retain the wider directional term.
breakend_admission_cap <- dbGetQuery(
  con,
  paste(
    "WITH events(ord, local_position) AS (VALUES",
    "(1, 5249::UBIGINT), (2, 5250::UBIGINT))",
    "SELECT ord, a.consequence, a.region FROM events,",
    "LATERAL unnest(_duckvep_annotate_breakend_rich(",
    "'r-test', 1::UINTEGER, local_position,",
    "1::UINTEGER, 124::UBIGINT, 10000::UBIGINT)) u(a)",
    "ORDER BY ord"
  )
)
expect_equal(breakend_admission_cap$ord, c(1, 2))
expect_identical(
  breakend_admission_cap$consequence,
  rep("feature_truncation&downstream_gene_variant", 2)
)
expect_identical(breakend_admission_cap$region, rep("downstream", 2))

# The fixed StructuralVariationOverlap endpoint admission is independent of
# the caller-managed directional window in the other direction too. With a
# zero transcript window, the exactly-5000 local allele contributes VEP's
# default intergenic term and the mate allele contributes feature_truncation.
breakend_zero_window <- dbGetQuery(
  con,
  paste(
    "WITH events(ord, local_position) AS (VALUES",
    "(1, 5249::UBIGINT), (2, 5250::UBIGINT))",
    "SELECT ord, a.consequence, a.region FROM events,",
    "LATERAL unnest(_duckvep_annotate_breakend_rich(",
    "'r-test', 1::UINTEGER, local_position,",
    "1::UINTEGER, 124::UBIGINT, 0::UBIGINT)) u(a)",
    "ORDER BY ord"
  )
)
expect_equal(breakend_zero_window$ord, c(1, 2))
expect_identical(
  breakend_zero_window$consequence,
  c("feature_truncation&intergenic_variant", "feature_truncation")
)
expect_true(all(is.na(breakend_zero_window$region)))

expect_error(
  dbGetQuery(
    con,
    paste(
      "SELECT _duckvep_annotate_structural_rich(",
      "'r-test', 1::UINTEGER, 100::UBIGINT, 200::UBIGINT,",
      "'DEL', 'GAIN')"
    )
  ),
  pattern = "structural type contradicts copy change"
)

expect_error(
  dbGetQuery(
    con,
    paste(
      "SELECT _duckvep_annotate_structural_rich(",
      "'r-test', 1::UINTEGER, 100::UBIGINT, 100::UBIGINT,",
      "'BND', 'UNKNOWN')"
    )
  ),
  pattern = "provide both mate coordinates for a breakend"
)

expect_error(
  dbGetQuery(
    con,
    paste(
      "SELECT _duckvep_annotate_breakend_rich(",
      "'r-test', 1::UINTEGER, 0::UBIGINT,",
      "1::UINTEGER, 170::UBIGINT)"
    )
  ),
  pattern = "invalid one-based endpoint coordinate"
)

reference_mismatch <- dbGetQuery(
  con,
  paste(
    "SELECT a.status, a.reason FROM unnest(_duckvep_annotate_small_rich(",
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
    "a.alternate_amino_acid FROM unnest(_duckvep_annotate_small_rich(",
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
    "LATERAL unnest(_duckvep_annotate_small_rich(",
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
    "LATERAL unnest(_duckvep_annotate_small_rich(",
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
    "LATERAL unnest(_duckvep_annotate_small_rich(",
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
    "LATERAL unnest(_duckvep_annotate_small_rich(",
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
    "LATERAL unnest(_duckvep_annotate_small_rich(",
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
    "SELECT a.consequence, a.status, a.reason FROM unnest(_duckvep_annotate_small_rich(",
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
    "SELECT a.consequence, a.status, a.reason FROM unnest(_duckvep_annotate_small_rich(",
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
    "LATERAL unnest(_duckvep_annotate_small_rich(",
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
    "SELECT a.consequence, a.status FROM unnest(_duckvep_annotate_small_rich(",
    "'r-complete', 1::UINTEGER, 1::UBIGINT, 'A', 'C', 0::UBIGINT",
    ")) u(a)"
  )
)
expect_identical(complete_intergenic$consequence, "intergenic_variant")
expect_identical(complete_intergenic$status, "supported")
asymmetric_upstream <- dbGetQuery(
  con,
  paste(
    "SELECT a.consequence, a.region FROM unnest(_duckvep_annotate_small_rich(",
    "'r-complete', 1::UINTEGER, 90::UBIGINT, 'A', 'C',",
    "10::UBIGINT, 0::UBIGINT)) u(a)"
  )
)
expect_identical(asymmetric_upstream$consequence, "upstream_gene_variant")
expect_identical(asymmetric_upstream$region, "upstream")
asymmetric_downstream <- dbGetQuery(
  con,
  paste(
    "SELECT a.region_mask, a.status_code",
    "FROM unnest(_duckvep_annotate_small_compact(",
    "'r-complete', 1::UINTEGER, 260::UBIGINT, 'A', 'C',",
    "0::UBIGINT, 10::UBIGINT)) u(a)"
  )
)
expect_equal(asymmetric_downstream$region_mask, 2)
expect_equal(asymmetric_downstream$status_code, 0)
expect_error(
  dbGetQuery(
    con,
    paste(
      "SELECT _duckvep_annotate_small_rich(",
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
    "SELECT a.consequence FROM unnest(_duckvep_annotate_small_rich(",
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
    "LATERAL unnest(_duckvep_annotate_small_rich(",
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
  "seq_region_attrib",
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
  regulation <- paste0("duckvep_r_", directory, "_regulation")
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
  if (identical(directory, "grch38")) {
    funcgen_schema <- "duckvep_r_grch38_funcgen"
    dbExecute(con, paste("CREATE SCHEMA", funcgen_schema))
    for (table in c("feature_type", "regulatory_feature", "motif_feature")) {
      path <- normalizePath(
        file.path(fixture_dir, paste0("funcgen_", table, ".parquet")),
        mustWork = TRUE
      )
      dbExecute(
        con,
        paste0(
          "CREATE VIEW ", funcgen_schema, ".",
          as.character(dbQuoteIdentifier(con, table)),
          " AS SELECT * FROM read_parquet(",
          as.character(dbQuoteString(con, path)), ")"
        )
      )
    }
    dbExecute(
      con,
      paste0(
        "CREATE TABLE ", regulation,
        " AS SELECT * FROM duckvep_ensembl_regulation_features(",
        as.character(dbQuoteString(con, funcgen_schema)), ", ",
        as.character(dbQuoteString(con, regions)), ")"
      )
    )
  } else {
    dbExecute(
      con,
      paste0(
        "CREATE TABLE ", regulation, "(",
        "regulation_feature_index UINTEGER, seq_region UINTEGER, ",
        "feature_start UINTEGER, feature_end UINTEGER, feature_kind UTINYINT)"
      )
    )
  }
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
  peptide_edit_query <- paste(
    "SELECT transcript_index, edit.protein_position,",
    "edit.alternate_amino_acid FROM", transcripts, ",",
    "LATERAL unnest(peptide_edits) AS u(edit)",
    "ORDER BY transcript_index, edit.protein_position"
  )
  interval_feature_query <- paste(
    "SELECT regulation_feature_index, seq_region, feature_start,",
    "feature_end, feature_kind FROM", regulation,
    "ORDER BY regulation_feature_index"
  )
  expect_true(
    load_model(
      model,
      model_queries,
      transcript_coverage_complete = TRUE,
      peptide_edit_query = peptide_edit_query,
      interval_feature_query = interval_feature_query
    )$loaded
  )
  list(
    model = model,
    regions = regions,
    transcripts = transcripts,
    regulation = regulation,
    summary = summary
  )
}

grch38_fixture <- prepare_ensembl_fixture("grch38", "GRCh38")
grch37_fixture <- prepare_ensembl_fixture("grch37", "GRCh37")
expect_equal(
  unlist(grch38_fixture$summary[1, ], use.names = FALSE),
  c(39, 14, 0, 40, 1)
)
expect_equal(
  unlist(grch37_fixture$summary[1, ], use.names = FALSE),
  c(39, 15, 0, 48, 0)
)
fixture_regulation_summary <- dbGetQuery(
  con,
  paste(
    "SELECT (SELECT count(*) FROM", grch38_fixture$regions,
    ") region_count, count(*) feature_count,",
    "count(*) FILTER (WHERE feature_kind = 2) motif_count FROM",
    grch38_fixture$regulation
  )
)
expect_equal(
  unlist(fixture_regulation_summary[1, ], use.names = FALSE),
  c(3, 3, 3)
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
    "FROM unnest(_duckvep_annotate_small_rich(",
    as.character(dbQuoteString(con, grch38_fixture$model)), ",",
    "2::UINTEGER, 32522::UBIGINT, 'C', 'T', 0::UBIGINT)) u(a)"
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
    "FROM unnest(_duckvep_annotate_small_rich(",
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
    "FROM unnest(_duckvep_annotate_small_rich(",
    as.character(dbQuoteString(con, grch38_fixture$model)), ",",
    "0::UINTEGER, 3307::UBIGINT,",
    "'A', 'C', 0::UBIGINT)) u(a)"
  )
)
expect_equal(fixture_mitochondrial$transcript_index, 5)
expect_identical(
  fixture_mitochondrial$consequence,
  "start_lost"
)
expect_identical(fixture_mitochondrial$status, "supported")
expect_true(is.na(fixture_mitochondrial$reason))

fixture_motifs <- dbGetQuery(
  con,
  paste(
    "SELECT a.regulation_feature_index, a.overlap_object,",
    "a.consequence, a.status FROM unnest(_duckvep_annotate_small_rich(",
    as.character(dbQuoteString(con, grch38_fixture$model)), ",",
    "1::UINTEGER, 75::UBIGINT, 'A', 'C', 0::UBIGINT)) u(a)",
    "WHERE a.overlap_object = 'transcription_factor_binding_site'",
    "ORDER BY a.regulation_feature_index"
  )
)
expect_equal(fixture_motifs$regulation_feature_index, c(0, 1))
expect_true(all(
  fixture_motifs$overlap_object == "transcription_factor_binding_site"
))
expect_true(all(
  fixture_motifs$consequence == "TF_binding_site_variant"
))
expect_true(all(fixture_motifs$status == "supported"))

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
    "SELECT duckvep_model_drop('r-public-terminal-partial') AS dropped"
  )$dropped
)
expect_true(
  dbGetQuery(
    con,
    "SELECT duckvep_model_drop('r-clinvar-terminal-partial') AS dropped"
  )$dropped
)
expect_true(
  dbGetQuery(
    con,
    "SELECT duckvep_model_drop('r-complete-span') AS dropped"
  )$dropped
)
expect_true(
  dbGetQuery(
    con,
    "SELECT duckvep_model_drop('r-noncoding-boundary') AS dropped"
  )$dropped
)
expect_true(
  dbGetQuery(
    con,
    "SELECT duckvep_model_drop('r-ensembl-mirna') AS dropped"
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
