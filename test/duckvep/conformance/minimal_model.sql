-- Resident model matching test/data/duckvep/minimal.{gff3,fa}.
-- The corpus differential uses the same four ordinary DuckDB relations for a
-- small fixture or a production-sized, separately prepared model database.

CREATE OR REPLACE TABLE duckvep_sequence_regions AS
SELECT * FROM (VALUES
  (0::UINTEGER, 'chrEmpty'::VARCHAR),
  (1::UINTEGER, 'chrDuck'::VARCHAR)
) t(seq_region, name);

CREATE OR REPLACE TABLE duckvep_transcripts AS
SELECT
  0::UINTEGER AS transcript_index,
  1::UINTEGER AS seq_region,
  100::UBIGINT AS transcript_start,
  250::UBIGINT AS transcript_end,
  1::TINYINT AS strand,
  0::UINTEGER AS gene_index,
  3::UBIGINT AS transcript_flags,
  120::UBIGINT AS cds_start,
  240::UBIGINT AS cds_end,
  'ATGGTACGTACGTACGTACGTACGTACGTACTACGTACGTACGTACGTACGTACGTACGTACGTACTGGTAA'::BLOB AS cds_sequence,
  1::UTINYINT AS codon_table,
  'ACG'::BLOB AS post_cds_bases,
  'TACGTACGTACGTACGTACG'::BLOB AS pre_cds_sequence,
  'ACGTACGTAC'::BLOB AS post_cds_sequence;

CREATE OR REPLACE TABLE duckvep_exons AS
SELECT * FROM (VALUES
  (0::UINTEGER, 100::UBIGINT, 150::UBIGINT, 1::UBIGINT, 51::UBIGINT, 0::TINYINT, 0::TINYINT),
  (0::UINTEGER, 200::UBIGINT, 250::UBIGINT, 52::UBIGINT, 102::UBIGINT, 0::TINYINT, 0::TINYINT)
) t(transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end, phase, end_phase);

CREATE OR REPLACE TABLE duckvep_transcript_names AS
SELECT 0::UINTEGER AS transcript_index, 'DUCK1-201'::VARCHAR AS transcript_id;
