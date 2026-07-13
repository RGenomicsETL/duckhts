-- Regenerate the offline DuckVEP Ensembl fixtures from the public release-116
-- databases. Run from the repository root after building the release extension:
--
--   mkdir -p test/data/duckvep/ensembl_core/{grch38,grch37}
--   duckdb -unsigned < test/scripts/prepare_duckvep_ensembl_fixture.sql
--   ./test/scripts/prepare_test_data.sh
--
-- This is an explicit networked staging step. Extension builds and tests read only
-- the checked-in Parquet files produced below.

INSTALL mysql;
LOAD mysql;
LOAD 'build/release/duckhts.duckdb_extension';

-- The archive server predates START TRANSACTION READ ONLY. These are immutable
-- anonymous source catalogs, so staging does not need remote transactions.
SET mysql_enable_transactions = false;
SET threads = 1;

ATTACH
  'host=ensembldb.ensembl.org port=5306 user=anonymous database=homo_sapiens_core_116_38'
  AS ensembl_grch38 (TYPE mysql, READ_ONLY);
ATTACH
  'host=ensembldb.ensembl.org port=3337 user=anonymous database=homo_sapiens_core_116_37'
  AS ensembl_grch37 (TYPE mysql, READ_ONLY);

-- These staging macros preserve the source table shapes. mysql_query keeps every
-- source join on the Ensembl server; query_table would pull large join inputs into
-- DuckDB before applying the two-region selection. region_ids_sql is pinned text in
-- this script, never user input.
CREATE OR REPLACE MACRO fixture_seq_regions(
  source_catalog, assembly_name, region_ids_sql
) AS TABLE
FROM mysql_query(source_catalog,
  'SELECT sr.* FROM seq_region sr '
  'JOIN coord_system cs USING (coord_system_id) '
  'WHERE cs.species_id = 1 AND cs.version = ''' || assembly_name || ''' '
  'AND sr.seq_region_id IN (' || region_ids_sql || ')'
);

CREATE OR REPLACE MACRO fixture_coord_systems(
  source_catalog, assembly_name, region_ids_sql
) AS TABLE
FROM mysql_query(source_catalog,
  'SELECT DISTINCT cs.* FROM coord_system cs '
  'JOIN seq_region sr USING (coord_system_id) '
  'WHERE cs.species_id = 1 AND cs.version = ''' || assembly_name || ''' '
  'AND sr.seq_region_id IN (' || region_ids_sql || ')'
);

CREATE OR REPLACE MACRO fixture_transcripts(
  source_catalog, assembly_name, region_ids_sql
) AS TABLE
FROM mysql_query(source_catalog,
  'SELECT t.* FROM transcript t '
  'JOIN seq_region sr USING (seq_region_id) '
  'JOIN coord_system cs USING (coord_system_id) '
  'WHERE cs.species_id = 1 AND cs.version = ''' || assembly_name || ''' '
  'AND sr.seq_region_id IN (' || region_ids_sql || ')'
);

CREATE OR REPLACE MACRO fixture_genes(
  source_catalog, assembly_name, region_ids_sql
) AS TABLE
FROM mysql_query(source_catalog,
  'SELECT DISTINCT g.* FROM gene g '
  'JOIN transcript t USING (gene_id) '
  'JOIN seq_region sr ON sr.seq_region_id = t.seq_region_id '
  'JOIN coord_system cs USING (coord_system_id) '
  'WHERE cs.species_id = 1 AND cs.version = ''' || assembly_name || ''' '
  'AND sr.seq_region_id IN (' || region_ids_sql || ')'
);

CREATE OR REPLACE MACRO fixture_translations(
  source_catalog, assembly_name, region_ids_sql
) AS TABLE
FROM mysql_query(source_catalog,
  'SELECT tl.* FROM translation tl '
  'JOIN transcript t USING (transcript_id) '
  'JOIN seq_region sr ON sr.seq_region_id = t.seq_region_id '
  'JOIN coord_system cs USING (coord_system_id) '
  'WHERE cs.species_id = 1 AND cs.version = ''' || assembly_name || ''' '
  'AND sr.seq_region_id IN (' || region_ids_sql || ')'
);

CREATE OR REPLACE MACRO fixture_exon_transcripts(
  source_catalog, assembly_name, region_ids_sql
) AS TABLE
FROM mysql_query(source_catalog,
  'SELECT et.* FROM exon_transcript et '
  'JOIN transcript t USING (transcript_id) '
  'JOIN seq_region sr ON sr.seq_region_id = t.seq_region_id '
  'JOIN coord_system cs USING (coord_system_id) '
  'WHERE cs.species_id = 1 AND cs.version = ''' || assembly_name || ''' '
  'AND sr.seq_region_id IN (' || region_ids_sql || ')'
);

CREATE OR REPLACE MACRO fixture_exons(
  source_catalog, assembly_name, region_ids_sql
) AS TABLE
FROM mysql_query(source_catalog,
  'SELECT DISTINCT e.* FROM exon e '
  'JOIN exon_transcript et USING (exon_id) '
  'JOIN transcript t USING (transcript_id) '
  'JOIN seq_region sr ON sr.seq_region_id = t.seq_region_id '
  'JOIN coord_system cs USING (coord_system_id) '
  'WHERE cs.species_id = 1 AND cs.version = ''' || assembly_name || ''' '
  'AND sr.seq_region_id IN (' || region_ids_sql || ')'
);

CREATE OR REPLACE MACRO fixture_transcript_attributes(
  source_catalog, assembly_name, region_ids_sql
) AS TABLE
FROM mysql_query(source_catalog,
  'SELECT ta.* FROM transcript_attrib ta '
  'JOIN transcript t USING (transcript_id) '
  'JOIN seq_region sr ON sr.seq_region_id = t.seq_region_id '
  'JOIN coord_system cs USING (coord_system_id) '
  'WHERE cs.species_id = 1 AND cs.version = ''' || assembly_name || ''' '
  'AND sr.seq_region_id IN (' || region_ids_sql || ')'
);

CREATE OR REPLACE MACRO fixture_translation_attributes(
  source_catalog, assembly_name, region_ids_sql
) AS TABLE
FROM mysql_query(source_catalog,
  'SELECT tla.* FROM translation_attrib tla '
  'JOIN translation tl USING (translation_id) '
  'JOIN transcript t USING (transcript_id) '
  'JOIN seq_region sr ON sr.seq_region_id = t.seq_region_id '
  'JOIN coord_system cs USING (coord_system_id) '
  'WHERE cs.species_id = 1 AND cs.version = ''' || assembly_name || ''' '
  'AND sr.seq_region_id IN (' || region_ids_sql || ')'
);

CREATE OR REPLACE MACRO fixture_attribute_types(
  source_catalog, assembly_name, region_ids_sql
) AS TABLE
FROM mysql_query(source_catalog, 'SELECT * FROM attrib_type');

CREATE OR REPLACE MACRO fixture_reference_chunks(
  source_catalog, assembly_name, region_ids_sql
) AS TABLE
WITH components AS MATERIALIZED (
  FROM mysql_query(source_catalog,
    'SELECT sr.seq_region_id, sr.name AS chrom, sr.length AS sequence_length, '
    'a.asm_start, a.asm_end, a.cmp_start, a.cmp_end, a.ori, d.sequence '
    'FROM seq_region sr '
    'JOIN coord_system cs USING (coord_system_id) '
    'JOIN assembly a ON a.asm_seq_region_id = sr.seq_region_id '
    'JOIN dna d ON d.seq_region_id = a.cmp_seq_region_id '
    'WHERE cs.species_id = 1 AND cs.version = ''' || assembly_name || ''' '
    'AND sr.seq_region_id IN (' || region_ids_sql || ') '
    'ORDER BY sr.seq_region_id, a.asm_start'
  )
), ordered AS MATERIALIZED (
  SELECT *, lag(asm_end) OVER (
    PARTITION BY seq_region_id ORDER BY asm_start
  ) AS previous_end
  FROM components
), validation AS MATERIALIZED (
  SELECT CASE
    WHEN (SELECT count(DISTINCT seq_region_id) FROM ordered) != 2 THEN
      error('DuckVEP fixture: every selected region must have assembly components')
    WHEN EXISTS (
      SELECT 1 FROM ordered
      WHERE ori NOT IN (-1, 1)
         OR asm_end - asm_start != cmp_end - cmp_start
         OR cmp_start < 1
         OR cmp_end > length(sequence)
         OR (previous_end IS NULL AND asm_start != 1)
         OR (previous_end IS NOT NULL AND asm_start != previous_end + 1)
    ) THEN error('DuckVEP fixture: invalid or discontinuous assembly components')
    WHEN EXISTS (
      SELECT 1
      FROM ordered
      GROUP BY seq_region_id, sequence_length
      HAVING max(asm_end) != sequence_length
    ) THEN error('DuckVEP fixture: assembly components do not cover a selected region')
    ELSE true
  END AS valid
), pieces AS MATERIALIZED (
  SELECT seq_region_id, chrom, sequence_length, asm_start,
         CASE ori
           WHEN 1 THEN substring(sequence, cmp_start, cmp_end - cmp_start + 1)
           ELSE seq_revcomp(
             substring(sequence, cmp_start, cmp_end - cmp_start + 1)
           )
         END AS piece
  FROM ordered
), assembled AS MATERIALIZED (
  SELECT seq_region_id, chrom, sequence_length,
         string_agg(piece, '' ORDER BY asm_start) AS seq
  FROM pieces
  GROUP BY seq_region_id, chrom, sequence_length
)
SELECT chrom, 0::BIGINT AS start, sequence_length::BIGINT AS "end", seq
FROM assembled
CROSS JOIN validation
WHERE validation.valid AND length(seq) = sequence_length
ORDER BY chrom;

CREATE OR REPLACE SCHEMA duckvep_grch38_core;
CREATE TABLE duckvep_grch38_core.coord_system AS
FROM fixture_coord_systems(
  'ensembl_grch38', 'GRCh38', '132907,2006446537'
);
CREATE TABLE duckvep_grch38_core.seq_region AS
FROM fixture_seq_regions(
  'ensembl_grch38', 'GRCh38', '132907,2006446537'
);
CREATE TABLE duckvep_grch38_core.transcript AS
FROM fixture_transcripts(
  'ensembl_grch38', 'GRCh38', '132907,2006446537'
);
CREATE TABLE duckvep_grch38_core.gene AS
FROM fixture_genes(
  'ensembl_grch38', 'GRCh38', '132907,2006446537'
);
CREATE TABLE duckvep_grch38_core.translation AS
FROM fixture_translations(
  'ensembl_grch38', 'GRCh38', '132907,2006446537'
);
CREATE TABLE duckvep_grch38_core.exon_transcript AS
FROM fixture_exon_transcripts(
  'ensembl_grch38', 'GRCh38', '132907,2006446537'
);
CREATE TABLE duckvep_grch38_core.exon AS
FROM fixture_exons(
  'ensembl_grch38', 'GRCh38', '132907,2006446537'
);
CREATE TABLE duckvep_grch38_core.transcript_attrib AS
FROM fixture_transcript_attributes(
  'ensembl_grch38', 'GRCh38', '132907,2006446537'
);
CREATE TABLE duckvep_grch38_core.translation_attrib AS
FROM fixture_translation_attributes(
  'ensembl_grch38', 'GRCh38', '132907,2006446537'
);
CREATE TABLE duckvep_grch38_core.attrib_type AS
FROM fixture_attribute_types(
  'ensembl_grch38', 'GRCh38', '132907,2006446537'
);
CREATE TABLE duckvep_grch38_reference AS
FROM fixture_reference_chunks(
  'ensembl_grch38', 'GRCh38', '132907,2006446537'
);

CREATE OR REPLACE SCHEMA duckvep_grch37_core;
CREATE TABLE duckvep_grch37_core.coord_system AS
FROM fixture_coord_systems(
  'ensembl_grch37', 'GRCh37', '100965601,27742'
);
CREATE TABLE duckvep_grch37_core.seq_region AS
FROM fixture_seq_regions(
  'ensembl_grch37', 'GRCh37', '100965601,27742'
);
CREATE TABLE duckvep_grch37_core.transcript AS
FROM fixture_transcripts(
  'ensembl_grch37', 'GRCh37', '100965601,27742'
);
CREATE TABLE duckvep_grch37_core.gene AS
FROM fixture_genes(
  'ensembl_grch37', 'GRCh37', '100965601,27742'
);
CREATE TABLE duckvep_grch37_core.translation AS
FROM fixture_translations(
  'ensembl_grch37', 'GRCh37', '100965601,27742'
);
CREATE TABLE duckvep_grch37_core.exon_transcript AS
FROM fixture_exon_transcripts(
  'ensembl_grch37', 'GRCh37', '100965601,27742'
);
CREATE TABLE duckvep_grch37_core.exon AS
FROM fixture_exons(
  'ensembl_grch37', 'GRCh37', '100965601,27742'
);
CREATE TABLE duckvep_grch37_core.transcript_attrib AS
FROM fixture_transcript_attributes(
  'ensembl_grch37', 'GRCh37', '100965601,27742'
);
CREATE TABLE duckvep_grch37_core.translation_attrib AS
FROM fixture_translation_attributes(
  'ensembl_grch37', 'GRCh37', '100965601,27742'
);
CREATE TABLE duckvep_grch37_core.attrib_type AS
FROM fixture_attribute_types(
  'ensembl_grch37', 'GRCh37', '100965601,27742'
);
CREATE TABLE duckvep_grch37_reference AS
FROM fixture_reference_chunks(
  'ensembl_grch37', 'GRCh37', '100965601,27742'
);

CREATE TEMP TABLE expected_reference_hashes(
  assembly VARCHAR,
  chrom VARCHAR,
  sequence_sha256 VARCHAR
);
INSERT INTO expected_reference_hashes VALUES
  ('GRCh38', 'HG2047_PATCH', 'b152b01e05646a9885bdc54b27a26f8670af6ac9c3a2661557e55d645d8ef0cc'),
  ('GRCh38', 'MT', 'f156ff3f65bbcc80c7ebb9936dceb96b1477b4f8f535c4e1dbe7baea225cbc66'),
  ('GRCh37', 'GL000201.1', '6420ef160c59c951482006f60f0864df79baa83ff63d9669ca9abf3ee415c197'),
  ('GRCh37', 'MT', 'f156ff3f65bbcc80c7ebb9936dceb96b1477b4f8f535c4e1dbe7baea225cbc66');

CREATE TEMP TABLE actual_reference_hashes AS
SELECT 'GRCh38'::VARCHAR AS assembly, chrom, sha256(seq) AS sequence_sha256
FROM duckvep_grch38_reference
UNION ALL
SELECT 'GRCh37'::VARCHAR AS assembly, chrom, sha256(seq) AS sequence_sha256
FROM duckvep_grch37_reference;

CREATE TEMP TABLE reference_validation AS
WITH reference_difference AS (
  SELECT *
  FROM actual_reference_hashes
  FULL OUTER JOIN expected_reference_hashes
    USING (assembly, chrom, sequence_sha256)
  WHERE actual_reference_hashes.chrom IS NULL
     OR expected_reference_hashes.chrom IS NULL
)
SELECT CASE
  WHEN EXISTS (SELECT 1 FROM reference_difference) THEN
    error('DuckVEP fixture: assembled reference does not match the pinned sequence hashes')
  ELSE true
END AS valid;

CREATE TEMP TABLE source_manifests AS
SELECT 'GRCh38'::VARCHAR AS assembly, sha256(content) AS source_manifest_sha256
FROM read_text(
  'https://ftp.ensembl.org/pub/release-116/mysql/homo_sapiens_core_116_38/CHECKSUMS'
)
UNION ALL
SELECT 'GRCh37'::VARCHAR AS assembly, sha256(content) AS source_manifest_sha256
FROM read_text(
  'https://ftp.ensembl.org/pub/grch37/release-116/mysql/homo_sapiens_core_116_37/CHECKSUMS'
);

CREATE TEMP TABLE source_validation AS
SELECT CASE
  WHEN max(source_manifest_sha256) FILTER (WHERE assembly = 'GRCh38') !=
       'c5afd2ec87ad7aa402035a070e81068e7bd8f4c5bdae9f3c698f8ff5b2fd63d1' OR
       max(source_manifest_sha256) FILTER (WHERE assembly = 'GRCh37') !=
       '46b8011c01b2628aa2678be753aa854528cf026c8d396afa54e99799d38ff8f0' THEN
    error('DuckVEP fixture: an Ensembl dump manifest does not match the pinned hash')
  ELSE true
END AS valid
FROM source_manifests;

CREATE TEMP TABLE grch38_regions AS
SELECT * FROM duckvep_ensembl_regions(
  'duckvep_grch38_core', 'duckvep_grch38_reference', 'GRCh38'
);
CREATE TEMP TABLE grch38_transcripts AS
SELECT * FROM duckvep_ensembl_transcripts(
  'duckvep_grch38_core', 'duckvep_grch38_reference', 'GRCh38'
);
CREATE TEMP TABLE grch37_regions AS
SELECT * FROM duckvep_ensembl_regions(
  'duckvep_grch37_core', 'duckvep_grch37_reference', 'GRCh37'
);
CREATE TEMP TABLE grch37_transcripts AS
SELECT * FROM duckvep_ensembl_transcripts(
  'duckvep_grch37_core', 'duckvep_grch37_reference', 'GRCh37'
);

CREATE TEMP TABLE fixture_receipts AS
SELECT * FROM duckvep_model_receipt(
  'grch38_regions', 'grch38_transcripts', 'Ensembl', '116', 'GRCh38',
  'c5afd2ec87ad7aa402035a070e81068e7bd8f4c5bdae9f3c698f8ff5b2fd63d1',
  '1ccef37f83c0784ebe74d89e4e0db8171e11625cd4b46ff47a9911e0ff471763',
  'all current transcripts on MT and HG2047_PATCH'
)
UNION ALL
SELECT * FROM duckvep_model_receipt(
  'grch37_regions', 'grch37_transcripts', 'Ensembl', '116', 'GRCh37',
  '46b8011c01b2628aa2678be753aa854528cf026c8d396afa54e99799d38ff8f0',
  '09de6c8bdf434f9fe29e9803c7403b145452fdeefb5387292af67811f1523d80',
  'all current transcripts on MT and GL000201.1'
);

CREATE TEMP TABLE fixture_acceptance AS
SELECT CASE
  WHEN count(*) != 2
    OR max(model_sha256) FILTER (WHERE assembly = 'GRCh38') !=
       '68bdd72a53027af9bc09d1853848b2609c34dee3b533fdcfa4eec5881dd12db0'
    OR max(model_sha256) FILTER (WHERE assembly = 'GRCh37') !=
       '88e3e65c1055c5a38ac821205c3f7a246042a781d73c487607d0c32c845db5cf'
    OR max(reference_base_count) FILTER (WHERE assembly = 'GRCh38') != 81258
    OR max(transcript_count) FILTER (WHERE assembly = 'GRCh38') != 39
    OR max(gene_count) FILTER (WHERE assembly = 'GRCh38') != 39
    OR max(coding_transcript_count) FILTER (WHERE assembly = 'GRCh38') != 14
    OR max(sequence_backed_transcript_count) FILTER (WHERE assembly = 'GRCh38') != 1
    OR max(sequence_withheld_transcript_count) FILTER (WHERE assembly = 'GRCh38') != 13
    OR max(exon_membership_count) FILTER (WHERE assembly = 'GRCh38') != 40
    OR max(cds_base_count) FILTER (WHERE assembly = 'GRCh38') != 906
    OR max(reference_base_count) FILTER (WHERE assembly = 'GRCh37') != 52717
    OR max(transcript_count) FILTER (WHERE assembly = 'GRCh37') != 39
    OR max(gene_count) FILTER (WHERE assembly = 'GRCh37') != 39
    OR max(coding_transcript_count) FILTER (WHERE assembly = 'GRCh37') != 15
    OR max(sequence_backed_transcript_count) FILTER (WHERE assembly = 'GRCh37') != 2
    OR max(sequence_withheld_transcript_count) FILTER (WHERE assembly = 'GRCh37') != 13
    OR max(exon_membership_count) FILTER (WHERE assembly = 'GRCh37') != 48
    OR max(cds_base_count) FILTER (WHERE assembly = 'GRCh37') != 1665
    OR min(region_count) != 2 OR max(region_count) != 2 THEN
    error('DuckVEP fixture: prepared models do not match the declared assemblies')
  ELSE true
END AS valid
FROM fixture_receipts;

CREATE TEMP TABLE fixture_manifest AS
SELECT
  r.source_name,
  r.source_version,
  CASE r.assembly
    WHEN 'GRCh38' THEN 'homo_sapiens_core_116_38'
    ELSE 'homo_sapiens_core_116_37'
  END AS source_database,
  r.assembly,
  CASE r.assembly
    WHEN 'GRCh38' THEN 'MT,HG2047_PATCH'
    ELSE 'MT,GL000201.1'
  END AS sequence_regions,
  r.source_manifest_sha256,
  r.reference_sha256,
  r.model_sha256,
  r.region_count,
  r.reference_base_count,
  r.transcript_count,
  r.gene_count,
  r.coding_transcript_count,
  r.sequence_backed_transcript_count,
  r.sequence_withheld_transcript_count,
  r.exon_membership_count,
  r.cds_base_count,
  r.transcript_filter
FROM fixture_receipts r
CROSS JOIN reference_validation
CROSS JOIN source_validation
CROSS JOIN fixture_acceptance
WHERE reference_validation.valid
  AND source_validation.valid
  AND fixture_acceptance.valid;

COPY duckvep_grch38_core.attrib_type TO
  'test/data/duckvep/ensembl_core/grch38/attrib_type.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch38_core.coord_system TO
  'test/data/duckvep/ensembl_core/grch38/coord_system.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch38_core.exon TO
  'test/data/duckvep/ensembl_core/grch38/exon.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch38_core.exon_transcript TO
  'test/data/duckvep/ensembl_core/grch38/exon_transcript.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch38_core.gene TO
  'test/data/duckvep/ensembl_core/grch38/gene.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch38_reference TO
  'test/data/duckvep/ensembl_core/grch38/reference_chunks.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch38_core.seq_region TO
  'test/data/duckvep/ensembl_core/grch38/seq_region.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch38_core.transcript TO
  'test/data/duckvep/ensembl_core/grch38/transcript.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch38_core.transcript_attrib TO
  'test/data/duckvep/ensembl_core/grch38/transcript_attrib.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch38_core.translation TO
  'test/data/duckvep/ensembl_core/grch38/translation.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch38_core.translation_attrib TO
  'test/data/duckvep/ensembl_core/grch38/translation_attrib.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT * FROM fixture_manifest WHERE assembly = 'GRCh38') TO
  'test/data/duckvep/ensembl_core/grch38/manifest.tsv'
  (FORMAT CSV, DELIMITER '\t', HEADER);

COPY duckvep_grch37_core.attrib_type TO
  'test/data/duckvep/ensembl_core/grch37/attrib_type.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch37_core.coord_system TO
  'test/data/duckvep/ensembl_core/grch37/coord_system.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch37_core.exon TO
  'test/data/duckvep/ensembl_core/grch37/exon.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch37_core.exon_transcript TO
  'test/data/duckvep/ensembl_core/grch37/exon_transcript.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch37_core.gene TO
  'test/data/duckvep/ensembl_core/grch37/gene.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch37_reference TO
  'test/data/duckvep/ensembl_core/grch37/reference_chunks.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch37_core.seq_region TO
  'test/data/duckvep/ensembl_core/grch37/seq_region.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch37_core.transcript TO
  'test/data/duckvep/ensembl_core/grch37/transcript.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch37_core.transcript_attrib TO
  'test/data/duckvep/ensembl_core/grch37/transcript_attrib.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch37_core.translation TO
  'test/data/duckvep/ensembl_core/grch37/translation.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY duckvep_grch37_core.translation_attrib TO
  'test/data/duckvep/ensembl_core/grch37/translation_attrib.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT * FROM fixture_manifest WHERE assembly = 'GRCh37') TO
  'test/data/duckvep/ensembl_core/grch37/manifest.tsv'
  (FORMAT CSV, DELIMITER '\t', HEADER);
