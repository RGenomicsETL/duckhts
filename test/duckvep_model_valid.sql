-- Tiny deterministic annotation-model candidate.
-- Loaded after scripts/duckvep_model_schema.sql.

INSERT INTO duckvep_model.model_manifest VALUES (
    repeat('1', 64), 2, 'ensembl-vep/116.0',
    '57ea5c52340acc1f156267f810ad162e26597082',
    'homo_sapiens', 9606, 'GRCh38.p14', 'GCA_000001405.29',
    'fixture-selection-v1', 'draft-unavailable',
    TIMESTAMP '2026-07-10 00:00:00'
);

INSERT INTO duckvep_model.model_source VALUES (
    repeat('1', 64), 'ensembl_core_116', 'Ensembl', 'core', 'primary',
    'GENCODE 50', '116', 'homo_sapiens_core_116_38',
    'https://github.com/Ensembl/ensembl',
    'c0cf13daa961d80584bad797b2eb0ff3a7500ef3',
    'ec7bb8dd2fcd6a7012bfd67aff8065b912915d541450fa4f3f05334a92944e8c', 0
);

INSERT INTO duckvep_model.model_selection_audit VALUES
    (repeat('1', 64), 'ensembl_core_116', 'seq_region', 'included', 'top_level', 1, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'gene', 'included', 'fixture_gene', 1, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'transcript', 'included', 'fixture_transcript', 1, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'exon', 'included', 'fixture_exons', 2, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'transcript_exon', 'included', 'fixture_memberships', 2, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'translation', 'included', 'fixture_translation', 1, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'attribute_type', 'included', 'fixture_attribute_types', 6, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'seq_region_attribute', 'included', 'fixture_seq_region_attribute', 1, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'gene_attribute', 'included', 'fixture_gene_attribute', 1, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'transcript_attribute', 'included', 'fixture_transcript_attributes', 3, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'translation_attribute', 'included', 'fixture_translation_attribute', 1, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'external_db', 'included', 'fixture_external_db', 1, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'xref', 'included', 'fixture_xref', 1, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'object_xref', 'included', 'fixture_object_xref', 1, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'seq_region_synonym', 'included', 'fixture_synonym', 1, NULL),
    (repeat('1', 64), 'ensembl_core_116', 'transcript_attribute', 'included',
     'translation_boundary_overrides_absent', 0, NULL);

INSERT INTO duckvep_model.model_build VALUES (
    repeat('1', 64), repeat('2', 64), 'https://example.invalid/duckvep-importer',
    repeat('a', 40), repeat('3', 64), 'https://example.invalid/duckvep-exporter',
    repeat('b', 40), repeat('4', 64), repeat('5', 64), repeat('6', 64),
    TIMESTAMP '2026-07-10 00:00:00'
);

INSERT INTO duckvep_model.model_artifact VALUES
    (repeat('1',64),repeat('2',64),'source_input','core_ddl','ensembl_core_116',
     'ensembl-116-table.sql','fixture://ensembl-116-table.sql',
     'ec7bb8dd2fcd6a7012bfd67aff8065b912915d541450fa4f3f05334a92944e8c',1,1,true),
    (repeat('1',64),repeat('2',64),'source_input','reference_fasta','ensembl_core_116',
     'fixture.fa','fixture://fixture.fa',repeat('8',64),17,1,true),
    (repeat('1',64),repeat('2',64),'source_input','reference_fasta_index','ensembl_core_116',
     'fixture.fa.fai','fixture://fixture.fa.fai',repeat('9',64),10,1,true),
    (repeat('1',64),repeat('2',64),'software_input','importer',NULL,
     'duckvep-model-importer','fixture://importer',repeat('3',64),1,NULL,true),
    (repeat('1',64),repeat('2',64),'software_input','sequence_state_exporter',NULL,
     'sequence-state-exporter','fixture://sequence-state-exporter',repeat('4',64),1,NULL,true),
    (repeat('1',64),repeat('2',64),'software_input','vep_api_modules',NULL,
     'vep-116-api-modules','fixture://vep-api-modules',repeat('c',64),1,NULL,true);

INSERT INTO duckvep_model.model_seq_region VALUES (
    repeat('1',64),'ensembl_core_116',1,'1','1',100,'chromosome','GRCh38',
    'GCA_000001405.29','NC_000001','11',repeat('e',64),false
);

INSERT INTO duckvep_model.model_gene VALUES (
    repeat('1',64),'ensembl_core_116',1,'1','ENSG00000000001',1,1,
    1,17,1,'ensembl','protein_coding',true,NULL,1
);

INSERT INTO duckvep_model.model_transcript VALUES (
    repeat('1',64),'ensembl_core_116',1,'1','ENST00000000001',1,1,1,
    1,17,1,'ensembl','protein_coding',true,1,1,true
);

INSERT INTO duckvep_model.model_exon VALUES
    (repeat('1',64),'ensembl_core_116',1,'1','ENSE00000000001',1,1,1,9,1,0,0),
    (repeat('1',64),'ensembl_core_116',2,'2','ENSE00000000002',1,1,12,17,1,0,0);

INSERT INTO duckvep_model.model_transcript_exon VALUES
    (repeat('1',64),'ensembl_core_116',1,1,1,1,9),
    (repeat('1',64),'ensembl_core_116',1,2,2,10,15);

INSERT INTO duckvep_model.model_translation VALUES (
    repeat('1',64),'ensembl_core_116',1,'1',1,'ENSP00000000001',1,
    1,1,2,6,1,15,1,15,0,1,true
);

INSERT INTO duckvep_model.model_attribute_type VALUES
    (repeat('1',64),'ensembl_core_116',1,'1','_rna_edit','RNA edit','RNA edit'),
    (repeat('1',64),'ensembl_core_116',2,'2','_selenocysteine','Selenocysteine','Peptide edit'),
    (repeat('1',64),'ensembl_core_116',3,'3','miRNA','Mature miRNA','Mature miRNA span'),
    (repeat('1',64),'ensembl_core_116',4,'4','MANE_Select','MANE Select',NULL),
    (repeat('1',64),'ensembl_core_116',5,'5','codon_table','Codon table',NULL),
    (repeat('1',64),'ensembl_core_116',6,'6','description','Description',NULL);

INSERT INTO duckvep_model.model_seq_region_attribute VALUES (
    repeat('1',64),'ensembl_core_116',1,1,5,1,'fixture:seq_region_attrib:1','1'
);

INSERT INTO duckvep_model.model_gene_attribute VALUES (
    repeat('1',64),'ensembl_core_116',1,1,6,1,'fixture:gene_attrib:1','fixture gene'
);

INSERT INTO duckvep_model.model_transcript_attribute VALUES
    (repeat('1',64),'ensembl_core_116',1,1,1,1,'fixture:transcript_attrib:1','4 4 G'),
    (repeat('1',64),'ensembl_core_116',2,1,3,1,'fixture:transcript_attrib:2','2 5'),
    (repeat('1',64),'ensembl_core_116',3,1,4,1,'fixture:transcript_attrib:3','NM_000001.1');

INSERT INTO duckvep_model.model_translation_attribute VALUES (
    repeat('1',64),'ensembl_core_116',1,1,2,1,'fixture:translation_attrib:1','3 3 U'
);

INSERT INTO duckvep_model.model_transcript_edit VALUES (
    repeat('1',64),'ensembl_core_116',1,1,1,'_rna_edit','raw_spliced',4,4,
    CAST('A' AS BLOB),CAST('A' AS BLOB),CAST('G' AS BLOB),1,'applied'
);

INSERT INTO duckvep_model.model_translation_edit VALUES (
    repeat('1',64),'ensembl_core_116',1,1,1,'_selenocysteine','peptide_unedited',3,3,
    CAST('*' AS BLOB),CAST('*' AS BLOB),CAST('U' AS BLOB),1,'applied'
);

INSERT INTO duckvep_model.model_mature_mirna VALUES (
    repeat('1',64),'ensembl_core_116',1,1,2,5,'raw_spliced',2
);

INSERT INTO duckvep_model.model_external_db VALUES (
    repeat('1',64),'ensembl_core_116',1,'1','RefSeq','116','KNOWN',1,
    'RefSeq','MISC',NULL,NULL,'RefSeq transcript'
);

INSERT INTO duckvep_model.model_xref VALUES (
    repeat('1',64),'ensembl_core_116',1,'1',1,'NM_000001','1','NM_000001.1',
    'fixture transcript','DIRECT','fixture'
);

INSERT INTO duckvep_model.model_object_xref VALUES (
    repeat('1',64),'ensembl_core_116','ensembl_core_116',1,'1',
    'transcript',1,1,true,NULL,NULL
);

INSERT INTO duckvep_model.model_seq_region_synonym VALUES (
    repeat('1',64),'ensembl_core_116',1,'1',1,'chr1',NULL
);

INSERT INTO duckvep_model.model_sequence_blob
SELECT sha256(CAST(sequence AS BLOB)), alphabet, length(sequence), CAST(sequence AS BLOB)
FROM (VALUES
    ('ATGAAATGAGGGTTT','dna'),
    ('ATGGAATGAGGGTTT','dna'),
    ('ME*GF','peptide'),
    ('MEUGF','peptide')
) AS fixture(sequence,alphabet);

INSERT INTO duckvep_model.model_sequence_state VALUES
    (repeat('1',64),'ensembl_core_116','transcript',1,'raw_spliced','dna',
     sha256(CAST('ATGAAATGAGGGTTT' AS BLOB)),15,'present'),
    (repeat('1',64),'ensembl_core_116','transcript',1,'edited_spliced','dna',
     sha256(CAST('ATGGAATGAGGGTTT' AS BLOB)),15,'present'),
    (repeat('1',64),'ensembl_core_116','translation',1,'translatable_cds','dna',
     sha256(CAST('ATGGAATGAGGGTTT' AS BLOB)),15,'present'),
    (repeat('1',64),'ensembl_core_116','translation',1,'peptide_unedited','peptide',
     sha256(CAST('ME*GF' AS BLOB)),5,'present'),
    (repeat('1',64),'ensembl_core_116','translation',1,'peptide_final','peptide',
     sha256(CAST('MEUGF' AS BLOB)),5,'present');

INSERT INTO duckvep_model.model_xref_identity VALUES (
    repeat('1',64),'ensembl_core_116','ensembl_core_116',1,'exact',
    'raw_spliced',sha256(CAST('ATGAAATGAGGGTTT' AS BLOB))
);
