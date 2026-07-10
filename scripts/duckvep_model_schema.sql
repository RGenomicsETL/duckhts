-- DuckVEP normalized annotation-model relation-pack schema (format version 2).
--
-- This private annotation-model schema implements the persistent relations in
-- design/duckvep_model.md.  It is an interchange contract for DuckDB
-- tables or equivalently-shaped Parquet relations, not a runtime ABI and not
-- a new cache format.  Source IDs remain opaque strings, stable IDs are
-- nullable metadata, coordinates are source-faithful 1-based inclusive
-- BIGINTs, and sequence/edit payloads have no fixed-size transport limit.
--
-- SQL constraints make an initialized DuckDB convenient to populate.  They
-- are not the conformance gate: Parquet bypasses them, so the companion
-- validator rechecks shape, enums, keys, cardinalities, references, identity
-- spelling/linkage, edit replay, and sequence-state coverage explicitly.

CREATE SCHEMA duckvep_model;

CREATE TABLE duckvep_model.model_manifest (
    model_id                    VARCHAR PRIMARY KEY,
    format_version              INTEGER NOT NULL CHECK (format_version = 2),
    compatibility_target        VARCHAR NOT NULL,
    compatibility_commit        VARCHAR NOT NULL,
    species                     VARCHAR NOT NULL,
    taxonomy_id                 UBIGINT NOT NULL CHECK (taxonomy_id > 0),
    assembly_name               VARCHAR NOT NULL,
    assembly_accession          VARCHAR NOT NULL,
    selection_policy_version    VARCHAR NOT NULL,
    canonicalization_version    VARCHAR NOT NULL,
    created_at                  TIMESTAMP NOT NULL
);

CREATE TABLE duckvep_model.model_source (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    provider                    VARCHAR NOT NULL,
    source_kind                 VARCHAR NOT NULL,
    selection_role              VARCHAR NOT NULL,
    annotation_set              VARCHAR NOT NULL,
    annotation_source_release   VARCHAR NOT NULL,
    source_database             VARCHAR,
    schema_repository           VARCHAR,
    schema_commit               VARCHAR,
    schema_sha256               VARCHAR,
    source_priority             INTEGER NOT NULL CHECK (source_priority >= 0),
    PRIMARY KEY (model_id, source_namespace)
);

CREATE TABLE duckvep_model.model_selection_audit (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    object_kind                 VARCHAR NOT NULL,
    decision                    VARCHAR NOT NULL CHECK (
        decision IN ('included', 'excluded', 'quarantined')
    ),
    reason_code                 VARCHAR NOT NULL,
    row_count                   UBIGINT NOT NULL,
    rejected_rows_artifact      VARCHAR,
    PRIMARY KEY (
        model_id, source_namespace, object_kind, decision, reason_code
    )
);

CREATE TABLE duckvep_model.model_build (
    model_id                    VARCHAR NOT NULL,
    model_build_id              VARCHAR NOT NULL,
    importer_repository         VARCHAR NOT NULL,
    importer_commit             VARCHAR NOT NULL,
    importer_sha256             VARCHAR NOT NULL,
    exporter_repository         VARCHAR NOT NULL,
    exporter_commit             VARCHAR NOT NULL,
    exporter_sha256             VARCHAR NOT NULL,
    invocation_sha256           VARCHAR NOT NULL,
    environment_sha256          VARCHAR NOT NULL,
    created_at                  TIMESTAMP NOT NULL,
    PRIMARY KEY (model_id, model_build_id)
);

CREATE TABLE duckvep_model.model_artifact (
    model_id                    VARCHAR NOT NULL,
    model_build_id              VARCHAR NOT NULL,
    role_class                  VARCHAR NOT NULL CHECK (
        role_class IN (
            'source_input', 'software_input', 'generated_output',
            'audit_output', 'oracle_output'
        )
    ),
    role                        VARCHAR NOT NULL,
    source_namespace            VARCHAR,
    logical_name                VARCHAR NOT NULL,
    locator                     VARCHAR NOT NULL,
    sha256                      VARCHAR NOT NULL,
    byte_count                  BIGINT NOT NULL CHECK (byte_count >= 0),
    row_count                   UBIGINT,
    required                    BOOLEAN NOT NULL,
    PRIMARY KEY (
        model_id, model_build_id, role_class, role, logical_name
    )
);

-- The physical DuckDB file or Parquet bundle is attested outside the pack.
-- Recording its checksum in model_artifact would make the enclosing artifact
-- hash itself.  generated_output rows are therefore for non-enclosing exports,
-- never for the relation pack that contains this table.

CREATE TABLE duckvep_model.model_seq_region (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    seq_region_key              UBIGINT NOT NULL CHECK (seq_region_key > 0),
    source_internal_id          VARCHAR NOT NULL,
    name                        VARCHAR NOT NULL,
    length                      BIGINT NOT NULL CHECK (length > 0),
    coord_system_name           VARCHAR NOT NULL,
    coord_system_version        VARCHAR,
    assembly_accession          VARCHAR NOT NULL,
    reference_accession         VARCHAR NOT NULL,
    reference_version           VARCHAR NOT NULL,
    sequence_sha256             VARCHAR NOT NULL,
    is_circular                 BOOLEAN NOT NULL,
    PRIMARY KEY (model_id, source_namespace, seq_region_key),
    UNIQUE (model_id, source_namespace, source_internal_id)
);

CREATE TABLE duckvep_model.model_gene (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    gene_key                    UBIGINT NOT NULL CHECK (gene_key > 0),
    source_internal_id          VARCHAR NOT NULL,
    stable_id                   VARCHAR,
    stable_id_version           INTEGER,
    seq_region_key              UBIGINT NOT NULL,
    start1                      BIGINT NOT NULL,
    end1                        BIGINT NOT NULL,
    strand                      INTEGER NOT NULL CHECK (strand IN (-1, 1)),
    source                      VARCHAR NOT NULL,
    biotype                     VARCHAR NOT NULL,
    is_current                  BOOLEAN NOT NULL,
    display_xref_key            UBIGINT,
    canonical_transcript_key    UBIGINT,
    PRIMARY KEY (model_id, source_namespace, gene_key),
    UNIQUE (model_id, source_namespace, source_internal_id),
    CHECK (start1 >= 1 AND end1 >= start1)
);

CREATE TABLE duckvep_model.model_transcript (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    transcript_key              UBIGINT NOT NULL CHECK (transcript_key > 0),
    source_internal_id          VARCHAR NOT NULL,
    stable_id                   VARCHAR,
    stable_id_version           INTEGER,
    gene_key                    UBIGINT,
    seq_region_key              UBIGINT NOT NULL,
    start1                      BIGINT NOT NULL,
    end1                        BIGINT NOT NULL,
    strand                      INTEGER NOT NULL CHECK (strand IN (-1, 1)),
    source                      VARCHAR NOT NULL,
    biotype                     VARCHAR NOT NULL,
    is_current                  BOOLEAN NOT NULL,
    display_xref_key            UBIGINT,
    canonical_translation_key   UBIGINT,
    is_canonical                BOOLEAN NOT NULL,
    PRIMARY KEY (model_id, source_namespace, transcript_key),
    UNIQUE (model_id, source_namespace, source_internal_id),
    CHECK (start1 >= 1 AND end1 >= start1)
);

CREATE TABLE duckvep_model.model_exon (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    exon_key                    UBIGINT NOT NULL CHECK (exon_key > 0),
    source_internal_id          VARCHAR NOT NULL,
    stable_id                   VARCHAR,
    stable_id_version           INTEGER,
    seq_region_key              UBIGINT NOT NULL,
    start1                      BIGINT NOT NULL,
    end1                        BIGINT NOT NULL,
    strand                      INTEGER NOT NULL CHECK (strand IN (-1, 1)),
    phase                       INTEGER NOT NULL CHECK (phase IN (-1, 0, 1, 2)),
    end_phase                   INTEGER NOT NULL CHECK (end_phase IN (-1, 0, 1, 2)),
    PRIMARY KEY (model_id, source_namespace, exon_key),
    UNIQUE (model_id, source_namespace, source_internal_id),
    CHECK (start1 >= 1 AND end1 >= start1)
);

CREATE TABLE duckvep_model.model_transcript_exon (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    transcript_key              UBIGINT NOT NULL,
    exon_key                    UBIGINT NOT NULL,
    exon_rank                   INTEGER NOT NULL CHECK (exon_rank > 0),
    raw_cdna_start1             BIGINT NOT NULL,
    raw_cdna_end1               BIGINT NOT NULL,
    PRIMARY KEY (
        model_id, source_namespace, transcript_key, exon_key, exon_rank
    ),
    UNIQUE (
        model_id, source_namespace, transcript_key, exon_rank
    ),
    CHECK (raw_cdna_start1 >= 1 AND raw_cdna_end1 >= raw_cdna_start1)
);

CREATE TABLE duckvep_model.model_translation (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    translation_key             UBIGINT NOT NULL CHECK (translation_key > 0),
    source_internal_id          VARCHAR NOT NULL,
    transcript_key              UBIGINT NOT NULL,
    stable_id                   VARCHAR,
    stable_id_version           INTEGER,
    start_exon_key              UBIGINT NOT NULL,
    start_offset1               BIGINT NOT NULL,
    end_exon_key                UBIGINT NOT NULL,
    end_offset1                 BIGINT NOT NULL,
    raw_cdna_coding_start1      BIGINT NOT NULL,
    raw_cdna_coding_end1        BIGINT NOT NULL,
    edited_cdna_coding_start1   BIGINT NOT NULL,
    edited_cdna_coding_end1     BIGINT NOT NULL,
    start_phase_padding         INTEGER NOT NULL,
    codon_table                 INTEGER NOT NULL,
    is_canonical                BOOLEAN NOT NULL,
    PRIMARY KEY (model_id, source_namespace, translation_key),
    UNIQUE (model_id, source_namespace, source_internal_id),
    CHECK (start_offset1 >= 1 AND end_offset1 >= 1),
    CHECK (
        raw_cdna_coding_start1 >= 1 AND
        raw_cdna_coding_end1 >= raw_cdna_coding_start1 AND
        edited_cdna_coding_start1 >= 1 AND
        edited_cdna_coding_end1 >= edited_cdna_coding_start1
    ),
    CHECK (start_phase_padding BETWEEN 0 AND 2),
    CHECK (codon_table IN (1, 2))
);

CREATE TABLE duckvep_model.model_attribute_type (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    attrib_type_key             UBIGINT NOT NULL CHECK (attrib_type_key > 0),
    source_internal_id          VARCHAR NOT NULL,
    code                        VARCHAR NOT NULL,
    name                        VARCHAR NOT NULL,
    description                 VARCHAR,
    PRIMARY KEY (model_id, source_namespace, attrib_type_key),
    UNIQUE (model_id, source_namespace, source_internal_id)
);

CREATE TABLE duckvep_model.model_seq_region_attribute (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    attribute_key               UBIGINT NOT NULL CHECK (attribute_key > 0),
    seq_region_key              UBIGINT NOT NULL,
    attrib_type_key             UBIGINT NOT NULL,
    duplicate_ordinal           INTEGER NOT NULL CHECK (duplicate_ordinal > 0),
    source_row_locator          VARCHAR NOT NULL,
    value                       VARCHAR NOT NULL,
    PRIMARY KEY (
        model_id, source_namespace, seq_region_key, attrib_type_key,
        value, duplicate_ordinal
    ),
    UNIQUE (model_id, source_namespace, attribute_key)
);

CREATE TABLE duckvep_model.model_gene_attribute (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    attribute_key               UBIGINT NOT NULL CHECK (attribute_key > 0),
    gene_key                    UBIGINT NOT NULL,
    attrib_type_key             UBIGINT NOT NULL,
    duplicate_ordinal           INTEGER NOT NULL CHECK (duplicate_ordinal > 0),
    source_row_locator          VARCHAR NOT NULL,
    value                       VARCHAR NOT NULL,
    PRIMARY KEY (
        model_id, source_namespace, gene_key, attrib_type_key,
        value, duplicate_ordinal
    ),
    UNIQUE (model_id, source_namespace, attribute_key)
);

CREATE TABLE duckvep_model.model_transcript_attribute (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    attribute_key               UBIGINT NOT NULL CHECK (attribute_key > 0),
    transcript_key              UBIGINT NOT NULL,
    attrib_type_key             UBIGINT NOT NULL,
    duplicate_ordinal           INTEGER NOT NULL CHECK (duplicate_ordinal > 0),
    source_row_locator          VARCHAR NOT NULL,
    value                       VARCHAR NOT NULL,
    PRIMARY KEY (
        model_id, source_namespace, transcript_key, attrib_type_key,
        value, duplicate_ordinal
    ),
    UNIQUE (model_id, source_namespace, attribute_key)
);

CREATE TABLE duckvep_model.model_translation_attribute (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    attribute_key               UBIGINT NOT NULL CHECK (attribute_key > 0),
    translation_key             UBIGINT NOT NULL,
    attrib_type_key             UBIGINT NOT NULL,
    duplicate_ordinal           INTEGER NOT NULL CHECK (duplicate_ordinal > 0),
    source_row_locator          VARCHAR NOT NULL,
    value                       VARCHAR NOT NULL,
    PRIMARY KEY (
        model_id, source_namespace, translation_key, attrib_type_key,
        value, duplicate_ordinal
    ),
    UNIQUE (model_id, source_namespace, attribute_key)
);

CREATE TABLE duckvep_model.model_transcript_edit (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    edit_key                    UBIGINT NOT NULL CHECK (edit_key > 0),
    transcript_key              UBIGINT NOT NULL,
    attribute_key               UBIGINT NOT NULL,
    code                        VARCHAR NOT NULL,
    coordinate_basis            VARCHAR NOT NULL,
    start1                      BIGINT NOT NULL,
    end1                        BIGINT NOT NULL,
    basis_ref_seq               BLOB NOT NULL,
    preapply_ref_seq            BLOB NOT NULL,
    alt_seq                     BLOB NOT NULL,
    apply_ordinal               INTEGER,
    status                      VARCHAR NOT NULL CHECK (
        status IN (
            'applied', 'source_ignored',
            'quarantined_ambiguous_order', 'unsupported_code'
        )
    ),
    PRIMARY KEY (model_id, source_namespace, edit_key),
    CHECK (coordinate_basis = 'raw_spliced'),
    CHECK (start1 >= 1 AND end1 >= 0 AND start1 - 1 <= end1),
    CHECK (octet_length(basis_ref_seq) = greatest(end1 - start1 + 1, 0)),
    CHECK (
        (status = 'applied' AND apply_ordinal IS NOT NULL AND apply_ordinal > 0) OR
        (status <> 'applied' AND apply_ordinal IS NULL)
    )
);

CREATE TABLE duckvep_model.model_translation_edit (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    edit_key                    UBIGINT NOT NULL CHECK (edit_key > 0),
    translation_key             UBIGINT NOT NULL,
    attribute_key               UBIGINT NOT NULL,
    code                        VARCHAR NOT NULL,
    coordinate_basis            VARCHAR NOT NULL,
    start1                      BIGINT NOT NULL,
    end1                        BIGINT NOT NULL,
    basis_ref_seq               BLOB NOT NULL,
    preapply_ref_seq            BLOB NOT NULL,
    alt_seq                     BLOB NOT NULL,
    apply_ordinal               INTEGER,
    status                      VARCHAR NOT NULL CHECK (
        status IN (
            'applied', 'source_ignored',
            'quarantined_ambiguous_order', 'unsupported_code'
        )
    ),
    PRIMARY KEY (model_id, source_namespace, edit_key),
    CHECK (coordinate_basis = 'peptide_unedited'),
    CHECK (start1 >= 1 AND end1 >= 0 AND start1 - 1 <= end1),
    CHECK (octet_length(basis_ref_seq) = greatest(end1 - start1 + 1, 0)),
    CHECK (
        (status = 'applied' AND apply_ordinal IS NOT NULL AND apply_ordinal > 0) OR
        (status <> 'applied' AND apply_ordinal IS NULL)
    )
);

CREATE TABLE duckvep_model.model_mature_mirna (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    feature_key                 UBIGINT NOT NULL CHECK (feature_key > 0),
    transcript_key              UBIGINT NOT NULL,
    cdna_start1                 BIGINT NOT NULL,
    cdna_end1                   BIGINT NOT NULL,
    coordinate_basis            VARCHAR NOT NULL,
    source_attribute_key        UBIGINT NOT NULL,
    PRIMARY KEY (model_id, source_namespace, feature_key),
    CHECK (coordinate_basis = 'raw_spliced'),
    CHECK (cdna_start1 >= 1 AND cdna_end1 >= cdna_start1)
);

-- Xref provenance is normalized deliberately.  external_db, xref,
-- object_xref, and seq_region_synonym source row IDs are different domains;
-- no flattened alias row may substitute one for another.
CREATE TABLE duckvep_model.model_external_db (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    external_db_key             UBIGINT NOT NULL CHECK (external_db_key > 0),
    source_internal_id          VARCHAR NOT NULL,
    db_name                     VARCHAR NOT NULL,
    db_release                  VARCHAR,
    status                      VARCHAR NOT NULL CHECK (
        status IN ('KNOWNXREF', 'KNOWN', 'XREF', 'PRED', 'ORTH', 'PSEUDO')
    ),
    priority                    INTEGER NOT NULL CHECK (priority >= 0),
    db_display_name             VARCHAR,
    type                        VARCHAR NOT NULL CHECK (
        type IN (
            'ARRAY', 'ALT_TRANS', 'ALT_GENE', 'MISC', 'LIT',
            'PRIMARY_DB_SYNONYM', 'ENSEMBL'
        )
    ),
    secondary_db_name           VARCHAR,
    secondary_db_table          VARCHAR,
    description                 VARCHAR,
    PRIMARY KEY (model_id, source_namespace, external_db_key),
    UNIQUE (model_id, source_namespace, source_internal_id)
);

CREATE TABLE duckvep_model.model_xref (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    xref_key                    UBIGINT NOT NULL CHECK (xref_key > 0),
    source_internal_id          VARCHAR NOT NULL,
    external_db_key             UBIGINT NOT NULL,
    accession                   VARCHAR NOT NULL,
    accession_version           VARCHAR,
    display_label               VARCHAR NOT NULL,
    description                 VARCHAR,
    info_type                   VARCHAR NOT NULL CHECK (
        info_type IN (
            'NONE', 'PROJECTION', 'MISC', 'DEPENDENT', 'DIRECT',
            'SEQUENCE_MATCH', 'INFERRED_PAIR', 'PROBE', 'UNMAPPED',
            'COORDINATE_OVERLAP', 'CHECKSUM'
        )
    ),
    info_text                   VARCHAR NOT NULL,
    PRIMARY KEY (model_id, source_namespace, xref_key),
    UNIQUE (model_id, source_namespace, source_internal_id)
);

CREATE TABLE duckvep_model.model_object_xref (
    model_id                    VARCHAR NOT NULL,
    object_source_namespace     VARCHAR NOT NULL,
    xref_source_namespace       VARCHAR NOT NULL,
    object_xref_key             UBIGINT NOT NULL CHECK (object_xref_key > 0),
    source_internal_id          VARCHAR NOT NULL,
    object_kind                 VARCHAR NOT NULL CHECK (
        object_kind IN ('gene', 'transcript', 'translation')
    ),
    object_key                  UBIGINT NOT NULL,
    xref_key                    UBIGINT NOT NULL,
    is_display_xref             BOOLEAN NOT NULL,
    linkage_annotation          VARCHAR,
    source_analysis_internal_id VARCHAR,
    PRIMARY KEY (
        model_id, object_source_namespace,
        xref_source_namespace, object_xref_key
    ),
    UNIQUE (
        model_id, object_source_namespace,
        xref_source_namespace, source_internal_id
    )
);

CREATE TABLE duckvep_model.model_seq_region_synonym (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    synonym_key                 UBIGINT NOT NULL CHECK (synonym_key > 0),
    source_internal_id          VARCHAR NOT NULL,
    seq_region_key              UBIGINT NOT NULL,
    synonym                     VARCHAR NOT NULL,
    external_db_key             UBIGINT,
    PRIMARY KEY (
        model_id, source_namespace, synonym_key
    ),
    UNIQUE (
        model_id, source_namespace, source_internal_id
    )
);

CREATE TABLE duckvep_model.model_xref_identity (
    model_id                    VARCHAR NOT NULL,
    object_source_namespace     VARCHAR NOT NULL,
    xref_source_namespace       VARCHAR NOT NULL,
    object_xref_key             UBIGINT NOT NULL,
    identity_status             VARCHAR NOT NULL CHECK (
        identity_status IN ('unverified', 'exact', 'mismatch')
    ),
    verified_sequence_role      VARCHAR,
    verified_sequence_sha256    VARCHAR,
    PRIMARY KEY (
        model_id, object_source_namespace,
        xref_source_namespace, object_xref_key
    ),
    CHECK (
        (identity_status = 'unverified' AND verified_sequence_role IS NULL AND
         verified_sequence_sha256 IS NULL) OR
        (identity_status IN ('exact', 'mismatch') AND
         verified_sequence_role IS NOT NULL AND
         verified_sequence_sha256 IS NOT NULL)
    )
);

CREATE TABLE duckvep_model.model_sequence_blob (
    sequence_sha256             VARCHAR NOT NULL,
    alphabet                    VARCHAR NOT NULL CHECK (alphabet IN ('dna', 'peptide')),
    byte_count                  BIGINT NOT NULL CHECK (byte_count >= 0),
    sequence_bytes              BLOB NOT NULL,
    PRIMARY KEY (sequence_sha256, alphabet)
);

CREATE TABLE duckvep_model.model_sequence_state (
    model_id                    VARCHAR NOT NULL,
    source_namespace            VARCHAR NOT NULL,
    owner_kind                  VARCHAR NOT NULL CHECK (
        owner_kind IN ('transcript', 'translation')
    ),
    owner_key                   UBIGINT NOT NULL CHECK (owner_key > 0),
    role                        VARCHAR NOT NULL CHECK (
        role IN (
            'raw_spliced', 'edited_spliced', 'translatable_cds',
            'peptide_unedited', 'peptide_final'
        )
    ),
    alphabet                    VARCHAR NOT NULL CHECK (alphabet IN ('dna', 'peptide')),
    sequence_sha256             VARCHAR,
    byte_count                  BIGINT CHECK (byte_count >= 0),
    status                      VARCHAR NOT NULL CHECK (
        status IN ('present', 'absent', 'provider_error')
    ),
    PRIMARY KEY (
        model_id, source_namespace, owner_kind, owner_key, role
    ),
    CHECK (
        (owner_kind = 'transcript' AND
         role IN ('raw_spliced', 'edited_spliced') AND alphabet = 'dna') OR
        (owner_kind = 'translation' AND role = 'translatable_cds' AND alphabet = 'dna') OR
        (owner_kind = 'translation' AND
         role IN ('peptide_unedited', 'peptide_final') AND alphabet = 'peptide')
    ),
    CHECK (
        (status = 'present' AND sequence_sha256 IS NOT NULL AND byte_count IS NOT NULL) OR
        (status <> 'present' AND sequence_sha256 IS NULL AND byte_count IS NULL)
    )
);
