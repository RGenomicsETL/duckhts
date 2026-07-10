#!/usr/bin/env Rscript

stopifnot(
  requireNamespace("DBI", quietly = TRUE),
  requireNamespace("duckdb", quietly = TRUE),
  requireNamespace("jsonlite", quietly = TRUE)
)

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
test_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else
  "test/duckvep_model.R"
root <- normalizePath(file.path(dirname(test_path), ".."), mustWork = TRUE)
cli <- file.path(root, "scripts", "duckvep_model.R")
valid_sql <- file.path(root, "test", "duckvep_model_valid.sql")
invalid_rank_sql <- file.path(root, "test", "duckvep_model_invalid_rank.sql")
rscript <- file.path(R.home("bin"), "Rscript")
scratch <- tempfile("duckvep-model-test-")
dir.create(scratch)
on.exit(unlink(scratch, recursive = TRUE, force = TRUE), add = TRUE)

run_cli <- function(...) {
  output <- suppressWarnings(system2(
    rscript, c(cli, ...), stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  text <- paste(output, collapse = "\n")
  list(
    status = as.integer(status),
    text = text,
    json = jsonlite::fromJSON(text, simplifyVector = FALSE)
  )
}

execute_sql_file <- function(database, path) {
  connection <- DBI::dbConnect(duckdb::duckdb(), dbdir = database)
  on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
  sql <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  invisible(DBI::dbExecute(connection, sql))
}

execute_sql <- function(database, sql) {
  connection <- DBI::dbConnect(duckdb::duckdb(), dbdir = database)
  on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
  invisible(DBI::dbExecute(connection, sql))
}

issue_codes <- function(result) {
  if (!length(result$json$issues)) return(character())
  vapply(result$json$issues, `[[`, character(1), "code")
}

empty_database <- file.path(scratch, "empty.duckdb")
initialized <- run_cli("init", empty_database)
stopifnot(initialized$status == 0L, identical(initialized$json$ok, TRUE))

structural <- run_cli("validate", empty_database, "--profile", "structural")
stopifnot(structural$status == 0L, identical(structural$json$ok, TRUE))

composite_blob_database <- file.path(scratch, "composite-blob.duckdb")
stopifnot(run_cli("init", composite_blob_database)$status == 0L)
execute_sql(
  composite_blob_database,
  paste(
    "INSERT INTO duckvep_model.model_sequence_blob VALUES",
    "(sha256(CAST('AC' AS BLOB)),'dna',2,CAST('AC' AS BLOB)),",
    "(sha256(CAST('AC' AS BLOB)),'peptide',2,CAST('AC' AS BLOB)),",
    "(sha256(CAST(repeat('A',5000) AS BLOB)),'dna',5000,",
    "CAST(repeat('A',5000) AS BLOB))"
  )
)
composite_blob <- run_cli(
  "validate", composite_blob_database, "--profile", "structural"
)
stopifnot(composite_blob$status == 0L, identical(composite_blob$json$ok, TRUE))

empty_candidate <- run_cli(
  "validate", empty_database, "--profile", "model_candidate", "--max-samples", "0"
)
stopifnot(empty_candidate$status == 1L, identical(empty_candidate$json$ok, FALSE))

valid_database <- file.path(scratch, "valid.duckdb")
stopifnot(run_cli("init", valid_database)$status == 0L)
execute_sql_file(valid_database, valid_sql)

candidate <- run_cli("validate", valid_database, "--profile", "model_candidate")
stopifnot(
  candidate$status == 0L,
  identical(candidate$json$ok, TRUE),
  identical(candidate$json$identity_authentication, "unavailable"),
  identical(candidate$json$physical_attestation, "external_required"),
  length(candidate$json$relations) == 26L,
  all(unlist(candidate$json$relations, use.names = FALSE) > 0),
  candidate$json$relations$model_sequence_state == 5
)

expect_mutation_code <- function(name, sql, code) {
  database <- file.path(scratch, paste0(name, ".duckdb"))
  stopifnot(run_cli("init", database)$status == 0L)
  execute_sql_file(database, valid_sql)
  execute_sql(database, sql)
  result <- run_cli(
    "validate", database, "--profile", "model_candidate", "--max-samples", "1"
  )
  stopifnot(result$status == 1L, code %in% issue_codes(result))
  invisible(result)
}

expect_mutation_code(
  "bad-manifest-pin",
  "UPDATE duckvep_model.model_manifest SET compatibility_target='ensembl-vep/115.0'",
  "DVMA2002"
)
expect_mutation_code(
  "bad-core-ddl-receipt",
  paste0(
    "UPDATE duckvep_model.model_artifact SET sha256=repeat('f',64) ",
    "WHERE role='core_ddl'"
  ),
  "DVMA2045"
)
expect_mutation_code(
  "bad-selection-count",
  paste0(
    "UPDATE duckvep_model.model_selection_audit SET row_count=2 ",
    "WHERE object_kind='transcript' AND reason_code='fixture_transcript'"
  ),
  "DVMA2032"
)
expect_mutation_code(
  "wrong-rejection-ledger-class",
  paste0(
    "UPDATE duckvep_model.model_selection_audit ",
    "SET rejected_rows_artifact='duckvep-model-importer' ",
    "WHERE object_kind='transcript' AND reason_code='fixture_transcript'"
  ),
  "DVMA2044"
)
expect_mutation_code(
  "bad-raw-attribute",
  paste0(
    "UPDATE duckvep_model.model_gene_attribute SET source_row_locator='' ",
    "WHERE attribute_key=1"
  ),
  "DVMA3310"
)
expect_mutation_code(
  "bad-transcript-envelope",
  paste(
    "UPDATE duckvep_model.model_gene SET end1=18;",
    "UPDATE duckvep_model.model_transcript SET end1=18;"
  ),
  "DVMA3105"
)
expect_mutation_code(
  "bad-sequence-region-assembly",
  "UPDATE duckvep_model.model_seq_region SET assembly_accession='GCA_000001405.13'",
  "DVMA3002"
)
expect_mutation_code(
  "bad-raw-coding-bound",
  "UPDATE duckvep_model.model_translation SET raw_cdna_coding_start1=2",
  "DVMA3201"
)
expect_mutation_code(
  "bad-edited-coding-map",
  "UPDATE duckvep_model.model_translation SET edited_cdna_coding_start1=2",
  "DVMA3323"
)
expect_mutation_code(
  "missing-source-input-receipt",
  "DELETE FROM duckvep_model.model_artifact WHERE role_class='source_input'",
  "DVMA2046"
)
expect_mutation_code(
  "mismatch-without-accession-version",
  paste(
    "UPDATE duckvep_model.model_xref SET accession_version=NULL;",
    "UPDATE duckvep_model.model_xref_identity",
    "SET identity_status='mismatch',verified_sequence_sha256=repeat('f',64);"
  ),
  "DVMA3542"
)
expect_mutation_code(
  "bad-exact-sequence-identity",
  paste0(
    "UPDATE duckvep_model.model_xref_identity ",
    "SET verified_sequence_sha256=repeat('f',64)"
  ),
  "DVMA3543"
)
expect_mutation_code(
  "cross-namespace-display-key-collision",
  paste(
    "INSERT INTO duckvep_model.model_source VALUES",
    "(repeat('1',64),'overlay','Fixture','overlay','crosswalk','fixture overlay','1',",
    "NULL,NULL,NULL,NULL,1);",
    "INSERT INTO duckvep_model.model_selection_audit VALUES",
    "(repeat('1',64),'overlay','external_db','included','fixture_external_db',1,NULL),",
    "(repeat('1',64),'overlay','xref','included','fixture_xref',1,NULL);",
    "INSERT INTO duckvep_model.model_artifact VALUES",
    "(repeat('1',64),repeat('2',64),'source_input','overlay_dump','overlay',",
    "'overlay.tsv','fixture://overlay.tsv',repeat('f',64),1,1,true);",
    "INSERT INTO duckvep_model.model_external_db VALUES",
    "(repeat('1',64),'overlay',1,'1','RefSeq','1','KNOWN',1,'RefSeq','MISC',",
    "NULL,NULL,'overlay');",
    "INSERT INTO duckvep_model.model_xref VALUES",
    "(repeat('1',64),'overlay',1,'1',1,'NM_000001','1','NM_000001.1',",
    "'overlay transcript','DIRECT','fixture');",
    "UPDATE duckvep_model.model_object_xref SET xref_source_namespace='overlay';",
    "UPDATE duckvep_model.model_xref_identity SET xref_source_namespace='overlay';"
  ),
  "DVMA3521"
)
expect_mutation_code(
  "mutable-database-view",
  paste(
    "CREATE SCHEMA mutable_backing;",
    "CREATE TABLE mutable_backing.model_xref AS",
    "SELECT * FROM duckvep_model.model_xref;",
    "DROP TABLE duckvep_model.model_xref;",
    "CREATE VIEW duckvep_model.model_xref AS",
    "SELECT * FROM mutable_backing.model_xref;"
  ),
  "DVMA1006"
)
expect_mutation_code(
  "bad-edit-replay",
  paste0(
    "UPDATE duckvep_model.model_transcript_edit SET preapply_ref_seq=CAST('T' AS BLOB) ",
    "WHERE edit_key=1"
  ),
  "DVMA3404"
)
expect_mutation_code(
  "bad-display-xref",
  "UPDATE duckvep_model.model_object_xref SET is_display_xref=false",
  "DVMA3522"
)
expect_mutation_code(
  "bad-source-key-order",
  paste(
    "UPDATE duckvep_model.model_external_db SET external_db_key=2;",
    "UPDATE duckvep_model.model_xref SET external_db_key=2;"
  ),
  "DVMA3700"
)

gate <- run_cli("validate", valid_database)
stopifnot(
  gate$status == 1L,
  identical(gate$json$ok, FALSE),
  identical(issue_codes(gate), c("DVMA2901", "DVMA2902"))
)

invalid_database <- file.path(scratch, "invalid-rank.duckdb")
stopifnot(run_cli("init", invalid_database)$status == 0L)
execute_sql_file(invalid_database, valid_sql)
execute_sql_file(invalid_database, invalid_rank_sql)
invalid_rank <- run_cli(
  "validate", invalid_database, "--profile", "model_candidate", "--max-samples", "2"
)
invalid_rank_again <- run_cli(
  "validate", invalid_database, "--profile", "model_candidate", "--max-samples", "2"
)
stopifnot(
  invalid_rank$status == 1L,
  identical(issue_codes(invalid_rank), "DVMA3101"),
  identical(invalid_rank$text, invalid_rank_again$text)
)

no_transcript_edit_database <- file.path(scratch, "no-transcript-edit.duckdb")
stopifnot(run_cli("init", no_transcript_edit_database)$status == 0L)
execute_sql_file(no_transcript_edit_database, valid_sql)
execute_sql(no_transcript_edit_database, "DELETE FROM duckvep_model.model_transcript_edit")
no_transcript_edit <- run_cli(
  "validate", no_transcript_edit_database, "--profile", "model_candidate"
)
stopifnot(
  no_transcript_edit$status == 1L,
  "DVMA3614" %in% issue_codes(no_transcript_edit)
)

no_peptide_edit_database <- file.path(scratch, "no-peptide-edit.duckdb")
stopifnot(run_cli("init", no_peptide_edit_database)$status == 0L)
execute_sql_file(no_peptide_edit_database, valid_sql)
execute_sql(no_peptide_edit_database, "DELETE FROM duckvep_model.model_translation_edit")
no_peptide_edit <- run_cli(
  "validate", no_peptide_edit_database, "--profile", "model_candidate"
)
stopifnot(
  no_peptide_edit$status == 1L,
  "DVMA3617" %in% issue_codes(no_peptide_edit)
)

overlap_database <- file.path(scratch, "overlap-replay.duckdb")
stopifnot(run_cli("init", overlap_database)$status == 0L)
execute_sql_file(overlap_database, valid_sql)
execute_sql(
  overlap_database,
  paste(
    "UPDATE duckvep_model.model_transcript_edit SET apply_ordinal=3 WHERE edit_key=1;",
    "INSERT INTO duckvep_model.model_attribute_type VALUES",
    "(repeat('1',64),'ensembl_core_116',7,'7','_transl_end','Translation end',NULL);",
    "INSERT INTO duckvep_model.model_transcript_attribute VALUES",
    "(repeat('1',64),'ensembl_core_116',4,1,1,1,'fixture:transcript_attrib:4','14 15'),",
    "(repeat('1',64),'ensembl_core_116',5,1,1,1,'fixture:transcript_attrib:5','10 15 C'),",
    "(repeat('1',64),'ensembl_core_116',6,1,7,1,'fixture:transcript_attrib:6','10');",
    "INSERT INTO duckvep_model.model_transcript_edit VALUES",
    "(repeat('1',64),'ensembl_core_116',2,1,4,'_rna_edit','raw_spliced',14,15,",
    "CAST('TT' AS BLOB),CAST('TT' AS BLOB),CAST('' AS BLOB),1,'applied'),",
    "(repeat('1',64),'ensembl_core_116',3,1,5,'_rna_edit','raw_spliced',10,15,",
    "CAST('GGGTTT' AS BLOB),CAST('GGGT' AS BLOB),CAST('C' AS BLOB),2,'applied');",
    "INSERT INTO duckvep_model.model_sequence_blob VALUES",
    "(sha256(CAST('ATGGAATGAC' AS BLOB)),'dna',10,CAST('ATGGAATGAC' AS BLOB));",
    "UPDATE duckvep_model.model_translation SET edited_cdna_coding_end1=10;",
    "UPDATE duckvep_model.model_sequence_state",
    "SET sequence_sha256=sha256(CAST('ATGGAATGAC' AS BLOB)),byte_count=10",
    "WHERE owner_kind='transcript' AND role='edited_spliced';"
  )
)
overlap_replay <- run_cli(
  "validate", overlap_database, "--profile", "structural", "--max-samples", "1"
)
stopifnot(
  !"DVMA3404" %in% issue_codes(overlap_replay),
  !"DVMA3323" %in% issue_codes(overlap_replay)
)

parquet_directory <- file.path(scratch, "parquet")
dir.create(parquet_directory)
connection <- DBI::dbConnect(duckdb::duckdb(), dbdir = valid_database, read_only = TRUE)
tables <- DBI::dbGetQuery(
  connection,
  paste(
    "SELECT table_name FROM information_schema.tables",
    "WHERE table_schema='duckvep_model' ORDER BY table_name"
  )
)$table_name
for (table in tables) {
  output <- file.path(parquet_directory, paste0(table, ".parquet"))
  DBI::dbExecute(
    connection,
    sprintf(
      "COPY duckvep_model.\"%s\" TO '%s' (FORMAT PARQUET)",
      table, gsub("'", "''", output, fixed = TRUE)
    )
  )
}
invisible(DBI::dbDisconnect(connection, shutdown = TRUE))

parquet_candidate <- run_cli(
  "validate", parquet_directory, "--profile", "model_candidate"
)
stopifnot(parquet_candidate$status == 0L, identical(parquet_candidate$json$ok, TRUE))

extra_parquet <- file.path(scratch, "extra-parquet")
dir.create(extra_parquet)
invisible(file.copy(list.files(parquet_directory, full.names = TRUE), extra_parquet))
invisible(file.copy(
  file.path(parquet_directory, "model_xref.parquet"),
  file.path(extra_parquet, "unexpected_relation.parquet")
))
unexpected_relation <- run_cli(
  "validate", extra_parquet, "--profile", "model_candidate", "--max-samples", "1"
)
stopifnot(
  unexpected_relation$status == 1L,
  "DVMA1000" %in% issue_codes(unexpected_relation)
)

artifact_only <- file.path(scratch, "artifact-only")
dir.create(artifact_only)
invisible(file.copy(
  file.path(parquet_directory, "model_artifact.parquet"),
  file.path(artifact_only, "model_artifact.parquet")
))
missing_manifest <- run_cli(
  "validate", artifact_only, "--profile", "model_candidate", "--max-samples", "0"
)
stopifnot(
  missing_manifest$status == 1L,
  identical(missing_manifest$json$ok, FALSE),
  all(issue_codes(missing_manifest) == "DVMA1001")
)

invalid_parquet <- file.path(scratch, "invalid-parquet")
dir.create(invalid_parquet)
invisible(file.copy(list.files(parquet_directory, full.names = TRUE), invalid_parquet))
translation_path <- file.path(invalid_parquet, "model_translation.parquet")
replacement_path <- file.path(invalid_parquet, "model_translation.new.parquet")
connection <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
invisible(DBI::dbExecute(
  connection,
  sprintf(
    paste(
      "COPY (SELECT * REPLACE (CAST(3 AS INTEGER) AS codon_table)",
      "FROM read_parquet('%s')) TO '%s' (FORMAT PARQUET)"
    ),
    gsub("'", "''", translation_path, fixed = TRUE),
    gsub("'", "''", replacement_path, fixed = TRUE)
  )
))
invisible(DBI::dbDisconnect(connection, shutdown = TRUE))
invisible(unlink(translation_path))
invisible(file.rename(replacement_path, translation_path))
external_db_path <- file.path(invalid_parquet, "model_external_db.parquet")
external_db_replacement <- file.path(
  invalid_parquet, "model_external_db.new.parquet"
)
connection <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
invisible(DBI::dbExecute(
  connection,
  sprintf(
    paste(
      "COPY (SELECT * REPLACE ('BROKEN' AS status)",
      "FROM read_parquet('%s')) TO '%s' (FORMAT PARQUET)"
    ),
    gsub("'", "''", external_db_path, fixed = TRUE),
    gsub("'", "''", external_db_replacement, fixed = TRUE)
  )
))
invisible(DBI::dbDisconnect(connection, shutdown = TRUE))
invisible(unlink(external_db_path))
invisible(file.rename(external_db_replacement, external_db_path))
bad_enum <- run_cli(
  "validate", invalid_parquet, "--profile", "model_candidate", "--max-samples", "1"
)
stopifnot(
  bad_enum$status == 1L,
  "DVMA3200" %in% issue_codes(bad_enum),
  "DVMA3500" %in% issue_codes(bad_enum)
)

invalid_sequence <- file.path(scratch, "invalid-sequence")
dir.create(invalid_sequence)
invisible(file.copy(list.files(parquet_directory, full.names = TRUE), invalid_sequence))
blob_path <- file.path(invalid_sequence, "model_sequence_blob.parquet")
blob_replacement <- file.path(invalid_sequence, "model_sequence_blob.new.parquet")
connection <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
invisible(DBI::dbExecute(
  connection,
  sprintf(
    paste(
      "COPY (SELECT * FROM read_parquet('%s') UNION ALL",
      "SELECT sha256(from_hex('ff')),'dna',CAST(1 AS BIGINT),from_hex('ff'))",
      "TO '%s' (FORMAT PARQUET)"
    ),
    gsub("'", "''", blob_path, fixed = TRUE),
    gsub("'", "''", blob_replacement, fixed = TRUE)
  )
))
invisible(DBI::dbDisconnect(connection, shutdown = TRUE))
invisible(unlink(blob_path))
invisible(file.rename(blob_replacement, blob_path))
bad_sequence <- run_cli(
  "validate", invalid_sequence, "--profile", "model_candidate", "--max-samples", "1"
)
stopifnot(
  bad_sequence$status == 1L,
  "DVMA3600" %in% issue_codes(bad_sequence)
)

schema_output <- suppressWarnings(system2(rscript, c(cli, "schema"), stdout = TRUE))
schema_text <- paste0(paste(schema_output, collapse = "\n"), "\n")
expected_schema <- paste0(paste(
  readLines(file.path(root, "scripts", "duckvep_model_schema.sql"), warn = FALSE),
  collapse = "\n"
), "\n")
stopifnot(identical(schema_text, expected_schema))

cat("duckvep_model validator tests: OK\n")
