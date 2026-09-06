#!/usr/bin/env Rscript
# Exercise the real nested-to-flat mapping without an extension or VEP process.
suppressMessages({ library(DBI); library(duckdb) })
source("r/duckhtsbench/R/duckvep_relations.R")
local({
  directory <- tempfile("duckvep-model-relations-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  path <- file.path(directory, "packed model.duckdb")
  con <- dbConnect(duckdb(), dbdir = path)
  on.exit(if (dbIsValid(con)) dbDisconnect(con, shutdown = TRUE), add = TRUE, after = FALSE)
  dbExecute(con, paste(readLines("test/duckvep/conformance/minimal_model.sql"), collapse = "\n"))
  flat_names <- c("duckvep_sequence_regions", "duckvep_transcripts", "duckvep_exons", "duckvep_transcript_names")
  expected <- setNames(lapply(flat_names, function(name) {
    dbGetQuery(con, paste("SELECT * FROM", name, "ORDER BY ALL"))
  }), flat_names)
  catalog <- dbGetQuery(con, "SELECT current_database() AS name")$name
  flat <- duckhts_bench_duckvep_relations(con, catalog)
  for (name in flat_names) {
    stopifnot(identical(dbGetQuery(con, paste("SELECT * FROM", flat[[name]], "ORDER BY ALL")),
      expected[[name]]))
  }
  dbExecute(con, "CREATE TABLE model_regions AS
    SELECT seq_region, name AS seq_region_name, 1000::UBIGINT AS sequence_length
    FROM duckvep_sequence_regions")
  error <- tryCatch({ duckhts_bench_duckvep_relations(con, catalog); "" }, error = conditionMessage)
  stopifnot(grepl("both model_regions and model_transcripts", error, fixed = TRUE))
  dbExecute(con, "CREATE TABLE model_transcripts AS
    SELECT t.*, n.transcript_id AS transcript_stable_id,
      (SELECT list(struct_pack(exon_start := e.exon_start, exon_end := e.exon_end,
          exon_cdna_start := e.exon_cdna_start, exon_cdna_end := e.exon_cdna_end,
          phase := e.phase, end_phase := e.end_phase) ORDER BY e.exon_cdna_start)
       FROM duckvep_exons e WHERE e.transcript_index = t.transcript_index) AS exons,
      [struct_pack(mature_mirna_start := 3::UBIGINT, mature_mirna_end := 8::UBIGINT)] AS mature_mirna_regions,
      [struct_pack(protein_position := 2::UINTEGER, alternate_amino_acid := 'U'::VARCHAR)] AS peptide_edits
    FROM duckvep_transcripts t JOIN duckvep_transcript_names n USING(transcript_index)")
  # Stale persisted projections must never override canonical nested model data.
  dbExecute(con, "UPDATE duckvep_transcripts SET cds_start = 999")
  dbExecute(con, "CHECKPOINT")
  dbDisconnect(con, shutdown = TRUE)
  before <- unname(tools::md5sum(path))

  driver <- duckdb()
  first <- dbConnect(driver)
  on.exit(dbDisconnect(first, shutdown = TRUE), add = TRUE, after = FALSE)
  attached <- 'quoted "model'
  dbExecute(first, paste("ATTACH", dbQuoteString(first, path), "AS",
    dbQuoteIdentifier(first, attached), "(READ_ONLY)"))
  relations <- duckhts_bench_duckvep_relations(first, attached)
  second <- dbConnect(driver)
  on.exit(dbDisconnect(second), add = TRUE, after = FALSE)
  for (name in flat_names) {
    columns <- paste(dbQuoteIdentifier(second, names(expected[[name]])), collapse = ", ")
    actual <- dbGetQuery(second, paste("SELECT", columns, "FROM", relations[[name]], "ORDER BY ALL"))
    stopifnot(identical(actual, expected[[name]]))
  }
  mirna <- dbGetQuery(second, paste("SELECT * FROM", relations[["duckvep_mature_mirna"]]))
  edits <- dbGetQuery(second, paste("SELECT * FROM", relations[["duckvep_peptide_edits"]]))
  stopifnot(nrow(mirna) == 1L, mirna$transcript_index == 0, mirna$mature_mirna_start == 3,
    mirna$mature_mirna_end == 8, nrow(edits) == 1L, edits$transcript_index == 0,
    edits$protein_position == 2, edits$alternate_amino_acid == "U",
    identical(unname(tools::md5sum(path)), before))
  cat("DuckVEP model relations: flat witnesses, nested projections, quoted catalogs, separate connections and unchanged model bytes: OK\n")
})
