# Small derived models for the SQL presentation differential and benchmark.
# Raw bytes belong to duckhtsbench's duckvep-projection registry workload:
# projection_reference, projection_reference_fai and projection_model_gff.
# These transformations are test fixtures, not a second general GFF importer.
duckvep_projection_cases <- c("forward", "reverse", "three_exon_phase2",
  "partial_cds_end", "noncoding_first_exon", "noncoding")

# The unchanged generator may emit the same allele/ID more than once. Give each
# physical record a comparison key; retain its original file alongside this copy.
# Only the ID column changes, never alleles, counts, ordering or annotations.
duckvep_projection_label_records <- function(input, output) {
  lines <- readLines(input)
  records <- which(!startsWith(lines, "#"))
  fields <- strsplit(lines[records], "\t", fixed = TRUE)
  for (i in seq_along(fields)) {
    stopifnot(length(fields[[i]]) >= 8L)
    fields[[i]][[3L]] <- sprintf("projection_%08d", i)
  }
  lines[records] <- vapply(fields, paste, character(1L), collapse = "\t")
  writeLines(lines, output)
  invisible(length(records))
}

duckvep_projection_fixture <- function(con, root, inputs, case, directory) {
  stopifnot(case %in% duckvep_projection_cases)
  quote <- function(x) as.character(DBI::dbQuoteString(con, x))
  execute <- function(sql) invisible(DBI::dbExecute(con, sql))
  execute(paste(readLines(file.path(root,
    "test/duckvep/conformance/minimal_model.sql")), collapse = "\n"))
  gff <- utils::read.delim(inputs[["projection_model_gff"]], header = FALSE,
    comment.char = "#", colClasses = "character", quote = "")
  if (case == "reverse") {
    gff$V7 <- "-"
    gff$V8[gff$V3 == "CDS"] <- c("1", "0")
    gff$V9[gff$V3 == "exon"] <- paste0(sub("rank=[12]", "",
      gff$V9[gff$V3 == "exon"]), c("rank=2", "rank=1"))
    execute("UPDATE duckvep_transcripts SET strand = -1,
      cds_sequence = seq_revcomp(decode(cds_sequence))::BLOB,
      pre_cds_sequence = seq_revcomp(decode(post_cds_sequence))::BLOB,
      post_cds_sequence = seq_revcomp(decode(pre_cds_sequence))::BLOB")
    execute("UPDATE duckvep_exons SET exon_cdna_start = 103 - exon_cdna_end,
      exon_cdna_end = 103 - exon_cdna_start")
  } else if (case %in% c("three_exon_phase2", "partial_cds_end", "noncoding_first_exon")) {
    starts <- c(100L, 150L, 220L)
    ends <- c(125L, 180L, 250L)
    coding_start <- if (case == "noncoding_first_exon") 158L else 108L
    coding_end <- if (case == "partial_cds_end") 239L else 240L
    # Ensembl phase 2 is GFF phase 1. The first, noncoding exon uses -1.
    phases <- if (case == "noncoding_first_exon") c(-1L, 0L, 1L) else c(2L, 2L, 0L)
    end_phases <- if (case == "noncoding_first_exon") c(-1L, 1L, 1L) else
      c(2L, 0L, if (case == "partial_cds_end") 2L else 0L)
    gff <- gff[1:2, ]
    cdna <- 1L
    execute("DELETE FROM duckvep_exons")
    reference <- DBI::dbGetQuery(con, paste0("SELECT SEQUENCE FROM read_fasta(",
      quote(inputs[["projection_reference"]]), ", region := 'chrDuck')"))$SEQUENCE
    stopifnot(length(reference) == 1L)
    cds <- pre <- post <- character()
    for (i in seq_along(starts)) {
      exon <- gff[2L, ]
      exon$V3 <- "exon"
      exon$V4 <- starts[[i]]
      exon$V5 <- ends[[i]]
      exon$V9 <- paste0("ID=exon:DUCK1-", i, ";Parent=transcript:DUCK1-201;rank=", i)
      gff <- rbind(gff, exon)
      lo <- max(starts[[i]], coding_start)
      hi <- min(ends[[i]], coding_end)
      if (lo <= hi) {
        coding <- exon
        coding$V3 <- "CDS"
        coding$V4 <- lo
        coding$V5 <- hi
        coding$V8 <- if (phases[[i]] == 0L) "0" else as.character(3L - phases[[i]])
        coding$V9 <- paste0("ID=CDS:DUCK1-", i, ";Parent=transcript:DUCK1-201")
        gff <- rbind(gff, coding)
        cds <- c(cds, substr(reference, lo, hi))
      }
      pre <- c(pre, substr(reference, starts[[i]], min(ends[[i]], coding_start - 1L)))
      post <- c(post, substr(reference, max(starts[[i]], coding_end + 1L), ends[[i]]))
      n <- ends[[i]] - starts[[i]] + 1L
      execute(sprintf("INSERT INTO duckvep_exons VALUES (0,%d,%d,%d,%d,%d,%d)",
        starts[[i]], ends[[i]], cdna, cdna+n-1L, phases[[i]], end_phases[[i]]))
      cdna <- cdna + n
    }
    padding <- if (case == "noncoding_first_exon") "" else "NN"
    execute(paste0("UPDATE duckvep_transcripts SET cds_start=", coding_start,
      ", cds_end=", coding_end, ", cds_sequence=", quote(paste0(padding, paste(cds, collapse = ""))),
      "::BLOB, pre_cds_sequence=", quote(paste(pre, collapse = "")),
      "::BLOB, post_cds_sequence=", quote(paste(post, collapse = "")), "::BLOB"))
  } else if (case == "noncoding") {
    gff <- gff[gff$V3 != "CDS", ]
    gff$V3[gff$V3 == "mRNA"] <- "ncRNA"
    gff$V9 <- sub("biotype=protein_coding", "biotype=lncRNA", gff$V9)
    execute("UPDATE duckvep_transcripts SET transcript_flags=0, cds_start=NULL, cds_end=NULL,
      codon_table=NULL, cds_sequence=NULL::BLOB,
      pre_cds_sequence=NULL::BLOB, post_cds_sequence=NULL::BLOB")
    execute("UPDATE duckvep_exons SET phase=-1, end_phase=-1")
  }
  gff_path <- file.path(directory, paste0(case, ".gff3"))
  writeLines("##gff-version 3", gff_path)
  utils::write.table(gff, gff_path, append = TRUE, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE)
  execute("CREATE OR REPLACE TABLE projection_transcripts AS SELECT t.*, exons,
    []::STRUCT(protein_position UINTEGER, alternate_amino_acid VARCHAR, edit_code VARCHAR)[] peptide_edits
    FROM duckvep_transcripts t JOIN (
      SELECT transcript_index, list(struct_pack(exon_start := exon_start, exon_end := exon_end,
        exon_cdna_start := exon_cdna_start, exon_cdna_end := exon_cdna_end, phase := phase,
        end_phase := end_phase) ORDER BY exon_cdna_start) exons
      FROM duckvep_exons GROUP BY transcript_index) e USING(transcript_index)")
  loaded <- DBI::dbGetQuery(con, "SELECT * FROM duckvep_model_load('projection',
    'SELECT seq_region FROM duckvep_sequence_regions ORDER BY seq_region',
    'SELECT * EXCLUDE(post_cds_bases) FROM duckvep_transcripts ORDER BY seq_region, transcript_start',
    'SELECT * FROM duckvep_exons ORDER BY transcript_index, exon_cdna_start')")
  stopifnot(isTRUE(loaded[[1L]]))
  gff_path
}

duckvep_projection_input <- function(con, vcf) {
  DBI::dbExecute(con, paste0("CREATE OR REPLACE TABLE projection_events AS SELECT
    (row_number() OVER (ORDER BY POS, ID)-1)::UBIGINT AS event_index,
    1::UINTEGER AS seq_region, POS::UBIGINT AS position, REF AS reference, ALT[1] AS alternate, ID,
    NULL::UBIGINT end_position, NULL::VARCHAR structural_type, NULL::VARCHAR copy_change,
    NULL::UINTEGER mate_seq_region, NULL::UBIGINT mate_position
    FROM read_bcf(", DBI::dbQuoteString(con, vcf), ") ORDER BY position, ID"))
  DBI::dbExecute(con, "CREATE OR REPLACE TABLE projection_annotations AS
    SELECT * FROM duckvep_annotate('projection_events', 'projection')")
  DBI::dbExecute(con, "CREATE OR REPLACE TABLE projection_actual AS
    SELECT p.*, e.ID, n.transcript_id FROM duckvep_transcript_projection(
      'projection_events', 'projection_annotations', 'projection_transcripts') p
    JOIN projection_events e USING(event_index)
    LEFT JOIN duckvep_transcript_names n USING(transcript_index)")
  invisible(NULL)
}
