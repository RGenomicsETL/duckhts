#' Project a prepared DuckVEP model as read-only flat relations.
#'
#' Nested model_regions/model_transcripts are the production storage authority.
#' Flat duckvep_* witness tables remain usable without repacking. Returned SQL
#' refers to the attached catalog, so model-load connections need no TEMP views.
#' This function neither writes the model nor changes any biological field.
#'
#' @param connection DuckDB connection with the model already attached.
#' @param catalog Attached model catalog name.
#' @return Named SQL table expressions for the flat model relations.
#' @export
duckhts_bench_duckvep_relations <- function(connection, catalog) {
  tables <- DBI::dbGetQuery(connection, paste(
    "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'",
    "AND table_catalog =", DBI::dbQuoteString(connection, catalog), "ORDER BY table_name"
  ))$table_name
  qualify <- function(table) as.character(DBI::dbQuoteIdentifier(connection,
    DBI::Id(catalog = catalog, schema = "main", table = table)))
  flat <- tables[startsWith(tables, "duckvep_")]
  result <- stats::setNames(vapply(flat, qualify, ""), flat)
  nested <- c("model_regions", "model_transcripts") %in% tables
  if (!any(nested)) return(result)
  if (!all(nested)) stop("prepared model requires both model_regions and model_transcripts", call. = FALSE)
  regions <- qualify("model_regions")
  transcripts <- qualify("model_transcripts")
  projections <- c(
    duckvep_sequence_regions = paste("SELECT seq_region, sequence_length, seq_region_name AS name FROM", regions),
    duckvep_transcripts = paste("SELECT * FROM", transcripts),
    duckvep_transcript_names = paste("SELECT transcript_index, transcript_stable_id AS transcript_id FROM", transcripts),
    duckvep_exons = paste("SELECT transcript_index, exon.* FROM", transcripts,
      "CROSS JOIN UNNEST(exons) AS u(exon)"),
    duckvep_mature_mirna = paste("SELECT transcript_index, region.* FROM", transcripts,
      "CROSS JOIN UNNEST(mature_mirna_regions) AS u(region)"),
    duckvep_peptide_edits = paste("SELECT transcript_index, edit.* FROM", transcripts,
      "CROSS JOIN UNNEST(peptide_edits) AS u(edit)")
  )
  result[names(projections)] <- paste0("(", projections, ") AS duckvep_prepared")
  result
}
