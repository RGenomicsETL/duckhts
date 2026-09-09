#' Replay Phased Transcript Haplotypes
#'
#' Consume a SELECT query of flat event-by-transcript-by-sample calls through
#' the bundled native replay stream. DuckDB derives complete phase-set domains
#' and sorts the calls. Each output row is one occupied shared path, with CDS,
#' protein, carrier keys and all contributing events. Incomplete calls and
#' projection/edit failures retain provenance without inventing sequence.
#'
#' `coding_blocks` groups physical edits that share an alternate codon or
#' displace then restore the reading frame. Each block gives one-based reference
#' `cds_start`, transcript-oriented `reference`/`alternate` spans, zero-based
#' `alt_start0` in the rebuilt CDS, `length_change`, `sequence_flags` and
#' `event_indices`: one source event ID per physical edit in ascending CDS order.
#' An ID can repeat for several differing islands, within or across blocks;
#' the list length is the physical edit count.
#' Join IDs to `contributors` for raw alleles and input provenance.
#' Spans include retained bases between edits; they are not aligned
#' differences or HGVS normalization. An insertion has empty reference; a deletion
#' has empty alternate. Unknown sequences have NULL blocks, not an empty known
#' result. Blocks reuse `max_leaf_edits` capacity and the per-call workspace limit.
#'
#' `cds_differences` instead contains aligned differing runs: zero-based
#' `ref_start0`, `alt_start0`, and `alignment_start0`, with borrowed sequence spans
#' materialized as `reference` and `alternate`. An empty span denotes a gap.
#' The reference uses replay's uppercase DNA spelling; model bytes stay immutable.
#' Runs join only when adjacent columns have the same gap/non-gap type on both
#' sides. Indel-bearing paths use the VEP-116 pure-Perl global-alignment score and
#' tie order; substitution-only paths compare corresponding positions. Repeated
#' sequence can therefore place a difference away from its contributing event.
#' Differences are not HGVS normalization or event-provenance reassignment.
#' Unknown sequences have NULL differences. `max_alignment_cells` bounds the
#' exact traceback band and `max_leaf_differences` bounds output runs; exhaustion
#' is an error, not approximate alignment or discarded differences.
#'
#' `protein_differences` has the same span fields and alignment mode, with
#' positions in amino acids. Its reference follows Ensembl-116 start-methionine,
#' terminal-stop and curated single-residue peptide-edit rules. Haplosaurus then
#' appends `*` only for an exact uppercase TAA/TAG/TGA raw-CDS suffix, including
#' nonstandard-table and partial-CDS cases. The alternate is the displayed
#' first-stop prefix without reference peptide edits. Unknown paths and reference
#' CDS shorter than one codon have NULL protein differences; identical known
#' proteins have an empty list. Both difference axes reuse the same native
#' scratch and are separately subject to the alignment-cell and run limits.
#'
#' `hgvs = TRUE` requests a VEP-116-derived protein HGVS suffix in `hgvsp`: equality,
#' one operation, or a cis allele such as `p.[(Gly2del;Ala4CysfsTer2)]`.
#' Protein operations use the completed path and may combine several physical
#' edits; source contributors and coding blocks remain unchanged. `hgvsp_status`
#' is `ok` only after the complete operation set has been rendered. Other statuses
#' retain NULL text without discarding the sequence or provenance. Conditional or
#' unavailable sequence, ordered overlapping replacements, missing peptide data,
#' unsupported coding contexts, and unrepresentable protein ends remain explicit.
#' Phase policy controls allele assignment, not HGVS nomenclature. Protein HGVS
#' retains VEP's local-peptide and unknown-residue presentation. An `ok` status
#' reports a supported computation, not independent HGVS-rule certification.
#' A path with one original ALT source uses independent-event VEP-116 HGVS,
#' including genomic shifting and absent results. Multiple differing islands
#' within one MNV retain that single source identity. Placement needing an
#' unavailable genomic FASTA returns `missing_reference`.
#' Prepared references retain their own residues and length without changing raw
#' CDS/frame facts or inventing source edits. Loss of only a reference stop marker
#' supplies no alternate extension, and insertions require reference flanks.
#' `max_hgvs_operations` bounds the working operation stack and final operations;
#' `max_hgvs_bytes` bounds text bytes excluding NUL; `max_hgvs_reference_bytes`
#' bounds the query-local FASTA result buffer, including NUL and line-ending
#' scratch. Sequence/edit scratch derives from `max_sequence_bases` and
#' `max_leaf_edits`, and allele scratch from `max_allele_bytes` and the literal
#' allele width. All DuckVEP-owned buffers count toward `workspace_limit`;
#' HTSlib handle/transport storage is separate. Exhaustion errors instead of
#' truncating text. Disabled HGVS allocates no HGVS buffers or reference handle
#' and returns `not_requested`. Identifiers are joined
#' through the model transcript ordinal; the suffix contains no accession.
#'
#' `stop_in_displaced_frame` reports whether any of the first translated stop
#' codon's three bases overlaps a frame-displaced span of the rebuilt CDS.
#' Displacement starts at a frame-changing edit and ends after the alternate
#' bases of the restoring edit, or continues downstream if unrestored.
#' It is FALSE for no stop or a stop after frame restoration, and NA when
#' sequence is unavailable. It is a sequence fact, not a combined SO consequence
#' or a claim that a restored DNA frame rescues the protein.
#'
#' Whole-haplotype SO, DNA HGVS, complete protein HGVS and structural-event
#' composition remain unfinished. Input must contain one row per
#' `event_index`, `transcript_index`, `sample_index`, with columns `seq_region`,
#' `position`, `reference`, `alternate`, `alt_index`, `alleles`, `phase_before`
#' and nullable `phase_set`. Event indices identify individual ALT events; retain
#' their source-record mapping. Transcript ordinals belong to the named model.
#' Candidate selection is explicit in this input relation.
#'
#' With `input_mode = "source_records"`, the required columns are `event_index`,
#' `seq_region`, `position`, `reference`, `alternates` (a character list),
#' `transcript_index`, `sample_index`, and `gt` (original VCF text).
#' Here `event_index` identifies the whole source record. This mode requires
#' `vep116_compat`: the pinned file profile consumes two slots and ignores PS.
#' Missing calls and undefined slots can yield `conditional` sequence with
#' evidence bit 8; this is not known phase or biological rescue. Contributor
#' `alt_index` is 0 for REF, a positive source ALT ordinal, or NA for an undefined
#' slot's full-REF deletion. Source ALT strings must be nonempty and nonmissing.
#' Projection failures still withhold sequence. Raw spelling must be retained at
#' ingestion; it cannot be reconstructed losslessly from decoded GT arrays.
#'
#' Preparation reads committed objects on the registry's retained connection;
#' caller-local temporary objects and uncommitted changes are not visible. One
#' preparation may run per registry at a time. Nested or concurrent preparation
#' returns a busy error; completed scans use independent native workspaces.
#'
#' @inheritParams rduckhts_geno
#' @param calls_query One nonempty SELECT query supplying the call relation.
#' @param model_name Name of an already loaded DuckVEP model.
#' @param phase_policy Strict GT/PS interpretation or VEP-116 called-slot order.
#'   Decoded missing calls remain incomplete; source-record input uses the
#'   pinned raw parser and explicitly conditional missing-slot interpretation.
#' @param input_mode `alt_events` for decoded per-ALT calls, or `source_records`
#'   for raw GT and complete source ALT lists under `vep116_compat`.
#' @param hgvs Whether to request bounded protein HGVS for supported completed paths.
#' @param ... Named positive integer workspace capacities accepted by
#'   `duckvep_haplotypes`, such as `max_active_events`, `max_active_carriers`,
#'   `max_sequence_bases`, `max_ploidy`, `max_phase_sets`, and `workspace_limit`.
#' @return A data frame, or invisible `TRUE` when creating `table_name`.
#' @export
rduckhts_haplotypes <- function(con, calls_query, model_name,
                               phase_policy = c("strict", "vep116_compat"),
                               ..., input_mode = c("alt_events", "source_records"),
                               hgvs = FALSE, table_name = NULL, overwrite = FALSE) {
  for (name in c("calls_query", "model_name")) {
    value <- get(name)
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop(name, " must be one nonempty string", call. = FALSE)
    }
  }
  phase_policy <- match.arg(phase_policy)
  input_mode <- match.arg(input_mode)
  if (input_mode == "source_records" && phase_policy != "vep116_compat") {
    stop("source_records requires phase_policy='vep116_compat'", call. = FALSE)
  }
  if (!is.logical(hgvs) || length(hgvs) != 1L || is.na(hgvs)) {
    stop("hgvs must be TRUE or FALSE", call. = FALSE)
  }
  limits <- list(...)
  if (length(limits) && (is.null(names(limits)) || anyDuplicated(names(limits)) ||
      any(!grepl("^[a-z][a-z_]*$", names(limits))))) {
    stop("capacities must have unique SQL parameter names", call. = FALSE)
  }
  params <- list(phase_policy = sql_quote_string(con, phase_policy),
                 input_mode = sql_quote_string(con, input_mode), hgvs = if (hgvs) "TRUE" else "FALSE")
  for (name in names(limits)) {
    value <- limits[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0 || value > 2^53 || value != floor(value)) {
      stop(name, " must be one positive, exactly representable integer", call. = FALSE)
    }
    params[[name]] <- format(value, scientific = FALSE, trim = TRUE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE", call. = FALSE)
  }
  query <- paste0("SELECT * FROM duckvep_haplotypes(", sql_quote_string(con, calls_query),
                  ",", sql_quote_string(con, model_name), build_param_str(params), ")")
  if (is.null(table_name)) return(DBI::dbGetQuery(con, query))
  prefix <- if (overwrite) "CREATE OR REPLACE TABLE " else "CREATE TABLE "
  DBI::dbExecute(con, paste0(prefix, sql_quote_identifier(con, table_name), " AS ", query))
  invisible(TRUE)
}
