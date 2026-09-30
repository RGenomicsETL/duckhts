#' Find Runs of Homozygosity
#'
#' Detect runs of homozygosity in a VCF or BCF with the two-state hidden Markov
#' model of `bcftools roh`: an autozygous state and a Hardy-Weinberg state,
#' decoded with the Viterbi algorithm. Segments, marker counts and start and end
#' positions equal `bcftools roh`; quality is the mean forward-backward phred
#' score, equal to the one decimal `bcftools` prints. The model is a port of
#' `vcfroh.c` and `HMM.c` (MIT, Genome Research Ltd), run by the native
#' `duckhts_roh_segments()` kernel over one list per sample and chromosome.
#'
#' Allele frequencies come from an INFO tag (`af_tag`, like `--AF-tag`, which
#' must be declared `Type=Float,Number=A`) or from a relation (`af_table`, like
#' `--AF-file`) with columns `chrom`, `pos`, `ref`, `alt` and `af`, matched to
#' each record by chromosome, position, REF and its ALT alleles joined by commas.
#' Sites without a usable frequency (absent, missing or exactly 0) are skipped,
#' as are records with more than one ALT or none.
#'
#' Genotype evidence is FORMAT/PL unless `gt_error` is given, which uses the
#' diploid GT calls with that phred error (`-G`) and needs no PL. With no
#' recombination map, transitions are the per-base-pair probabilities
#' `hw_to_az` and `az_to_hw` compounded over the physical distance between
#' sites, exactly as `bcftools roh` does without `-m` or `-M`. `rec_rate` is a
#' constant rate per base pair (`-M`). `genetic_map` names a relation with
#' columns `chrom`, `pos` and `cm` (cumulative centimorgans, the IMPUTE2
#' format), interpolated as `-m` does; chromosomes it does not cover are
#' skipped.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param path One VCF/BCF path or URI.
#' @param af_tag INFO tag holding the ALT allele frequency, for example `"AF"`.
#' @param af_table Name of a table or view of allele frequencies, instead of
#'   `af_tag`. Exactly one of the two is required.
#' @param genetic_map Optional name of a genetic-map table or view.
#' @param hw_to_az,az_to_hw Per-base-pair transition probabilities
#'   Hardy-Weinberg to autozygous (`-a`) and back (`-H`), each in `[0, 1]`.
#' @param gt_error Optional phred error for GT-only emissions (`-G`), at least 0.
#' @param rec_rate Optional constant recombination rate per base pair (`-M`).
#' @param samples Optional HTSlib sample selector: comma-separated inclusion,
#'   leading `^` exclusion, `"-"` for all, or `""` for none.
#' @param table_name Optional output table. `NULL` returns a data frame.
#' @param overwrite Whether an existing output table may be replaced.
#' @return A data frame (or invisible `TRUE` when `table_name` is given) with one
#'   row per run: `sample`, `chrom`, `start` and `end` (one-based, inclusive,
#'   the first and last marker), `length` (`end - start + 1`), `n_markers` and
#'   `quality`. Rows are unordered.
#' @examples
#' con <- rduckhts_connect()
#' path <- system.file("extdata", "roh_fixture.vcf.gz", package = "Rduckhts")
#' roh <- rduckhts_roh(con, path, af_tag = "AF")
#' roh[order(roh$sample, roh$chrom, roh$start), ]
#' DBI::dbDisconnect(con, shutdown = TRUE)
#' @export
rduckhts_roh <- function(
  con, path, af_tag = NULL, af_table = NULL, genetic_map = NULL,
  hw_to_az = 6.7e-8, az_to_hw = 5e-9, gt_error = NULL, rec_rate = NULL,
  samples = NULL, table_name = NULL, overwrite = FALSE
) {
  .somalier_validate_output(con, table_name, overwrite)
  .somalier_scalar_text(path, "path")
  if (is.null(af_tag) == is.null(af_table)) {
    stop("exactly one of af_tag and af_table is required", call. = FALSE)
  }
  if (!is.null(af_tag)) .somalier_scalar_text(af_tag, "af_tag")
  if (!is.null(af_table)) .somalier_validate_name(af_table, "af_table")
  if (!is.null(genetic_map)) .somalier_validate_name(genetic_map, "genetic_map")
  if (!is.null(samples)) .somalier_scalar_text(samples, "samples", allow_empty = TRUE)
  hw_to_az <- .roh_number(hw_to_az, "hw_to_az", maximum = 1)
  az_to_hw <- .roh_number(az_to_hw, "az_to_hw", maximum = 1)
  if (!is.null(gt_error)) gt_error <- .roh_number(gt_error, "gt_error")
  if (!is.null(rec_rate)) rec_rate <- .roh_number(rec_rate, "rec_rate")

  arguments <- c(
    sql_quote_string(con, path),
    sql_quote_string(con, if (is.null(af_tag)) af_table else af_tag),
    if (!is.null(genetic_map)) sql_quote_string(con, genetic_map),
    paste0("hw_to_az := ", hw_to_az),
    paste0("az_to_hw := ", az_to_hw),
    if (!is.null(gt_error)) paste0("gt_error := ", gt_error),
    if (!is.null(rec_rate)) paste0("rec_rate := ", rec_rate),
    if (!is.null(samples)) paste0("samples := ", sql_quote_string(con, samples))
  )
  macro <- if (is.null(af_tag)) "duckhts_roh_af_table" else "duckhts_roh"
  query <- paste0("SELECT * FROM ", macro, "(", paste(arguments, collapse = ", "), ")")
  .somalier_publish_query(con, query, table_name, overwrite)
}

# A finite number in [0, maximum], formatted with full precision for SQL.
.roh_number <- function(value, name, maximum = Inf) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0 || value > maximum) {
    stop(name, " must be one finite number in [0, ",
         if (is.finite(maximum)) format(maximum) else "Inf", "]", call. = FALSE)
  }
  sprintf("%.17g", value)
}
