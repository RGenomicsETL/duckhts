#' Estimate a Conditional Second-Genome Share from Low-Pass Counts
#'
#' Uses first-read and discordant-pair moments with fixed external calibration.
#' This estimator is separate from [rduckhts_count_error_fit()]. It is suitable
#' for low-pass, cell-free DNA or low-input DNA counts only under the supplied
#' population frequencies and calibration; the share is not an unconditional
#' biological fraction. Matched-frequency unrelated and IBD1 hypotheses can
#' have identical moments and are reported as unresolved.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param counts_table Table or view containing `sample_id`, `candidate`,
#'   `chrom`, one-based `pos`, `ref_count`, `alt_count`, `other_count`,
#'   `sample_af` and `candidate_af`. Counts and `sample_af` must agree across
#'   candidates; site keys must be unique within each sample and candidate.
#'   Frequencies refer to the counted alternate allele. Sites missing any
#'   candidate frequency are unavailable to every candidate of that sample.
#' @param seq_error Required fixed sequencing error in `[0, 0.5)`.
#' @param artefact_weight Required fixed uniform-count artefact weight in `[0, 1)`.
#' @param homozygosity_excess Required fixed excess of homozygosity in `[0, 1]`.
#' @param block_bases Bases per chromosome-local block, at least 1.
#' @param max_depth Maximum ref plus alt depth, 1 to 65535. Deeper sites are
#'   excluded; `other_count` contributes only to QC totals.
#' @param min_blocks Informative blocks required for 95 percent
#'   limits, at least 2, default 30. Matched frequencies require pair-bearing
#'   blocks; otherwise read-bearing blocks count. Lower values support small-fixture checks,
#'   not well-calibrated inference. Below this floor the share is reported but
#'   limits are NA and status is `insufficient_blocks`.
#' @param max_blocks Maximum physical blocks per sample, 1 to 4096.
#' @param max_candidates Maximum candidates per sample, 1 to 64.
#' @param max_state_bytes Maximum bytes of concurrently owned native moment
#'   buffers and fit copies. DuckDB owns SQL sorting, materialization and spill.
#' @param table_name Optional output table; otherwise return a data frame.
#' @param overwrite Whether to replace an existing output table.
#' @return One row per sample, candidate and relationship model. IBD1 rows
#'   require exact `candidate_af == sample_af` at every input site, including
#'   zero-read and excluded sites. Columns include the bounded share (0 to 0.5),
#'   block-jackknife standard error and score-inverted limits, singles share,
#'   pair residual, relationship diagnostic, QC read/pair totals, calibration,
#'   score rank and status. The discordance ratio divides discordant pairs by
#'   `SUM(choose(depth, 2) * sample_af * (1 - sample_af))`.
#'   When `table_name` is supplied, returns invisible `TRUE`.
#' @md
#' @export
rduckhts_count_mixture_lowpass <- function(
  con, counts_table, seq_error, artefact_weight, homozygosity_excess,
  block_bases = 1e7, max_depth = 8, min_blocks = 30, max_blocks = 512,
  max_candidates = 16, max_state_bytes = 2^28,
  table_name = NULL, overwrite = FALSE
) {
  .somalier_validate_output(con, table_name, overwrite)
  .somalier_validate_name(counts_table, "counts_table")
  calibration <- function(value, name, upper, closed = FALSE) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || (if (closed) value > upper else value >= upper)) {
      stop(sprintf("%s must be one finite calibration value in [0, %s%s",
                   name, upper, if (closed) "]" else ")"), call. = FALSE)
    }
    sprintf("%.17g", value)
  }
  whole <- function(value, name, upper, lower = 1) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != floor(value) || value < lower || value > upper) {
      stop(sprintf("%s must be one whole number in [%s, %s]", name, lower, upper), call. = FALSE)
    }
    sprintf("%.0f", value)
  }
  arguments <- c(
    sql_quote_string(con, counts_table),
    calibration(seq_error, "seq_error", 0.5),
    calibration(artefact_weight, "artefact_weight", 1),
    calibration(homozygosity_excess, "homozygosity_excess", 1, TRUE),
    paste0("block_bases := ", whole(block_bases, "block_bases", 2^62)),
    paste0("max_depth := ", whole(max_depth, "max_depth", 65535)),
    paste0("min_blocks := ", whole(min_blocks, "min_blocks", 2^31 - 1, 2)),
    paste0("max_blocks := ", whole(max_blocks, "max_blocks", 4096)),
    paste0("max_candidates := ", whole(max_candidates, "max_candidates", 64)),
    paste0("max_state_bytes := ", whole(max_state_bytes, "max_state_bytes", 2^62))
  )
  query <- paste0("SELECT * FROM duckhts_count_mixture_lowpass(", paste(arguments, collapse = ", "), ")")
  .somalier_publish_query(con, query, table_name, overwrite)
}
