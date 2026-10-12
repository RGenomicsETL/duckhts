.depth_number <- function(value, name, check, requirement) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || !check(value)) {
    stop(name, " must be ", requirement, call. = FALSE)
  }
  sprintf("%.17g", value)
}

#' Correct Bin Read Counts for GC Content
#'
#' Gives each bin of each sample an expected read count from the median read
#' rate of the sample's autosomal training bins in the same GC stratum, so that
#' `ratio = reads / expected` is about 1 for two copies, 0.5 for one and 1.5
#' for three. The output feeds [rduckhts_depth_region_fraction()] and
#' [rduckhts_depth_sex_fraction()].
#'
#' The GC stratum of a bin is `floor(gc / gc_step)`. Training bins are
#' autosomal bins that pass the bin rules, hold at least one read and have a
#' read rate inside the `outlier_quantile` and `1 - outlier_quantile`
#' quantiles of the sample. A bin is kept unless it lacks a count or a GC
#' value, has a mappability below `min_mappability` (when the relation has a
#' `mappability` column), falls in a stratum with fewer than
#' `min_stratum_bins` training bins, or has a ratio above `max_ratio`; `reason`
#' names the rule. Only `max_ratio` looks at a bin's count, so a bin with no
#' reads stays kept. Autosomal means a contig whose `duckhts_contig_key()` is
#' not X, Y or MT. Restrict contigs and remove problem bins in SQL first.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param bins_table Name of a table or view with `sample_id`, `chrom`,
#'   `start`, `end` (zero-based, half-open), `reads` and `gc` (the G+C share of
#'   the A, C, G and T reference bases, from 0 to 1), and optionally
#'   `mappability` (0 to 1).
#' @param min_mappability Bins with a lower mappability are not kept.
#' @param gc_step Width of a GC stratum, in (0, 1].
#' @param min_stratum_bins Training bins a stratum needs, a whole number.
#' @param outlier_quantile Training bins outside this quantile of the read
#'   rate, and its complement, set no stratum median; in \[0, 0.5).
#' @param max_ratio Bins with a higher ratio are not kept; above 1.
#' @param table_name Optional table to create instead of returning a data
#'   frame.
#' @param overwrite Whether `table_name` may replace an existing table.
#' @return A data frame (or invisible `TRUE` when `table_name` is given) with
#'   one row per input row: `sample_id`, `chrom`, `start`, `end`, `reads`,
#'   `gc`, `mappability`, `gc_stratum`, `stratum_bins`, `expected`, `ratio`,
#'   `log2_ratio`, `reason`, `kept` and `method`.
#' @seealso [rduckhts_bam_bin_counts()], [rduckhts_fasta_nuc()]
#' @examples
#' con <- rduckhts_connect()
#' DBI::dbExecute(con, paste(
#'   "CREATE TEMP TABLE bins AS SELECT 'S1' AS sample_id, CAST(c AS VARCHAR) AS chrom,",
#'   "b * 1000000 AS start, (b + 1) * 1000000 AS \"end\",",
#'   "0.35 + 0.25 * ((b * 37 + c * 11) % 100) / 100.0 AS gc,",
#'   "CAST(600 + 4 * ((b * 13 + c * 7) % 11) AS BIGINT) AS reads",
#'   "FROM range(1, 23) t(c), range(0, 100) r(b)"
#' ))
#' head(rduckhts_depth_bins(con, "bins")[, c("chrom", "start", "reads", "ratio", "kept")])
#' DBI::dbDisconnect(con, shutdown = TRUE)
#' @export
rduckhts_depth_bins <- function(con, bins_table, min_mappability = 0.9, gc_step = 0.01,
                                min_stratum_bins = 10, outlier_quantile = 0.01,
                                max_ratio = 4, table_name = NULL, overwrite = FALSE) {
  .somalier_validate_output(con, table_name, overwrite)
  .somalier_validate_name(bins_table, "bins_table")
  arguments <- c(
    sql_quote_string(con, bins_table),
    paste0("min_mappability := ", .depth_number(
      min_mappability, "min_mappability", function(x) x >= 0 && x <= 1, "a number in [0, 1]")),
    paste0("gc_step := ", .depth_number(
      gc_step, "gc_step", function(x) x > 0 && x <= 1, "a number in (0, 1]")),
    paste0("min_stratum_bins := ", .depth_number(
      min_stratum_bins, "min_stratum_bins", function(x) x == floor(x) && x >= 1 && x <= 1e9,
      "a whole number in [1, 1e9]")),
    paste0("outlier_quantile := ", .depth_number(
      outlier_quantile, "outlier_quantile", function(x) x >= 0 && x < 0.5, "a number in [0, 0.5)")),
    paste0("max_ratio := ", .depth_number(
      max_ratio, "max_ratio", function(x) x > 1, "a finite number above 1"))
  )
  query <- paste0("SELECT * FROM duckhts_depth_bins(", paste(arguments, collapse = ", "), ")")
  .somalier_publish_query(con, query, table_name, overwrite)
}

#' Estimate the Share of Cells with a One-Copy Change of a Region
#'
#' For each sample and named region, compares the depth of the region's kept
#' bins with the depth of the sample's other kept autosomal bins:
#' `ratio = (C_r / E_r) / (C_b / E_b)` with `C` and `E` the sums of reads and
#' expected counts, and `shift = 2 * (ratio - 1)`, the share of cells with a
#' one-copy gain (positive) or loss (negative) relative to two copies when the
#' baseline is unchanged. A one-copy change in every cell gives 1 or -1.
#'
#' Without `reference_table`, the interval treats bin counts as independent
#' and quasi-Poisson, with the dispersion of the baseline bins, and is
#' `shift +- 1.959964 shift_se`. With `reference_table`, corrected bins of
#' reference samples, the reported `shift` is the sample's shift minus the
#' mean shift of the references (leaving out the tested sample if it is a
#' reference), and the interval uses Student's t with `reference_samples - 1`
#' degrees of freedom. Other changes in the sample move the baseline.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param bins_table Name of a table or view of corrected bins with
#'   `sample_id`, `chrom`, `start`, `end`, `reads`, `expected` and `kept`, such
#'   as the output of [rduckhts_depth_bins()].
#' @param regions_table Name of a table or view with `region`, `chrom`, `start`
#'   and `end` (zero-based, half-open); rows with one name form one region, and
#'   a bin belongs to it when its midpoint lies in one of its intervals.
#' @param reference_table Optional name of a table or view of corrected bins
#'   of reference samples, in the layout of `bins_table`.
#' @inheritParams rduckhts_depth_bins
#' @return A data frame (or invisible `TRUE` when `table_name` is given) with
#'   one row per sample and region: `sample_id`, `region`, `bins`,
#'   `baseline_bins`, `ratio`, `shift`, `shift_se`, `shift_low`, `shift_high`,
#'   `dispersion`, `reference_shift`, `reference_samples`, `status` (`ok`,
#'   `no_bins`, `no_baseline` or `few_reference_samples`) and `method`.
#' @examples
#' con <- rduckhts_connect()
#' DBI::dbExecute(con, paste(
#'   "CREATE TEMP TABLE bins AS SELECT 'S1' AS sample_id, CAST(c AS VARCHAR) AS chrom,",
#'   "b * 1000000 AS start, (b + 1) * 1000000 AS \"end\",",
#'   "0.35 + 0.25 * ((b * 37 + c * 11) % 100) / 100.0 AS gc,",
#'   "CAST((600 + 4 * ((b * 13 + c * 7) % 11)) * CASE WHEN c = 8 THEN 1.1 ELSE 1 END AS BIGINT) AS reads",
#'   "FROM range(1, 23) t(c), range(0, 100) r(b)"
#' ))
#' rduckhts_depth_bins(con, "bins", table_name = "corrected")
#' DBI::dbExecute(con, paste(
#'   "CREATE TEMP TABLE chromosomes AS SELECT CAST(c AS VARCHAR) AS region,",
#'   "CAST(c AS VARCHAR) AS chrom, 0 AS start, 100000000 AS \"end\" FROM range(1, 23) t(c)"
#' ))
#' fractions <- rduckhts_depth_region_fraction(con, "corrected", "chromosomes")
#' fractions[fractions$region == "8", c("region", "shift", "shift_low", "shift_high")]
#' DBI::dbDisconnect(con, shutdown = TRUE)
#' @export
rduckhts_depth_region_fraction <- function(con, bins_table, regions_table, reference_table = NULL,
                                           table_name = NULL, overwrite = FALSE) {
  .somalier_validate_output(con, table_name, overwrite)
  .somalier_validate_name(bins_table, "bins_table")
  .somalier_validate_name(regions_table, "regions_table")
  arguments <- c(sql_quote_string(con, bins_table), sql_quote_string(con, regions_table))
  if (!is.null(reference_table)) {
    .somalier_validate_name(reference_table, "reference_table")
    arguments <- c(arguments, sql_quote_string(con, reference_table))
  }
  query <- paste0("SELECT * FROM duckhts_depth_region_fraction(",
                  paste(arguments, collapse = ", "), ")")
  .somalier_publish_query(con, query, table_name, overwrite)
}

#' Estimate the Share of a Second Genome of the Other Sex from X and Y Depth
#'
#' Treats each sample as a major genome and a second genome of the other sex
#' with a DNA share `f` below 0.5. `x_copies` and `y_copies` are twice the
#' depth ratio of the kept X or Y bins to the kept autosomal bins. For an XX
#' major genome `x_copies = 2 - f` and `y_copies = f`; for an XY one
#' `x_copies = 1 + f` and `y_copies = 1 - f`. The major genome is XX when
#' `x_copies >= 1.5`. `fraction_x` and `fraction_y` are the two estimates,
#' unclipped; `fraction` combines them by inverse variance and `agreement_z`
#' compares them, since an X or Y copy change of the major genome moves them
#' differently.
#'
#' Without `reference_table`, one X copy is 1, two are 2, one Y copy is 1 and
#' a genome without Y has no Y depth. With `reference_table`, its XY samples
#' set the depth of one X and one Y copy, and its XX samples, when present,
#' set two X copies and the Y background; the tested sample is left out of
#' its own reference. Remove the pseudo-autosomal regions of X, and Y bins that
#' are not single-copy, in SQL first.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param bins_table Name of a table or view of corrected bins, as for
#'   [rduckhts_depth_region_fraction()].
#' @param reference_table Optional name of a table or view of corrected bins
#'   of reference samples with a `sex` column, `XX` or `XY`.
#' @inheritParams rduckhts_depth_bins
#' @return A data frame (or invisible `TRUE` when `table_name` is given) with
#'   one row per sample: `sample_id`, `autosomal_bins`, `x_bins`, `y_bins`,
#'   `dispersion`, `x_copies`, `x_copies_se`, `y_copies`, `y_copies_se`,
#'   `major_sex`, `fraction_x`, `fraction_x_se`, `fraction_y`,
#'   `fraction_y_se`, `fraction`, `fraction_se`, `fraction_low`,
#'   `fraction_high`, `agreement_z`, `reference_samples`, `status` (`ok`,
#'   `no_autosomal_bins`, `no_xy_reference`, `bad_reference`, `no_sex_bins`,
#'   `no_x_bins`, `no_y_bins` or `x_y_disagree`) and `method`.
#' @examples
#' con <- rduckhts_connect()
#' DBI::dbExecute(con, paste(
#'   "CREATE TEMP TABLE bins AS SELECT 'S1' AS sample_id, chrom,",
#'   "b * 1000000 AS start, (b + 1) * 1000000 AS \"end\",",
#'   "0.35 + 0.25 * ((b * 37 + length(chrom) * 11) % 100) / 100.0 AS gc,",
#'   "CAST((600 + 4 * ((b * 13) % 11)) * CASE chrom WHEN 'X' THEN 0.95 WHEN 'Y' THEN 0.05",
#'   "ELSE 1 END AS BIGINT) AS reads",
#'   "FROM (SELECT CAST(c AS VARCHAR) AS chrom FROM range(1, 23) t(c)",
#'   "UNION ALL SELECT 'X' UNION ALL SELECT 'Y') g, range(0, 60) r(b)"
#' ))
#' rduckhts_depth_bins(con, "bins", table_name = "corrected")
#' rduckhts_depth_sex_fraction(con, "corrected")[, c("major_sex", "fraction", "status")]
#' DBI::dbDisconnect(con, shutdown = TRUE)
#' @export
rduckhts_depth_sex_fraction <- function(con, bins_table, reference_table = NULL,
                                        table_name = NULL, overwrite = FALSE) {
  .somalier_validate_output(con, table_name, overwrite)
  .somalier_validate_name(bins_table, "bins_table")
  arguments <- sql_quote_string(con, bins_table)
  if (!is.null(reference_table)) {
    .somalier_validate_name(reference_table, "reference_table")
    arguments <- c(arguments, sql_quote_string(con, reference_table))
  }
  query <- paste0("SELECT * FROM duckhts_depth_sex_fraction(",
                  paste(arguments, collapse = ", "), ")")
  .somalier_publish_query(con, query, table_name, overwrite)
}
