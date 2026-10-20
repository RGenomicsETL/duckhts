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
#' Normalization uses eligible autosomal bins with positive reads; empty bins
#' retain zero ratios when their bin and stratum rules pass.
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
#' @param gc_step Width of a GC stratum, in \[1/2147483647, 1\], so the
#'   stratum number fits the integer output column.
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
      gc_step, "gc_step", function(x) x >= 1 / .Machine$integer.max && x <= 1,
      "a number in [1/2147483647, 1]")),
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
#'   "CAST((600 + 4 * ((b * 13 + c * 7) % 11)) *",
#'   "CASE WHEN c = 8 THEN 1.1 ELSE 1 END AS BIGINT) AS reads",
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
  macro <- if (is.null(reference_table)) "duckhts_depth_region_fraction" else
    "duckhts_depth_region_fraction_ref"
  query <- paste0("SELECT * FROM ", macro, "(", paste(arguments, collapse = ", "), ")")
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
  macro <- if (is.null(reference_table)) "duckhts_depth_sex_fraction" else
    "duckhts_depth_sex_fraction_ref"
  query <- paste0("SELECT * FROM ", macro, "(", paste(arguments, collapse = ", "), ")")
  .somalier_publish_query(con, query, table_name, overwrite)
}

.depth_number_list <- function(values, name, check, requirement) {
  if (!is.numeric(values) || length(values) < 1L || length(values) > 100L ||
      !all(is.finite(values)) || !all(check(values))) {
    stop(name, " must be 1 to 100 numbers, each ", requirement, call. = FALSE)
  }
  paste0("[", paste(sprintf("%.17g", values), collapse = ", "), "]")
}

.depth_state_query <- function(con, macro, bins_table, normal_fractions, ploidies, max_cn, txn_e,
                               txn_strength, max_iter, min_segment_bins, alt_fraction_threshold,
                               close_margin, max_bins) {
  .somalier_validate_name(bins_table, "bins_table")
  arguments <- c(
    sql_quote_string(con, bins_table),
    paste0("normal_fractions := ", .depth_number_list(
      normal_fractions, "normal_fractions", function(x) x > 0 & x < 1, "in (0, 1)")),
    paste0("ploidies := ", .depth_number_list(
      ploidies, "ploidies", function(x) x >= 1 & x <= 10, "in [1, 10]")),
    paste0("max_cn := ", .depth_number(
      max_cn, "max_cn", function(x) x == floor(x) && x >= 2 && x <= 20, "a whole number in [2, 20]")),
    paste0("txn_e := ", .depth_number(
      txn_e, "txn_e", function(x) x > 0 && x < 1, "a number in (0, 1)")),
    paste0("txn_strength := ", .depth_number(
      txn_strength, "txn_strength", function(x) x > 0, "a positive number")),
    paste0("max_iter := ", .depth_number(
      max_iter, "max_iter", function(x) x == floor(x) && x >= 2 && x <= 1000,
      "a whole number in [2, 1000]")),
    paste0("min_segment_bins := ", .depth_number(
      min_segment_bins, "min_segment_bins", function(x) x == floor(x) && x >= 0 && x <= 1e9,
      "a whole number in [0, 1e9]")),
    paste0("alt_fraction_threshold := ", .depth_number(
      alt_fraction_threshold, "alt_fraction_threshold", function(x) x >= 0 && x <= 1,
      "a number in [0, 1]")),
    paste0("close_margin := ", .depth_number(
      close_margin, "close_margin", function(x) x >= 0, "a number of at least 0")),
    paste0("max_bins := ", .depth_number(
      max_bins, "max_bins", function(x) x == floor(x) && x >= 1 && x <= 1e9,
      "a whole number in [1, 1e9]"))
  )
  paste0("SELECT * FROM ", macro, "(", paste(arguments, collapse = ", "), ")")
}

#' Decode Tumour Fraction, Ploidy and Copy-Number Segments from Corrected Bins
#'
#' Fits the copy-number state model of ichorCNA (Adalsteinsson, Ha, Freeman
#' et al. 2017; GavinHaLab/ichorCNA tag v0.6.1) to each sample's kept
#' corrected bins: a hidden Markov model over copy numbers 1 to `max_cn` with
#' a Student's t emission on the log ratio, whose state means follow from the
#' normal fraction `n` and tumour ploidy `phi`, fitted by expectation
#' maximisation with the upstream priors from every start of
#' `ploidies` x `normal_fractions`. Each start gives one solution, ranked by
#' its objective (data log-likelihood plus log prior).
#' [rduckhts_depth_state_solutions()] returns the solutions;
#' [rduckhts_depth_state_segments()] returns their Viterbi segments.
#'
#' The model learns from autosomal bins (contigs whose `duckhts_contig_key()`
#' is not X, Y or MT) and decodes every kept bin. A kept bin without a ratio
#' takes every state with likelihood 1, as upstream carries it. Upstream's
#' guard sets the fraction to zero (`status` `guard_zero`) when the longest
#' altered training segment has at most `min_segment_bins` bins and at most
#' `alt_fraction_threshold` of the training bins are altered; a solution ranked
#' below the best with an objective within `close_margin` of it is
#' `close_to_best`, the rest are `ok`. A sample without copy-number change
#' therefore reports no fraction, by the guard or by a fit whose normal
#' fraction is near 1; in that fit the state means lie within the noise of
#' each other and whole chromosomes can take any state, so its copy numbers
#' are not evidence. Fraction and ploidy trade against each other,
#' so all solutions are returned; a ploidy-3 explanation of a ploidy-2 genome
#' usually ranks second. A genome with one altered region does not identify
#' the fraction either: a one-copy gain in a share `m` of cells and a two-copy
#' gain in `m / 2` give the same ratio, and the prior decides, as in upstream.
#' The port and its kept upstream quirks are documented in the DuckHTS catalog
#' entry and `design/depth_fraction_lane.md`.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param bins_table Name of a table or view of corrected bins with
#'   `sample_id`, `chrom`, `start`, `end`, `log2_ratio` and `kept`, such as
#'   the output of [rduckhts_depth_bins()]. Restrict contigs and remove problem
#'   bins in SQL first.
#' @param normal_fractions Starts of the normal fraction, each in (0, 1).
#' @param ploidies Starts of the tumour ploidy, each in \[1, 10\]; the EM
#'   starts from their floor. A normal start of 0.95 runs only at ploidy 2,
#'   as upstream.
#' @param max_cn Highest copy-number state, a whole number in \[2, 20\].
#' @param txn_e Self-transition probability in (0, 1).
#' @param txn_strength Transition pseudo-counts, positive.
#' @param max_iter EM iteration cap, a whole number in \[2, 1000\].
#' @param min_segment_bins,alt_fraction_threshold The guard's two thresholds.
#' @param close_margin Objective difference below which a lower-ranked
#'   solution is `close_to_best`.
#' @param max_bins Kept bins one sample may hold; the model holds about
#'   `max_cn * 28 + 40` bytes a bin.
#' @inheritParams rduckhts_depth_bins
#' @return A data frame (or invisible `TRUE` when `table_name` is given).
#'   Solutions have one row per sample and start: `sample_id`, `solution`
#'   (rank, 1 is best), `start_normal`, `start_ploidy`, `normal_fraction`,
#'   `tumour_fraction`, `ploidy`, `log_likelihood`, `altered_fraction`,
#'   `iterations`, `status`, `segments`, `bins`, `upstream` and `method`.
#'   Segments have one row per sample, solution and segment: the solution
#'   columns through `status`, then `chrom`, `start`, `end`, `bins`, `state`
#'   (the copy number of the Viterbi state), `copy_number` (after upstream's
#'   integer correction of high-level amplifications), `median_log2_ratio`,
#'   `upstream` and `method`.
#' @examples
#' con <- rduckhts_connect()
#' DBI::dbExecute(con, paste(
#'   "CREATE TEMP TABLE bins AS SELECT 'S1' AS sample_id, CAST(c AS VARCHAR) AS chrom,",
#'   "b * 1000000 AS start, (b + 1) * 1000000 AS \"end\", true AS kept,",
#'   "log2(CASE WHEN c = 8 THEN 1.15 ELSE 1 END) +",
#'   "0.05 * (((b * 37 + c * 11) % 100) / 50.0 - 1) AS log2_ratio",
#'   "FROM range(1, 23) t(c), range(0, CASE WHEN c = 8 THEN 60 ELSE 20 END) r(b)"
#' ))
#' # Chromosome 8 has 60 bins so the gain passes the 50-bin guard. One start
#' # keeps the example short; the default grid has ten.
#' best <- rduckhts_depth_state_solutions(con, "bins", normal_fractions = 0.7, ploidies = 2)
#' best[, c("tumour_fraction", "ploidy", "status")]
#' segments <- rduckhts_depth_state_segments(con, "bins", normal_fractions = 0.7, ploidies = 2)
#' segments[segments$chrom == "8", c("chrom", "copy_number", "bins")]
#' DBI::dbDisconnect(con, shutdown = TRUE)
#' @export
rduckhts_depth_state_solutions <- function(con, bins_table,
                                           normal_fractions = c(0.5, 0.6, 0.7, 0.8, 0.9),
                                           ploidies = c(2, 3), max_cn = 7, txn_e = 0.9999999,
                                           txn_strength = 1e7, max_iter = 20, min_segment_bins = 50,
                                           alt_fraction_threshold = 0.05, close_margin = 10,
                                           max_bins = 1e6, table_name = NULL, overwrite = FALSE) {
  .somalier_validate_output(con, table_name, overwrite)
  query <- .depth_state_query(con, "duckhts_depth_state_solutions", bins_table, normal_fractions,
                              ploidies, max_cn, txn_e, txn_strength, max_iter, min_segment_bins,
                              alt_fraction_threshold, close_margin, max_bins)
  .somalier_publish_query(con, query, table_name, overwrite)
}

#' @rdname rduckhts_depth_state_solutions
#' @export
rduckhts_depth_state_segments <- function(con, bins_table,
                                          normal_fractions = c(0.5, 0.6, 0.7, 0.8, 0.9),
                                          ploidies = c(2, 3), max_cn = 7, txn_e = 0.9999999,
                                          txn_strength = 1e7, max_iter = 20, min_segment_bins = 50,
                                          alt_fraction_threshold = 0.05, close_margin = 10,
                                          max_bins = 1e6, table_name = NULL, overwrite = FALSE) {
  .somalier_validate_output(con, table_name, overwrite)
  query <- .depth_state_query(con, "duckhts_depth_state_segments", bins_table, normal_fractions,
                              ploidies, max_cn, txn_e, txn_strength, max_iter, min_segment_bins,
                              alt_fraction_threshold, close_margin, max_bins)
  .somalier_publish_query(con, query, table_name, overwrite)
}
