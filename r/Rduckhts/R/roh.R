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
#' Alternatively, `reference_table` and `proportions_table` select the
#' ancestry-tuned path. The reference is long-format with columns `chromosome`,
#' `position`, `allele_a`, `allele_b`, `group_id` and `frequency` (frequency of
#' `allele_b`); proportions have `sample_id`, `group_id` and `proportion`.
#' Their group sets must match exactly, and each VCF sample must have a
#' proportions row for every group. Per-site AF is the sum of each group's
#' frequency weighted by the sample's proportion, with the proportions divided
#' by their sum (which must be above 0 and at most 1, so `sum_to_one = FALSE`
#' results and rounded proportions are accepted). REF=`allele_a`, ALT=`allele_b`
#' uses that AF; reversed alleles use `1 - AF`. Reference alleles must be on
#' the forward strand of the VCF's assembly, as in a FASTA-anchored panel; no
#' strand flip is attempted, so palindromic (A/T, C/G) sites are oriented by REF
#' like any other, and a site whose alleles match neither order is not used.
#' Contig names are matched with `duckhts_contig_key()` once per distinct name
#' (one leading `chr` removed, M/MT written as MT, X and Y uppercased), so
#' `chr1` and `1` match, as do `chrX` and `X`; accessions, patches and numeric
#' sex chromosomes are not mapped. `af_clamp` limits nonzero-clamp frequencies to
#' `[af_clamp, 1-af_clamp]`; the default keeps zero population frequencies from
#' being treated as impossible, and zero disables clamping. Only called sites
#' are frequency-weighted. Supply exactly one frequency source: `af_tag`,
#' `af_table`, or both ancestry relations.
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
#' @param af_table Name of a table or view of site allele frequencies, instead
#'   of `af_tag` or the ancestry inputs.
#' @param reference_table Long-format reference with chromosome, position,
#'   allele_a, allele_b, group_id and frequency columns.
#' @param proportions_table Long-format sample ancestry proportions with
#'   sample_id, group_id and proportion columns. Its group set must exactly
#'   match `reference_table`.
#' @param af_clamp Clamp ancestry-tuned frequencies to `[af_clamp, 1-af_clamp]`;
#'   zero disables clamping. Must be in `[0, 0.5)`.
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
  samples = NULL, table_name = NULL, overwrite = FALSE,
  reference_table = NULL, proportions_table = NULL, af_clamp = 1e-3
) {
  .somalier_validate_output(con, table_name, overwrite)
  .somalier_scalar_text(path, "path")
  source <- .roh_frequency_source(con, af_tag, af_table, reference_table,
                                  proportions_table, af_clamp)
  options <- .roh_validated_options(genetic_map, hw_to_az, az_to_hw, gt_error,
                                   rec_rate, samples)
  arguments <- c(
    sql_quote_string(con, path), source$arguments,
    if (!is.null(options$genetic_map)) sql_quote_string(con, options$genetic_map),
    paste0("hw_to_az := ", options$hw_to_az),
    paste0("az_to_hw := ", options$az_to_hw),
    if (!is.null(options$gt_error)) paste0("gt_error := ", options$gt_error),
    if (!is.null(options$rec_rate)) paste0("rec_rate := ", options$rec_rate),
    if (source$ancestry) paste0("af_clamp := ", source$af_clamp),
    if (!is.null(options$samples)) {
      paste0("samples := ", sql_quote_string(con, options$samples))
    }
  )
  query <- paste0("SELECT * FROM ", source$macro, "(",
                  paste(arguments, collapse = ", "), ")")
  .somalier_publish_query(con, query, table_name, overwrite)
}

# A finite number in [0, maximum], formatted with full precision for SQL.
.roh_frequency_source <- function(con, af_tag, af_table, reference_table,
                                proportions_table, af_clamp) {
  ancestry <- !is.null(reference_table) || !is.null(proportions_table)
  if (ancestry && (is.null(reference_table) || is.null(proportions_table))) {
    stop("reference_table and proportions_table must be supplied together", call. = FALSE)
  }
  sources <- sum(!is.null(af_tag), !is.null(af_table), ancestry)
  if (sources != 1L) {
    stop("supply exactly one of af_tag, af_table, or ancestry relations", call. = FALSE)
  }
  if (!is.null(af_tag)) .somalier_scalar_text(af_tag, "af_tag")
  if (!is.null(af_table)) .somalier_validate_name(af_table, "af_table")
  if (ancestry) {
    .somalier_validate_name(reference_table, "reference_table")
    .somalier_validate_name(proportions_table, "proportions_table")
    af_clamp <- .roh_number(af_clamp, "af_clamp", maximum = 0.49999999999999994)
    arguments <- c(sql_quote_string(con, reference_table),
                   sql_quote_string(con, proportions_table))
    macro <- "duckhts_roh_ancestry"
  } else {
    arguments <- sql_quote_string(con, if (is.null(af_tag)) af_table else af_tag)
    macro <- if (is.null(af_tag)) "duckhts_roh_af_table" else "duckhts_roh"
  }
  list(ancestry = ancestry, af_clamp = af_clamp, arguments = arguments, macro = macro)
}

.roh_validated_options <- function(genetic_map, hw_to_az, az_to_hw, gt_error,
                                   rec_rate, samples) {
  if (!is.null(genetic_map)) .somalier_validate_name(genetic_map, "genetic_map")
  if (!is.null(samples)) .somalier_scalar_text(samples, "samples", allow_empty = TRUE)
  list(
    genetic_map = genetic_map,
    hw_to_az = .roh_number(hw_to_az, "hw_to_az", maximum = 1),
    az_to_hw = .roh_number(az_to_hw, "az_to_hw", maximum = 1),
    gt_error = if (is.null(gt_error)) NULL else .roh_number(gt_error, "gt_error"),
    rec_rate = if (is.null(rec_rate)) NULL else .roh_number(rec_rate, "rec_rate"),
    samples = samples
  )
}

.roh_number <- function(value, name, maximum = Inf) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0 || value > maximum) {
    stop(name, " must be one finite number in [0, ",
         if (is.finite(maximum)) format(maximum) else "Inf", "]", call. = FALSE)
  }
  sprintf("%.17g", value)
}
