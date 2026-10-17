library(tinytest)
library(DBI)

# Expected segments are the RG lines of `bcftools roh` 1.23.1-70-g6dbd8fef on
# roh_fixture.vcf.gz (see test/sql/roh.test and test/scripts/roh_bcftools_expected.sh).
# Quality is compared at the one decimal bcftools prints.
.roh_expected <- function(rows) {
  out <- do.call(rbind, lapply(rows, function(row) {
    data.frame(sample = row[[1]], chrom = row[[2]], start = row[[3]], end = row[[4]],
               length = row[[5]], n_markers = row[[6]], quality = row[[7]],
               stringsAsFactors = FALSE)
  }))
  out[order(out$sample, out$chrom, out$start), ]
}

.roh_observed <- function(roh) {
  roh$quality <- round(roh$quality, 1)
  roh <- roh[order(roh$sample, roh$chrom, roh$start), ]
  rownames(roh) <- NULL
  roh
}

.roh_compare <- function(observed, expected) {
  rownames(expected) <- NULL
  expect_equal(nrow(observed), nrow(expected))
  expect_equal(observed$sample, expected$sample)
  expect_equal(observed$chrom, expected$chrom)
  expect_equal(as.numeric(observed$start), as.numeric(expected$start))
  expect_equal(as.numeric(observed$end), as.numeric(expected$end))
  expect_equal(as.numeric(observed$length), as.numeric(expected$length))
  expect_equal(as.integer(observed$n_markers), as.integer(expected$n_markers))
  expect_equal(observed$quality, expected$quality, tolerance = 0, scale = 1)
}

test_roh_parity <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  vcf <- system.file("extdata", "roh_fixture.vcf.gz", package = "Rduckhts")
  af_file <- system.file("extdata", "roh_af.tsv.gz", package = "Rduckhts")
  map1 <- system.file("extdata", "roh_map_chr1.txt", package = "Rduckhts")
  map2 <- system.file("extdata", "roh_map_chr2.txt", package = "Rduckhts")

  dbExecute(con, sprintf(paste(
    "CREATE TABLE roh_af AS SELECT chrom, pos, split_part(alleles, ',', 1) AS ref,",
    "split_part(alleles, ',', 2) AS alt, af FROM read_csv(%s, delim = '\\t', header = false,",
    "columns = {'chrom': 'VARCHAR', 'pos': 'BIGINT', 'alleles': 'VARCHAR', 'af': 'DOUBLE'})"),
    dbQuoteString(con, af_file)))
  dbExecute(con, paste(
    "CREATE TABLE roh_ancestry_reference AS SELECT chrom AS chromosome, pos AS position,",
    "alt AS allele_a, ref AS allele_b, 'A' AS group_id, 0.2::DOUBLE AS frequency FROM roh_af",
    "WHERE NOT ((ref='A' AND alt='T') OR (ref='T' AND alt='A') OR (ref='C' AND alt='G') OR (ref='G' AND alt='C'))",
    "UNION ALL SELECT chrom, pos, alt, ref, 'B', 0.4::DOUBLE FROM roh_af",
    "WHERE NOT ((ref='A' AND alt='T') OR (ref='T' AND alt='A') OR (ref='C' AND alt='G') OR (ref='G' AND alt='C'))"))
  dbExecute(con, paste(
    "CREATE TABLE roh_ancestry_proportions AS SELECT * FROM (VALUES",
    "('S1', 'A', 1.0), ('S1', 'B', 0.0), ('S2', 'A', 1.0), ('S2', 'B', 0.0),",
    "('S3', 'A', 1.0), ('S3', 'B', 0.0), ('S4', 'A', 1.0), ('S4', 'B', 0.0)",
    ") AS p(sample_id, group_id, proportion)"))
  dbExecute(con, paste(
    "CREATE TABLE roh_ancestry_af AS SELECT chrom, pos, ref, alt, 0.8::DOUBLE AS af FROM roh_af",
    "WHERE NOT ((ref='A' AND alt='T') OR (ref='T' AND alt='A') OR (ref='C' AND alt='G') OR (ref='G' AND alt='C'))"))
  dbExecute(con, sprintf(paste(
    "CREATE TABLE roh_map AS SELECT 'chr1' AS chrom, position AS pos, \"Genetic_Map(cM)\" AS cm",
    "FROM read_csv(%s, delim = ' ') UNION ALL SELECT 'chr2', position, \"Genetic_Map(cM)\"",
    "FROM read_csv(%s, delim = ' ')"), dbQuoteString(con, map1), dbQuoteString(con, map2)))

  # PL emissions, frequencies from INFO/AF.
  .roh_compare(.roh_observed(rduckhts_roh(con, vcf, af_tag = "AF")), .roh_expected(list(
    list("S1", "chr1", 1042655, 2396807, 1354153, 113, 32.0),
    list("S4", "chr1", 1876100, 3319712, 1443613, 131, 44.8),
    list("S2", "chr2", 693881, 1722963, 1029083, 114, 32.9),
    list("S4", "chr2", 23420, 2400223, 2376804, 239, 46.6))))

  # GT emissions (-G30).
  .roh_compare(.roh_observed(rduckhts_roh(con, vcf, af_tag = "AF", gt_error = 30)), .roh_expected(list(
    list("S1", "chr1", 1042655, 2396807, 1354153, 112, 40.6),
    list("S4", "chr1", 1876100, 3319712, 1443613, 131, 56.0),
    list("S2", "chr2", 693881, 1722963, 1029083, 114, 39.6),
    list("S4", "chr2", 23420, 2400223, 2376804, 239, 51.7))))

  # Frequencies from a relation (--AF-file).
  .roh_compare(.roh_observed(rduckhts_roh(con, vcf, af_table = "roh_af")), .roh_expected(list(
    list("S1", "chr1", 1042655, 2396807, 1354153, 113, 32.0),
    list("S4", "chr1", 1876100, 3319712, 1443613, 131, 44.9),
    list("S2", "chr2", 693881, 1722963, 1029083, 114, 33.2),
    list("S4", "chr2", 23420, 2400223, 2376804, 238, 46.5))))

  # Ancestry-weighted AF is group A's 0.2, reversed to 0.8 for these sites.
  ancestry <- rduckhts_roh(con, vcf, reference_table = "roh_ancestry_reference",
                           proportions_table = "roh_ancestry_proportions", gt_error = 30,
                           samples = "S1")
  equivalent <- rduckhts_roh(con, vcf, af_table = "roh_ancestry_af", gt_error = 30,
                             samples = "S1")
  .roh_compare(.roh_observed(ancestry), .roh_observed(equivalent))
  ancestry_map <- rduckhts_roh(con, vcf, reference_table = "roh_ancestry_reference",
                               proportions_table = "roh_ancestry_proportions",
                               genetic_map = "roh_map", gt_error = 30, samples = "S1")
  equivalent_map <- rduckhts_roh(con, vcf, af_table = "roh_ancestry_af", genetic_map = "roh_map",
                                 gt_error = 30, samples = "S1")
  .roh_compare(.roh_observed(ancestry_map), .roh_observed(equivalent_map))

  # bcftools 1.24: roh --AF-file roh_af.tsv.gz -m roh_map_{CHROM}.txt.
  .roh_compare(.roh_observed(rduckhts_roh(
    con, vcf, af_table = "roh_af", site_rule = "bcftools", genetic_map = "roh_map"
  )), .roh_expected(list(
    list("S1", "chr1", 1432927, 2396807, 963881, 95, 42.4),
    list("S4", "chr1", 2246017, 3319712, 1073696, 121, 77.4),
    list("S2", "chr2", 693881, 1722963, 1029083, 114, 18.6),
    list("S4", "chr2", 23420, 2400223, 2376804, 238, 94.0))))

  # Constant recombination rate (-M 1e-6).
  .roh_compare(.roh_observed(rduckhts_roh(con, vcf, af_tag = "AF", rec_rate = 1e-6)), .roh_expected(list(
    list("S1", "chr1", 1432927, 2396807, 963881, 95, 49.7),
    list("S4", "chr1", 1876100, 3319712, 1443613, 131, 66.6),
    list("S2", "chr2", 693881, 1722963, 1029083, 114, 40.7),
    list("S4", "chr2", 23420, 2400223, 2376804, 239, 79.2))))

  # Genetic map (-m).
  .roh_compare(.roh_observed(rduckhts_roh(con, vcf, af_tag = "AF", genetic_map = "roh_map")), .roh_expected(list(
    list("S1", "chr1", 1432927, 2396807, 963881, 95, 41.4),
    list("S4", "chr1", 2246017, 3319712, 1073696, 121, 77.5),
    list("S2", "chr2", 693881, 1722963, 1029083, 114, 17.7),
    list("S4", "chr2", 23420, 2400223, 2376804, 239, 94.4))))

  # Sample selection, and a table target.
  .roh_compare(.roh_observed(rduckhts_roh(con, vcf, af_tag = "AF", samples = "S2,S4")), .roh_expected(list(
    list("S4", "chr1", 1876100, 3319712, 1443613, 131, 44.8),
    list("S2", "chr2", 693881, 1722963, 1029083, 114, 32.9),
    list("S4", "chr2", 23420, 2400223, 2376804, 239, 46.6))))
  expect_true(rduckhts_roh(con, vcf, af_tag = "AF", table_name = "roh_out"))
  expect_equal(dbGetQuery(con, "SELECT count(*)::INTEGER AS n FROM roh_out")$n, 4L)
  expect_error(rduckhts_roh(con, vcf, af_tag = "AF", table_name = "roh_out"))
  expect_true(rduckhts_roh(con, vcf, af_tag = "AF", table_name = "roh_out", overwrite = TRUE))
}
test_roh_parity()

test_roh_degenerate_inputs <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  empty <- system.file("extdata", "roh_empty.vcf", package = "Rduckhts")
  single <- system.file("extdata", "roh_single_site.vcf", package = "Rduckhts")
  gt_only <- system.file("extdata", "roh_gt_only.vcf", package = "Rduckhts")

  expect_equal(nrow(rduckhts_roh(con, empty, af_tag = "AF")), 0L)
  one <- rduckhts_roh(con, single, af_tag = "AF")
  expect_equal(nrow(one), 1L)
  expect_equal(c(one$start, one$end, one$length, one$n_markers), c(5000, 5000, 1, 1))
  expect_equal(round(one$quality, 1), 13.1)

  # GT-only VCF (no FORMAT/PL): bcftools roh -G30 gives S1 20000-70000, 51 markers, 43.9.
  gt <- rduckhts_roh(con, gt_only, af_tag = "AF", gt_error = 30)
  expect_equal(nrow(gt), 1L)
  expect_equal(c(gt$start, gt$end, gt$length, gt$n_markers), c(20000, 70000, 50001, 51))
  expect_equal(round(gt$quality, 1), 43.9)
}
test_roh_degenerate_inputs()

test_roh_validation <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  vcf <- system.file("extdata", "roh_fixture.vcf.gz", package = "Rduckhts")
  expect_error(rduckhts_roh(con, vcf), "exactly one of af_tag, af_table, or ancestry relations")
  expect_error(rduckhts_roh(con, vcf, af_tag = "AF", af_table = "x"), "exactly one of af_tag, af_table, or ancestry relations")
  expect_error(rduckhts_roh(con, vcf, reference_table = "x"), "must be supplied together")
  expect_error(rduckhts_roh(con, vcf, reference_table = "x", proportions_table = "y",
                            af_clamp = 0.5), "af_clamp must be")
  expect_error(rduckhts_roh(con, vcf, af_tag = "AF", hw_to_az = 2), "hw_to_az must be")
  expect_error(rduckhts_roh(con, vcf, af_tag = "AF", az_to_hw = -1), "az_to_hw must be")
  expect_error(rduckhts_roh(con, vcf, af_tag = "AF", gt_error = -1), "gt_error must be")
  expect_error(rduckhts_roh(con, vcf, af_tag = "AF", rec_rate = NA_real_), "rec_rate must be")
  expect_error(rduckhts_roh(con, NA_character_, af_tag = "AF"), "path must be")
  # The memory limits are whole numbers of at least 1, checked in R, and the
  # extension reports a group or a total that exceeds them.
  expect_error(rduckhts_roh(con, vcf, af_tag = "AF", max_sites = 0), "max_sites must be")
  expect_error(rduckhts_roh(con, vcf, af_tag = "AF", max_sites = 1e9), "max_sites must be at most")
  expect_error(rduckhts_roh(con, vcf, af_tag = "AF", max_site_bytes = 1.5),
               "max_site_bytes must be")
  expect_error(rduckhts_roh(con, vcf, af_tag = "AF", max_sites = 3), "max_sites")
  expect_error(rduckhts_roh(con, vcf, af_tag = "AF", max_site_bytes = 100), "max_site_bytes")
  expect_equal(.roh_observed(rduckhts_roh(con, vcf, af_tag = "AF", max_sites = 1e6,
                                          max_site_bytes = 1e6)),
               .roh_observed(rduckhts_roh(con, vcf, af_tag = "AF")))
  # The kernel rejects unsorted positions and mismatched lists by name.
  expect_error(dbGetQuery(con, paste(
    "SELECT duckhts_roh_segments([200, 100], [0.3, 0.3], [[0, 30, 60], [0, 30, 60]],",
    "NULL::BIGINT[], NULL::DOUBLE[], NULL::DOUBLE, 6.7e-8, 5e-9)")),
    "positions must be sorted ascending")
  expect_error(dbGetQuery(con, paste(
    "SELECT duckhts_roh_segments([100, 200], [0.3], [[0, 30, 60], [0, 30, 60]],",
    "NULL::BIGINT[], NULL::DOUBLE[], NULL::DOUBLE, 6.7e-8, 5e-9)")),
    "differ in length")
}
test_roh_validation()

# Read-count emissions: a synthetic contamination titration with sampled reads.
# One chromosome of 400 sites 12.5 kb apart, homozygous over sites 121-280 and drawn
# from Hardy-Weinberg proportions elsewhere; depths are Poisson(30) and each read comes
# from a contaminating individual of the same population with probability alpha.
# Ignoring contamination must shorten the runs; modelling it must recover them.
.roh_titration_counts <- function(alpha, seed = 318L) {
  set.seed(seed)
  n <- 400L
  af <- stats::runif(n, 0.1, 0.9)
  draw <- function() stats::rbinom(n, 1L, af) + stats::rbinom(n, 1L, af)
  genotype <- draw()
  inside <- seq_len(n) %in% 121:280
  genotype[inside] <- 2L * stats::rbinom(sum(inside), 1L, af[inside])
  contaminant <- draw()
  depth <- stats::rpois(n, 30)
  error <- 1e-3
  q <- c(error, 0.5, 1 - error)
  from_contaminant <- stats::rbinom(n, depth, alpha)
  counted <- stats::rbinom(n, depth - from_contaminant, q[genotype + 1L]) +
    stats::rbinom(n, from_contaminant, q[contaminant + 1L])
  data.frame(sample_id = "S1", chrom = "1", pos = seq_len(n) * 12500L,
             ref_count = depth - counted, alt_count = counted, af = af)
}

test_roh_counts <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  truth_start <- 121 * 12500
  truth_end <- 280 * 12500
  covered <- function(roh) {
    if (!nrow(roh)) return(0)
    sum(pmax(0, pmin(roh$end, truth_end) - pmax(roh$start, truth_start) + 1))
  }
  alphas <- c(0, 0.02, 0.05, 0.10)
  aware <- unaware <- aware_outside <- numeric(length(alphas))
  for (k in seq_along(alphas)) {
    dbWriteTable(con, "roh_counts", .roh_titration_counts(alphas[[k]]), overwrite = TRUE)
    calls <- rduckhts_roh_counts(con, "roh_counts", contamination = alphas[[k]])
    aware[[k]] <- covered(calls)
    aware_outside[[k]] <- sum(calls$length) - aware[[k]]
    unaware[[k]] <- covered(rduckhts_roh_counts(con, "roh_counts"))
  }
  # The truth run is found without contamination.
  expect_true(aware[[1]] > 0.8 * (truth_end - truth_start + 1))
  expect_equal(unaware[[1]], aware[[1]])
  # Ignoring 10% contamination erodes the run; modelling it keeps at least 90% of it.
  expect_true(unaware[[4]] < aware[[4]])
  expect_true(unaware[[4]] < unaware[[1]])
  expect_true(all(aware >= 0.9 * aware[[1]]))
  # Recovery is not bought by overcalling: calls outside the truth run stay under
  # 5% of its length at every contamination level.
  expect_true(all(aware_outside <= 0.05 * (truth_end - truth_start + 1)))

  # The wrapper is the SQL macro, and NULL-free arguments are validated in R.
  sql <- dbGetQuery(con, "SELECT * FROM duckhts_roh_counts('roh_counts', contamination := 0.1)")
  wrapped <- rduckhts_roh_counts(con, "roh_counts", contamination = 0.1)
  expect_equal(.roh_observed(wrapped), .roh_observed(sql))
  expect_error(rduckhts_roh_counts(con, "roh_counts", seq_error = 0), "seq_error must be")
  expect_error(rduckhts_roh_counts(con, "roh_counts", seq_error = 0.5), "seq_error must be")
  expect_error(rduckhts_roh_counts(con, "roh_counts", contamination = 1), "contamination must be")
  expect_error(rduckhts_roh_counts(con, "roh_counts", contamination = -0.1), "contamination must be")
  expect_error(rduckhts_roh_counts(con, "roh_counts", contamination = NA_real_), "contamination must be")
  # allele_balance and artefact_weight pass through to the macro and are validated in R.
  sql <- dbGetQuery(con, paste("SELECT * FROM duckhts_roh_counts('roh_counts', contamination := 0.1,",
                               "allele_balance := 0.48, artefact_weight := 0.003)"))
  wrapped <- rduckhts_roh_counts(con, "roh_counts", contamination = 0.1,
                                 allele_balance = 0.48, artefact_weight = 0.003)
  expect_true(nrow(wrapped) > 0)
  expect_equal(.roh_observed(wrapped), .roh_observed(sql))
  # The defaults are the decode without the two terms, quality bits included.
  sorted <- function(rows) {
    rows <- rows[order(rows$sample, rows$chrom, rows$start), ]
    rownames(rows) <- NULL
    rows
  }
  expect_identical(sorted(rduckhts_roh_counts(con, "roh_counts", allele_balance = 0.5,
                                              artefact_weight = 0)),
                   sorted(dbGetQuery(con, paste("SELECT * FROM duckhts_roh_counts('roh_counts',",
                                                "seq_error := 1e-3, contamination := 0.0)"))))
  for (bad in list(0.349, 0.651, NA_real_, c(0.5, 0.5), "0.5")) {
    expect_error(rduckhts_roh_counts(con, "roh_counts", allele_balance = bad), "allele_balance must be")
  }
  for (bad in list(-1e-9, 0.0501, NA_real_, c(0, 0), "0")) {
    expect_error(rduckhts_roh_counts(con, "roh_counts", artefact_weight = bad), "artefact_weight must be")
  }
  expect_error(rduckhts_roh_counts(con, NA_character_), "counts_table")
  expect_error(rduckhts_roh_counts(con, "roh_counts", max_sites = 0), "max_sites must be")
  expect_error(rduckhts_roh_counts(con, "roh_counts", max_sites = 5), "max_sites")
  expect_error(rduckhts_roh_counts(con, "roh_counts", max_site_bytes = 100), "max_site_bytes")
}
test_roh_counts()

test_roh_counts_fitted <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  counts <- .roh_titration_counts(0.1)
  counts$sample_id <- "S1"
  second <- counts
  second$sample_id <- "S2"
  dbWriteTable(con, "counts", rbind(counts, second))
  # Columns other than the four fitted numbers and artefact_weight_sd (status
  # here) are ignored.
  fits <- data.frame(sample = c("S1", "S2"), seq_error = 0.001, contamination = 0.1,
                     status = "caller_selected", allele_balance = 0.5, artefact_weight = 0,
                     artefact_weight_sd = 0)
  dbWriteTable(con, "fits", fits)
  # The fitted decode adds stable_fraction and gained_bases to the columns of the scalar one;
  # the run columns are compared without it.
  ordered <- function(rows) {
    rows <- rows[order(rows$sample, rows$chrom, rows$start), setdiff(names(rows), c("stable_fraction", "gained_bases"))]
    rownames(rows) <- NULL
    rows
  }
  scalar <- rduckhts_roh_counts(con, "counts", contamination = 0.1)
  fitted <- rduckhts_roh_counts(con, "counts", fit_table = "fits")
  expect_identical(names(fitted), c(names(scalar), "stable_fraction", "gained_bases"))
  expect_identical(ordered(fitted), ordered(scalar))
  # With artefact_weight_sd 0 the three decodes are the same: every run has 1.
  expect_true(nrow(fitted) > 0)
  expect_true(all(fitted$stable_fraction == 1))
  expect_true(all(fitted$gained_bases == 0))
  # Unused scalar arguments do not constrain fitted parameters.
  expect_identical(ordered(rduckhts_roh_counts(
    con, "counts", fit_table = "fits", seq_error = NA, contamination = -1,
    allele_balance = 2, artefact_weight = -1
  )), ordered(scalar))

  fits$seq_error[[2]] <- 0.02
  fits$contamination[[2]] <- 0.2
  fits$allele_balance[[2]] <- 0.48
  fits$artefact_weight[[2]] <- 0.003
  fits$artefact_weight_sd[[2]] <- 0.001
  dbWriteTable(con, "fits", fits, overwrite = TRUE)
  dbWriteTable(con, "map", data.frame(chrom = "1", pos = c(1, 5000000), cm = c(0, 5)))
  for (map in list(NULL, "map")) {
    fitted <- rduckhts_roh_counts(con, "counts", fit_table = "fits", genetic_map = map)
    scalar <- lapply(seq_len(nrow(fits)), function(i) {
      rows <- rduckhts_roh_counts(con, "counts", genetic_map = map,
                                  seq_error = fits$seq_error[[i]],
                                  contamination = fits$contamination[[i]],
                                  allele_balance = fits$allele_balance[[i]],
                                  artefact_weight = fits$artefact_weight[[i]])
      rows[rows$sample == fits$sample[[i]], ]
    })
    expect_equal(length(unique(fitted$sample)), 2L)
    expect_identical(ordered(fitted), ordered(do.call(rbind, scalar)))
    expect_true(all(fitted$stable_fraction >= 0 & fitted$stable_fraction <= 1))
    expect_true(all(fitted$gained_bases >= 0))
  }
  # A NULL artefact_weight_sd counts as 0; a negative one is an error.
  no_interval <- fits
  no_interval$artefact_weight_sd <- NA_real_
  dbWriteTable(con, "no_interval_fits", no_interval)
  no_interval_runs <- rduckhts_roh_counts(con, "counts", fit_table = "no_interval_fits")
  expect_true(all(no_interval_runs$stable_fraction == 1))
  expect_true(all(no_interval_runs$gained_bases == 0))
  no_interval$artefact_weight_sd[[2]] <- -0.001
  dbWriteTable(con, "no_interval_fits", no_interval, overwrite = TRUE)
  expect_error(rduckhts_roh_counts(con, "counts", fit_table = "no_interval_fits"),
               "artefact_weight_sd must be finite and at least 0")
  expect_error(rduckhts_roh_counts(con, "counts", fit_table = NA_character_), "fit_table")
  expect_error(rduckhts_roh_counts(con, "counts", fit_table = c("fits", "fits")), "fit_table")
  dbWriteTable(con, "bad_fits", fits[1, ], overwrite = TRUE)
  expect_error(rduckhts_roh_counts(con, "counts", fit_table = "bad_fits"),
               "missing fit row for sample S2")
  out_of_range <- c(seq_error = 0.5, contamination = 1, allele_balance = 0.34, artefact_weight = 0.06)
  for (parameter in names(out_of_range)) {
    bad <- fits
    bad[[parameter]][[2]] <- NA_real_
    dbWriteTable(con, "bad_fits", bad, overwrite = TRUE)
    expect_error(rduckhts_roh_counts(con, "counts", fit_table = "bad_fits"),
                 paste("NULL", parameter, "for sample S2"))
    bad[[parameter]][[2]] <- out_of_range[[parameter]]
    dbWriteTable(con, "bad_fits", bad, overwrite = TRUE)
    expect_error(rduckhts_roh_counts(con, "counts", fit_table = "bad_fits"),
                 paste(parameter, "must be in"))
  }
  dbWriteTable(con, "bad_fits", rbind(fits, fits), overwrite = TRUE)
  # Duplicates are counted, so the error does not depend on a DuckDB setting.
  expect_error(rduckhts_roh_counts(con, "counts", fit_table = "bad_fits"),
               "duplicate fit rows for sample")
}
test_roh_counts_fitted()

# Equal 0.5 transitions make sites independent: homozygous pairs form AZ runs,
# and every third (balanced) site is HW. Lengths 1..9 cover partial traceback bytes.
test_roh_short_traceback <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  for (n in seq_len(9L)) {
    pos <- seq_len(n)
    counts <- data.frame(sample_id = "S1", chrom = "1", pos = pos,
                         ref_count = ifelse(pos %% 3L == 0L, 15L, 30L),
                         alt_count = ifelse(pos %% 3L == 0L, 15L, 0L), af = 0.5)
    dbWriteTable(con, "short_counts", counts, overwrite = TRUE)
    runs <- rduckhts_roh_counts(con, "short_counts", hw_to_az = 0.5, az_to_hw = 0.5)
    runs <- runs[order(runs$start), ]
    starts <- seq.int(1L, n, by = 3L)
    ends <- pmin(starts + 1L, n)
    expect_equal(as.integer(runs$start), starts)
    expect_equal(as.integer(runs$end), ends)
    expect_equal(runs$n_markers, ends - starts + 1L)
  }
}
test_roh_short_traceback()

# site_rule on roh_repeated.vcf.gz: each used test position is a run of one
# marker. "bcftools" gives the runs of bcftools roh 1.24 --AF-file -a 0.01
# -H 0.01 (see test/sql/roh.test); the default keeps matching records.
test_roh_site_rule <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  path <- system.file("extdata", "roh_repeated.vcf.gz", package = "Rduckhts")
  af_file <- system.file("extdata", "roh_repeated_af.tsv.gz", package = "Rduckhts")
  expect_true(nzchar(path) && nzchar(af_file))
  dbExecute(con, sprintf(paste(
    "CREATE TEMP TABLE repeated_af AS SELECT chrom, pos, split_part(alleles, ',', 1) AS ref,",
    "split_part(alleles, ',', 2) AS alt, af FROM read_csv(%s, delim = '\t', header = false,",
    "columns = {'chrom': 'VARCHAR', 'pos': 'BIGINT', 'alleles': 'VARCHAR', 'af': 'DOUBLE'})"),
    dbQuoteString(con, af_file)))
  dbExecute(con, "CREATE TEMP TABLE repeated_af_first AS SELECT * FROM repeated_af WHERE NOT (pos = 4000 AND alt = 'G')")

  starts <- function(...) sort(as.integer(rduckhts_roh(con, path, hw_to_az = 0.01, az_to_hw = 0.01, ...)$start))
  expect_equal(starts(af_table = "repeated_af_first", site_rule = "bcftools"),
               c(500L, 1000L, 5000L, 6000L, 9000L))
  expect_equal(starts(af_table = "repeated_af"),
               c(500L, 1000L, 2000L, 4000L, 5000L, 6000L, 7000L, 8000L, 9000L, 10000L))
  expect_error(starts(af_table = "repeated_af", site_rule = "bcftools"),
               "needs one af_table row per position, but chr1:4000 has 2 rows")
  for (bad in list("first", NA_character_, c("matching", "bcftools"), 1)) {
    expect_error(starts(af_table = "repeated_af_first", site_rule = bad),
                 "site_rule must be \"matching\" or \"bcftools\"")
  }
  expect_error(starts(af_tag = "AF", site_rule = "bcftools"), "site_rule = \"bcftools\" needs af_table")

  # af_order: the least value at a position stands for the first line of the
  # AF file. `line` is the order of roh_repeated_af.tsv.gz (4000 A,T before
  # A,G); `line_swapped` swaps its two 4000 lines. bcftools 1.24 gives a run
  # at 4000 only for the swapped file (see test/sql/roh.test).
  dbExecute(con, paste(
    "CREATE TEMP TABLE repeated_af_lines AS SELECT *,",
    "CASE WHEN pos = 4000 THEN 19 - line ELSE line END AS line_swapped",
    "FROM (SELECT *, row_number() OVER () AS line FROM repeated_af)"))
  expect_equal(dbGetQuery(con, "SELECT alt FROM repeated_af_lines WHERE pos = 4000 ORDER BY line")$alt,
               c("T", "G"))
  expect_equal(starts(af_table = "repeated_af_lines", site_rule = "bcftools", af_order = "line"),
               c(500L, 1000L, 5000L, 6000L, 9000L))
  expect_equal(starts(af_table = "repeated_af_lines", site_rule = "bcftools", af_order = "line_swapped"),
               c(500L, 1000L, 4000L, 5000L, 6000L, 9000L))
  expect_error(starts(af_table = "repeated_af_lines", site_rule = "bcftools", af_order = "line_number"),
               "af_order 'line_number' is not a column of af_table")
  dbExecute(con, "CREATE TEMP TABLE repeated_af_tied AS SELECT * REPLACE (9 AS line) FROM repeated_af_lines")
  expect_error(starts(af_table = "repeated_af_tied", site_rule = "bcftools", af_order = "line"),
               "chr1:4000 repeats a value of line")
  for (bad in list(NA_character_, "", c("line", "line_swapped"), 1)) {
    expect_error(starts(af_table = "repeated_af_lines", site_rule = "bcftools", af_order = bad),
                 "af_order must be one nonempty string")
  }
  expect_error(starts(af_table = "repeated_af_lines", af_order = "line"),
               "af_order needs site_rule = \"bcftools\"")
}
test_roh_site_rule()
