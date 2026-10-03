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
  aware <- unaware <- numeric(length(alphas))
  for (k in seq_along(alphas)) {
    dbWriteTable(con, "roh_counts", .roh_titration_counts(alphas[[k]]), overwrite = TRUE)
    aware[[k]] <- covered(rduckhts_roh_counts(con, "roh_counts", contamination = alphas[[k]]))
    unaware[[k]] <- covered(rduckhts_roh_counts(con, "roh_counts"))
  }
  # The truth run is found without contamination.
  expect_true(aware[[1]] > 0.8 * (truth_end - truth_start + 1))
  expect_equal(unaware[[1]], aware[[1]])
  # Ignoring 10% contamination erodes the run; modelling it keeps at least 90% of it.
  expect_true(unaware[[4]] < aware[[4]])
  expect_true(unaware[[4]] < unaware[[1]])
  expect_true(all(aware >= 0.9 * aware[[1]]))

  # The wrapper is the SQL macro, and NULL-free arguments are validated in R.
  sql <- dbGetQuery(con, "SELECT * FROM duckhts_roh_counts('roh_counts', contamination := 0.1)")
  wrapped <- rduckhts_roh_counts(con, "roh_counts", contamination = 0.1)
  expect_equal(.roh_observed(wrapped), .roh_observed(sql))
  expect_error(rduckhts_roh_counts(con, "roh_counts", seq_error = 0), "seq_error must be")
  expect_error(rduckhts_roh_counts(con, "roh_counts", seq_error = 0.5), "seq_error must be")
  expect_error(rduckhts_roh_counts(con, "roh_counts", contamination = 1), "contamination must be")
  expect_error(rduckhts_roh_counts(con, "roh_counts", contamination = -0.1), "contamination must be")
  expect_error(rduckhts_roh_counts(con, "roh_counts", contamination = NA_real_), "contamination must be")
  expect_error(rduckhts_roh_counts(con, NA_character_), "counts_table")
}
test_roh_counts()
