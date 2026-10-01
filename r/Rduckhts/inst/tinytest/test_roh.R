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
  expect_error(rduckhts_roh(con, vcf), "exactly one of af_tag and af_table")
  expect_error(rduckhts_roh(con, vcf, af_tag = "AF", af_table = "x"), "exactly one of af_tag and af_table")
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
