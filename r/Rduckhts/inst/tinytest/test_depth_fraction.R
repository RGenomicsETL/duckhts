library(tinytest)
library(DBI)

# An independent oracle in base R: the corrected bins, the region fraction and
# the sex-chromosome fraction computed from the formulas of the catalog
# entries, with R's quantile(type = 7), median() and arithmetic on vectors
# rather than SQL. The macros must agree with it to 1e-9, and on synthetic
# data with a known truth both must recover it within the declared tolerance.

set.seed(343)
contigs <- c(as.character(1:22), "X", "Y")
bins_per_contig <- c(rep(60L, 22L), 60L, 20L)
genome <- do.call(rbind, lapply(seq_along(contigs), function(i) {
  b <- seq_len(bins_per_contig[[i]]) - 1L
  data.frame(chrom = contigs[[i]], start = b * 1e6, end = (b + 1) * 1e6,
             gc = 0.35 + 0.25 * stats::runif(length(b)), stringsAsFactors = FALSE)
}))
# The last bin of chromosome 2 is half as long.
genome$end[genome$chrom == "2" & genome$start == 59e6] <- 59.5e6
# Copies of X, Y and chromosome 8: an XX genome, an XY genome, an XX genome
# with an XY second genome at 10%, and an XX genome with chromosome 8 gained in
# 20% of cells; then three XX and three XY reference genomes.
design <- data.frame(
  sample_id = c("xx", "xy", "xx_with_xy_10", "gain8_20", paste0("ref_xx_", 1:3), paste0("ref_xy_", 1:3)),
  x = c(2, 1, 1.9, 2, 2, 2, 2, 1, 1, 1), y = c(0, 1, 0.1, 0, 0, 0, 0, 1, 1, 1),
  c8 = c(2, 2, 2, 2.2, 2, 2, 2, 2, 2, 2), stringsAsFactors = FALSE)
bias <- function(gc) 1 + 3 * (gc - 0.45) - 20 * (gc - 0.45)^2
raw <- do.call(rbind, lapply(seq_len(nrow(design)), function(i) {
  copies <- ifelse(genome$chrom == "X", design$x[[i]], ifelse(genome$chrom == "Y", design$y[[i]],
            ifelse(genome$chrom == "8", design$c8[[i]], 2)))
  mean <- 700 * bias(genome$gc) * copies / 2 * (genome$end - genome$start) / 1e6
  data.frame(sample_id = design$sample_id[[i]], genome, reads = stats::rpois(nrow(genome), mean),
             stringsAsFactors = FALSE)
}))
# One bin of the first sample in each kind of rejection.
extra <- data.frame(sample_id = "xx", chrom = "1", start = c(100, 101, 102, 103) * 1e6,
                    end = c(101, 102, 103, 104) * 1e6, gc = c(NA, 0.45, 0.95, 0.45),
                    reads = c(700L, NA, 700L, 7000L), stringsAsFactors = FALSE)
raw <- rbind(raw, extra)

oracle_bins <- function(raw, gc_step = 0.01, min_stratum_bins = 10, q = 0.01, max_ratio = 4) {
  raw$autosomal <- !raw$chrom %in% c("X", "Y", "MT")
  raw$stratum <- floor(raw$gc / gc_step)
  raw$rate <- raw$reads / (raw$end - raw$start)
  raw$reason <- ifelse(is.na(raw$reads), "no_reads", ifelse(is.na(raw$gc), "no_gc", NA))
  out <- lapply(split(raw, raw$sample_id), function(s) {
    pool <- s[is.na(s$reason) & s$autosomal & s$reads > 0, ]
    limits <- stats::quantile(pool$rate, c(q, 1 - q), type = 7, names = FALSE)
    training <- pool[pool$reads > 0 & pool$rate >= limits[[1L]] & pool$rate <= limits[[2L]], ]
    strata <- tapply(training$rate, training$stratum, stats::median)
    sizes <- tapply(training$rate, training$stratum, length)
    key <- as.character(s$stratum)
    s$stratum_bins <- as.numeric(ifelse(is.na(s$stratum) | is.na(sizes[key]), 0, sizes[key]))
    s$reason <- ifelse(is.na(s$reason) & s$stratum_bins < min_stratum_bins, "sparse_stratum", s$reason)
    stratum_rate <- as.numeric(strata[key])
    centred <- is.na(s$reason) & s$autosomal & s$reads > 0
    centre <- stats::median(s$rate[centred] / stratum_rate[centred])
    s$expected <- ifelse(is.na(s$reason), stratum_rate * centre * (s$end - s$start), NA)
    s$ratio <- s$reads / s$expected
    s$reason <- ifelse(is.na(s$reason) & s$ratio > max_ratio, "high_ratio", s$reason)
    s$kept <- is.na(s$reason)
    s
  })
  do.call(rbind, out)
}

con <- rduckhts_connect()
dbWriteTable(con, "depth_raw", raw)

corrected <- rduckhts_depth_bins(con, "depth_raw")
expect_equal(nrow(corrected), nrow(raw))
expect_equal(names(corrected), c("sample_id", "chrom", "start", "end", "reads", "gc", "mappability",
                                 "gc_stratum", "stratum_bins", "expected", "ratio", "log2_ratio",
                                 "reason", "kept", "method"))
oracle <- oracle_bins(raw)
key <- function(x) paste(x$sample_id, x$chrom, x$start)
oracle <- oracle[match(key(corrected), key(oracle)), ]
expect_equal(corrected$reason, oracle$reason)
expect_equal(corrected$kept, oracle$kept)
expect_equal(as.numeric(corrected$stratum_bins), as.numeric(oracle$stratum_bins))
expect_equal(corrected$expected, oracle$expected, tolerance = 1e-9)
expect_equal(corrected$ratio, oracle$ratio, tolerance = 1e-9)
rejected <- corrected[corrected$sample_id == "xx" & corrected$start >= 100e6, ]
expect_equal(rejected$reason[order(rejected$start)], c("no_gc", "no_reads", "sparse_stratum", "high_ratio"))

# Hand-derived sparse samples: every positive bin has 100 reads and the
# same length and GC, so its expected count is 100 and its ratio is 1.
# Zero-read bins still have ratio 0 and remain available to estimators.
for (empty_bins in c(49L, 50L, 80L)) {
  sparse <- data.frame(sample_id = "sparse", chrom = "1", start = (0:99) * 1000,
                       end = (1:100) * 1000, gc = 0.4,
                       reads = c(rep(0L, empty_bins), rep(100L, 100L - empty_bins)))
  dbWriteTable(con, "depth_sparse", sparse, overwrite = TRUE)
  actual <- rduckhts_depth_bins(con, "depth_sparse")
  expected <- oracle_bins(sparse)
  expected <- expected[match(key(actual), key(expected)), ]
  expect_equal(nrow(actual), 100L)
  expect_true(all(actual$kept))
  expect_equal(actual$expected, rep(100, 100))
  expect_equal(actual$ratio, actual$reads / 100)
  expect_true(all(is.na(actual$log2_ratio[actual$reads == 0])))
  expect_equal(actual$expected, expected$expected, tolerance = 1e-9)
  expect_equal(actual$ratio, expected$ratio, tolerance = 1e-9)
}

# The same rows with a mappability column, and the table_name output.
dbWriteTable(con, "depth_raw_mappable", transform(raw, mappability = 1))
expect_true(rduckhts_depth_bins(con, "depth_raw_mappable", table_name = "depth_corrected"))
from_table <- dbGetQuery(con, "SELECT * FROM depth_corrected")
from_table <- from_table[match(key(corrected), key(from_table)), ]
expect_equal(from_table$expected, corrected$expected, tolerance = 1e-12)
expect_true(all(from_table$mappability == 1, na.rm = TRUE))

# The region fraction against the oracle, without and with the reference set.
dbExecute(con, "CREATE TEMP VIEW depth_tests AS SELECT * FROM depth_corrected WHERE sample_id NOT LIKE 'ref_%'")
dbExecute(con, paste(
  "CREATE TEMP VIEW depth_references AS SELECT *, CASE WHEN sample_id LIKE 'ref_xy_%'",
  "THEN 'XY' ELSE 'XX' END AS sex FROM depth_corrected WHERE sample_id LIKE 'ref_%'"))
dbExecute(con, paste(
  "CREATE TEMP TABLE depth_chromosomes AS SELECT CAST(c AS VARCHAR) AS region,",
  "CAST(c AS VARCHAR) AS chrom, 0 AS start, 99000000 AS \"end\" FROM range(1, 23) t(c)"))

oracle_shift <- function(s, region) {
  kept <- s[s$kept, ]
  inside <- kept$chrom == region & (kept$start + kept$end) / 2 < 99e6
  base <- kept[!inside & kept$autosomal, ]
  r <- kept[inside, ]
  rb <- sum(base$reads) / sum(base$expected)
  ratio <- (sum(r$reads) / sum(r$expected)) / rb
  phi <- max(1, sum((base$reads - rb * base$expected)^2 / (rb * base$expected)) / (nrow(base) - 1))
  se <- 2 * sqrt(phi * max(sum(r$reads), 1) / (sum(r$expected) * rb)^2 + ratio^2 * phi / max(sum(base$reads), 1))
  c(shift = 2 * (ratio - 1), se = se)
}
oracle_corrected <- oracle_bins(raw)
oracle_corrected$kept[is.na(oracle_corrected$kept)] <- FALSE

regions <- rduckhts_depth_region_fraction(con, "depth_tests", "depth_chromosomes")
expect_equal(nrow(regions), 4L * 22L)
for (sample in c("xx", "gain8_20")) {
  for (region in c("1", "8")) {
    row <- regions[regions$sample_id == sample & regions$region == region, ]
    expected <- oracle_shift(oracle_corrected[oracle_corrected$sample_id == sample, ], region)
    expect_equal(row$shift, expected[["shift"]], tolerance = 1e-9)
    expect_equal(row$shift_se, expected[["se"]], tolerance = 1e-9)
    expect_equal(row$shift_high - row$shift, 1.959963985 * row$shift_se, tolerance = 1e-9)
  }
}
# The planted gain within 0.04 (about four standard errors at 60 bins).
gain <- regions[regions$sample_id == "gain8_20" & regions$region == "8", ]
expect_true(abs(gain$shift - 0.2) < 0.04)
expect_true(gain$shift_low <= 0.2 && 0.2 <= gain$shift_high)

with_reference <- rduckhts_depth_region_fraction(con, "depth_tests", "depth_chromosomes",
                                                 reference_table = "depth_references")
row <- with_reference[with_reference$sample_id == "gain8_20" & with_reference$region == "8", ]
references <- vapply(paste0("ref_", c("xx_1", "xx_2", "xx_3", "xy_1", "xy_2", "xy_3")), function(id) {
  oracle_shift(oracle_corrected[oracle_corrected$sample_id == id, ], "8")
}, numeric(2))
own <- oracle_shift(oracle_corrected[oracle_corrected$sample_id == "gain8_20", ], "8")
expect_equal(row$reference_samples, 6)
expect_equal(row$reference_shift, mean(references["shift", ]), tolerance = 1e-9)
expect_equal(row$shift, own[["shift"]] - mean(references["shift", ]), tolerance = 1e-9)
expected_se <- sqrt(own[["se"]]^2 + max(0, stats::var(references["shift", ]) - mean(references["se", ]^2)) +
                      stats::var(references["shift", ]) / 6)
expect_equal(row$shift_se, expected_se, tolerance = 1e-9)
expect_equal(row$shift_high - row$shift, stats::qt(0.975, 5) * row$shift_se, tolerance = 1e-6)

# The sex-chromosome fraction against the oracle and the truth.
oracle_sex <- function(s) {
  kept <- s[s$kept, ]
  a <- kept[kept$autosomal, ]
  depth <- sum(a$reads) / sum(a$expected)
  phi <- max(1, sum((a$reads - depth * a$expected)^2 / (depth * a$expected)) / (nrow(a) - 1))
  copies <- function(chrom) {
    k <- kept[kept$chrom == chrom, ]
    value <- 2 * (sum(k$reads) / sum(k$expected)) / depth
    se <- sqrt((2 / (sum(k$expected) * depth))^2 * phi * max(sum(k$reads), 1) +
               value^2 * phi / max(sum(a$reads), 1))
    c(value = value, se = se)
  }
  x <- copies("X")
  y <- copies("Y")
  xx <- x[["value"]] >= 1.5
  fx <- if (xx) 2 - x[["value"]] else x[["value"]] - 1
  fy <- if (xx) y[["value"]] else 1 - y[["value"]]
  w <- 1 / x[["se"]]^2 + 1 / y[["se"]]^2
  c(x = x[["value"]], y = y[["value"]], fraction = (fx / x[["se"]]^2 + fy / y[["se"]]^2) / w,
    se = 1 / sqrt(w))
}
sex <- rduckhts_depth_sex_fraction(con, "depth_tests")
expect_equal(nrow(sex), 4L)
for (sample in c("xx", "xy", "xx_with_xy_10")) {
  row <- sex[sex$sample_id == sample, ]
  expected <- oracle_sex(oracle_corrected[oracle_corrected$sample_id == sample, ])
  expect_equal(row$x_copies, expected[["x"]], tolerance = 1e-9)
  expect_equal(row$y_copies, expected[["y"]], tolerance = 1e-9)
  expect_equal(row$fraction, expected[["fraction"]], tolerance = 1e-9)
  expect_equal(row$fraction_se, expected[["se"]], tolerance = 1e-9)
}
truth <- c(xx = 0, xy = 0, xx_with_xy_10 = 0.1)
for (sample in names(truth)) {
  expect_true(abs(sex$fraction[sex$sample_id == sample] - truth[[sample]]) < 0.03)
}
expect_equal(sex$major_sex[match(c("xx", "xy", "xx_with_xy_10"), sex$sample_id)], c("XX", "XY", "XX"))

with_sex_reference <- rduckhts_depth_sex_fraction(con, "depth_tests", reference_table = "depth_references")
expect_true(all(with_sex_reference$reference_samples == 6))
expect_true(abs(with_sex_reference$fraction[with_sex_reference$sample_id == "xx_with_xy_10"] - 0.1) < 0.03)

DBI::dbExecute(con, paste(
  "CREATE TEMP VIEW depth_mixed_references AS SELECT * FROM depth_references UNION ALL",
  "SELECT d.* REPLACE(s.id AS sample_id, kept AND CASE s.id",
  "WHEN 'unusable_no_autosomes' THEN chrom IN ('X', 'Y')",
  "WHEN 'unusable_no_sex' THEN chrom NOT IN ('X', 'Y')",
  "WHEN 'x_only' THEN chrom != 'Y' WHEN 'y_only' THEN chrom != 'X' END AS kept)",
  "FROM depth_references d CROSS JOIN (VALUES ('unusable_no_autosomes'),",
  "('unusable_no_sex'), ('x_only'), ('y_only')) s(id) WHERE d.sample_id = 'ref_xy_1'"
))
mixed_reference <- rduckhts_depth_sex_fraction(con, "depth_tests", reference_table = "depth_mixed_references")
expect_true(all(mixed_reference$reference_samples == 8))
self_reference <- rduckhts_depth_sex_fraction(con, "depth_mixed_references",
                                            reference_table = "depth_mixed_references")
expect_equal(nrow(self_reference), 10L)
expect_equal(as.numeric(self_reference$reference_samples),
             ifelse(startsWith(self_reference$sample_id, "unusable_"), 8, 7))

DBI::dbExecute(con, paste(
  "CREATE TEMP TABLE depth_large_positions AS SELECT 'S1' AS sample_id, chrom, start, \"end\",",
  "reads, 100.0 AS expected, true AS kept FROM (VALUES ('1', 0::BIGINT, 2::BIGINT, 100.0),",
  "('1', 2, 4, 100.0), ('2', 9223372036854775804, 9223372036854775807, 110.0))",
  "t(chrom, start, \"end\", reads)"
))
DBI::dbExecute(con, paste(
  "CREATE TEMP TABLE depth_large_regions AS SELECT * FROM (VALUES",
  "('inside', '2', 9223372036854775805::BIGINT, 9223372036854775806::BIGINT),",
  "('outside', '2', 9223372036854775806, 9223372036854775807)) t(region, chrom, start, \"end\")"
))
large <- rduckhts_depth_region_fraction(con, "depth_large_positions", "depth_large_regions")
expect_equal(as.numeric(large$bins[large$region == "inside"]), 1)
expect_equal(large$shift[large$region == "inside"], 0.2)
expect_equal(as.numeric(large$bins[large$region == "outside"]), 0)

# Argument checks at the wrapper boundary, and DuckDB's errors preserved.
expect_error(rduckhts_depth_bins(con, "depth_raw", gc_step = 0), "gc_step")
expect_error(rduckhts_depth_bins(con, "depth_raw", gc_step = 1e-10), "gc_step")
DBI::dbExecute(con, paste(
  "CREATE TEMP TABLE depth_gc_limit AS SELECT 'S1' AS sample_id, '1' AS chrom,",
  "0::BIGINT AS start, 100::BIGINT AS \"end\", 10::BIGINT AS reads, 1.0 AS gc"
))
expect_equal(rduckhts_depth_bins(con, "depth_gc_limit", gc_step = 1 / .Machine$integer.max,
                               min_stratum_bins = 1)$gc_stratum, .Machine$integer.max)
expect_error(rduckhts_depth_bins(con, "depth_raw", min_stratum_bins = 2.5), "min_stratum_bins")
expect_error(rduckhts_depth_bins(con, "depth_raw", outlier_quantile = 0.5), "outlier_quantile")
expect_error(rduckhts_depth_bins(con, "depth_raw", max_ratio = 1), "max_ratio")
expect_error(rduckhts_depth_bins(con, "depth_raw", min_mappability = NA), "min_mappability")
expect_error(rduckhts_depth_bins(con, c("a", "b")), "bins_table")
expect_error(rduckhts_depth_region_fraction(con, "depth_tests", NA_character_), "regions_table")
expect_error(rduckhts_depth_sex_fraction(con, "depth_tests", reference_table = 1), "reference_table")
dbExecute(con, paste(
  "CREATE TEMP VIEW depth_null_chrom AS SELECT * REPLACE",
  "(CASE WHEN chrom = '8' THEN NULL ELSE chrom END AS chrom) FROM depth_tests"
))
expect_error(rduckhts_depth_sex_fraction(con, "depth_null_chrom"), "chrom must not be NULL")
dbExecute(con, paste(
  "CREATE TEMP VIEW depth_null_reference_chrom AS SELECT * REPLACE",
  "(CASE WHEN chrom = '8' THEN NULL ELSE chrom END AS chrom) FROM depth_references"
))
expect_error(rduckhts_depth_sex_fraction(con, "depth_tests",
                                        reference_table = "depth_null_reference_chrom"),
             "chrom must not be NULL")
dbExecute(con, "CREATE TEMP VIEW depth_percent AS SELECT sample_id, chrom, start, \"end\", reads, gc * 100 AS gc FROM depth_raw")
expect_error(rduckhts_depth_bins(con, "depth_percent"), "gc must be a share")
dbExecute(con, "CREATE TEMP VIEW depth_bad_sex AS SELECT *, 'XXY' AS sex FROM depth_corrected WHERE sample_id LIKE 'ref_%'")
expect_error(rduckhts_depth_sex_fraction(con, "depth_tests", reference_table = "depth_bad_sex"),
             "reference sex must be XX or XY")

dbDisconnect(con, shutdown = TRUE)
