library(tinytest)
library(DBI)

test_record_major_genotypes <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, "SET threads=4")
  fixture <- function(name) {
    path <- system.file("extdata", name, package = "Rduckhts")
    stopifnot(nzchar(path), file.exists(path))
    path
  }
  quote <- function(text) as.character(dbQuoteString(con, text))
  expected <- data.frame(
    record_index = rep(0:5, each = 2), sample_index = rep(0:1, 6),
    alleles = c("[0, 1, 2, NULL]", "[NULL, 1]", "[NULL, NULL]", "[0, 0]",
                "[NULL, NULL]", "[0, 0]", NA, NA, NA, NA, "[200, 1]", "[1]"),
    phase_before = c("[false, true, false, true]", "[false, false]", "[true, true]", "[false, false]",
                     "[true, true]", "[false, false]", NA, NA, NA, NA, "[true, true]", "[true]"),
    phase_set = c("10", NA, NA, NA, NA, NA, "30", NA, NA, NA, "2147483647", "-1")
  )
  flatten <- function(source) dbGetQuery(con, paste0(
    "SELECT record_index::INTEGER AS record_index, c.sample_index::INTEGER AS sample_index, ",
    "c.alleles::VARCHAR AS alleles, c.phase_before::VARCHAR AS phase_before, ",
    "c.phase_set::VARCHAR AS phase_set FROM (SELECT record_index, unnest(calls) AS c FROM ",
    source, ") ORDER BY record_index, sample_index"))
  for (extension in c("vcf", "bcf", "vcf.gz")) {
    path <- fixture(paste0("geno_calls.", extension))
    expect_equal(rduckhts_bcf_samples(con, path),
                 data.frame(sample_index = 0:1, sample_name = c("S1", "S2")))
    expect_equal(rduckhts_bcf_samples(con, path, "S2,S1,S2"), rduckhts_bcf_samples(con, path))
    expect_equal(rduckhts_bcf_samples(con, path, "^S1"),
                 data.frame(sample_index = 1L, sample_name = "S2"))
    expect_equal(nrow(rduckhts_bcf_samples(con, path, "")), 0L)
    expect_equal(nrow(rduckhts_bcf_samples(con, path, "^S1,S2")), 0L)
    expect_true(rduckhts_geno(con, "geno_rows", path, overwrite = TRUE, decode_error_policy = "error"))
    expect_equal(flatten("geno_rows"), expected)
    expect_equal(nrow(rduckhts_geno(con, path = path)), 6L)
    for (selector in c("", "^S1,S2")) {
      expect_true(rduckhts_geno(con, "empty_calls", path, samples = selector, overwrite = TRUE))
      expect_equal(dbGetQuery(con, "SELECT count(*) AS n, sum(len(calls))::INTEGER AS calls FROM empty_calls"),
                   data.frame(n = 6, calls = 0L))
    }
    expect_true(rduckhts_geno(con, "selected", path, samples = "S2", overwrite = TRUE))
    expected_selected <- expected[expected$sample_index == 1L, ]
    rownames(expected_selected) <- NULL
    expect_equal(flatten("selected"), expected_selected)
    expect_true(rduckhts_geno(con, "sparse", path, non_reference_only = TRUE, overwrite = TRUE))
    expect_equal(dbGetQuery(con, "SELECT count(*) AS n, sum(len(calls))::INTEGER AS calls FROM sparse"),
                 data.frame(n = 6, calls = 4L))
    expect_equal(flatten("sparse"), expected[c(1L, 2L, 11L, 12L), ], check.attributes = FALSE)
    expect_true(rduckhts_geno(con, "sequential", path, scan_mode = "sequential", overwrite = TRUE))
    expect_equal(flatten("sequential"), expected)
    expect_silent(rduckhts_bcf(con, "bcf_selected", path, samples = "S2", overwrite = TRUE))
    actual <- dbGetQuery(con, "SELECT POS, FORMAT_GT_S2 FROM bcf_selected ORDER BY ALL")
    control <- dbGetQuery(con, paste0("SELECT POS, FORMAT_GT_S2 FROM read_bcf(", quote(path), ") ORDER BY ALL"))
    expect_equal(actual, control)
    if (extension != "vcf") {
      expect_true(rduckhts_geno(con, "regions", path, region = "chrG:20-30,chrG:20-20",
                                samples = "S2", decompression_threads = 2L, overwrite = TRUE))
      expect_equal(dbGetQuery(con, "SELECT record_index::INTEGER AS i, ID FROM regions ORDER BY i"),
                   data.frame(i = 0:2, ID = c("missing", "missing", "absent_gt")))
      expect_true(rduckhts_geno(con, "unknown", path, region = "unknown:1-20", overwrite = TRUE))
      expect_equal(dbGetQuery(con, "SELECT count(*) AS n FROM unknown")$n, 0)
    }
  }

  # Indexing sample cells and allele slots are independent: one record may
  # contain more samples than a DuckDB output chunk contains record rows.
  rduckhts_geno(con, "many_samples", fixture("tidy_chunk_boundary.vcf"))
  expect_equal(dbGetQuery(con, paste(
    "SELECT len(calls)::INTEGER AS n, calls[2053].sample_index AS last_sample,",
    "calls[2053].alleles::VARCHAR AS gt FROM many_samples")),
    data.frame(n = 2053L, last_sample = 2052L, gt = "[0, 1]"))
  rduckhts_geno(con, "many_records", fixture("bcf_filter_list_regression.vcf"))
  expect_equal(dbGetQuery(con, paste(
    "SELECT count(*) AS n, min(record_index)::INTEGER AS first_record,",
    "max(record_index)::INTEGER AS last_record, sum(len(calls))::INTEGER AS calls FROM many_records")),
    data.frame(n = 5000, first_record = 0L, last_record = 4999L, calls = 5000L))

  for (extension in c("vcf", "bcf")) {
    for (kind in c("ps_type", "ps_number", "ps_width", "gt_allele")) {
      path <- fixture(paste0("geno_", kind, ".", extension))
      expect_error(rduckhts_geno(con, path = path, decode_error_policy = "error"), pattern = "FORMAT/")
      rduckhts_geno(con, "bad_field", path, decode_error_policy = "null", overwrite = TRUE)
      if (kind == "gt_allele") {
        expect_equal(dbGetQuery(con, "SELECT calls[1].alleles IS NULL AS missing, calls[1].phase_set::INTEGER AS ps FROM bad_field"),
                     data.frame(missing = TRUE, ps = 10L))
      } else {
        expect_equal(dbGetQuery(con, "SELECT calls[1].alleles::VARCHAR AS gt, calls[1].phase_set IS NULL AS missing FROM bad_field"),
                     data.frame(gt = "[0, 1]", missing = TRUE))
      }
    }
  }
  path <- fixture("geno_calls.bcf")
  payload <- fixture("geno_ps_payload.bcf")
  expect_error(rduckhts_geno(con, path = payload, decode_error_policy = "error"), pattern = "encoded BCF type CHAR")
  for (policy in c("null", "warn")) {
    rduckhts_geno(con, "bad_payload", payload, decode_error_policy = policy, overwrite = TRUE)
    expect_equal(dbGetQuery(con, "SELECT calls[1].alleles::VARCHAR AS gt, calls[1].phase_set IS NULL AS missing FROM bad_payload"),
                 data.frame(gt = "[0, 1]", missing = TRUE))
  }
  phase_path <- fixture("geno_vcf44.vcf")
  explicit <- dbGetQuery(con, paste0(
    "SELECT calls[1].phase_before::VARCHAR AS first, calls[2].phase_before::VARCHAR AS second ",
    "FROM read_geno(", quote(phase_path), ") ORDER BY record_index"))
  expect_equal(explicit, data.frame(
    first = c("[false, true]", "[false]", "[true, false]", "[true, false, true]"),
    second = c("[true, false]", "[true]", "[false, true]", "[false, true, true]")))
  formatted <- dbGetQuery(con, paste0("SELECT FORMAT_GT_S1, FORMAT_GT_S2 FROM read_bcf(", quote(phase_path), ") ORDER BY POS"))
  expect_equal(formatted, data.frame(FORMAT_GT_S1 = c("/0|1", "/0", "|./1", "|0/1|."),
                                     FORMAT_GT_S2 = c("|0/1", "1", "/.|1", "/1|0|1")))
  # Per-file sample selection reuses the same public reader; it is not a new
  # native multi-file parser or a sample-name rebasing layer.
  rduckhts_bcf_multi(con, "selected_files", path, samples = "S2")
  expect_true("FORMAT_GT_S2" %in% dbListFields(con, "selected_files"))
  expect_false("FORMAT_GT_S1" %in% dbListFields(con, "selected_files"))
  expect_error(rduckhts_geno(con, path = path, samples = "unknown"), pattern = "not in the header")
  expect_error(rduckhts_bcf_samples(con, path, "unknown"), pattern = "not in the header")
  expect_error(rduckhts_geno(con, path = path, samples = c("S1", "S2")), pattern = "one non-missing")
  expect_error(rduckhts_geno(con, path = path, decompression_threads = -1), pattern = "whole number")
  expect_error(rduckhts_geno(con, path = path, non_reference_only = NA), pattern = "TRUE or FALSE")
  expect_error(rduckhts_geno(con, path = path, scan_mode = "random"), pattern = "scan_mode")
  expect_error(rduckhts_geno(con, path = path, region = "chrG", scan_mode = "sequential"), pattern = "incompatible")
  expect_error(rduckhts_geno(con, path = fixture("malformed_bad_pos.vcf")), pattern = "failed to read")
  # Failed replacement must leave the existing table usable.
  expect_error(rduckhts_geno(con, "geno_rows", path, samples = "unknown", overwrite = TRUE), pattern = "not in the header")
  expect_equal(flatten("geno_rows"), expected)
  expect_equal(dbGetQuery(con, "SELECT 42 AS n")$n, 42L)
}

test_record_major_genotypes()
