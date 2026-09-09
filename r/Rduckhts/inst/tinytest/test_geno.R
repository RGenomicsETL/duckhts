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

  # Missing fields and changing ploidy recur across several output chunks.
  repeated <- tempfile("geno-chunks-", fileext = ".vcf")
  on.exit(unlink(repeated), add = TRUE)
  seed <- readLines(fixture("geno_calls.vcf"))
  writeLines(c(seed[startsWith(seed, "#")], rep(seed[!startsWith(seed, "#")], 1000L)), repeated)
  rduckhts_geno(con, "repeated_calls", repeated, decode_error_policy = "error")
  expected_repeated <- expected[rep(seq_len(nrow(expected)), 1000L), ]
  expected_repeated$record_index <- rep(0:5999, each = 2L)
  rownames(expected_repeated) <- NULL
  expect_equal(flatten("repeated_calls"), expected_repeated)
  # One record can contain more sample cells than an output chunk has rows.
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
    for (policy in c("null", "warn", "error")) {
      rduckhts_geno(con, "selected_ps", fixture(paste0("geno_ps_width.", extension)),
                    samples = "S2", decode_error_policy = policy, overwrite = TRUE)
      expect_equal(dbGetQuery(con, paste(
        "SELECT calls[1].sample_index AS sample, calls[1].phase_set::INTEGER AS ps FROM selected_ps")),
        data.frame(sample = 1, ps = 20L))
    }
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

test_selected_genotype_format <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  fields <- c("AD", "DP", "GQ", "GL", "VI", "VF", "ST")
  fixture <- function(name) {
    path <- system.file("extdata", name, package = "Rduckhts")
    stopifnot(nzchar(path), file.exists(path))
    path
  }
  flatten <- function(source) paste0(
    "SELECT CHROM, POS, ID, REF, ALT, s.sample_name AS sample, ",
    paste0("c.format.", fields, " AS ", fields, collapse = ", "),
    " FROM ", source, " g, unnest(g.calls) u(c), samples s WHERE c.sample_index=s.sample_index")
  expected <- function(source) paste0(
    "SELECT CHROM, POS, ID, REF, ALT, SAMPLE_ID AS sample, ",
    paste0("FORMAT_", fields, " AS ", fields, collapse = ", "), " FROM ", source)
  compare <- function(left, right) {
    # Bag comparison counts every physical duplicate; NULL items stay in place.
    sql <- sprintf(paste(
      "SELECT (SELECT count(*) FROM ((%s) EXCEPT ALL (%s))) AS extra,",
      "(SELECT count(*) FROM ((%s) EXCEPT ALL (%s))) AS missing"), left, right, right, left)
    expect_equal(dbGetQuery(con, sql), data.frame(extra = 0, missing = 0))
  }
  for (extension in c("vcf", "bcf", "vcf.gz")) {
    path <- fixture(paste0("geno_format.", extension))
    quoted <- dbQuoteString(con, path)
    dbWriteTable(con, "samples", rduckhts_bcf_samples(con, path), overwrite = TRUE)
    for (selection in c("-", "S2", "^S2", "")) {
      rduckhts_geno(con, "selected", path, samples = selection, format_fields = fields,
                    decode_error_policy = "error", overwrite = TRUE)
      baseline <- sprintf("read_bcf(%s, tidy_format := true, samples := %s, scan_mode := 'sequential')",
                          quoted, dbQuoteString(con, selection))
      expect_equal(dbGetQuery(con, "SELECT count(*) AS n FROM selected")$n, 6)
      if (nzchar(selection)) compare(flatten("selected"), expected(baseline))
      else expect_equal(dbGetQuery(con, "SELECT sum(len(calls)) AS n FROM selected")$n, 0)
    }
    rduckhts_geno(con, "selected", path, format_fields = fields, overwrite = TRUE)
    first <- dbGetQuery(con, paste(
      "SELECT calls[1].format.AD::VARCHAR AS ad, calls[1].format.GL::VARCHAR AS gl,",
      "calls[1].format.ST::VARCHAR AS st, calls[1].alleles::VARCHAR AS gt",
      "FROM selected WHERE record_index=0"))
    expect_equal(first, data.frame(ad = "[10, NULL, 5]",
      gl = "[0.0, -1.0, NULL, -3.0, -4.0, -5.0]", st = "[a, NULL, b]", gt = "[NULL, NULL]"))
    expect_equal(dbGetQuery(con, paste(
      "SELECT calls[1].alleles IS NULL AS missing_gt, calls[1].format.AD::VARCHAR AS ad",
      "FROM selected WHERE record_index=2")), data.frame(missing_gt = TRUE, ad = "[8, 9]"))
    expect_equal(dbGetQuery(con, sprintf(paste(
      "SELECT INFO_MI::VARCHAR AS i, INFO_MF::VARCHAR AS f, INFO_MS::VARCHAR AS s",
      "FROM read_bcf(%s) WHERE POS=10"), quoted)),
      data.frame(i = "[1, NULL, 3]", f = "[1.5, NULL, 2.5]", s = "[a, NULL, b]"))
    if (extension != "vcf") {
      region <- "chrG:20-30,chrG:30-30"
      rduckhts_geno(con, "selected", path, format_fields = fields, region = region, overwrite = TRUE)
      compare(flatten("selected"), expected(sprintf("read_bcf(%s, tidy_format := true, region := %s)",
                                                   quoted, dbQuoteString(con, region))))
      expect_equal(dbGetQuery(con, "SELECT count(*) AS n FROM selected")$n, 3)
    }
    for (empty in list(NULL, character())) {
      default <- rduckhts_geno(con, path = path)
      expect_equal(rduckhts_geno(con, path = path, format_fields = empty), default)
    }
  }
  # Repeated growing/missing fields exercise worker-cache reset across output chunks.
  repeated <- tempfile("geno-format-chunks-", fileext = ".vcf")
  on.exit(unlink(repeated), add = TRUE)
  seed <- readLines(fixture("geno_format.vcf"))
  writeLines(c(seed[startsWith(seed, "#")], rep(seed[!startsWith(seed, "#")], 1000L)), repeated)
  rduckhts_geno(con, "selected", repeated, format_fields = fields, overwrite = TRUE)
  expect_equal(dbGetQuery(con, "SELECT count(*) AS n, sum(len(calls)) AS calls FROM selected"),
               data.frame(n = 6000, calls = 12000))
  compare(flatten("selected"), expected(sprintf("read_bcf(%s, tidy_format := true, scan_mode := 'sequential')",
                                               dbQuoteString(con, repeated))))
  # A single call list and the tidy counterpart span more than 2,048 samples.
  many_samples <- tempfile("geno-format-samples-", fileext = ".vcf")
  on.exit(unlink(many_samples), add = TRUE)
  row <- strsplit(seed[!startsWith(seed, "#")][1L], "\t", fixed = TRUE)[[1L]]
  header <- strsplit(seed[startsWith(seed, "#CHROM")], "\t", fixed = TRUE)[[1L]]
  writeLines(c(seed[startsWith(seed, "##")],
    paste(c(header[1:9], paste0("S", seq_len(2053L))), collapse = "\t"),
    paste(c(row[1:9], rep(row[10:11], length.out = 2053L)), collapse = "\t")), many_samples)
  dbWriteTable(con, "samples", rduckhts_bcf_samples(con, many_samples), overwrite = TRUE)
  rduckhts_geno(con, "selected", many_samples, format_fields = fields, overwrite = TRUE)
  expect_equal(dbGetQuery(con, paste(
    "SELECT len(calls)::INTEGER AS n, calls[2053].sample_index AS last_sample,",
    "calls[2053].format.AD::VARCHAR AS ad FROM selected")),
    data.frame(n = 2053L, last_sample = 2052L, ad = "[10, NULL, 5]"))
  compare(flatten("selected"), expected(sprintf("read_bcf(%s, tidy_format := true, scan_mode := 'sequential')",
                                               dbQuoteString(con, many_samples))))
  path <- fixture("geno_format.bcf")
  for (invalid in list(1, NA_character_, "", c("AD", NA_character_))) {
    expect_error(rduckhts_geno(con, path = path, format_fields = invalid), pattern = "format_fields")
  }
  expect_error(rduckhts_geno(con, path = path, format_fields = c("AD", "ad")), pattern = "duplicate")
  expect_error(rduckhts_geno(con, path = path, format_fields = "absent"), pattern = "not declared")
  expect_error(rduckhts_geno(con, path = path, format_fields = "GT"), pattern = "already exposed")
  expect_error(rduckhts_geno(con, path = path, format_fields = "PS"), pattern = "already exposed")
  clash <- fixture("bcf_format_type_clash.bcf")
  expect_error(rduckhts_geno(con, path = clash, format_fields = "XX", decode_error_policy = "error"),
               pattern = "FORMAT/XX")
  for (policy in c("null", "warn")) {
    rduckhts_geno(con, "clash", clash, format_fields = "XX", decode_error_policy = policy, overwrite = TRUE)
    expect_true(dbGetQuery(con, "SELECT calls[1].format.XX IS NULL AS missing FROM clash")$missing)
  }
  expect_equal(dbGetQuery(con, sprintf(paste(
    "SELECT count(*) AS n FROM read_geno(%s, format_fields := ['XX'], decode_error_policy := 'error')"),
    dbQuoteString(con, clash)))$n, 1)
  expect_equal(dbGetQuery(con, "SELECT 42 AS n")$n, 42L)
}

test_record_major_genotypes()
test_genotype_numeric_scalars <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  expected <- data.frame(pos = seq(10, 70, 10), si = c(8, 10, 13, NA, 20, NA, 23),
                         sf = c(2.5, 4.5, 6.5, 9.5, NA, NA, 15.5),
                         si2 = c(9, 11, 14, NA, 21, NA, 24),
                         sf2 = c(3.5, 5.5, 7.5, 10.5, NA, NA, 16.5))
  for (extension in c("vcf", "bcf", "vcf.gz")) {
    path <- system.file("extdata", paste0("bcf_scalar_counts.", extension), package = "Rduckhts")
    expect_true(nzchar(path))
    for (policy in c("null", "warn")) {
      rduckhts_geno(con, "scalar_calls", path, format_fields = c("SI", "SF"),
                    decode_error_policy = policy, overwrite = TRUE)
      expect_equal(dbGetQuery(con, paste(
        "SELECT POS AS pos, calls[1].format.SI AS si, calls[1].format.SF AS sf,",
        "calls[2].format.SI AS si2, calls[2].format.SF AS sf2 FROM scalar_calls ORDER BY POS")), expected)
    }
    for (field in c("SI", "SF")) expect_error(
      rduckhts_geno(con, path = path, format_fields = field, decode_error_policy = "error"),
      pattern = "has 2 values for header Number=1 at chrS:")
    rduckhts_geno(con, "scalar_calls", path, format_fields = c("SI", "SF"), samples = "",
                  decode_error_policy = "error", overwrite = TRUE)
    expect_equal(dbGetQuery(con, "SELECT count(*) AS n, sum(len(calls)) AS calls FROM scalar_calls"),
                  data.frame(n = 7, calls = 0))
    if (extension != "vcf") {
      rduckhts_geno(con, "scalar_calls", path, region = "chrS:40-40", format_fields = "SI",
                    samples = "S1", decode_error_policy = "error", overwrite = TRUE)
      expect_equal(dbGetQuery(con, paste(
        "SELECT calls[1].sample_index AS sample, calls[1].format.SI AS si FROM scalar_calls")),
        data.frame(sample = 0, si = 16))
    }
  }
  expect_equal(dbGetQuery(con, "SELECT 42 AS n")$n, 42L)
}

test_selected_genotype_format()
test_genotype_numeric_scalars()
