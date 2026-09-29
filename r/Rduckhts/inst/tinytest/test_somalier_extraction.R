library(tinytest)
library(DBI)

test_somalier_sites_import <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  source <- system.file("extdata", "somalier_sites.vcf", package = "Rduckhts")
  expect_true(nzchar(source) && file.exists(source))

  sites <- rduckhts_somalier_import_sites(con, source, "GRCh38")
  expect_equal(sites$site_index, 0:2)
  expect_equal(sites$region, c("chr1", "chr1", "chr2"))
  expect_equal(sites$allele_a, c("A", "A", "C"))
  expect_equal(sites$allele_b, c("G", "G", "T"))
  expect_equal(round(sites$population_b_af, 2), c(0.20, 0.25, 0.70))
  expect_equal(sites$source_ref, c("G", "A", "T"))
  expect_equal(sites$source_alt, c("A", "G", "C"))

  expect_true(rduckhts_somalier_import_sites(
    con, source, "GRCh38", table_name = "imported_sites"
  ))
  identity <- dbGetQuery(con, paste(
    "SELECT length(duckhts_somalier_panel_sha256('imported_sites')) panel,",
    "length(duckhts_somalier_frequency_sha256(",
    "'imported_sites', 'imported_sites')) frequency"
  ))
  expect_equal(identity$panel, 64)
  expect_equal(identity$frequency, 64)
  expect_error(
    rduckhts_somalier_import_sites(con, source, "GRCh38", max_sites = 2),
    pattern = "exceeds max_sites"
  )
  expect_error(
    rduckhts_somalier_import_sites(con, source, character()),
    pattern = "assembly must be one nonempty string"
  )
}

test_somalier_vcf_count_extraction <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  source <- system.file(
    "extdata", "mapping_number_families.vcf", package = "Rduckhts"
  )
  expect_true(nzchar(source) && file.exists(source))
  dbExecute(con, paste(
    "CREATE TEMP TABLE extraction_panel AS SELECT * FROM (VALUES",
    "('GRCh38', 0::UBIGINT, 'chr1', 100::UBIGINT, 'A', 'C'),",
    "('GRCh38', 1::UBIGINT, 'chr1', 200::UBIGINT, 'A', 'G'),",
    "('GRCh38', 2::UBIGINT, 'chr1', 300::UBIGINT, 'C', 'T'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))

  expect_true(rduckhts_somalier_vcf_counts(
    con, source, panel_table = "extraction_panel",
    samples = "S2", table_name = "extracted_counts"
  ))
  observed <- dbGetQuery(con, paste(
    "SELECT sample_id, site_index::INTEGER AS site_index,",
    "a::INTEGER AS a, b::INTEGER AS b, other::INTEGER AS other, status",
    "FROM extracted_counts ORDER BY site_index"
  ))
  expected <- data.frame(
    sample_id = rep("S2", 3), site_index = 0:2,
    a = c(0L, 0L, NA), b = c(20L, 10L, NA), other = c(0L, 5L, NA),
    status = c("measured", "measured", "unavailable_no_record")
  )
  expect_equal(observed, expected)

  bcf <- system.file("extdata", "geno_format.bcf", package = "Rduckhts")
  expect_true(nzchar(bcf) && file.exists(bcf))
  dbExecute(con, paste(
    "CREATE TEMP TABLE bcf_extraction_panel AS SELECT * FROM (VALUES",
    "('GRCh38', 0::UBIGINT, 'chrG', 20::UBIGINT, 'C', 'T'),",
    "('GRCh38', 1::UBIGINT, 'chrG', 60::UBIGINT, 'A', 'G'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))
  bcf_counts <- rduckhts_somalier_vcf_counts(
    con, bcf, panel_table = "bcf_extraction_panel", samples = "S2"
  )
  bcf_counts <- bcf_counts[order(bcf_counts$site_index), , drop = FALSE]
  rownames(bcf_counts) <- NULL
  expect_equal(bcf_counts$site_index, 0:1)
  expect_equal(bcf_counts$a, c(4, NA))
  expect_equal(bcf_counts$b, c(5, NA))
  expect_equal(bcf_counts$status, c("measured", "unavailable_no_record"))

  expect_error(
    rduckhts_somalier_vcf_counts(con, source, panel_table = "extraction_panel",
                                 filter_policy = "drop"),
    pattern = "should be one of"
  )
  expect_error(
    rduckhts_somalier_vcf_counts(con, character(), panel_table = "extraction_panel"),
    pattern = "path must be one nonempty string"
  )
}

test_somalier_vcf_count_contract_edges <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  source <- system.file(
    "extdata", "mapping_number_families.vcf", package = "Rduckhts"
  )
  expect_true(nzchar(source) && file.exists(source))
  lines <- readLines(source, warn = FALSE)
  temporary_sources <- c(
    tempfile(fileext = ".vcf"), tempfile(fileext = ".vcf"),
    tempfile(fileext = ".vcf"), tempfile(fileext = ".vcf"),
    tempfile(fileext = ".vcf"), tempfile(fileext = ".vcf")
  )
  on.exit(unlink(temporary_sources, force = TRUE), add = TRUE)

  dbExecute(con, paste(
    "CREATE TEMP TABLE count_edge_panel AS SELECT * FROM (VALUES",
    "('GRCh38', 0::UBIGINT, 'chr1', 100::UBIGINT, 'A', 'C'),",
    "('GRCh38', 1::UBIGINT, 'chr1', 200::UBIGINT, 'A', 'G'),",
    "('GRCh38', 2::UBIGINT, 'chr1', 300::UBIGINT, 'C', 'T'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))
  dbExecute(con, paste(
    "CREATE TEMP TABLE count_edge_offpanel AS",
    "SELECT * FROM count_edge_panel WHERE position = 100"
  ))

  reordered <- rduckhts_somalier_vcf_counts(
    con, source, panel_table = "count_edge_panel", samples = "S2,S1"
  )
  reordered <- reordered[order(
    as.integer(reordered$source_sample_index),
    as.integer(reordered$site_index)
  ), , drop = FALSE]
  rownames(reordered) <- NULL
  expect_equal(reordered[c(
    "sample_id", "source_sample_index", "site_index", "a", "b", "other",
    "status"
  )], data.frame(
    sample_id = c("S1", "S1", "S1", "S2", "S2", "S2"),
    source_sample_index = c(0L, 0L, 0L, 1L, 1L, 1L),
    site_index = c(0L, 1L, 2L, 0L, 1L, 2L),
    a = c(9L, 5L, NA, 0L, 0L, NA),
    b = c(3L, 9L, NA, 20L, 10L, NA),
    other = c(0L, 4L, NA, 0L, 5L, NA),
    status = c(
      "measured", "measured", "unavailable_no_record",
      "measured", "measured", "unavailable_no_record"
    )
  ))

  no_samples <- rduckhts_somalier_vcf_counts(
    con, source, panel_table = "count_edge_panel", samples = ""
  )
  expect_equal(nrow(no_samples), 0L)
  expect_true(all(c("sample_id", "site_index", "status") %in% names(no_samples)))

  empty_source <- temporary_sources[[1L]]
  writeLines(lines[startsWith(lines, "#")], empty_source, useBytes = TRUE)
  empty_counts <- rduckhts_somalier_vcf_counts(
    con, empty_source, panel_table = "count_edge_panel"
  )
  empty_counts <- empty_counts[order(
    as.integer(empty_counts$source_sample_index),
    as.integer(empty_counts$site_index)
  ), , drop = FALSE]
  rownames(empty_counts) <- NULL
  expect_equal(empty_counts[c("sample_id", "site_index", "status")], data.frame(
    sample_id = rep(c("S1", "S2"), each = 3L),
    site_index = rep(0:2, 2L),
    status = rep("unavailable_no_record", 6L)
  ))

  invalid_gt_lines <- lines
  invalid_gt_row <- grep("^chr1\t200\t", invalid_gt_lines)
  expect_equal(length(invalid_gt_row), 1L)
  invalid_gt_lines[[invalid_gt_row]] <- sub(
    "1/2:18:9,5,4:", "0/3:18:9,5,4:", invalid_gt_lines[[invalid_gt_row]],
    fixed = TRUE
  )
  invalid_gt <- temporary_sources[[2L]]
  writeLines(invalid_gt_lines, invalid_gt, useBytes = TRUE)
  expect_error(
    rduckhts_somalier_vcf_counts(
      con, invalid_gt, panel_table = "count_edge_offpanel", samples = "S1"
    ),
    pattern = "invalid FORMAT/GT allele"
  )
  unselected_invalid_gt <- rduckhts_somalier_vcf_counts(
    con, invalid_gt, panel_table = "count_edge_offpanel", samples = "S2"
  )
  expect_equal(unselected_invalid_gt[c(
    "sample_id", "source_sample_index", "a", "b", "status"
  )], data.frame(
    sample_id = "S2", source_sample_index = 1L, a = 0L, b = 20L,
    status = "measured"
  ))

  irregular_ad_lines <- lines
  multiallelic_row <- grep("^chr1\t200\t", irregular_ad_lines)
  irregular_ad_lines[[multiallelic_row]] <- sub(
    "1/2:18:9,5,4:", "1/2:18:9,5:",
    irregular_ad_lines[[multiallelic_row]], fixed = TRUE
  )
  irregular_ad_lines[[multiallelic_row]] <- sub(
    "0/2:15:10,0,5:", "0/2:15:10,-1,5,2:",
    irregular_ad_lines[[multiallelic_row]], fixed = TRUE
  )
  irregular_ad <- temporary_sources[[3L]]
  writeLines(irregular_ad_lines, irregular_ad, useBytes = TRUE)
  offpanel_ad_counts <- rduckhts_somalier_vcf_counts(
    con, irregular_ad, panel_table = "count_edge_offpanel", samples = "S1,S2"
  )
  offpanel_ad_counts <- offpanel_ad_counts[order(
    as.integer(offpanel_ad_counts$source_sample_index)
  ), , drop = FALSE]
  rownames(offpanel_ad_counts) <- NULL
  expect_equal(offpanel_ad_counts[c(
    "sample_id", "source_sample_index", "a", "b", "status"
  )], data.frame(
    sample_id = c("S1", "S2"), source_sample_index = 0:1,
    a = c(9L, 0L), b = c(3L, 20L), status = c("measured", "measured")
  ))

  filter_lines <- lines
  pos_100 <- grep("^chr1\t100\t", filter_lines)
  pos_200 <- grep("^chr1\t200\t", filter_lines)
  filter_lines[[pos_100]] <- sub("\tPASS\t", "\t.\t", filter_lines[[pos_100]])
  filter_lines[[pos_200]] <- sub("\tPASS\t", "\tq10\t", filter_lines[[pos_200]])
  chrom_header <- match(TRUE, startsWith(filter_lines, "#CHROM"))
  filter_lines <- append(
    filter_lines, '##FILTER=<ID=q10,Description="Low quality">',
    after = chrom_header - 1L
  )
  filter_source <- temporary_sources[[4L]]
  writeLines(filter_lines, filter_source, useBytes = TRUE)

  default_filter <- rduckhts_somalier_vcf_counts(
    con, filter_source, panel_table = "count_edge_panel", samples = "S1",
    filter_policy = "pass_or_unapplied"
  )
  default_filter <- default_filter[order(as.integer(default_filter$site_index)), ]
  expect_equal(default_filter$status, c(
    "measured", "unavailable_filtered", "unavailable_no_record"
  ))
  include_filtered <- rduckhts_somalier_vcf_counts(
    con, filter_source, panel_table = "count_edge_panel", samples = "S1",
    filter_policy = "include_all"
  )
  include_filtered <- include_filtered[order(as.integer(include_filtered$site_index)), ]
  expect_equal(include_filtered[c("a", "b", "other", "status")], data.frame(
    a = c(9L, 5L, NA), b = c(3L, 9L, NA), other = c(0L, 4L, NA),
    status = c("measured", "measured", "unavailable_no_record")
  ))
  expect_error(
    rduckhts_somalier_vcf_counts(
      con, filter_source, panel_table = "count_edge_panel", samples = "S1",
      filter_policy = "error"
    ),
    pattern = "panel record has a failing FILTER"
  )

  chrom_header <- match(TRUE, startsWith(lines, "#CHROM"))
  additions <- c(
    paste(
      "chr1", "300", "first", "C", "T", ".", "PASS", ".",
      "GT:DP:AD:PL:FT:LAA", "0/1:10:4,6:0,10,20:PASS:1",
      "1/1:20:0,10:50,5,0:PASS:1", sep = "\t"
    ),
    paste(
      "chr1", "400", "between", "T", "G", ".", "PASS", ".",
      "GT:DP:AD:PL:FT:LAA", "0/1:10:4,6:0,10,20:PASS:1",
      "1/1:20:0,10:50,5,0:PASS:1", sep = "\t"
    ),
    paste(
      "chr1", "100", "bad_ad", "A", "G", ".", "PASS", ".",
      "GT:DP:AD:PL:FT:LAA", "0/1:12:9:0,10,20:PASS:1",
      "1/1:20:0,20:50,5,0:PASS:1", sep = "\t"
    ),
    paste(
      "chr1", "200", "duplicate_alt", "G", "A,A", ".", "PASS", ".",
      "GT:DP:AD:PL:FT:LAA", "1/2:18:9,5,4:0,10,20,30,40,50:PASS:1,2",
      "0/2:15:10,0,5:0,15,25,35,45,55:PASS:2", sep = "\t"
    ),
    paste(
      "chr1", "300", "failing_filter", "C", "T", ".", "q10", ".",
      "GT:DP:AD:PL:FT:LAA", "0/1:10:4,6:0,10,20:PASS:1",
      "1/1:20:0,10:50,5,0:PASS:1", sep = "\t"
    )
  )
  duplicate_lines <- c(
    lines[seq_len(chrom_header - 1L)],
    '##FILTER=<ID=q10,Description="Low quality">',
    lines[chrom_header:length(lines)], additions
  )
  duplicate_source <- temporary_sources[[5L]]
  writeLines(duplicate_lines, duplicate_source, useBytes = TRUE)
  dbExecute(con, paste(
    "CREATE TEMP TABLE duplicate_panel_100 AS SELECT * FROM (VALUES",
    "('GRCh38', 0::UBIGINT, 'chr1', 100::UBIGINT, 'A', 'C'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))
  dbExecute(con, paste(
    "CREATE TEMP TABLE duplicate_panel_200 AS SELECT * FROM (VALUES",
    "('GRCh38', 0::UBIGINT, 'chr1', 200::UBIGINT, 'A', 'G'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))
  dbExecute(con, paste(
    "CREATE TEMP TABLE duplicate_panel_300 AS SELECT * FROM (VALUES",
    "('GRCh38', 0::UBIGINT, 'chr1', 300::UBIGINT, 'C', 'T'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))

  for (panel_name in c(
    "duplicate_panel_100", "duplicate_panel_200", "duplicate_panel_300"
  )) {
    expect_error(
      rduckhts_somalier_vcf_counts(
        con, duplicate_source, panel_table = panel_name, samples = "S1"
      ),
      pattern = "multiple source records occur at a panel site"
    )
  }
  quoted_duplicate_source <- as.character(dbQuoteString(con, duplicate_source))
  duplicate_counts_sql <- paste0(
    "SELECT count(*) FROM duckhts_somalier_vcf_counts(",
    quoted_duplicate_source,
    ", 'duplicate_panel_100', samples := 'S1')"
  )
  expect_error(
    dbGetQuery(con, duplicate_counts_sql),
    pattern = "multiple source records occur at a panel site"
  )
  expect_error(
    dbGetQuery(con, sub("count\\(\\*\\)", "status", duplicate_counts_sql)),
    pattern = "multiple source records occur at a panel site"
  )

  stress_source <- temporary_sources[[6L]]
  stress_connection <- file(stress_source, open = "wt")
  on.exit(try(close(stress_connection), silent = TRUE), add = TRUE)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    '##FILTER=<ID=PASS,Description="All filters passed">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"
  ), stress_connection)
  stress_record <- paste(
    "chr1", "100", ".", "A", "C", ".", "PASS", ".", "GT:AD", "0/1:9,3",
    sep = "\t"
  )
  for (batch in seq_len(100L)) {
    writeLines(rep(stress_record, 10000L), stress_connection)
  }
  close(stress_connection)
  spill_directory <- tempfile("duckhts-vcf-counts-spill-")
  dir.create(spill_directory)
  on.exit(unlink(spill_directory, recursive = TRUE), add = TRUE)
  dbExecute(con, paste0(
    "SET temp_directory=", dbQuoteString(con, spill_directory)
  ))
  dbExecute(con, "SET memory_limit='128MB'")
  expect_error(
    rduckhts_somalier_vcf_counts(
      con, stress_source, panel_table = "count_edge_offpanel", samples = "S1"
    ),
    pattern = "multiple source records occur at a panel site"
  )
  dbExecute(con, "SET memory_limit='1GB'")
}

test_somalier_bam_count_extraction <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  extdata <- function(name) system.file("extdata", name, package = "Rduckhts")
  paths <- vapply(
    c("range.bam", "range.bam.bai", "range.cram", "range.cram.crai",
      "ce.fa", "ce.fa.fai"),
    extdata, character(1)
  )
  expect_true(all(nzchar(paths)) && all(file.exists(paths)))
  dbExecute(con, paste(
    "CREATE TABLE bam_extraction_panel AS SELECT * FROM (VALUES",
    "('WBcel235', 0::UBIGINT, 'CHROMOSOME_I', 1::UBIGINT, 'A', 'G'),",
    "('WBcel235', 1::UBIGINT, 'CHROMOSOME_I', 914::UBIGINT, 'A', 'C'),",
    "('WBcel235', 2::UBIGINT, 'CHROMOSOME_I', 2::UBIGINT, 'A', 'G'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))

  extract <- function(source, index, worker_count = 1L) {
    rduckhts_somalier_bam_counts(
      con, source, "sample-1", paths[["ce.fa"]],
      panel_table = "bam_extraction_panel", index_path = index,
      reference_index_path = paths[["ce.fa.fai"]],
      worker_count = worker_count
    )
  }
  bam <- extract(paths[["range.bam"]], paths[["range.bam.bai"]])
  bam_workers <- extract(
    paths[["range.bam"]], paths[["range.bam.bai"]], worker_count = 4L
  )
  cram <- extract(paths[["range.cram"]], paths[["range.cram.crai"]])
  panel_file <- tempfile(fileext = ".parquet")
  on.exit(unlink(panel_file), add = TRUE)
  dbExecute(con, sprintf(
    "COPY bam_extraction_panel TO %s (FORMAT PARQUET)",
    as.character(dbQuoteString(con, panel_file))
  ))
  parquet <- rduckhts_somalier_bam_counts(
    con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
    panel_parquet = panel_file, index_path = paths[["range.bam.bai"]],
    reference_index_path = paths[["ce.fa.fai"]], worker_count = 2L
  )
  inferred <- rduckhts_somalier_bam_counts(
    con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
    panel_table = "bam_extraction_panel"
  )
  canonical <- function(value) {
    value <- value[order(value$site_index), , drop = FALSE]
    rownames(value) <- NULL
    value
  }
  bam <- canonical(bam)
  bam_workers <- canonical(bam_workers)
  cram <- canonical(cram)
  parquet <- canonical(parquet)
  inferred <- canonical(inferred)
  comparison_columns <- setdiff(names(bam), "source_path")
  expect_equal(bam, bam_workers)
  expect_equal(bam[comparison_columns], cram[comparison_columns])
  expect_equal(bam, parquet)
  expect_equal(bam[comparison_columns], inferred[comparison_columns])
  expect_equal(bam$a, c(0, 1, NA))
  expect_equal(bam$b, c(0, 0, NA))
  expect_equal(bam$other, c(0, 0, NA))
  expect_equal(
    bam$status,
    c("measured", "measured", "unavailable_reference_mismatch")
  )

  brace_paths <- vapply(
    c("somalier_brace.bam", "somalier_brace.bam.bai",
      "somalier_brace.fa", "somalier_brace.fa.fai"),
    extdata, character(1)
  )
  expect_true(all(nzchar(brace_paths)) && all(file.exists(brace_paths)))
  dbExecute(con, paste(
    "CREATE TABLE somalier_brace_panel AS SELECT * FROM (VALUES",
    "('synthetic', 0::UBIGINT, 'ctg}part', 10::UBIGINT, 'A', 'G'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))
  brace <- rduckhts_somalier_bam_counts(
    con, brace_paths[["somalier_brace.bam"]], "brace-sample",
    brace_paths[["somalier_brace.fa"]], panel_table = "somalier_brace_panel",
    index_path = brace_paths[["somalier_brace.bam.bai"]],
    reference_index_path = brace_paths[["somalier_brace.fa.fai"]]
  )
  expect_equal(brace[c("region", "a", "b", "other", "status")], data.frame(
    region = "ctg}part", a = 1, b = 0, other = 0, status = "measured"
  ))

  stricter <- rduckhts_somalier_bam_counts(
    con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
    panel_table = "bam_extraction_panel",
    index_path = paths[["range.bam.bai"]],
    reference_index_path = paths[["ce.fa.fai"]], min_mapq = 24
  )
  expect_equal(stricter$a[2], 0)
  expect_error(
    rduckhts_somalier_bam_counts(
      con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
      panel_table = "bam_extraction_panel", overlap_policy = "unknown"
    ),
    pattern = "should be one of"
  )
  expect_error(
    rduckhts_somalier_bam_counts(
      con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
      panel_table = "bam_extraction_panel", min_baseq = 256
    ),
    pattern = "min_baseq"
  )
  expect_error(
    rduckhts_somalier_bam_counts(
      con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
      panel_table = "bam_extraction_panel", worker_count = 0
    ),
    pattern = "worker_count"
  )
  dbExecute(con, "CREATE TEMP TABLE bam_temp_panel AS SELECT * FROM bam_extraction_panel")
  expect_error(
    rduckhts_somalier_bam_counts(
      con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
      panel_table = "bam_temp_panel"
    ),
    pattern = "does not exist"
  )
}

test_somalier_sites_import()
test_somalier_vcf_count_extraction()
test_somalier_vcf_count_contract_edges()
test_somalier_bam_count_extraction()
