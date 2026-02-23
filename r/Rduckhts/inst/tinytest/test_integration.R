# Integration tests for Rduckhts package - shows RBCFTools-like table creation
library(tinytest)
library(DBI)

# Test table creation pattern (like RBCFTools)
test_table_creation <- function() {
  # Create an in-memory DuckDB connection
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)

  # Ensure bundled extension can be loaded
  expect_silent(rduckhts_load(con))

  bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
  bcf_index_path <- system.file("extdata", "vcf_file.bcf.csi", package = "Rduckhts")
  formatcols_vcf_path <- system.file("extdata", "formatcols.vcf.gz", package = "Rduckhts")
  formatcols_vcf_index_path <- system.file("extdata", "formatcols.vcf.gz.csi", package = "Rduckhts")
  bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
  bam_index_path <- system.file("extdata", "range.bam.bai", package = "Rduckhts")
  fasta_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
  fasta_index_path <- system.file("extdata", "ce.fa.fai", package = "Rduckhts")
  fastq_r1 <- system.file("extdata", "r1.fq", package = "Rduckhts")
  fastq_r2 <- system.file("extdata", "r2.fq", package = "Rduckhts")
  gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
  gff_index_path <- system.file("extdata", "gff_file.gff.gz.tbi", package = "Rduckhts")
  tabix_path <- system.file("extdata", "rg.sam.gz", package = "Rduckhts")
  header_tabix_path <- system.file(
    "extdata",
    "header_tabix.tsv.gz",
    package = "Rduckhts"
  )
  meta_tabix_path <- system.file(
    "extdata",
    "meta_tabix.tsv.gz",
    package = "Rduckhts"
  )
  vep_path <- system.file("extdata", "test_vep.vcf", package = "Rduckhts")

  expect_true(file.exists(bcf_path))
  expect_true(file.exists(bcf_index_path))
  expect_true(file.exists(formatcols_vcf_path))
  expect_true(file.exists(formatcols_vcf_index_path))
  expect_true(file.exists(bam_path))
  expect_true(file.exists(bam_index_path))
  expect_true(file.exists(fasta_path))
  expect_true(file.exists(fasta_index_path))
  expect_true(file.exists(fastq_r1))
  expect_true(file.exists(fastq_r2))
  expect_true(file.exists(gff_path))
  expect_true(file.exists(gff_index_path))
  expect_true(file.exists(tabix_path))
  expect_true(file.exists(header_tabix_path))
  expect_true(file.exists(meta_tabix_path))
  expect_true(file.exists(vep_path))

  expect_silent(rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE))
  expect_silent(rduckhts_bcf(
    con,
    "variants_idx",
    bcf_path,
    region = "1:3000150-3000151",
    index_path = bcf_index_path,
    overwrite = TRUE
  ))
  expect_silent(rduckhts_bam(con, "reads", bam_path, overwrite = TRUE))
  expect_silent(rduckhts_bam(
    con,
    "reads_idx",
    bam_path,
    region = "CHROMOSOME_I:1-1000",
    index_path = bam_index_path,
    overwrite = TRUE
  ))
  expect_silent(rduckhts_fasta(con, "sequences", fasta_path, overwrite = TRUE))
  expect_silent(rduckhts_fasta(
    con,
    "sequences_region",
    fasta_path,
    region = "CHROMOSOME_I:1-10",
    index_path = fasta_index_path,
    overwrite = TRUE
  ))
  expect_silent(rduckhts_fasta_index(con, fasta_path))
  expect_silent(rduckhts_fastq(
    con,
    "fastq_reads",
    fastq_r1,
    mate_path = fastq_r2,
    overwrite = TRUE
  ))
  expect_silent(rduckhts_gff(
    con,
    "annotations",
    gff_path,
    attributes_map = TRUE,
    overwrite = TRUE
  ))
  expect_silent(rduckhts_tabix(con, "tabix_data", tabix_path, overwrite = TRUE))
  expect_silent(rduckhts_tabix(
    con,
    "tabix_idx",
    gff_path,
    region = "X:2934816-2935190",
    index_path = gff_index_path,
    overwrite = TRUE
  ))
  expect_silent(rduckhts_bcf(con, "vep_variants", vep_path, overwrite = TRUE))
  expect_true(DBI::dbGetQuery(con, "SELECT count(*) AS n FROM variants_idx")$n[1] == 2)
  expect_true(DBI::dbGetQuery(con, "SELECT count(*) AS n FROM reads_idx")$n[1] == 2)
  expect_true(DBI::dbGetQuery(con, "SELECT count(*) AS n FROM tabix_idx")$n[1] == 4)
  expect_true(DBI::dbGetQuery(con, "SELECT count(*) AS n FROM sequences_region")$n[1] == 1)
  expect_true(DBI::dbGetQuery(con, "SELECT length(SEQUENCE) AS n FROM sequences_region")$n[1] == 10)

  header_meta <- rduckhts_hts_header(con, bcf_path)
  expect_true(nrow(header_meta) > 0)
  expect_true("record_type" %in% names(header_meta))
  header_raw <- rduckhts_hts_header(con, bcf_path, mode = "raw")
  expect_true(nrow(header_raw) > 0)
  expect_true("raw" %in% names(header_raw))

  index_meta <- rduckhts_hts_index(con, bcf_path, index_path = bcf_index_path)
  expect_true(nrow(index_meta) > 0)
  expect_true("index_type" %in% names(index_meta))
  index_spans <- rduckhts_hts_index_spans(con, bcf_path, index_path = bcf_index_path)
  expect_true(nrow(index_spans) > 0)
  expect_true("chunk_beg_vo" %in% names(index_spans))
  index_raw <- rduckhts_hts_index_raw(
    con,
    formatcols_vcf_path,
    index_path = formatcols_vcf_index_path
  )
  expect_true(nrow(index_raw) == 1)
  expect_true("raw" %in% names(index_raw))
  expect_true(DBI::dbExistsTable(con, "annotations"))
  if (DBI::dbExistsTable(con, "annotations")) {
    expect_silent(DBI::dbGetQuery(con, "SELECT * FROM annotations LIMIT 1"))

    # Check MAP types specifically
    schema <- DBI::dbGetQuery(con, "DESCRIBE annotations")
    message("GFF schema with attributes_map=TRUE:\n", paste(capture.output(schema), collapse = "\n"))

    # Check what type the attributes column has
    attr_type <- DBI::dbGetQuery(
      con,
      "SELECT typeof(attributes) as attr_type FROM annotations LIMIT 1"
    )$attr_type[1]
    message("Attributes column type: ", attr_type)

    # Get sample to see MAP structure
    sample <- DBI::dbGetQuery(con, "SELECT attributes FROM annotations LIMIT 1")
    message("Sample MAP content:\n", paste(capture.output(sample), collapse = "\n"))
  }

  expect_true(DBI::dbExistsTable(con, "reads"))
  if (DBI::dbExistsTable(con, "reads")) {
    expect_silent(DBI::dbGetQuery(con, "SELECT * FROM reads LIMIT 1"))
  }

  expect_true(DBI::dbExistsTable(con, "sequences"))
  if (DBI::dbExistsTable(con, "sequences")) {
    expect_silent(DBI::dbGetQuery(con, "SELECT * FROM sequences LIMIT 1"))
  }

  expect_true(DBI::dbExistsTable(con, "fastq_reads"))
  if (DBI::dbExistsTable(con, "fastq_reads")) {
    expect_silent(DBI::dbGetQuery(con, "SELECT * FROM fastq_reads LIMIT 1"))
  }

  expect_true(DBI::dbExistsTable(con, "annotations"))
  if (DBI::dbExistsTable(con, "annotations")) {
    expect_silent(DBI::dbGetQuery(con, "SELECT * FROM annotations LIMIT 1"))
  }

  expect_true(DBI::dbExistsTable(con, "tabix_data"))
  if (DBI::dbExistsTable(con, "tabix_data")) {
    expect_silent(DBI::dbGetQuery(con, "SELECT * FROM tabix_data LIMIT 1"))
  }

  expect_silent(rduckhts_tabix(
    con,
    "tabix_header",
    header_tabix_path,
    header = TRUE,
    overwrite = TRUE
  ))
  tabix_header_cols <- DBI::dbGetQuery(con, "PRAGMA table_info('tabix_header')")
  expect_equal(tabix_header_cols$name, c("chrom", "pos", "value"))

  expect_silent(rduckhts_tabix(
    con,
    "tabix_named",
    meta_tabix_path,
    header_names = c("chr", "pos", "val"),
    overwrite = TRUE
  ))
  tabix_named_cols <- DBI::dbGetQuery(con, "PRAGMA table_info('tabix_named')")
  expect_equal(tabix_named_cols$name, c("chr", "pos", "val"))

  expect_silent(rduckhts_tabix(
    con,
    "tabix_auto",
    meta_tabix_path,
    auto_detect = TRUE,
    overwrite = TRUE
  ))
  tabix_auto_type <- DBI::dbGetQuery(
    con,
    "SELECT typeof(column1) AS t FROM tabix_auto LIMIT 1"
  )$t[1]
  expect_equal(tabix_auto_type, "BIGINT")

  expect_silent(rduckhts_tabix(
    con,
    "tabix_types",
    meta_tabix_path,
    column_types = c("VARCHAR", "BIGINT", "VARCHAR"),
    overwrite = TRUE
  ))
  tabix_types_type <- DBI::dbGetQuery(
    con,
    "SELECT typeof(column1) AS t FROM tabix_types LIMIT 1"
  )$t[1]
  expect_equal(tabix_types_type, "BIGINT")

  expect_true(DBI::dbExistsTable(con, "vep_variants"))
  if (DBI::dbExistsTable(con, "vep_variants")) {
    expect_silent(DBI::dbGetQuery(
      con,
      "SELECT VEP_Allele FROM vep_variants LIMIT 1"
    ))
  }

  # Test overwrite parameter validation
  if (DBI::dbExistsTable(con, "test_table")) {
    DBI::dbRemoveTable(con, "test_table")
  }

  expect_silent(rduckhts_bcf(con, "test_table", bcf_path, overwrite = TRUE))
  expect_error(rduckhts_bcf(con, "test_table", bcf_path, overwrite = FALSE))

  dbDisconnect(con, shutdown = TRUE)
  message("Table creation pattern tests passed!")
}

# Run the test
test_table_creation()

message("Integration tests completed!")
