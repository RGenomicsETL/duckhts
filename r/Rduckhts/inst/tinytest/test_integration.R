# Integration tests for Rduckhts package - shows RBCFTools-like table creation
library(tinytest)
library(DBI)

# Test table creation pattern (like RBCFTools)
test_table_creation <- function() {
  # Create an in-memory DuckDB connection
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)

  # Ensure bundled extension exists and can be loaded
  ext_path <- system.file(
    "extdata",
    "duckhts.duckdb_extension",
    package = "Rduckhts"
  )
  expect_true(file.exists(ext_path))
  expect_silent(rduckhts_load(con, ext_path))

  bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
  bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
  fasta_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
  fastq_r1 <- system.file("extdata", "r1.fq", package = "Rduckhts")
  fastq_r2 <- system.file("extdata", "r2.fq", package = "Rduckhts")
  gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
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
  expect_true(file.exists(bam_path))
  expect_true(file.exists(fasta_path))
  expect_true(file.exists(fastq_r1))
  expect_true(file.exists(fastq_r2))
  expect_true(file.exists(gff_path))
  expect_true(file.exists(tabix_path))
  expect_true(file.exists(header_tabix_path))
  expect_true(file.exists(meta_tabix_path))
  expect_true(file.exists(vep_path))

  expect_silent(rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE))
  expect_silent(rduckhts_bam(con, "reads", bam_path, overwrite = TRUE))
  expect_silent(rduckhts_fasta(con, "sequences", fasta_path, overwrite = TRUE))
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
  expect_silent(rduckhts_bcf(con, "vep_variants", vep_path, overwrite = TRUE))
  expect_true(DBI::dbExistsTable(con, "annotations"))
  if (DBI::dbExistsTable(con, "annotations")) {
    expect_silent(DBI::dbGetQuery(con, "SELECT * FROM annotations LIMIT 1"))

    # Check MAP types specifically
    schema <- DBI::dbGetQuery(con, "DESCRIBE annotations")
    cat("GFF schema with attributes_map=TRUE:\n")
    print(schema)

    # Check what type the attributes column has
    attr_type <- DBI::dbGetQuery(
      con,
      "SELECT typeof(attributes) as attr_type FROM annotations LIMIT 1"
    )$attr_type[1]
    cat("Attributes column type:", attr_type, "\n")

    # Get sample to see MAP structure
    sample <- DBI::dbGetQuery(con, "SELECT attributes FROM annotations LIMIT 1")
    cat("Sample MAP content:\n")
    print(sample)
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
