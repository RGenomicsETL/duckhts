# Integration tests for Rduckhts package - shows RBCFTools-like table creation
library(tinytest)
library(DBI)

# Test table creation pattern (like RBCFTools)
test_table_creation <- function() {
  # Create an in-memory DuckDB connection
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)

  # Ensure bundled extension exists and can be loaded
  ext_path <- system.file("extdata", "duckhts.duckdb_extension", package = "Rduckhts")
  expect_true(file.exists(ext_path))
  expect_silent(rduckhts_load(con, ext_path))

  bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
  bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
  fasta_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
  fastq_r1 <- system.file("extdata", "r1.fq", package = "Rduckhts")
  fastq_r2 <- system.file("extdata", "r2.fq", package = "Rduckhts")
  gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
  tabix_path <- system.file("extdata", "rg.sam.gz", package = "Rduckhts")

  expect_true(file.exists(bcf_path))
  expect_true(file.exists(bam_path))
  expect_true(file.exists(fasta_path))
  expect_true(file.exists(fastq_r1))
  expect_true(file.exists(fastq_r2))
  expect_true(file.exists(gff_path))
  expect_true(file.exists(tabix_path))

  expect_silent(rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE))
  expect_silent(rduckhts_bam(con, "reads", bam_path, overwrite = TRUE))
  expect_silent(rduckhts_fasta(con, "sequences", fasta_path, overwrite = TRUE))
  expect_silent(rduckhts_fastq(con, "fastq_reads", fastq_r1, mate_path = fastq_r2, overwrite = TRUE))
  expect_silent(rduckhts_gff(con, "annotations", gff_path, attributes_map = TRUE, overwrite = TRUE))
  expect_silent(rduckhts_tabix(con, "tabix_data", tabix_path, overwrite = TRUE))

  expect_true(nrow(DBI::dbGetQuery(con, "SELECT * FROM variants LIMIT 1")) >= 0)
  expect_true(nrow(DBI::dbGetQuery(con, "SELECT * FROM reads LIMIT 1")) >= 0)
  expect_true(nrow(DBI::dbGetQuery(con, "SELECT * FROM sequences LIMIT 1")) >= 0)
  expect_true(nrow(DBI::dbGetQuery(con, "SELECT * FROM fastq_reads LIMIT 1")) >= 0)
  expect_true(nrow(DBI::dbGetQuery(con, "SELECT * FROM annotations LIMIT 1")) >= 0)
  expect_true(nrow(DBI::dbGetQuery(con, "SELECT * FROM tabix_data LIMIT 1")) >= 0)

  # Test overwrite parameter validation
  if (DBI::dbExistsTable(con, "test_table")) {
    DBI::dbRemoveTable(con, "test_table")
  }

  expect_silent(rduckhts_bcf(con, "test_table", bcf_path, overwrite = TRUE))
  expect_error(rduckhts_bcf(con, "test_table", bcf_path, overwrite = FALSE))

  dbDisconnect(con, shutdown = TRUE)
  cat("Table creation pattern tests passed!\n")
}

# Run the test
test_table_creation()

cat("Integration tests completed!\n")
