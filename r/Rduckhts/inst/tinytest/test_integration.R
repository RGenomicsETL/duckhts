# Integration tests for Rduckhts package - shows RBCFTools-like table creation
library(tinytest)
library(DBI)

# Skip if duckdb not available
if (!requireNamespace("duckdb", quietly = TRUE)) {
  cat("Skipping integration tests - duckdb not available\n")
  quit(status = 0)
}

# Test table creation pattern (like RBCFTools)
test_table_creation <- function() {
  # Create an in-memory DuckDB connection
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)

  # Ensure bundled extension exists and can be loaded
  ext_path <- system.file("extdata", "duckhts.duckdb_extension", package = "Rduckhts")
  expect_true(file.exists(ext_path))
  expect_silent(rduckhts_load(con, ext_path))

  # Test table creation functions without extension (should fail gracefully)
  expect_error(rduckhts_bcf(con, "variants", "nonexistent.vcf"))
  expect_error(rduckhts_bam(con, "reads", "nonexistent.bam"))
  expect_error(rduckhts_fasta(con, "sequences", "nonexistent.fa"))

  # Test overwrite parameter validation
  if (dbExistsTable(con, "test_table")) {
    dbRemoveTable(con, "test_table")
  }

  expect_error(rduckhts_bcf(
    con,
    "test_table",
    "nonexistent.vcf",
    overwrite = FALSE
  ))
  expect_error(rduckhts_bcf(
    con,
    "test_table",
    "nonexistent.vcf",
    overwrite = TRUE
  ))

  dbDisconnect(con, shutdown = TRUE)
  cat("Table creation pattern tests passed!\n")
}

# Run the test
test_table_creation()

cat("Integration tests completed!\n")
