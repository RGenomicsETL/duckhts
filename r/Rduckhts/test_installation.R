#!/usr/bin/env Rscript
# Simple test script to verify Rduckhts installation

library(Rduckhts)
library(DBI)
library(duckdb)

cat("Testing Rduckhts installation...\n")

# Create connection
con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))

# Test extension loading
cat("Loading DuckHTS extension...\n")
tryCatch({
  rduckhts_load(con)
  cat("✓ DuckHTS extension loaded successfully\n")
}, error = function(e) {
  cat("✗ Failed to load DuckHTS extension:", e$message, "\n")
  cat("This indicates the extension file may not have been built during installation.\n")
  quit(status = 1)
})

# Test that functions are available
functions_to_test <- c("rduckhts_bcf", "rduckhts_bam", "rduckhts_fasta", 
                     "rduckhts_fastq", "rduckhts_gff", "rduckhts_gtf", "rduckhts_tabix")

cat("Checking function availability:\n")
for (func in functions_to_test) {
  if (exists(func)) {
    cat("✓", func, "available\n")
  } else {
    cat("✗", func, "missing\n")
  }
}

# Clean up
dbDisconnect(con, shutdown = TRUE)
cat("\nRduckhts installation test completed!\n")