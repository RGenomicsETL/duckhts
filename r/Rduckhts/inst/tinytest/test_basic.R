# Basic functionality tests for Rduckhts package
library(tinytest)

# Test package loading
expect_true(requireNamespace("Rduckhts", quietly = TRUE))

# Test basic functions exist
expect_true(exists("rduckhts_load"))
expect_true(exists("rduckhts_bcf"))
expect_true(exists("rduckhts_bam"))
expect_true(exists("rduckhts_fasta"))
expect_true(exists("rduckhts_fastq"))
expect_true(exists("rduckhts_gff"))
expect_true(exists("rduckhts_gtf"))
expect_true(exists("rduckhts_tabix"))

# Test function signatures
expect_equal(length(formals(rduckhts_load)), 2)
expect_equal(length(formals(rduckhts_bcf)), 5)
expect_equal(length(formals(rduckhts_bam)), 5)
expect_equal(length(formals(rduckhts_fasta)), 3)
expect_equal(length(formals(rduckhts_fastq)), 5)
expect_equal(length(formals(rduckhts_gff)), 5)
expect_equal(length(formals(rduckhts_gtf)), 5)
expect_equal(length(formals(rduckhts_tabix)), 4)

# Test that DBI is available
expect_true(requireNamespace("DBI", quietly = TRUE))

# Test parameter validation - these should fail gracefully without a connection
expect_error(rduckhts_bcf(NULL, "test", "nonexistent.vcf"))
expect_error(rduckhts_bam(NULL, "test", "nonexistent.bam"))
expect_error(rduckhts_fasta(NULL, "test", "nonexistent.fa"))
expect_error(rduckhts_fastq(NULL, "test", "nonexistent.fq"))
expect_error(rduckhts_gff(NULL, "test", "nonexistent.gff"))
expect_error(rduckhts_gtf(NULL, "test", "nonexistent.gtf"))
expect_error(rduckhts_tabix(NULL, "test", "nonexistent.bed.gz"))

cat("All basic tests passed!\n")
