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
expect_true(exists("rduckhts_hts_header"))
expect_true(exists("rduckhts_hts_index"))
expect_true(exists("rduckhts_hts_index_spans"))
expect_true(exists("rduckhts_hts_index_raw"))

# Test function signatures
expect_equal(length(formals(rduckhts_load)), 2)
expect_equal(length(formals(rduckhts_bcf)), 7)
expect_equal(length(formals(rduckhts_bam)), 9)
expect_equal(length(formals(rduckhts_fasta)), 4)
expect_equal(length(formals(rduckhts_fastq)), 6)
expect_equal(length(formals(rduckhts_gff)), 6)
expect_equal(length(formals(rduckhts_gtf)), 6)
expect_equal(length(formals(rduckhts_tabix)), 10)
expect_equal(length(formals(rduckhts_hts_header)), 4)
expect_equal(length(formals(rduckhts_hts_index)), 4)
expect_equal(length(formals(rduckhts_hts_index_spans)), 4)
expect_equal(length(formals(rduckhts_hts_index_raw)), 4)

# Test that DBI is available
expect_true(requireNamespace("DBI", quietly = TRUE))

# Test type mapping functions
expect_true(exists("duckdb_type_mappings"))
expect_true(exists("normalize_tabix_types"))

# Test duckdb_type_mappings function
mappings <- duckdb_type_mappings()
expect_true(is.list(mappings))
expect_true("duckdb_to_r" %in% names(mappings))
expect_true("r_to_duckdb" %in% names(mappings))
expect_true(is.character(mappings$duckdb_to_r))
expect_true(is.character(mappings$r_to_duckdb))

# Test specific mappings (corrected based on actual DuckDB behavior)
# Note: mappings return named character vectors, so we need to extract the values
expect_equal(as.character(mappings$duckdb_to_r["BIGINT"]), "double")
expect_equal(as.character(mappings$duckdb_to_r["DOUBLE"]), "double")
expect_equal(as.character(mappings$duckdb_to_r["VARCHAR"]), "character")
expect_equal(as.character(mappings$duckdb_to_r["BOOLEAN"]), "logical")

expect_equal(as.character(mappings$r_to_duckdb["integer"]), "INTEGER")
expect_equal(as.character(mappings$r_to_duckdb["numeric"]), "DOUBLE")
# Test the actual value works (named vector)
logical_mapping <- mappings$r_to_duckdb["logical"]
expect_equal(as.character(logical_mapping), "BOOLEAN")
expect_equal(names(logical_mapping), "logical")

expect_equal(as.character(mappings$r_to_duckdb["character"]), "VARCHAR")

# Test normalize_tabix_types function
test_types <- c("integer", "numeric", "character", "logical", "unknown")
expected <- c("BIGINT", "DOUBLE", "VARCHAR", "BOOLEAN", "UNKNOWN")
expect_equal(normalize_tabix_types(test_types), expected)

# Test complex type helper functions
expect_true(exists("detect_complex_types"))
expect_true(exists("extract_array_element"))
expect_true(exists("extract_map_data"))

# Test that type mappings include complex types
mappings <- duckdb_type_mappings()
expect_true("MAP" %in% names(mappings$duckdb_to_r))
expect_equal(as.character(mappings$duckdb_to_r["MAP"]), "data.frame")
expect_true(any(grepl("\\[", names(mappings$duckdb_to_r)))) # Array types

# Test example files are bundled
expect_true(file.exists(system.file("extdata", "ce.fa", package = "Rduckhts")))
expect_true(file.exists(system.file("extdata", "r1.fq", package = "Rduckhts")))
expect_true(file.exists(system.file(
  "extdata",
  "vcf_file.bcf",
  package = "Rduckhts"
)))

# Test parameter validation - these should fail gracefully without a connection
expect_error(rduckhts_bcf(NULL, "test", "nonexistent.vcf"))
expect_error(rduckhts_bam(NULL, "test", "nonexistent.bam"))
expect_error(rduckhts_fasta(NULL, "test", "nonexistent.fa"))
expect_error(rduckhts_fastq(NULL, "test", "nonexistent.fq"))
expect_error(rduckhts_gff(NULL, "test", "nonexistent.gff"))
expect_error(rduckhts_gtf(NULL, "test", "nonexistent.gtf"))
expect_error(rduckhts_tabix(NULL, "test", "nonexistent.bed.gz"))
expect_error(rduckhts_hts_header(NULL, "nonexistent.bcf"))
expect_error(rduckhts_hts_index(NULL, "nonexistent.bcf"))
expect_error(rduckhts_hts_index_spans(NULL, "nonexistent.bcf"))
expect_error(rduckhts_hts_index_raw(NULL, "nonexistent.bcf"))

message("All basic tests passed!")
