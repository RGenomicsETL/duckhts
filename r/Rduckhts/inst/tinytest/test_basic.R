# Basic functionality tests for Rduckhts package
library(tinytest)

# Test package loading
expect_true(requireNamespace("Rduckhts", quietly = TRUE))

# Test basic functions exist
expect_true(exists("duckhts_build"))
expect_true(exists("rduckhts_load"))
expect_true(exists("rduckhts_bcf"))
expect_true(exists("rduckhts_bam"))
expect_true(exists("rduckhts_pileup"))
expect_true(exists("rduckhts_bam_index"))
expect_true(exists("rduckhts_bcf_index"))
expect_true(exists("rduckhts_bgzip"))
expect_true(exists("rduckhts_bgunzip"))
expect_true(exists("rduckhts_bam_bin_counts"))
expect_true(exists("rduckhts_bam_bed_coverage"))
expect_true(exists("rduckhts_mosdepth"))
expect_true(exists("rduckhts_samtools_idxstats"))
expect_true(exists("rduckhts_fasta"))
expect_true(exists("rduckhts_fasta_index"))
expect_true(exists("rduckhts_fastq"))
expect_true(exists("rduckhts_bigwig"))
expect_true(exists("rduckhts_htslib_config"))
expect_true(exists("rduckhts_htslib_info"))
expect_true(exists("rduckhts_htslib_version"))
expect_true(exists("rduckhts_detect_quality_encoding"))
expect_true(exists("rduckhts_gff"))
expect_true(exists("rduckhts_gtf"))
expect_true(exists("rduckhts_tabix"))
expect_true(exists("rduckhts_tabix_index"))
expect_true(exists("rduckhts_functions"))
expect_true(exists("rduckhts_hts_header"))
expect_true(exists("rduckhts_hts_index"))
expect_true(exists("rduckhts_hts_index_spans"))
expect_true(exists("rduckhts_hts_index_raw"))
expect_true(exists("rduckhts_liftover"))
expect_true(exists("rduckhts_bcf_convert_parquet"))
expect_true(exists("rduckhts_bam_convert_parquet"))
expect_true(exists("rduckhts_gff_convert_parquet"))
expect_true(exists("rduckhts_tabix_convert_parquet"))

# Test function signatures
expect_identical(
  names(formals(rduckhts_load)),
  c("con", "extension_path")
)
expect_identical(
  names(formals(rduckhts_bcf)),
  c("con", "table_name", "path", "region", "index_path", "tidy_format",
    "additional_csq_column_types", "scan_mode", "decompression_threads",
    "decode_error_policy", "overwrite")
)
expect_identical(
  names(formals(rduckhts_bam)),
  c("con", "table_name", "path", "region", "index_path", "reference",
    "standard_tags", "auxiliary_tags", "sequence_encoding",
    "quality_representation", "cigar_representation", "scan_mode", "decompression_threads", "overwrite")
)
expect_identical(
  names(formals(rduckhts_pileup)),
  c("con", "table_name", "path", "region", "index_path", "min_mapq", "flag_mask", "overwrite")
)
expect_identical(
  names(formals(rduckhts_bam_index)),
  c("con", "path", "index_path", "min_shift", "threads")
)
expect_identical(
  names(formals(rduckhts_bcf_index)),
  c("con", "path", "index_path", "min_shift", "threads")
)
expect_identical(
  names(formals(rduckhts_bgzip)),
  c("con", "path", "output_path", "threads", "level", "keep", "overwrite")
)
expect_identical(
  names(formals(rduckhts_bgunzip)),
  c("con", "path", "output_path", "threads", "keep", "overwrite")
)
expect_identical(
  names(formals(rduckhts_bam_bin_counts)),
  c("con", "path", "bin_width", "chrom", "include_unmapped", "reference", "index_path",
    "mapq", "require_flags", "exclude_flags", "rmdup", "stats")
)
expect_identical(
  names(formals(rduckhts_bam_bed_coverage)),
  c("con", "path", "bed_path", "reference", "index_path", "bed_index_path", "mapq",
    "min_baseq", "min_read_len", "require_flags", "exclude_flags", "min_depth",
    "max_depth", "decompression_threads", "fragment_mode", "strand_outputs", "processing_threads")
)
expect_identical(
  names(formals(rduckhts_mosdepth)),
  c("con", "prefix", "path", "chrom", "by", "fasta", "read_groups", "no_per_base",
    "threads", "processing_threads", "flag", "include_flag", "fast_mode", "fragment_mode", "use_median", "mapq", "min_frag_len",
    "max_frag_len", "precision_digits", "quantize", "thresholds", "index_path",
    "overwrite")
)
expect_identical(
  names(formals(rduckhts_samtools_idxstats)),
  c("con", "path", "output", "index_path", "threads", "overwrite")
)
expect_identical(
  names(formals(rduckhts_fasta)),
  c("con", "table_name", "path", "region", "index_path",
    "gzi_path", "sequence_encoding", "scan_mode", "overwrite")
)
expect_identical(
  names(formals(rduckhts_fasta_index)),
  c("con", "path", "index_path")
)
expect_identical(
  names(formals(rduckhts_fastq)),
  c("con", "table_name", "path", "mate_path", "interleaved",
    "sequence_encoding", "quality_representation", "input_quality_encoding",
    "scan_mode", "overwrite")
)
expect_identical(
  names(formals(rduckhts_bigwig)),
  c("con", "table_name", "path", "region", "blocks_per_iteration", "overwrite")
)
expect_identical(
  names(formals(rduckhts_htslib_config)),
  c("link", "validate")
)
expect_identical(names(formals(rduckhts_htslib_info)), "con")
expect_identical(names(formals(rduckhts_htslib_version)), "con")
expect_identical(
  names(formals(rduckhts_detect_quality_encoding)),
  c("con", "path", "max_records")
)
expect_identical(
  names(formals(rduckhts_gff)),
  c("con", "table_name", "path", "region", "index_path", "header",
    "header_names", "auto_detect", "column_types", "scan_mode", "attributes_map",
    "attributes_list", "attributes_pairs", "strict", "overwrite")
)
expect_identical(
  names(formals(rduckhts_gtf)),
  c("con", "table_name", "path", "region", "index_path", "header",
    "header_names", "auto_detect", "column_types", "scan_mode", "attributes_map",
    "attributes_list", "attributes_pairs", "overwrite")
)
expect_identical(
  names(formals(rduckhts_tabix)),
  c("con", "table_name", "path", "region", "index_path", "header",
    "header_names", "auto_detect", "column_types", "scan_mode", "overwrite")
)
expect_identical(
  names(formals(rduckhts_tabix_index)),
  c("con", "path", "preset", "index_path", "min_shift", "threads",
    "seq_col", "start_col", "end_col", "comment_char", "skip_lines")
)
expect_identical(
  names(formals(rduckhts_functions)),
  c("category", "kind")
)
expect_identical(
  names(formals(rduckhts_hts_header)),
  c("con", "path", "format", "mode")
)
expect_identical(
  names(formals(rduckhts_hts_index)),
  c("con", "path", "format", "index_path")
)
expect_identical(
  names(formals(rduckhts_hts_index_spans)),
  c("con", "path", "format", "index_path")
)
expect_identical(
  names(formals(rduckhts_hts_index_raw)),
  c("con", "path", "format", "index_path")
)
expect_identical(
  names(formals(rduckhts_liftover)),
  c("con", "query", "chain_path", "dst_fasta_ref", "chrom_col", "pos_col",
    "ref_col", "alt_col", "src_fasta_ref", "max_snp_gap", "max_indel_inc",
    "lift_mt", "end_pos_col", "no_left_align")
)
expect_identical(
  names(formals(rduckhts_bcf_convert_parquet)),
  c("con", "path", "output", "columns", "region", "index_path", "tidy_format",
    "additional_csq_column_types", "decompression_threads", "where", "compression",
    "row_group_size", "partition_by", "include_metadata", "header_text", "metadata",
    "metadata_json_file", "write_format_version", "overwrite")
)
expect_identical(
  names(formals(rduckhts_bam_convert_parquet)),
  c("con", "path", "output", "columns", "region", "index_path", "reference",
    "standard_tags", "auxiliary_tags", "sequence_encoding", "quality_representation",
    "cigar_representation", "decompression_threads", "where", "compression",
    "row_group_size", "partition_by", "include_metadata", "header_text", "metadata",
    "metadata_json_file", "write_format_version", "overwrite")
)
expect_identical(
  names(formals(rduckhts_gff_convert_parquet)),
  c("con", "path", "output", "columns", "region", "index_path", "header",
    "header_names", "auto_detect", "column_types", "attributes_map", "attributes_list",
    "attributes_pairs", "strict", "where", "compression", "row_group_size", "partition_by",
    "include_metadata", "header_text", "metadata", "metadata_json_file", "write_format_version", "overwrite")
)
expect_identical(
  names(formals(rduckhts_tabix_convert_parquet)),
  c("con", "path", "output", "columns", "region", "index_path", "header",
    "header_names", "auto_detect", "column_types", "where", "compression", "row_group_size",
    "partition_by", "include_metadata", "header_text", "metadata", "metadata_json_file", "write_format_version",
    "overwrite")
)
expect_identical(
  names(formals(rduckhts_munge)),
  c("con", "query", "fasta_ref", "preset", "column_map", "column_map_file",
    "iffy_tag", "mismatch_tag", "ns", "nc", "ne")
)

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
expect_true(file.exists(system.file("extdata", "legacy_phred64.fq", package = "Rduckhts")))
expect_true(file.exists(system.file(
  "extdata",
  "vcf_file.bcf",
  package = "Rduckhts"
)))
expect_true(file.exists(system.file(
  "function_catalog",
  "functions.tsv",
  package = "Rduckhts"
)))

# The installed rebuild path and package bootstrap share one DuckVEP source
# inventory. Every listed source must be present in the bundled tree.
duckvep_kernel_sources <- getFromNamespace(
  "duckhts_duckvep_kernel_source_files",
  "Rduckhts"
)()
expect_true(length(duckvep_kernel_sources) > 0L)
expect_true(all(file.exists(file.path(
  system.file("duckhts_extension", "duckvep", "kernel", "src", package = "Rduckhts"),
  duckvep_kernel_sources
))))

catalog <- rduckhts_functions()
expect_true(is.data.frame(catalog))
expect_true(all(c("name", "kind", "category", "signature", "description") %in% names(catalog)))
expect_true("seq_revcomp" %in% catalog$name)
expect_true("seq_encode_4bit" %in% catalog$name)
expect_true("seq_decode_4bit" %in% catalog$name)
expect_true("cigar_has_soft_clip" %in% catalog$name)
expect_true("cigar_reference_length" %in% catalog$name)
expect_true("sam_flag_bits" %in% catalog$name)
expect_true("read_bcf" %in% catalog$name)
expect_true("bgzip" %in% catalog$name)
expect_true("bam_index" %in% catalog$name)
expect_true("detect_quality_encoding" %in% catalog$name)
expect_true("duckhts_cgranges_create" %in% catalog$name)
expect_true("duckhts_cgranges_has_overlap" %in% catalog$name)
expect_true("duckhts_cgranges_count_overlaps" %in% catalog$name)
expect_true("duckhts_cgranges_overlaps" %in% catalog$name)
expect_true("bcftools_liftover" %in% catalog$name)
expect_true("duckdb_liftover" %in% catalog$name)
expect_true("bcftools_norm_row" %in% catalog$name)
expect_true("duckhts_bcftools_norm" %in% catalog$name)
expect_equal(unique(rduckhts_functions(kind = "scalar")$kind), "scalar")
expect_equal(unique(rduckhts_functions(category = "Readers")$category), "Readers")
expect_equal(unique(rduckhts_functions(category = "CIGAR Utils")$category), "CIGAR Utils")

# Test parameter validation - these should fail gracefully without a connection
expect_error(rduckhts_bcf(NULL, "test", "nonexistent.vcf"))
expect_error(rduckhts_bam(NULL, "test", "nonexistent.bam"))
expect_error(rduckhts_bam(NULL, "test", "nonexistent.bam", decompression_threads = -1))
expect_error(rduckhts_bam(NULL, "test", "nonexistent.bam", decompression_threads = 1.5))
expect_error(rduckhts_bam(NULL, "test", "nonexistent.bam", scan_mode = "parallel"))
expect_error(rduckhts_bam_index(NULL, "nonexistent.bam"))
expect_error(rduckhts_bcf_index(NULL, "nonexistent.vcf.gz"))
expect_error(rduckhts_bgzip(NULL, "nonexistent.txt"))
expect_error(rduckhts_bgunzip(NULL, "nonexistent.txt.gz"))
expect_error(rduckhts_fasta(NULL, "test", "nonexistent.fa"))
expect_error(rduckhts_fasta_index(NULL, "nonexistent.fa"))
expect_error(rduckhts_fastq(NULL, "test", "nonexistent.fq"))
expect_error(rduckhts_detect_quality_encoding(NULL, "nonexistent.fq"))
expect_error(rduckhts_gff(NULL, "test", "nonexistent.gff"))
expect_error(rduckhts_gtf(NULL, "test", "nonexistent.gtf"))
expect_error(rduckhts_tabix(NULL, "test", "nonexistent.bed.gz"))
expect_error(rduckhts_tabix_index(NULL, "nonexistent.bed.gz"))
expect_error(rduckhts_hts_header(NULL, "nonexistent.bcf"))
expect_error(rduckhts_hts_index(NULL, "nonexistent.bcf"))
expect_error(rduckhts_hts_index_spans(NULL, "nonexistent.bcf"))
expect_error(rduckhts_hts_index_raw(NULL, "nonexistent.bcf"))
expect_error(rduckhts_liftover(NULL, "SELECT 1 AS chrom, 1 AS pos", "x.chain", "y.fa"))

# Test multi-file reading functions exist
expect_true(exists("rduckhts_bam_multi"))
expect_true(exists("rduckhts_bcf_multi"))
expect_true(exists("rduckhts_fastq_multi"))
expect_true(exists("rduckhts_fasta_multi"))
expect_true(exists("rduckhts_bed_multi"))
expect_true(exists("rduckhts_tabix_multi"))
expect_true(exists("rduckhts_gff_multi"))
expect_true(exists("rduckhts_gtf_multi"))

# Test multi-file reading function signatures
expect_identical(
  names(formals(rduckhts_bam_multi)),
  c("con", "table_name", "files", "region", "index_path", "reference",
    "standard_tags", "auxiliary_tags", "sequence_encoding",
    "quality_representation", "cigar_representation", "scan_mode", "decompression_threads", ".params", "overwrite")
)
expect_identical(
  names(formals(rduckhts_bcf_multi)),
  c("con", "table_name", "files", "region", "index_path", "tidy_format",
    "additional_csq_column_types", "scan_mode", "decompression_threads", ".params", "overwrite")
)
expect_identical(
  names(formals(rduckhts_fastq_multi)),
  c("con", "table_name", "files", "mate_path", "interleaved",
    "sequence_encoding", "quality_representation",
    "input_quality_encoding", "scan_mode", ".params", "overwrite")
)
expect_identical(
  names(formals(rduckhts_fasta_multi)),
  c("con", "table_name", "files", "region", "index_path",
    "gzi_path", "sequence_encoding", "scan_mode", ".params", "overwrite")
)
expect_identical(
  names(formals(rduckhts_bed_multi)),
  c("con", "table_name", "files", "region", "index_path", "scan_mode", ".params",
    "overwrite")
)
expect_identical(
  names(formals(rduckhts_tabix_multi)),
  c("con", "table_name", "files", "region", "index_path", "header",
    "header_names", "auto_detect", "column_types", "scan_mode", ".params", "overwrite")
)
expect_identical(
  names(formals(rduckhts_gff_multi)),
  c("con", "table_name", "files", "region", "index_path", "header",
    "header_names", "auto_detect", "column_types", "scan_mode", "attributes_map",
    "attributes_list", "attributes_pairs", "strict", ".params", "overwrite")
)
expect_identical(
  names(formals(rduckhts_gtf_multi)),
  c("con", "table_name", "files", "region", "index_path", "header",
    "header_names", "auto_detect", "column_types", "scan_mode", "attributes_map",
    "attributes_list", "attributes_pairs", ".params", "overwrite")
)

# Test .params validation (these should fail with bad input, no connection needed)
expect_error(rduckhts_bam_multi(NULL, "t", "*.bam", .params = "not a data.frame"))
expect_error(rduckhts_bam_multi(NULL, "t", "*.bam", .params = data.frame(x = 1)))
expect_error(rduckhts_bam_multi(NULL, "t", "*.bam", decompression_threads = -1))
expect_error(rduckhts_bam_multi(NULL, "t", "*.bam", scan_mode = "parallel"))

message("All basic tests passed!")
