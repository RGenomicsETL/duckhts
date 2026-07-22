/**
 * DuckHTS Extension Entry Point
 *
 * Registers all HTS reader table functions with DuckDB.
 */

#define DUCKDB_EXTENSION_NAME duckhts

#include "duckdb_extension.h"
#include <htslib/hts.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "duckhts_simd.h"
#include "wasm_http_hfile.h"

DUCKDB_EXTENSION_EXTERN

/* bcf_reader.c */
extern void register_read_bcf_function(duckdb_connection connection);
extern void register_read_bcf_v2_function(duckdb_connection connection);
extern void register_read_bcf_appender_function(duckdb_connection connection, duckdb_database database);
/* bam_reader.c */
extern void register_read_bam_function(duckdb_connection connection);
/* bam_pileup.c */
extern void register_read_pileup_function(duckdb_connection connection);
/* seq_reader.c */
extern void register_read_fasta_function(duckdb_connection connection);
extern void register_read_fastq_function(duckdb_connection connection);
extern void register_fasta_index_function(duckdb_connection connection);
/* fastq_qc.c */
extern void register_duckhts_fastq_qc_function(duckdb_connection connection);
/* interval_udf.c */
extern void register_read_bed_function(duckdb_connection connection);
extern void register_fasta_nuc_function(duckdb_connection connection);
/* bigwig_reader.c */
extern void register_read_bigwig_function(duckdb_connection connection);
/* bgzip.c */
extern void register_bgzip_function(duckdb_connection connection);
extern void register_bgunzip_function(duckdb_connection connection);
/* hts_index_builder.c */
extern void register_bam_index_function(duckdb_connection connection);
extern void register_bcf_index_function(duckdb_connection connection);
extern void register_tabix_index_function(duckdb_connection connection);
/* liftover_udf.c */
extern void register_liftover_functions(duckdb_connection connection);
/* bcftools_norm_udf.c */
extern void register_bcftools_norm_functions(duckdb_connection connection);
/* munge_udf.c */
extern void register_munge_functions(duckdb_connection connection);
/* kmer_udf.c */
extern void register_kmer_udf_functions(duckdb_connection connection);
/* tabix_reader.c */
extern void register_read_tabix_function(duckdb_connection connection);
extern void register_read_gtf_function(duckdb_connection connection);
extern void register_read_gff_function(duckdb_connection connection);
/* hts_meta_reader.c */
extern void register_read_hts_header_function(duckdb_connection connection);
extern void register_read_hts_index_function(duckdb_connection connection);
extern void register_read_hts_index_spans_function(duckdb_connection connection);
/* quality_encoding_reader.c */
extern void register_detect_quality_encoding_function(duckdb_connection connection);
/* score_udf.c */
extern void register_bcftools_score_function(duckdb_connection connection);
/* mosdepth_table.c */
extern void register_duckhts_mosdepth_function(duckdb_connection connection);
/* bam_bin_counts.c */
extern void register_bam_bin_counts_function(duckdb_connection connection);
/* samtools_idxstats_table.c */
extern void register_duckhts_samtools_idxstats_function(duckdb_connection connection);
/* bam_bed_coverage.c */
extern void register_duckhts_bam_bed_coverage_function(duckdb_connection connection);
/* cgranges_api.c */
extern void register_duckhts_cgranges_functions(duckdb_connection connection, duckdb_database database);
/* duckvep resident model and annotation adapter */
extern void register_duckvep_functions(duckdb_connection connection, duckdb_database database);
extern bool register_duckvep_ensembl_functions(duckdb_connection connection);
extern bool register_duckvep_sql_functions(duckdb_connection connection);
/* variantkey_udf.c */
extern void register_variantkey_functions(duckdb_connection connection);
/* simd/duckhts_simd_dispatch.c */
extern void register_duckhts_simd_functions(duckdb_connection connection);

static const char *duckhts_loaded_htslib_version = "";
static unsigned int duckhts_loaded_htslib_features = 0;
static char duckhts_loaded_htslib_feature_string[1200];
static pthread_once_t duckhts_htslib_snapshot_once = PTHREAD_ONCE_INIT;

static void duckhts_snapshot_htslib(void) {
    const char *version = hts_version();
    const char *feature_string = hts_feature_string();

    duckhts_loaded_htslib_version = version ? version : "";
    duckhts_loaded_htslib_features = hts_features();
    snprintf(duckhts_loaded_htslib_feature_string,
             sizeof(duckhts_loaded_htslib_feature_string), "%s",
             feature_string ? feature_string : "");
}

static void duckhts_htslib_version_scalar(duckdb_function_info info,
                                           duckdb_data_chunk input,
                                           duckdb_vector output) {
    idx_t row_count = duckdb_data_chunk_get_size(input);
    idx_t row;
    (void)info;

    for (row = 0; row < row_count; row++) {
        duckdb_vector_assign_string_element(output, row, duckhts_loaded_htslib_version);
    }
}

static void duckhts_htslib_feature_string_scalar(duckdb_function_info info,
                                                  duckdb_data_chunk input,
                                                  duckdb_vector output) {
    idx_t row_count = duckdb_data_chunk_get_size(input);
    idx_t row;
    (void)info;

    for (row = 0; row < row_count; row++) {
        duckdb_vector_assign_string_element(output, row, duckhts_loaded_htslib_feature_string);
    }
}

static void duckhts_htslib_features_scalar(duckdb_function_info info,
                                            duckdb_data_chunk input,
                                            duckdb_vector output) {
    uint32_t *values = (uint32_t *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);
    idx_t row;
    (void)info;

    for (row = 0; row < row_count; row++) {
        values[row] = (uint32_t)duckhts_loaded_htslib_features;
    }
}

static void register_duckhts_htslib_string_function(duckdb_connection connection,
                                                     const char *name,
                                                     duckdb_scalar_function_t callback) {
    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);

    duckdb_scalar_function_set_name(function, name);
    duckdb_scalar_function_set_return_type(function, varchar_type);
    duckdb_scalar_function_set_function(function, callback);
    duckdb_register_scalar_function(connection, function);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_scalar_function(&function);
}

static void register_duckhts_htslib_functions(duckdb_connection connection) {
    duckdb_scalar_function function;
    duckdb_logical_type uinteger_type;

    pthread_once(&duckhts_htslib_snapshot_once, duckhts_snapshot_htslib);

    register_duckhts_htslib_string_function(connection, "duckhts_htslib_version",
                                             duckhts_htslib_version_scalar);
    register_duckhts_htslib_string_function(connection, "duckhts_htslib_feature_string",
                                             duckhts_htslib_feature_string_scalar);

    function = duckdb_create_scalar_function();
    uinteger_type = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    duckdb_scalar_function_set_name(function, "duckhts_htslib_features");
    duckdb_scalar_function_set_return_type(function, uinteger_type);
    duckdb_scalar_function_set_function(function, duckhts_htslib_features_scalar);
    duckdb_register_scalar_function(connection, function);
    duckdb_destroy_logical_type(&uinteger_type);
    duckdb_destroy_scalar_function(&function);
}

static bool run_sql_or_fail(duckdb_connection connection, const char *sql) {
    duckdb_result result;
    duckdb_state state = duckdb_query(connection, sql, &result);
    if (state != DuckDBSuccess) {
        const char *err = duckdb_result_error(&result);
        if (err && *err) {
            fprintf(stderr, "[duckhts] failed SQL registration: %s\n", err);
        }
        duckdb_destroy_result(&result);
        return false;
    }
    duckdb_destroy_result(&result);
    return true;
}

static bool run_sql_parts_or_fail(duckdb_connection connection, const char *const *parts, size_t count) {
    size_t total_len = 0;
    size_t offset = 0;
    char *sql = NULL;
    size_t i;
    bool ok;

    for (i = 0; i < count; i++) {
        total_len += strlen(parts[i]);
    }

    sql = (char *)malloc(total_len + 1);
    if (!sql) {
        fprintf(stderr, "[duckhts] failed SQL registration: out of memory\n");
        return false;
    }

    for (i = 0; i < count; i++) {
        size_t len = strlen(parts[i]);
        memcpy(sql + offset, parts[i], len);
        offset += len;
    }
    sql[offset] = '\0';

    ok = run_sql_or_fail(connection, sql);
    free(sql);
    return ok;
}

DUCKDB_EXTENSION_ENTRYPOINT(duckdb_connection connection,
                            duckdb_extension_info info,
                            struct duckdb_extension_access* access) {
    (void)info;
    (void)access;
    register_wasm_http_hfile_backend();
    duckhts_simd_init();

    register_duckhts_htslib_functions(connection);
    register_duckhts_simd_functions(connection);
    register_read_bcf_function(connection);
    register_read_bcf_v2_function(connection);
    register_read_bcf_appender_function(connection, *access->get_database(info));
    register_read_bam_function(connection);
    register_read_pileup_function(connection);
    register_bcftools_norm_functions(connection);
    register_read_fasta_function(connection);
    register_read_fastq_function(connection);
    register_duckhts_fastq_qc_function(connection);
    register_fasta_index_function(connection);
    register_read_bed_function(connection);
    register_read_bigwig_function(connection);
    register_fasta_nuc_function(connection);
    register_bgzip_function(connection);
    register_bgunzip_function(connection);
    register_bam_index_function(connection);
    register_bcf_index_function(connection);
    register_tabix_index_function(connection);
    register_liftover_functions(connection);
    register_munge_functions(connection);
    register_kmer_udf_functions(connection);
    register_read_tabix_function(connection);
    register_read_gtf_function(connection);
    register_read_gff_function(connection);
    register_read_hts_header_function(connection);
    register_read_hts_index_function(connection);
    register_read_hts_index_spans_function(connection);
    register_detect_quality_encoding_function(connection);
    register_bcftools_score_function(connection);
    register_duckhts_mosdepth_function(connection);
    register_bam_bin_counts_function(connection);
    register_duckhts_samtools_idxstats_function(connection);
    register_duckhts_bam_bed_coverage_function(connection);
    register_variantkey_functions(connection);
    register_duckhts_cgranges_functions(connection, *access->get_database(info));
    register_duckvep_functions(connection, *access->get_database(info));
    if (!register_duckvep_ensembl_functions(connection)) {
        return false;
    }
    if (!register_duckvep_sql_functions(connection)) {
        return false;
    }
    if (!run_sql_or_fail(connection,
        "CREATE OR REPLACE MACRO duckhts_quote_ident(x) AS "
        "CASE WHEN x IS NULL THEN NULL ELSE '\"' || replace(x, '\"', '\"\"') || '\"' END")) {
        return false;
    }
    if (!run_sql_or_fail(connection,
        "CREATE OR REPLACE MACRO duckhts_quote_string(x) AS "
        "CASE WHEN x IS NULL THEN NULL ELSE '''' || replace(x, '''', '''''') || '''' END")) {
        return false;
    }
    if (!run_sql_or_fail(connection,
        "CREATE OR REPLACE MACRO duckhts_json_path_key(x) AS "
        "'$.\"' || replace(replace(x, '\\', '\\\\'), '\"', '\\\"') || '\"'")) {
        return false;
    }
    if (!run_sql_or_fail(connection,
        "CREATE OR REPLACE MACRO duckhts_duckdb_type_supported(candidate_type_name) AS ("
        "EXISTS (SELECT 1 FROM duckdb_types() AS dt WHERE lower(dt.type_name) = lower(candidate_type_name)))")) {
        return false;
    }
    if (!run_sql_or_fail(connection,
        "CREATE OR REPLACE MACRO duckhts_duckdb_supports_variant() AS ("
        "duckhts_duckdb_type_supported('VARIANT'))")) {
        return false;
    }
    if (!run_sql_or_fail(connection,
        "CREATE OR REPLACE MACRO duckhts_duckdb_supports_geometry() AS ("
        "duckhts_duckdb_type_supported('GEOMETRY'))")) {
        return false;
    }
    if (!run_sql_or_fail(connection,
        "CREATE OR REPLACE MACRO duckhts_parquet_ident_list(cols) AS ("
        "CASE WHEN len(cols) = 0 THEN '*' "
        "ELSE (SELECT string_agg(duckhts_quote_ident(x), ', ') FROM unnest(cols) AS t(x)) END)")) {
        return false;
    }
    if (!run_sql_or_fail(connection,
        "CREATE OR REPLACE MACRO duckhts_raw_header_text(path, format_hint) AS ("
        "(SELECT coalesce(string_agg(raw, '\\n' ORDER BY idx), '') "
        "FROM read_hts_header(path, format := format_hint, mode := 'raw')))")) {
        return false;
    }
    {
        static const char *const parquet_metadata_sql[] = {
            "CREATE OR REPLACE MACRO duckhts_parquet_metadata_args(",
            "source_format, reader, source_path, header_key, header_value, header_corrected, ",
            "columns, partition_by, filter_sql, layout_metadata := map([]::VARCHAR[], []::VARCHAR[]), ",
            "metadata := map([]::VARCHAR[], []::VARCHAR[]), metadata_json_file := NULL, write_format_version := '1') AS (",
            "SELECT string_agg(duckhts_quote_ident(key) || ' := ' || duckhts_quote_string(value), ', ' ORDER BY key) ",
            "FROM (SELECT key, value, row_number() OVER (PARTITION BY key ORDER BY priority DESC) AS rn ",
            "FROM (",
            "SELECT * FROM (VALUES ",
            "('duckhts_write_format_version', write_format_version, 10), ",
            "('duckhts_source_format', source_format, 10), ",
            "('duckhts_reader', reader, 10), ",
            "('duckhts_source_path', source_path, 10), ",
            "('duckhts_header_corrected', header_corrected, 10), ",
            "('duckhts_filter', coalesce(filter_sql, ''), 10), ",
            "('duckhts_columns', CASE WHEN len(columns) = 0 THEN '*' ELSE array_to_string(columns, ',') END, 10), ",
            "('duckhts_partition_by', CASE WHEN len(partition_by) = 0 THEN '' ELSE array_to_string(partition_by, ',') END, 10), ",
            "(header_key, coalesce(header_value, ''), 10)) AS defaults(key, value, priority) ",
            "UNION ALL SELECT e.key, e.value, 20 FROM unnest(map_entries(layout_metadata)) AS t(e) ",
            "UNION ALL SELECT key, value, 30 FROM query(CASE WHEN metadata_json_file IS NULL THEN ",
            "'SELECT NULL::VARCHAR AS key, NULL::VARCHAR AS value WHERE false' ELSE ",
            "'SELECT key, CASE ",
            "WHEN json_type(json_extract(content, duckhts_json_path_key(key))) IN (''OBJECT'', ''ARRAY'') THEN CAST(json_extract(content, duckhts_json_path_key(key)) AS VARCHAR) ",
            "WHEN json_type(json_extract(content, duckhts_json_path_key(key))) = ''NULL'' THEN '''' ",
            "ELSE json_extract_string(content, duckhts_json_path_key(key)) END AS value ",
            "FROM read_text(' || duckhts_quote_string(metadata_json_file) || '), ",
            "unnest(CASE WHEN json_type(content) = ''OBJECT'' THEN json_keys(content) ",
            "ELSE error(''metadata_json_file must contain a top-level JSON object'') END) AS q(key)' END) ",
            "UNION ALL SELECT e.key, e.value, 40 FROM unnest(map_entries(metadata)) AS t(e)",
            ") all_metadata) ranked WHERE rn = 1)"
        };
        if (!run_sql_parts_or_fail(connection, parquet_metadata_sql,
                                   sizeof(parquet_metadata_sql) / sizeof(parquet_metadata_sql[0]))) {
            return false;
        }
    }
    {
        static const char *const parquet_copy_sql[] = {
            "CREATE OR REPLACE MACRO duckhts_parquet_copy_sql(reader, path, output, reader_args := '', ",
            "source_format := '', header_key := '', header_value := '', header_corrected := 'false', ",
            "columns := []::VARCHAR[], partition_by := []::VARCHAR[], filter_sql := NULL, ",
            "layout_metadata := map([]::VARCHAR[], []::VARCHAR[]), metadata := map([]::VARCHAR[], []::VARCHAR[]), ",
            "metadata_json_file := NULL, compression := 'zstd', row_group_size := 100000, include_metadata := true, overwrite := false, ",
            "write_format_version := '1') AS (",
            "'COPY (SELECT ' || duckhts_parquet_ident_list(columns) || ' FROM ' || reader || '(' || duckhts_quote_string(path) || ",
            "CASE WHEN length(reader_args) > 0 THEN ', ' || reader_args ELSE '' END || ')' || ",
            "CASE WHEN filter_sql IS NULL OR length(filter_sql) = 0 THEN '' ELSE ' WHERE (' || filter_sql || ')' END || ",
            "') TO ' || duckhts_quote_string(output) || ' (FORMAT PARQUET, COMPRESSION ' || duckhts_quote_string(upper(compression)) || ",
            "', ROW_GROUP_SIZE ' || CAST(row_group_size AS VARCHAR) || ', OVERWRITE ' || CASE WHEN overwrite THEN 'true' ELSE 'false' END || ",
            "CASE WHEN len(partition_by) = 0 THEN '' ELSE ', PARTITION_BY (' || duckhts_parquet_ident_list(partition_by) || ')' END || ",
            "CASE WHEN include_metadata THEN ', KV_METADATA struct_pack(' || ",
            "duckhts_parquet_metadata_args(source_format, reader, path, header_key, header_value, header_corrected, columns, partition_by, filter_sql, layout_metadata, metadata, metadata_json_file, write_format_version) || ')' ",
            "ELSE '' END || ')' )"
        };
        if (!run_sql_parts_or_fail(connection, parquet_copy_sql,
                                   sizeof(parquet_copy_sql) / sizeof(parquet_copy_sql[0]))) {
            return false;
        }
    }
    {
        static const char *const bcf_convert_sql[] = {
            "CREATE OR REPLACE MACRO duckhts_bcf_convert_parquet_sql(path, output, columns := []::VARCHAR[], region := NULL, index_path := NULL, ",
            "tidy_format := false, additional_csq_column_types := NULL, decompression_threads := 0, where_sql := NULL, compression := 'zstd', ",
            "row_group_size := 100000, partition_by := []::VARCHAR[], include_metadata := true, header_text := NULL, ",
            "metadata := map([]::VARCHAR[], []::VARCHAR[]), metadata_json_file := NULL, overwrite := false, write_format_version := '1') AS (",
            "duckhts_parquet_copy_sql('read_bcf', path, output, reader_args := ltrim(",
            "(CASE WHEN region IS NULL THEN '' ELSE ', region := ' || duckhts_quote_string(region) END) || ",
            "(CASE WHEN index_path IS NULL THEN '' ELSE ', index_path := ' || duckhts_quote_string(index_path) END) || ",
            "', tidy_format := ' || CASE WHEN tidy_format THEN 'true' ELSE 'false' END || ",
            "(CASE WHEN additional_csq_column_types IS NULL THEN '' ELSE ', additional_csq_column_types := ' || duckhts_quote_string(additional_csq_column_types) END) || ",
            "', decompression_threads := ' || CAST(decompression_threads AS VARCHAR), ', '), ",
            "source_format := 'vcf_bcf', header_key := 'vcf_header', ",
            "header_value := coalesce(header_text, duckhts_raw_header_text(path, 'bcf')), ",
            "header_corrected := CASE WHEN header_text IS NULL THEN 'false' ELSE 'true' END, ",
            "columns := columns, partition_by := partition_by, filter_sql := where_sql, ",
            "layout_metadata := map(['duckhts_tidy_format', 'duckhts_region'], [CASE WHEN tidy_format THEN 'true' ELSE 'false' END, coalesce(region, '')]), ",
            "metadata := metadata, metadata_json_file := metadata_json_file, compression := compression, row_group_size := row_group_size, ",
            "include_metadata := include_metadata, overwrite := overwrite, write_format_version := write_format_version))"
        };
        if (!run_sql_parts_or_fail(connection, bcf_convert_sql,
                                   sizeof(bcf_convert_sql) / sizeof(bcf_convert_sql[0]))) {
            return false;
        }
    }
    {
        static const char *const bam_convert_sql[] = {
            "CREATE OR REPLACE MACRO duckhts_bam_convert_parquet_sql(path, output, columns := []::VARCHAR[], region := NULL, index_path := NULL, ",
            "reference := NULL, standard_tags := false, auxiliary_tags := false, sequence_encoding := NULL, quality_representation := NULL, ",
            "cigar_representation := NULL, decompression_threads := 2, where_sql := NULL, compression := 'zstd', row_group_size := 100000, ",
            "partition_by := []::VARCHAR[], include_metadata := true, header_text := NULL, metadata := map([]::VARCHAR[], []::VARCHAR[]), ",
            "metadata_json_file := NULL, overwrite := false, write_format_version := '1') AS (",
            "duckhts_parquet_copy_sql('read_bam', path, output, reader_args := ltrim(",
            "(CASE WHEN region IS NULL THEN '' ELSE ', region := ' || duckhts_quote_string(region) END) || ",
            "(CASE WHEN index_path IS NULL THEN '' ELSE ', index_path := ' || duckhts_quote_string(index_path) END) || ",
            "(CASE WHEN reference IS NULL THEN '' ELSE ', reference := ' || duckhts_quote_string(reference) END) || ",
            "', standard_tags := ' || CASE WHEN standard_tags THEN 'true' ELSE 'false' END || ",
            "', auxiliary_tags := ' || CASE WHEN auxiliary_tags THEN 'true' ELSE 'false' END || ",
            "(CASE WHEN sequence_encoding IS NULL THEN '' ELSE ', sequence_encoding := ' || duckhts_quote_string(sequence_encoding) END) || ",
            "(CASE WHEN quality_representation IS NULL THEN '' ELSE ', quality_representation := ' || duckhts_quote_string(quality_representation) END) || ",
            "(CASE WHEN cigar_representation IS NULL THEN '' ELSE ', cigar_representation := ' || duckhts_quote_string(cigar_representation) END) || ",
            "', decompression_threads := ' || CAST(decompression_threads AS VARCHAR), ', '), ",
            "source_format := 'sam_bam_cram', header_key := 'sam_header', ",
            "header_value := coalesce(header_text, duckhts_raw_header_text(path, 'bam')), ",
            "header_corrected := CASE WHEN header_text IS NULL THEN 'false' ELSE 'true' END, ",
            "columns := columns, partition_by := partition_by, filter_sql := where_sql, ",
            "layout_metadata := map(['duckhts_region', 'duckhts_standard_tags', 'duckhts_auxiliary_tags', 'duckhts_sequence_encoding', 'duckhts_quality_representation', 'duckhts_cigar_representation'], ",
            "[coalesce(region, ''), CASE WHEN standard_tags THEN 'true' ELSE 'false' END, CASE WHEN auxiliary_tags THEN 'true' ELSE 'false' END, coalesce(sequence_encoding, 'string'), coalesce(quality_representation, 'string'), coalesce(cigar_representation, 'string')]), ",
            "metadata := metadata, metadata_json_file := metadata_json_file, compression := compression, row_group_size := row_group_size, ",
            "include_metadata := include_metadata, overwrite := overwrite, write_format_version := write_format_version))"
        };
        if (!run_sql_parts_or_fail(connection, bam_convert_sql,
                                   sizeof(bam_convert_sql) / sizeof(bam_convert_sql[0]))) {
            return false;
        }
    }
    {
        static const char *const gff_convert_sql[] = {
            "CREATE OR REPLACE MACRO duckhts_gff_convert_parquet_sql(path, output, columns := []::VARCHAR[], region := NULL, index_path := NULL, ",
            "header := NULL, header_names := []::VARCHAR[], auto_detect := NULL, column_types := []::VARCHAR[], attributes_map := false, ",
            "attributes_list := false, attributes_pairs := false, strict := false, where_sql := NULL, compression := 'zstd', row_group_size := 100000, ",
            "partition_by := []::VARCHAR[], include_metadata := true, header_text := NULL, metadata := map([]::VARCHAR[], []::VARCHAR[]), metadata_json_file := NULL, overwrite := false, write_format_version := '1') AS (",
            "duckhts_parquet_copy_sql('read_gff', path, output, reader_args := ltrim(",
            "(CASE WHEN region IS NULL THEN '' ELSE ', region := ' || duckhts_quote_string(region) END) || ",
            "(CASE WHEN index_path IS NULL THEN '' ELSE ', index_path := ' || duckhts_quote_string(index_path) END) || ",
            "(CASE WHEN header IS NULL THEN '' ELSE ', header := ' || CASE WHEN header THEN 'true' ELSE 'false' END END) || ",
            "(CASE WHEN len(header_names) = 0 THEN '' ELSE ', header_names := [' || (SELECT string_agg(duckhts_quote_string(x), ', ') FROM unnest(header_names) AS t(x)) || ']' END) || ",
            "(CASE WHEN auto_detect IS NULL THEN '' ELSE ', auto_detect := ' || CASE WHEN auto_detect THEN 'true' ELSE 'false' END END) || ",
            "(CASE WHEN len(column_types) = 0 THEN '' ELSE ', column_types := [' || (SELECT string_agg(duckhts_quote_string(x), ', ') FROM unnest(column_types) AS t(x)) || ']' END) || ",
            "', attributes_map := ' || CASE WHEN attributes_map THEN 'true' ELSE 'false' END || ",
            "', attributes_list := ' || CASE WHEN attributes_list THEN 'true' ELSE 'false' END || ",
            "', attributes_pairs := ' || CASE WHEN attributes_pairs THEN 'true' ELSE 'false' END || ",
            "', strict := ' || CASE WHEN strict THEN 'true' ELSE 'false' END, ', '), ",
            "source_format := 'gff', header_key := 'gff_header', header_value := coalesce(header_text, duckhts_raw_header_text(path, 'tabix')), ",
            "header_corrected := CASE WHEN header_text IS NULL THEN 'false' ELSE 'true' END, columns := columns, partition_by := partition_by, filter_sql := where_sql, ",
            "layout_metadata := map(['duckhts_region', 'duckhts_header', 'duckhts_header_names', 'duckhts_auto_detect', 'duckhts_column_types', 'duckhts_attributes_map', 'duckhts_attributes_list', 'duckhts_attributes_pairs', 'duckhts_strict'], ",
            "[coalesce(region, ''), CASE WHEN header IS NULL THEN '' WHEN header THEN 'true' ELSE 'false' END, array_to_string(header_names, ','), CASE WHEN auto_detect IS NULL THEN '' WHEN auto_detect THEN 'true' ELSE 'false' END, array_to_string(column_types, ','), CASE WHEN attributes_map THEN 'true' ELSE 'false' END, CASE WHEN attributes_list THEN 'true' ELSE 'false' END, CASE WHEN attributes_pairs THEN 'true' ELSE 'false' END, CASE WHEN strict THEN 'true' ELSE 'false' END]), ",
            "metadata := metadata, metadata_json_file := metadata_json_file, compression := compression, row_group_size := row_group_size, include_metadata := include_metadata, overwrite := overwrite, write_format_version := write_format_version))"
        };
        if (!run_sql_parts_or_fail(connection, gff_convert_sql,
                                   sizeof(gff_convert_sql) / sizeof(gff_convert_sql[0]))) {
            return false;
        }
    }
    {
        static const char *const tabix_convert_sql[] = {
            "CREATE OR REPLACE MACRO duckhts_tabix_convert_parquet_sql(path, output, columns := []::VARCHAR[], region := NULL, index_path := NULL, ",
            "header := NULL, header_names := []::VARCHAR[], auto_detect := NULL, column_types := []::VARCHAR[], where_sql := NULL, compression := 'zstd', row_group_size := 100000, ",
            "partition_by := []::VARCHAR[], include_metadata := true, header_text := NULL, metadata := map([]::VARCHAR[], []::VARCHAR[]), metadata_json_file := NULL, overwrite := false, write_format_version := '1') AS (",
            "duckhts_parquet_copy_sql('read_tabix', path, output, reader_args := ltrim(",
            "(CASE WHEN region IS NULL THEN '' ELSE ', region := ' || duckhts_quote_string(region) END) || ",
            "(CASE WHEN index_path IS NULL THEN '' ELSE ', index_path := ' || duckhts_quote_string(index_path) END) || ",
            "(CASE WHEN header IS NULL THEN '' ELSE ', header := ' || CASE WHEN header THEN 'true' ELSE 'false' END END) || ",
            "(CASE WHEN len(header_names) = 0 THEN '' ELSE ', header_names := [' || (SELECT string_agg(duckhts_quote_string(x), ', ') FROM unnest(header_names) AS t(x)) || ']' END) || ",
            "(CASE WHEN auto_detect IS NULL THEN '' ELSE ', auto_detect := ' || CASE WHEN auto_detect THEN 'true' ELSE 'false' END END) || ",
            "(CASE WHEN len(column_types) = 0 THEN '' ELSE ', column_types := [' || (SELECT string_agg(duckhts_quote_string(x), ', ') FROM unnest(column_types) AS t(x)) || ']' END), ', '), ",
            "source_format := 'tabix', header_key := 'tabix_header', header_value := coalesce(header_text, duckhts_raw_header_text(path, 'tabix')), ",
            "header_corrected := CASE WHEN header_text IS NULL THEN 'false' ELSE 'true' END, columns := columns, partition_by := partition_by, filter_sql := where_sql, ",
            "layout_metadata := map(['duckhts_region', 'duckhts_header', 'duckhts_header_names', 'duckhts_auto_detect', 'duckhts_column_types'], ",
            "[coalesce(region, ''), CASE WHEN header IS NULL THEN '' WHEN header THEN 'true' ELSE 'false' END, array_to_string(header_names, ','), CASE WHEN auto_detect IS NULL THEN '' WHEN auto_detect THEN 'true' ELSE 'false' END, array_to_string(column_types, ',')]), ",
            "metadata := metadata, metadata_json_file := metadata_json_file, compression := compression, row_group_size := row_group_size, include_metadata := include_metadata, overwrite := overwrite, write_format_version := write_format_version))"
        };
        if (!run_sql_parts_or_fail(connection, tabix_convert_sql,
                                   sizeof(tabix_convert_sql) / sizeof(tabix_convert_sql[0]))) {
            return false;
        }
    }
    {
        static const char *const hts_union_query_sql[] = {
            "CREATE OR REPLACE MACRO hts_union_query(reader, pattern, params := '') AS (",
            "(SELECT string_agg(",
            "'SELECT ''' || replace(file, '''', '''''') || ''' AS filename, * FROM ' ",
            "|| reader || '(''' || replace(file, '''', '''''') || '''' ",
            "|| CASE WHEN length(params) > 0 THEN ', ' || params ELSE '' END ",
            "|| ')', ",
            "' UNION ALL BY NAME '",
            ") FROM glob(pattern) g(file))",
            ")"
        };
        if (!run_sql_parts_or_fail(connection, hts_union_query_sql,
                                    sizeof(hts_union_query_sql) / sizeof(hts_union_query_sql[0]))) {
            return false;
        }
    }
    {
        static const char *const hts_region_union_query_sql[] = {
            "CREATE OR REPLACE MACRO hts_region_union_query(reader, path, regions, params := '') AS (",
            "(SELECT string_agg(",
            "'SELECT ' || CAST(ord AS VARCHAR) || ' AS duckhts_region_shard_id, ' ",
            "|| duckhts_quote_string(region) || ' AS duckhts_region_shard, ' ",
            "|| duckhts_quote_string(path) || ' AS filename, * FROM ' ",
            "|| reader || '(' || duckhts_quote_string(path) ",
            "|| ', region := ' || duckhts_quote_string(region) ",
            "|| CASE WHEN length(params) > 0 THEN ', ' || params ELSE '' END ",
            "|| ')', ",
            "' UNION ALL BY NAME ' ORDER BY ord",
            ") FROM unnest(regions) WITH ORDINALITY AS t(region, ord))",
            ")"
        };
        if (!run_sql_parts_or_fail(connection, hts_region_union_query_sql,
                                    sizeof(hts_region_union_query_sql) / sizeof(hts_region_union_query_sql[0]))) {
            return false;
        }
    }
    if (!run_sql_or_fail(connection,
        "CREATE OR REPLACE MACRO duckdb_munge_preset_map(preset) AS "
        "CASE "
        "WHEN upper(preset) = 'PLINK' THEN map(['SNP','BP','CHR','A1','A2','P','OR','BETA','INFO','FRQ','SE'], ['SNP','BP','CHR','A1','A2','P','OR','BETA','INFO','FRQ','SE']) "
        "WHEN upper(preset) = 'PLINK2' THEN map(['SNP','BP','CHR','A1','A2','P','Z','OR','BETA','INFO','FRQ','SE','LP'], ['ID','POS','CHROM','A1','AX','P','Z_STAT','OR','BETA','MACH_R2','A1_FREQ','SE','LOG10_P']) "
        "WHEN upper(preset) = 'REGENIE' THEN map(['SNP','BP','CHR','A1','A2','BETA','N','N_CAS','N_CON','INFO','FRQ','SE','LP','AC'], ['ID','GENPOS','CHROM','ALLELE1','ALLELE0','BETA','N','N_CASES','N_CONTROLS','INFO','A1FREQ','SE','LOG10P','AC_ALLELE1']) "
        "WHEN upper(preset) = 'SAIGE' THEN map(['SNP','BP','CHR','A1','A2','P','BETA','N','N_CAS','N_CON','INFO','FRQ','SE','AC'], ['SNPID','POS','CHR','Allele2','Allele1','p.value','BETA','N','N_case','N_ctrl','imputationInfo','AF_Allele2','SE','AC_Allele2']) "
        "WHEN upper(preset) = 'BOLT' THEN map(['SNP','BP','CHR','A1','A2','P','BETA','FRQ','SE'], ['SNP','BP','CHR','ALLELE1','ALLELE0','P_LINREG','BETA','A1FREQ','SE']) "
        "WHEN upper(preset) = 'METAL' THEN map(['SNP','BP','CHR','A1','A2','P','Z','BETA','N','FRQ','SE','LP','NEFF','HET_I2','HET_P','HET_LP','DIRE'], ['MarkerName','Position','Chromosome','Allele1','Allele2','P-value','Zscore','Effect','N','Freq1','StdErr','log(P)','Weight','HetISq','HetPVal','logHetP','Direction']) "
        "WHEN upper(preset) = 'PGS' THEN map(['CHR','BP','SNP','A1','A2','OR','BETA','FRQ'], ['chr_name','chr_position','rsID','effect_allele','other_allele','OR','effect_weight','allelefrequency_effect']) "
        "WHEN upper(preset) = 'SSF' OR upper(preset) = 'GWAS-SSF' THEN map(['CHR','BP','SNP','A1','A2','P','OR','BETA','N','INFO','FRQ','SE'], ['chromosome','base_pair_location','variant_id','effect_allele','other_allele','p_value','odds_ratio','beta','n','info','effect_allele_frequency','standard_error']) "
        "ELSE error('duckdb_munge: unknown preset') END")) {
        return false;
    }
    if (!run_sql_or_fail(connection,
        "CREATE OR REPLACE MACRO duckdb_munge_resolved_map(preset, column_map, column_map_file) AS "
        "CASE "
        "WHEN preset <> '' AND map_extract_value(column_map, '') IS NULL THEN error('duckdb_munge: specify only one of preset, column_map, or column_map_file') "
        "WHEN preset <> '' AND column_map_file <> '' THEN error('duckdb_munge: specify only one of preset, column_map, or column_map_file') "
        "WHEN map_extract_value(column_map, '') IS NULL AND column_map_file <> '' THEN error('duckdb_munge: specify only one of preset, column_map, or column_map_file') "
        "WHEN map_extract_value(column_map, '') IS NULL THEN column_map "
        "WHEN column_map_file <> '' THEN error('duckdb_munge: column_map_file is not supported directly in SQL; use column_map or the R wrapper') "
        "WHEN preset <> '' THEN duckdb_munge_preset_map(preset) "
        "ELSE error('duckdb_munge: one of preset, column_map, or column_map_file is required') "
        "END")) {
        return false;
    }
    {
        static const char *const duckdb_munge_sql[] = {
            "CREATE OR REPLACE MACRO duckdb_munge(table_name, preset := '', column_map := map([''], ['']), column_map_file := '', ",
            "fasta_ref := NULL, iffy_tag := 'IFFY', mismatch_tag := 'REF_MISMATCH', ns := NULL, nc := NULL, ne := NULL) AS TABLE ",
            "SELECT mu.chrom, mu.pos, mu.id, mu.ref, mu.alt, mu.alleles_swapped, mu.filter, ",
            "mu.ns, mu.ez, mu.nc, mu.es, mu.se, mu.lp, mu.af, mu.ac, mu.ne ",
            "FROM ",
            "query(",
            "  'SELECT ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'CHR')), 'NULL') || ' AS __chrom, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'BP')), 'NULL') || ' AS __pos, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'A1')), 'NULL') || ' AS __a1, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'A2')), 'NULL') || ' AS __a2, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'SNP')), 'NULL') || ' AS __id, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'P')), 'NULL') || ' AS __p, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'Z')), 'NULL') || ' AS __z, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'OR')), 'NULL') || ' AS __or, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'BETA')), 'NULL') || ' AS __beta, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'N')), 'NULL') || ' AS __n, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'N_CAS')), 'NULL') || ' AS __n_cas, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'N_CON')), 'NULL') || ' AS __n_con, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'INFO')), 'NULL') || ' AS __info, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'FRQ')), 'NULL') || ' AS __frq, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'SE')), 'NULL') || ' AS __se, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'LP')), 'NULL') || ' AS __lp, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'AC')), 'NULL') || ' AS __ac, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'NEFF')), 'NULL') || ' AS __neff, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'NEFFDIV2')), 'NULL') || ' AS __neffdiv2, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'HET_I2')), 'NULL') || ' AS __het_i2, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'HET_P')), 'NULL') || ' AS __het_p, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'HET_LP')), 'NULL') || ' AS __het_lp, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'DIRE')), 'NULL') || ' AS __dire FROM ' || table_name",
            ") src, ",
            "LATERAL (SELECT bcftools_munge_row(",
            "  src.__chrom, src.__pos, src.__a1, src.__a2, src.__id, src.__p, src.__z, src.__or, src.__beta, src.__n, src.__n_cas, src.__n_con, ",
            "  src.__info, src.__frq, src.__se, src.__lp, src.__ac, src.__neff, src.__neffdiv2, src.__het_i2, src.__het_p, src.__het_lp, src.__dire, ",
            "  fasta_ref, iffy_tag, mismatch_tag, ns, nc, ne",
            ") AS mu) q"
        };
        if (!run_sql_parts_or_fail(connection, duckdb_munge_sql, sizeof(duckdb_munge_sql) / sizeof(duckdb_munge_sql[0]))) {
            return false;
        }
    }
    /* Variant with METAL/heterogeneity output columns (SI, I2, CQ, ED) */
    {
        static const char *const duckdb_munge_metal_sql[] = {
            "CREATE OR REPLACE MACRO duckdb_munge_metal(table_name, preset := '', column_map := map([''], ['']), column_map_file := '', ",
            "fasta_ref := NULL, iffy_tag := 'IFFY', mismatch_tag := 'REF_MISMATCH', ns := NULL, nc := NULL, ne := NULL) AS TABLE ",
            "SELECT mu.chrom, mu.pos, mu.id, mu.ref, mu.alt, mu.alleles_swapped, mu.filter, ",
            "mu.ns, mu.ez, mu.nc, mu.es, mu.se, mu.lp, mu.af, mu.ac, mu.ne, mu.si, mu.i2, mu.cq, mu.ed ",
            "FROM ",
            "query(",
            "  'SELECT ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'CHR')), 'NULL') || ' AS __chrom, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'BP')), 'NULL') || ' AS __pos, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'A1')), 'NULL') || ' AS __a1, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'A2')), 'NULL') || ' AS __a2, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'SNP')), 'NULL') || ' AS __id, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'P')), 'NULL') || ' AS __p, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'Z')), 'NULL') || ' AS __z, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'OR')), 'NULL') || ' AS __or, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'BETA')), 'NULL') || ' AS __beta, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'N')), 'NULL') || ' AS __n, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'N_CAS')), 'NULL') || ' AS __n_cas, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'N_CON')), 'NULL') || ' AS __n_con, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'INFO')), 'NULL') || ' AS __info, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'FRQ')), 'NULL') || ' AS __frq, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'SE')), 'NULL') || ' AS __se, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'LP')), 'NULL') || ' AS __lp, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'AC')), 'NULL') || ' AS __ac, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'NEFF')), 'NULL') || ' AS __neff, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'NEFFDIV2')), 'NULL') || ' AS __neffdiv2, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'HET_I2')), 'NULL') || ' AS __het_i2, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'HET_P')), 'NULL') || ' AS __het_p, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'HET_LP')), 'NULL') || ' AS __het_lp, ' || ",
            "  coalesce(duckhts_quote_ident(map_extract_value(duckdb_munge_resolved_map(preset, column_map, column_map_file), 'DIRE')), 'NULL') || ' AS __dire FROM ' || table_name",
            ") src, ",
            "LATERAL (SELECT bcftools_munge_row(",
            "  src.__chrom, src.__pos, src.__a1, src.__a2, src.__id, src.__p, src.__z, src.__or, src.__beta, src.__n, src.__n_cas, src.__n_con, ",
            "  src.__info, src.__frq, src.__se, src.__lp, src.__ac, src.__neff, src.__neffdiv2, src.__het_i2, src.__het_p, src.__het_lp, src.__dire, ",
            "  fasta_ref, iffy_tag, mismatch_tag, ns, nc, ne",
            ") AS mu) q"
        };
        if (!run_sql_parts_or_fail(connection, duckdb_munge_metal_sql,
                                   sizeof(duckdb_munge_metal_sql) / sizeof(duckdb_munge_metal_sql[0]))) {
            return false;
        }
    }
    if (!run_sql_or_fail(connection,
        "CREATE OR REPLACE MACRO read_hts_index_raw(path, format := NULL, index_path := NULL) AS TABLE "
        "SELECT "
        "index_type, index_path, meta AS raw "
        "FROM read_hts_index(path, format := format, index_path := index_path) "
        "WHERE meta IS NOT NULL "
        "LIMIT 1")) {
        return false;
    }
    if (!run_sql_or_fail(connection,
        "CREATE OR REPLACE MACRO duckdb_liftover(table_name, chrom_col, pos_col, ref_col := NULL, alt_col := NULL, "
        "chain_path := NULL, dst_fasta_ref := NULL, src_fasta_ref := NULL, max_snp_gap := 1, max_indel_inc := 250, "
        "lift_mt := false, end_pos_col := NULL, no_left_align := false) AS TABLE "
        "SELECT lo.* "
        "FROM query('SELECT ' || chrom_col || ' AS __duckhts_chrom, ' || pos_col || ' AS __duckhts_pos, ' || "
        "coalesce(ref_col, 'NULL') || ' AS __duckhts_ref, ' || coalesce(alt_col, 'NULL') || ' AS __duckhts_alt, ' || "
        "coalesce(end_pos_col, 'NULL::BIGINT') || ' AS __duckhts_end_pos FROM ' || table_name) src, "
        "LATERAL (SELECT CASE "
        "WHEN src.__duckhts_chrom IS NULL THEN error('bcftools_liftover: chrom must be non-null') "
        "WHEN length(src.__duckhts_chrom) = 0 THEN error('bcftools_liftover: chrom must be non-empty') "
        "WHEN src.__duckhts_pos IS NULL THEN error('bcftools_liftover: pos must be non-null') "
        "WHEN src.__duckhts_pos < 1 THEN error('bcftools_liftover: pos must be >= 1') "
        "ELSE bcftools_liftover(src.__duckhts_chrom, src.__duckhts_pos, src.__duckhts_ref, src.__duckhts_alt, "
        "chain_path, dst_fasta_ref, src_fasta_ref, max_snp_gap, max_indel_inc, lift_mt, src.__duckhts_end_pos, no_left_align) "
        "END AS lo) q")) {
        return false;
    }
    {
        static const char *const duckhts_bcftools_norm_sql[] = {
            "CREATE OR REPLACE MACRO duckhts_bcftools_norm(table_name, fasta_ref, chrom_col := 'chrom', pos_col := 'pos', ref_col := 'ref', alt_col := 'alt', ",
            "split_multiallelic := false, end_pos_col := NULL, svlen_col := NULL, fasta_index_path := NULL, gzi_path := NULL) AS TABLE ",
            "FROM query(CASE WHEN split_multiallelic THEN ",
            "'SELECT src.*, ' || ",
            "'(bcftools_norm_row(src.' || duckhts_quote_ident(chrom_col) || ', src.' || duckhts_quote_ident(pos_col) || ', src.' || duckhts_quote_ident(ref_col) || ', CASE WHEN u0.alt_index IS NULL THEN NULL::VARCHAR[] ELSE list_value(u0.alt_item) END, ' || ",
            "duckhts_quote_string(fasta_ref) || ', ' || coalesce('src.' || duckhts_quote_ident(end_pos_col), 'NULL::BIGINT') || ', ' || coalesce('src.' || duckhts_quote_ident(svlen_col), 'NULL::BIGINT') || ', ' || ",
            "coalesce(duckhts_quote_string(fasta_index_path), 'NULL') || ', ' || coalesce(duckhts_quote_string(gzi_path), 'NULL') || ')).pos_normed AS pos_normed, ' || ",
            "'(bcftools_norm_row(src.' || duckhts_quote_ident(chrom_col) || ', src.' || duckhts_quote_ident(pos_col) || ', src.' || duckhts_quote_ident(ref_col) || ', CASE WHEN u0.alt_index IS NULL THEN NULL::VARCHAR[] ELSE list_value(u0.alt_item) END, ' || ",
            "duckhts_quote_string(fasta_ref) || ', ' || coalesce('src.' || duckhts_quote_ident(end_pos_col), 'NULL::BIGINT') || ', ' || coalesce('src.' || duckhts_quote_ident(svlen_col), 'NULL::BIGINT') || ', ' || ",
            "coalesce(duckhts_quote_string(fasta_index_path), 'NULL') || ', ' || coalesce(duckhts_quote_string(gzi_path), 'NULL') || ')).end_pos_normed AS end_pos_normed, ' || ",
            "'(bcftools_norm_row(src.' || duckhts_quote_ident(chrom_col) || ', src.' || duckhts_quote_ident(pos_col) || ', src.' || duckhts_quote_ident(ref_col) || ', CASE WHEN u0.alt_index IS NULL THEN NULL::VARCHAR[] ELSE list_value(u0.alt_item) END, ' || ",
            "duckhts_quote_string(fasta_ref) || ', ' || coalesce('src.' || duckhts_quote_ident(end_pos_col), 'NULL::BIGINT') || ', ' || coalesce('src.' || duckhts_quote_ident(svlen_col), 'NULL::BIGINT') || ', ' || ",
            "coalesce(duckhts_quote_string(fasta_index_path), 'NULL') || ', ' || coalesce(duckhts_quote_string(gzi_path), 'NULL') || ')).ref_normed AS ref_normed, ' || ",
            "'(bcftools_norm_row(src.' || duckhts_quote_ident(chrom_col) || ', src.' || duckhts_quote_ident(pos_col) || ', src.' || duckhts_quote_ident(ref_col) || ', CASE WHEN u0.alt_index IS NULL THEN NULL::VARCHAR[] ELSE list_value(u0.alt_item) END, ' || ",
            "duckhts_quote_string(fasta_ref) || ', ' || coalesce('src.' || duckhts_quote_ident(end_pos_col), 'NULL::BIGINT') || ', ' || coalesce('src.' || duckhts_quote_ident(svlen_col), 'NULL::BIGINT') || ', ' || ",
            "coalesce(duckhts_quote_string(fasta_index_path), 'NULL') || ', ' || coalesce(duckhts_quote_string(gzi_path), 'NULL') || ')).alt_normed[1] AS alt_normed, ' || ",
            "'(bcftools_norm_row(src.' || duckhts_quote_ident(chrom_col) || ', src.' || duckhts_quote_ident(pos_col) || ', src.' || duckhts_quote_ident(ref_col) || ', CASE WHEN u0.alt_index IS NULL THEN NULL::VARCHAR[] ELSE list_value(u0.alt_item) END, ' || ",
            "duckhts_quote_string(fasta_ref) || ', ' || coalesce('src.' || duckhts_quote_ident(end_pos_col), 'NULL::BIGINT') || ', ' || coalesce('src.' || duckhts_quote_ident(svlen_col), 'NULL::BIGINT') || ', ' || ",
            "coalesce(duckhts_quote_string(fasta_index_path), 'NULL') || ', ' || coalesce(duckhts_quote_string(gzi_path), 'NULL') || ')).normed AS normed, ' || ",
            "'(bcftools_norm_row(src.' || duckhts_quote_ident(chrom_col) || ', src.' || duckhts_quote_ident(pos_col) || ', src.' || duckhts_quote_ident(ref_col) || ', CASE WHEN u0.alt_index IS NULL THEN NULL::VARCHAR[] ELSE list_value(u0.alt_item) END, ' || ",
            "duckhts_quote_string(fasta_ref) || ', ' || coalesce('src.' || duckhts_quote_ident(end_pos_col), 'NULL::BIGINT') || ', ' || coalesce('src.' || duckhts_quote_ident(svlen_col), 'NULL::BIGINT') || ', ' || ",
            "coalesce(duckhts_quote_string(fasta_index_path), 'NULL') || ', ' || coalesce(duckhts_quote_string(gzi_path), 'NULL') || ')).norm_status AS norm_status, u0.alt_index AS alt_index ' || ",
            "'FROM (SELECT * FROM ' || table_name || ') src ' || ",
            "'LEFT JOIN LATERAL UNNEST(duckhts_alt_to_list(src.' || duckhts_quote_ident(alt_col) || ')) WITH ORDINALITY AS u0(alt_item, alt_index) ON TRUE' ",
            "ELSE ",
            "'SELECT src.*, unnest(bcftools_norm_row(src.' || duckhts_quote_ident(chrom_col) || ', src.' || duckhts_quote_ident(pos_col) || ', src.' || duckhts_quote_ident(ref_col) || ', src.' || duckhts_quote_ident(alt_col) || ', ' || ",
            "duckhts_quote_string(fasta_ref) || ', ' || coalesce('src.' || duckhts_quote_ident(end_pos_col), 'NULL::BIGINT') || ', ' || coalesce('src.' || duckhts_quote_ident(svlen_col), 'NULL::BIGINT') || ', ' || ",
            "coalesce(duckhts_quote_string(fasta_index_path), 'NULL') || ', ' || coalesce(duckhts_quote_string(gzi_path), 'NULL') || ')) ' || ",
            "'FROM (SELECT * FROM ' || table_name || ') src' ",
            "END)"
        };
        if (!run_sql_parts_or_fail(connection, duckhts_bcftools_norm_sql,
                                   sizeof(duckhts_bcftools_norm_sql) / sizeof(duckhts_bcftools_norm_sql[0]))) {
            return false;
        }
    }

    return true;
}
