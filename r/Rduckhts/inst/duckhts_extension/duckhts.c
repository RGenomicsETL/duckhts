/**
 * DuckHTS Extension Entry Point
 *
 * Registers all HTS reader table functions with DuckDB.
 */

#define DUCKDB_EXTENSION_NAME duckhts

#include "duckdb_extension.h"

DUCKDB_EXTENSION_EXTERN

/* bcf_reader.c */
extern void register_read_bcf_function(duckdb_connection connection);
/* bam_reader.c */
extern void register_read_bam_function(duckdb_connection connection);
/* seq_reader.c */
extern void register_read_fasta_function(duckdb_connection connection);
extern void register_read_fastq_function(duckdb_connection connection);
/* tabix_reader.c */
extern void register_read_tabix_function(duckdb_connection connection);
extern void register_read_gtf_function(duckdb_connection connection);
extern void register_read_gff_function(duckdb_connection connection);
/* hts_meta_reader.c */
extern void register_read_hts_header_function(duckdb_connection connection);
extern void register_read_hts_index_function(duckdb_connection connection);

static void run_sql_no_fail(duckdb_connection connection, const char *sql) {
    duckdb_result result;
    if (duckdb_query(connection, sql, &result) == DuckDBSuccess) {
        duckdb_destroy_result(&result);
    }
}

DUCKDB_EXTENSION_ENTRYPOINT(duckdb_connection connection,
                            duckdb_extension_info info,
                            struct duckdb_extension_access* access) {
    (void)info;
    (void)access;

    register_read_bcf_function(connection);
    register_read_bam_function(connection);
    register_read_fasta_function(connection);
    register_read_fastq_function(connection);
    register_read_tabix_function(connection);
    register_read_gtf_function(connection);
    register_read_gff_function(connection);
    register_read_hts_header_function(connection);
    register_read_hts_index_function(connection);
    run_sql_no_fail(connection,
        "CREATE OR REPLACE MACRO read_hts_index_spans(path, format := NULL, index_path := NULL) AS TABLE "
        "SELECT "
        "file_format, seqname, tid, "
        "CAST(NULL AS BIGINT) AS bin, "
        "CAST(NULL AS UBIGINT) AS chunk_beg_vo, "
        "CAST(NULL AS UBIGINT) AS chunk_end_vo, "
        "CAST(NULL AS UBIGINT) AS chunk_bytes, "
        "CAST(NULL AS BIGINT) AS seq_start, "
        "length AS seq_end, "
        "mapped, unmapped, n_no_coor, index_type, index_path, meta "
        "FROM read_hts_index(path, format := format, index_path := index_path)");
    run_sql_no_fail(connection,
        "CREATE OR REPLACE MACRO read_hts_index_raw(path, format := NULL, index_path := NULL) AS TABLE "
        "SELECT "
        "index_type, index_path, meta AS raw "
        "FROM read_hts_index(path, format := format, index_path := index_path) "
        "WHERE meta IS NOT NULL "
        "LIMIT 1");

    return true;
}
