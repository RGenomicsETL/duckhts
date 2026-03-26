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
extern void register_fasta_index_function(duckdb_connection connection);
/* interval_udf.c */
extern void register_read_bed_function(duckdb_connection connection);
extern void register_fasta_nuc_function(duckdb_connection connection);
/* bgzip.c */
extern void register_bgzip_function(duckdb_connection connection);
extern void register_bgunzip_function(duckdb_connection connection);
/* hts_index_builder.c */
extern void register_bam_index_function(duckdb_connection connection);
extern void register_bcf_index_function(duckdb_connection connection);
extern void register_tabix_index_function(duckdb_connection connection);
/* liftover_udf.c */
extern void register_liftover_functions(duckdb_connection connection);
/* kmer_udf.c */
extern void register_kmer_udf_functions(duckdb_connection connection);
/* tabix_reader.c */
extern void register_read_tabix_function(duckdb_connection connection);
extern void register_read_gtf_function(duckdb_connection connection);
extern void register_read_gff_function(duckdb_connection connection);
/* hts_meta_reader.c */
extern void register_read_hts_header_function(duckdb_connection connection);
extern void register_read_hts_index_function(duckdb_connection connection);
/* quality_encoding_reader.c */
extern void register_detect_quality_encoding_function(duckdb_connection connection);

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
    register_fasta_index_function(connection);
    register_read_bed_function(connection);
    register_fasta_nuc_function(connection);
    register_bgzip_function(connection);
    register_bgunzip_function(connection);
    register_bam_index_function(connection);
    register_bcf_index_function(connection);
    register_tabix_index_function(connection);
    register_liftover_functions(connection);
    register_kmer_udf_functions(connection);
    register_read_tabix_function(connection);
    register_read_gtf_function(connection);
    register_read_gff_function(connection);
    register_read_hts_header_function(connection);
    register_read_hts_index_function(connection);
    register_detect_quality_encoding_function(connection);
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
    run_sql_no_fail(connection,
        "CREATE OR REPLACE MACRO duckdb_liftover(table_name, chrom_col, pos_col, ref_col := NULL, alt_col := NULL, "
        "chain_path := NULL, dst_fasta_ref := NULL, src_fasta_ref := NULL, max_snp_gap := 1, max_indel_inc := 250) AS TABLE "
        "SELECT lo.* "
        "FROM query('SELECT ' || chrom_col || ' AS __duckhts_chrom, ' || pos_col || ' AS __duckhts_pos, ' || "
        "coalesce(ref_col, 'NULL') || ' AS __duckhts_ref, ' || coalesce(alt_col, 'NULL') || ' AS __duckhts_alt FROM ' || table_name) src, "
        "LATERAL (SELECT bcftools_liftover(src.__duckhts_chrom, src.__duckhts_pos, src.__duckhts_ref, src.__duckhts_alt, "
        "chain_path, dst_fasta_ref, src_fasta_ref, max_snp_gap, max_indel_inc) AS lo) q");

    return true;
}
