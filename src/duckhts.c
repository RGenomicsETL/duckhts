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

    return true;
}
