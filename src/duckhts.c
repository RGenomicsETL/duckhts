/**
 * DuckHTS Extension Entry Point
 *
 * Registers all HTS reader table functions with DuckDB.
 */

#define DUCKDB_EXTENSION_NAME duckhts

#include "duckdb_extension.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
/* quality_encoding_reader.c */
extern void register_detect_quality_encoding_function(duckdb_connection connection);
/* score_udf.c */
extern void register_bcftools_score_function(duckdb_connection connection);

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
    register_munge_functions(connection);
    register_kmer_udf_functions(connection);
    register_read_tabix_function(connection);
    register_read_gtf_function(connection);
    register_read_gff_function(connection);
    register_read_hts_header_function(connection);
    register_read_hts_index_function(connection);
    register_detect_quality_encoding_function(connection);
    register_bcftools_score_function(connection);
    if (!run_sql_or_fail(connection,
        "CREATE OR REPLACE MACRO duckhts_quote_ident(x) AS "
        "CASE WHEN x IS NULL THEN NULL ELSE '\"' || replace(x, '\"', '\"\"') || '\"' END")) {
        return false;
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
        "FROM read_hts_index(path, format := format, index_path := index_path)")) {
        return false;
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

    return true;
}
