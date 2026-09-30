/* The one definition of duckhts_somalier_panel_sha256, installed by LOAD and
   by the private panel-reading instance of duckhts_somalier_bam_counts. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "duckhts_somalier.h"
#include "duckhts_somalier_panel_sql.h"

static const char panel_sha256_sql_a[] =
"CREATE OR REPLACE MACRO duckhts_somalier_panel_sha256(panel_table) AS ("
"WITH __dht_panel_rows AS MATERIALIZED ("
"SELECT CAST(assembly AS VARCHAR) AS assembly, "
"CAST(site_index AS UBIGINT) AS site_index, CAST(region AS VARCHAR) AS region, "
"CAST(position AS UBIGINT) AS position, CAST(allele_a AS VARCHAR) AS allele_a, "
"CAST(allele_b AS VARCHAR) AS allele_b FROM query_table(panel_table)"
"), __dht_string_validation AS MATERIALIZED (SELECT CASE "
"WHEN count(*) FILTER (WHERE strlen(assembly) > ";
static const char panel_sha256_sql_b[] =
" OR strlen(region) > ";
static const char panel_sha256_sql_c[] =
") != 0 THEN error('duckhts_somalier_panel_sha256: panel assembly and region "
"must be at most ";
static const char panel_sha256_sql_d[] =
" bytes') ELSE true END AS valid FROM __dht_panel_rows), "
"__dht_bounded_panel_rows AS MATERIALIZED (SELECT p.* FROM __dht_panel_rows p "
"CROSS JOIN __dht_string_validation sv WHERE sv.valid), "
"__dht_validation AS MATERIALIZED (SELECT CASE "
"WHEN count(*) = 0 THEN error('duckhts_somalier_panel_sha256: panel is empty') "
"WHEN count(*) FILTER (WHERE assembly IS NULL OR site_index IS NULL "
"OR region IS NULL OR position IS NULL "
"OR allele_a IS NULL OR allele_b IS NULL) != 0 THEN "
"error('duckhts_somalier_panel_sha256: panel identity fields cannot be NULL') "
"WHEN count(DISTINCT assembly) != 1 OR min(assembly) = '' THEN "
"error('duckhts_somalier_panel_sha256: panel must use one assembly') "
"WHEN count(*) > 100000000 THEN "
"error('duckhts_somalier_panel_sha256: panel exceeds 100000000 sites') "
"WHEN count(DISTINCT site_index) != count(*) OR min(site_index) != 0 "
"OR max(site_index) != count(*) - 1 THEN "
"error('duckhts_somalier_panel_sha256: site_index must be unique and dense from zero') "
"WHEN count(*) FILTER (WHERE region = '' OR position = 0 "
"OR NOT regexp_matches(allele_a, '^[ACGT]$') "
"OR NOT regexp_matches(allele_b, '^[ACGT]$') "
"OR allele_a >= allele_b) != 0 THEN "
"error('duckhts_somalier_panel_sha256: panel sites require a positive "
"position and canonical distinct uppercase biallelic SNP alleles A < B') "
"WHEN coalesce(max(CAST(site_index AS HUGEINT)) FILTER (WHERE region NOT IN "
DUCKHTS_SOMALIER_SEX_REGIONS_SQL "), -1) >= "
"coalesce(min(CAST(site_index AS HUGEINT)) FILTER (WHERE region IN "
DUCKHTS_SOMALIER_SEX_REGIONS_SQL "), 9223372036854775807) THEN "
"error('duckhts_somalier_panel_sha256: autosomal sites must precede X/Y sites "
"in site_index') "
"WHEN count(DISTINCT struct_pack(region := region, pos1 := position)) "
"!= count(*) THEN "
"error('duckhts_somalier_panel_sha256: duplicate physical region and position') "
"ELSE true END AS valid FROM __dht_bounded_panel_rows), "
"__dht_validated_panel_rows AS MATERIALIZED (SELECT p.* "
"FROM __dht_bounded_panel_rows p CROSS JOIN __dht_validation v WHERE v.valid), "
"__dht_row_hashes AS MATERIALIZED (SELECT site_index, "
"sha256(site_index::VARCHAR || ':' || hex(encode(region)) || ':' || "
"position::VARCHAR || ':' || hex(encode(allele_a)) || ':' || "
"hex(encode(allele_b))) AS row_hash FROM __dht_validated_panel_rows), "
"__dht_block_hashes AS MATERIALIZED (SELECT site_index // 4096 AS block_index, "
"sha256(string_agg(row_hash, '' ORDER BY site_index)) AS block_hash "
"FROM __dht_row_hashes GROUP BY block_index), "
"__dht_panel_summary AS MATERIALIZED (SELECT min(assembly) AS assembly, "
"count(*) AS site_count FROM __dht_validated_panel_rows) "
"SELECT sha256('duckhts-somalier-panel-v3;biallelic-snp;autosomal-then-sex-chromosomes;' || "
"hex(encode(ps.assembly)) || ';' || ps.site_count::VARCHAR || ';' || "
"string_agg(b.block_hash, '' ORDER BY b.block_index)) "
"FROM __dht_block_hashes b CROSS JOIN __dht_panel_summary ps "
"GROUP BY ps.assembly, ps.site_count)";

char *duckhts_somalier_panel_sha256_macro_sql(void) {
    char max_identity_bytes[16];
    const char *const parts[] = {
        panel_sha256_sql_a, max_identity_bytes,
        panel_sha256_sql_b, max_identity_bytes,
        panel_sha256_sql_c, max_identity_bytes,
        panel_sha256_sql_d
    };
    size_t length = 0;
    char *sql;
    size_t offset = 0;

    snprintf(max_identity_bytes, sizeof(max_identity_bytes), "%u",
             (unsigned int)DUCKHTS_SOMALIER_MAX_IDENTITY_BYTES);
    for (size_t i = 0; i < sizeof(parts) / sizeof(parts[0]); i++) {
        length += strlen(parts[i]);
    }
    sql = malloc(length + 1);
    if (sql == NULL) return NULL;
    for (size_t i = 0; i < sizeof(parts) / sizeof(parts[0]); i++) {
        size_t part_length = strlen(parts[i]);
        memcpy(sql + offset, parts[i], part_length);
        offset += part_length;
    }
    sql[offset] = '\0';
    return sql;
}
