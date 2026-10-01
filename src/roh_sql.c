/* duckhts_roh and duckhts_roh_af_table: runs of homozygosity from a VCF/BCF.
 *
 * The macros read records with read_bcf(tidy_format := true), build one
 * list(... ORDER BY pos) per sample and chromosome, and pass it to the native
 * duckhts_roh_segments kernel. Relations (frequencies, genetic map) are passed
 * by table name through query_table, as in duckhts_cgranges_from_table. A
 * macro cannot take an absent relation, so the genetic map is a separate
 * overload and the frequency source is chosen by macro name.
 */
#include "duckhts_registration.h"

#include <stdbool.h>
#include <stddef.h>

#define ROH_OPTIONS \
    "hw_to_az := 6.7e-8, az_to_hw := 5e-9, gt_error := NULL, rec_rate := NULL, samples := NULL"

/* The evidence column is chosen at bind time so that GT mode never needs
 * FORMAT/PL, and projection pushdown decodes only the columns used. */
#define ROH_EVIDENCE \
    "COLUMNS(lambda c: c = CASE WHEN gt_error IS NULL THEN 'FORMAT_PL' ELSE 'FORMAT_GT' END)"

/* Records as bcftools roh sees them: one ALT other than the symbolic unseen
 * allele, and it is the first ALT. */
#define ROH_RECORDS \
    "FROM read_bcf(path, tidy_format := true, samples := samples, " \
    "decode_error_policy := 'error') " \
    "WHERE len(list_filter(ALT, lambda a: a NOT IN ('<*>', '<NON_REF>'))) = 1 " \
    "AND ALT[1] NOT IN ('<*>', '<NON_REF>')"

static const char sites_from_info_tag[] =
    "WITH __dht_sites AS NOT MATERIALIZED (SELECT CHROM AS chrom, POS AS pos, "
    "SAMPLE_ID AS smp, "
    "CAST(list_extract(COLUMNS(lambda c: c = 'INFO_' || af_tag), 1) AS DOUBLE) AS af, "
    ROH_EVIDENCE " AS evidence " ROH_RECORDS "), ";

static const char sites_from_relation[] =
    "WITH __dht_freq AS (SELECT CAST(chrom AS VARCHAR) AS f_chrom, CAST(pos AS BIGINT) AS f_pos, "
    "CAST(ref AS VARCHAR) AS f_ref, CAST(alt AS VARCHAR) AS f_alt, CAST(af AS DOUBLE) AS f_af "
    "FROM query_table(af_table)), "
    "__dht_calls AS NOT MATERIALIZED (SELECT CHROM AS chrom, POS AS pos, REF AS ref, "
    "array_to_string(ALT, ',') AS alt, SAMPLE_ID AS smp, " ROH_EVIDENCE " AS evidence "
    ROH_RECORDS "), "
    "__dht_sites AS NOT MATERIALIZED (SELECT c.chrom, c.pos, c.smp, f.f_af AS af, c.evidence "
    "FROM __dht_calls AS c LEFT JOIN __dht_freq AS f ON f.f_chrom = c.chrom "
    "AND f.f_pos = c.pos AND f.f_ref = c.ref AND f.f_alt = c.alt), ";

/* One sorted list per sample and chromosome. The columns travel in a single
 * struct list so that all of them share one ordering, including ties. */
static const char site_lists[] =
    "__dht_lists AS (SELECT smp, chrom, list(struct_pack(pos := pos, af := af, "
    "pl := CASE WHEN gt_error IS NULL THEN (CAST(evidence AS INTEGER[]))[1:3] END, "
    "dosage := CASE WHEN gt_error IS NOT NULL AND "
    "regexp_matches(CAST(evidence AS VARCHAR), '^[01][/|][01]$') THEN "
    "CAST(substr(CAST(evidence AS VARCHAR), 1, 1) AS INTEGER) + "
    "CAST(substr(CAST(evidence AS VARCHAR), 3, 1) AS INTEGER) END) ORDER BY pos) AS s "
    "FROM __dht_sites GROUP BY smp, chrom)";

static const char map_lists[] =
    ", __dht_map AS (SELECT CAST(chrom AS VARCHAR) AS m_chrom, "
    "list(struct_pack(pos := CAST(pos AS BIGINT), cm := CAST(cm AS DOUBLE)) ORDER BY pos) AS m "
    "FROM query_table(genetic_map) GROUP BY chrom)";

static const char no_map_arguments[] = "NULL::BIGINT[], NULL::DOUBLE[]";
static const char map_arguments[] =
    "list_transform(m.m, lambda x: x.pos), list_transform(m.m, lambda x: x.cm)";

static const char kernel_pl[] =
    "CASE WHEN gt_error IS NULL THEN duckhts_roh_segments("
    "list_transform(l.s, lambda x: x.pos), list_transform(l.s, lambda x: x.af), "
    "list_transform(l.s, lambda x: x.pl), ";
static const char kernel_gt[] =
    ", CAST(rec_rate AS DOUBLE), CAST(hw_to_az AS DOUBLE), CAST(az_to_hw AS DOUBLE)) ELSE "
    "duckhts_roh_segments(list_transform(l.s, lambda x: x.pos), "
    "list_transform(l.s, lambda x: x.af), list_transform(l.s, lambda x: x.dosage), "
    "CAST(gt_error AS DOUBLE), ";
static const char kernel_end[] =
    ", CAST(rec_rate AS DOUBLE), CAST(hw_to_az AS DOUBLE), CAST(az_to_hw AS DOUBLE)) END";

static const char select_head[] =
    " SELECT seg.smp AS \"sample\", seg.chrom AS chrom, seg.r.start AS start, "
    "seg.r.\"end\" AS \"end\", seg.r.\"end\" - seg.r.start + 1 AS length, "
    "seg.r.n_markers AS n_markers, seg.r.quality AS quality "
    "FROM (SELECT l.smp, l.chrom, unnest(";
static const char from_lists[] = ") AS r FROM __dht_lists AS l) AS seg";
static const char from_lists_with_map[] =
    ") AS r FROM __dht_lists AS l JOIN __dht_map AS m ON m.m_chrom = l.chrom) AS seg";

static const char ancestry_sites[] =
    "WITH __dht_ref_rows AS (SELECT "
    "coalesce(try_cast(regexp_replace(CAST(chromosome AS VARCHAR), '^chr', '') AS INTEGER)::VARCHAR, "
    "CAST(chromosome AS VARCHAR)) AS chrom, CAST(position AS BIGINT) AS pos, "
    "CAST(allele_a AS VARCHAR) AS allele_a, CAST(allele_b AS VARCHAR) AS allele_b, "
    "CAST(group_id AS VARCHAR) AS group_id, CAST(frequency AS DOUBLE) AS frequency "
    "FROM query_table(reference_table)), "
    "__dht_prop_rows AS (SELECT CAST(sample_id AS VARCHAR) AS smp, "
    "CAST(group_id AS VARCHAR) AS group_id, CAST(proportion AS DOUBLE) AS proportion "
    "FROM query_table(proportions_table)), "
    "__dht_calls AS MATERIALIZED (SELECT CHROM AS chrom, POS AS pos, REF AS ref, "
    "ALT[1] AS alt, SAMPLE_ID AS smp, " ROH_EVIDENCE " AS evidence " ROH_RECORDS "), ";

static const char ancestry_validations[] =
    "__dht_ref_matches AS (SELECT r.* FROM __dht_ref_rows AS r JOIN "
    "(SELECT DISTINCT regexp_replace(chrom, '^chr', '') AS chrom, pos, ref, alt "
    "FROM __dht_calls) AS c ON c.chrom = r.chrom AND c.pos = r.pos AND "
    "((c.ref = r.allele_a AND c.alt = r.allele_b) OR "
    "(c.ref = r.allele_b AND c.alt = r.allele_a))), "
    "__dht_called_samples AS (SELECT DISTINCT smp FROM __dht_calls), "
    "__dht_group_guard AS (SELECT CASE "
    "WHEN (SELECT list_sort(list_distinct(list(group_id))) FROM __dht_ref_rows) != "
    "(SELECT list_sort(list_distinct(list(group_id))) FROM __dht_prop_rows) "
    "THEN error('reference and proportions group sets must match exactly') "
    "WHEN EXISTS (SELECT 1 FROM __dht_ref_matches GROUP BY chrom, pos, allele_a, allele_b "
    "HAVING count(*) != (SELECT count(DISTINCT group_id) FROM __dht_ref_rows) OR "
    "count(DISTINCT group_id) != count(*)) "
    "THEN error('reference must have exactly one frequency per site and group') "
    "WHEN EXISTS (SELECT 1 FROM __dht_prop_rows AS p JOIN __dht_called_samples AS s "
    "USING (smp) GROUP BY smp HAVING count(*) != "
    "(SELECT count(DISTINCT group_id) FROM __dht_ref_rows) OR count(DISTINCT group_id) != count(*)) "
    "THEN error('each called sample must have exactly one proportion per group') "
    "WHEN EXISTS (SELECT 1 FROM __dht_ref_matches WHERE frequency IS NULL OR "
    "NOT isfinite(frequency) OR frequency < 0 OR frequency > 1) OR "
    "EXISTS (SELECT 1 FROM __dht_prop_rows AS p JOIN __dht_called_samples AS s "
    "USING (smp) WHERE proportion IS NULL OR NOT isfinite(proportion) OR "
    "proportion < 0 OR proportion > 1) "
    "THEN error('frequencies and proportions must be finite values in [0, 1]') "
    "WHEN EXISTS (SELECT 1 FROM __dht_prop_rows AS p JOIN __dht_called_samples AS s "
    "USING (smp) GROUP BY smp HAVING abs(sum(proportion) - 1.0) > 1e-6) "
    "THEN error('ancestry proportions must sum to one per sample') "
    "ELSE true END AS valid), "
    "__dht_ref AS (SELECT chrom, pos, allele_a, allele_b, "
    "list(frequency ORDER BY group_id) AS frequencies FROM __dht_ref_matches "
    "GROUP BY chrom, pos, allele_a, allele_b), "
    "__dht_props AS (SELECT smp, list(proportion ORDER BY group_id) AS proportions "
    "FROM __dht_prop_rows GROUP BY smp), "
    "__dht_sites AS NOT MATERIALIZED (SELECT c.chrom, c.pos, c.smp, "
    "CASE WHEN p.proportions IS NULL THEN error('every VCF sample requires ancestry proportions') "
    "WHEN NOT g.valid THEN error('reference and proportions group sets must match exactly') "
    "WHEN r.frequencies IS NULL THEN NULL "
    "WHEN af_clamp < 0 OR af_clamp >= 0.5 THEN error('af_clamp must be 0 or in (0, 0.5)') "
    "WHEN c.ref = r.allele_a AND c.alt = r.allele_b THEN "
    "greatest(af_clamp, least(1.0 - af_clamp, "
    "list_inner_product(r.frequencies, p.proportions))) "
    "WHEN c.ref = r.allele_b AND c.alt = r.allele_a THEN "
    "greatest(af_clamp, least(1.0 - af_clamp, "
    "1.0 - list_inner_product(r.frequencies, p.proportions))) END AS af, c.evidence "
    "FROM __dht_calls AS c CROSS JOIN __dht_group_guard AS g "
    "LEFT JOIN __dht_ref AS r ON regexp_replace(c.chrom, '^chr', '') = r.chrom "
    "AND c.pos = r.pos AND ((c.ref = r.allele_a AND c.alt = r.allele_b) OR "
    "(c.ref = r.allele_b AND c.alt = r.allele_a)) "
    "AND NOT ((r.allele_a = 'A' AND r.allele_b = 'T') OR "
    "(r.allele_a = 'T' AND r.allele_b = 'A') OR "
    "(r.allele_a = 'C' AND r.allele_b = 'G') OR "
    "(r.allele_a = 'G' AND r.allele_b = 'C')) "
    "LEFT JOIN __dht_props AS p ON p.smp = c.smp), ";

#define ROH_MAX_PARTS 32

typedef struct {
    const char *name;
    const char *source_parameter;
    const char *sites;
} roh_source_t;

/* Both overloads of one macro: the plain body and the body with a genetic map. */
static bool register_roh_macro(duckhts_registration_t *registration, const roh_source_t *source) {
    const char *parts[ROH_MAX_PARTS];
    size_t count = 0;
    parts[count++] = "CREATE OR REPLACE MACRO ";
    parts[count++] = source->name;
    parts[count++] = "(path, ";
    parts[count++] = source->source_parameter;
    parts[count++] = ", " ROH_OPTIONS ") AS TABLE (";
    parts[count++] = source->sites;
    parts[count++] = site_lists;
    parts[count++] = select_head;
    parts[count++] = kernel_pl;
    parts[count++] = no_map_arguments;
    parts[count++] = kernel_gt;
    parts[count++] = no_map_arguments;
    parts[count++] = kernel_end;
    parts[count++] = from_lists;
    parts[count++] = "), (path, ";
    parts[count++] = source->source_parameter;
    parts[count++] = ", genetic_map VARCHAR, " ROH_OPTIONS ") AS TABLE (";
    parts[count++] = source->sites;
    parts[count++] = site_lists;
    parts[count++] = map_lists;
    parts[count++] = select_head;
    parts[count++] = kernel_pl;
    parts[count++] = map_arguments;
    parts[count++] = kernel_gt;
    parts[count++] = map_arguments;
    parts[count++] = kernel_end;
    parts[count++] = from_lists_with_map;
    parts[count++] = ")";
    return duckhts_register_sql_parts(registration, parts, count);
}

static bool register_roh_ancestry_macro(duckhts_registration_t *registration) {
    const char *parts[ROH_MAX_PARTS];
    size_t count = 0;
    parts[count++] = "CREATE OR REPLACE MACRO duckhts_roh_ancestry(path, reference_table, "
                     "proportions_table, af_clamp := 1e-3, " ROH_OPTIONS ") AS TABLE (";
    parts[count++] = ancestry_sites;
    parts[count++] = ancestry_validations;
    parts[count++] = site_lists;
    parts[count++] = select_head;
    parts[count++] = kernel_pl;
    parts[count++] = no_map_arguments;
    parts[count++] = kernel_gt;
    parts[count++] = no_map_arguments;
    parts[count++] = kernel_end;
    parts[count++] = from_lists;
    parts[count++] = "), (path, reference_table, proportions_table, genetic_map VARCHAR, "
                     "af_clamp := 1e-3, " ROH_OPTIONS ") AS TABLE (";
    parts[count++] = ancestry_sites;
    parts[count++] = ancestry_validations;
    parts[count++] = site_lists;
    parts[count++] = map_lists;
    parts[count++] = select_head;
    parts[count++] = kernel_pl;
    parts[count++] = map_arguments;
    parts[count++] = kernel_gt;
    parts[count++] = map_arguments;
    parts[count++] = kernel_end;
    parts[count++] = from_lists_with_map;
    parts[count++] = ")";
    return duckhts_register_sql_parts(registration, parts, count);
}

bool register_duckhts_roh_sql(duckhts_registration_t *registration) {
    static const roh_source_t from_info_tag = {"duckhts_roh", "af_tag", sites_from_info_tag};
    static const roh_source_t from_relation = {"duckhts_roh_af_table", "af_table",
                                               sites_from_relation};
    return register_roh_macro(registration, &from_info_tag) &&
           register_roh_macro(registration, &from_relation) &&
           register_roh_ancestry_macro(registration);
}
