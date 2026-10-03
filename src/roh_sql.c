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
#define ROH_ALT_FILTER \
    "WHERE len(list_filter(ALT, lambda a: a NOT IN ('<*>', '<NON_REF>'))) = 1 " \
    "AND ALT[1] NOT IN ('<*>', '<NON_REF>')"
#define ROH_RECORDS \
    "FROM read_bcf(path, tidy_format := true, samples := samples, " \
    "decode_error_policy := 'error') " ROH_ALT_FILTER

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
    /* The long reference (sites x groups rows) is read where it is used rather
     * than held: rescanning the caller's relation is cheaper than keeping
     * millions of rows with string group names live for the whole query. */
    "WITH __dht_ref_rows AS NOT MATERIALIZED (SELECT chromosome AS contig, CAST(position AS BIGINT) AS pos, "
    "CAST(allele_a AS VARCHAR) AS allele_a, CAST(allele_b AS VARCHAR) AS allele_b, "
    "CAST(group_id AS VARCHAR) AS group_id, CAST(frequency AS DOUBLE) AS frequency "
    "FROM query_table(reference_table)), "
    "__dht_prop_rows AS (SELECT CAST(sample_id AS VARCHAR) AS smp, "
    "CAST(group_id AS VARCHAR) AS group_id, CAST(proportion AS DOUBLE) AS proportion "
    "FROM query_table(proportions_table)), "
    /* The reference group set, indexed once in group_id order; every per-site
     * and per-sample list below is ordered by this index. */
    "__dht_groups AS (SELECT group_id, row_number() OVER (ORDER BY group_id) AS g "
    "FROM (SELECT DISTINCT group_id FROM __dht_ref_rows WHERE group_id IS NOT NULL)), "
    /* A sites-only scan finds the reference sites to pivot, so the per-sample
     * calls stream into the lists instead of being held for reuse. */
    "__dht_vcf_sites AS (SELECT DISTINCT CHROM AS chrom, "
    "POS AS pos, REF AS ref, ALT[1] AS alt "
    "FROM read_bcf(path, decode_error_policy := 'error') " ROH_ALT_FILTER "), "
    "__dht_calls AS NOT MATERIALIZED (SELECT CHROM AS chrom, "
    "POS AS pos, REF AS ref, "
    "ALT[1] AS alt, SAMPLE_ID AS smp, " ROH_EVIDENCE " AS evidence " ROH_RECORDS "), ";

/* Contig names are resolved once per distinct name, never per row: the
 * reference's contigs and the VCF's are matched with duckhts_contig_key(), the
 * shared conservative contig key, into a small map, and every row-level join
 * is on raw contig equality. Frequencies therefore carry the VCF's own CHROM.
 * A reference that spells one contig two ways duplicates its rows per site and
 * fails the one-row-per-group check below.
 *
 * Reference alleles are taken to be on the forward strand of the VCF's
 * assembly, as for a FASTA-anchored panel; no strand flip is attempted. Each
 * site, palindromic (A/T, C/G) or not, is oriented by REF: a reference row
 * whose alleles match REF/ALT in either order is used, any other is not.
 *
 * Every join below is on column equality. An OR over the two allele
 * orientations in a join condition is planned as a nested-loop join, so
 * alleles are compared as an unordered pair, and the reference is joined to
 * the calls once per orientation.
 *
 * The matched reference rows are aggregated once per called site (unordered
 * allele pair): frequencies and group indexes in group order, and whether the
 * frequencies are valid and the rows share one orientation. Every per-site
 * check reads that site-level relation instead of the long rows. */
static const char ancestry_validations[] =
    "__dht_contig_map AS (SELECT r.contig, v.chrom FROM "
    "(SELECT DISTINCT contig FROM __dht_ref_rows) AS r JOIN "
    "(SELECT DISTINCT chrom FROM __dht_vcf_sites) AS v "
    "ON duckhts_contig_key(CAST(r.contig AS VARCHAR)) = duckhts_contig_key(v.chrom)), "
    "__dht_ref AS (SELECT c.chrom, r.pos, any_value(r.allele_a) AS allele_a, "
    "any_value(r.allele_b) AS allele_b, count(*) AS n_rows, "
    "list(g.g ORDER BY g.g) AS groups, list(r.frequency ORDER BY g.g) AS frequencies, "
    "bool_and(r.frequency IS NOT NULL AND isfinite(r.frequency) AND r.frequency >= 0 "
    "AND r.frequency <= 1) AS valid_frequencies, "
    "min(r.allele_a) = max(r.allele_a) AS one_orientation "
    "FROM __dht_ref_rows AS r JOIN __dht_contig_map AS m ON m.contig = r.contig "
    "JOIN __dht_vcf_sites AS c ON c.chrom = m.chrom AND c.pos = r.pos "
    "AND least(c.ref, c.alt) = least(r.allele_a, r.allele_b) "
    "AND greatest(c.ref, c.alt) = greatest(r.allele_a, r.allele_b) "
    "LEFT JOIN __dht_groups AS g ON g.group_id = r.group_id "
    "GROUP BY c.chrom, r.pos, least(r.allele_a, r.allele_b), greatest(r.allele_a, r.allele_b)), "
    /* The argument checks come first, so they never depend on which
     * reference rows match a VCF site. */
    "__dht_group_guard AS (SELECT CASE "
    "WHEN af_clamp IS NULL OR NOT isfinite(af_clamp) OR af_clamp < 0 OR af_clamp >= 0.5 "
    "THEN error('af_clamp must be 0 or in (0, 0.5)') "
    "WHEN (SELECT list(group_id ORDER BY group_id) FROM __dht_groups) IS DISTINCT FROM "
    "(SELECT list(group_id ORDER BY group_id) FROM "
    "(SELECT DISTINCT group_id FROM __dht_prop_rows WHERE group_id IS NOT NULL)) "
    "THEN error('reference and proportions group sets must match exactly') "
    /* Exactly one row per group: the sorted group indexes are 1..G, and no
     * row lacks a group (a NULL index sorts last and breaks the equality). */
    "WHEN EXISTS (SELECT 1 FROM __dht_ref WHERE n_rows != (SELECT count(*) FROM __dht_groups) "
    "OR groups IS DISTINCT FROM range(1, (SELECT count(*) FROM __dht_groups) + 1)) "
    "THEN error('reference must have exactly one frequency per site and group') "
    "WHEN EXISTS (SELECT 1 FROM __dht_ref WHERE NOT one_orientation) "
    "THEN error('reference rows of one site must use one allele orientation') "
    "WHEN EXISTS (SELECT 1 FROM __dht_prop_rows GROUP BY smp HAVING count(*) != "
    "(SELECT count(*) FROM __dht_groups) OR count(DISTINCT group_id) != count(*)) "
    "THEN error('each sample must have exactly one proportion per group') "
    "WHEN EXISTS (SELECT 1 FROM __dht_ref WHERE NOT valid_frequencies) OR "
    "EXISTS (SELECT 1 FROM __dht_prop_rows WHERE proportion IS NULL OR "
    "NOT isfinite(proportion) OR proportion < 0 OR proportion > 1) "
    "THEN error('frequencies and proportions must be finite values in [0, 1]') "
    /* Proportions are normalised by their sum: the wrapper rounds each one to
     * seven decimals, and sum_to_one = FALSE allows a sum below one. */
    "WHEN EXISTS (SELECT 1 FROM __dht_prop_rows GROUP BY smp "
    "HAVING NOT (sum(proportion) > 0 AND sum(proportion) <= 1.0 + 1e-6)) "
    "THEN error('ancestry proportions must sum to more than 0 and at most 1 per sample') "
    "ELSE true END AS valid), "
    "__dht_ref_oriented AS (SELECT chrom, pos, allele_a AS ref, allele_b AS alt, "
    "false AS reversed, frequencies FROM __dht_ref UNION ALL "
    "SELECT chrom, pos, allele_b, allele_a, true, frequencies FROM __dht_ref), "
    "__dht_props AS (SELECT p.smp, list(p.proportion / p.total ORDER BY g.g) AS proportions "
    "FROM (SELECT *, sum(proportion) OVER (PARTITION BY smp) AS total FROM __dht_prop_rows) AS p "
    "JOIN __dht_groups AS g ON g.group_id = p.group_id GROUP BY p.smp), "
    "__dht_sites AS NOT MATERIALIZED (SELECT c.chrom, c.pos, c.smp, "
    "CASE WHEN p.proportions IS NULL THEN error('every VCF sample requires ancestry proportions') "
    "WHEN NOT g.valid THEN error('reference and proportions group sets must match exactly') "
    "WHEN r.frequencies IS NULL THEN NULL "
    "WHEN NOT r.reversed THEN greatest(af_clamp, least(1.0 - af_clamp, "
    "list_inner_product(r.frequencies, p.proportions))) "
    "ELSE greatest(af_clamp, least(1.0 - af_clamp, "
    "1.0 - list_inner_product(r.frequencies, p.proportions))) END AS af, c.evidence "
    "FROM __dht_calls AS c CROSS JOIN __dht_group_guard AS g "
    "LEFT JOIN __dht_ref_oriented AS r ON c.chrom = r.chrom AND c.pos = r.pos "
    "AND c.ref = r.ref AND c.alt = r.alt "
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
