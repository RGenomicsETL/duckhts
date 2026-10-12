/* duckhts_roh, duckhts_roh_af_table and duckhts_roh_ancestry: runs of homozygosity
 * from a VCF/BCF; duckhts_roh_counts: from read counts in a caller relation.
 *
 * The macros read records with read_bcf(tidy_format := true), collect the
 * sites of each sample and chromosome with the native __duckhts_roh_sites
 * aggregate, and decode each packed site list with __duckhts_roh_decode.
 * Relations (frequencies, genetic map) are passed by table name through
 * query_table, as in duckhts_cgranges_from_table. A macro cannot take an
 * absent relation, so the genetic map is a separate overload and the frequency
 * source is chosen by macro name.
 *
 * Memory. DuckDB manages and can spill the scans, joins and reference stages.
 * The site lists are native buffers that DuckHTS bounds: max_sites per sample
 * and chromosome, and max_site_bytes for all the lists held at once, each with
 * an explicit error (roh_sites.c).
 */
#include "count_sql.h"
#include "duckhts_registration.h"
#include "roh_sites.h"

#include <stdbool.h>
#include <stddef.h>

/* The SQL below passes the evidence kinds as literals. */
_Static_assert(DUCKHTS_ROH_SITES_PL == 1 && DUCKHTS_ROH_SITES_GT == 2 &&
               DUCKHTS_ROH_SITES_COUNTS == 3, "the ROH macro SQL names these kinds 1, 2 and 3");
_Static_assert(DUCKHTS_ROH_MAX_SITES >= 20000000, "the default max_sites is within the hard limit");

#define ROH_LIMITS "max_sites := 20000000, max_site_bytes := 4294967296"
#define ROH_OPTIONS \
    "hw_to_az := 6.7e-8, az_to_hw := 5e-9, gt_error := NULL, rec_rate := NULL, samples := NULL, " \
    ROH_LIMITS
#define ROH_LIMIT_ARGUMENTS "CAST(max_sites AS BIGINT), CAST(max_site_bytes AS BIGINT)"
/* PL evidence unless gt_error is given. */
#define ROH_KIND "CASE WHEN gt_error IS NULL THEN 1 ELSE 2 END"

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
#define ROH_SYMBOLIC "('<*>', '<NON_REF>')"

static const char sites_from_info_tag[] =
    "WITH __dht_sites AS NOT MATERIALIZED (SELECT CHROM AS chrom, POS AS pos, "
    "SAMPLE_ID AS smp, "
    "CAST(list_extract(COLUMNS(lambda c: c = 'INFO_' || af_tag), 1) AS DOUBLE) AS af, "
    ROH_EVIDENCE " AS evidence " ROH_RECORDS "), ";

/* duckhts_htslib_is_indel(ref, alts): whether HTSlib's bcf_get_variant_types()
 * sets VCF_INDEL for a record with these alleles, a port of
 * bcf_set_variant_type() (htslib vcf.c:5361-5428 in 1.23, unchanged in 1.24).
 * An allele is an indel when, after the common prefix (compared without case),
 * one allele is used up and the other is not and does not end in a breakend
 * bracket; or when both remain and, after the common suffix, one shrinks to
 * its first base and the base where the other stops equals it. '*', symbolic
 * and breakend alleles and a single-base pair are never indels, so an MNP such
 * as AC>GT and a complex substitution such as ACG>AT are not. The record is an
 * indel when any ALT is. ind_p is the common prefix length; ind_s the common
 * suffix length of the rest, leaving at least one base of each. */
static const char htslib_is_indel_macro[] =
    "CREATE OR REPLACE MACRO duckhts_htslib_is_indel(ref, alts) AS "
    "coalesce(list_bool_or(list_transform(alts, lambda ind_x: "
    "CASE WHEN ind_x = '*' OR starts_with(ind_x, '<') OR starts_with(ind_x, '[') "
    "OR starts_with(ind_x, ']') OR (length(ref) = 1 AND length(ind_x) = 1) THEN false "
    "ELSE list_transform([struct_pack(r := upper(ref), a := upper(ind_x))], lambda ind_v: "
    "list_transform([coalesce(list_position(list_transform("
    "range(1, least(length(ind_v.r), length(ind_v.a)) + 1), "
    "lambda ind_i: substr(ind_v.r, ind_i, 1) != substr(ind_v.a, ind_i, 1)), true) - 1, "
    "least(length(ind_v.r), length(ind_v.a)))], lambda ind_p: "
    "CASE WHEN ind_p = length(ind_v.r) AND ind_p = length(ind_v.a) THEN false "
    "WHEN ind_p = length(ind_v.r) THEN NOT (ends_with(ind_v.a, ']') OR ends_with(ind_v.a, '[')) "
    "WHEN ind_p = length(ind_v.a) THEN true "
    "WHEN ends_with(ind_v.a, ']') OR ends_with(ind_v.a, '[') THEN false "
    "ELSE list_transform([coalesce(list_position(list_transform("
    "range(0, least(length(ind_v.r), length(ind_v.a)) - ind_p - 1), "
    "lambda ind_j: substr(ind_v.r, length(ind_v.r) - ind_j, 1) != "
    "substr(ind_v.a, length(ind_v.a) - ind_j, 1)), true) - 1, "
    "least(length(ind_v.r), length(ind_v.a)) - ind_p - 1)], lambda ind_s: "
    "CASE WHEN ind_s = length(ind_v.a) - ind_p - 1 AND ind_s = length(ind_v.r) - ind_p - 1 THEN false "
    "WHEN ind_s = length(ind_v.a) - ind_p - 1 "
    "THEN substr(ind_v.r, length(ind_v.r) - ind_s, 1) = substr(ind_v.a, ind_p + 1, 1) "
    "WHEN ind_s = length(ind_v.r) - ind_p - 1 "
    "THEN substr(ind_v.r, ind_p + 1, 1) = substr(ind_v.a, length(ind_v.a) - ind_s, 1) "
    "ELSE false END)[1] END)[1])[1] END)), false)";

/* site_rule chooses the record that gives a listed position its site. Both
 * arms are written out; the macro argument is a constant, so DuckDB removes
 * the arm whose site_rule test is false before running it.
 *
 * 'matching': a record without a matching frequency row is dropped before the
 * sites are collected, so it never claims a repeated position. With a left
 * join it would arrive before or after the matching record of its position as
 * the hash join happens to emit them, and the record kept would depend on the
 * plan.
 *
 * 'bcftools': bcftools roh --AF-file (bcftools 1.23 vcfroh.c and htslib
 * synced_bcf_reader.c). The frequency file is a target list with alleles
 * (vcfroh.c:1272-1276), and bcftools reads only the first line of a position
 * (synced_bcf_reader.c:1525-1551), so a relation must give one row per
 * position (the guard below). At a listed position, in file order (FILE_OFFSET):
 * a record whose indel class differs from the frequency row's is passed over
 * while another record of the position follows; the last record is always
 * seen (synced_bcf_reader.c:865-886). The row is an indel when any of its
 * alleles is longer than one character (1466-1518). A record with no
 * non-symbolic ALT, or several, claims nothing (vcfroh.c:1010-1027); the first
 * remaining record claims the position (1036-1065), and the position is
 * skipped unless that record's alleles equal the row's (read_AF, 661-686). A
 * claiming record whose first ALT is symbolic is skipped as well: bcftools
 * would decode its second allele, and the evidence here is that of the first
 * ALT. File order needs BGZF input. */
static const char sites_from_relation[] =
    "WITH __dht_freq AS NOT MATERIALIZED (SELECT CAST(chrom AS VARCHAR) AS f_chrom, "
    "CAST(pos AS BIGINT) AS f_pos, CAST(ref AS VARCHAR) AS f_ref, CAST(alt AS VARCHAR) AS f_alt, "
    "CAST(af AS DOUBLE) AS f_af FROM query_table(af_table)), "
    "__dht_calls AS NOT MATERIALIZED (SELECT CHROM AS chrom, POS AS pos, REF AS ref, "
    "array_to_string(ALT, ',') AS alt, FILE_OFFSET AS file_offset, SAMPLE_ID AS smp, "
    ROH_EVIDENCE " AS evidence FROM read_bcf(path, tidy_format := true, samples := samples, "
    "decode_error_policy := 'error', include_file_offset := true) " ROH_ALT_FILTER "), "
    "__dht_bcftools_targets AS NOT MATERIALIZED (SELECT f_chrom AS t_chrom, f_pos AS t_pos, "
    "any_value(f_ref) AS t_ref, any_value(f_alt) AS t_alt, any_value(f_af) AS t_af, "
    "coalesce(length(any_value(f_ref)) > 1 OR "
    "list_bool_or(list_transform(string_split(any_value(f_alt), ','), lambda a: length(a) > 1)), "
    "false) AS t_indel "
    "FROM __dht_freq WHERE site_rule = 'bcftools' AND f_chrom IS NOT NULL AND f_pos IS NOT NULL "
    "GROUP BY f_chrom, f_pos), "
    /* Every record at a listed position, from a scan without samples. */
    "__dht_bcftools_records AS NOT MATERIALIZED (SELECT r.CHROM AS chrom, r.POS AS pos, "
    "r.REF AS ref, r.ALT AS alts, CASE WHEN r.FILE_OFFSET IS NULL THEN "
    "error('site_rule := ''bcftools'' needs file order: use a BGZF-compressed VCF or BCF') "
    "ELSE r.FILE_OFFSET END AS file_offset, "
    "max(r.FILE_OFFSET) OVER (PARTITION BY r.CHROM, r.POS) AS last_offset, "
    "t.t_ref, t.t_alt, t.t_af, t.t_indel "
    "FROM read_bcf(path, samples := '', decode_error_policy := 'error', include_file_offset := true) AS r "
    "JOIN __dht_bcftools_targets AS t ON t.t_chrom = r.CHROM AND t.t_pos = r.POS), "
    "__dht_bcftools_claims AS NOT MATERIALIZED (SELECT chrom, pos, file_offset, t_af, "
    "coalesce(ref = t_ref AND array_to_string(alts, ',') = t_alt AND alts[1] NOT IN " ROH_SYMBOLIC ", "
    "false) AS usable FROM __dht_bcftools_records "
    "WHERE (file_offset = last_offset OR duckhts_htslib_is_indel(ref, alts) = t_indel) "
    "AND len(list_filter(alts, lambda a: a NOT IN " ROH_SYMBOLIC ")) = 1 "
    "QUALIFY row_number() OVER (PARTITION BY chrom, pos ORDER BY file_offset) = 1), "
    "__dht_sites AS NOT MATERIALIZED (SELECT c.chrom, c.pos, c.smp, f.f_af AS af, c.evidence "
    "FROM __dht_calls AS c JOIN __dht_freq AS f ON f.f_chrom = c.chrom "
    "AND f.f_pos = c.pos AND f.f_ref = c.ref AND f.f_alt = c.alt WHERE site_rule = 'matching' "
    "UNION ALL SELECT c.chrom, c.pos, c.smp, k.t_af, c.evidence "
    "FROM __dht_calls AS c JOIN __dht_bcftools_claims AS k ON k.chrom = c.chrom "
    "AND k.pos = c.pos AND k.file_offset = c.file_offset "
    "WHERE k.usable AND site_rule = 'bcftools'), ";

/* One packed site list per sample and chromosome. The evidence column is
 * FORMAT_PL or FORMAT_GT, chosen at bind time, and only the branch of that
 * kind is evaluated: the RR, RA and AA likelihoods (a NULL where FORMAT/PL has
 * fewer than three values), or the dosage of a diploid biallelic call. */
static const char site_lists[] =
    "__dht_lists AS (SELECT smp, chrom, __duckhts_roh_sites(pos, af, "
    "CASE WHEN gt_error IS NULL THEN (CAST(evidence AS INTEGER[]))[1] "
    "WHEN regexp_matches(CAST(evidence AS VARCHAR), '^[01][/|][01]$') THEN "
    "CAST(substr(CAST(evidence AS VARCHAR), 1, 1) AS INTEGER) + "
    "CAST(substr(CAST(evidence AS VARCHAR), 3, 1) AS INTEGER) END, "
    "CASE WHEN gt_error IS NULL THEN (CAST(evidence AS INTEGER[]))[2] END, "
    "CASE WHEN gt_error IS NULL THEN (CAST(evidence AS INTEGER[]))[3] END, "
    ROH_KIND ", " ROH_LIMIT_ARGUMENTS ") AS s "
    "FROM __dht_sites GROUP BY smp, chrom)";

static const char map_lists[] =
    ", __dht_map AS (SELECT CAST(chrom AS VARCHAR) AS m_chrom, "
    "list(struct_pack(pos := CAST(pos AS BIGINT), cm := CAST(cm AS DOUBLE)) ORDER BY pos) AS m "
    "FROM query_table(genetic_map) GROUP BY chrom)";

static const char no_map_arguments[] = "NULL::BIGINT[], NULL::DOUBLE[]";
static const char map_arguments[] =
    "list_transform(m.m, lambda x: x.pos), list_transform(m.m, lambda x: x.cm)";

/* __duckhts_roh_decode(sites, gt_error, seq_error, contamination, map_pos,
 * map_cm, rec_rate, hw_to_az, az_to_hw): the model parameters of the evidence
 * kind, then the map arguments, then the transitions. */
static const char select_head[] =
    " SELECT seg.smp AS \"sample\", seg.chrom AS chrom, seg.r.start AS start, "
    "seg.r.\"end\" AS \"end\", seg.r.\"end\" - seg.r.start + 1 AS length, "
    "seg.r.n_markers AS n_markers, seg.r.quality AS quality "
    "FROM (SELECT l.smp, l.chrom, unnest(__duckhts_roh_decode(l.s, ";
static const char genotype_model[] = "CAST(gt_error AS DOUBLE), NULL::DOUBLE, NULL::DOUBLE, ";
static const char counts_model[] =
    "NULL::DOUBLE, CAST(seq_error AS DOUBLE), CAST(contamination AS DOUBLE), ";
static const char decode_end[] =
    ", CAST(rec_rate AS DOUBLE), CAST(hw_to_az AS DOUBLE), CAST(az_to_hw AS DOUBLE))";
static const char from_lists[] = ") AS r FROM __dht_lists AS l) AS seg";
static const char from_lists_with_map[] =
    ") AS r FROM __dht_lists AS l JOIN __dht_map AS m ON m.m_chrom = l.chrom) AS seg";

/* A second arm that returns no row but always runs, so that an invalid model
 * parameter or limit fails even when there is no site to decode.
 * __duckhts_roh_valid_args() applies the checks of the decode itself. */
#define ROH_GUARD_COLUMNS \
    " UNION ALL SELECT NULL::VARCHAR, NULL::VARCHAR, NULL::BIGINT, NULL::BIGINT, NULL::BIGINT, " \
    "NULL::INTEGER, NULL::DOUBLE "
#define ROH_GENOTYPE_ARGUMENTS_VALID \
    "__duckhts_roh_valid_args(" ROH_KIND ", CAST(gt_error AS DOUBLE), NULL::DOUBLE, NULL::DOUBLE, " \
    "CAST(rec_rate AS DOUBLE), CAST(hw_to_az AS DOUBLE), CAST(az_to_hw AS DOUBLE), " \
    ROH_LIMIT_ARGUMENTS ")"
static const char genotype_guard[] = ROH_GUARD_COLUMNS "WHERE NOT " ROH_GENOTYPE_ARGUMENTS_VALID;
/* duckhts_roh_af_table also checks site_rule, and under 'bcftools' that the
 * relation gives at most one row per position, naming the first such position. */
static const char af_table_guard[] =
    ROH_GUARD_COLUMNS "WHERE CASE WHEN site_rule IN ('matching', 'bcftools') THEN false "
    "ELSE error('site_rule must be ''matching'' or ''bcftools''') IS NULL END"
    ROH_GUARD_COLUMNS "FROM (SELECT f_chrom, f_pos, count(*) AS n FROM __dht_freq "
    "WHERE site_rule = 'bcftools' AND f_chrom IS NOT NULL AND f_pos IS NOT NULL "
    "GROUP BY f_chrom, f_pos HAVING count(*) > 1 ORDER BY f_chrom, f_pos LIMIT 1) AS d "
    "WHERE error('site_rule := ''bcftools'' needs one af_table row per position, but ' || "
    "d.f_chrom || ':' || d.f_pos || ' has ' || d.n || ' rows; bcftools reads the first line "
    "of a position, and a relation has no line order') IS NULL";
/* The ancestry macro also runs its relation checks (__dht_group_guard) here,
 * so they do not depend on a call matching a reference site. */
static const char ancestry_guard[] =
    ROH_GUARD_COLUMNS "FROM __dht_group_guard AS g WHERE NOT (g.valid AND "
    ROH_GENOTYPE_ARGUMENTS_VALID ")";
#define ROH_COUNTS_ARGUMENTS_VALID \
    "__duckhts_roh_valid_args(3, NULL::DOUBLE, CAST(seq_error AS DOUBLE), " \
    "CAST(contamination AS DOUBLE), CAST(rec_rate AS DOUBLE), CAST(hw_to_az AS DOUBLE), " \
    "CAST(az_to_hw AS DOUBLE), " ROH_LIMIT_ARGUMENTS ")"
static const char counts_guard[] = ROH_GUARD_COLUMNS "WHERE NOT " ROH_COUNTS_ARGUMENTS_VALID;

/* One fitted parameter of the sample `sample`, read from the relation `from`
 * that carries __dht_fit's columns: an error when the sample has no fit row,
 * several, or a NULL value. */
#define ROH_FIT_VALUE(from, column, sample) \
    "CASE WHEN " from ".fit_rows IS NULL THEN error('missing fit row for sample ' || coalesce(" sample ", 'NULL')) " \
    "WHEN " from ".fit_rows > 1 THEN error('duplicate fit rows for sample ' || coalesce(" sample ", 'NULL')) " \
    "WHEN " from "." column " IS NULL THEN error('NULL " column " for sample ' || coalesce(" sample ", 'NULL')) " \
    "ELSE " from "." column " END"
static const char fitted_counts_model[] =
    "NULL::DOUBLE, " ROH_FIT_VALUE("l", "seq_error", "l.smp") ", " ROH_FIT_VALUE("l", "contamination", "l.smp") ", ";
/* Every sample of the counts relation is checked against its fit, also one
 * whose chromosomes the genetic map drops from the decode; the second arm
 * checks the other arguments when there is no sample. */
static const char fitted_counts_guard[] =
    ROH_GUARD_COLUMNS "FROM (SELECT DISTINCT CAST(sample_id AS VARCHAR) AS smp FROM query_table(counts_table)) AS s "
    "LEFT JOIN __dht_fit AS fit ON fit.smp IS NOT DISTINCT FROM s.smp "
    "WHERE NOT __duckhts_roh_valid_args(3, NULL::DOUBLE, " ROH_FIT_VALUE("fit", "seq_error", "s.smp") ", "
    ROH_FIT_VALUE("fit", "contamination", "s.smp") ", CAST(rec_rate AS DOUBLE), CAST(hw_to_az AS DOUBLE), "
    "CAST(az_to_hw AS DOUBLE), " ROH_LIMIT_ARGUMENTS ")"
    ROH_GUARD_COLUMNS "FROM (SELECT 0.001 AS seq_error, 0.0 AS contamination) AS defaults "
    "WHERE NOT " ROH_COUNTS_ARGUMENTS_VALID;

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
 * check reads that site-level relation instead of the long rows.
 *
 * The (group index, frequency) pairs of a site are collected unordered and
 * sorted per site. list(... ORDER BY ...) is not used: DuckDB's ordered
 * aggregate keeps a 2048-row buffer for every group with more than 16 rows,
 * which for the chr20 reference (111,058 sites, 21 groups) costs 903 MB where
 * this form costs 159 MB. */
static const char ancestry_validations[] =
    "__dht_contig_map AS (SELECT r.contig, v.chrom FROM "
    "(SELECT DISTINCT contig FROM __dht_ref_rows) AS r JOIN "
    "(SELECT DISTINCT chrom FROM __dht_vcf_sites) AS v "
    "ON duckhts_contig_key(CAST(r.contig AS VARCHAR)) = duckhts_contig_key(v.chrom)), "
    "__dht_ref AS (SELECT chrom, pos, allele_a, allele_b, n_rows, "
    "list_transform(pairs, lambda x: x.g) AS groups, "
    "list_transform(pairs, lambda x: x.f) AS frequencies, valid_frequencies, one_orientation "
    "FROM (SELECT c.chrom, r.pos, any_value(r.allele_a) AS allele_a, "
    "any_value(r.allele_b) AS allele_b, count(*) AS n_rows, "
    "list_sort(list(struct_pack(g := g.g, f := r.frequency))) AS pairs, "
    "bool_and(r.frequency IS NOT NULL AND isfinite(r.frequency) AND r.frequency >= 0 "
    "AND r.frequency <= 1) AS valid_frequencies, "
    "min(r.allele_a) = max(r.allele_a) AS one_orientation "
    "FROM __dht_ref_rows AS r JOIN __dht_contig_map AS m ON m.contig = r.contig "
    "JOIN __dht_vcf_sites AS c ON c.chrom = m.chrom AND c.pos = r.pos "
    "AND least(c.ref, c.alt) = least(r.allele_a, r.allele_b) "
    "AND greatest(c.ref, c.alt) = greatest(r.allele_a, r.allele_b) "
    "LEFT JOIN __dht_groups AS g ON g.group_id = r.group_id "
    "GROUP BY c.chrom, r.pos, least(r.allele_a, r.allele_b), "
    "greatest(r.allele_a, r.allele_b))), "
    /* The argument checks come first, so they never depend on which
     * reference rows match a VCF site. */
    "__dht_group_guard AS MATERIALIZED (SELECT CASE "
    "WHEN af_clamp IS NULL OR NOT isfinite(af_clamp) OR af_clamp < 0 OR af_clamp >= 0.5 "
    "THEN error('af_clamp must be 0 or in (0, 0.5)') "
    /* duckhts_contig_key() is not injective (1 and chr1 share a key), so each key
     * must name one contig on each side, or one contig would be decoded twice. */
    "WHEN EXISTS (SELECT 1 FROM __dht_contig_map GROUP BY chrom HAVING count(*) > 1) "
    "THEN error('two reference contig names match one VCF contig; use one name per contig') "
    "WHEN EXISTS (SELECT 1 FROM __dht_contig_map GROUP BY contig HAVING count(*) > 1) "
    "THEN error('two VCF contig names match one reference contig; use one name per contig') "
    "WHEN (SELECT list(group_id ORDER BY group_id) FROM __dht_groups) IS DISTINCT FROM "
    "(SELECT list(group_id ORDER BY group_id) FROM "
    "(SELECT DISTINCT group_id FROM __dht_prop_rows WHERE group_id IS NOT NULL)) "
    "THEN error('reference and proportions group sets must match exactly') "
    /* Exactly one row per group: the sorted group indexes are 1..G, and no
     * row lacks a group (a NULL index breaks the equality wherever it sorts). */
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
    /* The selected samples come from the header, so this does not depend on
     * a call matching a reference site. */
    "WHEN EXISTS (SELECT 1 FROM read_bcf_samples(path, samples := samples) AS s "
    "WHERE s.sample_name NOT IN (SELECT smp FROM __dht_prop_rows WHERE smp IS NOT NULL)) "
    "THEN error('every VCF sample requires ancestry proportions') "
    "ELSE true END AS valid), ";

/* Per-sample frequencies of the called sites: the proportion-weighted sum of
 * the group frequencies, for REF = allele_a, or one minus it when reversed. */
static const char ancestry_frequencies[] =
    "__dht_ref_oriented AS (SELECT chrom, pos, allele_a AS ref, allele_b AS alt, "
    "false AS reversed, frequencies FROM __dht_ref UNION ALL "
    "SELECT chrom, pos, allele_b, allele_a, true, frequencies FROM __dht_ref), "
    "__dht_props AS (SELECT p.smp, list(p.proportion / p.total ORDER BY g.g) AS proportions "
    "FROM (SELECT *, sum(proportion) OVER (PARTITION BY smp) AS total FROM __dht_prop_rows) AS p "
    "JOIN __dht_groups AS g ON g.group_id = p.group_id GROUP BY p.smp), "
    "__dht_sites AS NOT MATERIALIZED (SELECT c.chrom, c.pos, c.smp, "
    "CASE WHEN NOT g.valid THEN error('reference and proportions group sets must match exactly') "
    "WHEN NOT r.reversed THEN greatest(af_clamp, least(1.0 - af_clamp, "
    "list_inner_product(r.frequencies, p.proportions))) "
    "ELSE greatest(af_clamp, least(1.0 - af_clamp, "
    "1.0 - list_inner_product(r.frequencies, p.proportions))) END AS af, c.evidence "
    /* As in duckhts_roh_af_table, a call without a reference site is dropped
     * here, so it never claims a repeated position. */
    "FROM __dht_calls AS c CROSS JOIN __dht_group_guard AS g "
    "JOIN __dht_ref_oriented AS r ON c.chrom = r.chrom AND c.pos = r.pos "
    "AND c.ref = r.ref AND c.alt = r.alt "
    "JOIN __dht_props AS p ON p.smp = c.smp), ";

/* Read counts from a caller relation: one row per sample and site, with the
 * frequency of the counted allele. Any frequency source composes in SQL
 * before the macro, for example counts at panel sites joined to INFO/AF, an
 * af_table or ancestry-weighted frequencies. */
#define ROH_COUNTS_FITTED_OPTIONS \
    "hw_to_az := 6.7e-8, az_to_hw := 5e-9, rec_rate := NULL, " ROH_LIMITS
#define ROH_COUNTS_OPTIONS \
    "seq_error := 1e-3, contamination := 0.0, " ROH_COUNTS_FITTED_OPTIONS

#define ROH_COUNTS_INPUT \
    "WITH __dht_counts AS (SELECT CAST(sample_id AS VARCHAR) AS smp, " \
    "CAST(chrom AS VARCHAR) AS chrom, CAST(pos AS BIGINT) AS pos, " \
    DUCKHTS_WHOLE_COUNT("ref_count") " AS ref_count, " DUCKHTS_WHOLE_COUNT("alt_count") " AS alt_count, " \
    "CAST(af AS DOUBLE) AS af FROM query_table(counts_table)), "

static const char counts_lists[] =
    ROH_COUNTS_INPUT
    "__dht_lists AS (SELECT smp, chrom, __duckhts_roh_sites(pos, af, ref_count, alt_count, "
    "NULL::INTEGER, 3, " ROH_LIMIT_ARGUMENTS ") AS s "
    "FROM __dht_counts GROUP BY smp, chrom)";

/* One row per sample of the fit relation, with its number of rows, so that a
 * duplicate is an error whatever DuckDB's scalar-subquery setting, and a join
 * cannot repeat a site. A NULL sample is a sample, as in the scalar decode and
 * the count-error fit, so the joins are NULL-safe. The fit joins the counts
 * rows before they are collected and its values become keys of the site
 * lists: joined after collection, DuckDB built its hash table on the site
 * lists, a second copy of every site. */
static const char fitted_counts_lists[] =
    ROH_COUNTS_INPUT
    "__dht_fit AS (SELECT CAST(sample AS VARCHAR) AS smp, CAST(seq_error AS DOUBLE) AS seq_error, "
    "CAST(contamination AS DOUBLE) AS contamination, "
    "count(*) OVER (PARTITION BY CAST(sample AS VARCHAR)) AS fit_rows FROM query_table(fit_table) "
    "QUALIFY row_number() OVER (PARTITION BY CAST(sample AS VARCHAR)) = 1), "
    "__dht_lists AS (SELECT c.smp, c.chrom, fit.fit_rows, fit.seq_error, fit.contamination, "
    "__duckhts_roh_sites(c.pos, c.af, c.ref_count, c.alt_count, NULL::INTEGER, 3, " ROH_LIMIT_ARGUMENTS ") AS s "
    "FROM __dht_counts AS c LEFT JOIN __dht_fit AS fit ON fit.smp IS NOT DISTINCT FROM c.smp "
    "GROUP BY c.smp, c.chrom, fit.fit_rows, fit.seq_error, fit.contamination)";

#define ROH_MAX_PARTS 32

typedef struct {
    const char *name;
    const char *source_parameter;
    const char *options;      /* the named parameters after the source, with their defaults */
    const char *sites;
    const char *source_guard; /* checks of the source's own options; may be empty */
} roh_source_t;

/* Both overloads of one macro: the plain body and the body with a genetic map. */
static bool register_roh_macro(duckhts_registration_t *registration, const roh_source_t *source) {
    const char *parts[ROH_MAX_PARTS];
    size_t count = 0;
    parts[count++] = "CREATE OR REPLACE MACRO ";
    parts[count++] = source->name;
    parts[count++] = "(path, ";
    parts[count++] = source->source_parameter;
    parts[count++] = ", ";
    parts[count++] = source->options;
    parts[count++] = ") AS TABLE (";
    parts[count++] = source->sites;
    parts[count++] = site_lists;
    parts[count++] = select_head;
    parts[count++] = genotype_model;
    parts[count++] = no_map_arguments;
    parts[count++] = decode_end;
    parts[count++] = from_lists;
    parts[count++] = genotype_guard;
    parts[count++] = source->source_guard;
    parts[count++] = "), (path, ";
    parts[count++] = source->source_parameter;
    parts[count++] = ", genetic_map VARCHAR, ";
    parts[count++] = source->options;
    parts[count++] = ") AS TABLE (";
    parts[count++] = source->sites;
    parts[count++] = site_lists;
    parts[count++] = map_lists;
    parts[count++] = select_head;
    parts[count++] = genotype_model;
    parts[count++] = map_arguments;
    parts[count++] = decode_end;
    parts[count++] = from_lists_with_map;
    parts[count++] = genotype_guard;
    parts[count++] = source->source_guard;
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
    parts[count++] = ancestry_frequencies;
    parts[count++] = site_lists;
    parts[count++] = select_head;
    parts[count++] = genotype_model;
    parts[count++] = no_map_arguments;
    parts[count++] = decode_end;
    parts[count++] = from_lists;
    parts[count++] = ancestry_guard;
    parts[count++] = "), (path, reference_table, proportions_table, genetic_map VARCHAR, "
                     "af_clamp := 1e-3, " ROH_OPTIONS ") AS TABLE (";
    parts[count++] = ancestry_sites;
    parts[count++] = ancestry_validations;
    parts[count++] = ancestry_frequencies;
    parts[count++] = site_lists;
    parts[count++] = map_lists;
    parts[count++] = select_head;
    parts[count++] = genotype_model;
    parts[count++] = map_arguments;
    parts[count++] = decode_end;
    parts[count++] = from_lists_with_map;
    parts[count++] = ancestry_guard;
    parts[count++] = ")";
    return duckhts_register_sql_parts(registration, parts, count);
}

static bool register_roh_counts_macro(duckhts_registration_t *registration, bool fitted) {
    const char *parts[ROH_MAX_PARTS];
    size_t count = 0;
    parts[count++] = fitted
        ? "CREATE OR REPLACE MACRO duckhts_roh_counts_fitted(counts_table, fit_table, "
          ROH_COUNTS_FITTED_OPTIONS ") AS TABLE ("
        : "CREATE OR REPLACE MACRO duckhts_roh_counts(counts_table, "
          ROH_COUNTS_OPTIONS ") AS TABLE (";
    parts[count++] = fitted ? fitted_counts_lists : counts_lists;
    parts[count++] = select_head;
    parts[count++] = fitted ? fitted_counts_model : counts_model;
    parts[count++] = no_map_arguments;
    parts[count++] = decode_end;
    parts[count++] = from_lists;
    parts[count++] = fitted ? fitted_counts_guard : counts_guard;
    parts[count++] = fitted
        ? "), (counts_table, fit_table, genetic_map VARCHAR, " ROH_COUNTS_FITTED_OPTIONS ") AS TABLE ("
        : "), (counts_table, genetic_map VARCHAR, " ROH_COUNTS_OPTIONS ") AS TABLE (";
    parts[count++] = fitted ? fitted_counts_lists : counts_lists;
    parts[count++] = map_lists;
    parts[count++] = select_head;
    parts[count++] = fitted ? fitted_counts_model : counts_model;
    parts[count++] = map_arguments;
    parts[count++] = decode_end;
    parts[count++] = from_lists_with_map;
    parts[count++] = fitted ? fitted_counts_guard : counts_guard;
    parts[count++] = ")";
    return duckhts_register_sql_parts(registration, parts, count);
}

bool register_duckhts_roh_sql(duckhts_registration_t *registration) {
    static const roh_source_t from_info_tag = {"duckhts_roh", "af_tag", ROH_OPTIONS,
                                               sites_from_info_tag, ""};
    static const roh_source_t from_relation = {"duckhts_roh_af_table", "af_table",
                                               ROH_OPTIONS ", site_rule := 'matching'",
                                               sites_from_relation, af_table_guard};
    const char *classifier[] = {htslib_is_indel_macro};
    return duckhts_register_sql_parts(registration, classifier, 1) &&
           register_roh_macro(registration, &from_info_tag) &&
           register_roh_macro(registration, &from_relation) &&
           register_roh_ancestry_macro(registration) &&
           register_roh_counts_macro(registration, false) &&
           register_roh_counts_macro(registration, true);
}
