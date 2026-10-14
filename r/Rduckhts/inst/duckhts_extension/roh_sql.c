/* duckhts_roh, duckhts_roh_af_table and duckhts_roh_ancestry: runs of homozygosity
 * from a VCF/BCF; duckhts_roh_counts: from read counts in a caller relation.
 *
 * The macros read records with read_bcf(tidy_format := true), collect the
 * sites of each sample and chromosome with the native __duckhts_roh_sites
 * aggregate, and decode each packed site list with __duckhts_roh_decode.
 * Relations (frequencies, genetic map) are passed by table name through
 * query_table, as in duckhts_cgranges_from_table. A macro cannot take an
 * absent relation, so the genetic map uses a separate _map macro and the frequency
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

/* duckhts__htslib_is_indel(ref, alts): whether HTSlib's bcf_get_variant_types()
 * sets VCF_INDEL for a record with these alleles, a port of
 * bcf_set_variant_type() (third_party/htslib vcf.c:5415, HTSlib 1.24).
 * An allele is an indel when, after the common prefix (compared without case),
 * one allele is used up and the other is not and does not end in a breakend
 * bracket; or when both remain and, after the common suffix, one shrinks to
 * its first base and the base where the other stops equals it. '*', symbolic
 * and breakend alleles and a single-base pair are never indels, so an MNP such
 * as AC>GT and a complex substitution such as ACG>AT are not. The record is an
 * indel when any ALT is. ind_p is the common prefix length; ind_s the common
 * suffix length of the rest, leaving at least one base of each. */
static const char htslib_is_indel_macro[] =
    "CREATE OR REPLACE MACRO duckhts__htslib_is_indel(ref, alts) AS "
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

/* duckhts_roh_af_table: a record without a matching frequency row is dropped
 * here, before the sites are collected, so it never claims a repeated position.
 * With a left join it would arrive before or after the matching record of its
 * position as the hash join happens to emit them, and the record kept would
 * depend on the plan. (bcftools roh --AF-file chooses differently at a repeated
 * position; duckhts_roh_af_table_bcftools below follows it.) */
static const char sites_from_relation[] =
    "WITH __dht_freq AS (SELECT CAST(chrom AS VARCHAR) AS f_chrom, CAST(pos AS BIGINT) AS f_pos, "
    "CAST(ref AS VARCHAR) AS f_ref, CAST(alt AS VARCHAR) AS f_alt, CAST(af AS DOUBLE) AS f_af "
    "FROM query_table(af_table)), "
    "__dht_calls AS NOT MATERIALIZED (SELECT CHROM AS chrom, POS AS pos, REF AS ref, "
    "array_to_string(ALT, ',') AS alt, SAMPLE_ID AS smp, " ROH_EVIDENCE " AS evidence "
    ROH_RECORDS "), "
    "__dht_sites AS NOT MATERIALIZED (SELECT c.chrom, c.pos, c.smp, f.f_af AS af, c.evidence "
    "FROM __dht_calls AS c JOIN __dht_freq AS f ON f.f_chrom = c.chrom "
    "AND f.f_pos = c.pos AND f.f_ref = c.ref AND f.f_alt = c.alt), ";

/* af_order names a column of af_table that gives the line order of a position.
 * A macro can select a column only by a pattern, and a pattern that matches no
 * column is a bind error, so the column is matched together with a NULL
 * sentinel column added last: the list holds [order value, NULL] when af_order
 * names a column (in any case, as SQL identifiers), and [NULL] when it does not
 * or is NULL. The list length tells the two apart, also on an empty relation. */
#define ROH_AF_ORDER_SENTINEL "__dht_no_af_order"
#define ROH_AF_ORDER_VALUES \
    "list_value(*COLUMNS(lambda c: lower(c) = lower(af_order) OR c = '" ROH_AF_ORDER_SENTINEL "'))"
#define ROH_AF_ORDER_UNKNOWN \
    "error('duckhts_roh_af_table_bcftools: af_order ''' || af_order || ''' is not a column of af_table')"

/* The caller's frequency relation, typed, with f_order the af_order value of
 * each row (NULL when af_order is NULL). Each use scans it again, rather than
 * holding it for several readers (a CTE read twice is materialized). When
 * af_order is NULL the order columns fold to constants and are pruned, so the
 * plan reads only chrom, pos, ref, alt and af. */
#define ROH_FREQ_ROWS \
    "(SELECT CAST(chrom AS VARCHAR) AS f_chrom, CAST(pos AS BIGINT) AS f_pos, " \
    "CAST(ref AS VARCHAR) AS f_ref, CAST(alt AS VARCHAR) AS f_alt, CAST(af AS DOUBLE) AS f_af, " \
    "CASE WHEN af_order IS NOT NULL AND len(__dht_af_order) = 1 THEN " ROH_AF_ORDER_UNKNOWN " " \
    "ELSE __dht_af_order[1] END AS f_order " \
    "FROM (SELECT *, " ROH_AF_ORDER_VALUES " AS __dht_af_order " \
    "FROM (SELECT *, NULL AS " ROH_AF_ORDER_SENTINEL " FROM query_table(af_table))))"

/* The value of `column` in the row that bcftools reads at a position: any row
 * when af_order is NULL (the guard allows one), else the row with the least
 * f_order, as bcftools keeps the first line. arg_min_null, unlike arg_min,
 * returns that row's value even when it is NULL. */
#define ROH_FREQ_FIRST(column) \
    "CASE WHEN af_order IS NULL THEN any_value(" column ") ELSE arg_min_null(" column ", f_order) END"

/* duckhts_roh_af_table_bcftools: the record of a listed position as bcftools
 * roh --AF-file chooses it. HTSlib functions are cited with lines of the
 * vendored third_party/htslib/synced_bcf_reader.c (HTSlib 1.24); bcftools is
 * not vendored, so vcfroh.c lines are those of the bcftools 1.23 mirror under
 * .sync/RBCFTools, and the tests compare with bcftools 1.24. Checked against
 * htslib develop 613169c1 (2026-10-08), where these targets paths are
 * unchanged; its new region merging (bcf_sr_regions_next_ with aux,
 * multi-region iterators) applies to -r/-R regions only.
 *
 * vcfroh's main() makes the frequency file a target list with alleles in its
 * third column (bcf_sr_set_targets(..., 1, 3), vcfroh.c:1272-1276). The
 * first target line at a position stays current, since
 * _bcf_sr_regions_overlap() (1793) advances only past lines that end before
 * the record (`while ( iseq==reg->iseq && reg->end < start )`, 1815), and
 * _regions_match_alleles() (1725) parses the alleles of that current line only
 * (reg->line.s, 1738). So bcftools reads the first line of a position. A
 * relation has no line order: without af_order it must give one row per
 * position, and with af_order the row with the least af_order value stands for
 * the first line (ROH_FREQ_FIRST; the guard below checks both). The row is an
 * indel when any of its alleles is longer than one character. At a listed
 * position, in file order (FILE_OFFSET), bcf_sr_next_line() (1001) passes over
 * a record whose indel class differs from the row's (the allele match, 1014)
 * when the next record of the file has the same POS
 * (`buffer[1]->pos!=buffer[0]->pos`, 1022); otherwise the record is seen. In
 * vcfroh(), a record with no non-symbolic ALT, or several, claims nothing
 * (vcfroh.c:1010-1027); the first remaining record claims the position
 * (prev_pos, 1036-1065), and process_line() skips the position unless that
 * record's alleles equal the row's (read_AF(), 661-686). A claiming record
 * whose first ALT is symbolic is skipped as well: bcftools would decode its
 * second allele, and the evidence here is that of the first ALT.
 *
 * The next record of the file. In a sorted file the records of one position
 * are contiguous, so the record after a record at (chrom, pos) has the same POS
 * when another record of (chrom, pos) follows it, or, since the test compares
 * POS only, when it is the last record of its contig and the next contig's
 * first record has that POS. One scan without samples, projecting CHROM, POS,
 * REF, ALT and FILE_OFFSET, feeds one aggregate with two grouping sets: one
 * row per contig (its first and last record) and one row per record whose POS
 * the frequency relation lists (all other records of a contig fall into one
 * group). Windows then run over that small result, never over the file.
 * File order needs BGZF input. */
static const char sites_from_relation_bcftools[] =
    "WITH __dht_calls AS NOT MATERIALIZED (SELECT CHROM AS chrom, POS AS pos, "
    "FILE_OFFSET AS file_offset, SAMPLE_ID AS smp, " ROH_EVIDENCE " AS evidence "
    "FROM read_bcf(path, tidy_format := true, samples := samples, "
    "decode_error_policy := 'error', include_file_offset := true) " ROH_ALT_FILTER "), "
    "__dht_listed_pos AS NOT MATERIALIZED (SELECT DISTINCT f_pos AS pos FROM " ROH_FREQ_ROWS " AS f "
    "WHERE f_pos IS NOT NULL), "
    "__dht_scan AS NOT MATERIALIZED (SELECT r.CHROM AS chrom, r.POS AS pos, r.REF AS ref, "
    "r.ALT AS alts, CASE WHEN r.FILE_OFFSET IS NULL THEN "
    "error('duckhts_roh_af_table_bcftools needs file order: use a BGZF-compressed VCF or BCF') "
    "ELSE r.FILE_OFFSET END AS file_offset, "
    "CASE WHEN p.pos IS NOT NULL THEN r.FILE_OFFSET END AS listed_offset "
    "FROM read_bcf(path, samples := '', decode_error_policy := 'error', include_file_offset := true) AS r "
    "LEFT JOIN __dht_listed_pos AS p ON p.pos = r.POS), "
    "__dht_groups AS NOT MATERIALIZED (SELECT chrom, listed_offset, "
    "GROUPING(listed_offset) AS per_contig, min(file_offset) AS first_offset, "
    "arg_min(pos, file_offset) AS first_pos, max(file_offset) AS last_offset, "
    "any_value(pos) AS pos, any_value(ref) AS ref, any_value(alts) AS alts "
    "FROM __dht_scan GROUP BY GROUPING SETS ((chrom), (chrom, listed_offset))), "
    /* Per contig, the first POS of the next contig in file order (lead over the
     * contig rows), then spread with the contig's last offset to its records. */
    "__dht_contig_next AS NOT MATERIALIZED (SELECT *, CASE WHEN per_contig = 1 THEN "
    "lead(first_pos) OVER (PARTITION BY per_contig ORDER BY first_offset) END AS next_first_pos "
    "FROM __dht_groups WHERE per_contig = 1 OR listed_offset IS NOT NULL), "
    "__dht_marked AS NOT MATERIALIZED (SELECT chrom, listed_offset, per_contig, pos, ref, alts, "
    "max(next_first_pos) OVER (PARTITION BY chrom) AS next_first_pos, "
    "max(CASE WHEN per_contig = 1 THEN last_offset END) OVER (PARTITION BY chrom) AS contig_last_offset, "
    "max(listed_offset) OVER (PARTITION BY chrom, pos) AS position_last_offset "
    "FROM __dht_contig_next), ";
/* The second part of sites_from_relation_bcftools (C99 bounds one string literal). */
static const char sites_from_relation_bcftools_claims[] =
    "__dht_bcftools_records AS NOT MATERIALIZED (SELECT m.chrom, m.pos, m.ref, m.alts, "
    "m.listed_offset AS file_offset, t.t_ref, t.t_alt, t.t_af, t.t_indel, "
    "m.listed_offset < m.position_last_offset OR (m.listed_offset = m.contig_last_offset "
    "AND m.next_first_pos = m.pos) AS followed_at_pos "
    "FROM __dht_marked AS m "
    "JOIN (SELECT f_chrom AS t_chrom, f_pos AS t_pos, "
    ROH_FREQ_FIRST("f_ref") " AS t_ref, " ROH_FREQ_FIRST("f_alt") " AS t_alt, "
    ROH_FREQ_FIRST("f_af") " AS t_af, "
    "coalesce(length(" ROH_FREQ_FIRST("f_ref") ") > 1 OR "
    "list_bool_or(list_transform(string_split(" ROH_FREQ_FIRST("f_alt") ", ','), lambda a: length(a) > 1)), "
    "false) AS t_indel FROM " ROH_FREQ_ROWS " AS f WHERE f_chrom IS NOT NULL AND f_pos IS NOT NULL "
    "GROUP BY f_chrom, f_pos) AS t ON t.t_chrom = m.chrom AND t.t_pos = m.pos "
    "WHERE m.per_contig = 0), "
    "__dht_bcftools_claims AS NOT MATERIALIZED (SELECT chrom, pos, file_offset, t_af, "
    "coalesce(ref = t_ref AND array_to_string(alts, ',') = t_alt AND alts[1] NOT IN " ROH_SYMBOLIC ", "
    "false) AS usable FROM __dht_bcftools_records "
    "WHERE CASE WHEN followed_at_pos THEN duckhts__htslib_is_indel(ref, alts) = t_indel ELSE true END "
    "AND len(list_filter(alts, lambda a: a NOT IN " ROH_SYMBOLIC ")) = 1 "
    "QUALIFY row_number() OVER (PARTITION BY chrom, pos ORDER BY file_offset) = 1), "
    "__dht_sites AS NOT MATERIALIZED (SELECT c.chrom, c.pos, c.smp, k.t_af AS af, c.evidence "
    "FROM __dht_calls AS c JOIN __dht_bcftools_claims AS k ON k.chrom = c.chrom "
    "AND k.pos = c.pos AND k.file_offset = c.file_offset WHERE k.usable), ";

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

/* __duckhts_roh_decode(sites, gt_error, seq_error, contamination,
 * allele_balance, artefact_weight, map_pos, map_cm, rec_rate, hw_to_az,
 * az_to_hw): the model parameters of the evidence kind, then the map
 * arguments, then the transitions. duckhts_roh_counts_fitted calls
 * __duckhts_roh_decode_stable, which takes artefact_weight_sd after
 * artefact_weight and gives each run a stable_fraction and gained_bases. */
#define ROH_RUN_COLUMNS \
    " SELECT seg.smp AS \"sample\", seg.chrom AS chrom, seg.r.start AS start, " \
    "seg.r.\"end\" AS \"end\", seg.r.\"end\" - seg.r.start + 1 AS length, " \
    "seg.r.n_markers AS n_markers, seg.r.quality AS quality"
static const char select_head[] =
    ROH_RUN_COLUMNS " FROM (SELECT l.smp, l.chrom, unnest(__duckhts_roh_decode(l.s, ";
#define ROH_NO_COUNT_MODEL "NULL::DOUBLE, NULL::DOUBLE, NULL::DOUBLE, NULL::DOUBLE"
#define ROH_COUNT_MODEL \
    "CAST(seq_error AS DOUBLE), CAST(contamination AS DOUBLE), " \
    "CAST(allele_balance AS DOUBLE), CAST(artefact_weight AS DOUBLE)"
/* The artefact_weight_sd of __duckhts_roh_valid_args() for a macro without
 * an interval. */
#define ROH_NO_INTERVAL "NULL::DOUBLE"
static const char genotype_model[] = "CAST(gt_error AS DOUBLE), " ROH_NO_COUNT_MODEL ", ";
static const char counts_model[] = "NULL::DOUBLE, " ROH_COUNT_MODEL ", ";
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
#define ROH_FITTED_GUARD_COLUMNS ROH_GUARD_COLUMNS ", NULL::DOUBLE, NULL::BIGINT "
#define ROH_GENOTYPE_ARGUMENTS_VALID \
    "__duckhts_roh_valid_args(" ROH_KIND ", CAST(gt_error AS DOUBLE), " ROH_NO_COUNT_MODEL ", " \
    ROH_NO_INTERVAL ", CAST(rec_rate AS DOUBLE), CAST(hw_to_az AS DOUBLE), CAST(az_to_hw AS DOUBLE), " \
    ROH_LIMIT_ARGUMENTS ")"
static const char genotype_guard[] = ROH_GUARD_COLUMNS "WHERE NOT " ROH_GENOTYPE_ARGUMENTS_VALID;
/* duckhts_roh_af_table_bcftools also checks the rows of each position, naming
 * the first failing position, whether or not a record is listed: without
 * af_order, at most one row; with it, a non-NULL af_order value in every row and
 * no value twice within a position. When af_order is NULL the CASE expressions
 * fold to the one-row check. A second arm checks that af_order names a column
 * even when af_table is empty: it joins a one-row sentinel relation to no row
 * of af_table, so the list of matched columns has its length without data. */
static const char af_table_bcftools_guard[] =
    ROH_GUARD_COLUMNS "FROM (SELECT f_chrom, f_pos, count(*) AS n, count(f_order) AS n_ordered "
    "FROM " ROH_FREQ_ROWS " AS f WHERE f_chrom IS NOT NULL AND f_pos IS NOT NULL "
    "GROUP BY f_chrom, f_pos HAVING CASE WHEN af_order IS NULL THEN count(*) > 1 "
    "ELSE count(f_order) < count(*) OR count(DISTINCT f_order) < count(*) END "
    "ORDER BY f_chrom, f_pos LIMIT 1) AS d "
    "WHERE error(CASE WHEN af_order IS NULL THEN "
    "'duckhts_roh_af_table_bcftools needs one af_table row per position, but ' || "
    "d.f_chrom || ':' || d.f_pos || ' has ' || d.n || ' rows; bcftools reads the first line "
    "of a position, and a relation has no line order; give af_order to choose it' "
    "WHEN d.n_ordered < d.n THEN 'duckhts_roh_af_table_bcftools needs an af_order value in every "
    "af_table row, but ' || d.f_chrom || ':' || d.f_pos || ' has a row with NULL ' || af_order "
    "ELSE 'duckhts_roh_af_table_bcftools needs distinct af_order values within a position, but ' || "
    "d.f_chrom || ':' || d.f_pos || ' repeats a value of ' || af_order END) IS NULL"
    ROH_GUARD_COLUMNS "FROM (SELECT len(" ROH_AF_ORDER_VALUES ") AS n "
    "FROM (SELECT NULL AS " ROH_AF_ORDER_SENTINEL ") AS s "
    "LEFT JOIN (SELECT * FROM query_table(af_table) LIMIT 0) AS f ON true) AS k "
    "WHERE CASE WHEN af_order IS NOT NULL AND k.n = 1 THEN " ROH_AF_ORDER_UNKNOWN " END IS NOT NULL";
/* The ancestry macro also runs its relation checks (__dht_group_guard) here,
 * so they do not depend on a call matching a reference site. */
static const char ancestry_guard[] =
    ROH_GUARD_COLUMNS "FROM __dht_group_guard AS g WHERE NOT (g.valid AND "
    ROH_GENOTYPE_ARGUMENTS_VALID ")";
#define ROH_COUNTS_ARGUMENTS_VALID \
    "__duckhts_roh_valid_args(3, NULL::DOUBLE, " ROH_COUNT_MODEL ", " ROH_NO_INTERVAL ", " \
    "CAST(rec_rate AS DOUBLE), CAST(hw_to_az AS DOUBLE), " \
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
/* The four fitted numbers of the sample `sample` and the artefact weight's
 * interval, in the order of __duckhts_roh_decode_stable() and
 * __duckhts_roh_valid_args(). A NULL artefact_weight_sd (fewer than three
 * blocks in the fit) counts as 0: the perturbed decodes equal the central
 * one, so every stable_fraction is 1 and every gained_bases 0. */
#define ROH_FIT_MODEL(from, sample) \
    ROH_FIT_VALUE(from, "seq_error", sample) ", " ROH_FIT_VALUE(from, "contamination", sample) ", " \
    ROH_FIT_VALUE(from, "allele_balance", sample) ", " ROH_FIT_VALUE(from, "artefact_weight", sample) ", " \
    "coalesce(" from ".artefact_weight_sd, 0.0)"
static const char fitted_counts_model[] = "NULL::DOUBLE, " ROH_FIT_MODEL("l", "l.smp") ", ";
/* Every sample of the counts relation is checked against its fit, also one
 * whose sites are all skipped or whose chromosomes the genetic map drops from
 * the decode. The check runs on every site list, before the decode, so the
 * counts relation is read once: read a second time for the check, a relation
 * that holds a join (a view filtered by an IN list, or counts joined to
 * frequencies) is a common subplan that DuckDB materializes with every row
 * and column. __duckhts_roh_valid_args() is TRUE or raises the error, and a
 * NULL fitted number raises in ROH_FIT_VALUE. The _map macro joins the map
 * with a LEFT JOIN and decodes only the lists with map nodes, so a list whose
 * chromosome has none is checked and gives no run. */
#define ROH_FITTED_VALID \
    "__duckhts_roh_valid_args(3, NULL::DOUBLE, " ROH_FIT_MODEL("l", "l.smp") ", " \
    "CAST(rec_rate AS DOUBLE), CAST(hw_to_az AS DOUBLE), CAST(az_to_hw AS DOUBLE), " \
    ROH_LIMIT_ARGUMENTS ")"
static const char fitted_select_head[] =
    ROH_RUN_COLUMNS ", seg.r.stable_fraction AS stable_fraction, seg.r.gained_bases AS gained_bases "
    "FROM (SELECT l.smp, l.chrom, unnest(CASE WHEN NOT " ROH_FITTED_VALID " THEN NULL "
    "ELSE __duckhts_roh_decode_stable(l.s, ";
static const char fitted_select_head_with_map[] =
    ROH_RUN_COLUMNS ", seg.r.stable_fraction AS stable_fraction, seg.r.gained_bases AS gained_bases "
    "FROM (SELECT l.smp, l.chrom, unnest(CASE WHEN NOT " ROH_FITTED_VALID " THEN NULL "
    "WHEN m.m_chrom IS NULL THEN NULL ELSE __duckhts_roh_decode_stable(l.s, ";
static const char fitted_decode_end[] =
    ", CAST(rec_rate AS DOUBLE), CAST(hw_to_az AS DOUBLE), CAST(az_to_hw AS DOUBLE)) END";
static const char fitted_from_lists_with_map[] =
    ") AS r FROM __dht_lists AS l LEFT JOIN __dht_map AS m ON m.m_chrom = l.chrom) AS seg";
/* The other arguments, checked also when there is no sample. */
static const char fitted_counts_guard[] =
    ROH_FITTED_GUARD_COLUMNS "FROM (SELECT 0.001 AS seq_error, 0.0 AS contamination, 0.5 AS allele_balance, "
    "0.0 AS artefact_weight) AS defaults WHERE NOT " ROH_COUNTS_ARGUMENTS_VALID;

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
/* New options go after the existing ones: a macro's options can be given by
 * position, so inserting one would shift every later positional value. */
#define ROH_COUNTS_OPTIONS \
    "seq_error := 1e-3, contamination := 0.0, " ROH_COUNTS_FITTED_OPTIONS \
    ", allele_balance := 0.5, artefact_weight := 0.0"

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
    "CAST(contamination AS DOUBLE) AS contamination, CAST(allele_balance AS DOUBLE) AS allele_balance, "
    "CAST(artefact_weight AS DOUBLE) AS artefact_weight, "
    "CAST(artefact_weight_sd AS DOUBLE) AS artefact_weight_sd, "
    "count(*) OVER (PARTITION BY CAST(sample AS VARCHAR)) AS fit_rows FROM query_table(fit_table) "
    "QUALIFY row_number() OVER (PARTITION BY CAST(sample AS VARCHAR)) = 1), "
    "__dht_lists AS (SELECT c.smp, c.chrom, fit.fit_rows, fit.seq_error, fit.contamination, "
    "fit.allele_balance, fit.artefact_weight, fit.artefact_weight_sd, "
    "__duckhts_roh_sites(c.pos, c.af, c.ref_count, c.alt_count, NULL::INTEGER, 3, " ROH_LIMIT_ARGUMENTS ") AS s "
    "FROM __dht_counts AS c LEFT JOIN __dht_fit AS fit ON fit.smp IS NOT DISTINCT FROM c.smp "
    "GROUP BY c.smp, c.chrom, fit.fit_rows, fit.seq_error, fit.contamination, fit.allele_balance, "
    "fit.artefact_weight, fit.artefact_weight_sd)";

#define ROH_MAX_PARTS 40

typedef struct {
    const char *name;
    const char *source_parameter;
    const char *options;      /* the named parameters after the source, with their defaults */
    const char *sites;
    const char *sites_continued; /* the rest of sites when one literal is too long; may be empty */
    const char *source_guard; /* checks of the source's own options; may be empty */
} roh_source_t;

/* Register the base macro and its genetic-map form under separate names. */
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
    parts[count++] = source->sites_continued;
    parts[count++] = site_lists;
    parts[count++] = select_head;
    parts[count++] = genotype_model;
    parts[count++] = no_map_arguments;
    parts[count++] = decode_end;
    parts[count++] = from_lists;
    parts[count++] = genotype_guard;
    parts[count++] = source->source_guard;
    parts[count++] = ")";
    if (!duckhts_register_sql_parts(registration, parts, count)) {
        return false;
    }
    count = 0;
    parts[count++] = "CREATE OR REPLACE MACRO ";
    parts[count++] = source->name;
    parts[count++] = "_map(path, ";
    parts[count++] = source->source_parameter;
    parts[count++] = ", genetic_map VARCHAR, ";
    parts[count++] = source->options;
    parts[count++] = ") AS TABLE (";
    parts[count++] = source->sites;
    parts[count++] = source->sites_continued;
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
    parts[count++] = ")";
    if (!duckhts_register_sql_parts(registration, parts, count)) {
        return false;
    }
    count = 0;
    parts[count++] = "CREATE OR REPLACE MACRO duckhts_roh_ancestry_map("
                     "path, reference_table, proportions_table, genetic_map VARCHAR, "
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
    parts[count++] = fitted ? fitted_select_head : select_head;
    parts[count++] = fitted ? fitted_counts_model : counts_model;
    parts[count++] = no_map_arguments;
    parts[count++] = fitted ? fitted_decode_end : decode_end;
    parts[count++] = from_lists;
    parts[count++] = fitted ? fitted_counts_guard : counts_guard;
    parts[count++] = ")";
    if (!duckhts_register_sql_parts(registration, parts, count)) {
        return false;
    }
    count = 0;
    parts[count++] = fitted
        ? "CREATE OR REPLACE MACRO duckhts_roh_counts_fitted_map("
          "counts_table, fit_table, genetic_map VARCHAR, " ROH_COUNTS_FITTED_OPTIONS ") AS TABLE ("
        : "CREATE OR REPLACE MACRO duckhts_roh_counts_map("
          "counts_table, genetic_map VARCHAR, " ROH_COUNTS_OPTIONS ") AS TABLE (";
    parts[count++] = fitted ? fitted_counts_lists : counts_lists;
    parts[count++] = map_lists;
    parts[count++] = fitted ? fitted_select_head_with_map : select_head;
    parts[count++] = fitted ? fitted_counts_model : counts_model;
    parts[count++] = map_arguments;
    parts[count++] = fitted ? fitted_decode_end : decode_end;
    parts[count++] = fitted ? fitted_from_lists_with_map : from_lists_with_map;
    parts[count++] = fitted ? fitted_counts_guard : counts_guard;
    parts[count++] = ")";
    return duckhts_register_sql_parts(registration, parts, count);
}

bool register_duckhts_roh_sql(duckhts_registration_t *registration) {
    static const roh_source_t from_info_tag = {"duckhts_roh", "af_tag", ROH_OPTIONS,
                                               sites_from_info_tag, "", ""};
    static const roh_source_t from_relation = {"duckhts_roh_af_table", "af_table", ROH_OPTIONS,
                                               sites_from_relation, "", ""};
    static const roh_source_t from_relation_bcftools = {
        "duckhts_roh_af_table_bcftools", "af_table", ROH_OPTIONS ", af_order := NULL",
        sites_from_relation_bcftools, sites_from_relation_bcftools_claims,
        af_table_bcftools_guard};
    const char *classifier[] = {htslib_is_indel_macro};
    return duckhts_register_sql_parts(registration, classifier, 1) &&
           register_roh_macro(registration, &from_info_tag) &&
           register_roh_macro(registration, &from_relation) &&
           register_roh_macro(registration, &from_relation_bcftools) &&
           register_roh_ancestry_macro(registration) &&
           register_roh_counts_macro(registration, false) &&
           register_roh_counts_macro(registration, true);
}
