#include "count_sql.h"
#include "duckhts_registration.h"
#include <stdbool.h>
#include <stddef.h>

#define CAL "CAST(seq_error AS DOUBLE), CAST(artefact_weight AS DOUBLE), CAST(homozygosity_excess AS DOUBLE)"
#define LIMITS DUCKHTS_WHOLE_OPTION("block_bases") ", " DUCKHTS_WHOLE_OPTION("max_depth") ", " \
    DUCKHTS_WHOLE_OPTION("min_blocks") ", " DUCKHTS_WHOLE_OPTION("max_blocks") ", " \
    DUCKHTS_WHOLE_OPTION("max_candidates") ", " DUCKHTS_WHOLE_OPTION("max_state_bytes")
#define FIT_ARGS CAL ", " DUCKHTS_WHOLE_OPTION("min_blocks") ", model = 'ibd1', " DUCKHTS_WHOLE_OPTION("max_state_bytes")

/* Exact binary sums keep fixed-size state per moment and block. Native
 * fitting receives one row per block, independent of panel site count. */
#define SUM(expr) "CAST(coalesce(__duckhts_mixture_sum(" expr "), 0) AS DOUBLE)"
#define S(a,b) SUM("n*(" a ")*(" b ")")
#define T(a,b) SUM("p*(" a ")*(" b ")")
#define TROW(a) T(a,"1") "," T(a,"ff") "," T(a,"gg") "," T(a,"ff*ff") "," T(a,"ff*gg") "," T(a,"gg*gg")
static const char *const mixture_sql[] = {
    "CREATE OR REPLACE MACRO duckhts_count_mixture_lowpass(counts_table, seq_error, artefact_weight, "
    "homozygosity_excess, block_bases := 10000000, max_depth := 8, min_blocks := 30, "
    "max_blocks := 512, max_candidates := 16, max_state_bytes := 268435456) AS TABLE (",
    "WITH raw AS NOT MATERIALIZED (SELECT "
    "CASE WHEN sample_id IS NULL OR CAST(sample_id AS VARCHAR) = '' THEN error('sample_id must not be empty') "
    "ELSE CAST(sample_id AS VARCHAR) END AS smp, "
    "CASE WHEN candidate IS NULL OR CAST(candidate AS VARCHAR) = '' THEN error('candidate must not be empty') "
    "ELSE CAST(candidate AS VARCHAR) END AS cand, "
    "CASE WHEN chrom IS NULL OR CAST(chrom AS VARCHAR) = '' THEN error('chrom must not be empty') "
    "ELSE CAST(chrom AS VARCHAR) END AS chr, "
    "CASE WHEN pos IS NULL OR pos < 1 THEN error('pos must be one-based') ELSE "
    DUCKHTS_WHOLE_AS("pos", "BIGINT", "pos must be whole") " END AS bp, "
    DUCKHTS_WHOLE_COUNT_AS("ref_count", "BIGINT") " AS r, "
    DUCKHTS_WHOLE_COUNT_AS("alt_count", "BIGINT") " AS a, "
    DUCKHTS_WHOLE_COUNT_AS("other_count", "BIGINT") " AS o, "
    "CAST(sample_af AS DOUBLE) AS f, CAST(candidate_af AS DOUBLE) AS g FROM query_table(counts_table)), ",
    "candidates AS (SELECT DISTINCT smp, cand FROM raw), "
    "candidate_sizes AS (SELECT smp, CASE WHEN count(*) > " DUCKHTS_WHOLE_OPTION("max_candidates")
    " THEN error('max_candidates exceeded') ELSE count(*) END AS nc FROM candidates GROUP BY smp), "
    "site_data AS (SELECT smp, chr, bp, first(r) AS r, first(a) AS a, first(o) AS o, first(f) AS f, "
    "CASE WHEN count(*) != count(DISTINCT cand) THEN error('duplicate site keys') "
    "WHEN count(DISTINCT struct_pack(r := r, a := a, o := o, f := f)) != 1 "
    "THEN error('counts and sample_af must agree across candidates') "
    "WHEN min(r) < 0 OR min(a) < 0 OR min(o) < 0 THEN error('read counts must not be negative') "
    "WHEN bool_or((f IS NOT NULL AND (NOT isfinite(f) OR f < 0 OR f > 1)) OR "
    "(g IS NOT NULL AND (NOT isfinite(g) OR g < 0 OR g > 1))) "
    "THEN error('frequencies must be finite and in [0,1]') "
    "ELSE count(g) END AS nf FROM raw GROUP BY smp, chr, bp), ",
    "prepared AS (SELECT s.smp, c.cand, s.chr, s.bp, s.r, s.a, s.o, s.f, r.g, "
    "s.nf = n.nc AND s.f IS NOT NULL AND s.r IS NOT NULL AND s.a IS NOT NULL AND s.o IS NOT NULL AS usable, "
    "coalesce(s.r::HUGEINT+s.a::HUGEINT > " DUCKHTS_WHOLE_OPTION("max_depth") ", false) AS excluded, "
    "dense_rank() OVER (PARTITION BY s.smp ORDER BY s.chr, (s.bp-1)//"
    DUCKHTS_WHOLE_OPTION("block_bases") ")-1 AS block "
    "FROM site_data s JOIN candidates c ON s.smp = c.smp JOIN candidate_sizes n ON s.smp = n.smp "
    "LEFT JOIN raw r ON s.smp = r.smp AND c.cand = r.cand AND s.chr = r.chr AND s.bp = r.bp), ",
    "effective AS (SELECT *, CASE WHEN usable AND NOT excluded THEN r+a ELSE 0 END AS n, "
    "CASE WHEN usable AND NOT excluded THEN a ELSE 0 END AS alt, "
    "coalesce(f,0) AS ff, coalesce(g,0) AS gg FROM prepared), "
    "pair_data AS (SELECT *, n*(n-1)//2 AS p, alt*(n-alt) AS discordant FROM effective), "
    "block_rows AS (SELECT smp, cand, block, struct_pack(counts := CAST([count(*), "
    "count(*) FILTER (WHERE usable AND NOT excluded), count(*) FILTER (WHERE usable AND NOT excluded AND n=0), "
    "count(*) FILTER (WHERE n=1), sum(n), coalesce(sum(o),0), count(*) FILTER (WHERE n>1), "
    "count(*) FILTER (WHERE excluded), sum((n-alt)*(n-alt-1)//2), sum(discordant), sum(alt*(alt-1)//2), "
    "count(*) FILTER (WHERE NOT usable AND NOT excluded)] AS BIGINT[]), "
    "s := [" S("1","1") "," S("1","ff") "," S("1","gg") ","
    S("ff","1") "," S("ff","ff") "," S("ff","gg") ","
    S("gg","1") "," S("gg","ff") "," S("gg","gg") "], "
    "y := [" SUM("alt") "," SUM("alt*ff") "," SUM("alt*gg") "], "
    "t := [" TROW("1") "," TROW("ff") ",",
    TROW("gg") "," TROW("ff*ff") ","
    TROW("ff*gg") "," TROW("gg*gg") "], "
    "z := [" SUM("discordant") "," SUM("discordant*ff") "," SUM("discordant*gg") ","
    SUM("discordant*ff*ff") "," SUM("discordant*ff*gg") "," SUM("discordant*gg*gg") "], "
    "denominator := " SUM("p*ff*(1-ff)") ", "
    "matched := bool_and(coalesce(f = g, false))) AS moments FROM pair_data GROUP BY smp, cand, block), ",
    "packed AS (SELECT smp, cand, __duckhts_mixture_blocks(block, moments, "
    DUCKHTS_WHOLE_OPTION("max_blocks") ", " DUCKHTS_WHOLE_OPTION("max_state_bytes")
    ") AS moments FROM block_rows GROUP BY smp, cand), "
    "fits AS (SELECT smp, cand, model, __duckhts_mixture_fit(moments, " FIT_ARGS ") AS fit "
    "FROM packed CROSS JOIN (VALUES ('unrelated'), ('ibd1')) models(model)) "
    "SELECT smp AS sample, cand AS candidate, model AS relationship_model, fit.* FROM fits "
    "WHERE model = 'unrelated' OR fit.matched_frequencies "
    "UNION ALL BY NAME SELECT NULL::VARCHAR AS sample, NULL::VARCHAR AS candidate, "
    "NULL::VARCHAR AS relationship_model, fit.* FROM (SELECT __duckhts_mixture_fit(NULL::BLOB, "
    CAL ", " DUCKHTS_WHOLE_OPTION("min_blocks") ", false, " DUCKHTS_WHOLE_OPTION("max_state_bytes")
    ") AS fit WHERE NOT __duckhts_mixture_valid_args(" CAL ", " LIMITS ")) guard",
    ")"
};
bool register_duckhts_count_mixture_sql(duckhts_registration_t *registration) {
    return duckhts_register_sql_parts(registration, mixture_sql, sizeof(mixture_sql)/sizeof(mixture_sql[0]));
}
