/* Copy-number fractions from read depth in fixed bins: design/depth_fraction_lane.md.
 *
 * duckhts_depth_bins(bins_table, ...) corrects bin counts for GC content with
 * a median read rate for each GC stratum. duckhts_depth_region_fraction and
 * duckhts_depth_sex_fraction turn corrected bins into the share of cells with
 * a one-copy change of a region, and the share of a second genome of the
 * other sex. Each estimator has an overload with a reference set.
 *
 * All three are SQL over caller-named relations. duckhts_depth_bins
 * materializes its input once and its holistic aggregates (median, quantile)
 * hold one value per autosomal bin and sample. The region estimator
 * materializes the narrow input columns once; the sex estimator streams its
 * input. Every other state is a few sums per sample, per sample and region,
 * or per GC stratum.
 */
#include "duckhts_registration.h"

#include <stdbool.h>
#include <stddef.h>

/* The 97.5% normal quantile, for intervals of +-1.96 standard errors. */
#define DEPTH_Z975 "1.959963985"

/* The 97.5% quantile of Student's t with `df` degrees of freedom: exact to six
 * decimals for 1 to 30 (R qt), and the Cornish-Fisher expansion of Abramowitz
 * and Stegun 26.7.5 above 30, which differs from qt by less than 3e-8 there. */
#define DEPTH_T975(df) \
    "CASE WHEN " df " < 1 THEN NULL WHEN " df " <= 30 THEN [12.706205, 4.302653, 3.182446, " \
    "2.776445, 2.570582, 2.446912, 2.364624, 2.306004, 2.262157, 2.228139, 2.200985, " \
    "2.178813, 2.160369, 2.144787, 2.131450, 2.119905, 2.109816, 2.100922, 2.093024, " \
    "2.085963, 2.079614, 2.073873, 2.068658, 2.063899, 2.059539, 2.055529, 2.051831, " \
    "2.048407, 2.045230, 2.042272][" df "] ELSE " DEPTH_Z975 \
    " + (" DEPTH_Z975 " ^ 3 + " DEPTH_Z975 ") / (4 * " df ")" \
    " + (5 * " DEPTH_Z975 " ^ 5 + 16 * " DEPTH_Z975 " ^ 3 + 3 * " DEPTH_Z975 ") / (96 * " df " ^ 2)" \
    " + (3 * " DEPTH_Z975 " ^ 7 + 19 * " DEPTH_Z975 " ^ 5 + 17 * " DEPTH_Z975 " ^ 3 - 15 * " \
    DEPTH_Z975 ") / (384 * " df " ^ 3)" \
    " + (79 * " DEPTH_Z975 " ^ 9 + 776 * " DEPTH_Z975 " ^ 7 + 1482 * " DEPTH_Z975 " ^ 5 - 1920 * " \
    DEPTH_Z975 " ^ 3 - 945 * " DEPTH_Z975 ") / (92160 * " df " ^ 4) END"

/* Autosomal means a contig whose join key is not X, Y or MT. */
#define DEPTH_AUTOSOMAL(chrom) "duckhts_contig_key(" chrom ") NOT IN ('X', 'Y', 'MT')"

/* floor(gc / gc_step) must fit the public INTEGER stratum column at gc = 1. */
#define DEPTH_GC_STEP_MIN "(1.0 / 2147483647)"

#define DEPTH_BINS_OPTIONS \
    "min_mappability := 0.9, gc_step := 0.01, min_stratum_bins := 10, " \
    "outlier_quantile := 0.01, max_ratio := 4.0"

/* ---- duckhts_depth_bins ------------------------------------------------- */

static const char *const depth_bins_sql[] = {
    "CREATE OR REPLACE MACRO duckhts_depth_bins(bins_table, " DEPTH_BINS_OPTIONS ") AS TABLE ",
    /* Options are checked once, before any row, and also for an empty input:
     * the guard arm at the end reads __dht_args. */
    "WITH __dht_args AS MATERIALIZED (SELECT CASE "
    "WHEN min_mappability IS NULL OR NOT (CAST(min_mappability AS DOUBLE) BETWEEN 0 AND 1) THEN "
    "error('duckhts_depth_bins: min_mappability must be in [0, 1]') "
    "WHEN gc_step IS NULL OR NOT (CAST(gc_step AS DOUBLE) BETWEEN " DEPTH_GC_STEP_MIN " AND 1) THEN "
    "error('duckhts_depth_bins: gc_step must be in [1/2147483647, 1]') "
    "WHEN min_stratum_bins IS NULL OR min_stratum_bins != trunc(min_stratum_bins) "
    "OR min_stratum_bins < 1 OR min_stratum_bins > 1000000000 THEN "
    "error('duckhts_depth_bins: min_stratum_bins must be a whole number in [1, 1000000000]') "
    "WHEN outlier_quantile IS NULL OR NOT (CAST(outlier_quantile AS DOUBLE) >= 0 "
    "AND CAST(outlier_quantile AS DOUBLE) < 0.5) THEN "
    "error('duckhts_depth_bins: outlier_quantile must be in [0, 0.5)') "
    "WHEN max_ratio IS NULL OR NOT (CAST(max_ratio AS DOUBLE) > 1 "
    "AND isfinite(CAST(max_ratio AS DOUBLE))) THEN "
    "error('duckhts_depth_bins: max_ratio must be finite and above 1') "
    "ELSE true END AS valid), ",
    /* mappability is optional: UNION ALL BY NAME adds it as NULL when the
     * relation lacks it. The other columns are required, so a missing one is
     * a binder error. */
    "__dht_in AS (SELECT CAST(sample_id AS VARCHAR) AS sample_id, CAST(chrom AS VARCHAR) AS chrom, "
    "CAST(start AS BIGINT) AS start, CAST(\"end\" AS BIGINT) AS \"end\", "
    "CASE WHEN reads != trunc(reads) THEN error('duckhts_depth_bins: reads must be whole numbers') "
    "ELSE CAST(reads AS BIGINT) END AS reads, "
    "CAST(gc AS DOUBLE) AS gc, CAST(mappability AS DOUBLE) AS mappability "
    "FROM (SELECT * FROM query_table(bins_table) "
    "UNION ALL BY NAME SELECT NULL::DOUBLE AS mappability WHERE false)), ",
    /* One materialized pass over the input: the limits, the strata and the
     * output all read it. reason holds the first rule a bin fails that does
     * not depend on its count. */
    "__dht_bins AS MATERIALIZED (SELECT "
    "CASE WHEN sample_id IS NULL THEN error('duckhts_depth_bins: sample_id must not be NULL') "
    "ELSE sample_id END AS sample_id, "
    "CASE WHEN chrom IS NULL THEN error('duckhts_depth_bins: chrom must not be NULL') "
    "ELSE chrom END AS chrom, "
    "CASE WHEN start IS NULL OR \"end\" IS NULL OR start < 0 OR \"end\" <= start THEN "
    "error('duckhts_depth_bins: bins need 0 <= start < end') ELSE start END AS start, \"end\", "
    "CASE WHEN reads < 0 THEN error('duckhts_depth_bins: reads must not be negative') "
    "ELSE reads END AS reads, "
    "CASE WHEN gc IS NOT NULL AND NOT (gc BETWEEN 0 AND 1) THEN "
    "error('duckhts_depth_bins: gc must be a share in [0, 1]') ELSE gc END AS gc, "
    "CASE WHEN mappability IS NOT NULL AND NOT (mappability BETWEEN 0 AND 1) THEN "
    "error('duckhts_depth_bins: mappability must be in [0, 1]') ELSE mappability END AS mappability, "
    DEPTH_AUTOSOMAL("chrom") " AS autosomal, "
    "CAST(floor(gc / CAST(gc_step AS DOUBLE)) AS INTEGER) AS gc_stratum, "
    "CAST(reads AS DOUBLE) / CAST(\"end\" - start AS DOUBLE) AS rate, "
    "CASE WHEN reads IS NULL THEN 'no_reads' WHEN gc IS NULL THEN 'no_gc' "
    "WHEN mappability < CAST(min_mappability AS DOUBLE) THEN 'low_mappability' END AS reason "
    "FROM __dht_in CROSS JOIN __dht_args WHERE __dht_args.valid), ",
    /* Training bins: autosomal, passing the bin rules, with at least one read
     * and a read rate inside the sample's outlier quantiles.
     * Continuous quantiles lie within the observed range. Floating-point
     * interpolation can exceed a constant rate by one ULP; clamping keeps
     * such a sample's positive bins eligible without a comparison tolerance. */
    "__dht_limits AS MATERIALIZED (SELECT sample_id, "
    "least(max(rate), greatest(min(rate), "
    "quantile_cont(rate, CAST(outlier_quantile AS DOUBLE)))) AS rate_low, "
    "least(max(rate), greatest(min(rate), "
    "quantile_cont(rate, 1 - CAST(outlier_quantile AS DOUBLE)))) AS rate_high "
    "FROM __dht_bins WHERE reason IS NULL AND autosomal AND reads > 0 GROUP BY sample_id), "
    "__dht_strata AS MATERIALIZED (SELECT b.sample_id, b.gc_stratum, median(b.rate) AS stratum_rate, "
    "count(*) AS stratum_bins FROM __dht_bins b JOIN __dht_limits l ON b.sample_id = l.sample_id "
    "WHERE b.reason IS NULL AND b.autosomal AND b.reads > 0 "
    "AND b.rate >= l.rate_low AND b.rate <= l.rate_high GROUP BY b.sample_id, b.gc_stratum), ",
    /* Read twice, by the centre and by the output: joining the small strata
     * table again costs less than a second copy of the bins. */
    "__dht_stratified AS NOT MATERIALIZED (SELECT b.*, s.stratum_rate, "
    "coalesce(s.stratum_bins, 0) AS stratum_bins, "
    "coalesce(b.reason, CASE WHEN coalesce(s.stratum_bins, 0) < min_stratum_bins "
    "THEN 'sparse_stratum' END) AS stratum_reason "
    "FROM __dht_bins b LEFT JOIN __dht_strata s "
    "ON b.sample_id = s.sample_id AND b.gc_stratum = s.gc_stratum), ",
    /* The ratio is scaled so that the median of the sample's positive-read
     * autosomal bins that pass the bin and stratum rules is 1. Empty bins
     * retain their output rows and zero ratios without setting the scale. */
    "__dht_centre AS MATERIALIZED (SELECT sample_id, median(rate / stratum_rate) AS centre "
    "FROM __dht_stratified WHERE stratum_reason IS NULL AND autosomal AND reads > 0 GROUP BY sample_id), "
    "__dht_expected AS (SELECT s.*, CASE WHEN s.stratum_reason IS NULL AND c.centre > 0 THEN "
    "s.stratum_rate * c.centre * CAST(s.\"end\" - s.start AS DOUBLE) END AS expected, "
    "coalesce(s.stratum_reason, CASE WHEN c.centre IS NULL OR c.centre <= 0 "
    "THEN 'no_autosomal_bins' END) AS centre_reason "
    "FROM __dht_stratified s LEFT JOIN __dht_centre c ON s.sample_id = c.sample_id), "
    "__dht_ratio AS (SELECT *, reads / expected AS ratio, coalesce(centre_reason, "
    "CASE WHEN reads / expected > CAST(max_ratio AS DOUBLE) THEN 'high_ratio' END) AS final_reason "
    "FROM __dht_expected) ",
    "SELECT sample_id, chrom, start, \"end\", reads, gc, mappability, gc_stratum, "
    "CAST(stratum_bins AS BIGINT) AS stratum_bins, expected, ratio, "
    "CASE WHEN reads > 0 THEN log2(ratio) END AS log2_ratio, final_reason AS reason, "
    "final_reason IS NULL AS kept, 'duckhts_depth_bins gc_stratum_median v1' AS method "
    "FROM __dht_ratio "
    "UNION ALL SELECT NULL::VARCHAR, NULL::VARCHAR, NULL::BIGINT, NULL::BIGINT, NULL::BIGINT, "
    "NULL::DOUBLE, NULL::DOUBLE, NULL::INTEGER, NULL::BIGINT, NULL::DOUBLE, NULL::DOUBLE, "
    "NULL::DOUBLE, NULL::VARCHAR, NULL::BOOLEAN, NULL::VARCHAR FROM __dht_args WHERE NOT valid",
};

/* ---- duckhts_depth_region_fraction --------------------------------------- */

/* Corrected bins of the tested samples (role 'sample') and of the reference
 * samples (role 'reference') go through one chain. */
#define REGION_INPUT(role, table) \
    "SELECT '" role "' AS role, CAST(sample_id AS VARCHAR) AS sample_id, " \
    "CAST(chrom AS VARCHAR) AS chrom, CAST(start AS BIGINT) AS start, " \
    "CAST(\"end\" AS BIGINT) AS \"end\", CAST(reads AS DOUBLE) AS reads, " \
    "CAST(expected AS DOUBLE) AS expected, CAST(kept AS BOOLEAN) AS kept FROM query_table(" table ")"

/* One materialized copy of the input's narrow columns: the per-sample sums
 * and the bin-region join both read it. A NULL kept counts as not kept. */
static const char region_rows[] =
    "__dht_rows AS MATERIALIZED (SELECT role, "
    "CASE WHEN sample_id IS NULL THEN error('duckhts_depth_region_fraction: sample_id must not be NULL') "
    "ELSE sample_id END AS sample_id, chrom, start, \"end\", coalesce(kept, false) AS kept, "
    "CASE WHEN kept AND NOT (reads >= 0 AND expected > 0 AND isfinite(expected)) THEN "
    "error('duckhts_depth_region_fraction: kept bins need reads >= 0 and a finite expected > 0') "
    "ELSE reads END AS reads, expected, " DEPTH_AUTOSOMAL("chrom") " AS autosomal, "
    /* Region endpoints are integers: floor(midpoint) gives the same half-open
     * membership without overflowing start + end or rounding large positions. */
    "start + (\"end\" - start) // 2 AS midpoint FROM __dht_input), "
    /* Sums over each sample's kept autosomal bins; a sample with no kept bin
     * keeps a row of zeros. */
    "__dht_samples AS (SELECT role, sample_id, count(*) FILTER (WHERE kept AND autosomal) AS n_a, "
    "coalesce(sum(reads) FILTER (WHERE kept AND autosomal), 0) AS c_a, "
    "coalesce(sum(expected) FILTER (WHERE kept AND autosomal), 0) AS e_a, "
    "coalesce(sum(reads * reads / expected) FILTER (WHERE kept AND autosomal), 0) AS q_a "
    "FROM __dht_rows GROUP BY role, sample_id), ";

/* The intervals of a region are merged first, so that a bin whose midpoint
 * lies in two overlapping intervals of one region joins it once. */
static const char region_members[] =
    "__dht_region_rows AS (SELECT "
    "CASE WHEN region IS NULL THEN error('duckhts_depth_region_fraction: region must not be NULL') "
    "ELSE CAST(region AS VARCHAR) END AS region, CAST(chrom AS VARCHAR) AS chrom, "
    "CASE WHEN start IS NULL OR \"end\" IS NULL OR start < 0 OR \"end\" <= start THEN "
    "error('duckhts_depth_region_fraction: regions need 0 <= start < end') "
    "ELSE CAST(start AS BIGINT) END AS start, CAST(\"end\" AS BIGINT) AS \"end\" "
    "FROM query_table(regions_table)), "
    "__dht_region_pieces AS (SELECT *, sum(CASE WHEN start <= reached THEN 0 ELSE 1 END) OVER ("
    "PARTITION BY region, chrom ORDER BY start, \"end\" ROWS UNBOUNDED PRECEDING) AS piece FROM ("
    "SELECT *, max(\"end\") OVER (PARTITION BY region, chrom ORDER BY start, \"end\" "
    "ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS reached FROM __dht_region_rows)), "
    "__dht_regions AS MATERIALIZED (SELECT region, chrom, min(start) AS start, max(\"end\") AS \"end\" "
    "FROM __dht_region_pieces GROUP BY region, chrom, piece), "
    /* A bin belongs to a region when its midpoint lies in one of the region's
     * intervals. */
    "__dht_region_sums AS (SELECT k.role, k.sample_id, r.region, count(*) AS n_r, "
    "sum(k.reads) AS c_r, sum(k.expected) AS e_r, count(*) FILTER (WHERE k.autosomal) AS n_ra, "
    "coalesce(sum(k.reads) FILTER (WHERE k.autosomal), 0) AS c_ra, "
    "coalesce(sum(k.expected) FILTER (WHERE k.autosomal), 0) AS e_ra, "
    "coalesce(sum(k.reads * k.reads / k.expected) FILTER (WHERE k.autosomal), 0) AS q_ra "
    "FROM __dht_rows k JOIN __dht_regions r "
    "ON k.chrom = r.chrom AND k.midpoint >= r.start AND k.midpoint < r.\"end\" "
    "WHERE k.kept GROUP BY k.role, k.sample_id, r.region), ";

static const char region_estimate[] =
    /* The baseline of a region is the sample's kept autosomal bins outside it. */
    "__dht_baseline AS (SELECT s.role, s.sample_id, n.region, coalesce(r.n_r, 0) AS bins, "
    "r.c_r, r.e_r, s.n_a - coalesce(r.n_ra, 0) AS baseline_bins, "
    "s.c_a - coalesce(r.c_ra, 0) AS c_b, s.e_a - coalesce(r.e_ra, 0) AS e_b, "
    "s.q_a - coalesce(r.q_ra, 0) AS q_b "
    "FROM __dht_samples s CROSS JOIN (SELECT DISTINCT region FROM __dht_regions) n "
    "LEFT JOIN __dht_region_sums r "
    "ON s.role = r.role AND s.sample_id = r.sample_id AND n.region = r.region), "
    /* Quasi-Poisson dispersion of the baseline bins: the sum of
     * (reads - R expected)^2 / (R expected) over n - 1, with R = c_b / e_b,
     * written with the sums; at least 1. */
    "__dht_dispersion AS (SELECT *, CASE WHEN baseline_bins >= 2 AND c_b > 0 THEN "
    "greatest(1.0, ((e_b / c_b) * q_b - c_b) / (baseline_bins - 1)) END AS dispersion, "
    "CASE WHEN bins > 0 AND baseline_bins >= 2 AND c_b > 0 THEN (c_r / e_r) / (c_b / e_b) END AS ratio "
    "FROM __dht_baseline), "
    /* The variance of a sum of counts is the dispersion times the sum, at
     * least one read, so that a region without reads keeps an interval. */
    "__dht_estimate AS (SELECT *, 2 * (ratio - 1) AS bin_shift, "
    "2 * sqrt(dispersion * greatest(c_r, 1) / (e_r * (c_b / e_b)) ^ 2 "
    "+ ratio ^ 2 * dispersion / greatest(c_b, 1)) AS bin_shift_se FROM __dht_dispersion) ";

static const char region_select_bins[] =
    "SELECT sample_id, region, CAST(bins AS BIGINT) AS bins, "
    "CAST(baseline_bins AS BIGINT) AS baseline_bins, ratio, bin_shift AS shift, "
    "bin_shift_se AS shift_se, bin_shift - " DEPTH_Z975 " * bin_shift_se AS shift_low, "
    "bin_shift + " DEPTH_Z975 " * bin_shift_se AS shift_high, dispersion, "
    "NULL::DOUBLE AS reference_shift, NULL::BIGINT AS reference_samples, "
    "CASE WHEN bins = 0 THEN 'no_bins' WHEN ratio IS NULL THEN 'no_baseline' ELSE 'ok' END AS status, "
    "'duckhts_depth_region_fraction bins v1' AS method FROM __dht_estimate";

/* The reference: each reference sample's shift for the region, leaving out
 * the tested sample when it is also a reference. The shift is the sample's
 * minus the reference mean. Its variance is the sample's own bin variance,
 * plus the reference variance beyond the reference samples' bin variance,
 * plus the variance of the reference mean. */
static const char region_select_reference[] =
    ", __dht_reference AS (SELECT region, count(bin_shift) AS n, sum(bin_shift) AS s, "
    "sum(bin_shift * bin_shift) AS ss, sum(bin_shift_se * bin_shift_se) AS vs "
    "FROM __dht_estimate WHERE role = 'reference' AND bin_shift IS NOT NULL GROUP BY region), "
    "__dht_own AS (SELECT sample_id, region, bin_shift, bin_shift_se FROM __dht_estimate "
    "WHERE role = 'reference' AND bin_shift IS NOT NULL), "
    "__dht_leave_one_out AS (SELECT t.*, "
    "coalesce(r.n, 0) - CAST(o.bin_shift IS NOT NULL AS BIGINT) AS ref_n, "
    "coalesce(r.s, 0) - coalesce(o.bin_shift, 0) AS ref_s, "
    "coalesce(r.ss, 0) - coalesce(o.bin_shift * o.bin_shift, 0) AS ref_ss, "
    "coalesce(r.vs, 0) - coalesce(o.bin_shift_se * o.bin_shift_se, 0) AS ref_vs "
    "FROM __dht_estimate t LEFT JOIN __dht_reference r ON t.region = r.region "
    "LEFT JOIN __dht_own o ON t.sample_id = o.sample_id AND t.region = o.region "
    "WHERE t.role = 'sample'), "
    "__dht_reference_moments AS (SELECT *, CASE WHEN ref_n >= 1 THEN ref_s / ref_n END AS ref_mean, "
    "CASE WHEN ref_n >= 2 THEN greatest(0.0, (ref_ss - ref_s * ref_s / ref_n) / (ref_n - 1)) END AS ref_var "
    "FROM __dht_leave_one_out), "
    "__dht_adjusted AS (SELECT *, bin_shift - ref_mean AS shift, sqrt(bin_shift_se ^ 2 "
    "+ greatest(0.0, ref_var - ref_vs / ref_n) + ref_var / ref_n) AS shift_se, "
    DEPTH_T975("CAST(ref_n - 1 AS BIGINT)") " AS t975 FROM __dht_reference_moments) "
    "SELECT sample_id, region, CAST(bins AS BIGINT) AS bins, "
    "CAST(baseline_bins AS BIGINT) AS baseline_bins, ratio, shift, shift_se, "
    "shift - t975 * shift_se AS shift_low, shift + t975 * shift_se AS shift_high, dispersion, "
    "ref_mean AS reference_shift, CAST(ref_n AS BIGINT) AS reference_samples, "
    "CASE WHEN bins = 0 THEN 'no_bins' WHEN ratio IS NULL THEN 'no_baseline' "
    "WHEN ref_n < 2 THEN 'few_reference_samples' ELSE 'ok' END AS status, "
    "'duckhts_depth_region_fraction reference v1' AS method FROM __dht_adjusted";

static const char *const depth_region_sql[] = {
    "CREATE OR REPLACE MACRO duckhts_depth_region_fraction(bins_table, regions_table) AS TABLE ("
    "WITH __dht_input AS (" REGION_INPUT("sample", "bins_table") "), ",
    region_rows,
    region_members,
    region_estimate,
    region_select_bins,
    ")",
};

static const char *const depth_region_ref_sql[] = {
    "CREATE OR REPLACE MACRO duckhts_depth_region_fraction_ref(bins_table, regions_table, reference_table) AS TABLE ("
    "WITH __dht_input AS (" REGION_INPUT("sample", "bins_table") " UNION ALL "
    REGION_INPUT("reference", "reference_table") "), ",
    region_rows,
    region_members,
    region_estimate,
    region_select_reference,
    ")",
};

/* ---- duckhts_depth_sex_fraction ------------------------------------------ */

#define SEX_INPUT(role, table, sex) \
    "SELECT '" role "' AS role, CAST(sample_id AS VARCHAR) AS sample_id, " sex " AS sex, " \
    "CAST(chrom AS VARCHAR) AS chrom, CAST(reads AS DOUBLE) AS reads, " \
    "CAST(expected AS DOUBLE) AS expected, CAST(kept AS BOOLEAN) AS kept FROM query_table(" table ")"

/* One pass over the input: each sample's sums of kept autosomal, X and Y
 * bins. A NULL kept counts as not kept. */
static const char sex_copies[] =
    "__dht_rows AS (SELECT role, CASE WHEN sample_id IS NULL THEN "
    "error('duckhts_depth_sex_fraction: sample_id must not be NULL') ELSE sample_id END AS sample_id, "
    "sex, coalesce(kept, false) AS kept, CASE WHEN kept AND NOT (reads >= 0 AND expected > 0 "
    "AND isfinite(expected)) THEN "
    "error('duckhts_depth_sex_fraction: kept bins need reads >= 0 and a finite expected > 0') "
    "ELSE reads END AS reads, expected, "
    "CASE WHEN chrom IS NULL THEN "
    "error('duckhts_depth_sex_fraction: chrom must not be NULL') "
    "ELSE CASE duckhts_contig_key(chrom) WHEN 'X' THEN 'X' WHEN 'Y' THEN 'Y' WHEN 'MT' THEN NULL "
    "ELSE 'A' END END AS kind FROM __dht_input), "
    "__dht_sums AS (SELECT role, sample_id, CASE WHEN min(sex) IS DISTINCT FROM max(sex) THEN "
    "error('duckhts_depth_sex_fraction: a reference sample must have one sex') ELSE min(sex) END AS sex, "
    "count(*) FILTER (WHERE kept AND kind = 'A') AS n_a, "
    "sum(reads) FILTER (WHERE kept AND kind = 'A') AS c_a, "
    "sum(expected) FILTER (WHERE kept AND kind = 'A') AS e_a, "
    "sum(reads * reads / expected) FILTER (WHERE kept AND kind = 'A') AS q_a, "
    "count(*) FILTER (WHERE kept AND kind = 'X') AS n_x, "
    "sum(reads) FILTER (WHERE kept AND kind = 'X') AS c_x, "
    "sum(expected) FILTER (WHERE kept AND kind = 'X') AS e_x, "
    "count(*) FILTER (WHERE kept AND kind = 'Y') AS n_y, "
    "sum(reads) FILTER (WHERE kept AND kind = 'Y') AS c_y, "
    "sum(expected) FILTER (WHERE kept AND kind = 'Y') AS e_y "
    "FROM __dht_rows GROUP BY role, sample_id), "
    /* Dispersion and depth of the autosomal bins, as in the region estimator. */
    "__dht_dispersion AS (SELECT *, CASE WHEN n_a >= 2 AND c_a > 0 THEN "
    "greatest(1.0, ((e_a / c_a) * q_a - c_a) / (n_a - 1)) END AS dispersion, "
    "CASE WHEN n_a >= 2 AND c_a > 0 THEN c_a / e_a END AS depth_a FROM __dht_sums), "
    /* Copies implied by depth: twice the depth ratio to the autosomes. */
    "__dht_copies AS (SELECT *, "
    "CASE WHEN n_x > 0 THEN 2 * (c_x / e_x) / depth_a END AS x_raw, "
    "CASE WHEN n_x > 0 THEN sqrt((2 / (e_x * depth_a)) ^ 2 * dispersion * greatest(c_x, 1) "
    "+ (2 * (c_x / e_x) / depth_a) ^ 2 * dispersion / greatest(c_a, 1)) END AS x_raw_se, "
    "CASE WHEN n_y > 0 THEN 2 * (c_y / e_y) / depth_a END AS y_raw, "
    "CASE WHEN n_y > 0 THEN sqrt((2 / (e_y * depth_a)) ^ 2 * dispersion * greatest(c_y, 1) "
    "+ (2 * (c_y / e_y) / depth_a) ^ 2 * dispersion / greatest(c_a, 1)) END AS y_raw_se "
    "FROM __dht_dispersion), ";

static const char sex_levels_ideal[] =
    "__dht_levels AS (SELECT *, 1.0 AS x_one, 2.0 AS x_two, 0.0 AS y_none, 1.0 AS y_one, "
    "NULL::BIGINT AS reference_samples, NULL::VARCHAR AS reference_status "
    "FROM __dht_copies WHERE role = 'sample'), ";

/* XY references set the depth of one X and one Y copy; XX references, when
 * there are any, set two X copies and the Y background. The tested sample is
 * left out of its own reference. */
static const char sex_levels_reference[] =
    "__dht_reference AS (SELECT sex, count(x_raw) AS nx, sum(x_raw) AS sx, count(y_raw) AS ny, "
    "sum(y_raw) AS sy, count(*) FILTER (WHERE x_raw IS NOT NULL OR y_raw IS NOT NULL) AS samples "
    "FROM __dht_copies WHERE role = 'reference' GROUP BY sex), "
    "__dht_reference_sums AS (SELECT t.sample_id, "
    "coalesce(xy.nx, 0) - CAST(o.sex = 'XY' AND o.x_raw IS NOT NULL AS BIGINT) AS xy_nx, "
    "coalesce(xy.sx, 0) - CASE WHEN o.sex = 'XY' THEN coalesce(o.x_raw, 0) ELSE 0 END AS xy_sx, "
    "coalesce(xy.ny, 0) - CAST(o.sex = 'XY' AND o.y_raw IS NOT NULL AS BIGINT) AS xy_ny, "
    "coalesce(xy.sy, 0) - CASE WHEN o.sex = 'XY' THEN coalesce(o.y_raw, 0) ELSE 0 END AS xy_sy, "
    "coalesce(xx.nx, 0) - CAST(o.sex = 'XX' AND o.x_raw IS NOT NULL AS BIGINT) AS xx_nx, "
    "coalesce(xx.sx, 0) - CASE WHEN o.sex = 'XX' THEN coalesce(o.x_raw, 0) ELSE 0 END AS xx_sx, "
    "coalesce(xx.ny, 0) - CAST(o.sex = 'XX' AND o.y_raw IS NOT NULL AS BIGINT) AS xx_ny, "
    "coalesce(xx.sy, 0) - CASE WHEN o.sex = 'XX' THEN coalesce(o.y_raw, 0) ELSE 0 END AS xx_sy, "
    "coalesce(xy.samples, 0) + coalesce(xx.samples, 0) "
    "- CAST(o.x_raw IS NOT NULL OR o.y_raw IS NOT NULL AS BIGINT) "
    "AS reference_samples "
    "FROM __dht_copies t LEFT JOIN (SELECT * FROM __dht_reference WHERE sex = 'XY') xy ON true "
    "LEFT JOIN (SELECT * FROM __dht_reference WHERE sex = 'XX') xx ON true "
    "LEFT JOIN (SELECT * FROM __dht_copies WHERE role = 'reference') o ON t.sample_id = o.sample_id "
    "WHERE t.role = 'sample'), ";

static const char sex_levels_reference_apply[] =
    "__dht_levels AS (SELECT c.*, CASE WHEN r.xy_nx > 0 THEN r.xy_sx / r.xy_nx END AS x_one, "
    "CASE WHEN r.xx_nx > 0 THEN r.xx_sx / r.xx_nx WHEN r.xy_nx > 0 THEN 2 * r.xy_sx / r.xy_nx END "
    "AS x_two, CASE WHEN r.xx_ny > 0 THEN r.xx_sy / r.xx_ny ELSE 0.0 END AS y_none, "
    "CASE WHEN r.xy_ny > 0 THEN r.xy_sy / r.xy_ny END AS y_one, r.reference_samples, "
    "CASE WHEN r.xy_nx = 0 AND r.xy_ny = 0 THEN 'no_xy_reference' END AS reference_status "
    "FROM __dht_copies c JOIN __dht_reference_sums r ON c.sample_id = r.sample_id "
    "WHERE c.role = 'sample'), ";

/* Calibrated copies, the major genome and the two estimates of f. A level
 * pair that does not separate (x_two <= x_one, y_one <= y_none) gives NULL. */
static const char sex_estimate[] =
    "__dht_calibrated AS (SELECT *, "
    "CASE WHEN x_two > x_one THEN 1 + (x_raw - x_one) / (x_two - x_one) END AS x_copies, "
    "CASE WHEN x_two > x_one THEN x_raw_se / (x_two - x_one) END AS x_copies_se, "
    "CASE WHEN y_one > y_none THEN (y_raw - y_none) / (y_one - y_none) END AS y_copies, "
    "CASE WHEN y_one > y_none THEN y_raw_se / (y_one - y_none) END AS y_copies_se "
    "FROM __dht_levels), "
    "__dht_major AS (SELECT *, CASE WHEN x_copies IS NOT NULL THEN "
    "CASE WHEN x_copies >= 1.5 THEN 'XX' ELSE 'XY' END "
    "WHEN y_copies IS NOT NULL THEN CASE WHEN y_copies >= 0.5 THEN 'XY' ELSE 'XX' END END AS major_sex "
    "FROM __dht_calibrated), "
    "__dht_parts AS (SELECT *, "
    "CASE major_sex WHEN 'XX' THEN 2 - x_copies WHEN 'XY' THEN x_copies - 1 END AS fraction_x, "
    "CASE major_sex WHEN 'XX' THEN y_copies WHEN 'XY' THEN 1 - y_copies END AS fraction_y, "
    "coalesce(1 / x_copies_se ^ 2, 0) + coalesce(1 / y_copies_se ^ 2, 0) AS weight "
    "FROM __dht_major), "
    "__dht_combined AS (SELECT *, CASE WHEN weight > 0 THEN (coalesce(fraction_x / x_copies_se ^ 2, 0) "
    "+ coalesce(fraction_y / y_copies_se ^ 2, 0)) / weight END AS fraction, "
    "CASE WHEN weight > 0 THEN 1 / sqrt(weight) END AS fraction_se, "
    "(fraction_x - fraction_y) / sqrt(x_copies_se ^ 2 + y_copies_se ^ 2) AS agreement_z "
    "FROM __dht_parts) ";

static const char sex_select[] =
    "SELECT sample_id, CAST(n_a AS BIGINT) AS autosomal_bins, CAST(n_x AS BIGINT) AS x_bins, "
    "CAST(n_y AS BIGINT) AS y_bins, dispersion, x_copies, x_copies_se, y_copies, y_copies_se, "
    "major_sex, fraction_x, x_copies_se AS fraction_x_se, fraction_y, y_copies_se AS fraction_y_se, "
    "fraction, fraction_se, fraction - " DEPTH_Z975 " * fraction_se AS fraction_low, "
    "fraction + " DEPTH_Z975 " * fraction_se AS fraction_high, agreement_z, reference_samples, "
    "CASE WHEN depth_a IS NULL THEN 'no_autosomal_bins' "
    "WHEN reference_status IS NOT NULL THEN reference_status "
    "WHEN x_copies IS NULL AND y_copies IS NULL AND (n_x > 0 OR n_y > 0) THEN 'bad_reference' "
    "WHEN n_x = 0 AND n_y = 0 THEN 'no_sex_bins' WHEN n_x = 0 THEN 'no_x_bins' "
    "WHEN n_y = 0 THEN 'no_y_bins' WHEN abs(agreement_z) > 3 THEN 'x_y_disagree' ELSE 'ok' END AS status, ";

static const char *const depth_sex_sql[] = {
    "CREATE OR REPLACE MACRO duckhts_depth_sex_fraction(bins_table) AS TABLE ("
    "WITH __dht_input AS (" SEX_INPUT("sample", "bins_table", "NULL::VARCHAR") "), ",
    sex_copies,
    sex_levels_ideal,
    sex_estimate,
    sex_select,
    "'duckhts_depth_sex_fraction bins v1' AS method FROM __dht_combined)",
};

static const char *const depth_sex_ref_sql[] = {
    "CREATE OR REPLACE MACRO duckhts_depth_sex_fraction_ref(bins_table, reference_table) AS TABLE ("
    "WITH __dht_input AS (" SEX_INPUT("sample", "bins_table", "NULL::VARCHAR") " UNION ALL "
    SEX_INPUT("reference", "reference_table",
              "CASE WHEN CAST(sex AS VARCHAR) IN ('XX', 'XY') THEN CAST(sex AS VARCHAR) "
              "ELSE error('duckhts_depth_sex_fraction: reference sex must be XX or XY') END") "), ",
    sex_copies,
    sex_levels_reference,
    sex_levels_reference_apply,
    sex_estimate,
    sex_select,
    "'duckhts_depth_sex_fraction reference v1' AS method FROM __dht_combined)",
};

/* ---- duckhts_depth_state_segments / duckhts_depth_state_solutions ------- */

/* The state model of depth_state_hmm.h (an ichorCNA v0.6.1 port) over kept
 * corrected bins. Bins are listed by chrom then start, as upstream orders
 * its chromosome factor; the model learns from autosomal bins and decodes
 * every bin. One decode runs per sample and macro call. */
#define DEPTH_STATE_OPTIONS \
    "normal_fractions := [0.5, 0.6, 0.7, 0.8, 0.9], ploidies := [2.0, 3.0], max_cn := 7, " \
    "txn_e := 0.9999999, txn_strength := 10000000.0, max_iter := 20, min_segment_bins := 50, " \
    "alt_fraction_threshold := 0.05, close_margin := 10.0, max_bins := 1000000"

#define DEPTH_STATE_ARGUMENTS \
    "normal_fractions, ploidies, max_cn, txn_e, txn_strength, max_iter, min_segment_bins, " \
    "alt_fraction_threshold, close_margin, max_bins"

static const char depth_state_lists[] =
    "WITH __dht_bins AS MATERIALIZED (SELECT "
    "CASE WHEN sample_id IS NULL THEN error('duckhts_depth_state: sample_id must not be NULL') "
    "ELSE CAST(sample_id AS VARCHAR) END AS sample_id, "
    "CASE WHEN chrom IS NULL THEN error('duckhts_depth_state: chrom must not be NULL') "
    "ELSE CAST(chrom AS VARCHAR) END AS chrom, "
    "CAST(start AS BIGINT) AS start, CAST(\"end\" AS BIGINT) AS \"end\", "
    "CAST(log2_ratio AS DOUBLE) AS log2_ratio, " DEPTH_AUTOSOMAL("chrom") " AS train "
    "FROM query_table(bins_table) WHERE coalesce(kept, false)), "
    "__dht_lists AS (SELECT sample_id, "
    "list(chrom ORDER BY chrom, start) AS chroms, "
    "list(start ORDER BY chrom, start) AS starts, "
    "list(\"end\" ORDER BY chrom, start) AS ends, "
    "list(log2_ratio ORDER BY chrom, start) AS ratios, "
    "list(train ORDER BY chrom, start) AS trains "
    "FROM __dht_bins GROUP BY sample_id), "
    "__dht_decoded AS (SELECT sample_id, unnest(__duckhts_depth_state_decode("
    "chroms, starts, ends, ratios, trains, "
    "CAST(normal_fractions AS DOUBLE[]), CAST(ploidies AS DOUBLE[]), CAST(max_cn AS INTEGER), "
    "CAST(txn_e AS DOUBLE), CAST(txn_strength AS DOUBLE), CAST(max_iter AS INTEGER), "
    "CAST(min_segment_bins AS INTEGER), CAST(alt_fraction_threshold AS DOUBLE), "
    "CAST(close_margin AS DOUBLE), CAST(max_bins AS BIGINT))) AS r FROM __dht_lists) ";

static const char *const depth_state_segments_sql[] = {
    "CREATE OR REPLACE MACRO duckhts_depth_state_segments(bins_table, " DEPTH_STATE_OPTIONS ") AS TABLE ",
    depth_state_lists,
    "SELECT sample_id, r.solution, r.start_normal, r.start_ploidy, r.normal_fraction, "
    "1 - r.normal_fraction AS tumour_fraction, r.ploidy, r.log_likelihood, r.altered_fraction, "
    "r.iterations, r.status, r.chrom, r.start, r.\"end\", r.bins, r.state, r.copy_number, "
    "r.median_log2_ratio, 'ichorCNA v0.6.1' AS upstream, 'duckhts_depth_state v1' AS method "
    "FROM __dht_decoded",
};

static const char *const depth_state_solutions_sql[] = {
    "CREATE OR REPLACE MACRO duckhts_depth_state_solutions(bins_table, " DEPTH_STATE_OPTIONS ") AS TABLE ",
    "SELECT sample_id, solution, start_normal, start_ploidy, normal_fraction, tumour_fraction, ploidy, "
    "log_likelihood, altered_fraction, iterations, status, count(*) AS segments, "
    "CAST(sum(bins) AS BIGINT) AS bins, upstream, method "
    "FROM duckhts_depth_state_segments(bins_table, " DEPTH_STATE_ARGUMENTS ") GROUP BY ALL",
};

bool register_duckhts_depth_sql(duckhts_registration_t *registration) {
    return duckhts_register_sql_parts(registration, depth_bins_sql,
                                      sizeof(depth_bins_sql) / sizeof(depth_bins_sql[0])) &&
           duckhts_register_sql_parts(registration, depth_state_segments_sql,
                                      sizeof(depth_state_segments_sql) / sizeof(depth_state_segments_sql[0])) &&
           duckhts_register_sql_parts(registration, depth_state_solutions_sql,
                                      sizeof(depth_state_solutions_sql) / sizeof(depth_state_solutions_sql[0])) &&
           duckhts_register_sql_parts(registration, depth_region_sql,
                                      sizeof(depth_region_sql) / sizeof(depth_region_sql[0])) &&
           duckhts_register_sql_parts(registration, depth_region_ref_sql,
                                      sizeof(depth_region_ref_sql) / sizeof(depth_region_ref_sql[0])) &&
           duckhts_register_sql_parts(registration, depth_sex_sql,
                                      sizeof(depth_sex_sql) / sizeof(depth_sex_sql[0])) &&
           duckhts_register_sql_parts(registration, depth_sex_ref_sql,
                                      sizeof(depth_sex_ref_sql) / sizeof(depth_sex_ref_sql[0]));
}
