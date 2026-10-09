/* duckhts_count_error_fit(counts_table, ...): read error, contamination and
 * related quantities fitted from allele read counts, one row per sample of
 * the input, including a sample with no usable site (status no_sites).
 *
 * The relation is the one duckhts_roh_counts takes: sample_id, chrom, pos,
 * ref_count, alt_count and af. DuckDB builds the histogram of (block,
 * frequency bin, depth, alt count) cells itself and can spill that work; the
 * native __duckhts_count_cells aggregate then packs the cells of a sample
 * into one bounded buffer (count_cells.c), and __duckhts_count_error_fit fits
 * the model on it (count_error_udf.c). The guard arm returns no row and
 * checks the arguments even when the relation is empty.
 */
#include "count_sql.h"
#include "duckhts_registration.h"

#include <stdbool.h>
#include <stddef.h>

#define FIT_OPTIONS \
    "block_bases := 10000000, freq_bins := 64, max_depth := 1000, min_block_sites := 1000, " \
    "max_cells := 5000000, max_cell_bytes := 1073741824"

/* The errors of a malformed row, raised wherever the row is seen. */
#define FIT_POS_ERROR "error('pos must be a one-based position, at least 1')"
#define FIT_NEGATIVE_ERROR "error('read counts must not be negative')"
/* A position is cast only when it is a whole number. */
#define FIT_WHOLE_POS DUCKHTS_WHOLE_AS("pos", "BIGINT", "pos must be a whole number")

/* The macro's arguments as the native functions take them: whole numbers,
 * checked before the cast so that 1.5 is an error and not 2. */
#define FIT_ARGUMENTS \
    DUCKHTS_WHOLE_OPTION("block_bases") ", " DUCKHTS_WHOLE_OPTION("freq_bins") ", " \
    DUCKHTS_WHOLE_OPTION("max_depth") ", " DUCKHTS_WHOLE_OPTION("min_block_sites") ", " \
    DUCKHTS_WHOLE_OPTION("max_cells") ", " DUCKHTS_WHOLE_OPTION("max_cell_bytes")

static const char *const count_error_fit_sql[] = {
    "CREATE OR REPLACE MACRO duckhts_count_error_fit(counts_table, " FIT_OPTIONS ") AS TABLE (",
    /* Rows as the decode reads them, with the counts kept wide until the
     * depth test, so that a count beyond INTEGER makes a skipped site and
     * not a conversion error. A negative count is an error. */
    "WITH __dht_counts AS (SELECT CAST(sample_id AS VARCHAR) AS smp, CAST(chrom AS VARCHAR) AS chrom, "
    FIT_WHOLE_POS " AS pos, " DUCKHTS_WHOLE_COUNT_AS("ref_count", "BIGINT") " AS ref_count, "
    DUCKHTS_WHOLE_COUNT_AS("alt_count", "BIGINT") " AS alt_count, CAST(af AS DOUBLE) AS af "
    "FROM query_table(counts_table)), ",
    /* pos is one-based, so block k holds k * block_bases + 1 to (k + 1) * block_bases. */
    "__dht_sites AS (SELECT smp, chrom, CASE WHEN pos IS NULL OR pos < 1 THEN " FIT_POS_ERROR " "
    "ELSE (pos - 1) // " DUCKHTS_WHOLE_OPTION("block_bases") " END AS block_index, af, "
    /* A count above max_depth makes a skipped site before any sum, so the
     * sum of two counts at most max_depth cannot overflow. */
    "CASE WHEN ref_count < 0 OR alt_count < 0 THEN " FIT_NEGATIVE_ERROR " "
    "WHEN ref_count > " DUCKHTS_WHOLE_OPTION("max_depth") " OR alt_count > " DUCKHTS_WHOLE_OPTION("max_depth") " "
    "THEN NULL ELSE ref_count + alt_count END AS depth, alt_count AS alt "
    "FROM __dht_counts WHERE ref_count IS NOT NULL AND alt_count IS NOT NULL "
    "AND af IS NOT NULL AND isfinite(af) AND af > 0 AND af < 1), ",
    /* The histogram: one row per cell, with the mean frequency of its sites.
     * Sites deeper than max_depth are left out, which bounds the cells. */
    "__dht_cells AS (SELECT smp, chrom, block_index, "
    "CAST(floor(af * " DUCKHTS_WHOLE_OPTION("freq_bins") ") AS INTEGER) AS bin, "
    "depth, alt, count(*) AS sites, avg(af) AS af FROM __dht_sites "
    "WHERE depth >= 1 AND depth <= " DUCKHTS_WHOLE_OPTION("max_depth") " GROUP BY ALL), ",
    /* Blocks of a sample are numbered from 0 in chromosome and position order. */
    "__dht_blocks AS (SELECT *, dense_rank() OVER (PARTITION BY smp ORDER BY chrom, block_index) - 1 AS block "
    "FROM __dht_cells), ",
    "__dht_hist AS (SELECT smp, __duckhts_count_cells(CAST(block AS INTEGER), af, CAST(depth AS INTEGER), "
    "CAST(alt AS INTEGER), sites, " DUCKHTS_WHOLE_OPTION("max_cells") ", "
    DUCKHTS_WHOLE_OPTION("max_cell_bytes") ") AS cells "
    "FROM __dht_blocks GROUP BY smp), ",
    /* Every sample of the input keeps its row: one with no usable site gets a
     * NULL histogram, which the fit reports as no_sites. The sample set reads
     * the input again rather than __dht_counts, so that DuckDB streams the
     * rows once into the histogram instead of materialising the CTE it
     * would otherwise reference twice; the same scan checks every row for a
     * bad position, a negative count or a fractional count, rows the
     * frequency and count filters of __dht_sites would otherwise drop before
     * their own checks ran. */
    "__dht_samples AS (SELECT CASE WHEN bad_pos > 0 THEN " FIT_POS_ERROR " "
    "WHEN bad_counts > 0 THEN " FIT_NEGATIVE_ERROR " ELSE smp END AS smp FROM ("
    "SELECT CAST(sample_id AS VARCHAR) AS smp, "
    "count(*) FILTER (WHERE " FIT_WHOLE_POS " IS NULL OR " FIT_WHOLE_POS " < 1) AS bad_pos, "
    "count(*) FILTER (WHERE " DUCKHTS_WHOLE_COUNT_AS("ref_count", "BIGINT") " < 0 OR "
    DUCKHTS_WHOLE_COUNT_AS("alt_count", "BIGINT") " < 0) AS bad_counts "
    "FROM query_table(counts_table) GROUP BY smp)), ",
    /* The join is NULL-safe: a NULL sample_id is one sample. */
    "__dht_fit AS (SELECT s.smp, __duckhts_count_error_fit(h.cells, " DUCKHTS_WHOLE_OPTION("min_block_sites") ", "
    DUCKHTS_WHOLE_OPTION("max_cell_bytes") ") AS f "
    "FROM __dht_samples s LEFT JOIN __dht_hist h ON s.smp IS NOT DISTINCT FROM h.smp) ",
    "SELECT smp AS \"sample\", f.sites, f.reads, f.mean_depth, f.blocks, f.seq_error, f.contamination, "
    "f.contamination_sd, f.homozygosity_excess, f.allele_balance, f.spread_hom, f.spread_het, "
    "f.artefact_weight, f.log_likelihood, f.contamination_relative, f.log_likelihood_relative, "
    "f.status, f.method FROM __dht_fit ",
    "UNION ALL SELECT NULL::VARCHAR, NULL::BIGINT, NULL::BIGINT, NULL::DOUBLE, NULL::INTEGER, NULL::DOUBLE, "
    "NULL::DOUBLE, NULL::DOUBLE, NULL::DOUBLE, NULL::DOUBLE, NULL::DOUBLE, NULL::DOUBLE, NULL::DOUBLE, "
    "NULL::DOUBLE, NULL::DOUBLE, NULL::DOUBLE, NULL::VARCHAR, NULL::VARCHAR "
    "WHERE NOT __duckhts_count_error_valid_args(" FIT_ARGUMENTS ")",
    ")",
};

bool register_duckhts_count_error_sql(duckhts_registration_t *registration) {
    return duckhts_register_sql_parts(registration, count_error_fit_sql,
                                      sizeof(count_error_fit_sql) / sizeof(count_error_fit_sql[0]));
}
