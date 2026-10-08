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

/* The macro's arguments as the native functions take them. */
#define FIT_ARGUMENTS \
    "CAST(block_bases AS BIGINT), CAST(freq_bins AS BIGINT), CAST(max_depth AS BIGINT), " \
    "CAST(min_block_sites AS BIGINT), CAST(max_cells AS BIGINT), CAST(max_cell_bytes AS BIGINT)"

static const char *const count_error_fit_sql[] = {
    "CREATE OR REPLACE MACRO duckhts_count_error_fit(counts_table, " FIT_OPTIONS ") AS TABLE (",
    /* Rows as the decode reads them. A negative count is an error. */
    "WITH __dht_counts AS (SELECT CAST(sample_id AS VARCHAR) AS smp, CAST(chrom AS VARCHAR) AS chrom, "
    "CAST(pos AS BIGINT) AS pos, " DUCKHTS_WHOLE_COUNT("ref_count") " AS ref_count, "
    DUCKHTS_WHOLE_COUNT("alt_count") " AS alt_count, CAST(af AS DOUBLE) AS af "
    "FROM query_table(counts_table)), ",
    /* pos is one-based, so block k holds k * block_bases + 1 to (k + 1) * block_bases. */
    "__dht_sites AS (SELECT smp, chrom, CASE WHEN pos IS NULL OR pos < 1 "
    "THEN error('pos must be a one-based position, at least 1') "
    "ELSE (pos - 1) // CAST(block_bases AS BIGINT) END AS block_index, af, "
    "CASE WHEN ref_count < 0 OR alt_count < 0 THEN error('read counts must not be negative') "
    "ELSE CAST(ref_count AS BIGINT) + CAST(alt_count AS BIGINT) END AS depth, alt_count AS alt "
    "FROM __dht_counts WHERE ref_count IS NOT NULL AND alt_count IS NOT NULL "
    "AND af IS NOT NULL AND isfinite(af) AND af > 0 AND af < 1), ",
    /* The histogram: one row per cell, with the mean frequency of its sites.
     * Sites deeper than max_depth are left out, which bounds the cells. */
    "__dht_cells AS (SELECT smp, chrom, block_index, CAST(floor(af * freq_bins) AS INTEGER) AS bin, "
    "depth, alt, count(*) AS sites, avg(af) AS af FROM __dht_sites "
    "WHERE depth >= 1 AND depth <= CAST(max_depth AS BIGINT) GROUP BY ALL), ",
    /* Blocks of a sample are numbered from 0 in chromosome and position order. */
    "__dht_blocks AS (SELECT *, dense_rank() OVER (PARTITION BY smp ORDER BY chrom, block_index) - 1 AS block "
    "FROM __dht_cells), ",
    "__dht_hist AS (SELECT smp, __duckhts_count_cells(CAST(block AS INTEGER), af, CAST(depth AS INTEGER), "
    "CAST(alt AS INTEGER), sites, CAST(max_cells AS BIGINT), CAST(max_cell_bytes AS BIGINT)) AS cells "
    "FROM __dht_blocks GROUP BY smp), ",
    /* Every sample of the input keeps its row: one with no usable site gets a
     * NULL histogram, which the fit reports as no_sites. The sample set reads
     * the input again rather than __dht_counts, so that DuckDB streams the
     * rows once into the histogram instead of materialising the CTE it
     * would otherwise reference twice. */
    "__dht_samples AS (SELECT DISTINCT CAST(sample_id AS VARCHAR) AS smp FROM query_table(counts_table)), "
    "__dht_fit AS (SELECT smp, __duckhts_count_error_fit(cells, CAST(min_block_sites AS BIGINT), "
    "CAST(max_cell_bytes AS BIGINT)) AS f "
    "FROM __dht_samples LEFT JOIN __dht_hist USING (smp)) ",
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
