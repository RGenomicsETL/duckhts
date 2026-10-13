/* Low-pass first and second fragment moments. No DuckDB dependency.
 * Owned block arrays belong to the caller; the fitting workspace is fixed-size.
 * Frequencies describe the allele counted as alt. */
#ifndef DUCKHTS_COUNT_MIXTURE_MODEL_H
#define DUCKHTS_COUNT_MIXTURE_MODEL_H
#include <stddef.h>
#include <stdint.h>

enum {
    DUCKHTS_MIXTURE_READ_BASIS = 3,
    DUCKHTS_MIXTURE_PAIR_BASIS = 6,
    DUCKHTS_MIXTURE_SCORE_DIM = 9
};
enum {
    MIX_SITES,
    MIX_USABLE,
    MIX_ZERO,
    MIX_SINGLETON,
    MIX_READS,
    MIX_OTHER,
    MIX_PAIR_SITES,
    MIX_EXCLUDED,
    MIX_RR,
    MIX_RA,
    MIX_AA,
    MIX_MISSING,
    MIX_COUNTERS
};
typedef struct {
    int64_t id;
    uint64_t counts[MIX_COUNTERS];
    double s[3][3], y[3], t[6][6], z[6];
    double denominator;
    int matched;
} duckhts_mixture_block_t;
typedef struct {
    double seq_error, artefact_weight, homozygosity_excess;
    int min_blocks;
} duckhts_mixture_calibration_t;
typedef struct {
    double share, lower, upper, se, singles_share, singles_se;
    double pair_residual, pair_residual_se, score, leverage;
    int rank, blocks, pair_blocks, failed_refits, at_bound, interval_hull;
    const char *status, *diagnostic;
} duckhts_mixture_fit_t;

/* Exact non-negative binary64 moment sums. Nineteen 64-bit limbs cover
 * subnormals through UINT64_MAX operands <= 2^31, without retained rows.
 * Finalization rounds once to nearest, ties to even. Add/merge return zero
 * for invalid operands or capacity overflow. Initialize the state to zero. */
#define DUCKHTS_MIXTURE_SUM_LIMBS 19
#define DUCKHTS_MIXTURE_SUM_MAX_OPERAND 2147483648.0
typedef struct {
    uint64_t limbs[DUCKHTS_MIXTURE_SUM_LIMBS];
} duckhts_mixture_sum_t;
int duckhts_mixture_sum_add(duckhts_mixture_sum_t *sum, double value);
int duckhts_mixture_sum_merge(duckhts_mixture_sum_t *sum, const duckhts_mixture_sum_t *from);
double duckhts_mixture_sum_value(const duckhts_mixture_sum_t *sum);

/* Returns zero on count overflow. Inputs are validated by the admission layer.
 * Missing sites and over-depth sites contribute QC counts but no moments. */
int duckhts_mixture_add(duckhts_mixture_block_t *block, uint64_t ref, uint64_t alt, uint64_t other,
                        double f, double g, int usable, int excluded);
int duckhts_mixture_merge(duckhts_mixture_block_t *to, const duckhts_mixture_block_t *from);
void duckhts_mixture_coefficients(double share, const duckhts_mixture_calibration_t *cal, int ibd1,
                                  double beta[3], double gamma[6]);
void duckhts_mixture_fit(const duckhts_mixture_block_t *blocks, size_t count,
                         const duckhts_mixture_calibration_t *cal, int ibd1,
                         duckhts_mixture_fit_t *fit);
#endif
