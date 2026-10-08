/* Maximum-likelihood fit of read error, contamination and related quantities
 * from allele read counts at sites of known population frequency.
 *
 * The model of one sample. A site has `depth` reads, `alt` of which show the
 * allele whose population frequency is `af`.
 *   sample genotype g in 0, 1, 2 copies: (1 - F) * HWE(af) + F * (1 - af, 0, af)
 *   second genome h in 0, 1, 2 copies: HWE(af) when unrelated; one allele of g
 *     and one draw at af when it is a parent or child of the sample
 *   P(read shows the allele | g, h) = (1 - c) * q[g] + c * q[h],
 *     q = (e, b, 1 - e)
 *   alt ~ beta-binomial(depth, p, rho): rho_hom at g in {0, 2}, rho_het at g = 1;
 *     rho = 0 is the binomial
 *   with probability w the site is an artefact whose alt count is uniform on
 *     0..depth, whatever its genotype
 * Parameters: e read error, c contamination, F excess homozygosity (a
 * nuisance: it also absorbs a mismatch between the sample's population and
 * the frequencies), b allele balance at heterozygotes, rho_hom and rho_het
 * spreads, w artefact weight.
 *
 * The likelihood depends on a site only through (depth, alt, af), so the fit
 * runs on histogram cells (count_cells.h). This header has no DuckDB
 * dependency.
 */
#ifndef DUCKHTS_COUNT_ERROR_MODEL_H
#define DUCKHTS_COUNT_ERROR_MODEL_H

#include <stddef.h>
#include <stdint.h>

#include "count_cells.h"

enum {
    DUCKHTS_COUNT_ERROR_PARAMETERS = 7,
    DUCKHTS_COUNT_ERROR_METHOD_VERSION = 1
};

typedef struct {
    double seq_error;
    double contamination;
    double homozygosity_excess;
    double allele_balance;
    double spread_hom;
    double spread_het;
    double artefact_weight;
} duckhts_count_error_params_t;

/* The second genome's relation to the sample. */
typedef enum {
    DUCKHTS_COUNT_ERROR_UNRELATED = 0,
    DUCKHTS_COUNT_ERROR_RELATIVE = 1
} duckhts_count_error_relation_t;

typedef struct {
    duckhts_count_error_params_t params;
    double log_likelihood;
    int converged;     /* the optimizer stopped on its tolerance, not its step limit */
    int at_bound;      /* contamination or the artefact weight is at its upper bound */
} duckhts_count_error_fit_t;

/* Both fits charge their working memory against max_cell_bytes
 * (duckhts_count_bytes_charge) and release it before returning: 8 bytes per
 * cell, 88 bytes per distinct (depth, alt) pair and 24 bytes per depth up to
 * the deepest. On a limit or an allocation failure they return 0 and write
 * the reason to `error`. */

/* Fits all parameters by maximum likelihood from a fixed set of starts and
 * returns the best. */
int duckhts_count_error_fit(const duckhts_count_cell_t *cells, size_t count,
                            duckhts_count_error_relation_t relation, uint64_t max_cell_bytes,
                            duckhts_count_error_fit_t *fit, char *error, size_t error_length);

/* Fits contamination only, the other parameters held at `fixed`, from `fixed`
 * as the start. For the per-block fits. */
int duckhts_count_error_fit_block(const duckhts_count_cell_t *cells, size_t count,
                                  const duckhts_count_error_params_t *fixed, uint64_t max_cell_bytes,
                                  duckhts_count_error_fit_t *fit, char *error, size_t error_length);

#endif
