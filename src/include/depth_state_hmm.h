/* depth_state_hmm.h -- copy-number state model over corrected bin depth, a
   port of the ichorCNA hidden Markov model (GavinHaLab/ichorCNA tag v0.6.1,
   R/EM.R, R/segmentation.R, R/utils.R and R/runIchorCNA.R) with the scaled
   forward-backward and Viterbi routines of HMMcopy 1.54.0 (src/forward_backward.c,
   src/viterbi.c). Method: Adalsteinsson, Ha, Freeman et al. 2017,
   https://doi.org/10.1038/s41467-017-00965-y. ichorCNA and HMMcopy are GPL-3.0;
   this is a re-implementation of the published model from the pinned sources
   under .sync/, with the upstream quirks named where they are kept.

   Scope: one sample, clonal states with copy numbers 1..max_cn (no
   homozygous-deletion state and no subclonal states), Student's t emission.
   The engine has no DuckDB dependency.

   Coordinates: bins arrive in the order the caller lists them, one contiguous
   block per contig (a contig changes when contig_id changes). Log ratios are
   log2(corrected ratio); NAN marks a bin without a ratio, which upstream
   carries with emission 1 in every state. */
#ifndef DUCKHTS_DEPTH_STATE_HMM_H
#define DUCKHTS_DEPTH_STATE_HMM_H

#include <stddef.h>
#include <stdint.h>

#define DUCKHTS_DEPTH_STATE_MAX_CN 20

typedef struct {
    int32_t max_cn;                 /* states are copy numbers 1..max_cn */
    double txn_e;                   /* self-transition probability (txnE) */
    double txn_strength;            /* transition pseudo-counts (txnStrength) */
    int32_t max_iter;               /* EM iteration cap (HMMsegment maxiter) */
    double lik_change;              /* relative change that ends EM */
    int32_t min_segment_bins;       /* guard: largest altered segment, in bins */
    double alt_fraction_threshold;  /* guard: altered share of training bins */
} duckhts_depth_state_params_t;

typedef struct {
    double start_normal;
    double start_ploidy;
    double normal;                  /* n: normal fraction; 1 after the guard */
    double ploidy;                  /* phi: tumour ploidy */
    double log_likelihood;          /* data log-likelihood plus log prior */
    double altered_fraction;        /* training bins whose state is not 2 */
    int32_t iterations;             /* EM iterations whose parameters are used */
    int32_t guard_zero;             /* the guard set the fraction to zero */
} duckhts_depth_state_solution_t;

typedef struct {
    size_t first_bin;               /* indices into the caller's bin arrays */
    size_t last_bin;
    int32_t bins;
    int32_t state;                  /* copy number of the Viterbi state */
    int32_t copy_number;            /* after the upstream integer correction */
    double median_log2_ratio;       /* median of the segment's finite ratios */
} duckhts_depth_state_segment_t;

typedef struct duckhts_depth_state duckhts_depth_state_t;

duckhts_depth_state_t *duckhts_depth_state_create(void);
void duckhts_depth_state_destroy(duckhts_depth_state_t *state);

/* Loads n_bins bins. train marks the bins EM learns from (upstream chrTrain:
 * the autosomes). Returns 0; -1 when storage cannot be allocated; -2 when
 * n_bins is 0 or above max_bins, or no training bin has a finite ratio. The
 * arrays are copied. */
int duckhts_depth_state_set_bins(duckhts_depth_state_t *state, size_t n_bins,
                                 const int32_t *contig_id, const int64_t *start,
                                 const int64_t *end, const double *log2_ratio,
                                 const uint8_t *train, size_t max_bins);

/* Fits one start (normal fraction start_normal, ploidy start_ploidy) with EM,
 * decodes the Viterbi path and builds the segments. Returns 0; -1 on
 * allocation failure; -2 on invalid parameters; -3 when a bin has zero
 * likelihood in every state (upstream stops with an error there too). */
int duckhts_depth_state_fit(duckhts_depth_state_t *state,
                            const duckhts_depth_state_params_t *params,
                            double start_normal, double start_ploidy,
                            duckhts_depth_state_solution_t *solution);

/* Segments of the last fit, valid until the next fit or set_bins. */
const duckhts_depth_state_segment_t *duckhts_depth_state_segments(
    const duckhts_depth_state_t *state, size_t *count);

#endif
