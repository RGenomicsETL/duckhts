/* Copy-number state model over corrected bin depth: see depth_state_hmm.h.
 *
 * The port follows ichorCNA v0.6.1 for one sample with the Student's t
 * emission and clonal states 1..max_cn (run_ichorCNA defaults: includeHOMD
 * FALSE, scStates NULL, likModel "t"). Upstream behaviour that is kept on
 * purpose, with the source line in R/EM.R, R/segmentation.R or R/utils.R:
 *   - the precision prior is set from the log2 ratios (getDefaultParameters is
 *     called before HMMsegment converts the data to natural log), while the
 *     emission sees natural-log ratios;
 *   - the initial-distribution update uses the first training bin only
 *     (estimateMixWeightsParamMap, rho[, 1]);
 *   - the transition prior term evaluates dirichletpdflog(A_prior, A), with the
 *     arguments in that order (priorProbs);
 *   - when the objective decreases, the previous iteration's n, phi and loglik
 *     are reported, but the Viterbi path uses the emissions and transitions of
 *     the rejected iteration (runEM, "Using previous iteration").
 * The optimizer differs: upstream calls L-BFGS-B on each parameter in turn;
 * this file climbs each one-dimensional objective from the same start with a
 * bounded adaptive step. Both find a local maximum near the start. */
#include "depth_state_hmm.h"

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define LN2 0.6931471805599453
#define DEPTH_PI 3.14159265358979323846
#define NU 2.1
#define ALPHA_LAMBDA 2.0
#define ALPHA_PHI 4.0
#define BETA_PHI 0.75
#define ALPHA_N 1.0
#define BETA_N 1.0
/* betapdflog(0.5, 2, 2): the subclone prevalence stays at its start without
 * subclonal states, and its prior is a constant of the objective. */
#define SP_PRIOR_LOG 0.4054651081081644
#define N_LO 1e-6
#define N_HI (1.0 - 1e-6)
#define PHI_LO DBL_EPSILON
#define PHI_HI 10.0
#define LAMBDA_LO 1e-5
#define MIN_PURITY_TO_CORRECT 0.05
#define CN_FLOOR (1.0 / 64.0)
#define CLIMB_EVALUATIONS 400

struct duckhts_depth_state {
    size_t n_bins;
    size_t capacity;
    int32_t *contig;
    int64_t *start;
    int64_t *end;
    double *x2;              /* log2 ratio, NAN when absent */
    double *x;               /* natural-log ratio, NAN when absent */
    uint8_t *train;
    size_t n_contigs;
    size_t *contig_begin;    /* n_contigs + 1 entries */
    size_t contig_capacity;
    size_t *index;           /* bin indices of one contig's training bins */
    int32_t K;
    size_t work_capacity;    /* K * capacity */
    double *py;              /* emission, k + K * t */
    double *rho;             /* responsibilities, training bins only */
    double *delta;           /* forward scratch, then Viterbi scores */
    int32_t *psi;
    int32_t *path;           /* Viterbi state per bin, -1 before a fit */
    double *median_scratch;
    duckhts_depth_state_segment_t *segments;
    size_t n_segments;
    size_t segment_capacity;
};

static double gammapdflog(double x, double a, double b) {
    return a * log(b) - lgamma(a) + (a - 1.0) * log(x) - b * x;
}

static double betapdflog(double x, double a, double b) {
    return -(lgamma(a) + lgamma(b) - lgamma(a + b)) + (a - 1.0) * log(x) + (b - 1.0) * log1p(-x);
}

/* dirichletpdflog(x, k) as upstream writes it: lgamma(sum k) - sum lgamma(k)
 * + sum (k - 1) log x. */
static double dirichletpdflog(const double *x, const double *k, int32_t K) {
    double sum = 0.0, lg = 0.0, l = 0.0;
    for (int32_t i = 0; i < K; i++) {
        sum += k[i];
        lg += lgamma(k[i]);
        l += (k[i] - 1.0) * log(x[i]);
    }
    return lgamma(sum) - lg + l;
}

static double t_constant(void) {
    return tgamma(NU / 2.0 + 0.5) / tgamma(NU / 2.0);
}

/* tdistPDF; an absent ratio has density 1 (upstream: p[is.na(p)] <- 1). */
static double tdist_pdf(double x, double mu, double lambda) {
    if (isnan(x)) return 1.0;
    double d = x - mu;
    return t_constant() * sqrt(lambda / (DEPTH_PI * NU)) * pow(1.0 + lambda * d * d / NU, -0.5 * NU - 0.5);
}

static double log_tdist_pdf(double x, double mu, double lambda) {
    if (isnan(x)) return 0.0;
    double d = x - mu;
    return log(t_constant()) + 0.5 * log(lambda / (DEPTH_PI * NU)) - (0.5 * NU + 0.5) * log1p(lambda * d * d / NU);
}

/* get2ComponentMixture: log(((1 - n) c + 2 n) / ((1 - n) phi + 2 n)). */
static double state_mean(int32_t copy_number, double n, double phi) {
    return log(((1.0 - n) * copy_number + 2.0 * n) / ((1.0 - n) * phi + 2.0 * n));
}

duckhts_depth_state_t *duckhts_depth_state_create(void) {
    duckhts_depth_state_t *state = calloc(1, sizeof(*state));
    return state;
}

void duckhts_depth_state_destroy(duckhts_depth_state_t *state) {
    if (!state) return;
    free(state->contig);
    free(state->start);
    free(state->end);
    free(state->x2);
    free(state->x);
    free(state->train);
    free(state->contig_begin);
    free(state->index);
    free(state->py);
    free(state->rho);
    free(state->delta);
    free(state->psi);
    free(state->path);
    free(state->median_scratch);
    free(state->segments);
    free(state);
}

static int reserve_bins(duckhts_depth_state_t *state, size_t n_bins) {
    if (n_bins <= state->capacity) return 0;
    if (n_bins > SIZE_MAX / sizeof(double)) return -1;
    int32_t *contig = realloc(state->contig, n_bins * sizeof(*contig));
    if (!contig) return -1;
    state->contig = contig;
    int64_t *start = realloc(state->start, n_bins * sizeof(*start));
    if (!start) return -1;
    state->start = start;
    int64_t *end = realloc(state->end, n_bins * sizeof(*end));
    if (!end) return -1;
    state->end = end;
    double *x2 = realloc(state->x2, n_bins * sizeof(*x2));
    if (!x2) return -1;
    state->x2 = x2;
    double *x = realloc(state->x, n_bins * sizeof(*x));
    if (!x) return -1;
    state->x = x;
    uint8_t *train = realloc(state->train, n_bins * sizeof(*train));
    if (!train) return -1;
    state->train = train;
    size_t *index = realloc(state->index, n_bins * sizeof(*index));
    if (!index) return -1;
    state->index = index;
    int32_t *path = realloc(state->path, n_bins * sizeof(*path));
    if (!path) return -1;
    state->path = path;
    double *median_scratch = realloc(state->median_scratch, n_bins * sizeof(*median_scratch));
    if (!median_scratch) return -1;
    state->median_scratch = median_scratch;
    state->capacity = n_bins;
    return 0;
}

static int reserve_work(duckhts_depth_state_t *state, int32_t K) {
    size_t need = (size_t)K * state->n_bins;
    if (need <= state->work_capacity) return 0;
    if (need > SIZE_MAX / sizeof(double)) return -1;
    double *py = realloc(state->py, need * sizeof(*py));
    if (!py) return -1;
    state->py = py;
    double *rho = realloc(state->rho, need * sizeof(*rho));
    if (!rho) return -1;
    state->rho = rho;
    double *delta = realloc(state->delta, need * sizeof(*delta));
    if (!delta) return -1;
    state->delta = delta;
    int32_t *psi = realloc(state->psi, need * sizeof(*psi));
    if (!psi) return -1;
    state->psi = psi;
    state->work_capacity = need;
    return 0;
}

static int reserve_contigs(duckhts_depth_state_t *state, size_t n_contigs) {
    if (n_contigs + 1 <= state->contig_capacity) return 0;
    size_t *begin = realloc(state->contig_begin, (n_contigs + 1) * sizeof(*begin));
    if (!begin) return -1;
    state->contig_begin = begin;
    state->contig_capacity = n_contigs + 1;
    return 0;
}

int duckhts_depth_state_set_bins(duckhts_depth_state_t *state, size_t n_bins,
                                 const int32_t *contig_id, const int64_t *start,
                                 const int64_t *end, const double *log2_ratio,
                                 const uint8_t *train, size_t max_bins) {
    if (n_bins == 0 || n_bins > max_bins) return -2;
    if (reserve_bins(state, n_bins) < 0) return -1;
    size_t n_contigs = 0;
    size_t finite_train = 0;
    for (size_t t = 0; t < n_bins; t++) {
        if (t == 0 || contig_id[t] != contig_id[t - 1]) n_contigs++;
        state->contig[t] = contig_id[t];
        state->start[t] = start[t];
        state->end[t] = end[t];
        state->x2[t] = isfinite(log2_ratio[t]) ? log2_ratio[t] : NAN;
        state->x[t] = isfinite(log2_ratio[t]) ? log2_ratio[t] * LN2 : NAN;
        state->train[t] = train[t] ? 1u : 0u;
        if (state->train[t] && isfinite(log2_ratio[t])) finite_train++;
        state->path[t] = -1;
    }
    if (finite_train < 2) return -2;
    if (reserve_contigs(state, n_contigs) < 0) return -1;
    size_t c = 0;
    for (size_t t = 0; t < n_bins; t++) {
        if (t == 0 || contig_id[t] != contig_id[t - 1]) state->contig_begin[c++] = t;
    }
    state->contig_begin[c] = n_bins;
    state->n_contigs = n_contigs;
    state->n_bins = n_bins;
    state->n_segments = 0;
    return 0;
}

/* ---- model parameters of one fit ---------------------------------------- */

typedef struct {
    int32_t K;
    double n, phi;
    double lambda[DUCKHTS_DEPTH_STATE_MAX_CN];
    double pi[DUCKHTS_DEPTH_STATE_MAX_CN];
} fit_point_t;

typedef struct {
    int32_t K;
    double A[DUCKHTS_DEPTH_STATE_MAX_CN * DUCKHTS_DEPTH_STATE_MAX_CN];
    double A_prior[DUCKHTS_DEPTH_STATE_MAX_CN * DUCKHTS_DEPTH_STATE_MAX_CN];
    double dir_prior[DUCKHTS_DEPTH_STATE_MAX_CN * DUCKHTS_DEPTH_STATE_MAX_CN];
    double kappa[DUCKHTS_DEPTH_STATE_MAX_CN];
    double beta_lambda[DUCKHTS_DEPTH_STATE_MAX_CN];
    double Z[DUCKHTS_DEPTH_STATE_MAX_CN * DUCKHTS_DEPTH_STATE_MAX_CN];
    double mu[DUCKHTS_DEPTH_STATE_MAX_CN];
    size_t first_train;      /* the first training bin (rho[, 1] upstream) */
} fit_model_t;

#define AT(m, i, j, K) ((m)[(i) + (size_t)(K) * (j)])

static void normalize_row(double *row, int32_t K) {
    double sum = 0.0;
    for (int32_t j = 0; j < K; j++) sum += row[j];
    if (sum == 0.0) sum = 1.0;
    for (int32_t j = 0; j < K; j++) row[j] /= sum;
}

static void set_means(fit_model_t *model, double n, double phi) {
    for (int32_t k = 0; k < model->K; k++) model->mu[k] = state_mean(k + 1, n, phi);
}

static void compute_emissions(duckhts_depth_state_t *state, const fit_model_t *model,
                              const double *lambda) {
    int32_t K = model->K;
    for (size_t t = 0; t < state->n_bins; t++) {
        for (int32_t k = 0; k < K; k++) {
            AT(state->py, k, t, K) = tdist_pdf(state->x[t], model->mu[k], lambda[k]);
        }
    }
}

/* HMMcopy forward_backward over the training bins of one contig, whose bin
 * indices are index[0..T). Writes rho for those bins, adds the log-likelihood
 * and the expected transition counts. Returns -3 when a bin has no state with
 * positive likelihood (normalizeInPlace stops there). */
static int forward_backward(duckhts_depth_state_t *state, const fit_model_t *model,
                            const double *pi, const size_t *index, size_t T,
                            double *loglik, double *Z) {
    int32_t K = model->K;
    double *alpha = state->delta;
    double m[DUCKHTS_DEPTH_STATE_MAX_CN], b[DUCKHTS_DEPTH_STATE_MAX_CN];
    double beta[DUCKHTS_DEPTH_STATE_MAX_CN], beta_next[DUCKHTS_DEPTH_STATE_MAX_CN];
    double eta[DUCKHTS_DEPTH_STATE_MAX_CN * DUCKHTS_DEPTH_STATE_MAX_CN];
    for (size_t t = 0; t < T; t++) {
        const double *py = state->py + K * index[t];
        double *a = alpha + K * t;
        if (t == 0) {
            for (int32_t k = 0; k < K; k++) a[k] = pi[k] * py[k];
        } else {
            const double *prev = alpha + K * (t - 1);
            for (int32_t j = 0; j < K; j++) {
                double s = 0.0;
                for (int32_t i = 0; i < K; i++) s += AT(model->A, i, j, K) * prev[i];
                a[j] = s * py[j];
            }
        }
        double scale = 0.0;
        for (int32_t k = 0; k < K; k++) scale += a[k];
        if (!(scale > 0.0)) return -3;
        for (int32_t k = 0; k < K; k++) a[k] /= scale;
        *loglik += log(scale);
    }
    for (int32_t k = 0; k < K; k++) {
        beta_next[k] = 1.0;
        AT(state->rho, k, index[T - 1], K) = alpha[K * (T - 1) + k];
    }
    for (size_t t = T - 1; t-- > 0;) {
        const double *py_next = state->py + K * index[t + 1];
        const double *a = alpha + K * t;
        double sum = 0.0;
        for (int32_t k = 0; k < K; k++) b[k] = beta_next[k] * py_next[k];
        for (int32_t i = 0; i < K; i++) {
            double s = 0.0;
            for (int32_t j = 0; j < K; j++) s += AT(model->A, i, j, K) * b[j];
            m[i] = s;
            sum += s;
        }
        if (!(sum > 0.0)) return -3;
        for (int32_t i = 0; i < K; i++) beta[i] = m[i] / sum;
        double eta_sum = 0.0;
        for (int32_t i = 0; i < K; i++) {
            for (int32_t j = 0; j < K; j++) {
                double v = AT(model->A, i, j, K) * a[i] * b[j];
                AT(eta, i, j, K) = v;
                eta_sum += v;
            }
        }
        if (!(eta_sum > 0.0)) return -3;
        for (int32_t i = 0; i < K * K; i++) Z[i] += eta[i] / eta_sum;
        double g_sum = 0.0;
        for (int32_t k = 0; k < K; k++) {
            m[k] = a[k] * beta[k];
            g_sum += m[k];
        }
        if (!(g_sum > 0.0)) return -3;
        for (int32_t k = 0; k < K; k++) AT(state->rho, k, index[t], K) = m[k] / g_sum;
        memcpy(beta_next, beta, (size_t)K * sizeof(double));
    }
    return 0;
}

/* HMMcopy viterbi over all bins of one contig, with log pi, log A and log py.
 * Returns -3 when no state has a finite score at some bin. */
static int viterbi(duckhts_depth_state_t *state, const fit_model_t *model,
                   const double *log_pi, const double *log_A, size_t begin, size_t end) {
    int32_t K = model->K;
    double *delta = state->delta;
    int32_t *psi = state->psi;
    for (size_t t = begin; t < end; t++) {
        const double *py = state->py + K * t;
        for (int32_t j = 0; j < K; j++) {
            double best;
            int32_t arg = 0;
            if (t == begin) {
                best = log_pi[j];
            } else {
                best = -INFINITY;
                for (int32_t i = 0; i < K; i++) {
                    double v = delta[K * (t - 1) + i] + AT(log_A, i, j, K);
                    if (v > best) {
                        best = v;
                        arg = i;
                    }
                }
            }
            delta[K * t + j] = best + log(py[j]);
            psi[K * t + j] = arg;
        }
        bool finite = false;
        for (int32_t j = 0; j < K; j++) finite = finite || delta[K * t + j] > -INFINITY;
        if (!finite) return -3;
    }
    double best = -INFINITY;
    int32_t arg = 0;
    for (int32_t j = 0; j < K; j++) {
        if (delta[K * (end - 1) + j] > best) {
            best = delta[K * (end - 1) + j];
            arg = j;
        }
    }
    state->path[end - 1] = arg;
    for (size_t t = end - 1; t-- > begin;) {
        state->path[t] = psi[K * (t + 1) + state->path[t + 1]];
    }
    return 0;
}

/* ---- M-step objectives ---------------------------------------------------- */

typedef struct {
    const duckhts_depth_state_t *state;
    const fit_model_t *model;
    double n, phi;
    const double *lambda;
    int32_t k;               /* the state whose precision moves */
} objective_ctx_t;

static double weighted_t_sum(const objective_ctx_t *ctx, const double *mu, const double *lambda,
                             int32_t only_k) {
    const duckhts_depth_state_t *state = ctx->state;
    int32_t K = ctx->model->K;
    double sum = 0.0;
    for (size_t t = 0; t < state->n_bins; t++) {
        if (!state->train[t]) continue;
        if (only_k >= 0) {
            sum += AT(state->rho, only_k, t, K) * log_tdist_pdf(state->x[t], mu[only_k], lambda[only_k]);
        } else {
            for (int32_t k = 0; k < K; k++) {
                sum += AT(state->rho, k, t, K) * log_tdist_pdf(state->x[t], mu[k], lambda[k]);
            }
        }
    }
    return sum;
}

static double objective_n(double n, void *pointer) {
    const objective_ctx_t *ctx = pointer;
    double mu[DUCKHTS_DEPTH_STATE_MAX_CN];
    for (int32_t k = 0; k < ctx->model->K; k++) mu[k] = state_mean(k + 1, n, ctx->phi);
    return weighted_t_sum(ctx, mu, ctx->lambda, -1) + betapdflog(n, ALPHA_N, BETA_N);
}

static double objective_phi(double phi, void *pointer) {
    const objective_ctx_t *ctx = pointer;
    double mu[DUCKHTS_DEPTH_STATE_MAX_CN];
    for (int32_t k = 0; k < ctx->model->K; k++) mu[k] = state_mean(k + 1, ctx->n, phi);
    return weighted_t_sum(ctx, mu, ctx->lambda, -1) + gammapdflog(phi, ALPHA_PHI, BETA_PHI);
}

static double objective_lambda(double lambda, void *pointer) {
    const objective_ctx_t *ctx = pointer;
    double lambdas[DUCKHTS_DEPTH_STATE_MAX_CN];
    memcpy(lambdas, ctx->lambda, (size_t)ctx->model->K * sizeof(double));
    lambdas[ctx->k] = lambda;
    return weighted_t_sum(ctx, ctx->model->mu, lambdas, ctx->k) +
           gammapdflog(lambda, ALPHA_LAMBDA, ctx->model->beta_lambda[ctx->k]);
}

/* A bounded climb from x0: steps that improve the objective double, steps
 * that do not halve, until the step is negligible. */
static double climb(double (*f)(double, void *), void *ctx, double x0, double lo, double hi,
                    double step) {
    double x = x0;
    double fx = f(x, ctx);
    for (int evaluations = 0; evaluations < CLIMB_EVALUATIONS; evaluations += 2) {
        if (step <= 1e-10 * (1.0 + fabs(x))) break;
        double up = fmin(x + step, hi), down = fmax(x - step, lo);
        double fu = up != x ? f(up, ctx) : -INFINITY;
        double fd = down != x ? f(down, ctx) : -INFINITY;
        if (fu > fx && fu >= fd) {
            x = up;
            fx = fu;
            step *= 2.0;
        } else if (fd > fx) {
            x = down;
            fx = fd;
            step *= 2.0;
        } else {
            step *= 0.5;
        }
    }
    return x;
}

static double prior_log(const fit_model_t *model, const fit_point_t *point) {
    int32_t K = model->K;
    double prior = 0.0;
    for (int32_t k = 0; k < K; k++) {
        /* dirichletpdflog(A_prior[k, ], A[k, ]): upstream's argument order. */
        double lg = 0.0, sum = 0.0, l = 0.0;
        for (int32_t j = 0; j < K; j++) {
            double a = AT(model->A, k, j, K);
            sum += a;
            lg += lgamma(a);
            l += (a - 1.0) * log(AT(model->A_prior, k, j, K));
        }
        prior += lgamma(sum) - lg + l;
    }
    prior += betapdflog(point->n, ALPHA_N, BETA_N);
    prior += SP_PRIOR_LOG;
    prior += gammapdflog(point->phi, ALPHA_PHI, BETA_PHI);
    for (int32_t k = 0; k < K; k++) {
        prior += gammapdflog(point->lambda[k], ALPHA_LAMBDA, model->beta_lambda[k]);
    }
    prior += dirichletpdflog(point->pi, model->kappa, K);
    return prior;
}

/* ---- segments -------------------------------------------------------------- */

static int compare_double(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

static double segment_median(duckhts_depth_state_t *state, size_t first, size_t last) {
    size_t n = 0;
    for (size_t t = first; t <= last; t++) {
        if (isfinite(state->x2[t])) state->median_scratch[n++] = state->x2[t];
    }
    if (n == 0) return NAN;
    qsort(state->median_scratch, n, sizeof(double), compare_double);
    return n % 2 ? state->median_scratch[n / 2]
                 : 0.5 * (state->median_scratch[n / 2 - 1] + state->median_scratch[n / 2]);
}

static int push_segment(duckhts_depth_state_t *state, size_t first, size_t last, int32_t copy_number) {
    if (state->n_segments == state->segment_capacity) {
        size_t capacity = state->segment_capacity ? state->segment_capacity * 2 : 64;
        duckhts_depth_state_segment_t *segments =
            realloc(state->segments, capacity * sizeof(*segments));
        if (!segments) return -1;
        state->segments = segments;
        state->segment_capacity = capacity;
    }
    duckhts_depth_state_segment_t *segment = &state->segments[state->n_segments++];
    segment->first_bin = first;
    segment->last_bin = last;
    segment->bins = (int32_t)(last - first + 1);
    segment->state = copy_number;
    segment->copy_number = copy_number;
    segment->median_log2_ratio = segment_median(state, first, last);
    return 0;
}

static int build_segments(duckhts_depth_state_t *state) {
    state->n_segments = 0;
    for (size_t c = 0; c < state->n_contigs; c++) {
        size_t begin = state->contig_begin[c], end = state->contig_begin[c + 1];
        size_t first = begin;
        for (size_t t = begin + 1; t <= end; t++) {
            if (t == end || state->path[t] != state->path[first]) {
                if (push_segment(state, first, t - 1, state->path[first] + 1) < 0) return -1;
                first = t;
            }
        }
    }
    return 0;
}

/* logRbasedCN with cellPrev 1 and cn 2 on the log2 segment median, which
 * HMMsegment converts back from natural log before correctIntegerCN. */
static double log_ratio_copy_number(double x, double purity, double ploidy) {
    double ct = (pow(2.0, x) * (2.0 * (1.0 - purity) + purity * ploidy) - 2.0 * (1.0 - purity)) / purity;
    return fmax(ct, CN_FLOOR);
}

/* correctIntegerCN on the segments: the autosomal rule applied to every
 * segment (no chrX sex model here), correctHOMD off. */
static void correct_copy_numbers(duckhts_depth_state_t *state, int32_t max_cn,
                                 double purity, double ploidy) {
    if (!(purity >= MIN_PURITY_TO_CORRECT)) return;
    for (size_t s = 0; s < state->n_segments; s++) {
        duckhts_depth_state_segment_t *segment = &state->segments[s];
        if (isnan(segment->median_log2_ratio)) continue;
        double lr_cn = log_ratio_copy_number(segment->median_log2_ratio, purity, ploidy);
        double rounded = nearbyint(lr_cn);
        if (segment->state >= max_cn || (isfinite(lr_cn) && lr_cn >= max_cn * 1.2)) {
            segment->copy_number = (int32_t)rounded;
        }
        int32_t corrected = segment->copy_number;
        if (((rounded < ploidy && corrected > ploidy) || (rounded > ploidy && corrected < ploidy)) &&
            fabs(rounded - corrected) > 2.0) {
            segment->copy_number = (int32_t)rounded;
        }
    }
}

/* The guard of run_ichorCNA: the fraction is zero when the longest altered
 * training segment (by bases) has at most min_segment_bins bins and at most
 * alt_fraction_threshold of the training bins are altered. Altered means the
 * Viterbi state is not two copies. */
static void apply_guard(duckhts_depth_state_t *state, const duckhts_depth_state_params_t *params,
                        duckhts_depth_state_solution_t *solution) {
    size_t train_bins = 0, altered_bins = 0;
    for (size_t t = 0; t < state->n_bins; t++) {
        if (!state->train[t]) continue;
        train_bins++;
        if (state->path[t] + 1 != 2) altered_bins++;
    }
    int64_t longest = -1;
    int32_t longest_bins = -1;
    for (size_t s = 0; s < state->n_segments; s++) {
        const duckhts_depth_state_segment_t *segment = &state->segments[s];
        if (!state->train[segment->first_bin] || segment->state == 2) continue;
        int64_t length = state->end[segment->last_bin] - state->start[segment->first_bin] + 1;
        if (length > longest) {
            longest = length;
            longest_bins = segment->bins;
        }
    }
    solution->altered_fraction = train_bins ? (double)altered_bins / (double)train_bins : 0.0;
    solution->guard_zero = 0;
    if (longest_bins <= params->min_segment_bins &&
        solution->altered_fraction <= params->alt_fraction_threshold) {
        solution->normal = 1.0;
        solution->guard_zero = 1;
    }
}

/* ---- the fit ---------------------------------------------------------------- */

static void copy_point(fit_point_t *dst, const fit_point_t *src) {
    memcpy(dst, src, sizeof(*dst));
}

int duckhts_depth_state_fit(duckhts_depth_state_t *state,
                            const duckhts_depth_state_params_t *params,
                            double start_normal, double start_ploidy,
                            duckhts_depth_state_solution_t *solution) {
    if (!state || !params || !solution || state->n_bins == 0) return -2;
    int32_t K = params->max_cn;
    if (K < 2 || K > DUCKHTS_DEPTH_STATE_MAX_CN) return -2;
    if (!(params->txn_e > 0.0 && params->txn_e < 1.0) || !(params->txn_strength > 0.0) ||
        params->max_iter < 2 || !(params->lik_change > 0.0) || params->min_segment_bins < 0 ||
        !(params->alt_fraction_threshold >= 0.0 && params->alt_fraction_threshold <= 1.0)) {
        return -2;
    }
    if (!(start_normal > 0.0 && start_normal < 1.0) || !(start_ploidy >= 1.0 && start_ploidy <= PHI_HI)) {
        return -2;
    }
    if (reserve_work(state, K) < 0) return -1;

    fit_model_t model;
    memset(&model, 0, sizeof(model));
    model.K = K;

    /* getDefaultParameters on the log2 ratios of the training bins. */
    double mean = 0.0, ss = 0.0;
    size_t count = 0;
    model.first_train = SIZE_MAX;
    for (size_t t = 0; t < state->n_bins; t++) {
        if (!state->train[t]) continue;
        if (model.first_train == SIZE_MAX) model.first_train = t;
        if (!isfinite(state->x2[t])) continue;
        count++;
        double delta = state->x2[t] - mean;
        mean += delta / (double)count;
        ss += delta * (state->x2[t] - mean);
    }
    double log_ratio_var = ss / (double)(count - 1) / (double)K;
    if (!(log_ratio_var > 0.0)) return -2;
    double base_lambda = 1.0 / log_ratio_var;

    fit_point_t current;
    memset(&current, 0, sizeof(current));
    current.K = K;
    current.n = start_normal;
    current.phi = floor(start_ploidy);
    for (int32_t k = 0; k < K; k++) {
        int32_t ct = k + 1;
        double lambda = base_lambda;
        if (ct == 2) lambda = base_lambda * 10.0;
        if (ct == K) lambda = base_lambda / 15.0;
        current.lambda[k] = lambda;
        model.beta_lambda[k] = ALPHA_LAMBDA / lambda;
        model.kappa[k] = ct == 2 ? 800.0 : 200.0;
    }
    double kappa_sum = 0.0;
    for (int32_t k = 0; k < K; k++) kappa_sum += model.kappa[k];
    for (int32_t k = 0; k < K; k++) current.pi[k] = model.kappa[k] / kappa_sum;
    /* getTransitionMatrix: e on the diagonal, the rest shared equally. */
    for (int32_t i = 0; i < K; i++) {
        double row[DUCKHTS_DEPTH_STATE_MAX_CN];
        for (int32_t j = 0; j < K; j++) {
            row[j] = i == j ? params->txn_e : (1.0 - params->txn_e) / (double)(K - 1);
        }
        normalize_row(row, K);
        for (int32_t j = 0; j < K; j++) {
            AT(model.A, i, j, K) = row[j];
            AT(model.A_prior, i, j, K) = row[j];
            AT(model.dir_prior, i, j, K) = row[j] * params->txn_strength;
        }
    }

    set_means(&model, current.n, current.phi);
    compute_emissions(state, &model, current.lambda);

    fit_point_t previous;
    copy_point(&previous, &current);
    double loglik_prev = -INFINITY;
    double loglik_cur = -INFINITY;
    int32_t i = 1;
    bool converged = false;
    bool reverted = false;
    while (!converged && i < params->max_iter) {
        i++;
        /* E-step over the training bins of each contig. */
        double data_loglik = 0.0;
        memset(model.Z, 0, sizeof(model.Z));
        for (size_t c = 0; c < state->n_contigs; c++) {
            size_t T = 0;
            for (size_t t = state->contig_begin[c]; t < state->contig_begin[c + 1]; t++) {
                if (state->train[t]) state->index[T++] = t;
            }
            if (T == 0) continue;
            int status = forward_backward(state, &model, current.pi, state->index, T,
                                          &data_loglik, model.Z);
            if (status < 0) return status;
        }
        /* M-step: transitions, initial distribution, then n, phi and the
         * precisions, each from the previous value with the others fixed. */
        fit_point_t next;
        copy_point(&next, &current);
        for (int32_t k = 0; k < K; k++) {
            double row[DUCKHTS_DEPTH_STATE_MAX_CN];
            for (int32_t j = 0; j < K; j++) row[j] = AT(model.Z, k, j, K) + AT(model.dir_prior, k, j, K);
            normalize_row(row, K);
            for (int32_t j = 0; j < K; j++) AT(model.A, k, j, K) = row[j];
        }
        {
            double pi_sum = 0.0;
            for (int32_t k = 0; k < K; k++) {
                next.pi[k] = AT(state->rho, k, model.first_train, K) + model.kappa[k] - 1.0;
                pi_sum += next.pi[k];
            }
            for (int32_t k = 0; k < K; k++) next.pi[k] /= pi_sum;
        }
        objective_ctx_t ctx = {state, &model, current.n, current.phi, current.lambda, -1};
        next.n = climb(objective_n, &ctx, current.n, N_LO, N_HI, 0.02);
        ctx.n = next.n;
        next.phi = climb(objective_phi, &ctx, current.phi, PHI_LO, PHI_HI, 0.1);
        ctx.phi = next.phi;
        set_means(&model, next.n, next.phi);
        for (int32_t k = 0; k < K; k++) {
            ctx.k = k;
            ctx.lambda = current.lambda;
            next.lambda[k] = climb(objective_lambda, &ctx, current.lambda[k], LAMBDA_LO, INFINITY,
                                   0.1 * current.lambda[k]);
        }
        compute_emissions(state, &model, next.lambda);
        double loglik = data_loglik + prior_log(&model, &next);
        copy_point(&previous, &current);
        copy_point(&current, &next);
        loglik_prev = loglik_cur;
        loglik_cur = loglik;
        if (fabs(loglik_cur - loglik_prev) / fabs(fmax(loglik_prev, -1e30)) < params->lik_change) {
            converged = true;
        }
        if (loglik_cur < loglik_prev) {
            /* "Using previous iteration": n, phi, lambda, pi and the objective
             * of the previous iteration are reported; the emissions and the
             * transitions just computed stay in use, as upstream leaves them. */
            copy_point(&current, &previous);
            loglik_cur = loglik_prev;
            i--;
            converged = true;
            reverted = true;
        }
    }
    (void)reverted;

    /* runViterbi over every bin of each contig with the final pi, A and
     * emissions. */
    double log_pi[DUCKHTS_DEPTH_STATE_MAX_CN];
    double log_A[DUCKHTS_DEPTH_STATE_MAX_CN * DUCKHTS_DEPTH_STATE_MAX_CN];
    for (int32_t k = 0; k < K; k++) log_pi[k] = log(current.pi[k]);
    for (int32_t k = 0; k < K * K; k++) log_A[k] = log(model.A[k]);
    for (size_t c = 0; c < state->n_contigs; c++) {
        int status = viterbi(state, &model, log_pi, log_A, state->contig_begin[c],
                             state->contig_begin[c + 1]);
        if (status < 0) return status;
    }
    if (build_segments(state) < 0) return -1;

    solution->start_normal = start_normal;
    solution->start_ploidy = start_ploidy;
    solution->normal = current.n;
    solution->ploidy = current.phi;
    solution->log_likelihood = loglik_cur;
    solution->iterations = i - 1;
    correct_copy_numbers(state, K, 1.0 - current.n, current.phi);
    apply_guard(state, params, solution);
    return 0;
}

const duckhts_depth_state_segment_t *duckhts_depth_state_segments(
    const duckhts_depth_state_t *state, size_t *count) {
    *count = state ? state->n_segments : 0;
    return state ? state->segments : NULL;
}
