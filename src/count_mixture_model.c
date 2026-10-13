#include "include/count_mixture_model.h"
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>

_Static_assert(FLT_RADIX == 2 && DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024 &&
               sizeof(double) == sizeof(uint64_t), "Mixture sums require binary64 doubles");

static int sum_word(duckhts_mixture_sum_t *sum, size_t limb, uint64_t word) {
    while (word) {
        if (limb >= DUCKHTS_MIXTURE_SUM_LIMBS)
            return 0;
        uint64_t previous = sum->limbs[limb];
        sum->limbs[limb] += word;
        word = sum->limbs[limb] < previous;
        limb++;
    }
    return 1;
}

int duckhts_mixture_sum_add(duckhts_mixture_sum_t *sum, double value) {
    if (!isfinite(value) || value < 0 || value > DUCKHTS_MIXTURE_SUM_MAX_OPERAND)
        return 0;
    uint64_t bits;
    memcpy(&bits, &value, sizeof(bits));
    uint64_t mantissa = bits & ((UINT64_C(1) << 52) - 1);
    unsigned exponent = (unsigned)((bits >> 52) & 2047);
    unsigned offset = 0;
    if (exponent) {
        mantissa |= UINT64_C(1) << 52;
        offset = exponent - 1;
    }
    unsigned shift = offset % 64;
    size_t limb = offset / 64;
    if (!sum_word(sum, limb, mantissa << shift))
        return 0;
    return !shift || sum_word(sum, limb + 1, mantissa >> (64 - shift));
}

int duckhts_mixture_sum_merge(duckhts_mixture_sum_t *sum, const duckhts_mixture_sum_t *from) {
    for (size_t i = 0; i < DUCKHTS_MIXTURE_SUM_LIMBS; i++) {
        if (!sum_word(sum, i, from->limbs[i]))
            return 0;
    }
    return 1;
}

double duckhts_mixture_sum_value(const duckhts_mixture_sum_t *sum) {
    size_t limb = DUCKHTS_MIXTURE_SUM_LIMBS;
    while (limb && !sum->limbs[limb - 1])
        limb--;
    if (!limb)
        return 0;
    unsigned high = 0;
    uint64_t word = sum->limbs[limb - 1];
    while (word >>= 1)
        high++;
    unsigned top = (unsigned)((limb - 1) * 64) + high;
    if (top <= 52)
        return ldexp((double)sum->limbs[0], -1074);
    unsigned offset = top - 52, shift = offset % 64;
    limb = offset / 64;
    uint64_t mantissa = sum->limbs[limb] >> shift;
    if (shift && limb + 1 < DUCKHTS_MIXTURE_SUM_LIMBS)
        mantissa |= sum->limbs[limb + 1] << (64 - shift);
    unsigned guard_shift = (offset - 1) % 64;
    size_t guard_limb = (offset - 1) / 64;
    int guard = (sum->limbs[guard_limb] >> guard_shift) & 1;
    int sticky = (sum->limbs[guard_limb] & ((UINT64_C(1) << guard_shift) - 1)) != 0;
    for (size_t i = 0; i < guard_limb; i++)
        sticky |= sum->limbs[i] != 0;
    if (guard && (sticky || (mantissa & 1)))
        mantissa++;
    return ldexp((double)mantissa, (int)offset - 1074);
}

#define DIM DUCKHTS_MIXTURE_SCORE_DIM
#define GRID 128

int duckhts_mixture_merge(duckhts_mixture_block_t *to, const duckhts_mixture_block_t *from) {
    for (int i = 0; i < MIX_COUNTERS; i++) {
        if (from->counts[i] > (uint64_t)INT64_MAX - to->counts[i])
            return 0;
    }
    for (int i = 0; i < MIX_COUNTERS; i++)
        to->counts[i] += from->counts[i];
    for (int i = 0; i < 3; i++) {
        to->y[i] += from->y[i];
        for (int j = 0; j < 3; j++)
            to->s[i][j] += from->s[i][j];
    }
    for (int i = 0; i < 6; i++) {
        to->z[i] += from->z[i];
        for (int j = 0; j < 6; j++)
            to->t[i][j] += from->t[i][j];
    }
    to->denominator += from->denominator;
    to->matched = to->matched && from->matched;
    return 1;
}

int duckhts_mixture_add(duckhts_mixture_block_t *block, uint64_t ref, uint64_t alt, uint64_t other,
                        double f, double g, int usable, int excluded) {
    duckhts_mixture_block_t site;
    memset(&site, 0, sizeof(site));
    site.matched = isfinite(f) && isfinite(g) && f == g;
    site.counts[MIX_SITES] = 1;
    site.counts[MIX_OTHER] = other;
    if (excluded)
        site.counts[MIX_EXCLUDED] = 1;
    else if (!usable)
        site.counts[MIX_MISSING] = 1;
    else {
        uint64_t n = ref + alt;
        double pairs = (double)n * (double)(n ? n - 1 : 0) / 2;
        double u[3] = {1, f, g}, x[6] = {1, f, g, f * f, f * g, g * g};
        site.counts[MIX_USABLE] = 1;
        site.counts[MIX_ZERO] = n == 0;
        site.counts[MIX_SINGLETON] = n == 1;
        site.counts[MIX_READS] = n;
        site.counts[MIX_PAIR_SITES] = n > 1;
        site.counts[MIX_RR] = ref * (ref ? ref - 1 : 0) / 2;
        site.counts[MIX_RA] = ref * alt;
        site.counts[MIX_AA] = alt * (alt ? alt - 1 : 0) / 2;
        site.denominator = pairs * f * (1 - f);
        for (int i = 0; i < 3; i++) {
            site.y[i] = (double)alt * u[i];
            for (int j = 0; j < 3; j++)
                site.s[i][j] = (double)n * u[i] * u[j];
        }
        for (int i = 0; i < 6; i++) {
            site.z[i] = (double)(ref * alt) * x[i];
            for (int j = 0; j < 6; j++)
                site.t[i][j] = pairs * x[i] * x[j];
        }
    }
    return duckhts_mixture_merge(block, &site);
}

void duckhts_mixture_coefficients(double m, const duckhts_mixture_calibration_t *c, int ibd1,
                                  double beta[3], double gamma[6]) {
    double e = c->seq_error, w = c->artefact_weight, F = c->homozygosity_excess;
    double a = 1 - 2 * e, h = 1 - m;
    double vx = h * h * (1 + F) / 2, vy = m * m / 2;
    /* The conditional transmitted-allele model has C = Vx/2. */
    if (ibd1)
        vx += m * h * (1 + F) / 2;
    beta[0] = (1 - w) * e + w / 2;
    beta[1] = (1 - w) * a * h;
    beta[2] = (1 - w) * a * m;
    gamma[0] = (1 - w) * 2 * e * (1 - e) + w / 3;
    gamma[1] = (1 - w) * 2 * a * a * (h - vx);
    gamma[2] = (1 - w) * 2 * a * a * (m - vy);
    gamma[3] = (1 - w) * 2 * a * a * (vx - h * h);
    gamma[4] = (1 - w) * (-4) * a * a * h * m;
    gamma[5] = (1 - w) * 2 * a * a * (vy - m * m);
}

typedef struct {
    duckhts_mixture_block_t total;
    const duckhts_mixture_block_t *blocks;
    size_t count, omitted;
    const duckhts_mixture_calibration_t *cal;
    int ibd1, rank;
    double weight[DIM][DIM], reads, pairs;
} problem_t;

static void residual(const duckhts_mixture_block_t *b, const problem_t *p, double m,
                     double r[DIM]) {
    double beta[3], gamma[6];
    duckhts_mixture_coefficients(m, p->cal, p->ibd1, beta, gamma);
    for (int i = 0; i < 3; i++) {
        r[i] = b->y[i];
        for (int j = 0; j < 3; j++)
            r[i] -= b->s[i][j] * beta[j];
        r[i] /= p->reads > 0 ? p->reads : 1;
    }
    for (int i = 0; i < 6; i++) {
        r[i + 3] = b->z[i];
        for (int j = 0; j < 6; j++)
            r[i + 3] -= b->t[i][j] * gamma[j];
        r[i + 3] /= p->pairs > 0 ? p->pairs : 1;
    }
}

static double score(const problem_t *p, double m) {
    double r[DIM], value = 0;
    residual(&p->total, p, m, r);
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++)
            value += r[i] * p->weight[i][j] * r[j];
    }
    return fmax(0, value);
}

/* Search every grid-local minimum; a fixed lower-share tie break makes the
 * constrained fit deterministic even when several minima have equal score. */
static double minimize(const problem_t *p) {
    double best = 0, value = score(p, 0);
    for (int k = 0; k <= GRID; k++) {
        double m = 0.5 * k / GRID, v = score(p, m);
        if (v < value) {
            best = m;
            value = v;
        }
        if ((k > 0 && v > score(p, m - 0.5 / GRID)) || (k < GRID && v > score(p, m + 0.5 / GRID)))
            continue;
        double lo = fmax(0, m - 0.5 / GRID), hi = fmin(0.5, m + 0.5 / GRID);
        double x = hi - (hi - lo) * 0.6180339887498949, y = lo + (hi - lo) * 0.6180339887498949;
        double fx = score(p, x), fy = score(p, y);
        for (int iter = 0; iter < 64; iter++) {
            if (fx < fy) {
                hi = y;
                y = x;
                fy = fx;
                x = hi - (hi - lo) * 0.6180339887498949;
                fx = score(p, x);
            } else {
                lo = x;
                x = y;
                fx = fy;
                y = lo + (hi - lo) * 0.6180339887498949;
                fy = score(p, y);
            }
        }
        m = (lo + hi) / 2;
        v = score(p, m);
        if (v < value) {
            best = m;
            value = v;
        }
    }
    return best;
}

/* Symmetric Jacobi eigensystem, with relative rank cutoff after exposure
 * scaling. The pseudoinverse drops duplicate frequency-basis directions. */
static int inverse(double matrix[DIM][DIM], double inverse_matrix[DIM][DIM]) {
    double vectors[DIM][DIM] = {{0}}, largest = 0;
    int rank = 0;
    for (int i = 0; i < DIM; i++)
        vectors[i][i] = 1;
    for (int sweep = 0; sweep < 2048; sweep++) {
        int a = 0, b = 1;
        double off = 0, diagonal = 0;
        for (int i = 0; i < DIM; i++) {
            diagonal = fmax(diagonal, fabs(matrix[i][i]));
            for (int j = i + 1; j < DIM; j++) {
                if (fabs(matrix[i][j]) > off) {
                    off = fabs(matrix[i][j]);
                    a = i;
                    b = j;
                }
            }
        }
        if (off <= fmax(1e-30, diagonal * 1e-14))
            break;
        double angle = 0.5 * atan2(2 * matrix[a][b], matrix[b][b] - matrix[a][a]);
        double c = cos(angle), s = sin(angle);
        double aa = matrix[a][a], bb = matrix[b][b], ab = matrix[a][b];
        for (int i = 0; i < DIM; i++) {
            if (i != a && i != b) {
                double ia = matrix[i][a], ib = matrix[i][b];
                matrix[i][a] = matrix[a][i] = c * ia - s * ib;
                matrix[i][b] = matrix[b][i] = s * ia + c * ib;
            }
            double va = vectors[i][a], vb = vectors[i][b];
            vectors[i][a] = c * va - s * vb;
            vectors[i][b] = s * va + c * vb;
        }
        matrix[a][a] = c * c * aa - 2 * s * c * ab + s * s * bb;
        matrix[b][b] = s * s * aa + 2 * s * c * ab + c * c * bb;
        matrix[a][b] = matrix[b][a] = 0;
    }
    memset(inverse_matrix, 0, DIM * DIM * sizeof(double));
    for (int i = 0; i < DIM; i++)
        largest = fmax(largest, matrix[i][i]);
    for (int k = 0; k < DIM; k++) {
        if (matrix[k][k] <= fmax(1e-28, largest * 1e-10))
            continue;
        rank++;
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                inverse_matrix[i][j] += vectors[i][k] * vectors[j][k] / matrix[k][k];
            }
        }
    }
    return rank;
}

static double fit_problem(problem_t *p) {
    double covariance[DIM][DIM] = {{0}}, mean[DIM];
    int B = 0;
    memset(p->weight, 0, sizeof(p->weight));
    for (int i = 0; i < DIM; i++)
        p->weight[i][i] = 1;
    p->reads = p->total.s[0][0];
    p->pairs = p->total.t[0][0];
    double m = minimize(p);
    residual(&p->total, p, m, mean);
    for (size_t k = 0; k < p->count; k++) {
        if (k != p->omitted && p->blocks[k].counts[MIX_READS])
            B++;
    }
    if (B < 2) {
        p->rank = 0;
        return m;
    }
    for (size_t k = 0; k < p->count; k++) {
        double r[DIM];
        if (k == p->omitted || !p->blocks[k].counts[MIX_READS])
            continue;
        residual(&p->blocks[k], p, m, r);
        for (int i = 0; i < DIM; i++)
            r[i] -= mean[i] / B;
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++)
                covariance[i][j] += r[i] * r[j] * B / (B - 1.0);
        }
    }
    p->rank = inverse(covariance, p->weight);
    if (!p->rank) {
        for (int i = 0; i < DIM; i++)
            p->weight[i][i] = 1;
        return m;
    }
    return minimize(p);
}

static double singles(const duckhts_mixture_block_t *b, const duckhts_mixture_calibration_t *c) {
    double a = (1 - c->artefact_weight) * (1 - 2 * c->seq_error);
    double contrast = b->s[1][1] - 2 * b->s[1][2] + b->s[2][2];
    if (contrast <= fmax(1e-12, b->s[0][0] * 1e-12))
        return NAN;
    double numerator = b->y[2] - b->y[1] -
                       ((1 - c->artefact_weight) * c->seq_error + c->artefact_weight / 2) *
                           (b->s[0][2] - b->s[0][1]) -
                       a * (b->s[1][2] - b->s[1][1]);
    return numerator / (a * contrast);
}

static double pair_difference(const duckhts_mixture_block_t *b,
                              const duckhts_mixture_calibration_t *c, double m) {
    double beta[3], gamma[6], r = b->z[0];
    if (!isfinite(m) || b->t[0][0] == 0)
        return NAN;
    duckhts_mixture_coefficients(m, c, 0, beta, gamma);
    for (int i = 0; i < 6; i++)
        r -= b->t[0][i] * gamma[i];
    return r / b->t[0][0];
}

static void subtract(duckhts_mixture_block_t *to, const duckhts_mixture_block_t *b) {
    for (int i = 0; i < MIX_COUNTERS; i++)
        to->counts[i] -= b->counts[i];
    for (int i = 0; i < 3; i++) {
        to->y[i] -= b->y[i];
        for (int j = 0; j < 3; j++)
            to->s[i][j] -= b->s[i][j];
    }
    for (int i = 0; i < 6; i++) {
        to->z[i] -= b->z[i];
        for (int j = 0; j < 6; j++)
            to->t[i][j] -= b->t[i][j];
    }
}

static double root(const problem_t *p, double lo, double hi, double cutoff) {
    int accepted = score(p, lo) <= cutoff;
    for (int i = 0; i < 48; i++) {
        double mid = (lo + hi) / 2;
        if ((score(p, mid) <= cutoff) == accepted)
            lo = mid;
        else
            hi = mid;
    }
    return (lo + hi) / 2;
}

void duckhts_mixture_fit(const duckhts_mixture_block_t *blocks, size_t count,
                         const duckhts_mixture_calibration_t *cal, int ibd1,
                         duckhts_mixture_fit_t *fit) {
    problem_t p;
    memset(&p, 0, sizeof(p));
    memset(fit, 0, sizeof(*fit));
    p.blocks = blocks;
    p.count = count;
    p.omitted = SIZE_MAX;
    p.cal = cal;
    p.ibd1 = ibd1;
    p.total.matched = 1;
    fit->share = fit->lower = fit->upper = fit->se = NAN;
    fit->singles_share = fit->singles_se = fit->pair_residual = fit->pair_residual_se = NAN;
    fit->score = NAN;
    fit->status = "no_sites";
    fit->diagnostic = "uninformative";
    for (size_t k = 0; k < count; k++) {
        if (!duckhts_mixture_merge(&p.total, &blocks[k])) {
            fit->status = "model_mismatch";
            return;
        }
        fit->blocks += blocks[k].counts[MIX_READS] > 0;
        fit->pair_blocks += blocks[k].counts[MIX_PAIR_SITES] > 0;
    }
    if (!p.total.counts[MIX_READS])
        return;
    int read_blocks = fit->blocks;
    if (p.total.matched)
        fit->blocks = fit->pair_blocks;
    fit->singles_share = singles(&p.total, cal);
    fit->share = fit_problem(&p);
    fit->score = score(&p, fit->share);
    fit->rank = p.rank;
    fit->at_bound = fit->share < 1e-8 || fit->share > 0.5 - 1e-8;
    fit->pair_residual = pair_difference(&p.total, cal, fit->singles_share);
    fit->diagnostic = p.total.matched ? "unresolved" : "uninformative";
    double r0[DIM], r1[DIM], r2[DIM], sensitivity = 0;
    residual(&p.total, &p, 0, r0);
    residual(&p.total, &p, .25, r1);
    residual(&p.total, &p, .5, r2);
    for (int i = 0; i < DIM; i++) {
        sensitivity += (r0[i] - r1[i]) * (r0[i] - r1[i]) + (r1[i] - r2[i]) * (r1[i] - r2[i]);
    }
    if (sensitivity < 1e-24) {
        fit->share = NAN;
        fit->status = !p.total.counts[MIX_PAIR_SITES] && !isfinite(fit->singles_share)
                          ? "no_population_contrast"
                          : "nonidentified";
        fit->diagnostic = "uninformative";
        return;
    }
    if (!p.total.counts[MIX_PAIR_SITES])
        fit->status = "no_pairs";
    else if (fit->blocks < cal->min_blocks)
        fit->status = "insufficient_blocks";
    else if (!p.rank)
        fit->status = "nonidentified";
    else
        fit->status = "ok_conditional";
    int mismatch = 0;
    if (p.rank > 1) {
        double df = p.rank - 1.0;
        double threshold = df + 2 * sqrt(df * log(1000.0)) + 2 * log(1000.0);
        mismatch = fit->score > threshold;
        if (mismatch && fit->blocks >= cal->min_blocks) {
            fit->status = "model_mismatch";
            if (!p.total.matched)
                fit->diagnostic = "model_mismatch";
        }
    }
    /* Welford spreads for the delete-one-block refits. Each deletion repeats
     * both GMM passes with fixed external calibration. */
    double mean = 0, spread = 0, sm = 0, ss = 0, rm = 0, rs = 0;
    int B = 0, SB = 0, RB = 0;
    for (size_t k = 0; k < count; k++) {
        if (!blocks[k].counts[MIX_READS])
            continue;
        fit->leverage = fmax(fit->leverage, blocks[k].s[0][0] / p.total.s[0][0]);
        if (p.total.t[0][0] > 0)
            fit->leverage = fmax(fit->leverage, blocks[k].t[0][0] / p.total.t[0][0]);
        problem_t deletion = p;
        deletion.omitted = k;
        subtract(&deletion.total, &blocks[k]);
        if (!deletion.total.counts[MIX_READS]) {
            fit->failed_refits++;
            continue;
        }
        double m = fit_problem(&deletion), delta = m - mean;
        if (!isfinite(m) || (p.rank > 0 && deletion.rank == 0)) {
            fit->failed_refits++;
            continue;
        }
        B++;
        mean += delta / B;
        spread += delta * (m - mean);
        double s = singles(&deletion.total, cal);
        if (isfinite(s)) {
            delta = s - sm;
            SB++;
            sm += delta / SB;
            ss += delta * (s - sm);
        }
        double r = pair_difference(&deletion.total, cal, s);
        if (isfinite(r)) {
            delta = r - rm;
            RB++;
            rm += delta / RB;
            rs += delta * (r - rm);
        }
    }
    if (B > 1)
        fit->se = sqrt((B - 1.0) / B * spread);
    if (SB > 1)
        fit->singles_se = sqrt((SB - 1.0) / SB * ss);
    if (RB > 1)
        fit->pair_residual_se = sqrt((RB - 1.0) / RB * rs);
    if (fit->blocks < cal->min_blocks || !p.rank || B != read_blocks || mismatch)
        return;
    if (!p.total.matched && isfinite(fit->pair_residual_se) && fit->pair_residual_se > 0) {
        double z = fit->pair_residual / fit->pair_residual_se;
        fit->diagnostic = fabs(z) <= 1.959963984540054 ? "compatible_unrelated"
                          : z < 0                      ? "shared_allele_signal"
                                                       : "model_mismatch";
    }
    if (fit->share < 1e-8)
        fit->diagnostic = "uninformative";
    /* Interior limits use the deletion variance to calibrate score curvature.
     * At a bound the block sandwich score cutoff avoids constrained-Wald
     * collapse; its reference cutoff is chi-square(1), 95 percent. */
    double h = 1e-4, cutoff = 3.841458820694124;
    if (!fit->at_bound && fit->se > 0) {
        double curvature =
            (score(&p, fit->share + h) + score(&p, fit->share - h) - 2 * fit->score) / (2 * h * h);
        if (curvature > 0)
            cutoff *= curvature * fit->se * fit->se;
    }
    cutoff += fit->score;
    fit->lower = NAN;
    fit->upper = NAN;
    int previous = score(&p, 0) <= cutoff, runs = previous;
    if (previous)
        fit->lower = 0;
    for (int k = 1; k <= 1024; k++) {
        double lo = 0.5 * (k - 1) / 1024, hi = 0.5 * k / 1024;
        int accepted = score(&p, hi) <= cutoff;
        if (accepted && !previous) {
            double edge = root(&p, lo, hi, cutoff);
            if (!isfinite(fit->lower))
                fit->lower = edge;
            runs++;
        }
        if (!accepted && previous)
            fit->upper = root(&p, lo, hi, cutoff);
        if (accepted && k == 1024)
            fit->upper = 0.5;
        previous = accepted;
    }
    /* A narrow accepted interval may lie wholly between grid points. */
    if (!isfinite(fit->lower)) {
        fit->lower = root(&p, 0, fit->share, cutoff);
        fit->upper = root(&p, fit->share, 0.5, cutoff);
    }
    fit->interval_hull = runs > 1;
    if (fit->lower == 0 && fit->upper == 0.5)
        fit->status = "uninformative";
}
