/* The count-error model: likelihood on histogram cells and its maximum by
 * Nelder-Mead from a fixed set of starts. No DuckDB dependency. */
#include "include/count_error_model.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define GENOTYPES 3
#define PAIRS (GENOTYPES * GENOTYPES)
#define BINOMIAL_SPREAD 1e-9
#define LOG_FLOOR 1e-300

/* Natural bounds of the parameters. Each is an open interval reached through
 * a logistic transform, so the optimizer works without constraints. */
static const double parameter_scale[DUCKHTS_COUNT_ERROR_PARAMETERS] = {0.45, 0.45, 1.0, 1.0, 0.5, 0.5, 0.5};

/* to_natural and to_working read the parameter struct as an array. */
_Static_assert(sizeof(duckhts_count_error_params_t) == DUCKHTS_COUNT_ERROR_PARAMETERS * sizeof(double),
               "the parameters are seven doubles in order");

static double logistic(double x) { return 1.0 / (1.0 + exp(-x)); }
static double logit(double p) { return log(p / (1.0 - p)); }

static void to_natural(const double *working, duckhts_count_error_params_t *params) {
    double *values = (double *)params;
    for (int i = 0; i < DUCKHTS_COUNT_ERROR_PARAMETERS; i++) {
        values[i] = parameter_scale[i] * logistic(working[i]);
    }
}

static void to_working(const duckhts_count_error_params_t *params, double *working) {
    const double *values = (const double *)params;
    for (int i = 0; i < DUCKHTS_COUNT_ERROR_PARAMETERS; i++) {
        double share = values[i] / parameter_scale[i];
        if (share < 1e-9) share = 1e-9;
        if (share > 1 - 1e-9) share = 1 - 1e-9;
        working[i] = logit(share);
    }
}

/* The distinct (depth, alt) pairs of the cells, with each cell's pair index:
 * the count likelihood of a pair is computed once per parameter vector. */
typedef struct {
    uint16_t depth;
    uint16_t alt;
    double log_choose;
} count_pair_t;

typedef struct {
    const duckhts_count_cell_t *cells;
    size_t count;
    count_pair_t *pairs;
    size_t pair_count;
    uint32_t *pair_of;          /* cell index -> pair index */
    double *pair_probability;   /* P(alt | depth, g, h): pair index * PAIRS + g * GENOTYPES + h */
    uint16_t max_depth;         /* the deepest pair */
    double *gamma_table;        /* 3 * (max_depth + 1): lgamma(x + a), lgamma(x + b), lgamma(x + a + b) */
    duckhts_count_error_relation_t relation;
    int fixed[DUCKHTS_COUNT_ERROR_PARAMETERS]; /* 1 where a parameter is held */
    double held[DUCKHTS_COUNT_ERROR_PARAMETERS];
    long evaluations;
} count_problem_t;

static int compare_pairs(const void *a, const void *b) {
    const count_pair_t *x = a, *y = b;
    if (x->depth != y->depth) return x->depth < y->depth ? -1 : 1;
    if (x->alt != y->alt) return x->alt < y->alt ? -1 : 1;
    return 0;
}

static void problem_release(count_problem_t *problem) {
    free(problem->pairs);
    free(problem->pair_of);
    free(problem->pair_probability);
    free(problem->gamma_table);
    memset(problem, 0, sizeof(*problem));
}

static int problem_init(count_problem_t *problem, const duckhts_count_cell_t *cells, size_t count,
                        duckhts_count_error_relation_t relation) {
    count_pair_t *sorted;
    memset(problem, 0, sizeof(*problem));
    problem->cells = cells;
    problem->count = count;
    problem->relation = relation;
    if (count == 0) return 1;
    sorted = malloc(count * sizeof(*sorted));
    problem->pair_of = malloc(count * sizeof(*problem->pair_of));
    if (!sorted || !problem->pair_of) {
        free(sorted);
        problem_release(problem);
        return 0;
    }
    for (size_t i = 0; i < count; i++) {
        sorted[i].depth = cells[i].depth;
        sorted[i].alt = cells[i].alt;
        sorted[i].log_choose = 0;
    }
    qsort(sorted, count, sizeof(*sorted), compare_pairs);
    problem->pairs = malloc(count * sizeof(*problem->pairs));
    if (!problem->pairs) {
        free(sorted);
        problem_release(problem);
        return 0;
    }
    for (size_t i = 0; i < count; i++) {
        if (problem->pair_count == 0 || compare_pairs(&sorted[i], &problem->pairs[problem->pair_count - 1]) != 0) {
            count_pair_t *pair = &problem->pairs[problem->pair_count++];
            pair->depth = sorted[i].depth;
            pair->alt = sorted[i].alt;
            pair->log_choose = lgamma(pair->depth + 1.0) - lgamma(pair->alt + 1.0) -
                               lgamma(pair->depth - pair->alt + 1.0);
        }
    }
    free(sorted);
    problem->max_depth = problem->pairs[problem->pair_count - 1].depth;
    problem->pair_probability = malloc(problem->pair_count * PAIRS * sizeof(double));
    problem->gamma_table = malloc(3 * (problem->max_depth + 1u) * sizeof(double));
    if (!problem->pair_probability || !problem->gamma_table) {
        problem_release(problem);
        return 0;
    }
    for (size_t i = 0; i < count; i++) {
        count_pair_t key = {cells[i].depth, cells[i].alt, 0};
        const count_pair_t *found = bsearch(&key, problem->pairs, problem->pair_count, sizeof(key), compare_pairs);
        problem->pair_of[i] = (uint32_t)(found - problem->pairs);
    }
    return 1;
}

/* Fills P(alt | depth, g, h) of every pair: a beta-binomial, or the binomial
 * when the spread is (near) zero. The beta-binomial needs lgamma at alt + a,
 * depth - alt + b and depth + a + b, where a and b depend on (g, h) only, so
 * each (g, h) fills three tables over 0..max_depth once and the pairs read
 * them. */
static void fill_pair_probabilities(count_problem_t *problem, const duckhts_count_error_params_t *params) {
    const double q[GENOTYPES] = {params->seq_error, params->allele_balance, 1 - params->seq_error};
    const double rho[GENOTYPES] = {params->spread_hom, params->spread_het, params->spread_hom};
    const size_t width = problem->max_depth + 1u;
    double *table_a = problem->gamma_table;
    double *table_b = table_a + width;
    double *table_ab = table_b + width;
    for (int g = 0; g < GENOTYPES; g++) {
        for (int h = 0; h < GENOTYPES; h++) {
            double p = (1 - params->contamination) * q[g] + params->contamination * q[h];
            const int at = g * GENOTYPES + h;
            if (p < 1e-12) p = 1e-12;
            if (p > 1 - 1e-12) p = 1 - 1e-12;
            if (rho[g] < BINOMIAL_SPREAD) {
                const double log_p = log(p), log_not_p = log1p(-p);
                for (size_t i = 0; i < problem->pair_count; i++) {
                    const count_pair_t *pair = &problem->pairs[i];
                    problem->pair_probability[i * PAIRS + at] =
                        exp(pair->log_choose + pair->alt * log_p + (pair->depth - pair->alt) * log_not_p);
                }
            } else {
                const double scale = (1 - rho[g]) / rho[g];
                const double a = p * scale, b = (1 - p) * scale;
                const double constant = lgamma(a + b) - lgamma(a) - lgamma(b);
                for (size_t x = 0; x < width; x++) {
                    table_a[x] = lgamma(x + a);
                    table_b[x] = lgamma(x + b);
                    table_ab[x] = lgamma(x + a + b);
                }
                for (size_t i = 0; i < problem->pair_count; i++) {
                    const count_pair_t *pair = &problem->pairs[i];
                    problem->pair_probability[i * PAIRS + at] =
                        exp(pair->log_choose + constant + table_a[pair->alt] + table_b[pair->depth - pair->alt] -
                            table_ab[pair->depth]);
                }
            }
        }
    }
}

static double problem_log_likelihood(count_problem_t *problem, const duckhts_count_error_params_t *params) {
    double total = 0;
    const double f_share = params->homozygosity_excess;
    problem->evaluations++;
    fill_pair_probabilities(problem, params);
    for (size_t i = 0; i < problem->count; i++) {
        const duckhts_count_cell_t *cell = &problem->cells[i];
        const double *pair = &problem->pair_probability[problem->pair_of[i] * PAIRS];
        const double f = cell->af;
        const double hwe[GENOTYPES] = {(1 - f) * (1 - f), 2 * f * (1 - f), f * f};
        const double prior[GENOTYPES] = {(1 - f_share) * hwe[0] + f_share * (1 - f), (1 - f_share) * hwe[1],
                                         (1 - f_share) * hwe[2] + f_share * f};
        double mixture = 0;
        if (cell->depth == 0) continue;
        for (int g = 0; g < GENOTYPES; g++) {
            double given_g[GENOTYPES];
            if (problem->relation == DUCKHTS_COUNT_ERROR_UNRELATED) {
                given_g[0] = hwe[0]; given_g[1] = hwe[1]; given_g[2] = hwe[2];
            } else if (g == 0) {
                /* one allele of g (reference), the other a draw at f */
                given_g[0] = 1 - f; given_g[1] = f; given_g[2] = 0;
            } else if (g == 1) {
                given_g[0] = (1 - f) / 2; given_g[1] = 0.5; given_g[2] = f / 2;
            } else {
                given_g[0] = 0; given_g[1] = 1 - f; given_g[2] = f;
            }
            for (int h = 0; h < GENOTYPES; h++) {
                if (given_g[h] == 0) continue;
                mixture += prior[g] * given_g[h] * pair[g * GENOTYPES + h];
            }
        }
        mixture = (1 - params->artefact_weight) * mixture + params->artefact_weight / (cell->depth + 1.0);
        total += cell->sites * log(mixture > LOG_FLOOR ? mixture : LOG_FLOOR);
    }
    return total;
}

double duckhts_count_error_log_likelihood(const duckhts_count_cell_t *cells, size_t count,
                                          const duckhts_count_error_params_t *params,
                                          duckhts_count_error_relation_t relation) {
    count_problem_t problem;
    double value;
    if (!problem_init(&problem, cells, count, relation)) return -HUGE_VAL;
    value = problem_log_likelihood(&problem, params);
    problem_release(&problem);
    return value;
}

/* The objective of the optimizer: minus the log-likelihood at a working
 * point, with the held parameters put back. */
static double objective(count_problem_t *problem, const double *working) {
    duckhts_count_error_params_t params;
    double full[DUCKHTS_COUNT_ERROR_PARAMETERS] = {0};
    double value;
    for (int i = 0, free_at = 0; i < DUCKHTS_COUNT_ERROR_PARAMETERS; i++) {
        full[i] = problem->fixed[i] ? problem->held[i] : working[free_at++];
    }
    to_natural(full, &params);
    value = -problem_log_likelihood(problem, &params);
    return isfinite(value) ? value : HUGE_VAL;
}

/* Nelder-Mead over `n` free working parameters. Returns 1 when the simplex
 * shrank below the tolerance within the step limit, 0 at the step limit. */
static int nelder_mead(count_problem_t *problem, double *point, int n, double step, int max_steps) {
    enum { MAX_FREE = DUCKHTS_COUNT_ERROR_PARAMETERS };
    double simplex[MAX_FREE + 1][MAX_FREE] = {{0}};
    double values[MAX_FREE + 1] = {0};
    double centroid[MAX_FREE] = {0}, reflected[MAX_FREE] = {0}, expanded[MAX_FREE] = {0}, contracted[MAX_FREE] = {0};
    const double alpha = 1.0, gamma = 2.0, rho = 0.5, sigma = 0.5;
    int converged = 0;

    for (int v = 0; v <= n; v++) {
        for (int i = 0; i < n; i++) simplex[v][i] = point[i] + (v == i + 1 ? step : 0);
        values[v] = objective(problem, simplex[v]);
    }
    for (int iteration = 0; iteration < max_steps; iteration++) {
        /* Order the vertices: best first. */
        for (int a = 1; a <= n; a++) {
            double value = values[a];
            double row[MAX_FREE];
            int b = a - 1;
            memcpy(row, simplex[a], sizeof(row));
            while (b >= 0 && values[b] > value) {
                values[b + 1] = values[b];
                memcpy(simplex[b + 1], simplex[b], sizeof(row));
                b--;
            }
            values[b + 1] = value;
            memcpy(simplex[b + 1], row, sizeof(row));
        }
        {
            /* The width of the simplex is measured on the natural scale: a
             * parameter at its bound saturates the logistic transform, so the
             * working coordinates there can stay wide while nothing moves. */
            double spread = 0;
            for (int i = 0; i < n; i++) {
                double range = fabs(logistic(simplex[n][i]) - logistic(simplex[0][i]));
                if (range > spread) spread = range;
            }
            /* Done when the simplex is small, or when the objective no longer
             * differs across it beyond what the summation resolves: a
             * log-likelihood flat to 1e-4 over the simplex is at its maximum
             * for any use of the fit. */
            if (isfinite(values[n]) && (spread < 1e-6 || fabs(values[n] - values[0]) < 1e-4)) {
                converged = 1;
                break;
            }
        }
        for (int i = 0; i < n; i++) {
            centroid[i] = 0;
            for (int v = 0; v < n; v++) centroid[i] += simplex[v][i];
            centroid[i] /= n;
            reflected[i] = centroid[i] + alpha * (centroid[i] - simplex[n][i]);
        }
        {
            double reflected_value = objective(problem, reflected);
            if (reflected_value < values[0]) {
                double expanded_value;
                for (int i = 0; i < n; i++) expanded[i] = centroid[i] + gamma * (reflected[i] - centroid[i]);
                expanded_value = objective(problem, expanded);
                if (expanded_value < reflected_value) {
                    memcpy(simplex[n], expanded, n * sizeof(double));
                    values[n] = expanded_value;
                } else {
                    memcpy(simplex[n], reflected, n * sizeof(double));
                    values[n] = reflected_value;
                }
                continue;
            }
            if (reflected_value < values[n - 1]) {
                memcpy(simplex[n], reflected, n * sizeof(double));
                values[n] = reflected_value;
                continue;
            }
            {
                double contracted_value;
                const double *toward = reflected_value < values[n] ? reflected : simplex[n];
                for (int i = 0; i < n; i++) contracted[i] = centroid[i] + rho * (toward[i] - centroid[i]);
                contracted_value = objective(problem, contracted);
                if (contracted_value < (reflected_value < values[n] ? reflected_value : values[n])) {
                    memcpy(simplex[n], contracted, n * sizeof(double));
                    values[n] = contracted_value;
                    continue;
                }
            }
            for (int v = 1; v <= n; v++) {
                for (int i = 0; i < n; i++) simplex[v][i] = simplex[0][i] + sigma * (simplex[v][i] - simplex[0][i]);
                values[v] = objective(problem, simplex[v]);
            }
        }
    }
    memcpy(point, simplex[0], n * sizeof(double));
    return converged;
}

/* Runs the optimizer from `start` on the free parameters and writes the fit. */
static void fit_from(count_problem_t *problem, const duckhts_count_error_params_t *start, int max_steps,
                     duckhts_count_error_fit_t *fit) {
    double working[DUCKHTS_COUNT_ERROR_PARAMETERS] = {0};
    double full[DUCKHTS_COUNT_ERROR_PARAMETERS] = {0};
    double free_point[DUCKHTS_COUNT_ERROR_PARAMETERS] = {0};
    int free_count = 0;
    to_working(start, full);
    for (int i = 0; i < DUCKHTS_COUNT_ERROR_PARAMETERS; i++) {
        if (problem->fixed[i]) problem->held[i] = full[i]; else free_point[free_count++] = full[i];
    }
    /* A second run from the first result, with a smaller simplex, confirms a
     * minimum that the first run reached at its step limit. */
    fit->converged = nelder_mead(problem, free_point, free_count, 0.5, max_steps) ||
                     nelder_mead(problem, free_point, free_count, 0.1, max_steps);
    for (int i = 0, free_at = 0; i < DUCKHTS_COUNT_ERROR_PARAMETERS; i++) {
        working[i] = problem->fixed[i] ? problem->held[i] : free_point[free_at++];
    }
    to_natural(working, &fit->params);
    fit->log_likelihood = problem_log_likelihood(problem, &fit->params);
    fit->at_bound = fit->params.contamination > 0.449 || fit->params.artefact_weight > 0.499 ||
                    fit->params.seq_error > 0.449;
}

int duckhts_count_error_fit(const duckhts_count_cell_t *cells, size_t count,
                            duckhts_count_error_relation_t relation, duckhts_count_error_fit_t *fit) {
    /* The starts differ in contamination, where the likelihood can have more
     * than one optimum; the other starts are usual values. */
    static const double contamination_starts[] = {0.005, 0.1};
    count_problem_t problem;
    duckhts_count_error_fit_t best;
    int have_best = 0;
    if (!problem_init(&problem, cells, count, relation)) return 0;
    for (size_t s = 0; s < sizeof(contamination_starts) / sizeof(contamination_starts[0]); s++) {
        duckhts_count_error_params_t start = {1e-3, contamination_starts[s], 0.05, 0.48, 1e-3, 1e-3, 1e-3};
        duckhts_count_error_fit_t candidate;
        fit_from(&problem, &start, 2000, &candidate);
        if (!have_best || candidate.log_likelihood > best.log_likelihood) {
            best = candidate;
            have_best = 1;
        }
    }
    problem_release(&problem);
    *fit = best;
    return 1;
}

int duckhts_count_error_fit_block(const duckhts_count_cell_t *cells, size_t count,
                                  const duckhts_count_error_params_t *fixed,
                                  duckhts_count_error_fit_t *fit) {
    count_problem_t problem;
    if (!problem_init(&problem, cells, count, DUCKHTS_COUNT_ERROR_UNRELATED)) return 0;
    for (int i = 0; i < DUCKHTS_COUNT_ERROR_PARAMETERS; i++) problem.fixed[i] = 1;
    problem.fixed[1] = 0; /* contamination */
    fit_from(&problem, fixed, 400, fit);
    problem_release(&problem);
    return 1;
}
