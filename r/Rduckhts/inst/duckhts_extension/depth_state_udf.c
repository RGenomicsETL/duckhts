/* __duckhts_depth_state_decode: the copy-number state model of
 * depth_state_hmm.h as a scalar over one sample's bin lists, behind the
 * duckhts_depth_state_solutions and duckhts_depth_state_segments macros of
 * depth_sql.c.
 *
 * __duckhts_depth_state_decode(chroms, starts, ends, log2_ratios, train,
 *     normal_starts, ploidy_starts, max_cn, txn_e, txn_strength, max_iter,
 *     min_segment_bins, alt_fraction_threshold, close_margin, max_bins)
 * returns one struct per solution and segment. Bins are listed in the order
 * the macro gives them: by contig, then start; a contig is the run of equal
 * chrom values. Every start of the grid (ploidy_starts x normal_starts, with
 * upstream's skip of normal 0.95 at ploidies other than 2) gives one solution,
 * ranked by its objective; the status is guard_zero when the guard set the
 * fraction to zero, close_to_best when a lower-ranked solution's objective is
 * within close_margin of the best, ok otherwise.
 *
 * Memory: the engine holds about max_cn * 28 bytes a bin plus 40 bytes a bin;
 * the solutions hold one segment record (48 bytes) per segment and start.
 * max_bins bounds the bins of one sample. */
#if defined(__MINGW32__) && !defined(__USE_MINGW_ANSI_STDIO)
#define __USE_MINGW_ANSI_STDIO 1
#endif
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "duckdb_list.h"
#include "depth_state_hmm.h"

#define DEPTH_STATE_NAME "duckhts_depth_state"
#define DEPTH_STATE_MAX_STARTS 100
#define DEPTH_STATE_ERRLEN 256

enum {
    OUT_SOLUTION,
    OUT_START_NORMAL,
    OUT_START_PLOIDY,
    OUT_NORMAL,
    OUT_PLOIDY,
    OUT_LOGLIK,
    OUT_ALTERED,
    OUT_ITERATIONS,
    OUT_STATUS,
    OUT_CHROM,
    OUT_START,
    OUT_END,
    OUT_BINS,
    OUT_STATE,
    OUT_COPY_NUMBER,
    OUT_MEDIAN,
    OUT_COUNT
};

typedef struct {
    bool valid;
    const duckdb_list_entry *entries;
    duckdb_vector child;
    const uint64_t *child_validity;
} list_t;

typedef struct {
    duckhts_depth_state_solution_t solution;
    duckhts_depth_state_segment_t *segments;
    size_t n_segments;
    size_t order;            /* position in the grid, for a stable rank */
    int32_t rank;
} fitted_t;

typedef struct {
    int32_t *contig_id;
    int64_t *start;
    int64_t *end;
    double *ratio;
    uint8_t *train;
    char **contig_name;      /* one per contig id */
    size_t n_contigs;
    size_t contig_capacity;
    size_t n_bins;
    fitted_t fits[DEPTH_STATE_MAX_STARTS];
    size_t n_fits;
} sample_t;

static bool row_valid(duckdb_vector vector, idx_t row) {
    uint64_t *mask = duckdb_vector_get_validity(vector);
    return !mask || duckdb_validity_row_is_valid(mask, row);
}

static bool element_valid(const uint64_t *mask, idx_t index) {
    return !mask || duckdb_validity_row_is_valid((uint64_t *)mask, index);
}

static list_t open_list(duckdb_vector vector, idx_t row) {
    list_t list = {0};
    list.valid = row_valid(vector, row);
    list.entries = duckdb_vector_get_data(vector);
    list.child = duckdb_list_vector_get_child(vector);
    list.child_validity = duckdb_vector_get_validity(list.child);
    return list;
}

static void set_error(duckdb_function_info info, const char *message) {
    char full[DEPTH_STATE_ERRLEN + sizeof(DEPTH_STATE_NAME) + 4];
    snprintf(full, sizeof(full), "%s: %s", DEPTH_STATE_NAME, message);
    duckdb_scalar_function_set_error(info, full);
}

static void sample_free(sample_t *sample) {
    free(sample->contig_id);
    free(sample->start);
    free(sample->end);
    free(sample->ratio);
    free(sample->train);
    for (size_t c = 0; c < sample->n_contigs; c++) free(sample->contig_name[c]);
    free(sample->contig_name);
    for (size_t f = 0; f < sample->n_fits; f++) free(sample->fits[f].segments);
    memset(sample, 0, sizeof(*sample));
}

static bool read_scalar_double(duckdb_vector vector, idx_t row, double *value) {
    if (!row_valid(vector, row)) return false;
    *value = ((const double *)duckdb_vector_get_data(vector))[row];
    return true;
}

static bool read_scalar_int32(duckdb_vector vector, idx_t row, int32_t *value) {
    if (!row_valid(vector, row)) return false;
    *value = ((const int32_t *)duckdb_vector_get_data(vector))[row];
    return true;
}

static bool read_scalar_int64(duckdb_vector vector, idx_t row, int64_t *value) {
    if (!row_valid(vector, row)) return false;
    *value = ((const int64_t *)duckdb_vector_get_data(vector))[row];
    return true;
}

/* Copies the bin lists of one row into the sample. Returns 1, 0 for an
 * error already reported, -1 for a NULL list (a NULL output row). */
static int load_bins(duckdb_function_info info, sample_t *sample, duckdb_vector *vectors,
                     idx_t row, size_t max_bins) {
    list_t chroms = open_list(vectors[0], row);
    list_t starts = open_list(vectors[1], row);
    list_t ends = open_list(vectors[2], row);
    list_t ratios = open_list(vectors[3], row);
    list_t trains = open_list(vectors[4], row);
    if (!chroms.valid || !starts.valid || !ends.valid || !ratios.valid || !trains.valid) return -1;
    duckdb_list_entry entry = chroms.entries[row];
    if (starts.entries[row].length != entry.length || ends.entries[row].length != entry.length ||
        ratios.entries[row].length != entry.length || trains.entries[row].length != entry.length) {
        set_error(info, "the bin lists of one sample differ in length");
        return 0;
    }
    if (entry.length == 0) {
        set_error(info, "a sample has no kept bins");
        return 0;
    }
    if (entry.length > max_bins) {
        set_error(info, "a sample has more kept bins than max_bins");
        return 0;
    }
    size_t n = entry.length;
    sample->contig_id = malloc(n * sizeof(*sample->contig_id));
    sample->start = malloc(n * sizeof(*sample->start));
    sample->end = malloc(n * sizeof(*sample->end));
    sample->ratio = malloc(n * sizeof(*sample->ratio));
    sample->train = malloc(n * sizeof(*sample->train));
    if (!sample->contig_id || !sample->start || !sample->end || !sample->ratio || !sample->train) {
        set_error(info, "out of memory");
        return 0;
    }
    const duckdb_string_t *chrom_data = duckdb_vector_get_data(chroms.child);
    const int64_t *start_data = duckdb_vector_get_data(starts.child);
    const int64_t *end_data = duckdb_vector_get_data(ends.child);
    const double *ratio_data = duckdb_vector_get_data(ratios.child);
    const bool *train_data = duckdb_vector_get_data(trains.child);
    for (size_t i = 0; i < n; i++) {
        idx_t at = entry.offset + i;
        if (!element_valid(chroms.child_validity, at) || !element_valid(starts.child_validity, at) ||
            !element_valid(ends.child_validity, at) || !element_valid(trains.child_validity, at)) {
            set_error(info, "chrom, start, end and the training flag of a bin must not be NULL");
            return 0;
        }
        const char *name = duckdb_string_t_data((duckdb_string_t *)&chrom_data[at]);
        uint32_t length = duckdb_string_t_length(chrom_data[at]);
        bool new_contig = sample->n_contigs == 0 ||
                          strlen(sample->contig_name[sample->n_contigs - 1]) != length ||
                          memcmp(sample->contig_name[sample->n_contigs - 1], name, length) != 0;
        if (new_contig) {
            if (sample->n_contigs == sample->contig_capacity) {
                size_t capacity = sample->contig_capacity ? sample->contig_capacity * 2 : 32;
                char **names = realloc(sample->contig_name, capacity * sizeof(*names));
                if (!names) {
                    set_error(info, "out of memory");
                    return 0;
                }
                sample->contig_name = names;
                sample->contig_capacity = capacity;
            }
            char *copy = malloc((size_t)length + 1);
            if (!copy) {
                set_error(info, "out of memory");
                return 0;
            }
            memcpy(copy, name, length);
            copy[length] = '\0';
            if (memchr(copy, '\0', length) != NULL) {
                free(copy);
                set_error(info, "a chrom value contains a NUL byte");
                return 0;
            }
            sample->contig_name[sample->n_contigs++] = copy;
        }
        sample->contig_id[i] = (int32_t)(sample->n_contigs - 1);
        sample->start[i] = start_data[at];
        sample->end[i] = end_data[at];
        if (sample->end[i] < sample->start[i]) {
            set_error(info, "a bin ends before it starts");
            return 0;
        }
        sample->ratio[i] = element_valid(ratios.child_validity, at) ? ratio_data[at] : NAN;
        sample->train[i] = train_data[at] ? 1u : 0u;
    }
    sample->n_bins = n;
    return 1;
}

static bool read_start_list(duckdb_function_info info, duckdb_vector vector, idx_t row,
                            const char *what, double *values, size_t *count) {
    list_t list = open_list(vector, row);
    if (!list.valid) {
        char message[DEPTH_STATE_ERRLEN];
        snprintf(message, sizeof(message), "%s must not be NULL", what);
        set_error(info, message);
        return false;
    }
    duckdb_list_entry entry = list.entries[row];
    if (entry.length == 0 || entry.length > DEPTH_STATE_MAX_STARTS) {
        char message[DEPTH_STATE_ERRLEN];
        snprintf(message, sizeof(message), "%s needs 1 to %d values", what, DEPTH_STATE_MAX_STARTS);
        set_error(info, message);
        return false;
    }
    const double *data = duckdb_vector_get_data(list.child);
    for (size_t i = 0; i < entry.length; i++) {
        if (!element_valid(list.child_validity, entry.offset + i)) {
            char message[DEPTH_STATE_ERRLEN];
            snprintf(message, sizeof(message), "%s must not contain NULL", what);
            set_error(info, message);
            return false;
        }
        values[i] = data[entry.offset + i];
    }
    *count = entry.length;
    return true;
}

static int compare_fits(const void *a, const void *b) {
    const fitted_t *x = a, *y = b;
    if (x->solution.log_likelihood != y->solution.log_likelihood) {
        return x->solution.log_likelihood > y->solution.log_likelihood ? -1 : 1;
    }
    return x->order < y->order ? -1 : x->order > y->order ? 1 : 0;
}

/* Fits every start of the grid. Returns 1, or 0 after reporting an error. */
static int fit_sample(duckdb_function_info info, sample_t *sample, duckhts_depth_state_t *state,
                      const duckhts_depth_state_params_t *params, const double *normals,
                      size_t n_normals, const double *ploidies, size_t n_ploidies, size_t max_bins,
                      double close_margin) {
    int status = duckhts_depth_state_set_bins(state, sample->n_bins, sample->contig_id, sample->start,
                                              sample->end, sample->ratio, sample->train, max_bins);
    if (status == -1) {
        set_error(info, "out of memory");
        return 0;
    }
    if (status == -2) {
        set_error(info, "a sample needs at least two training bins with a finite log2 ratio");
        return 0;
    }
    if (n_normals * n_ploidies > DEPTH_STATE_MAX_STARTS) {
        set_error(info, "the start grid has more than 100 points");
        return 0;
    }
    for (size_t p = 0; p < n_ploidies; p++) {
        for (size_t q = 0; q < n_normals; q++) {
            double n0 = normals[q], p0 = ploidies[p];
            if (!(n0 > 0.0 && n0 < 1.0)) {
                set_error(info, "normal_fractions must be in (0, 1)");
                return 0;
            }
            if (!(p0 >= 1.0 && p0 <= 10.0)) {
                set_error(info, "ploidies must be in [1, 10]");
                return 0;
            }
            /* run_ichorCNA skips normal 0.95 unless the ploidy start is 2. */
            if (n0 == 0.95 && p0 != 2.0) continue;
            fitted_t *fit = &sample->fits[sample->n_fits];
            memset(fit, 0, sizeof(*fit));
            status = duckhts_depth_state_fit(state, params, n0, p0, &fit->solution);
            if (status == -1) {
                set_error(info, "out of memory");
                return 0;
            }
            if (status == -2) {
                set_error(info, "invalid model parameters");
                return 0;
            }
            if (status == -3) {
                set_error(info, "a bin has zero likelihood in every state; its log2 ratio is too far from every state mean");
                return 0;
            }
            size_t count;
            const duckhts_depth_state_segment_t *segments = duckhts_depth_state_segments(state, &count);
            fit->segments = malloc(count * sizeof(*segments));
            if (!fit->segments) {
                set_error(info, "out of memory");
                return 0;
            }
            memcpy(fit->segments, segments, count * sizeof(*segments));
            fit->n_segments = count;
            fit->order = sample->n_fits;
            sample->n_fits++;
        }
    }
    if (sample->n_fits == 0) {
        set_error(info, "the start grid is empty after upstream's skip of normal 0.95 at ploidies other than 2");
        return 0;
    }
    qsort(sample->fits, sample->n_fits, sizeof(sample->fits[0]), compare_fits);
    for (size_t f = 0; f < sample->n_fits; f++) sample->fits[f].rank = (int32_t)f + 1;
    (void)close_margin;
    return 1;
}

static const char *status_name(const sample_t *sample, const fitted_t *fit, double close_margin) {
    if (fit->solution.guard_zero) return "guard_zero";
    if (fit->rank > 1 &&
        sample->fits[0].solution.log_likelihood - fit->solution.log_likelihood <= close_margin) {
        return "close_to_best";
    }
    return "ok";
}

static bool write_rows(duckdb_function_info info, duckdb_vector output, idx_t row,
                       const sample_t *sample, double close_margin) {
    size_t total = 0;
    for (size_t f = 0; f < sample->n_fits; f++) total += sample->fits[f].n_segments;
    duckdb_list_entry entry;
    if (!duckhts_list_extend(output, total, &entry)) {
        set_error(info, "output allocation failed");
        return false;
    }
    duckdb_vector child = duckdb_list_vector_get_child(output);
    duckdb_vector columns[OUT_COUNT];
    for (int i = 0; i < OUT_COUNT; i++) columns[i] = duckdb_struct_vector_get_child(child, i);
    int32_t *solution = duckdb_vector_get_data(columns[OUT_SOLUTION]);
    double *start_normal = duckdb_vector_get_data(columns[OUT_START_NORMAL]);
    double *start_ploidy = duckdb_vector_get_data(columns[OUT_START_PLOIDY]);
    double *normal = duckdb_vector_get_data(columns[OUT_NORMAL]);
    double *ploidy = duckdb_vector_get_data(columns[OUT_PLOIDY]);
    double *loglik = duckdb_vector_get_data(columns[OUT_LOGLIK]);
    double *altered = duckdb_vector_get_data(columns[OUT_ALTERED]);
    int32_t *iterations = duckdb_vector_get_data(columns[OUT_ITERATIONS]);
    int64_t *start = duckdb_vector_get_data(columns[OUT_START]);
    int64_t *end = duckdb_vector_get_data(columns[OUT_END]);
    int32_t *bins = duckdb_vector_get_data(columns[OUT_BINS]);
    int32_t *state = duckdb_vector_get_data(columns[OUT_STATE]);
    int32_t *copy_number = duckdb_vector_get_data(columns[OUT_COPY_NUMBER]);
    double *median = duckdb_vector_get_data(columns[OUT_MEDIAN]);
    duckdb_vector_ensure_validity_writable(columns[OUT_MEDIAN]);
    uint64_t *median_validity = duckdb_vector_get_validity(columns[OUT_MEDIAN]);
    idx_t at = entry.offset;
    for (size_t f = 0; f < sample->n_fits; f++) {
        const fitted_t *fit = &sample->fits[f];
        const char *status = status_name(sample, fit, close_margin);
        for (size_t s = 0; s < fit->n_segments; s++, at++) {
            const duckhts_depth_state_segment_t *segment = &fit->segments[s];
            solution[at] = fit->rank;
            start_normal[at] = fit->solution.start_normal;
            start_ploidy[at] = fit->solution.start_ploidy;
            normal[at] = fit->solution.normal;
            ploidy[at] = fit->solution.ploidy;
            loglik[at] = fit->solution.log_likelihood;
            altered[at] = fit->solution.altered_fraction;
            iterations[at] = fit->solution.iterations;
            duckdb_vector_assign_string_element(columns[OUT_STATUS], at, status);
            duckdb_vector_assign_string_element(columns[OUT_CHROM], at,
                                                sample->contig_name[sample->contig_id[segment->first_bin]]);
            start[at] = sample->start[segment->first_bin];
            end[at] = sample->end[segment->last_bin];
            bins[at] = segment->bins;
            state[at] = segment->state;
            copy_number[at] = segment->copy_number;
            if (isnan(segment->median_log2_ratio)) {
                median[at] = 0.0;
                duckdb_validity_set_row_invalid(median_validity, at);
            } else {
                median[at] = segment->median_log2_ratio;
                duckdb_validity_set_row_valid(median_validity, at);
            }
        }
    }
    ((duckdb_list_entry *)duckdb_vector_get_data(output))[row] = entry;
    return true;
}

static void set_row_null(duckdb_vector output, idx_t row) {
    duckdb_vector_ensure_validity_writable(output);
    duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
}

static void depth_state_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector vectors[15];
    for (idx_t i = 0; i < 15; i++) vectors[i] = duckdb_data_chunk_get_vector(input, i);
    duckhts_depth_state_t *state = duckhts_depth_state_create();
    if (!state) {
        set_error(info, "out of memory");
        return;
    }
    duckdb_list_vector_set_size(output, 0);
    idx_t rows = duckdb_data_chunk_get_size(input);
    for (idx_t row = 0; row < rows; row++) {
        duckhts_depth_state_params_t params;
        double normals[DEPTH_STATE_MAX_STARTS], ploidies[DEPTH_STATE_MAX_STARTS];
        size_t n_normals, n_ploidies;
        double close_margin;
        int64_t max_bins;
        int32_t max_cn, max_iter, min_segment_bins;
        if (!read_scalar_int32(vectors[7], row, &max_cn) || !read_scalar_double(vectors[8], row, &params.txn_e) ||
            !read_scalar_double(vectors[9], row, &params.txn_strength) ||
            !read_scalar_int32(vectors[10], row, &max_iter) ||
            !read_scalar_int32(vectors[11], row, &min_segment_bins) ||
            !read_scalar_double(vectors[12], row, &params.alt_fraction_threshold) ||
            !read_scalar_double(vectors[13], row, &close_margin) ||
            !read_scalar_int64(vectors[14], row, &max_bins)) {
            set_error(info, "options must not be NULL");
            break;
        }
        if (max_cn < 2 || max_cn > DUCKHTS_DEPTH_STATE_MAX_CN) {
            set_error(info, "max_cn must be in [2, 20]");
            break;
        }
        if (!(params.txn_e > 0.0 && params.txn_e < 1.0)) {
            set_error(info, "txn_e must be in (0, 1)");
            break;
        }
        if (!(params.txn_strength > 0.0 && isfinite(params.txn_strength))) {
            set_error(info, "txn_strength must be positive");
            break;
        }
        if (max_iter < 2 || max_iter > 1000) {
            set_error(info, "max_iter must be in [2, 1000]");
            break;
        }
        if (min_segment_bins < 0) {
            set_error(info, "min_segment_bins must not be negative");
            break;
        }
        if (!(params.alt_fraction_threshold >= 0.0 && params.alt_fraction_threshold <= 1.0)) {
            set_error(info, "alt_fraction_threshold must be in [0, 1]");
            break;
        }
        if (!(close_margin >= 0.0)) {
            set_error(info, "close_margin must not be negative");
            break;
        }
        if (max_bins < 1) {
            set_error(info, "max_bins must be at least 1");
            break;
        }
        params.max_cn = max_cn;
        params.max_iter = max_iter;
        params.lik_change = 1e-4;
        params.min_segment_bins = min_segment_bins;
        if (!read_start_list(info, vectors[5], row, "normal_fractions", normals, &n_normals) ||
            !read_start_list(info, vectors[6], row, "ploidies", ploidies, &n_ploidies)) {
            break;
        }
        sample_t sample;
        memset(&sample, 0, sizeof(sample));
        int loaded = load_bins(info, &sample, vectors, row, (size_t)max_bins);
        if (loaded < 0) {
            sample_free(&sample);
            set_row_null(output, row);
            continue;
        }
        bool ok = loaded == 1 &&
                  fit_sample(info, &sample, state, &params, normals, n_normals, ploidies, n_ploidies,
                             (size_t)max_bins, close_margin) &&
                  write_rows(info, output, row, &sample, close_margin);
        sample_free(&sample);
        if (!ok) break;
    }
    duckhts_depth_state_destroy(state);
}

bool register_duckhts_depth_state_functions(duckdb_connection connection) {
    duckdb_logical_type varchar = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type varchar_list = duckdb_create_list_type(varchar);
    duckdb_logical_type bigint_list = duckdb_create_list_type(bigint);
    duckdb_logical_type real_list = duckdb_create_list_type(real);
    duckdb_logical_type boolean_list = duckdb_create_list_type(boolean);
    duckdb_logical_type member_types[OUT_COUNT] = {
        integer, real, real, real, real, real, real, integer, varchar,
        varchar, bigint, bigint, integer, integer, integer, real};
    const char *member_names[OUT_COUNT] = {
        "solution", "start_normal", "start_ploidy", "normal_fraction", "ploidy",
        "log_likelihood", "altered_fraction", "iterations", "status",
        "chrom", "start", "end", "bins", "state", "copy_number", "median_log2_ratio"};
    duckdb_logical_type row = duckdb_create_struct_type(member_types, member_names, OUT_COUNT);
    duckdb_logical_type rows = duckdb_create_list_type(row);

    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, "__duckhts_depth_state_decode");
    duckdb_scalar_function_add_parameter(function, varchar_list);
    duckdb_scalar_function_add_parameter(function, bigint_list);
    duckdb_scalar_function_add_parameter(function, bigint_list);
    duckdb_scalar_function_add_parameter(function, real_list);
    duckdb_scalar_function_add_parameter(function, boolean_list);
    duckdb_scalar_function_add_parameter(function, real_list);
    duckdb_scalar_function_add_parameter(function, real_list);
    duckdb_scalar_function_add_parameter(function, integer);
    duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_add_parameter(function, integer);
    duckdb_scalar_function_add_parameter(function, integer);
    duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_add_parameter(function, bigint);
    duckdb_scalar_function_set_return_type(function, rows);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, depth_state_scalar);
    bool ok = duckdb_register_scalar_function(connection, function) == DuckDBSuccess;
    duckdb_destroy_scalar_function(&function);

    duckdb_destroy_logical_type(&rows);
    duckdb_destroy_logical_type(&row);
    duckdb_destroy_logical_type(&boolean_list);
    duckdb_destroy_logical_type(&real_list);
    duckdb_destroy_logical_type(&bigint_list);
    duckdb_destroy_logical_type(&varchar_list);
    duckdb_destroy_logical_type(&boolean);
    duckdb_destroy_logical_type(&real);
    duckdb_destroy_logical_type(&integer);
    duckdb_destroy_logical_type(&bigint);
    duckdb_destroy_logical_type(&varchar);
    return ok;
}
