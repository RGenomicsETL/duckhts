#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN
#include "include/count_mixture_model.h"
#include "include/count_cells.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ERRLEN 256
#define MAX_BLOCKS 4096
#define MAX_CANDIDATES 64
#define MAGIC UINT64_C(0x4448544d49583101)

typedef struct {
    uint64_t magic, count;
} mixture_header_t;
typedef struct {
    unsigned char *buffer;
    size_t count, capacity;
    uint64_t limit;
    int64_t max_blocks;
} mixture_state_t;

static void aggregate_error(duckdb_function_info info, const char *message) {
    duckdb_aggregate_function_set_error(info, message);
}
static bool valid(duckdb_vector v, idx_t row) {
    uint64_t *mask = duckdb_vector_get_validity(v);
    return !mask || duckdb_validity_row_is_valid(mask, row);
}
static duckhts_mixture_block_t *blocks(mixture_state_t *s) {
    return (duckhts_mixture_block_t *)(s->buffer + sizeof(mixture_header_t));
}
static size_t bytes(size_t capacity) {
    return sizeof(mixture_header_t) + capacity * sizeof(duckhts_mixture_block_t);
}
static bool reserve(mixture_state_t *s, size_t count, char *error) {
    if (count > (size_t)s->max_blocks) {
        snprintf(error, ERRLEN, "max_blocks exceeded; increase block_bases or max_blocks");
        return false;
    }
    if (count <= s->capacity)
        return true;
    size_t capacity = s->capacity ? s->capacity * 2 : 1;
    if (capacity < count)
        capacity = count;
    if (capacity > (size_t)s->max_blocks)
        capacity = (size_t)s->max_blocks;
    size_t old_bytes = s->buffer ? bytes(s->capacity) : 0;
    size_t new_bytes = bytes(capacity);
    uint64_t held;
    /* A moving realloc may temporarily own both allocations. */
    if (!duckhts_count_bytes_charge(new_bytes, s->limit, &held)) {
        snprintf(error, ERRLEN, "max_state_bytes exceeded by low-pass moment state");
        return false;
    }
    unsigned char *buffer = realloc(s->buffer, new_bytes);
    if (!buffer) {
        duckhts_count_bytes_release(new_bytes);
        snprintf(error, ERRLEN, "out of memory for low-pass moment blocks");
        return false;
    }
    duckhts_count_bytes_release(old_bytes);
    s->buffer = buffer;
    s->capacity = capacity;
    return true;
}
static void release(mixture_state_t *s) {
    if (s->buffer) {
        free(s->buffer);
        duckhts_count_bytes_release(bytes(s->capacity));
    }
    memset(s, 0, sizeof(*s));
}
static idx_t state_size(duckdb_function_info info) {
    (void)info;
    return sizeof(mixture_state_t);
}
static void state_init(duckdb_function_info info, duckdb_aggregate_state state) {
    (void)info;
    memset(state, 0, sizeof(mixture_state_t));
}
static void state_destroy(duckdb_aggregate_state *states, idx_t count) {
    for (idx_t i = 0; i < count; i++)
        release((mixture_state_t *)states[i]);
}
static duckhts_mixture_block_t *get_block(mixture_state_t *s, int64_t id, char *error) {
    for (size_t k = 0; k < s->count; k++)
        if (blocks(s)[k].id == id)
            return &blocks(s)[k];
    if (!reserve(s, s->count + 1, error))
        return NULL;
    duckhts_mixture_block_t *b = &blocks(s)[s->count++];
    memset(b, 0, sizeof(*b));
    b->id = id;
    b->matched = 1;
    return b;
}
enum { IN_COUNTS, IN_S, IN_Y, IN_T, IN_Z, IN_DENOMINATOR, IN_MATCHED, IN_FIELDS };
static const char *input_names[IN_FIELDS] = {"counts", "s",           "y",      "t",
                                             "z",      "denominator", "matched"};
static const size_t input_lengths[5] = {MIX_COUNTERS, 9, 3, 36, 6};
static duckdb_logical_type moments_type(void) {
    duckdb_logical_type fields[IN_FIELDS];
    for (int i = 0; i < IN_FIELDS; i++) {
        duckdb_logical_type child =
            duckdb_create_logical_type(i == IN_COUNTS    ? DUCKDB_TYPE_BIGINT
                                       : i == IN_MATCHED ? DUCKDB_TYPE_BOOLEAN
                                                         : DUCKDB_TYPE_DOUBLE);
        fields[i] = i < IN_DENOMINATOR ? duckdb_create_list_type(child) : child;
        if (i < IN_DENOMINATOR)
            duckdb_destroy_logical_type(&child);
    }
    duckdb_logical_type type = duckdb_create_struct_type(fields, input_names, IN_FIELDS);
    for (int i = 0; i < IN_FIELDS; i++)
        duckdb_destroy_logical_type(&fields[i]);
    return type;
}
static void update(duckdb_function_info info, duckdb_data_chunk input,
                   duckdb_aggregate_state *states) {
    duckdb_vector v[4], fields[IN_FIELDS];
    for (int i = 0; i < 4; i++)
        v[i] = duckdb_data_chunk_get_vector(input, i);
    for (int i = 0; i < IN_FIELDS; i++)
        fields[i] = duckdb_struct_vector_get_child(v[1], i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        char error[ERRLEN];
        mixture_state_t *s = (mixture_state_t *)states[row];
        if (!s || !valid(v[0], row) || !valid(v[1], row) || !valid(v[2], row) ||
            !valid(v[3], row)) {
            aggregate_error(info, "invalid low-pass moment row");
            return;
        }
        int64_t id = ((int64_t *)duckdb_vector_get_data(v[0]))[row];
        int64_t max_blocks = ((int64_t *)duckdb_vector_get_data(v[2]))[row];
        int64_t limit = ((int64_t *)duckdb_vector_get_data(v[3]))[row];
        if (id < 0 || max_blocks < 1 || max_blocks > MAX_BLOCKS || limit < 1) {
            aggregate_error(info, "invalid low-pass state limits or block");
            return;
        }
        if (s->limit && (s->limit != (uint64_t)limit || s->max_blocks != max_blocks)) {
            aggregate_error(info, "low-pass state limits must be constant");
            return;
        }
        s->limit = (uint64_t)limit;
        s->max_blocks = max_blocks;
        duckhts_mixture_block_t b;
        memset(&b, 0, sizeof(b));
        b.id = id;
        double *destinations[4] = {&b.s[0][0], b.y, &b.t[0][0], b.z};
        for (int i = 0; i < IN_FIELDS; i++) {
            if (!valid(fields[i], row)) {
                aggregate_error(info, "NULL low-pass moments");
                return;
            }
            if (i >= IN_DENOMINATOR)
                continue;
            duckdb_list_entry entry = ((duckdb_list_entry *)duckdb_vector_get_data(fields[i]))[row];
            if (entry.length != input_lengths[i] ||
                entry.offset > duckdb_list_vector_get_size(fields[i]) ||
                entry.length > duckdb_list_vector_get_size(fields[i]) - entry.offset) {
                aggregate_error(info, "invalid low-pass moment list length");
                return;
            }
            duckdb_vector child = duckdb_list_vector_get_child(fields[i]);
            for (size_t k = 0; k < entry.length; k++) {
                if (!valid(child, entry.offset + k)) {
                    aggregate_error(info, "NULL low-pass moment element");
                    return;
                }
                if (i == IN_COUNTS) {
                    int64_t value = ((int64_t *)duckdb_vector_get_data(child))[entry.offset + k];
                    if (value < 0) {
                        aggregate_error(info, "negative low-pass counter");
                        return;
                    }
                    b.counts[k] = (uint64_t)value;
                } else {
                    double value = ((double *)duckdb_vector_get_data(child))[entry.offset + k];
                    if (!isfinite(value)) {
                        aggregate_error(info, "nonfinite low-pass moment");
                        return;
                    }
                    destinations[i - 1][k] = value;
                }
            }
        }
        b.denominator = ((double *)duckdb_vector_get_data(fields[IN_DENOMINATOR]))[row];
        b.matched = ((bool *)duckdb_vector_get_data(fields[IN_MATCHED]))[row];
        duckhts_mixture_block_t *to = get_block(s, id, error);
        if (!to) {
            aggregate_error(info, error);
            return;
        }
        if (!duckhts_mixture_merge(to, &b)) {
            aggregate_error(info, "low-pass integer counters overflow BIGINT");
            return;
        }
    }
}
static void combine(duckdb_function_info info, duckdb_aggregate_state *source,
                    duckdb_aggregate_state *target, idx_t count) {
    for (idx_t i = 0; i < count; i++) {
        char error[ERRLEN];
        mixture_state_t *from = (mixture_state_t *)source[i], *to = (mixture_state_t *)target[i];
        if (!from->count)
            continue;
        if (to->limit && (to->limit != from->limit || to->max_blocks != from->max_blocks)) {
            aggregate_error(info, "low-pass state limits must be constant");
            return;
        }
        to->limit = from->limit;
        to->max_blocks = from->max_blocks;
        for (size_t k = 0; k < from->count; k++) {
            duckhts_mixture_block_t *b = get_block(to, blocks(from)[k].id, error);
            if (!b) {
                aggregate_error(info, error);
                return;
            }
            if (!duckhts_mixture_merge(b, &blocks(from)[k])) {
                aggregate_error(info, "low-pass integer counters overflow BIGINT");
                return;
            }
        }
    }
}
static void finalize(duckdb_function_info info, duckdb_aggregate_state *states,
                     duckdb_vector output, idx_t count, idx_t offset) {
    (void)info;
    duckdb_vector_ensure_validity_writable(output);
    for (idx_t i = 0; i < count; i++) {
        mixture_state_t *s = (mixture_state_t *)states[i];
        if (!s->count) {
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), i + offset);
            continue;
        }
        mixture_header_t header = {MAGIC, s->count};
        memcpy(s->buffer, &header, sizeof(header));
        duckdb_vector_assign_string_element_len(output, i + offset, (const char *)s->buffer,
                                                bytes(s->count));
        release(s);
    }
}

/* One named definition controls both the result schema and assignment order. */
enum {
    OUT_SHARE,
    OUT_LOWER,
    OUT_UPPER,
    OUT_SE,
    OUT_SINGLES,
    OUT_SINGLES_SE,
    OUT_RESIDUAL,
    OUT_RESIDUAL_SE,
    OUT_DENOMINATOR,
    OUT_RATIO,
    OUT_ERROR,
    OUT_ARTEFACT,
    OUT_F,
    OUT_SCORE,
    OUT_LEVERAGE,
    OUT_SITES,
    OUT_USABLE,
    OUT_ZERO,
    OUT_SINGLETON,
    OUT_READS,
    OUT_OTHER,
    OUT_PAIR_SITES,
    OUT_EXCLUDED,
    OUT_RR,
    OUT_RA,
    OUT_AA,
    OUT_MISSING,
    OUT_BLOCKS,
    OUT_PAIR_BLOCKS,
    OUT_RANK,
    OUT_FAILED,
    OUT_AT_BOUND,
    OUT_MATCHED,
    OUT_HULL,
    OUT_STATUS,
    OUT_DIAGNOSTIC,
    OUT_METHOD,
    OUT_COUNT
};
static const char *names[OUT_COUNT] = {"share",
                                       "share_lower",
                                       "share_upper",
                                       "share_se",
                                       "singles_share",
                                       "singles_se",
                                       "pair_residual",
                                       "pair_residual_se",
                                       "pair_denominator",
                                       "discordance_ratio",
                                       "seq_error",
                                       "artefact_weight",
                                       "homozygosity_excess",
                                       "moment_score",
                                       "max_block_leverage",
                                       "sites",
                                       "usable_sites",
                                       "zero_depth_sites",
                                       "singleton_sites",
                                       "reads",
                                       "other_reads",
                                       "positive_pair_sites",
                                       "excluded_sites",
                                       "concordant_ref_pairs",
                                       "discordant_pairs",
                                       "concordant_alt_pairs",
                                       "unavailable_sites",
                                       "blocks",
                                       "pair_blocks",
                                       "score_rank",
                                       "jackknife_failed_refits",
                                       "at_bound",
                                       "matched_frequencies",
                                       "interval_hull",
                                       "status",
                                       "relation_diagnostic",
                                       "method"};
static duckdb_type kind(int i) {
    if (i < OUT_SITES)
        return DUCKDB_TYPE_DOUBLE;
    if (i <= OUT_FAILED)
        return DUCKDB_TYPE_BIGINT;
    if (i <= OUT_HULL)
        return DUCKDB_TYPE_BOOLEAN;
    return DUCKDB_TYPE_VARCHAR;
}
static duckdb_logical_type result_type(void) {
    duckdb_logical_type fields[OUT_COUNT];
    for (int i = 0; i < OUT_COUNT; i++)
        fields[i] = duckdb_create_logical_type(kind(i));
    duckdb_logical_type result = duckdb_create_struct_type(fields, names, OUT_COUNT);
    for (int i = 0; i < OUT_COUNT; i++)
        duckdb_destroy_logical_type(&fields[i]);
    return result;
}
static const char *calibration_error(double e, double w, double F) {
    if (!isfinite(e) || e < 0 || e >= 0.5)
        return "seq_error must be finite and in [0,0.5)";
    if (!isfinite(w) || w < 0 || w >= 1)
        return "artefact_weight must be finite and in [0,1)";
    if (!isfinite(F) || F < 0 || F > 1)
        return "homozygosity_excess must be finite and in [0,1]";
    return NULL;
}
static void valid_args(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector v[9];
    for (int i = 0; i < 9; i++)
        v[i] = duckdb_data_chunk_get_vector(input, i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        for (int i = 0; i < 9; i++) {
            if (!valid(v[i], row)) {
                duckdb_scalar_function_set_error(info, "low-pass arguments must not be NULL");
                return;
            }
        }
        const char *error = calibration_error(((double *)duckdb_vector_get_data(v[0]))[row],
                                              ((double *)duckdb_vector_get_data(v[1]))[row],
                                              ((double *)duckdb_vector_get_data(v[2]))[row]);
        if (error) {
            duckdb_scalar_function_set_error(info, error);
            return;
        }
        for (int i = 3; i < 9; i++) {
            int64_t n = ((int64_t *)duckdb_vector_get_data(v[i]))[row];
            if (n < 1 || (i == 4 && n > 65535) || (i == 5 && (n < 2 || n > INT32_MAX)) ||
                (i == 6 && n > MAX_BLOCKS) || (i == 7 && n > MAX_CANDIDATES)) {
                duckdb_scalar_function_set_error(
                    info, "invalid low-pass resource limit (min_blocks >= 2, max_depth <= 65535, "
                          "max_blocks <= 4096, max_candidates <= 64)");
                return;
            }
        }
        ((bool *)duckdb_vector_get_data(output))[row] = true;
    }
}
static int compare_blocks(const void *a, const void *b) {
    const duckhts_mixture_block_t *x = a, *y = b;
    return x->id < y->id ? -1 : x->id > y->id;
}
static void fit_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector v[7], out[OUT_COUNT];
    for (int i = 0; i < 7; i++)
        v[i] = duckdb_data_chunk_get_vector(input, i);
    for (int i = 0; i < OUT_COUNT; i++)
        out[i] = duckdb_struct_vector_get_child(output, i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        for (int i = 1; i < 7; i++) {
            if (!valid(v[i], row)) {
                duckdb_scalar_function_set_error(info, "low-pass arguments must not be NULL");
                return;
            }
        }
        duckhts_mixture_calibration_t c = {((double *)duckdb_vector_get_data(v[1]))[row],
                                           ((double *)duckdb_vector_get_data(v[2]))[row],
                                           ((double *)duckdb_vector_get_data(v[3]))[row], 0};
        int64_t min_blocks = ((int64_t *)duckdb_vector_get_data(v[4]))[row];
        int64_t limit = ((int64_t *)duckdb_vector_get_data(v[6]))[row];
        const char *error =
            calibration_error(c.seq_error, c.artefact_weight, c.homozygosity_excess);
        if (!error && (min_blocks < 2 || min_blocks > INT32_MAX || limit < 1))
            error = "invalid low-pass fit limits";
        if (error) {
            duckdb_scalar_function_set_error(info, error);
            return;
        }
        c.min_blocks = (int)min_blocks;
        bool ibd1 = ((bool *)duckdb_vector_get_data(v[5]))[row];
        size_t count = 0, allocated = 0;
        duckhts_mixture_block_t *b = NULL, total;
        memset(&total, 0, sizeof(total));
        total.matched = 1;
        if (valid(v[0], row)) {
            duckdb_string_t blob = ((duckdb_string_t *)duckdb_vector_get_data(v[0]))[row];
            uint32_t length = duckdb_string_t_length(blob);
            const char *data = duckdb_string_t_data(&blob);
            mixture_header_t header;
            if (length < sizeof(header)) {
                duckdb_scalar_function_set_error(info, "invalid low-pass state header");
                return;
            }
            memcpy(&header, data, sizeof(header));
            if (header.magic != MAGIC || header.count > MAX_BLOCKS || header.count == 0 ||
                length != bytes(header.count)) {
                duckdb_scalar_function_set_error(info, "invalid low-pass state layout");
                return;
            }
            count = (size_t)header.count;
            allocated = count * sizeof(*b);
            uint64_t held;
            if (!duckhts_count_bytes_charge(allocated, (uint64_t)limit, &held)) {
                duckdb_scalar_function_set_error(info, "max_state_bytes exceeded by low-pass fit");
                return;
            }
            b = malloc(allocated);
            if (!b) {
                duckhts_count_bytes_release(allocated);
                duckdb_scalar_function_set_error(info, "out of memory for low-pass fit");
                return;
            }
            memcpy(b, data + sizeof(header), allocated);
            qsort(b, count, sizeof(*b), compare_blocks);
            for (size_t k = 0; k < count; k++) {
                if (!duckhts_mixture_merge(&total, &b[k])) {
                    free(b);
                    duckhts_count_bytes_release(allocated);
                    duckdb_scalar_function_set_error(info,
                                                     "low-pass integer counters overflow BIGINT");
                    return;
                }
            }
        }
        duckhts_mixture_fit_t f;
        duckhts_mixture_fit(b, count, &c, ibd1, &f);
        free(b);
        if (allocated)
            duckhts_count_bytes_release(allocated);
        double values[OUT_SITES] = {
            f.share,
            f.lower,
            f.upper,
            f.se,
            f.singles_share,
            f.singles_se,
            f.pair_residual,
            f.pair_residual_se,
            total.denominator,
            total.denominator > 0 ? (double)total.counts[MIX_RA] / total.denominator : NAN,
            c.seq_error,
            c.artefact_weight,
            c.homozygosity_excess,
            f.score,
            f.leverage};
        for (int i = 0; i < OUT_SITES; i++) {
            if (!isfinite(values[i])) {
                duckdb_vector_ensure_validity_writable(out[i]);
                duckdb_validity_set_row_invalid(duckdb_vector_get_validity(out[i]), row);
            } else
                ((double *)duckdb_vector_get_data(out[i]))[row] = values[i];
        }
        for (int i = 0; i < MIX_COUNTERS; i++)
            ((int64_t *)duckdb_vector_get_data(out[OUT_SITES + i]))[row] = (int64_t)total.counts[i];
        ((int64_t *)duckdb_vector_get_data(out[OUT_BLOCKS]))[row] = f.blocks;
        ((int64_t *)duckdb_vector_get_data(out[OUT_PAIR_BLOCKS]))[row] = f.pair_blocks;
        ((int64_t *)duckdb_vector_get_data(out[OUT_RANK]))[row] = f.rank;
        ((int64_t *)duckdb_vector_get_data(out[OUT_FAILED]))[row] = f.failed_refits;
        ((bool *)duckdb_vector_get_data(out[OUT_AT_BOUND]))[row] = f.at_bound;
        ((bool *)duckdb_vector_get_data(out[OUT_MATCHED]))[row] = total.matched;
        ((bool *)duckdb_vector_get_data(out[OUT_HULL]))[row] = f.interval_hull;
        duckdb_vector_assign_string_element(out[OUT_STATUS], row, f.status);
        duckdb_vector_assign_string_element(out[OUT_DIAGNOSTIC], row, f.diagnostic);
        duckdb_vector_assign_string_element(out[OUT_METHOD], row,
                                            "duckhts_count_mixture_lowpass v1");
    }
}

static idx_t exact_sum_size(duckdb_function_info info) {
    (void)info;
    return sizeof(duckhts_mixture_sum_t);
}

static void exact_sum_init(duckdb_function_info info, duckdb_aggregate_state state) {
    (void)info;
    memset(state, 0, sizeof(duckhts_mixture_sum_t));
}

static void exact_sum_update(duckdb_function_info info, duckdb_data_chunk input,
                             duckdb_aggregate_state *states) {
    duckdb_vector vector = duckdb_data_chunk_get_vector(input, 0);
    double *values = duckdb_vector_get_data(vector);
    idx_t count = duckdb_data_chunk_get_size(input);
    for (idx_t i = 0; i < count; i++) {
        if (valid(vector, i) &&
            !duckhts_mixture_sum_add((duckhts_mixture_sum_t *)states[i], values[i])) {
            aggregate_error(info, "Invalid or overflowing low-pass moment sum");
            return;
        }
    }
}

static void exact_sum_combine(duckdb_function_info info, duckdb_aggregate_state *source,
                              duckdb_aggregate_state *target, idx_t count) {
    for (idx_t i = 0; i < count; i++) {
        if (!duckhts_mixture_sum_merge((duckhts_mixture_sum_t *)target[i],
                                       (duckhts_mixture_sum_t *)source[i])) {
            aggregate_error(info, "Low-pass moment sum capacity exceeded");
            return;
        }
    }
}

static void exact_sum_finalize(duckdb_function_info info, duckdb_aggregate_state *states,
                               duckdb_vector result, idx_t count, idx_t offset) {
    (void)info;
    double *values = duckdb_vector_get_data(result);
    for (idx_t i = 0; i < count; i++)
        values[offset + i] = duckhts_mixture_sum_value((duckhts_mixture_sum_t *)states[i]);
}

bool register_duckhts_count_mixture_functions(duckdb_connection connection) {
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type blob = duckdb_create_logical_type(DUCKDB_TYPE_BLOB);
    duckdb_logical_type result = result_type();
    duckdb_aggregate_function aggregate = duckdb_create_aggregate_function();
    duckdb_aggregate_function_set_name(aggregate, "__duckhts_mixture_blocks");
    duckdb_logical_type moments = moments_type();
    duckdb_aggregate_function_add_parameter(aggregate, bigint);
    duckdb_aggregate_function_add_parameter(aggregate, moments);
    duckdb_aggregate_function_add_parameter(aggregate, bigint);
    duckdb_aggregate_function_add_parameter(aggregate, bigint);
    duckdb_destroy_logical_type(&moments);
    duckdb_aggregate_function_set_return_type(aggregate, blob);
    duckdb_aggregate_function_set_special_handling(aggregate);
    duckdb_aggregate_function_set_functions(aggregate, state_size, state_init, update, combine,
                                            finalize);
    duckdb_aggregate_function_set_destructor(aggregate, state_destroy);
    bool ok = duckdb_register_aggregate_function(connection, aggregate) == DuckDBSuccess;
    duckdb_destroy_aggregate_function(&aggregate);
    aggregate = duckdb_create_aggregate_function();
    duckdb_aggregate_function_set_name(aggregate, "__duckhts_mixture_sum");
    duckdb_aggregate_function_add_parameter(aggregate, real);
    duckdb_aggregate_function_set_return_type(aggregate, real);
    duckdb_aggregate_function_set_functions(aggregate, exact_sum_size, exact_sum_init,
                                            exact_sum_update, exact_sum_combine, exact_sum_finalize);
    ok = ok && duckdb_register_aggregate_function(connection, aggregate) == DuckDBSuccess;
    duckdb_destroy_aggregate_function(&aggregate);
    duckdb_scalar_function scalar = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(scalar, "__duckhts_mixture_fit");
    for (int i = 0; i < 7; i++)
        duckdb_scalar_function_add_parameter(scalar, i == 0   ? blob
                                                     : i < 4  ? real
                                                     : i == 5 ? boolean
                                                              : bigint);
    duckdb_scalar_function_set_return_type(scalar, result);
    duckdb_scalar_function_set_special_handling(scalar);
    duckdb_scalar_function_set_function(scalar, fit_scalar);
    ok = ok && duckdb_register_scalar_function(connection, scalar) == DuckDBSuccess;
    duckdb_destroy_scalar_function(&scalar);
    scalar = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(scalar, "__duckhts_mixture_valid_args");
    for (int i = 0; i < 9; i++)
        duckdb_scalar_function_add_parameter(scalar, i < 3 ? real : bigint);
    duckdb_scalar_function_set_return_type(scalar, boolean);
    duckdb_scalar_function_set_special_handling(scalar);
    duckdb_scalar_function_set_function(scalar, valid_args);
    ok = ok && duckdb_register_scalar_function(connection, scalar) == DuckDBSuccess;
    duckdb_destroy_scalar_function(&scalar);
    duckdb_destroy_logical_type(&result);
    duckdb_destroy_logical_type(&blob);
    duckdb_destroy_logical_type(&boolean);
    duckdb_destroy_logical_type(&real);
    duckdb_destroy_logical_type(&bigint);
    return ok;
}
