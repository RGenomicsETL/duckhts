/* duckhts_roh_segments: vectorised runs-of-homozygosity kernel over one list per row.
 *
 * The caller builds per-sample, per-chromosome lists with list(... ORDER BY pos). One
 * row is decoded at a time, so peak memory follows the longest list, not the chunk.
 */
#if defined(__MINGW32__) && !defined(__USE_MINGW_ANSI_STDIO)
#define __USE_MINGW_ANSI_STDIO 1
#endif
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "duckdb_list.h"
#include "roh_hmm.h"

#define ROH_NAME "duckhts_roh_segments"
#define ROH_ERRLEN 256
#define ROH_MAX_POS1 INT32_MAX

typedef struct {
    bool has_gt; /* dosage evidence with gt_error instead of PL */
    duckdb_vector positions, af, evidence, gt_error, map_pos, map_cm, rec_rate, hw_to_az, az_to_hw;
} roh_args_t;

typedef struct {
    int32_t *map_pos0;
    double *map_rate;
    size_t map_capacity;
} roh_map_t;

typedef struct {
    bool valid;
    const duckdb_list_entry *entries;
    duckdb_vector child;
    const uint64_t *child_validity;
} roh_list_t;

static bool row_valid(duckdb_vector vector, idx_t row) {
    uint64_t *mask = duckdb_vector_get_validity(vector);
    return !mask || duckdb_validity_row_is_valid(mask, row);
}

static bool element_valid(const uint64_t *mask, idx_t index) {
    return !mask || duckdb_validity_row_is_valid((uint64_t *)mask, index);
}

static roh_list_t open_list(duckdb_vector vector, idx_t row) {
    roh_list_t list = {0};
    list.valid = row_valid(vector, row);
    list.entries = duckdb_vector_get_data(vector);
    list.child = duckdb_list_vector_get_child(vector);
    list.child_validity = duckdb_vector_get_validity(list.child);
    return list;
}

/* MinGW's GCC treats the printf archetype as Microsoft's format rules, which
 * reject %lld/%llu. Use its C99-conforming stdio and check formats as such.
 * Clang (Windows ARM64) has no gnu_printf archetype and accepts %llu as printf. */
#if defined(__MINGW32__) && !defined(__clang__)
#define ROH_PRINTF_ARCHETYPE gnu_printf
#else
#define ROH_PRINTF_ARCHETYPE printf
#endif
static void roh_error(duckdb_function_info info, const char *format, ...)
    __attribute__((format(ROH_PRINTF_ARCHETYPE, 2, 3)));

static void roh_error(duckdb_function_info info, const char *format, ...) {
    char message[ROH_ERRLEN];
    char full[ROH_ERRLEN + sizeof(ROH_NAME) + 4];
    va_list args;
    va_start(args, format);
    vsnprintf(message, sizeof(message), format, args);
    va_end(args);
    snprintf(full, sizeof(full), "%s: %s", ROH_NAME, message);
    duckdb_scalar_function_set_error(info, full);
}

static void set_row_null(duckdb_vector output, idx_t row) {
    duckdb_vector_ensure_validity_writable(output);
    duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
}

static bool read_double_arg(duckdb_vector vector, idx_t row, double *value) {
    if (!row_valid(vector, row)) return false;
    *value = ((const double *)duckdb_vector_get_data(vector))[row];
    return true;
}

static void free_map(roh_map_t *map) {
    free(map->map_pos0);
    free(map->map_rate);
}

/* Loads the map nodes of one row. Returns 1 with *present set, or 0 after
 * reporting an error. */
static int load_map(duckdb_function_info info, const roh_args_t *args, idx_t row, roh_map_t *map,
                    size_t *count) {
    roh_list_t pos_list = open_list(args->map_pos, row);
    roh_list_t cm_list = open_list(args->map_cm, row);
    *count = 0;
    if (!pos_list.valid && !cm_list.valid) return 1;
    if (!pos_list.valid || !cm_list.valid) {
        roh_error(info, "genetic map positions and cM values must both be given or both NULL");
        return 0;
    }
    idx_t n = pos_list.entries[row].length;
    if (n != cm_list.entries[row].length) {
        roh_error(info, "genetic map positions and cM values differ in length (%llu vs %llu)",
                  (unsigned long long)n, (unsigned long long)cm_list.entries[row].length);
        return 0;
    }
    if (n == 0) {
        roh_error(info, "genetic map is empty");
        return 0;
    }
    if (n > map->map_capacity) {
        int32_t *pos0 = realloc(map->map_pos0, n * sizeof(*pos0));
        if (pos0 == NULL) goto oom;
        map->map_pos0 = pos0;
        double *rate = realloc(map->map_rate, n * sizeof(*rate));
        if (rate == NULL) goto oom;
        map->map_rate = rate;
        map->map_capacity = n;
    }
    const int64_t *positions = duckdb_vector_get_data(pos_list.child);
    const double *cm = duckdb_vector_get_data(cm_list.child);
    idx_t pos_offset = pos_list.entries[row].offset;
    idx_t cm_offset = cm_list.entries[row].offset;
    for (idx_t i = 0; i < n; i++) {
        if (!element_valid(pos_list.child_validity, pos_offset + i) ||
            !element_valid(cm_list.child_validity, cm_offset + i)) {
            roh_error(info, "genetic map entry %llu has a NULL position or cM value",
                      (unsigned long long)(i + 1));
            return 0;
        }
        int64_t pos1 = positions[pos_offset + i];
        double value = cm[cm_offset + i];
        if (pos1 < 1 || pos1 > ROH_MAX_POS1) {
            roh_error(info, "genetic map position %lld is outside 1..%d", (long long)pos1,
                      ROH_MAX_POS1);
            return 0;
        }
        if (!isfinite(value)) {
            roh_error(info, "genetic map cM values must be finite");
            return 0;
        }
        if (i > 0 && pos1 - 1 <= map->map_pos0[i - 1]) {
            roh_error(info, "genetic map positions must be strictly ascending");
            return 0;
        }
        map->map_pos0[i] = (int32_t)(pos1 - 1);
        map->map_rate[i] = value * 0.01;
    }
    *count = n;
    return 1;
oom:
    roh_error(info, "out of memory for the genetic map");
    return 0;
}

/* Decodes one row. Returns 1 when the row holds a result (segments are in roh),
 * 0 after an error, and -1 for a NULL result. */
static int decode_row(duckdb_function_info info, const roh_args_t *args, idx_t row,
                      duckhts_roh_t *roh, roh_map_t *map) {
    double hw_to_az, az_to_hw, gt_error = 0, rec_rate = 0;
    if (!read_double_arg(args->hw_to_az, row, &hw_to_az) ||
        !read_double_arg(args->az_to_hw, row, &az_to_hw)) return -1;
    if (args->has_gt && !read_double_arg(args->gt_error, row, &gt_error)) return -1;
    if (!isfinite(hw_to_az) || hw_to_az < 0 || hw_to_az > 1) {
        roh_error(info, "hw_to_az must be in [0, 1]");
        return 0;
    }
    if (!isfinite(az_to_hw) || az_to_hw < 0 || az_to_hw > 1) {
        roh_error(info, "az_to_hw must be in [0, 1]");
        return 0;
    }
    if (args->has_gt && (!isfinite(gt_error) || gt_error < 0)) {
        roh_error(info, "gt_error must be a finite phred value of at least 0");
        return 0;
    }
    if (read_double_arg(args->rec_rate, row, &rec_rate) && (!isfinite(rec_rate) || rec_rate < 0)) {
        roh_error(info, "rec_rate must be finite and at least 0");
        return 0;
    }

    roh_list_t positions = open_list(args->positions, row);
    roh_list_t af = open_list(args->af, row);
    roh_list_t evidence = open_list(args->evidence, row);
    if (!positions.valid || !af.valid || !evidence.valid) return -1;
    idx_t n = positions.entries[row].length;
    const char *evidence_name = args->has_gt ? "genotype" : "PL";
    if (af.entries[row].length != n) {
        roh_error(info, "positions and af lists differ in length (%llu vs %llu)",
                  (unsigned long long)n, (unsigned long long)af.entries[row].length);
        return 0;
    }
    if (evidence.entries[row].length != n) {
        roh_error(info, "positions and %s lists differ in length (%llu vs %llu)", evidence_name,
                  (unsigned long long)n, (unsigned long long)evidence.entries[row].length);
        return 0;
    }

    size_t map_count;
    if (!load_map(info, args, row, map, &map_count)) return 0;
    duckhts_roh_params_t params = {
        .hw_to_az = hw_to_az,
        .az_to_hw = az_to_hw,
        .rec_rate = rec_rate,
        .map_pos0 = map_count ? map->map_pos0 : NULL,
        .map_rate = map_count ? map->map_rate : NULL,
        .map_n = map_count,
    };
    if (duckhts_roh_begin(roh, &params, n) != 0) {
        roh_error(info, "out of memory for %llu sites", (unsigned long long)n);
        return 0;
    }

    const int64_t *pos = duckdb_vector_get_data(positions.child);
    const double *freq = duckdb_vector_get_data(af.child);
    idx_t pos_offset = positions.entries[row].offset;
    idx_t af_offset = af.entries[row].offset;
    idx_t ev_offset = evidence.entries[row].offset;

    /* PL evidence is a list of lists; GT evidence is a flat list. */
    const duckdb_list_entry *pl_entries = NULL;
    const int32_t *pl_values = NULL;
    const uint64_t *pl_validity = NULL;
    const int32_t *dosage = NULL;
    if (args->has_gt) {
        dosage = duckdb_vector_get_data(evidence.child);
    } else {
        pl_entries = duckdb_vector_get_data(evidence.child);
        duckdb_vector pl_child = duckdb_list_vector_get_child(evidence.child);
        pl_values = duckdb_vector_get_data(pl_child);
        pl_validity = duckdb_vector_get_validity(pl_child);
    }

    int64_t prev_pos1 = 0;
    for (idx_t i = 0; i < n; i++) {
        if (!element_valid(positions.child_validity, pos_offset + i)) {
            roh_error(info, "positions cannot contain NULL (entry %llu)", (unsigned long long)(i + 1));
            return 0;
        }
        int64_t pos1 = pos[pos_offset + i];
        if (pos1 < 1 || pos1 > ROH_MAX_POS1) {
            roh_error(info, "position %lld is outside 1..%d", (long long)pos1, ROH_MAX_POS1);
            return 0;
        }
        if (pos1 < prev_pos1) {
            roh_error(info, "positions must be sorted ascending (entry %llu is %lld after %lld)",
                      (unsigned long long)(i + 1), (long long)pos1, (long long)prev_pos1);
            return 0;
        }
        bool duplicate = pos1 == prev_pos1;
        prev_pos1 = pos1;

        /* A valid AF element outside [0, 1] is an error even on skipped sites. */
        bool has_af = element_valid(af.child_validity, af_offset + i);
        double alt_freq = has_af ? freq[af_offset + i] : 0;
        if (has_af && !isnan(alt_freq) && (alt_freq < 0 || alt_freq > 1)) {
            roh_error(info, "af must be in [0, 1] (entry %llu is %g)", (unsigned long long)(i + 1),
                      alt_freq);
            return 0;
        }

        double pdg[3];
        bool usable;
        if (args->has_gt) {
            idx_t at = ev_offset + i;
            usable = element_valid(evidence.child_validity, at);
            if (usable) {
                int32_t value = dosage[at];
                if (value < 0 || value > 2) {
                    roh_error(info, "genotype dosage must be 0, 1 or 2 (entry %llu is %d)",
                              (unsigned long long)(i + 1), value);
                    return 0;
                }
                usable = duckhts_roh_pdg_from_gt(gt_error, value, pdg);
            }
        } else {
            idx_t at = ev_offset + i;
            usable = element_valid(evidence.child_validity, at) && pl_entries[at].length == 3;
            if (usable) {
                idx_t base = pl_entries[at].offset;
                usable = element_valid(pl_validity, base) && element_valid(pl_validity, base + 1) &&
                         element_valid(pl_validity, base + 2);
                if (usable) {
                    usable = duckhts_roh_pdg_from_pl(roh, pl_values[base], pl_values[base + 1],
                                                     pl_values[base + 2], pdg);
                }
            }
        }

        /* bcftools skips repeated positions, sites without a usable AF (missing or
         * exactly 0) and sites without usable genotype evidence. */
        if (duplicate || !usable || !has_af || isnan(alt_freq) || alt_freq == 0.0) continue;
        duckhts_roh_push(roh, (int32_t)(pos1 - 1), pdg, alt_freq);
    }
    if (duckhts_roh_decode(roh) != 0) {
        roh_error(info, "out of memory while decoding");
        return 0;
    }
    return 1;
}

static bool write_segments(duckdb_function_info info, duckdb_vector output, idx_t row,
                           const duckhts_roh_t *roh) {
    size_t count;
    const duckhts_roh_segment_t *segments = duckhts_roh_segments(roh, &count);
    duckdb_list_entry entry;
    if (!duckhts_list_extend(output, count, &entry)) {
        roh_error(info, "output allocation failed");
        return false;
    }
    duckdb_vector child = duckdb_list_vector_get_child(output);
    int64_t *starts = duckdb_vector_get_data(duckdb_struct_vector_get_child(child, 0));
    int64_t *ends = duckdb_vector_get_data(duckdb_struct_vector_get_child(child, 1));
    int32_t *markers = duckdb_vector_get_data(duckdb_struct_vector_get_child(child, 2));
    double *quality = duckdb_vector_get_data(duckdb_struct_vector_get_child(child, 3));
    for (size_t i = 0; i < count; i++) {
        starts[entry.offset + i] = segments[i].start_pos1;
        ends[entry.offset + i] = segments[i].end_pos1;
        markers[entry.offset + i] = segments[i].n_markers;
        quality[entry.offset + i] = segments[i].quality;
    }
    ((duckdb_list_entry *)duckdb_vector_get_data(output))[row] = entry;
    return true;
}

static void roh_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    roh_args_t args = {0};
    idx_t columns = duckdb_data_chunk_get_column_count(input);
    args.has_gt = columns == 9;
    idx_t at = 0;
    args.positions = duckdb_data_chunk_get_vector(input, at++);
    args.af = duckdb_data_chunk_get_vector(input, at++);
    args.evidence = duckdb_data_chunk_get_vector(input, at++);
    if (args.has_gt) args.gt_error = duckdb_data_chunk_get_vector(input, at++);
    args.map_pos = duckdb_data_chunk_get_vector(input, at++);
    args.map_cm = duckdb_data_chunk_get_vector(input, at++);
    args.rec_rate = duckdb_data_chunk_get_vector(input, at++);
    args.hw_to_az = duckdb_data_chunk_get_vector(input, at++);
    args.az_to_hw = duckdb_data_chunk_get_vector(input, at++);

    duckhts_roh_t *roh = duckhts_roh_create();
    roh_map_t map = {0};
    if (roh == NULL) {
        roh_error(info, "out of memory");
        return;
    }
    duckdb_list_vector_set_size(output, 0);
    idx_t rows = duckdb_data_chunk_get_size(input);
    for (idx_t row = 0; row < rows; row++) {
        int status = decode_row(info, &args, row, roh, &map);
        if (status == 0) break;
        if (status < 0) {
            set_row_null(output, row);
            continue;
        }
        if (!write_segments(info, output, row, roh)) break;
    }
    free_map(&map);
    duckhts_roh_destroy(roh);
}

static bool register_overload(duckdb_connection connection, bool has_gt) {
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type bigint_list = duckdb_create_list_type(bigint);
    duckdb_logical_type real_list = duckdb_create_list_type(real);
    duckdb_logical_type integer_list = duckdb_create_list_type(integer);
    duckdb_logical_type pl_list = duckdb_create_list_type(integer_list);
    duckdb_logical_type member_types[4] = {bigint, bigint, integer, real};
    const char *member_names[4] = {"start", "end", "n_markers", "quality"};
    duckdb_logical_type segment = duckdb_create_struct_type(member_types, member_names, 4);
    duckdb_logical_type segments = duckdb_create_list_type(segment);

    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, ROH_NAME);
    duckdb_scalar_function_add_parameter(function, bigint_list);
    duckdb_scalar_function_add_parameter(function, real_list);
    duckdb_scalar_function_add_parameter(function, has_gt ? integer_list : pl_list);
    if (has_gt) duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_add_parameter(function, bigint_list);
    duckdb_scalar_function_add_parameter(function, real_list);
    duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_set_return_type(function, segments);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, roh_scalar);
    bool ok = duckdb_register_scalar_function(connection, function) == DuckDBSuccess;
    duckdb_destroy_scalar_function(&function);

    duckdb_destroy_logical_type(&segments);
    duckdb_destroy_logical_type(&segment);
    duckdb_destroy_logical_type(&pl_list);
    duckdb_destroy_logical_type(&integer_list);
    duckdb_destroy_logical_type(&real_list);
    duckdb_destroy_logical_type(&bigint_list);
    duckdb_destroy_logical_type(&real);
    duckdb_destroy_logical_type(&integer);
    duckdb_destroy_logical_type(&bigint);
    return ok;
}

bool register_duckhts_roh_functions(duckdb_connection connection) {
    return register_overload(connection, false) && register_overload(connection, true);
}
