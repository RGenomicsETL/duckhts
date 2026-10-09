/* __duckhts_count_cells: collects the histogram cells of one sample for a
 * count-error fit.
 *
 * Memory semantics. DuckDB builds the histogram itself, one row for each
 * (block, frequency bin, depth, alt count) cell of a sample, and can spill
 * that work. This aggregate then keeps 16 bytes per cell in native buffers
 * that DuckHTS bounds:
 *   - max_cells caps one sample, checked as rows arrive;
 *   - max_cell_bytes bounds the native histogram memory of the process: a
 *     group does not grow when the bytes held by all groups and all running
 *     fits would pass its own max_cell_bytes. Concurrent queries share the
 *     count of bytes held and each applies its own limit to it. This count
 *     is separate from the ROH site buffers' (roh_sites.c).
 * Both limits raise query errors. Finalize hands a group to DuckDB as one
 * BLOB (see count_cells.h), which is DuckDB's memory, and frees the buffer.
 * Every native buffer counts against max_cell_bytes: the partial groups of
 * each thread, the copy DuckDB's combine step makes of them, the old and new
 * block while a buffer grows, and the working memory of the fit that reads
 * the BLOB (count_error_udf.c, count_error_model.c). Cells need no order
 * here (the fit sorts them), so no sort and no ORDER BY: DuckDB's
 * ordered-aggregate wrapper does not initialise the states of a C API
 * aggregate (https://github.com/duckdb/duckdb/issues/26109).
 *
 * Measured with DuckDB 1.5.1: when a query fails inside this aggregate with
 * many threads, DuckDB destroys the states of some threads and not others
 * (14 of 14 states of a 20-thread run were never destroyed, 0 of 1 with one
 * thread), not even at process exit. The buffers of such states are leaked
 * and stay charged against max_cell_bytes, so the count stays truthful. The
 * ROH site buffers (roh_sites.c) have the same exposure.
 */
#if defined(__MINGW32__) && !defined(__USE_MINGW_ANSI_STDIO)
#define __USE_MINGW_ANSI_STDIO 1
#endif
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <math.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "include/count_cells.h"

#define CELLS_NAME "__duckhts_count_cells"
#define CELLS_ERRLEN 256
#define CELLS_MIN_CAPACITY 16u
/* A state whose cells finalize has already handed to DuckDB. */
#define CELLS_CONSUMED UINT64_MAX

/* Aggregate arguments, in order. */
enum {
    CELLS_IN_BLOCK = 0,
    CELLS_IN_AF,
    CELLS_IN_DEPTH,
    CELLS_IN_ALT,
    CELLS_IN_SITES,
    CELLS_IN_MAX_CELLS,
    CELLS_IN_MAX_CELL_BYTES,
    CELLS_IN_COUNT
};

typedef struct {
    unsigned char *buffer; /* a header, then `capacity` cells */
    uint64_t count;
    uint64_t capacity;
    uint64_t max_cells;     /* 0 before the first row; CELLS_CONSUMED after finalize */
    uint64_t max_cell_bytes;
} count_cells_state_t;

/* Bytes of native histogram memory held by this process. */
static atomic_uint_fast64_t count_cell_bytes_held;

int duckhts_count_bytes_charge(uint64_t bytes, uint64_t max_cell_bytes, uint64_t *would_hold) {
    uint64_t held = atomic_fetch_add(&count_cell_bytes_held, bytes) + bytes;
    if (held > max_cell_bytes) {
        atomic_fetch_sub(&count_cell_bytes_held, bytes);
        *would_hold = held;
        return 0;
    }
    return 1;
}

void duckhts_count_bytes_release(uint64_t bytes) {
    atomic_fetch_sub(&count_cell_bytes_held, bytes);
}

static void cells_error(duckdb_function_info info, const char *message) {
    char full[CELLS_ERRLEN + 32];
    snprintf(full, sizeof(full), "duckhts_count_error_fit: %s", message);
    duckdb_aggregate_function_set_error(info, full);
}

static bool input_valid(duckdb_vector vector, idx_t row) {
    uint64_t *mask = duckdb_vector_get_validity(vector);
    return !mask || duckdb_validity_row_is_valid(mask, row);
}

static uint64_t buffer_bytes(uint64_t capacity) {
    return sizeof(duckhts_count_cells_header_t) + capacity * sizeof(duckhts_count_cell_t);
}

static duckhts_count_cell_t *cells_of(unsigned char *buffer) {
    return (duckhts_count_cell_t *)(buffer + sizeof(duckhts_count_cells_header_t));
}

static void release_state(count_cells_state_t *state) {
    if (state->buffer != NULL) {
        duckhts_count_bytes_release(buffer_bytes(state->capacity));
        free(state->buffer);
    }
    state->buffer = NULL;
    state->count = 0;
    state->capacity = 0;
}

/* Makes room for `need` cells. Returns false with `error` written when a limit
 * is reached or the allocation fails; the state is then unchanged. */
static bool reserve_cells(count_cells_state_t *state, uint64_t need, char *error) {
    if (need <= state->capacity) return true;
    if (need > state->max_cells) {
        snprintf(error, CELLS_ERRLEN, "one sample holds more than %llu histogram cells (max_cells)",
                 (unsigned long long)state->max_cells);
        return false;
    }
    uint64_t capacity = state->capacity ? state->capacity * 2 : CELLS_MIN_CAPACITY;
    if (capacity < need) capacity = need;
    if (capacity > state->max_cells) capacity = state->max_cells;

    /* The new block is reserved in full before it exists and the old one is
     * released only after the copy, so the budget covers the moment both live. */
    uint64_t old_bytes = state->buffer ? buffer_bytes(state->capacity) : 0;
    uint64_t new_bytes = buffer_bytes(capacity);
    uint64_t would_hold = 0;
    if (!duckhts_count_bytes_charge(new_bytes, state->max_cell_bytes, &would_hold)) {
        snprintf(error, CELLS_ERRLEN,
                 "histogram buffers would hold %llu bytes, more than max_cell_bytes = %llu; fit "
                 "fewer samples per query or raise max_cell_bytes",
                 (unsigned long long)would_hold, (unsigned long long)state->max_cell_bytes);
        return false;
    }
    unsigned char *grown = new_bytes <= SIZE_MAX ? malloc((size_t)new_bytes) : NULL;
    if (grown == NULL) {
        duckhts_count_bytes_release(new_bytes);
        snprintf(error, CELLS_ERRLEN, "out of memory for %llu cells", (unsigned long long)capacity);
        return false;
    }
    if (state->buffer != NULL) {
        memcpy(grown, state->buffer, (size_t)buffer_bytes(state->count));
        free(state->buffer);
        duckhts_count_bytes_release(old_bytes);
    }
    state->buffer = grown;
    state->capacity = capacity;
    return true;
}

static idx_t cells_state_size(duckdb_function_info info) {
    (void)info;
    return (idx_t)sizeof(count_cells_state_t);
}

static void cells_state_init(duckdb_function_info info, duckdb_aggregate_state state) {
    (void)info;
    if (state != NULL) memset(state, 0, sizeof(count_cells_state_t));
}

static void cells_state_destroy(duckdb_aggregate_state *states, idx_t count) {
    if (states == NULL) return;
    for (idx_t i = 0; i < count; i++) {
        if (states[i] != NULL) release_state((count_cells_state_t *)states[i]);
    }
}

/* Sets the limits of a state from its first row or first combined source, and
 * checks that later ones agree. */
static bool adopt_limits(count_cells_state_t *state, uint64_t max_cells, uint64_t max_cell_bytes,
                         char *error) {
    if (state->max_cells == CELLS_CONSUMED) {
        snprintf(error, CELLS_ERRLEN, "internal error: histogram was already finalized");
        return false;
    }
    if (state->max_cells == 0) {
        state->max_cells = max_cells;
        state->max_cell_bytes = max_cell_bytes;
        return true;
    }
    if (state->max_cells != max_cells || state->max_cell_bytes != max_cell_bytes) {
        snprintf(error, CELLS_ERRLEN, "internal error: limits change within a group");
        return false;
    }
    return true;
}

static void cells_update(duckdb_function_info info, duckdb_data_chunk input,
                         duckdb_aggregate_state *states) {
    duckdb_vector vectors[CELLS_IN_COUNT];
    char error[CELLS_ERRLEN];
    if (states == NULL) {
        cells_error(info, "internal error: aggregate state is missing");
        return;
    }
    for (idx_t i = 0; i < CELLS_IN_COUNT; i++) vectors[i] = duckdb_data_chunk_get_vector(input, i);
    const int32_t *blocks = duckdb_vector_get_data(vectors[CELLS_IN_BLOCK]);
    const double *frequencies = duckdb_vector_get_data(vectors[CELLS_IN_AF]);
    const int32_t *depths = duckdb_vector_get_data(vectors[CELLS_IN_DEPTH]);
    const int32_t *alts = duckdb_vector_get_data(vectors[CELLS_IN_ALT]);
    const int64_t *sites = duckdb_vector_get_data(vectors[CELLS_IN_SITES]);
    const int64_t *max_cells = duckdb_vector_get_data(vectors[CELLS_IN_MAX_CELLS]);
    const int64_t *max_cell_bytes = duckdb_vector_get_data(vectors[CELLS_IN_MAX_CELL_BYTES]);

    idx_t rows = duckdb_data_chunk_get_size(input);
    for (idx_t row = 0; row < rows; row++) {
        count_cells_state_t *state = (count_cells_state_t *)states[row];
        duckhts_count_cell_t cell;
        if (state == NULL) {
            cells_error(info, "internal error: aggregate state is missing");
            return;
        }
        bool has_max_cells = input_valid(vectors[CELLS_IN_MAX_CELLS], row);
        bool has_max_cell_bytes = input_valid(vectors[CELLS_IN_MAX_CELL_BYTES], row);
        const char *limit_error = duckhts_count_check_limits(
            has_max_cells, has_max_cells ? max_cells[row] : 0,
            has_max_cell_bytes, has_max_cell_bytes ? max_cell_bytes[row] : 0);
        if (limit_error != NULL) {
            cells_error(info, limit_error);
            return;
        }
        if (!adopt_limits(state, (uint64_t)max_cells[row], (uint64_t)max_cell_bytes[row], error)) {
            cells_error(info, error);
            return;
        }
        /* The macro makes every cell from whole rows, so a NULL here is its error. */
        for (int field = CELLS_IN_BLOCK; field <= CELLS_IN_SITES; field++) {
            if (!input_valid(vectors[field], row)) {
                cells_error(info, "internal error: a histogram cell has a NULL field");
                return;
            }
        }
        if (blocks[row] < 0 || blocks[row] > DUCKHTS_COUNT_MAX_BLOCKS) {
            snprintf(error, sizeof(error), "a sample has more than %d blocks; raise block_bases",
                     DUCKHTS_COUNT_MAX_BLOCKS + 1);
            cells_error(info, error);
            return;
        }
        if (depths[row] < 1 || depths[row] > DUCKHTS_COUNT_MAX_DEPTH || alts[row] < 0 || alts[row] > depths[row]) {
            cells_error(info, "internal error: a histogram cell has a depth or count out of range");
            return;
        }
        if (!(frequencies[row] > 0 && frequencies[row] < 1)) {
            cells_error(info, "internal error: a histogram cell has a frequency outside (0, 1)");
            return;
        }
        if (sites[row] < 1 || sites[row] > UINT32_MAX) {
            cells_error(info, "internal error: a histogram cell has a site count out of range");
            return;
        }
        if (!reserve_cells(state, state->count + 1, error)) {
            cells_error(info, error);
            return;
        }
        cell.sites = (uint32_t)sites[row];
        cell.af = (float)frequencies[row];
        /* A double inside (0, 1) can round to 0 or 1 as a float; the cell
         * keeps the nearest float inside, as count_cells.h promises. */
        if (cell.af >= 1.0f) cell.af = nextafterf(1.0f, 0.0f);
        if (cell.af <= 0.0f) cell.af = nextafterf(0.0f, 1.0f);
        cell.block = (uint16_t)blocks[row];
        cell.depth = (uint16_t)depths[row];
        cell.alt = (uint16_t)alts[row];
        cell.reserved = 0;
        memcpy(&cells_of(state->buffer)[state->count], &cell, sizeof(cell));
        state->count++;
    }
}

static void cells_combine(duckdb_function_info info, duckdb_aggregate_state *source,
                          duckdb_aggregate_state *target, idx_t count) {
    char error[CELLS_ERRLEN];
    if (source == NULL || target == NULL) {
        cells_error(info, "internal error: aggregate state is missing during combine");
        return;
    }
    for (idx_t i = 0; i < count; i++) {
        const count_cells_state_t *from = (const count_cells_state_t *)source[i];
        count_cells_state_t *to = (count_cells_state_t *)target[i];
        if (from == NULL || to == NULL) {
            cells_error(info, "internal error: aggregate state is missing during combine");
            return;
        }
        if (from->max_cells == 0 || from->count == 0) continue;
        if (from->max_cells == CELLS_CONSUMED) {
            cells_error(info, "internal error: histogram was already finalized");
            return;
        }
        if (!adopt_limits(to, from->max_cells, from->max_cell_bytes, error) ||
            !reserve_cells(to, to->count + from->count, error)) {
            cells_error(info, error);
            return;
        }
        memcpy(cells_of(to->buffer) + to->count, cells_of(from->buffer),
               (size_t)(from->count * sizeof(duckhts_count_cell_t)));
        to->count += from->count;
    }
}

static void cells_finalize(duckdb_function_info info, duckdb_aggregate_state *source,
                           duckdb_vector result, idx_t count, idx_t offset) {
    if (source == NULL) {
        cells_error(info, "internal error: aggregate state is missing during finalize");
        return;
    }
    duckdb_vector_ensure_validity_writable(result);
    uint64_t *validity = duckdb_vector_get_validity(result);
    for (idx_t i = 0; i < count; i++) {
        count_cells_state_t *state = (count_cells_state_t *)source[i];
        idx_t row = offset + i;
        duckhts_count_cells_header_t header = {.reserved = 0};
        uint64_t bytes;
        if (state == NULL) {
            cells_error(info, "internal error: aggregate state is missing during finalize");
            return;
        }
        if (state->max_cells == CELLS_CONSUMED) {
            cells_error(info, "internal error: histogram was already finalized");
            return;
        }
        if (state->max_cells == 0 || state->count == 0) {
            duckdb_validity_set_row_invalid(validity, row);
            continue;
        }
        memcpy(header.magic, DUCKHTS_COUNT_CELLS_MAGIC, DUCKHTS_COUNT_CELLS_MAGIC_LENGTH);
        memcpy(state->buffer, &header, sizeof(header));
        bytes = buffer_bytes(state->count);
        if (bytes > UINT32_MAX) {
            cells_error(info, "internal error: histogram exceeds a BLOB");
            return;
        }
        duckdb_vector_assign_string_element_len(result, row, (const char *)state->buffer, (idx_t)bytes);
        /* DuckDB now owns a copy; the group is not held twice. */
        release_state(state);
        state->max_cells = CELLS_CONSUMED;
    }
}

bool register_duckhts_count_cells(duckdb_connection connection) {
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type blob = duckdb_create_logical_type(DUCKDB_TYPE_BLOB);
    bool ok;

    duckdb_aggregate_function function = duckdb_create_aggregate_function();
    duckdb_aggregate_function_set_name(function, CELLS_NAME);
    duckdb_aggregate_function_add_parameter(function, integer); /* block */
    duckdb_aggregate_function_add_parameter(function, real);    /* af */
    duckdb_aggregate_function_add_parameter(function, integer); /* depth */
    duckdb_aggregate_function_add_parameter(function, integer); /* alt */
    duckdb_aggregate_function_add_parameter(function, bigint);  /* sites */
    duckdb_aggregate_function_add_parameter(function, bigint);  /* max_cells */
    duckdb_aggregate_function_add_parameter(function, bigint);  /* max_cell_bytes */
    duckdb_aggregate_function_set_return_type(function, blob);
    duckdb_aggregate_function_set_special_handling(function);
    duckdb_aggregate_function_set_functions(function, cells_state_size, cells_state_init,
                                            cells_update, cells_combine, cells_finalize);
    duckdb_aggregate_function_set_destructor(function, cells_state_destroy);
    ok = duckdb_register_aggregate_function(connection, function) == DuckDBSuccess;
    duckdb_destroy_aggregate_function(&function);

    duckdb_destroy_logical_type(&blob);
    duckdb_destroy_logical_type(&real);
    duckdb_destroy_logical_type(&integer);
    duckdb_destroy_logical_type(&bigint);
    return ok;
}
