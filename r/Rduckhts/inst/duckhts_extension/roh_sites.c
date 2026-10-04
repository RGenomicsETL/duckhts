/* __duckhts_roh_sites: collects the sites of one sample and chromosome for a
 * runs-of-homozygosity decode.
 *
 * Memory semantics. The ROH macros used to build one DuckDB list per sample
 * and chromosome, about 240 bytes per site, which DuckDB could neither cap nor
 * spill. This aggregate keeps 16 bytes per site (24 for read counts) in native
 * buffers that DuckHTS bounds itself:
 *   - max_sites caps one sample and chromosome, checked as rows arrive;
 *   - max_site_bytes caps the buffers of every live group in the process, so
 *     concurrent decodes share it.
 * Both limits raise query errors. Finalize sorts a group by position, hands it
 * to DuckDB as one BLOB (see roh_sites.h) and frees the buffer. Every native
 * buffer counts against max_site_bytes: the partial groups of each thread, the
 * copy DuckDB's combine step makes of them, the old and new block while a
 * buffer grows, and the scratch copy that sorting a group that arrived out of
 * order needs.
 *
 * The sort is stable: sites of one position keep their arrival order, so the
 * first record of a repeated position stays first, as bcftools sees the file.
 *
 * Never call this aggregate with ORDER BY: DuckDB's ordered-aggregate wrapper
 * does not initialise the states of a C API aggregate
 * (https://github.com/duckdb/duckdb/issues/26109). The sort here makes it
 * unnecessary.
 */
#if defined(__MINGW32__) && !defined(__USE_MINGW_ANSI_STDIO)
#define __USE_MINGW_ANSI_STDIO 1
#endif
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <math.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "roh_hmm.h"
#include "roh_sites.h"

#define SITES_NAME "__duckhts_roh_sites"
#define SITES_ERRLEN 256
#define SITES_MIN_CAPACITY 16u
#define SITES_MAX_POS1 INT32_MAX
/* Kind of a state whose records finalize has already handed to DuckDB. */
#define SITES_CONSUMED UINT32_MAX

_Static_assert(offsetof(duckhts_roh_site_t, pos1) == offsetof(duckhts_roh_count_site_t, pos1),
               "both site records keep the position at one offset");

/* Aggregate arguments, in order. e1, e2 and e3 are the evidence of the kind:
 * the RR, RA and AA likelihoods; the dosage; or the other-allele and
 * counted-allele read counts. */
enum {
    SITES_IN_POS = 0,
    SITES_IN_AF,
    SITES_IN_E1,
    SITES_IN_E2,
    SITES_IN_E3,
    SITES_IN_KIND,
    SITES_IN_MAX_SITES,
    SITES_IN_MAX_SITE_BYTES,
    SITES_IN_COUNT
};

typedef struct {
    unsigned char *buffer; /* a header, then `capacity` records of the kind */
    uint64_t count;
    uint64_t capacity;
    uint64_t max_sites;
    uint64_t max_site_bytes;
    uint32_t kind; /* duckhts_roh_sites_kind_t; 0 before the first row */
} roh_sites_state_t;

/* Bytes of site buffers held by the live states of this process. */
static atomic_uint_fast64_t roh_site_bytes_held;

static void sites_error(duckdb_function_info info, const char *message) {
    char full[SITES_ERRLEN + 16];
    snprintf(full, sizeof(full), "duckhts_roh: %s", message);
    duckdb_aggregate_function_set_error(info, full);
}

static bool input_valid(duckdb_vector vector, idx_t row) {
    uint64_t *mask = duckdb_vector_get_validity(vector);
    return !mask || duckdb_validity_row_is_valid(mask, row);
}

static uint64_t buffer_bytes(uint32_t kind, uint64_t capacity) {
    return sizeof(duckhts_roh_sites_header_t) + capacity * duckhts_roh_site_size(kind);
}

static unsigned char *records_of(unsigned char *buffer) {
    return buffer + sizeof(duckhts_roh_sites_header_t);
}

static void release_state(roh_sites_state_t *state) {
    if (state->buffer != NULL) {
        atomic_fetch_sub(&roh_site_bytes_held, buffer_bytes(state->kind, state->capacity));
        free(state->buffer);
    }
    state->buffer = NULL;
    state->count = 0;
    state->capacity = 0;
}

/* Makes room for `need` records. Returns false with `error` written when a
 * limit is reached or the allocation fails; the state is then unchanged. */
static bool reserve_sites(roh_sites_state_t *state, uint64_t need, char *error) {
    if (need <= state->capacity) return true;
    if (need > state->max_sites) {
        snprintf(error, SITES_ERRLEN,
                 "one sample and chromosome holds more than %llu sites (max_sites)",
                 (unsigned long long)state->max_sites);
        return false;
    }
    uint64_t capacity = state->capacity ? state->capacity * 2 : SITES_MIN_CAPACITY;
    if (capacity < need) capacity = need;
    if (capacity > state->max_sites) capacity = state->max_sites;

    /* The new block is reserved in full before it exists and the old one is
     * released only after the copy, so the budget covers the moment both live. */
    uint64_t old_bytes = state->buffer ? buffer_bytes(state->kind, state->capacity) : 0;
    uint64_t new_bytes = buffer_bytes(state->kind, capacity);
    uint64_t held = atomic_fetch_add(&roh_site_bytes_held, new_bytes) + new_bytes;
    if (held > state->max_site_bytes) {
        atomic_fetch_sub(&roh_site_bytes_held, new_bytes);
        snprintf(error, SITES_ERRLEN,
                 "site buffers would hold %llu bytes, more than max_site_bytes = %llu; decode "
                 "fewer samples per query or raise max_site_bytes",
                 (unsigned long long)held, (unsigned long long)state->max_site_bytes);
        return false;
    }
    unsigned char *grown = new_bytes <= SIZE_MAX ? malloc((size_t)new_bytes) : NULL;
    if (grown == NULL) {
        atomic_fetch_sub(&roh_site_bytes_held, new_bytes);
        snprintf(error, SITES_ERRLEN, "out of memory for %llu sites", (unsigned long long)capacity);
        return false;
    }
    if (state->buffer != NULL) {
        memcpy(grown, state->buffer,
               (size_t)buffer_bytes(state->kind, state->count));
        free(state->buffer);
        atomic_fetch_sub(&roh_site_bytes_held, old_bytes);
    }
    state->buffer = grown;
    state->capacity = capacity;
    return true;
}

static idx_t sites_state_size(duckdb_function_info info) {
    (void)info;
    return (idx_t)sizeof(roh_sites_state_t);
}

static void sites_state_init(duckdb_function_info info, duckdb_aggregate_state state) {
    (void)info;
    if (state != NULL) memset(state, 0, sizeof(roh_sites_state_t));
}

static void sites_state_destroy(duckdb_aggregate_state *states, idx_t count) {
    if (states == NULL) return;
    for (idx_t i = 0; i < count; i++) {
        if (states[i] != NULL) release_state((roh_sites_state_t *)states[i]);
    }
}

/* Sets the kind and limits of a state from its first row or first combined
 * source, and checks that later ones agree. */
static bool adopt_settings(roh_sites_state_t *state, uint32_t kind, uint64_t max_sites,
                           uint64_t max_site_bytes, char *error) {
    if (state->kind == SITES_CONSUMED) {
        snprintf(error, SITES_ERRLEN, "internal error: site list was already finalized");
        return false;
    }
    if (state->kind == 0) {
        state->kind = kind;
        state->max_sites = max_sites;
        state->max_site_bytes = max_site_bytes;
        return true;
    }
    if (state->kind != kind || state->max_sites != max_sites ||
        state->max_site_bytes != max_site_bytes) {
        snprintf(error, SITES_ERRLEN, "internal error: evidence kind or limits change within a group");
        return false;
    }
    return true;
}

static void sites_update(duckdb_function_info info, duckdb_data_chunk input,
                         duckdb_aggregate_state *states) {
    duckdb_vector vectors[SITES_IN_COUNT];
    char error[SITES_ERRLEN];
    if (states == NULL) {
        sites_error(info, "internal error: aggregate state is missing");
        return;
    }
    for (idx_t i = 0; i < SITES_IN_COUNT; i++) vectors[i] = duckdb_data_chunk_get_vector(input, i);
    const int64_t *positions = duckdb_vector_get_data(vectors[SITES_IN_POS]);
    const double *frequencies = duckdb_vector_get_data(vectors[SITES_IN_AF]);
    const int32_t *e1 = duckdb_vector_get_data(vectors[SITES_IN_E1]);
    const int32_t *e2 = duckdb_vector_get_data(vectors[SITES_IN_E2]);
    const int32_t *e3 = duckdb_vector_get_data(vectors[SITES_IN_E3]);
    const int32_t *kinds = duckdb_vector_get_data(vectors[SITES_IN_KIND]);
    const int64_t *max_sites = duckdb_vector_get_data(vectors[SITES_IN_MAX_SITES]);
    const int64_t *max_site_bytes = duckdb_vector_get_data(vectors[SITES_IN_MAX_SITE_BYTES]);

    idx_t rows = duckdb_data_chunk_get_size(input);
    for (idx_t row = 0; row < rows; row++) {
        roh_sites_state_t *state = (roh_sites_state_t *)states[row];
        if (state == NULL) {
            sites_error(info, "internal error: aggregate state is missing");
            return;
        }
        bool has_max_sites = input_valid(vectors[SITES_IN_MAX_SITES], row);
        bool has_max_site_bytes = input_valid(vectors[SITES_IN_MAX_SITE_BYTES], row);
        const char *limit_error = duckhts_roh_check_limits(
            has_max_sites, has_max_sites ? max_sites[row] : 0,
            has_max_site_bytes, has_max_site_bytes ? max_site_bytes[row] : 0);
        if (limit_error != NULL) {
            sites_error(info, limit_error);
            return;
        }
        if (!input_valid(vectors[SITES_IN_KIND], row) || duckhts_roh_site_size((uint32_t)kinds[row]) == 0) {
            sites_error(info, "internal error: unknown evidence kind");
            return;
        }
        uint32_t kind = (uint32_t)kinds[row];
        if (!adopt_settings(state, kind, (uint64_t)max_sites[row], (uint64_t)max_site_bytes[row], error)) {
            sites_error(info, error);
            return;
        }

        if (!input_valid(vectors[SITES_IN_POS], row)) {
            sites_error(info, "positions cannot contain NULL");
            return;
        }
        int64_t pos1 = positions[row];
        if (pos1 < 1 || pos1 > SITES_MAX_POS1) {
            snprintf(error, sizeof(error), "position %lld is outside 1..%d", (long long)pos1,
                     SITES_MAX_POS1);
            sites_error(info, error);
            return;
        }
        /* A frequency outside [0, 1] is an error even on a site that is skipped. */
        double af = input_valid(vectors[SITES_IN_AF], row) ? frequencies[row] : NAN;
        if (!isnan(af) && (af < 0 || af > 1)) {
            snprintf(error, sizeof(error), "af must be in [0, 1] (position %lld is %g)",
                     (long long)pos1, af);
            sites_error(info, error);
            return;
        }
        bool has_e1 = input_valid(vectors[SITES_IN_E1], row);
        bool has_e2 = input_valid(vectors[SITES_IN_E2], row);
        bool has_e3 = input_valid(vectors[SITES_IN_E3], row);
        /* Each count present is checked on its own, so a negative count is an
         * error even when the other count of its site is NULL. */
        if (kind == DUCKHTS_ROH_SITES_COUNTS && ((has_e1 && e1[row] < 0) || (has_e2 && e2[row] < 0))) {
            snprintf(error, sizeof(error), "read counts must be at least 0 (position %lld)",
                     (long long)pos1);
            sites_error(info, error);
            return;
        }
        if (kind == DUCKHTS_ROH_SITES_GT && has_e1 && (e1[row] < 0 || e1[row] > 2)) {
            snprintf(error, sizeof(error), "genotype dosage must be 0, 1 or 2 (position %lld is %d)",
                     (long long)pos1, e1[row]);
            sites_error(info, error);
            return;
        }

        if (!reserve_sites(state, state->count + 1, error)) {
            sites_error(info, error);
            return;
        }
        unsigned char *slot = records_of(state->buffer) + state->count * duckhts_roh_site_size(kind);
        if (kind == DUCKHTS_ROH_SITES_COUNTS) {
            duckhts_roh_count_site_t site = {.af = af, .pos1 = (uint32_t)pos1};
            site.usable = has_e1 && has_e2;
            if (site.usable) {
                site.other = e1[row];
                site.counted = e2[row];
            }
            memcpy(slot, &site, sizeof(site));
        } else if (kind == DUCKHTS_ROH_SITES_GT) {
            duckhts_roh_site_t site = {.af = af, .pos1 = (uint32_t)pos1};
            site.usable = has_e1;
            if (site.usable) site.evidence[0] = (uint8_t)e1[row];
            memcpy(slot, &site, sizeof(site));
        } else {
            duckhts_roh_site_t site = {.af = af, .pos1 = (uint32_t)pos1};
            /* Usability is decided on the likelihoods as given; the cap comes after. */
            site.usable = has_e1 && has_e2 && has_e3 && duckhts_roh_pl_usable(e1[row], e2[row], e3[row]);
            if (site.usable) {
                int32_t cap = DUCKHTS_ROH_PHRED_TABLE - 1;
                site.evidence[0] = (uint8_t)(e1[row] < cap ? e1[row] : cap);
                site.evidence[1] = (uint8_t)(e2[row] < cap ? e2[row] : cap);
                site.evidence[2] = (uint8_t)(e3[row] < cap ? e3[row] : cap);
            }
            memcpy(slot, &site, sizeof(site));
        }
        state->count++;
    }
}

static void sites_combine(duckdb_function_info info, duckdb_aggregate_state *source,
                          duckdb_aggregate_state *target, idx_t count) {
    char error[SITES_ERRLEN];
    if (source == NULL || target == NULL) {
        sites_error(info, "internal error: aggregate state is missing during combine");
        return;
    }
    for (idx_t i = 0; i < count; i++) {
        const roh_sites_state_t *from = (const roh_sites_state_t *)source[i];
        roh_sites_state_t *to = (roh_sites_state_t *)target[i];
        if (from == NULL || to == NULL) {
            sites_error(info, "internal error: aggregate state is missing during combine");
            return;
        }
        if (from->kind == 0 || from->count == 0) continue;
        if (from->kind == SITES_CONSUMED) {
            sites_error(info, "internal error: site list was already finalized");
            return;
        }
        if (!adopt_settings(to, from->kind, from->max_sites, from->max_site_bytes, error)) {
            sites_error(info, error);
            return;
        }
        if (!reserve_sites(to, to->count + from->count, error)) {
            sites_error(info, error);
            return;
        }
        uint64_t size = duckhts_roh_site_size(to->kind);
        memcpy(records_of(to->buffer) + to->count * size, records_of(from->buffer),
               (size_t)(from->count * size));
        to->count += from->count;
    }
}

static uint32_t record_pos(const unsigned char *records, uint64_t index, uint64_t size) {
    uint32_t pos1;
    memcpy(&pos1, records + index * size + offsetof(duckhts_roh_site_t, pos1), sizeof(pos1));
    return pos1;
}

/* End of the run of non-decreasing positions that starts at `start`. */
static uint64_t run_end(const unsigned char *records, uint64_t start, uint64_t count, uint64_t size) {
    uint64_t at = start + 1;
    while (at < count && record_pos(records, at - 1, size) <= record_pos(records, at, size)) at++;
    return at;
}

/* Merges the runs of `from` pairwise into `to` and returns how many runs it
 * read. Records of one position keep their order: a tie takes the left run. */
static uint64_t merge_runs(const unsigned char *from, unsigned char *to, uint64_t count,
                           uint64_t size) {
    uint64_t runs = 0;
    uint64_t at = 0;
    uint64_t out = 0;
    while (at < count) {
        uint64_t middle = run_end(from, at, count, size);
        uint64_t end = middle < count ? run_end(from, middle, count, size) : middle;
        uint64_t left = at;
        uint64_t right = middle;
        while (left < middle && right < end) {
            if (record_pos(from, right, size) < record_pos(from, left, size)) {
                memcpy(to + out * size, from + right * size, (size_t)size);
                right++;
            } else {
                memcpy(to + out * size, from + left * size, (size_t)size);
                left++;
            }
            out++;
        }
        memcpy(to + out * size, from + left * size, (size_t)((middle - left) * size));
        out += middle - left;
        memcpy(to + out * size, from + right * size, (size_t)((end - right) * size));
        out += end - right;
        runs += middle < count ? 2 : 1;
        at = end;
    }
    return runs;
}

/* Sorts the records of a state by position, keeping the arrival order of equal
 * positions. The scratch copy counts against max_site_bytes. Returns false
 * with `error` written when it exceeds the limit or cannot be allocated. */
static bool sort_sites(roh_sites_state_t *state, char *error) {
    uint64_t size = duckhts_roh_site_size(state->kind);
    unsigned char *records = records_of(state->buffer);
    if (run_end(records, 0, state->count, size) == state->count) return true;

    uint64_t scratch_bytes = buffer_bytes(state->kind, state->count);
    uint64_t held = atomic_fetch_add(&roh_site_bytes_held, scratch_bytes) + scratch_bytes;
    if (held > state->max_site_bytes) {
        atomic_fetch_sub(&roh_site_bytes_held, scratch_bytes);
        snprintf(error, SITES_ERRLEN,
                 "sorting one sample and chromosome would bring the site buffers to %llu bytes, "
                 "more than max_site_bytes = %llu; decode fewer samples per query or raise "
                 "max_site_bytes",
                 (unsigned long long)held, (unsigned long long)state->max_site_bytes);
        return false;
    }
    unsigned char *scratch = scratch_bytes <= SIZE_MAX ? malloc((size_t)scratch_bytes) : NULL;
    if (scratch == NULL) {
        atomic_fetch_sub(&roh_site_bytes_held, scratch_bytes);
        snprintf(error, SITES_ERRLEN,
                 "out of memory while sorting the sites of one sample and chromosome");
        return false;
    }
    unsigned char *from = records;
    unsigned char *to = records_of(scratch);
    uint64_t runs;
    do {
        runs = merge_runs(from, to, state->count, size);
        unsigned char *swap = from;
        from = to;
        to = swap;
    } while (runs > 2);

    /* `from` now holds the sorted records. Keep the buffer that owns them. */
    if (from == records) {
        free(scratch);
        atomic_fetch_sub(&roh_site_bytes_held, scratch_bytes);
        return true;
    }
    /* The scratch copy is already counted; the old buffer goes. */
    atomic_fetch_sub(&roh_site_bytes_held, buffer_bytes(state->kind, state->capacity));
    free(state->buffer);
    state->buffer = scratch;
    state->capacity = state->count;
    return true;
}

static void sites_finalize(duckdb_function_info info, duckdb_aggregate_state *source,
                           duckdb_vector result, idx_t count, idx_t offset) {
    if (source == NULL) {
        sites_error(info, "internal error: aggregate state is missing during finalize");
        return;
    }
    duckdb_vector_ensure_validity_writable(result);
    uint64_t *validity = duckdb_vector_get_validity(result);
    for (idx_t i = 0; i < count; i++) {
        roh_sites_state_t *state = (roh_sites_state_t *)source[i];
        idx_t row = offset + i;
        if (state == NULL) {
            sites_error(info, "internal error: aggregate state is missing during finalize");
            return;
        }
        if (state->kind == SITES_CONSUMED) {
            sites_error(info, "internal error: site list was already finalized");
            return;
        }
        if (state->kind == 0 || state->count == 0) {
            duckdb_validity_set_row_invalid(validity, row);
            continue;
        }
        char error[SITES_ERRLEN];
        if (!sort_sites(state, error)) {
            sites_error(info, error);
            return;
        }
        duckhts_roh_sites_header_t header = {.kind = (uint8_t)state->kind};
        memcpy(header.magic, DUCKHTS_ROH_SITES_MAGIC, DUCKHTS_ROH_SITES_MAGIC_LENGTH);
        memcpy(state->buffer, &header, sizeof(header));
        uint64_t bytes = buffer_bytes(state->kind, state->count);
        if (bytes > UINT32_MAX) {
            sites_error(info, "internal error: site list exceeds a BLOB");
            return;
        }
        duckdb_vector_assign_string_element_len(result, row, (const char *)state->buffer, (idx_t)bytes);
        /* DuckDB now owns a copy; the group is not held twice. */
        release_state(state);
        state->kind = SITES_CONSUMED;
    }
}

bool register_duckhts_roh_sites(duckdb_connection connection) {
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type blob = duckdb_create_logical_type(DUCKDB_TYPE_BLOB);

    duckdb_aggregate_function function = duckdb_create_aggregate_function();
    duckdb_aggregate_function_set_name(function, SITES_NAME);
    duckdb_aggregate_function_add_parameter(function, bigint);  /* pos */
    duckdb_aggregate_function_add_parameter(function, real);    /* af */
    duckdb_aggregate_function_add_parameter(function, integer); /* e1 */
    duckdb_aggregate_function_add_parameter(function, integer); /* e2 */
    duckdb_aggregate_function_add_parameter(function, integer); /* e3 */
    duckdb_aggregate_function_add_parameter(function, integer); /* kind */
    duckdb_aggregate_function_add_parameter(function, bigint);  /* max_sites */
    duckdb_aggregate_function_add_parameter(function, bigint);  /* max_site_bytes */
    duckdb_aggregate_function_set_return_type(function, blob);
    duckdb_aggregate_function_set_special_handling(function);
    duckdb_aggregate_function_set_functions(function, sites_state_size, sites_state_init,
                                            sites_update, sites_combine, sites_finalize);
    duckdb_aggregate_function_set_destructor(function, sites_state_destroy);
    bool ok = duckdb_register_aggregate_function(connection, function) == DuckDBSuccess;
    duckdb_destroy_aggregate_function(&function);

    duckdb_destroy_logical_type(&blob);
    duckdb_destroy_logical_type(&real);
    duckdb_destroy_logical_type(&integer);
    duckdb_destroy_logical_type(&bigint);
    return ok;
}
