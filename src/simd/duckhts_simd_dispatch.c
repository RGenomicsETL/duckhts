#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && !defined(__STDC_NO_ATOMICS__) && \
    !defined(__EMSCRIPTEN__) && !defined(DUCKDB_WASM_EXTENSION)
#include <stdatomic.h>
#define DUCKHTS_SIMD_HAVE_C11_ATOMICS 1
#else
#define DUCKHTS_SIMD_HAVE_C11_ATOMICS 0
#endif

#include "duckhts_simd_internal.h"

#if DUCKHTS_SIMD_HAVE_C11_ATOMICS
static _Atomic(const duckhts_simd_ops_t *) duckhts_simd_ops;
static _Atomic(const char *) duckhts_simd_requested;
static const duckhts_simd_ops_t *duckhts_simd_load_ops(void) {
    return atomic_load_explicit(&duckhts_simd_ops, memory_order_acquire);
}
static void duckhts_simd_store_ops(const duckhts_simd_ops_t *ops) {
    atomic_store_explicit(&duckhts_simd_ops, ops, memory_order_release);
}
static const char *duckhts_simd_load_requested(void) {
    return atomic_load_explicit(&duckhts_simd_requested, memory_order_acquire);
}
static void duckhts_simd_store_requested(const char *name) {
    atomic_store_explicit(&duckhts_simd_requested, name, memory_order_release);
}
#else
static const duckhts_simd_ops_t *duckhts_simd_ops;
static const char *duckhts_simd_requested;
static const duckhts_simd_ops_t *duckhts_simd_load_ops(void) { return duckhts_simd_ops; }
static void duckhts_simd_store_ops(const duckhts_simd_ops_t *ops) { duckhts_simd_ops = ops; }
static const char *duckhts_simd_load_requested(void) { return duckhts_simd_requested; }
static void duckhts_simd_store_requested(const char *name) { duckhts_simd_requested = name; }
#endif

static int duckhts_simd_initialized = 0;

static int duckhts_backend_name_equal(const char *a, const char *b) {
    return a && b && strcmp(a, b) == 0;
}

static const char *duckhts_simd_backend_canonical(const char *backend) {
    if (!backend || duckhts_backend_name_equal(backend, "auto")) return "auto";
    if (duckhts_backend_name_equal(backend, "scalar")) return "scalar";
    if (duckhts_backend_name_equal(backend, "avx2")) return "avx2";
    return NULL;
}

static const duckhts_simd_ops_t *duckhts_simd_best_ops(void) {
    const duckhts_simd_ops_t *ops = duckhts_simd_avx2_ops_if_available();
    return ops ? ops : duckhts_simd_scalar_ops();
}

static const duckhts_simd_ops_t *duckhts_simd_ops_for_backend(const char *backend) {
    if (!backend || duckhts_backend_name_equal(backend, "auto")) {
        return duckhts_simd_best_ops();
    }
    if (duckhts_backend_name_equal(backend, "scalar")) {
        return duckhts_simd_scalar_ops();
    }
    if (duckhts_backend_name_equal(backend, "avx2")) {
        return duckhts_simd_avx2_ops_if_available();
    }
    return NULL;
}

void duckhts_simd_init(void) {
    if (duckhts_simd_initialized) return;
    duckhts_simd_store_requested("auto");
    duckhts_simd_store_ops(duckhts_simd_best_ops());
    duckhts_simd_initialized = 1;
}

int duckhts_simd_backend_available(const char *backend) {
    if (!backend) return 0;
    if (duckhts_backend_name_equal(backend, "auto")) return 1;
    if (duckhts_backend_name_equal(backend, "scalar")) return 1;
    if (duckhts_backend_name_equal(backend, "avx2")) return duckhts_simd_avx2_available();
    return 0;
}

int duckhts_simd_set_backend(const char *backend, char *err, size_t err_len) {
    const duckhts_simd_ops_t *ops;
    const char *canonical;

    if (!duckhts_simd_initialized) duckhts_simd_init();

    canonical = duckhts_simd_backend_canonical(backend);
    if (!canonical) {
        if (err && err_len) {
            snprintf(err, err_len,
                     "unknown SIMD backend '%s' (expected auto, scalar, or a compiled SIMD backend)",
                     backend ? backend : "");
        }
        return 0;
    }

    ops = duckhts_simd_ops_for_backend(canonical);
    if (!ops) {
        if (err && err_len) {
            snprintf(err, err_len, "SIMD backend '%s' is not available in this process", canonical);
        }
        return 0;
    }

    duckhts_simd_store_requested(canonical);
    duckhts_simd_store_ops(ops);
    return 1;
}

const char *duckhts_simd_selected_backend(void) {
    if (!duckhts_simd_initialized) duckhts_simd_init();
    return duckhts_simd_load_ops()->name;
}

const char *duckhts_simd_requested_backend(void) {
    if (!duckhts_simd_initialized) duckhts_simd_init();
    return duckhts_simd_load_requested();
}

void duckhts_simd_base_counts(const char *seq, size_t len,
                              duckhts_simd_base_counts_t *out) {
    const duckhts_simd_ops_t *ops;
    if (!duckhts_simd_initialized) duckhts_simd_init();
    ops = duckhts_simd_load_ops();
    ops->base_counts(seq, len, out);
}

static inline void set_null_at(duckdb_vector vector, idx_t row) {
    duckdb_vector_ensure_validity_writable(vector);
    duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vector), row);
}

static inline int row_is_valid(duckdb_vector vector, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    if (!validity) return 1;
    return duckdb_validity_row_is_valid(validity, row);
}

static char ascii_lower(char c) {
    if (c >= 'A' && c <= 'Z') return (char)(c - 'A' + 'a');
    return c;
}

static int ascii_space(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
}

static int normalize_backend_name(const char *src, idx_t src_len, char *dst, size_t dst_len) {
    idx_t start = 0;
    idx_t end = src_len;
    size_t out_len;

    while (start < end && ascii_space(src[start])) start++;
    while (end > start && ascii_space(src[end - 1])) end--;

    out_len = (size_t)(end - start);
    if (out_len == 0 || out_len >= dst_len) return 0;

    for (size_t i = 0; i < out_len; i++) {
        dst[i] = ascii_lower(src[start + (idx_t)i]);
    }
    dst[out_len] = '\0';
    return 1;
}

static int backend_arg_at(duckdb_vector vector, idx_t row, char *dst, size_t dst_len) {
    duckdb_string_t *data = (duckdb_string_t *)duckdb_vector_get_data(vector);
    duckdb_string_t *value = &data[row];
    return normalize_backend_name(duckdb_string_t_data(value), duckdb_string_t_length(*value), dst, dst_len);
}

static void duckhts_simd_backend_scalar(duckdb_function_info info,
                                        duckdb_data_chunk input,
                                        duckdb_vector output) {
    (void)info;
    idx_t row_count = duckdb_data_chunk_get_size(input);
    const char *backend = duckhts_simd_selected_backend();
    idx_t len = (idx_t)strlen(backend);
    for (idx_t row = 0; row < row_count; row++) {
        duckdb_vector_assign_string_element_len(output, row, backend, len);
    }
}

static void duckhts_simd_requested_backend_scalar(duckdb_function_info info,
                                                  duckdb_data_chunk input,
                                                  duckdb_vector output) {
    (void)info;
    idx_t row_count = duckdb_data_chunk_get_size(input);
    const char *backend = duckhts_simd_requested_backend();
    idx_t len = (idx_t)strlen(backend);
    for (idx_t row = 0; row < row_count; row++) {
        duckdb_vector_assign_string_element_len(output, row, backend, len);
    }
}

static void duckhts_simd_backend_available_scalar(duckdb_function_info info,
                                                  duckdb_data_chunk input,
                                                  duckdb_vector output) {
    (void)info;
    duckdb_vector backend_vec = duckdb_data_chunk_get_vector(input, 0);
    bool *out_data = (bool *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < row_count; row++) {
        char backend[32];
        if (!row_is_valid(backend_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        if (!backend_arg_at(backend_vec, row, backend, sizeof(backend))) {
            out_data[row] = false;
            continue;
        }
        out_data[row] = duckhts_simd_backend_available(backend) != 0;
    }
}

static void duckhts_simd_set_backend_scalar(duckdb_function_info info,
                                            duckdb_data_chunk input,
                                            duckdb_vector output) {
    duckdb_vector backend_vec = duckdb_data_chunk_get_vector(input, 0);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < row_count; row++) {
        char backend[32];
        char err[160];
        const char *selected;
        idx_t selected_len;

        if (!row_is_valid(backend_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        if (!backend_arg_at(backend_vec, row, backend, sizeof(backend))) {
            duckdb_scalar_function_set_error(info, "duckhts_simd_set_backend: backend must be a non-empty short string");
            return;
        }
        if (!duckhts_simd_set_backend(backend, err, sizeof(err))) {
            char msg[208];
            snprintf(msg, sizeof(msg), "duckhts_simd_set_backend: %s", err);
            duckdb_scalar_function_set_error(info, msg);
            return;
        }

        selected = duckhts_simd_selected_backend();
        selected_len = (idx_t)strlen(selected);
        duckdb_vector_assign_string_element_len(output, row, selected, selected_len);
    }
}

static void register_simd_noarg_string_function(duckdb_connection connection,
                                                const char *name,
                                                duckdb_scalar_function_t fn_ptr) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_scalar_function_set_name(fn, name);
    duckdb_scalar_function_set_return_type(fn, varchar_type);
    duckdb_scalar_function_set_volatile(fn);
    duckdb_scalar_function_set_function(fn, fn_ptr);
    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_simd_backend_available_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_scalar_function_set_name(fn, "duckhts_simd_backend_available");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_function(fn, duckhts_simd_backend_available_scalar);
    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_simd_set_backend_function(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_scalar_function_set_name(fn, "duckhts_simd_set_backend");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, varchar_type);
    duckdb_scalar_function_set_volatile(fn);
    duckdb_scalar_function_set_function(fn, duckhts_simd_set_backend_scalar);
    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_scalar_function(&fn);
}

void register_duckhts_simd_functions(duckdb_connection connection) {
    register_simd_noarg_string_function(connection, "duckhts_simd_backend", duckhts_simd_backend_scalar);
    register_simd_noarg_string_function(connection, "duckhts_simd_requested_backend", duckhts_simd_requested_backend_scalar);
    register_simd_backend_available_function(connection);
    register_simd_set_backend_function(connection);
}
