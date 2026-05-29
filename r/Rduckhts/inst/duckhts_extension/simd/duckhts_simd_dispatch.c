#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <config.h>

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

static const char *duckhts_simd_dispatch_mode_name = "single-shlib-ops-table";

typedef int (*duckhts_simd_status_fn)(void);
typedef const duckhts_simd_ops_t *(*duckhts_simd_ops_provider_fn)(void);

typedef struct duckhts_simd_backend_entry {
    const char *name;
    int selectable;
    duckhts_simd_status_fn compiled;
    duckhts_simd_status_fn cpu_supported;
    duckhts_simd_ops_provider_fn ops_if_available;
    int priority;
} duckhts_simd_backend_entry_t;

static int duckhts_simd_true(void) { return 1; }
static int duckhts_simd_false(void) { return 0; }
static const duckhts_simd_ops_t *duckhts_simd_no_ops(void) { return (const duckhts_simd_ops_t *)0; }

static int duckhts_backend_name_equal(const char *a, const char *b) {
    return a && b && strcmp(a, b) == 0;
}

static int duckhts_cpu_has_sse2(void) {
#if (defined(__x86_64__) || defined(__i386__)) && \
    defined(HAVE_BUILTIN_CPU_SUPPORT_SSSE3) && HAVE_BUILTIN_CPU_SUPPORT_SSSE3
    __builtin_cpu_init();
    return __builtin_cpu_supports("sse2") != 0;
#else
    return 0;
#endif
}

static int duckhts_cpu_has_sse41(void) {
#if (defined(__x86_64__) || defined(__i386__)) && \
    defined(HAVE_BUILTIN_CPU_SUPPORT_SSSE3) && HAVE_BUILTIN_CPU_SUPPORT_SSSE3
    __builtin_cpu_init();
    return __builtin_cpu_supports("sse4.1") != 0;
#else
    return 0;
#endif
}

static const duckhts_simd_backend_entry_t duckhts_simd_backend_entries[] = {
    {"scalar",       1, duckhts_simd_true,                  duckhts_simd_true,                         duckhts_simd_scalar_ops,                     70},
    {"sse2",         0, duckhts_simd_false,                 duckhts_cpu_has_sse2,                      duckhts_simd_no_ops,                         40},
    {"sse41",        0, duckhts_simd_false,                 duckhts_cpu_has_sse41,                     duckhts_simd_no_ops,                         30},
    {"avx2",         1, duckhts_simd_avx2_compiled,         duckhts_simd_avx2_cpu_supported,          duckhts_simd_avx2_ops_if_available,          20},
    {"avx512",       1, duckhts_simd_avx512_compiled,       duckhts_simd_avx512_cpu_supported,        duckhts_simd_avx512_ops_if_available,        10},
    {"neon",         1, duckhts_simd_neon_compiled,         duckhts_simd_neon_cpu_supported,          duckhts_simd_neon_ops_if_available,          50},
    {"wasm_simd128", 1, duckhts_simd_wasm_simd128_compiled, duckhts_simd_wasm_simd128_cpu_supported, duckhts_simd_wasm_simd128_ops_if_available, 60}
};

static size_t duckhts_simd_backend_count(void) {
    return sizeof(duckhts_simd_backend_entries) / sizeof(duckhts_simd_backend_entries[0]);
}

static const duckhts_simd_backend_entry_t *duckhts_simd_find_backend(const char *backend) {
    if (!backend) return (const duckhts_simd_backend_entry_t *)0;
    for (size_t i = 0; i < duckhts_simd_backend_count(); i++) {
        if (duckhts_backend_name_equal(backend, duckhts_simd_backend_entries[i].name)) {
            return &duckhts_simd_backend_entries[i];
        }
    }
    return (const duckhts_simd_backend_entry_t *)0;
}

static const char *duckhts_simd_backend_canonical(const char *backend) {
    const duckhts_simd_backend_entry_t *entry;
    if (!backend || duckhts_backend_name_equal(backend, "auto")) return "auto";
    entry = duckhts_simd_find_backend(backend);
    return entry ? entry->name : (const char *)0;
}

static int duckhts_simd_entry_compiled(const duckhts_simd_backend_entry_t *entry) {
    return entry && entry->compiled && entry->compiled() != 0;
}

static int duckhts_simd_entry_cpu_supported(const duckhts_simd_backend_entry_t *entry) {
    return entry && entry->cpu_supported && entry->cpu_supported() != 0;
}

static int duckhts_simd_entry_available(const duckhts_simd_backend_entry_t *entry) {
    return duckhts_simd_entry_compiled(entry) && duckhts_simd_entry_cpu_supported(entry);
}

static const duckhts_simd_ops_t *duckhts_simd_best_ops(void) {
    const duckhts_simd_backend_entry_t *best = (const duckhts_simd_backend_entry_t *)0;
    const duckhts_simd_ops_t *best_ops = (const duckhts_simd_ops_t *)0;

    for (size_t i = 0; i < duckhts_simd_backend_count(); i++) {
        const duckhts_simd_backend_entry_t *entry = &duckhts_simd_backend_entries[i];
        const duckhts_simd_ops_t *ops;
        if (!entry->selectable || !duckhts_simd_entry_available(entry)) continue;
        ops = entry->ops_if_available ? entry->ops_if_available() : (const duckhts_simd_ops_t *)0;
        if (!ops) continue;
        if (!best || entry->priority < best->priority) {
            best = entry;
            best_ops = ops;
        }
    }

    return best_ops ? best_ops : duckhts_simd_scalar_ops();
}

static const duckhts_simd_ops_t *duckhts_simd_ops_for_backend(const char *backend) {
    const duckhts_simd_backend_entry_t *entry;
    if (!backend || duckhts_backend_name_equal(backend, "auto")) {
        return duckhts_simd_best_ops();
    }
    entry = duckhts_simd_find_backend(backend);
    if (!entry || !entry->selectable || !duckhts_simd_entry_available(entry)) {
        return (const duckhts_simd_ops_t *)0;
    }
    return entry->ops_if_available ? entry->ops_if_available() : (const duckhts_simd_ops_t *)0;
}

void duckhts_simd_init(void) {
    if (duckhts_simd_initialized) return;
    duckhts_simd_store_requested("auto");
    duckhts_simd_store_ops(duckhts_simd_best_ops());
    duckhts_simd_initialized = 1;
}

int duckhts_simd_backend_compiled(const char *backend) {
    return duckhts_simd_entry_compiled(duckhts_simd_find_backend(backend));
}

int duckhts_simd_backend_cpu_supported(const char *backend) {
    return duckhts_simd_entry_cpu_supported(duckhts_simd_find_backend(backend));
}

int duckhts_simd_backend_available(const char *backend) {
    return duckhts_simd_entry_available(duckhts_simd_find_backend(backend));
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

static void duckhts_simd_backend_status_scalar(duckdb_function_info info,
                                                duckdb_data_chunk input,
                                                duckdb_vector output,
                                                int (*status_fn)(const char *)) {
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
        out_data[row] = status_fn(backend) != 0;
    }
}

static void duckhts_simd_backend_compiled_scalar(duckdb_function_info info,
                                                 duckdb_data_chunk input,
                                                 duckdb_vector output) {
    duckhts_simd_backend_status_scalar(info, input, output, duckhts_simd_backend_compiled);
}

static void duckhts_simd_backend_cpu_supported_scalar(duckdb_function_info info,
                                                      duckdb_data_chunk input,
                                                      duckdb_vector output) {
    duckhts_simd_backend_status_scalar(info, input, output, duckhts_simd_backend_cpu_supported);
}

static void duckhts_simd_backend_available_scalar(duckdb_function_info info,
                                                  duckdb_data_chunk input,
                                                  duckdb_vector output) {
    duckhts_simd_backend_status_scalar(info, input, output, duckhts_simd_backend_available);
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

static void register_simd_bool_backend_function(duckdb_connection connection,
                                                const char *name,
                                                duckdb_scalar_function_t fn_ptr) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_scalar_function_set_name(fn, name);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_function(fn, fn_ptr);
    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_scalar_function(&fn);
}

typedef struct duckhts_simd_info_init {
    idx_t offset;
} duckhts_simd_info_init_t;

static void destroy_simd_info_init(void *ptr) {
    duckdb_free(ptr);
}

static void duckhts_simd_info_bind(duckdb_bind_info info) {
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);

    duckdb_bind_add_result_column(info, "backend", varchar_type);
    duckdb_bind_add_result_column(info, "compiled", bool_type);
    duckdb_bind_add_result_column(info, "cpu_supported", bool_type);
    duckdb_bind_add_result_column(info, "available", bool_type);
    duckdb_bind_add_result_column(info, "selected", bool_type);
    duckdb_bind_add_result_column(info, "requested", bool_type);
    duckdb_bind_add_result_column(info, "dispatch_mode", varchar_type);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);
}

static void duckhts_simd_info_init(duckdb_init_info info) {
    duckhts_simd_info_init_t *init = (duckhts_simd_info_init_t *)duckdb_malloc(sizeof(*init));
    if (!init) {
        duckdb_init_set_error(info, "duckhts_simd_info: failed to allocate init state");
        return;
    }
    init->offset = 0;
    duckdb_init_set_init_data(info, init, destroy_simd_info_init);
}

static void duckhts_simd_info_scan(duckdb_function_info info, duckdb_data_chunk output) {
    duckhts_simd_info_init_t *init = (duckhts_simd_info_init_t *)duckdb_function_get_init_data(info);
    idx_t vector_size = duckdb_vector_size();
    idx_t row_count = 0;
    idx_t n_backends = (idx_t)duckhts_simd_backend_count();
    const char *selected = duckhts_simd_selected_backend();
    const char *requested = duckhts_simd_requested_backend();

    duckdb_vector backend_vec = duckdb_data_chunk_get_vector(output, 0);
    bool *compiled = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 1));
    bool *cpu_supported = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 2));
    bool *available = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 3));
    bool *selected_out = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 4));
    bool *requested_out = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 5));
    duckdb_vector dispatch_mode_vec = duckdb_data_chunk_get_vector(output, 6);

    if (!init) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    while (row_count < vector_size && init->offset < n_backends) {
        const char *backend = duckhts_simd_backend_entries[init->offset].name;
        int is_compiled = duckhts_simd_backend_compiled(backend) != 0;
        int is_cpu_supported = duckhts_simd_backend_cpu_supported(backend) != 0;

        duckdb_vector_assign_string_element(backend_vec, row_count, backend);
        compiled[row_count] = is_compiled;
        cpu_supported[row_count] = is_cpu_supported;
        available[row_count] = (is_compiled && is_cpu_supported) != 0;
        selected_out[row_count] = duckhts_backend_name_equal(backend, selected) != 0;
        requested_out[row_count] = duckhts_backend_name_equal(backend, requested) != 0;
        duckdb_vector_assign_string_element(dispatch_mode_vec, row_count, duckhts_simd_dispatch_mode_name);

        row_count++;
        init->offset++;
    }

    duckdb_data_chunk_set_size(output, row_count);
}

static void register_simd_info_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "duckhts_simd_info");
    duckdb_table_function_set_bind(tf, duckhts_simd_info_bind);
    duckdb_table_function_set_init(tf, duckhts_simd_info_init);
    duckdb_table_function_set_function(tf, duckhts_simd_info_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
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
    register_simd_bool_backend_function(connection, "duckhts_simd_backend_compiled", duckhts_simd_backend_compiled_scalar);
    register_simd_bool_backend_function(connection, "duckhts_simd_backend_cpu_supported", duckhts_simd_backend_cpu_supported_scalar);
    register_simd_bool_backend_function(connection, "duckhts_simd_backend_available", duckhts_simd_backend_available_scalar);
    register_simd_info_function(connection);
    register_simd_set_backend_function(connection);
}
