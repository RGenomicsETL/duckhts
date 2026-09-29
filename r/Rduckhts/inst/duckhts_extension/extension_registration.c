#include "duckhts_registration.h"
#include "duckhts_public_macros.h"
DUCKDB_EXTENSION_EXTERN

#include <stdatomic.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define DUCKHTS_MAX_MACROS 128
#define MACRO_PREFIX "CREATE OR REPLACE MACRO "

typedef struct {
    char *name;
    char *sql;
    bool public;
} duckhts_macro_t;

static duckhts_macro_t macros[DUCKHTS_MAX_MACROS];
static size_t macro_count;
static char definitions_sha256[65];
static atomic_flag macro_lock = ATOMIC_FLAG_INIT;
/* The definition list is built under macro_lock by the first successful LOAD in the
   process and never changes afterwards. Readers use it only once it is published. */
static atomic_bool definitions_published = false;

bool duckhts_registration_error(duckhts_registration_t *registration, const char *message) {
    registration->access->set_error(registration->info, message);
    return false;
}

static bool macro_is_public(const char *name) {
    for (size_t i = 0; i < sizeof(duckhts_public_macros) / sizeof(duckhts_public_macros[0]); i++) {
        if (strcmp(name, duckhts_public_macros[i]) == 0) {
            return true;
        }
    }
    return false;
}

bool duckhts_macro_registration_begin(duckhts_registration_t *registration) {
    while (atomic_flag_test_and_set_explicit(&macro_lock, memory_order_acquire)) {
    }
    if (!atomic_load_explicit(&definitions_published, memory_order_acquire)) {
        for (size_t i = 0; i < macro_count; i++) {
            free(macros[i].name);
            free(macros[i].sql);
            macros[i].name = NULL;
            macros[i].sql = NULL;
        }
        macro_count = 0;
    }

    /* The decision concerns the database instance's default catalog, which is what the
       extension's initialization connection sees; a caller's USE of another catalog
       before LOAD does not change it. Those callers install the TEMP definitions. */
    duckdb_result result;
    duckdb_state state = duckdb_query(registration->connection,
        "SELECT count(*) = 1 AND coalesce(bool_and(path IS NULL OR path = ''), false) "
        "AND coalesce(bool_and(NOT readonly AND NOT internal AND type = 'duckdb'), false) "
        "FROM duckdb_databases() WHERE database_name = current_database()", &result);
    if (state != DuckDBSuccess) {
        const char *error = duckdb_result_error(&result);
        duckhts_registration_error(registration, error != NULL ? error : "Cannot inspect default database");
        duckdb_destroy_result(&result);
        atomic_flag_clear_explicit(&macro_lock, memory_order_release);
        return false;
    }
    duckdb_data_chunk chunk = duckdb_fetch_chunk(result);
    bool *values = duckdb_vector_get_data(duckdb_data_chunk_get_vector(chunk, 0));
    registration->install_persistent_macros = values[0];
    duckdb_destroy_data_chunk(&chunk);
    registration->macro_index = 0;
    duckdb_destroy_result(&result);
    return true;
}

static bool macro_digest(duckhts_registration_t *registration) {
    size_t length = 0;
    for (size_t i = 0; i < macro_count; i++) {
        size_t sql_length = strlen(macros[i].sql);
        if (length > SIZE_MAX - 2 || sql_length > SIZE_MAX - length - 2) {
            return duckhts_registration_error(registration, "Macro definitions are too large");
        }
        length += sql_length + 2;
    }
    if (length > UINT32_MAX) {
        return duckhts_registration_error(registration, "Macro definitions are too large");
    }
    char *joined = malloc(length);
    if (joined == NULL) {
        return duckhts_registration_error(registration, "Cannot allocate macro digest input");
    }
    size_t offset = 0;
    for (size_t i = 0; i < macro_count; i++) {
        size_t sql_length = strlen(macros[i].sql);
        joined[offset++] = macros[i].public ? 1 : 0;
        memcpy(joined + offset, macros[i].sql, sql_length);
        joined[offset + sql_length] = '\0';
        offset += sql_length + 1;
    }
    duckdb_prepared_statement prepared = NULL;
    duckdb_result result;
    duckdb_state state = duckdb_prepare(registration->connection, "SELECT sha256(?)", &prepared);
    if (state == DuckDBSuccess) {
        state = duckdb_bind_varchar_length(prepared, 1, joined, length);
    }
    if (state == DuckDBSuccess) {
        state = duckdb_execute_prepared(prepared, &result);
        if (state == DuckDBSuccess) {
            duckdb_data_chunk chunk = duckdb_fetch_chunk(result);
            duckdb_string_t *values = duckdb_vector_get_data(duckdb_data_chunk_get_vector(chunk, 0));
            if (duckdb_string_t_length(values[0]) != 64) {
                state = DuckDBError;
            } else {
                memcpy(definitions_sha256, duckdb_string_t_data(&values[0]), 64);
                definitions_sha256[64] = '\0';
            }
            duckdb_destroy_data_chunk(&chunk);
            duckdb_destroy_result(&result);
        }
    }
    duckdb_destroy_prepare(&prepared);
    free(joined);
    if (state != DuckDBSuccess) {
        return duckhts_registration_error(registration, "Cannot hash macro definitions");
    }
    return true;
}

bool duckhts_macro_registration_end(duckhts_registration_t *registration) {
    bool ok = registration->macro_index == macro_count;
    if (!ok) {
        duckhts_registration_error(registration, "Inconsistent DuckHTS macro registration order");
    } else if (!atomic_load_explicit(&definitions_published, memory_order_acquire)) {
        ok = macro_digest(registration);
        if (ok) {
            atomic_store_explicit(&definitions_published, true, memory_order_release);
        }
    }
    atomic_flag_clear_explicit(&macro_lock, memory_order_release);
    return ok;
}

bool duckhts_register_sql(duckhts_registration_t *registration, const char *sql) {
    const size_t prefix_length = sizeof(MACRO_PREFIX) - 1;
    if (strncmp(sql, MACRO_PREFIX, prefix_length) != 0) {
        duckhts_registration_error(registration, "Expected CREATE OR REPLACE MACRO statement");
        atomic_flag_clear_explicit(&macro_lock, memory_order_release);
        return false;
    }
    const char *name_end = strchr(sql + prefix_length, '(');
    if (name_end == NULL || name_end == sql + prefix_length ||
        registration->macro_index >= DUCKHTS_MAX_MACROS) {
        duckhts_registration_error(registration, "Invalid DuckHTS macro name or count");
        atomic_flag_clear_explicit(&macro_lock, memory_order_release);
        return false;
    }
    size_t name_length = (size_t)(name_end - sql - prefix_length);
    size_t index = registration->macro_index++;
    if (index == macro_count) {
        macros[index].name = malloc(name_length + 1);
        macros[index].sql = malloc(strlen(sql) + 1);
        if (macros[index].name == NULL || macros[index].sql == NULL) {
            free(macros[index].name);
            free(macros[index].sql);
            macros[index].name = NULL;
            macros[index].sql = NULL;
            duckhts_registration_error(registration, "Cannot allocate macro definition");
            atomic_flag_clear_explicit(&macro_lock, memory_order_release);
            return false;
        }
        memcpy(macros[index].name, sql + prefix_length, name_length);
        macros[index].name[name_length] = '\0';
        strcpy(macros[index].sql, sql);
        macros[index].public = macro_is_public(macros[index].name);
        macro_count++;
    } else if (strcmp(macros[index].sql, sql) != 0) {
        duckhts_registration_error(registration, "Inconsistent DuckHTS macro definition");
        atomic_flag_clear_explicit(&macro_lock, memory_order_release);
        return false;
    }

    if (!registration->install_persistent_macros) {
        return true;
    }
    duckdb_result result;
    duckdb_state state = duckdb_query(registration->connection, sql, &result);
    if (state != DuckDBSuccess) {
        const char *error = duckdb_result_error(&result);
        duckhts_registration_error(registration,
            error != NULL && *error != '\0' ? error : "DuckHTS SQL registration failed");
    }
    duckdb_destroy_result(&result);
    if (state != DuckDBSuccess) {
        atomic_flag_clear_explicit(&macro_lock, memory_order_release);
    }
    return state == DuckDBSuccess;
}

bool duckhts_register_sql_parts(duckhts_registration_t *registration,
                               const char *const *parts, size_t count) {
    size_t length = 0;
    for (size_t i = 0; i < count; i++) {
        size_t part_length = strlen(parts[i]);
        if (part_length > SIZE_MAX - length - 1) {
            duckhts_registration_error(registration, "DuckHTS SQL registration text is too large");
            atomic_flag_clear_explicit(&macro_lock, memory_order_release);
            return false;
        }
        length += part_length;
    }

    char *sql = malloc(length + 1);
    if (sql == NULL) {
        duckhts_registration_error(registration, "DuckHTS could not allocate SQL registration text");
        atomic_flag_clear_explicit(&macro_lock, memory_order_release);
        return false;
    }
    size_t offset = 0;
    for (size_t i = 0; i < count; i++) {
        size_t part_length = strlen(parts[i]);
        memcpy(sql + offset, parts[i], part_length);
        offset += part_length;
    }
    sql[offset] = '\0';
    bool ok = duckhts_register_sql(registration, sql);
    free(sql);
    return ok;
}

typedef struct {
    size_t offset;
} macro_scan_t;

static void macro_definitions_bind(duckdb_bind_info info) {
    const char *names[] = {"install_order", "name", "public", "sql", "definitions_sha256"};
    duckdb_type types[] = {DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_BOOLEAN,
                           DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR};
    for (size_t i = 0; i < 5; i++) {
        duckdb_logical_type type = duckdb_create_logical_type(types[i]);
        duckdb_bind_add_result_column(info, names[i], type);
        duckdb_destroy_logical_type(&type);
    }
}

static void macro_scan_destroy(void *ptr) {
    duckdb_free(ptr);
}

static void macro_definitions_init(duckdb_init_info info) {
    if (!atomic_load_explicit(&definitions_published, memory_order_acquire)) {
        duckdb_init_set_error(info, "DuckHTS macro definitions are not available until DuckHTS has loaded");
        return;
    }
    macro_scan_t *scan = duckdb_malloc(sizeof(*scan));
    if (scan == NULL) {
        duckdb_init_set_error(info, "Cannot allocate macro scan");
        return;
    }
    scan->offset = 0;
    duckdb_init_set_init_data(info, scan, macro_scan_destroy);
}

static void macro_definitions_scan(duckdb_function_info info, duckdb_data_chunk output) {
    macro_scan_t *scan = duckdb_function_get_init_data(info);
    idx_t capacity = duckdb_vector_size();
    idx_t count = 0;
    uint32_t *orders = duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 0));
    bool *public_flags = duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 2));
    while (scan->offset < macro_count && count < capacity) {
        duckhts_macro_t *macro = &macros[scan->offset];
        orders[count] = (uint32_t)scan->offset + 1;
        public_flags[count] = macro->public;
        duckdb_vector_assign_string_element(duckdb_data_chunk_get_vector(output, 1), count, macro->name);
        const char *sql = macro->sql;
        size_t prefix = sizeof("CREATE OR REPLACE ") - 1;
        size_t sql_length = strlen(sql);
        char *temporary = malloc(sql_length + sizeof("TEMP "));
        if (temporary == NULL) {
            duckdb_function_set_error(info, "Cannot allocate TEMP macro statement");
            return;
        }
        memcpy(temporary, sql, prefix);
        memcpy(temporary + prefix, "TEMP ", sizeof("TEMP ") - 1);
        memcpy(temporary + prefix + sizeof("TEMP ") - 1, sql + prefix, sql_length - prefix + 1);
        duckdb_vector_assign_string_element(duckdb_data_chunk_get_vector(output, 3), count, temporary);
        free(temporary);
        duckdb_vector_assign_string_element(duckdb_data_chunk_get_vector(output, 4), count, definitions_sha256);
        count++;
        scan->offset++;
    }
    duckdb_data_chunk_set_size(output, count);
}

bool duckhts_register_macro_definitions(duckhts_registration_t *registration) {
    duckdb_table_function function = duckdb_create_table_function();
    duckdb_table_function_set_name(function, "duckhts_macro_definitions");
    duckdb_table_function_set_bind(function, macro_definitions_bind);
    duckdb_table_function_set_init(function, macro_definitions_init);
    duckdb_table_function_set_function(function, macro_definitions_scan);
    duckdb_state state = duckdb_register_table_function(registration->connection, function);
    duckdb_destroy_table_function(&function);
    if (state != DuckDBSuccess) {
        return duckhts_registration_error(registration, "Cannot register macro definitions table function");
    }
    return true;
}
