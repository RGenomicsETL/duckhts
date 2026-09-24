#include "duckhts_registration.h"
DUCKDB_EXTENSION_EXTERN

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

bool duckhts_registration_error(duckhts_registration_t *registration, const char *message) {
    registration->access->set_error(registration->info, message);
    return false;
}

bool duckhts_register_sql(duckhts_registration_t *registration, const char *sql) {
    duckdb_result result;
    duckdb_state state = duckdb_query(registration->connection, sql, &result);

    if (state != DuckDBSuccess) {
        const char *error = duckdb_result_error(&result);
        duckhts_registration_error(registration,
            error != NULL && *error != '\0' ? error : "DuckHTS SQL registration failed");
    }
    duckdb_destroy_result(&result);
    return state == DuckDBSuccess;
}

bool duckhts_register_sql_parts(duckhts_registration_t *registration,
                               const char *const *parts, size_t count) {
    size_t length = 0;
    for (size_t i = 0; i < count; i++) {
        size_t part_length = strlen(parts[i]);
        if (part_length > SIZE_MAX - length - 1) {
            return duckhts_registration_error(registration,
                "DuckHTS SQL registration text is too large");
        }
        length += part_length;
    }

    char *sql = malloc(length + 1);
    if (sql == NULL) {
        return duckhts_registration_error(registration,
            "DuckHTS could not allocate SQL registration text");
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
