#ifndef DUCKHTS_REGISTRATION_H
#define DUCKHTS_REGISTRATION_H

#include "duckdb_extension.h"

#include <stdbool.h>
#include <stddef.h>

/* Borrowed handles valid only during the extension initialization call. */
typedef struct {
    duckdb_connection connection;
    duckdb_extension_info info;
    struct duckdb_extension_access *access;
    size_t macro_index;
    bool install_persistent_macros;
} duckhts_registration_t;

bool duckhts_registration_error(duckhts_registration_t *registration, const char *message);
bool duckhts_macro_registration_begin(duckhts_registration_t *registration);
bool duckhts_macro_registration_end(duckhts_registration_t *registration);
bool duckhts_register_macro_definitions(duckhts_registration_t *registration);
bool duckhts_register_sql(duckhts_registration_t *registration, const char *sql);
bool duckhts_register_sql_parts(duckhts_registration_t *registration,
                               const char *const *parts, size_t count);

#endif
