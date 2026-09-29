#include "duckhts_registration.h"
DUCKDB_EXTENSION_EXTERN

DUCKDB_EXTENSION_ENTRYPOINT(duckdb_connection connection, duckdb_extension_info info,
                           struct duckdb_extension_access *access) {
    duckhts_registration_t registration = {
        .connection = connection,
        .info = info,
        .access = access
    };
    if (!duckhts_macro_registration_begin(&registration)) {
        return false;
    }
    const char *const sql[] = {
        "CREATE OR REPLACE MACRO duckhts_init_failure() AS ",
        "(SELECT * FROM __duckhts_missing_init_relation)"
    };
    if (!duckhts_register_sql_parts(&registration, sql, sizeof(sql) / sizeof(sql[0]))) {
        return false;
    }
    return duckhts_macro_registration_end(&registration);
}
