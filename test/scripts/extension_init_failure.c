#include "duckhts_registration.h"
DUCKDB_EXTENSION_EXTERN

DUCKDB_EXTENSION_ENTRYPOINT(duckdb_connection connection, duckdb_extension_info info,
                           struct duckdb_extension_access *access) {
    duckhts_registration_t registration = {
        .connection = connection,
        .info = info,
        .access = access
    };
    const char *const sql[] = {
        "SELECT * FROM ",
        "__duckhts_missing_init_relation"
    };
    return duckhts_register_sql_parts(&registration, sql, sizeof(sql) / sizeof(sql[0]));
}
