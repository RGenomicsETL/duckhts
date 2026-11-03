/* Native regression for errors reported by callers holding the registration lock.
 * Build: cc -std=c11 -Wall -Wextra -Werror -fPIC -shared -Iduckdb_capi -Isrc/include \
 *   test/scripts/extension_init_oom_probe.c -ldl -o build/extension_init_oom_probe.so
 */
#define _GNU_SOURCE
#include "duckhts_registration.h"
#include <dlfcn.h>
#include <stdio.h>
#include <string.h>

#undef duckdb_open
#undef duckdb_close
#undef duckdb_connect
#undef duckdb_disconnect

static char error_message[512];

static void set_error(duckdb_extension_info info, const char *message) {
    (void)info;
    snprintf(error_message, sizeof(error_message), "%s", message);
}

int extension_init_oom_check(const char *path) {
    void *handle = dlopen(path, RTLD_NOW | RTLD_NOLOAD);
    if (!handle) return 1;
    duckdb_ext_api_v1 *api = dlsym(handle, "duckdb_ext_api");
    bool (*begin)(duckhts_registration_t *) = dlsym(handle, "duckhts_macro_registration_begin");
    bool (*error)(duckhts_registration_t *, const char *) = dlsym(handle, "duckhts_registration_error");
    bool (*end)(duckhts_registration_t *) = dlsym(handle, "duckhts_macro_registration_end");
    if (!api || !begin || !error || !end) return 2;
    duckdb_database database;
    duckdb_connection connection;
    if (api->duckdb_open(NULL, &database) != DuckDBSuccess) return 3;
    if (api->duckdb_connect(database, &connection) != DuckDBSuccess) return 4;
    struct duckdb_extension_access access = {.set_error = set_error};
    duckhts_registration_t registration = {.connection = connection, .access = &access};
    if (!begin(&registration)) return 5;
    /* The panel-macro allocation failure in duckhts.c takes this exact exit. */
    const char *message = "DuckHTS could not build the Somalier panel macro";
    if (error(&registration, message) || strcmp(error_message, message) != 0) return 6;
    fprintf(stderr, "Injected registration error returned; reacquiring registration lock\n");
    if (!begin(&registration)) return 7;
    /* No macros were registered in this probe; end reports that and releases. */
    if (end(&registration)) return 8;
    api->duckdb_disconnect(&connection);
    api->duckdb_close(&database);
    dlclose(handle);
    return 0;
}
