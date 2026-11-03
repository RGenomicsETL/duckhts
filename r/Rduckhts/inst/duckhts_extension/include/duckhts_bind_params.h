#ifndef DUCKHTS_BIND_PARAMS_H
#define DUCKHTS_BIND_PARAMS_H

#include "duckdb_extension.h"

#include <stdbool.h>
#include <stdint.h>
#include <string.h>

/* Named-parameter readers for table-function bind. An absent or NULL parameter
 * gives the fallback. Include after DUCKDB_EXTENSION_EXTERN: the calls resolve
 * through the extension's C API table. */

static inline duckdb_value duckhts_bind_value(duckdb_bind_info info, const char *name) {
    duckdb_value value = duckdb_bind_get_named_parameter(info, name);
    if (value && duckdb_is_null_value(value)) duckdb_destroy_value(&value);
    return value;
}

static inline bool duckhts_bind_bool(duckdb_bind_info info, const char *name, bool fallback) {
    duckdb_value value = duckhts_bind_value(info, name);
    bool result = value ? duckdb_get_bool(value) : fallback;
    if (value) duckdb_destroy_value(&value);
    return result;
}

static inline int64_t duckhts_bind_int64(duckdb_bind_info info, const char *name, int64_t fallback) {
    duckdb_value value = duckhts_bind_value(info, name);
    int64_t result = value ? duckdb_get_int64(value) : fallback;
    if (value) duckdb_destroy_value(&value);
    return result;
}

/* Sets *out to an owned copy (duckdb_free), or NULL when absent or NULL. Returns
 * false only when a present value could not be copied. */
static inline bool duckhts_bind_varchar(duckdb_bind_info info, const char *name, char **out) {
    duckdb_value value = duckhts_bind_value(info, name);
    *out = value ? duckdb_get_varchar(value) : NULL;
    bool ok = !value || *out;
    if (value) duckdb_destroy_value(&value);
    return ok;
}

/* Index of the value in `choices` (NULL-terminated), `fallback` when absent, -1 when
 * the value is not one of the choices, -2 when it could not be read. */
static inline int duckhts_bind_choice(duckdb_bind_info info, const char *name,
                                      const char *const *choices, int fallback) {
    char *text = NULL;
    if (!duckhts_bind_varchar(info, name, &text)) return -2;
    int result = text ? -1 : fallback;
    for (int i = 0; text && choices[i]; i++) {
        if (strcmp(text, choices[i]) == 0) result = i;
    }
    if (text) duckdb_free(text);
    return result;
}

#endif
