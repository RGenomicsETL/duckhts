/* Test-only interposition on the loaded extension's copied DuckDB API table.
 * Never replace the host allocator or expose a production fault-injection API.
 * The driver must run one query at a time with SET threads=1. */
#define _GNU_SOURCE
#include "duckdb_extension.h"
#include <dlfcn.h>
#include <stdlib.h>

static duckdb_ext_api_v1 duckdb_ext_api;
#include "duckdb_alloc.h"
#undef duckdb_malloc
#undef duckdb_free

static duckdb_ext_api_v1 *extension_api;
static void *extension_handle;
static void *(*saved_malloc)(size_t);
static void (*saved_free)(void *);
static void *live[16384];
static long live_count, attempts, fail_at, failures;

static void *probe_malloc(size_t bytes) {
    attempts++;
    if (attempts == fail_at) {
        failures++;
        return NULL;
    }
    void *data = saved_malloc(bytes);
    if (!data) abort(); /* An uncontrolled OOM is not the injected failure. */
    for (size_t i = 0; i < sizeof(live) / sizeof(*live); i++) {
        if (!live[i]) {
            live[i] = data;
            live_count++;
            return data;
        }
    }
    abort(); /* Test ledger capacity is not a reader resource limit. */
}

static void probe_free(void *data) {
    if (data) {
        for (size_t i = 0; i < sizeof(live) / sizeof(*live); i++) {
            if (live[i] == data) {
                live[i] = NULL;
                live_count--;
                break;
            }
        }
    }
    /* Also release host-created values (e.g. duckdb_get_varchar results),
     * which are deliberately outside the first-party allocation ledger. */
    saved_free(data);
}

int reader_alloc_open(const char *path) {
    extension_handle = dlopen(path, RTLD_NOW | RTLD_NOLOAD);
    if (!extension_handle) return 1;
    extension_api = dlsym(extension_handle, "duckdb_ext_api");
    if (!extension_api) return 2;
    saved_malloc = extension_api->duckdb_malloc;
    saved_free = extension_api->duckdb_free;
    return !saved_malloc || !saved_free;
}

int reader_alloc_arm(long nth) {
    if (!extension_api || live_count) return 1;
    attempts = failures = 0;
    fail_at = nth;
    extension_api->duckdb_malloc = probe_malloc;
    extension_api->duckdb_free = probe_free;
    return 0;
}

long reader_alloc_attempts(void) { return attempts; }
long reader_alloc_failures(void) { return failures; }
long reader_alloc_live(void) { return live_count; }

void reader_alloc_disarm(void) {
    extension_api->duckdb_malloc = saved_malloc;
    extension_api->duckdb_free = saved_free;
}

void reader_alloc_close(void) {
    reader_alloc_disarm();
    dlclose(extension_handle);
    extension_handle = NULL;
    extension_api = NULL;
}

int reader_alloc_helper_checks(void) {
    duckdb_ext_api.duckdb_malloc = probe_malloc;
    long before = attempts;
    if (duckhts_alloc_array(UINT64_MAX, 2) || duckhts_alloc_array(0, 1) ||
        duckhts_alloc_array(1, 0) || attempts != before) return 1;
    unsigned char *bytes = duckhts_alloc_array(3, 7);
    if (!bytes) return 2;
    for (int i = 0; i < 21; i++) {
        if (bytes[i]) return 3;
    }
    probe_free(bytes);
    char *text = duckhts_copy_string("sample");
    if (!text || strcmp(text, "sample")) return 4;
    probe_free(text);
    fail_at = attempts + 1;
    if (duckhts_alloc_array(1, 7)) return 5;
    fail_at = attempts + 1;
    if (duckhts_copy_string("sample")) return 6;
    return live_count != 0;
}

/* R's .C passes addresses; these adapters exercise the installed package's
 * separately built extension without depending on R headers or changing it. */
void reader_alloc_r_open(char **path, int *status) {
    *status = reader_alloc_open(*path);
}

void reader_alloc_r_arm(int *nth, int *status) {
    *status = reader_alloc_arm(*nth);
}

void reader_alloc_r_stats(int *count, int *remaining, int *failed) {
    *count = (int)attempts;
    *remaining = (int)live_count;
    *failed = (int)failures;
}
