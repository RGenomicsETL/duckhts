#ifndef DUCKHTS_DUCKDB_ALLOC_H
#define DUCKHTS_DUCKDB_ALLOC_H

#include <stdint.h>
#include <string.h>

/* Include after duckdb_extension.h and DUCKDB_EXTENSION_EXTERN. Storage is
 * owned by the caller and released with duckdb_free. Empty arrays return NULL.
 * Zero each allocation immediately so partial owners are safe to destroy. */
static inline void *duckhts_alloc_array(idx_t count, size_t width) {
    if (!count || !width || count > SIZE_MAX / width) return NULL;
    size_t bytes = (size_t)count * width;
    void *data = duckdb_malloc(bytes);
    if (data) memset(data, 0, bytes);
    return data;
}

static inline char *duckhts_copy_string(const char *text) {
    if (!text) return NULL;
    size_t length = strlen(text);
    if (length == SIZE_MAX) return NULL;
    char *copy = duckdb_malloc(length + 1);
    if (copy) memcpy(copy, text, length + 1);
    return copy;
}

#endif
