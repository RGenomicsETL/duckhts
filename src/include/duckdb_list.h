#ifndef DUCKHTS_DUCKDB_LIST_H
#define DUCKHTS_DUCKDB_LIST_H

#include <stdint.h>

/* Include after duckdb_extension.h / DUCKDB_EXTENSION_EXTERN. DuckDB owns
 * child storage. Do not retain child data pointers across this operation or
 * publish the row entry until it succeeds. Failure leaves entry untouched. */
static inline int duckhts_list_extend(duckdb_vector vec, idx_t count,
                                      duckdb_list_entry *entry) {
    idx_t offset = duckdb_list_vector_get_size(vec);
    if (count > UINT64_MAX - offset) return 0;
    if (count && (duckdb_list_vector_reserve(vec, offset + count) != DuckDBSuccess ||
                  duckdb_list_vector_set_size(vec, offset + count) != DuckDBSuccess)) return 0;
    entry->offset = offset;
    entry->length = count;
    return 1;
}

#endif
