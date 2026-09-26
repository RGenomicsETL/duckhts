#ifndef DUCKHTS_NAMED_ATTRIBUTE_COLUMNS_H
#define DUCKHTS_NAMED_ATTRIBUTE_COLUMNS_H

#include "duckdb.h"
#include <stddef.h>

/* Bind-owned keys; projected entries borrow keys from the bind until scan completion. */
typedef struct {
    char *name;
    size_t length;
} duckhts_attribute_key;

typedef struct {
    duckhts_attribute_key *keys;
    idx_t count;
    idx_t first_column;
} duckhts_attribute_columns;

typedef struct {
    const duckhts_attribute_key *key;
    idx_t output_column;
    const char *value;
    size_t value_length;
} duckhts_projected_attribute;

/* reserved_names are the reader's output names, including optional columns. */
bool duckhts_attribute_columns_bind(duckdb_bind_info info, const char *parameter,
                                    const char *const *reserved_names, size_t reserved_count,
                                    duckhts_attribute_columns *columns);
void duckhts_attribute_columns_declare(duckdb_bind_info info,
                                       duckhts_attribute_columns *columns, idx_t first_column);
void duckhts_attribute_columns_destroy(duckhts_attribute_columns *columns);

/* Map projected DuckDB ids to keys. Returns NULL with count=0 when none are projected. */
duckhts_projected_attribute *duckhts_attribute_columns_project(
    const duckhts_attribute_columns *columns, const idx_t *column_ids,
    idx_t column_count, idx_t *projected_count);

/* Reader supplies borrowed value spans (NULL for absent keys) for each row. */
void duckhts_attribute_columns_write(const duckhts_projected_attribute *projected,
                                     idx_t projected_count, duckdb_vector *vectors, idx_t row);

#endif
