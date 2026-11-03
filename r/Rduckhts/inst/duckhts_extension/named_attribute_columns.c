#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "named_attribute_columns.h"
#include "duckhts_column_names.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void duckhts_attribute_columns_destroy(duckhts_attribute_columns *columns) {
    for (idx_t i = 0; i < columns->count; i++) duckdb_free(columns->keys[i].name);
    free(columns->keys);
    memset(columns, 0, sizeof(*columns));
}

bool duckhts_attribute_columns_bind(duckdb_bind_info info, const char *parameter,
                                    const char *const *reserved_names, size_t reserved_count,
                                    duckhts_attribute_columns *columns) {
    duckdb_value list = duckdb_bind_get_named_parameter(info, parameter);
    if (!list) return true;
    const char *error = NULL;
    char detail[512];
    char message[600];
    duckhts_column_names names;
    duckhts_column_names_init(&names, 0);
    for (size_t j = 0; j < reserved_count; j++) {
        if (duckhts_column_names_claim(&names, reserved_names[j], NULL, detail, sizeof(detail)) != 0) {
            error = detail;
            goto done;
        }
    }
    if (duckdb_is_null_value(list)) {
        error = "attributes must be a non-NULL VARCHAR[]";
        goto done;
    }
    idx_t count = duckdb_get_list_size(list);
    if (count > 256) {
        error = "attributes supports at most 256 keys";
        goto done;
    }
    if (count > 0) {
        columns->keys = calloc((size_t)count, sizeof(*columns->keys));
        if (!columns->keys) {
            error = "Out of memory allocating attribute keys";
            goto done;
        }
    }
    for (idx_t i = 0; i < count; i++) {
        duckdb_value value = duckdb_get_list_child(list, i);
        char *name = duckdb_is_null_value(value) ? NULL : duckdb_get_varchar(value);
        duckdb_destroy_value(&value);
        columns->keys[i].name = name;
        columns->count++;
        if (!name || !name[0] || strlen(name) > 1024) {
            error = "attributes keys must be non-NULL, nonempty strings of at most 1024 bytes";
            break;
        }
        columns->keys[i].length = strlen(name);
        if (duckhts_column_names_claim(&names, name, NULL, detail, sizeof(detail)) != 0) {
            snprintf(message, sizeof(message), "attributes: %s", detail);
            error = message;
            break;
        }
    }
done:
    duckhts_column_names_destroy(&names);
    duckdb_destroy_value(&list);
    if (error) duckdb_bind_set_error(info, error);
    return error == NULL;
}

void duckhts_attribute_columns_declare(duckdb_bind_info info,
                                       duckhts_attribute_columns *columns, idx_t first_column) {
    columns->first_column = first_column;
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    for (idx_t i = 0; i < columns->count; i++) {
        duckdb_bind_add_result_column(info, columns->keys[i].name, varchar_type);
    }
    duckdb_destroy_logical_type(&varchar_type);
}

duckhts_projected_attribute *duckhts_attribute_columns_project(
    const duckhts_attribute_columns *columns, const idx_t *column_ids,
    idx_t column_count, idx_t *projected_count) {
    *projected_count = 0;
    for (idx_t i = 0; i < column_count; i++) {
        if (column_ids[i] >= columns->first_column &&
            column_ids[i] - columns->first_column < columns->count) (*projected_count)++;
    }
    if (!*projected_count || *projected_count > SIZE_MAX / sizeof(duckhts_projected_attribute)) {
        return NULL;
    }
    duckhts_projected_attribute *projected = calloc((size_t)*projected_count, sizeof(*projected));
    if (!projected) return NULL;
    idx_t next = 0;
    for (idx_t i = 0; i < column_count; i++) {
        if (column_ids[i] >= columns->first_column &&
            column_ids[i] - columns->first_column < columns->count) {
            projected[next].key = &columns->keys[column_ids[i] - columns->first_column];
            projected[next].output_column = i;
            next++;
        }
    }
    return projected;
}

void duckhts_attribute_columns_write(const duckhts_projected_attribute *projected,
                                     idx_t projected_count, duckdb_vector *vectors, idx_t row) {
    for (idx_t i = 0; i < projected_count; i++) {
        duckdb_vector vector = vectors[projected[i].output_column];
        if (projected[i].value) {
            duckdb_vector_assign_string_element_len(vector, row, projected[i].value,
                                                     projected[i].value_length);
        } else {
            duckdb_vector_ensure_validity_writable(vector);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vector), row);
        }
    }
}
