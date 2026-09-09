/* INFO and FORMAT materialization preserves encoded value positions. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN
#include "include/bcf_field_vector.h"
#include "include/duckdb_list.h"
#include <htslib/vcf.h>
#include <string.h>

duckdb_logical_type duckhts_bcf_field_type(int type, int is_list) {
    duckdb_type id = type == BCF_HT_INT ? DUCKDB_TYPE_INTEGER :
        type == BCF_HT_REAL ? DUCKDB_TYPE_FLOAT :
        type == BCF_HT_FLAG ? DUCKDB_TYPE_BOOLEAN : DUCKDB_TYPE_VARCHAR;
    duckdb_logical_type element = duckdb_create_logical_type(id);
    if (!is_list) return element;
    duckdb_logical_type list = duckdb_create_list_type(element);
    duckdb_destroy_logical_type(&element);
    return list;
}

void duckhts_bcf_field_null(duckdb_vector vector, idx_t row, int is_list) {
    duckdb_vector_ensure_validity_writable(vector);
    duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vector), row);
    if (is_list) ((duckdb_list_entry *)duckdb_vector_get_data(vector))[row] =
        (duckdb_list_entry){duckdb_list_vector_get_size(vector), 0};
}

int duckhts_bcf_field_int32(duckdb_vector vector, idx_t row,
                           const int32_t *values, int count, int is_list) {
    if (!values || count <= 0) {
        duckhts_bcf_field_null(vector, row, is_list);
        return 1;
    }
    if (!is_list) {
        if (values[0] == bcf_int32_missing || values[0] == bcf_int32_vector_end)
            duckhts_bcf_field_null(vector, row, 0);
        else ((int32_t *)duckdb_vector_get_data(vector))[row] = values[0];
        return 1;
    }
    int length = 0;
    while (length < count && values[length] != bcf_int32_vector_end) length++;
    duckdb_list_entry entry;
    if (!duckhts_list_extend(vector, (idx_t)length, &entry)) return 0;
    duckdb_vector child = duckdb_list_vector_get_child(vector);
    int32_t *data = duckdb_vector_get_data(child);
    for (int i = 0; i < length; i++) {
        if (values[i] == bcf_int32_missing) duckhts_bcf_field_null(child, entry.offset + i, 0);
        else data[entry.offset + i] = values[i];
    }
    ((duckdb_list_entry *)duckdb_vector_get_data(vector))[row] = entry;
    return 1;
}

int duckhts_bcf_field_float(duckdb_vector vector, idx_t row,
                           const float *values, int count, int is_list) {
    if (!values || count <= 0) {
        duckhts_bcf_field_null(vector, row, is_list);
        return 1;
    }
    if (!is_list) {
        if (bcf_float_is_missing(values[0]) || bcf_float_is_vector_end(values[0]))
            duckhts_bcf_field_null(vector, row, 0);
        else ((float *)duckdb_vector_get_data(vector))[row] = values[0];
        return 1;
    }
    int length = 0;
    while (length < count && !bcf_float_is_vector_end(values[length])) length++;
    duckdb_list_entry entry;
    if (!duckhts_list_extend(vector, (idx_t)length, &entry)) return 0;
    duckdb_vector child = duckdb_list_vector_get_child(vector);
    float *data = duckdb_vector_get_data(child);
    for (int i = 0; i < length; i++) {
        if (bcf_float_is_missing(values[i])) duckhts_bcf_field_null(child, entry.offset + i, 0);
        else data[entry.offset + i] = values[i];
    }
    ((duckdb_list_entry *)duckdb_vector_get_data(vector))[row] = entry;
    return 1;
}

int duckhts_bcf_format_write(duckdb_vector vector, idx_t row,
                             const duckhts_bcf_format_t *values, int header_type,
                             int is_list, int sample) {
    if (!values->count) {
        duckhts_bcf_field_null(vector, row, is_list);
        return 1;
    }
    size_t offset = (size_t)sample * values->stride;
    if (header_type == BCF_HT_INT)
        return duckhts_bcf_field_int32(vector, row, (const int32_t *)values->data + offset,
                                       values->stride, is_list);
    if (header_type == BCF_HT_REAL)
        return duckhts_bcf_field_float(vector, row, (const float *)values->data + offset,
                                       values->stride, is_list);
    return duckhts_bcf_field_string(vector, row, values->strings[sample], is_list);
}

int duckhts_bcf_field_string(duckdb_vector vector, idx_t row,
                            const char *value, int is_list) {
    if (!value || strcmp(value, ".") == 0) {
        duckhts_bcf_field_null(vector, row, is_list);
        return 1;
    }
    if (!is_list) {
        duckdb_vector_assign_string_element(vector, row, value);
        return 1;
    }
    idx_t count = *value ? 1 : 0;
    for (const char *p = value; *p; p++) if (*p == ',') count++;
    duckdb_list_entry entry;
    if (!duckhts_list_extend(vector, count, &entry)) return 0;
    duckdb_vector child = duckdb_list_vector_get_child(vector);
    const char *start = value;
    for (idx_t i = 0; i < count; i++) {
        const char *end = strchr(start, ',');
        size_t length = end ? (size_t)(end - start) : strlen(start);
        if (length == 1 && start[0] == '.') duckhts_bcf_field_null(child, entry.offset + i, 0);
        else duckdb_vector_assign_string_element_len(child, entry.offset + i, start, length);
        start += length + (end != NULL);
    }
    ((duckdb_list_entry *)duckdb_vector_get_data(vector))[row] = entry;
    return 1;
}
