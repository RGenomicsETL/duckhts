#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <stdint.h>
#include <stddef.h>
#include <string.h>

#define FIND_SITES_MAX_CANDIDATES 1000000u

typedef struct {
    uint64_t cell;
    uint64_t position;
    bool occupied;
} spacing_cell_t;

static bool valid_row(duckdb_vector vector, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    return validity == NULL || duckdb_validity_row_is_valid(validity, row);
}

static spacing_cell_t *lookup_cell(spacing_cell_t *cells, size_t mask,
                                   uint64_t cell) {
    size_t index = (size_t)(cell * UINT64_C(11400714819323198485)) & mask;
    while (cells[index].occupied && cells[index].cell != cell) {
        index = (index + 1u) & mask;
    }
    return &cells[index];
}

static void spacing_scalar(duckdb_function_info info, duckdb_data_chunk input,
                           duckdb_vector output) {
    duckdb_vector positions = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector distance = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector child = duckdb_list_vector_get_child(positions);
    duckdb_list_entry *lists = duckdb_vector_get_data(positions);
    uint64_t *values = duckdb_vector_get_data(child);
    uint64_t *distances = duckdb_vector_get_data(distance);
    idx_t child_size = duckdb_list_vector_get_size(positions);
    idx_t rows = duckdb_data_chunk_get_size(input);
    idx_t total = 0u;

    for (idx_t row = 0u; row < rows; row++) {
        duckdb_list_entry list = lists[row];
        if (!valid_row(positions, row) || !valid_row(distance, row) ||
            distances[row] == 0u || list.offset > child_size ||
            list.length > child_size - list.offset ||
            list.length > FIND_SITES_MAX_CANDIDATES ||
            list.length > FIND_SITES_MAX_CANDIDATES - total) {
            duckdb_scalar_function_set_error(info,
                "duckhts_somalier_spacing: require non-NULL positions, positive distance, and at most 1000000 candidates per chunk");
            return;
        }
        for (idx_t i = 0u; i < list.length; i++) {
            if (!valid_row(child, list.offset + i) ||
                values[list.offset + i] == 0u) {
                duckdb_scalar_function_set_error(info,
                    "duckhts_somalier_spacing: positions must be positive and non-NULL");
                return;
            }
        }
        total += list.length;
    }
    if (duckdb_list_vector_reserve(output, total) != DuckDBSuccess ||
        duckdb_list_vector_set_size(output, total) != DuckDBSuccess) {
        duckdb_scalar_function_set_error(info,
            "duckhts_somalier_spacing: could not reserve output");
        return;
    }
    duckdb_list_entry *result = duckdb_vector_get_data(output);
    duckdb_vector result_child = duckdb_list_vector_get_child(output);
    bool *selected = duckdb_vector_get_data(result_child);
    idx_t offset = 0u;
    for (idx_t row = 0u; row < rows; row++) {
        duckdb_list_entry list = lists[row];
        uint64_t gap = distances[row];
        size_t capacity = 2u;
        spacing_cell_t *cells;
        while (capacity < (size_t)list.length * 2u) capacity *= 2u;
        cells = duckdb_malloc(capacity * sizeof(*cells));
        if (cells == NULL) {
            duckdb_scalar_function_set_error(info,
                "duckhts_somalier_spacing: allocation failed");
            return;
        }
        memset(cells, 0, capacity * sizeof(*cells));
        result[row].offset = offset;
        result[row].length = list.length;
        for (idx_t i = 0u; i < list.length; i++) {
            uint64_t position = values[list.offset + i];
            uint64_t cell = position / gap;
            bool keep = true;
            uint64_t first = cell == 0u ? 0u : cell - 1u;
            uint64_t last = cell == UINT64_MAX ? cell : cell + 1u;
            for (uint64_t neighbor = first; neighbor <= last; neighbor++) {
                spacing_cell_t *previous = lookup_cell(cells, capacity - 1u,
                                                       neighbor);
                if (previous->occupied) {
                    uint64_t p = previous->position;
                    if ((position >= p ? position - p : p - position) < gap) {
                        keep = false;
                        break;
                    }
                }
                if (neighbor == last) break;
            }
            selected[offset + i] = keep;
            if (keep) {
                spacing_cell_t *slot = lookup_cell(cells, capacity - 1u, cell);
                slot->cell = cell;
                slot->position = position;
                slot->occupied = true;
            }
        }
        duckdb_free(cells);
        offset += list.length;
    }
}

void register_duckhts_somalier_spacing(duckdb_connection connection) {
    duckdb_logical_type position = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type distance = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type positions = duckdb_create_list_type(position);
    duckdb_logical_type result = duckdb_create_list_type(boolean);
    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, "duckhts_somalier_spacing");
    duckdb_scalar_function_add_parameter(function, positions);
    duckdb_scalar_function_add_parameter(function, distance);
    duckdb_scalar_function_set_return_type(function, result);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, spacing_scalar);
    duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);
    duckdb_destroy_logical_type(&position);
    duckdb_destroy_logical_type(&distance);
    duckdb_destroy_logical_type(&boolean);
    duckdb_destroy_logical_type(&positions);
    duckdb_destroy_logical_type(&result);
}
