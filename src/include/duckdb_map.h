#ifndef DUCKHTS_DUCKDB_MAP_H
#define DUCKHTS_DUCKDB_MAP_H

#include "duckdb_list.h"

#include <stddef.h>

/* Row writer for MAP(VARCHAR, VARCHAR) vectors. Include after duckdb_extension.h /
 * DUCKDB_EXTENSION_EXTERN.
 *
 * A MAP vector is a LIST of STRUCT(key, value). begin reserves `capacity` child slots
 * at the end of the child list. Each put copies one key span and one value span into
 * the next reserved slot. finish publishes the pairs put so far as the row's entry.
 * A reserved slot that receives no pair stays in the child list, and no entry refers
 * to it.
 *
 * The caller owns the pair policy: which pairs to write and in what order, repeated
 * keys, escaping, empty values, and NULL versus empty rows. DuckDB stores a NULL key
 * or value for a span that is not valid UTF-8.
 *
 * The vectors are borrowed from the output chunk. Use a duckhts_map_pairs only
 * between its begin and its finish. */
typedef struct {
    duckdb_vector map;
    duckdb_vector keys;
    duckdb_vector values;
    duckdb_list_entry entry; /* offset: first reserved slot; length: pairs put */
    idx_t capacity;          /* reserved slots */
} duckhts_map_pairs;

/* 1 on success. 0 when the reservation fails; the child list keeps its size. */
static inline int duckhts_map_pairs_begin(duckdb_vector map, idx_t capacity,
                                          duckhts_map_pairs *pairs) {
    duckdb_list_entry reserved;
    if (!duckhts_list_extend(map, capacity, &reserved)) return 0;
    duckdb_vector child = duckdb_list_vector_get_child(map);
    pairs->map = map;
    pairs->keys = duckdb_struct_vector_get_child(child, 0);
    pairs->values = duckdb_struct_vector_get_child(child, 1);
    pairs->entry.offset = reserved.offset;
    pairs->entry.length = 0;
    pairs->capacity = capacity;
    return 1;
}

/* 1 on success. 0 without writing when every reserved slot already holds a pair. */
static inline int duckhts_map_pairs_put(duckhts_map_pairs *pairs, const char *key,
                                        size_t key_len, const char *value,
                                        size_t value_len) {
    if (pairs->entry.length == pairs->capacity) return 0;
    idx_t slot = pairs->entry.offset + pairs->entry.length;
    duckdb_vector_assign_string_element_len(pairs->keys, slot, key, (idx_t)key_len);
    duckdb_vector_assign_string_element_len(pairs->values, slot, value, (idx_t)value_len);
    pairs->entry.length++;
    return 1;
}

static inline void duckhts_map_pairs_finish(const duckhts_map_pairs *pairs, idx_t row) {
    ((duckdb_list_entry *)duckdb_vector_get_data(pairs->map))[row] = pairs->entry;
}

/* Marks a row of a LIST or MAP vector NULL. The row's entry has length 0 and starts
 * at the current child size, so it never points past the child list. */
static inline void duckhts_list_set_null(duckdb_vector list, idx_t row) {
    duckdb_vector_ensure_validity_writable(list);
    duckdb_validity_set_row_invalid(duckdb_vector_get_validity(list), row);
    duckdb_list_entry entry = {duckdb_list_vector_get_size(list), 0};
    ((duckdb_list_entry *)duckdb_vector_get_data(list))[row] = entry;
}

#endif
