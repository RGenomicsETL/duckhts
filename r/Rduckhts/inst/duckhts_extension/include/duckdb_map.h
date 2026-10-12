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
 * keys, escaping, empty values, and NULL versus empty rows. A MAP key may not be
 * NULL, so put leaves out a pair whose key is not valid UTF-8 (the slot stays
 * unused); DuckDB stores a NULL value for a value that is not valid UTF-8.
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

/* 1 when the bytes are UTF-8 as RFC 3629 defines it, which is what DuckDB accepts:
 * no overlong form, no surrogate, nothing above U+10FFFF. NUL bytes are valid. */
static inline int duckhts_utf8_valid(const char *text, size_t len) {
    const unsigned char *s = (const unsigned char *)text;
    size_t i = 0;
    while (i < len) {
        unsigned char c = s[i];
        size_t extra;
        unsigned char low = 0x80, high = 0xBF; /* range of the second byte */
        if (c < 0x80) {
            i++;
            continue;
        } else if (c >= 0xC2 && c <= 0xDF) {
            extra = 1;
        } else if (c >= 0xE0 && c <= 0xEF) {
            extra = 2;
            if (c == 0xE0) low = 0xA0;  /* overlong */
            if (c == 0xED) high = 0x9F; /* surrogates */
        } else if (c >= 0xF0 && c <= 0xF4) {
            extra = 3;
            if (c == 0xF0) low = 0x90;  /* overlong */
            if (c == 0xF4) high = 0x8F; /* above U+10FFFF */
        } else {
            return 0;
        }
        if (len - i <= extra) return 0;
        if (s[i + 1] < low || s[i + 1] > high) return 0;
        for (size_t k = 2; k <= extra; k++) {
            if (s[i + k] < 0x80 || s[i + k] > 0xBF) return 0;
        }
        i += extra + 1;
    }
    return 1;
}

/* 1 on success, including a pair left out because its key is not valid UTF-8.
 * 0 without writing when every reserved slot already holds a pair. */
static inline int duckhts_map_pairs_put(duckhts_map_pairs *pairs, const char *key,
                                        size_t key_len, const char *value,
                                        size_t value_len) {
    if (!duckhts_utf8_valid(key, key_len)) return 1;
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
