#include "duckhts_column_names.h"

#include <htslib/khash.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

KHASH_MAP_INIT_STR(duckhts_colname, char *)

/* DuckDB identifier comparisons fold ASCII letters, not the host locale. */
static char *fold_copy(const char *name) {
    size_t length = strlen(name);
    char *folded = malloc(length + 1);
    if (!folded) return NULL;
    for (size_t i = 0; i <= length; i++) {
        unsigned char c = (unsigned char)name[i];
        folded[i] = (char)(c >= 'A' && c <= 'Z' ? c + ('a' - 'A') : c);
    }
    return folded;
}

void duckhts_column_names_init(duckhts_column_names *set, int suffix_duplicates) {
    set->names = NULL;
    set->suffix_duplicates = suffix_duplicates;
}

void duckhts_column_names_destroy(duckhts_column_names *set) {
    khash_t(duckhts_colname) *names = set->names;
    if (names) {
        for (khiter_t k = kh_begin(names); k != kh_end(names); k++) {
            if (!kh_exist(names, k)) continue;
            free((char *)kh_key(names, k));
            free(kh_val(names, k));
        }
        kh_destroy(duckhts_colname, names);
    }
    set->names = NULL;
}

/* Insert `name` unless its folded form is present. Returns 1 when inserted, 0 when
 * present (with *existing set), -1 on allocation failure. */
static int insert_name(khash_t(duckhts_colname) *names, const char *name,
                       const char **inserted, const char **existing) {
    char *folded = fold_copy(name);
    int absent;
    if (!folded) return -1;
    khiter_t k = kh_put(duckhts_colname, names, folded, &absent);
    if (absent < 0) {
        free(folded);
        return -1;
    }
    if (!absent) {
        free(folded);
        *existing = kh_val(names, k);
        return 0;
    }
    size_t size = strlen(name) + 1;
    char *copy = malloc(size);
    if (copy) memcpy(copy, name, size);
    if (!copy) {
        kh_del(duckhts_colname, names, k);
        free(folded);
        return -1;
    }
    kh_val(names, k) = copy;
    *inserted = copy;
    return 1;
}

int duckhts_column_names_claim(duckhts_column_names *set, const char *name,
                               const char **claimed, char *err, size_t err_size) {
    const char *inserted = NULL;
    const char *existing = NULL;
    if (!set->names && !(set->names = kh_init(duckhts_colname))) goto out_of_memory;
    int rc = insert_name(set->names, name, &inserted, &existing);
    if (rc < 0) goto out_of_memory;
    if (rc == 0 && !set->suffix_duplicates) {
        snprintf(err, err_size, "%s and %s name the same column (column names ignore case)",
                 existing, name);
        return -1;
    }
    if (rc == 0) {
        size_t size = strlen(name) + sizeof("__4294967295");
        char *suffixed = malloc(size);
        if (!suffixed) goto out_of_memory;
        /* Every pass adds a name or finds one already claimed, so this ends. */
        for (unsigned n = 2; rc == 0; n++) {
            snprintf(suffixed, size, "%s__%u", name, n);
            rc = insert_name(set->names, suffixed, &inserted, &existing);
        }
        free(suffixed);
        if (rc < 0) goto out_of_memory;
    }
    if (claimed) *claimed = inserted;
    return 0;

out_of_memory:
    snprintf(err, err_size, "out of memory checking column names");
    return -2;
}
