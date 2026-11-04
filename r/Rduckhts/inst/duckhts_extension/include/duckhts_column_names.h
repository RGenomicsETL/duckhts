#ifndef DUCKHTS_COLUMN_NAMES_H
#define DUCKHTS_COLUMN_NAMES_H

#include <stddef.h>

/* Result column names as DuckDB compares identifiers: ASCII letters fold to lower
 * case. A bind-time set; every name claimed is copied and owned by the set. */
typedef struct {
    void *names;           /* khash: folded name -> claimed name, both owned */
    int suffix_duplicates; /* rename a colliding name to name__2, name__3, ... */
} duckhts_column_names;

void duckhts_column_names_init(duckhts_column_names *set, int suffix_duplicates);
void duckhts_column_names_destroy(duckhts_column_names *set);

/* Claim `name` for a result column. On success returns 0 and, when `claimed` is not
 * NULL, points it at the name to declare (borrowed from the set): `name` itself, or
 * a suffixed copy under the suffix policy. Returns -1 with `err` set when the name
 * collides and duplicates are not suffixed, and -2 with `err` set on allocation
 * failure. */
int duckhts_column_names_claim(duckhts_column_names *set, const char *name,
                               const char **claimed, char *err, size_t err_size);

#endif
