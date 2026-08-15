#ifndef DUCKHTS_LIFTOVER_NW_LIMIT_H
#define DUCKHTS_LIFTOVER_NW_LIMIT_H

#include <stddef.h>

#define LIFTOVER_NW_MAX_CELLS ((size_t)4 * 1024 * 1024)

static inline int liftover_nw_cell_count(size_t s_len, size_t t_len, size_t *cells) {
    size_t rows, columns;
    if (!cells || s_len >= LIFTOVER_NW_MAX_CELLS || t_len >= LIFTOVER_NW_MAX_CELLS) return 0;
    rows = s_len + 1;
    columns = t_len + 1;
    if (rows > LIFTOVER_NW_MAX_CELLS / columns) return 0;
    *cells = rows * columns;
    return 1;
}

#endif
