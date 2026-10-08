/* Packed histogram cells of allele read counts, for one count-error fit.
 *
 * __duckhts_count_cells() collects the cells of one sample into these records
 * and returns them as one BLOB: a header, then the cells. A cell is the number
 * of sites of one block that have one (depth, alt count) pair in one frequency
 * bin, with the mean frequency of those sites. __duckhts_count_error_fit()
 * reads that BLOB. The layout is private to those two functions and is never
 * persisted. Cells with one key may repeat; the fit sums them.
 *
 * This header has no DuckDB dependency.
 */
#ifndef DUCKHTS_COUNT_CELLS_H
#define DUCKHTS_COUNT_CELLS_H

#include <stdint.h>

#define DUCKHTS_COUNT_CELLS_MAGIC "DHCC"
#define DUCKHTS_COUNT_CELLS_MAGIC_LENGTH 4

typedef struct {
    char magic[DUCKHTS_COUNT_CELLS_MAGIC_LENGTH];
    uint32_t reserved;
} duckhts_count_cells_header_t;

typedef struct {
    uint32_t sites; /* sites in the cell */
    float af;       /* mean frequency of the counted allele over those sites, in (0, 1) */
    uint16_t block; /* block of the genome, from 0 */
    uint16_t depth; /* reads at each site */
    uint16_t alt;   /* reads of the counted allele at each site, at most depth */
    uint16_t reserved;
} duckhts_count_cell_t;

_Static_assert(sizeof(duckhts_count_cells_header_t) == 8, "cell header is 8 bytes");
_Static_assert(sizeof(duckhts_count_cell_t) == 16, "a cell is 16 bytes");

/* Limits of the packed histogram. A BLOB length is 32 bits. */
#define DUCKHTS_COUNT_MAX_CELLS 100000000
#define DUCKHTS_COUNT_MAX_BLOCKS 65535
#define DUCKHTS_COUNT_MAX_DEPTH 65535

/* Native histogram memory held by this process: the buffers of the aggregate
 * and the working memory of the fits (count_cells.c). Charge adds `bytes`
 * when the total stays within max_cell_bytes and returns 1; otherwise it
 * adds nothing, writes the total that would have been held, and returns 0.
 * Release subtracts what a charge added. */
int duckhts_count_bytes_charge(uint64_t bytes, uint64_t max_cell_bytes, uint64_t *would_hold);
void duckhts_count_bytes_release(uint64_t bytes);

/* Checks the two memory limits of a fit. Returns NULL when they are valid,
 * otherwise the error text. has_* is 0 for a SQL NULL. */
static inline const char *duckhts_count_check_limits(int has_max_cells, int64_t max_cells,
                                                     int has_max_cell_bytes, int64_t max_cell_bytes) {
    if (!has_max_cells || max_cells < 1 || max_cells > DUCKHTS_COUNT_MAX_CELLS) {
        return "max_cells must be between 1 and 100000000";
    }
    if (!has_max_cell_bytes || max_cell_bytes < 1) return "max_cell_bytes must be at least 1";
    return (const char *)0;
}

#endif
