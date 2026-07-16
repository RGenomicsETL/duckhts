/*
 * Reusable pure-C variant storage for sorted annotation sweeps. This file has no
 * DuckDB dependency and performs no sequence-region name lookup.
 *
 * Transcript models remain resident while variants arrive in bounded groups. The
 * tile owns the kernel's structure-of-arrays view, allele bytes, and variant IDs;
 * callers may therefore release their input chunk after appending it. The tile caps
 * only the number of variants. The annotation kernel separately bounds and resumes
 * result emission because transcript overlap determines the result count.
 *
 * Sequence region is an already-resolved `uint16` ordinal. The tile enforces:
 *   1. all rows for a contig are contiguous (one run per tile);
 *   2. pos1 is nondecreasing within a run;
 *   3. a contig must not reappear after its run has closed.
 * Run-validation state persists across reset_rows, so ordering is checked across tile
 * boundaries as well as within one tile.
 *
 * Current scope includes validated VCF-shaped small variants plus single-interval
 * structural events. Structural operation and copy direction are explicit fields;
 * breakends use the separate paired-locus annotation path.
 */
#ifndef DUCKVEP_VARIANT_TILE_H
#define DUCKVEP_VARIANT_TILE_H

#include <stddef.h>
#include <stdint.h>

#include "duckvep_kernel.h" /* duckvep_variant_batch_t, duckvep_variant_kind_t */

typedef struct duckvep_variant_tile duckvep_variant_tile_t;

typedef enum duckvep_tile_status {
    DUCKVEP_TILE_APPENDED = 0,   /* row stored */
    DUCKVEP_TILE_FULL_CAPACITY,  /* NOT stored: flush this tile (SAME chrom continues), then retry */
    DUCKVEP_TILE_FULL_CHROM,     /* NOT stored: flush + close_run (NEW chrom begins), then retry */
    DUCKVEP_TILE_INVALID         /* NOT stored: hard validation error (see _error) */
} duckvep_tile_status_t;

/* Tile holding up to `variant_capacity` rows over a contig-id space of `max_chroms`
 * (the interner count — bounds the closed-set bitset). NULL on OOM or zero args. */
duckvep_variant_tile_t *duckvep_variant_tile_create(size_t variant_capacity,
                                                    size_t max_chroms);
void duckvep_variant_tile_destroy(duckvep_variant_tile_t *tile);

/* Append one fully-read, already-interned ordinary variant row. SV rows passed
 * here receive UNKNOWN subtype/direction; adapters with structural metadata should
 * call duckvep_variant_tile_append_event. On FULL_* the row is NOT stored. */
duckvep_tile_status_t duckvep_variant_tile_append(
    duckvep_variant_tile_t *tile, uint16_t chrom_id, uint32_t pos1, uint32_t end1,
    uint8_t kind, const char *ref, size_t ref_len, const char *alt, size_t alt_len,
    const char *variant_id, size_t variant_id_len);

/* Extended append preserving structural operation and copy-number direction. */
duckvep_tile_status_t duckvep_variant_tile_append_event(
    duckvep_variant_tile_t *tile, uint16_t chrom_id, uint32_t pos1, uint32_t end1,
    uint8_t kind, uint8_t sv_type, uint8_t copy_change,
    const char *ref, size_t ref_len, const char *alt, size_t alt_len,
    const char *variant_id, size_t variant_id_len);

/* Clear stored rows for the next tile, PRESERVING run-validation state (current open
 * chrom, last pos1, closed-set). Use after a FULL_CAPACITY flush. */
void duckvep_variant_tile_reset_rows(duckvep_variant_tile_t *tile);

/* Clear rows AND close the current run: mark the current chrom closed and require the
 * next append to open a fresh run. Use after a FULL_CHROM flush. */
void duckvep_variant_tile_close_run(duckvep_variant_tile_t *tile);

/* Borrowed kernel view over the currently-stored rows (valid until the next mutating
 * call). Pass straight to duckvep_annotate_tile. */
void duckvep_variant_tile_batch(const duckvep_variant_tile_t *tile,
                                duckvep_variant_batch_t *out);

size_t duckvep_variant_tile_count(const duckvep_variant_tile_t *tile);

/* variant_id metadata for an emitted `variant_idx` (index into the batch). Returns the
 * owned id bytes (NOT NUL-terminated), sets *len; NULL if idx is out of range. */
const char *duckvep_variant_tile_variant_id(const duckvep_variant_tile_t *tile,
                                            size_t idx, size_t *len);

/* Last validation-error string (owned by the tile), or NULL. */
const char *duckvep_variant_tile_error(const duckvep_variant_tile_t *tile);

#endif /* DUCKVEP_VARIANT_TILE_H */
