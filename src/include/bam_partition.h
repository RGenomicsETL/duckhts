#ifndef DUCKHTS_BAM_PARTITION_H
#define DUCKHTS_BAM_PARTITION_H

/*
 * Coordinate partitions of an indexed BAM scan.
 *
 * A scan that claims whole references runs on one thread when the file holds
 * one reference, or one large one. These partitions cut a reference whose
 * indexed data exceeds `piece_bytes` into byte-balanced pieces, using the
 * file offsets the index reports. Each piece is a half-open position range
 * [beg, end) of one reference. A record belongs to the piece that holds its
 * start position (duckhts_bam_partition_owns), so a record that crosses a
 * boundary is read by both pieces' iterators and kept by exactly one.
 *
 * The three arrays share one allocation (beg points at its start), so the
 * claim loop reads flat arrays indexed by one atomic counter. Building them
 * queries the index only; it reads no record.
 */

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdint.h>

typedef struct {
    int n;
    int capacity;
    int64_t *beg;
    int64_t *end;
    int32_t *tid;
} duckhts_bam_partitions;

/*
 * Builds one piece per reference of `hdr`, in reference order. With
 * piece_bytes > 0, a reference with more indexed compressed bytes than
 * piece_bytes becomes up to max_pieces byte-balanced pieces. With align > 1,
 * each interior boundary is rounded down to a multiple of align, so a caller
 * that accumulates fixed-width bins never shares a bin between pieces.
 * piece_bytes must be 0 for an index whose offsets are not BGZF virtual
 * offsets (CRAM). Returns 1, or 0 when out of memory; on failure `out` is
 * freed. `out` must be zero-initialized or freed before the call.
 */
int duckhts_bam_partitions_build(hts_idx_t *idx, const sam_hdr_t *hdr,
                                 int64_t piece_bytes, int max_pieces, int64_t align,
                                 duckhts_bam_partitions *out);

void duckhts_bam_partitions_free(duckhts_bam_partitions *parts);

/* Whether a record read by the iterator of a piece starting at `beg` belongs
 * to that piece. Records without a reference belong to their own partition.
 * The first piece of a reference (beg 0) keeps every record its iterator
 * returns, including a record with a reference and no position (POS 0). */
static inline int duckhts_bam_partition_owns(const bam1_t *rec, int64_t beg) {
    return beg <= 0 || rec->core.tid < 0 || rec->core.pos >= beg;
}

#endif
