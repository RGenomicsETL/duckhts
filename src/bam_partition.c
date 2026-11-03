/* Coordinate partitions of an indexed BAM scan; see include/bam_partition.h. */

#include "include/bam_partition.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

void duckhts_bam_partitions_free(duckhts_bam_partitions *parts) {
    if (!parts) return;
    free(parts->beg);
    memset(parts, 0, sizeof(*parts));
}

/* Appends one piece, growing the three arrays together. Returns 0 when out of memory. */
static int partitions_push(duckhts_bam_partitions *parts, int32_t tid, int64_t beg, int64_t end) {
    if (parts->n == parts->capacity) {
        size_t each = 2 * sizeof(int64_t) + sizeof(int32_t);
        int next;
        int64_t *block;

        if (parts->capacity > INT_MAX / 2) return 0;
        next = parts->capacity ? parts->capacity * 2 : 64;
        if ((size_t)next > SIZE_MAX / each) return 0;
        block = malloc((size_t)next * each);
        if (!block) return 0;
        if (parts->n) {
            memcpy(block, parts->beg, (size_t)parts->n * sizeof(int64_t));
            memcpy(block + next, parts->end, (size_t)parts->n * sizeof(int64_t));
            memcpy((int32_t *)(block + 2 * (size_t)next), parts->tid, (size_t)parts->n * sizeof(int32_t));
        }
        free(parts->beg);
        parts->beg = block;
        parts->end = block + next;
        parts->tid = (int32_t *)(block + 2 * (size_t)next);
        parts->capacity = next;
    }
    parts->tid[parts->n] = tid;
    parts->beg[parts->n] = beg;
    parts->end[parts->n] = end;
    parts->n++;
    return 1;
}

/* Compressed file offsets (BGZF block addresses) of the first and last indexed
 * chunk overlapping [beg, end) on `tid`. Returns 1 when the index has data there. */
static int index_byte_span(hts_idx_t *idx, int tid, hts_pos_t beg, hts_pos_t end,
                           uint64_t *first, uint64_t *last) {
    hts_itr_t *itr = sam_itr_queryi(idx, tid, beg, end);
    int found = 0;

    if (!itr) return 0;
    for (int i = 0; i < itr->n_off; i++) {
        uint64_t u = itr->off[i].u >> 16;
        uint64_t v = itr->off[i].v >> 16;
        if (!found || u < *first) *first = u;
        if (!found || v > *last) *last = v;
        found = 1;
    }
    hts_itr_destroy(itr);
    return found;
}

/* The first indexed byte at or after `pos`: the smallest chunk start over a window
 * that begins at pos. In a coordinate-sorted file a record overlapping a later
 * position either starts after the window, at a larger offset, or also overlaps
 * the window, so a narrow window answers; it widens only across a stretch without
 * data. Returns 0 when no data lies at or after pos. */
static int index_first_byte_at(hts_idx_t *idx, int tid, int64_t pos, int64_t length,
                               uint64_t *at) {
    uint64_t unused = 0;

    for (int64_t width = 1 << 14;; width *= 8) {
        int64_t end = width >= length - pos ? length : pos + width;
        if (index_byte_span(idx, tid, pos, end, at, &unused)) return 1;
        if (end >= length) return 0;
    }
}

/*
 * Boundary j is the smallest position, to within BOUNDARY_GRANULE bases (or the
 * alignment), whose first indexed byte reaches the j-th equal share of the
 * reference's bytes, found by bisection. Bisection follows the data wherever it
 * lies, so a reference whose records occupy a short prefix is still split. A share
 * that falls inside one granule of a deep pileup adds no boundary; the reference
 * then gets fewer, larger pieces.
 */
#define BOUNDARY_GRANULE ((int64_t)1 << 14)

int duckhts_bam_partitions_build(hts_idx_t *idx, const sam_hdr_t *hdr,
                                 int64_t piece_bytes, int max_pieces, int64_t align,
                                 duckhts_bam_partitions *out) {
    int n_refs = sam_hdr_nref(hdr);
    int64_t granule = align > BOUNDARY_GRANULE ? align : BOUNDARY_GRANULE;

    duckhts_bam_partitions_free(out);
    if (max_pieces < 1) max_pieces = 1;
    for (int tid = 0; tid < n_refs; tid++) {
        hts_pos_t length = sam_hdr_tid2len(hdr, tid);
        uint64_t first = 0, last = 0;
        int64_t previous = 0;
        int pieces = 1;

        if (piece_bytes > 0 && length > 1 &&
            index_byte_span(idx, tid, 0, length, &first, &last) && last > first) {
            uint64_t wanted = (last - first) / (uint64_t)piece_bytes + 1;
            pieces = wanted > (uint64_t)max_pieces ? max_pieces : (int)wanted;
        }
        for (int j = 1; j < pieces; j++) {
            uint64_t share = first + (last - first) * (uint64_t)j / (uint64_t)pieces;
            int64_t lo = previous + 1;
            int64_t hi = length;
            int64_t pos;

            while (hi - lo > granule) {
                int64_t mid = lo + (hi - lo) / 2;
                uint64_t at = 0;
                if (index_first_byte_at(idx, tid, mid, length, &at) && at < share) lo = mid + 1;
                else hi = mid;
            }
            pos = hi;
            if (align > 1) pos -= pos % align;
            if (pos >= length) break; /* no data reaches this share */
            if (pos <= previous) continue;
            if (!partitions_push(out, tid, previous, pos)) goto oom;
            previous = pos;
        }
        if (!partitions_push(out, tid, previous, HTS_POS_MAX)) goto oom;
    }
    return 1;

oom:
    duckhts_bam_partitions_free(out);
    return 0;
}
