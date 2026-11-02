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

/*
 * A boundary is a sampled position (two samples per piece, evenly spaced over
 * the reference length) where the indexed data starts at or past the next
 * equal share of the reference's bytes. Deep regions can put several shares in
 * one sample interval; the reference then gets fewer, larger pieces.
 */
int duckhts_bam_partitions_build(hts_idx_t *idx, const sam_hdr_t *hdr,
                                 int64_t piece_bytes, int max_pieces, int64_t align,
                                 duckhts_bam_partitions *out) {
    int n_refs = sam_hdr_nref(hdr);

    duckhts_bam_partitions_free(out);
    if (max_pieces < 1) max_pieces = 1;
    for (int tid = 0; tid < n_refs; tid++) {
        hts_pos_t length = sam_hdr_tid2len(hdr, tid);
        uint64_t first = 0, last = 0;
        int64_t previous = 0;
        int pieces = 1;
        int made = 1;

        if (piece_bytes > 0 && length > 1 &&
            index_byte_span(idx, tid, 0, length, &first, &last) && last > first) {
            uint64_t wanted = (last - first) / (uint64_t)piece_bytes + 1;
            pieces = wanted > (uint64_t)max_pieces ? max_pieces : (int)wanted;
        }
        for (int k = 1; pieces > 1 && k < 2 * pieces && made < pieces; k++) {
            int64_t pos = (int64_t)length * k / (2 * pieces);
            uint64_t at = 0, unused = 0;
            uint64_t share = first + (last - first) * (uint64_t)made / (uint64_t)pieces;

            if (align > 1) pos -= pos % align;
            if (pos <= previous) continue;
            if (!index_byte_span(idx, tid, pos, length, &at, &unused)) break;
            if (at < share) continue;
            if (!partitions_push(out, tid, previous, pos)) goto oom;
            previous = pos;
            made++;
        }
        if (!partitions_push(out, tid, previous, HTS_POS_MAX)) goto oom;
    }
    return 1;

oom:
    duckhts_bam_partitions_free(out);
    return 0;
}
