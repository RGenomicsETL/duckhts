/* duckhts_bam_mate_overlap_counts: read 1 and read 2 of one fragment compared
 * at the reference positions where both mates have an aligned base.
 *
 * A real allele sits on the molecule, so both mates show it. A sequencing
 * error usually shows in one mate only. So the rate of disagreement is about
 * twice the instrument error rate, with no reference, mask or contamination
 * model.
 *
 * One pass over coordinate-sorted alignments, one worker. The first mate of a
 * pair whose mate starts inside it is held until its mate arrives, and dropped
 * when the scan passes the mate start without it. At most max_pending mates
 * are held at once. The counters are fixed. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include "include/bam_site_counts.h"
#include "include/hts_io_tuning.h"
#include "include/region_list.h"

enum {
    MATE_QUALITY_MAX = 93,     /* a higher base quality is counted at this quality */
    MATE_QUALITY_MISSING = 94, /* the read stores no base qualities */
    MATE_QUALITIES = 95,
    MATE_BASES = 4,            /* A, C, G, T */
    MATE_VECTOR_COLUMNS = 5,
    /* A BAM name has at most 254 characters and its NUL (BAM_MAX_QNAME_LEN). */
    MATE_QNAME_BYTES = 256,
    MATE_ENTRY_BLOCK = 1024,   /* held-mate entries allocated at once */
    MATE_TABLE_START = 1024,   /* first capacity of the name table, in entries */
    MATE_PENDING_MAX = 100000000
};

#define MATE_CELLS ((size_t)MATE_QUALITIES * MATE_QUALITIES * MATE_BASES * MATE_BASES)
/* Reference bases fetched at once to resolve '=' read bases. */
#define MATE_WINDOW_BASES ((hts_pos_t)1 << 16)

_Static_assert(MATE_QUALITY_MISSING == MATE_QUALITY_MAX + 1 && MATE_QUALITIES == MATE_QUALITY_MISSING + 1,
               "quality axis");
_Static_assert(MATE_QNAME_BYTES > BAM_MAX_QNAME_LEN, "a BAM name and its NUL fit the held name");

typedef struct {
    char *path;
    char *region;
    char *index_path;
    char *reference;
    char *reference_index_path;
    int min_mapq;
    int require_flags;
    int exclude_flags;
    int max_pending;
} mate_bind_t;

typedef struct {
    uint64_t *cells; /* MATE_CELLS counters, owned (malloc family) */
    size_t next_cell;
    int counted;
} mate_state_t;

/* One held mate. The name table borrows qname, so an entry must not move while
 * it is held: entries live in fixed blocks and are reused through a free list.
 * qname is the first member, so the table's qname pointer converts back to
 * the entry (C11 6.7.2.1p15). */
typedef struct {
    char qname[MATE_QNAME_BYTES];
    uint8_t *data;        /* CIGAR, packed sequence, qualities; owned (malloc family) */
    size_t data_capacity;
    hts_pos_t pos;        /* zero-based start of this mate */
    hts_pos_t mpos;       /* zero-based start of its mate */
    uint64_t hash;        /* of qname */
    size_t heap_index;
    uint32_t qname_len;
    uint32_t n_cigar;
    int32_t l_qseq;
    uint16_t flag;
} mate_held_t;

_Static_assert(offsetof(mate_held_t, qname) == 0, "a held qname converts back to its entry");

/* The held mates. Every entry in the blocks is either held (in the heap and
 * the name table) or free (in free_entries). entry_count never exceeds
 * max_pending, so held memory is at most max_pending entries, each with the
 * CIGAR, sequence and qualities of the longest alignment it has held. */
typedef struct {
    mate_held_t **blocks;          /* owned, as is each block (malloc family) */
    size_t block_count;
    size_t entry_count;
    mate_held_t **free_entries;    /* owned; capacity entry_count */
    size_t free_count;
    mate_held_t **heap;            /* owned; capacity entry_count; a min-heap on mpos */
    size_t held_count;
    duckhts_bam_site_overlap_scratch_t table; /* slots owned (calloc) */
    size_t max_pending;
} mate_pending_t;

/* Borrowed view of one mate: the record in hand or a held entry. */
typedef struct {
    const uint32_t *cigar;
    const uint8_t *sequence;  /* htslib 4-bit packed */
    const uint8_t *qualities;
    hts_pos_t pos;            /* zero-based start */
    uint32_t n_cigar;
    int32_t l_qseq;
    int has_qualities;
} mate_view_t;

/* The aligned blocks (CIGAR M, = or X) of one view in reference order. The
 * current block is [block_ref, block_ref + block_length) on the reference,
 * from query offset block_query. */
typedef struct {
    const mate_view_t *view;
    uint32_t next_op;
    hts_pos_t reference;
    int64_t query;
    hts_pos_t block_ref;
    int64_t block_query;
    int64_t block_length;
} mate_blocks_t;

/* Reference bases [beg0, end0) of contig tid, loaded only to resolve a '='
 * read base. fai is NULL without a reference; a '=' base is then left out. */
typedef struct {
    faidx_t *fai;
    int tid;
    hts_pos_t beg0;
    hts_pos_t end0;
    hts_pos_t contig_length;
    char *bases; /* owned by htslib's allocator: free() */
} mate_reference_t;

static size_t mate_cell(int quality1, int quality2, int base1, int base2) {
    return (((size_t)quality1 * MATE_QUALITIES + (size_t)quality2) * MATE_BASES + (size_t)base1) * MATE_BASES +
           (size_t)base2;
}

static int base_code(char base) {
    switch (base) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
    }
}

/* htslib 4-bit base (1, 2, 4, 8 for A, C, G, T) to 0..3, or -1. */
static int nt16_code(int nt16) {
    switch (nt16) {
    case 1: return 0;
    case 2: return 1;
    case 4: return 2;
    case 8: return 3;
    default: return -1;
    }
}

/* ---- held mates ---- */

static mate_held_t *held_from_qname(const char *qname) {
    return (mate_held_t *)(void *)(uintptr_t)qname;
}

static void heap_place(mate_pending_t *pending, size_t at, mate_held_t *entry) {
    pending->heap[at] = entry;
    entry->heap_index = at;
}

static void heap_sift_up(mate_pending_t *pending, size_t at) {
    mate_held_t *entry = pending->heap[at];
    while (at > 0) {
        size_t parent = (at - 1) / 2;
        if (pending->heap[parent]->mpos <= entry->mpos) break;
        heap_place(pending, at, pending->heap[parent]);
        at = parent;
    }
    heap_place(pending, at, entry);
}

static void heap_sift_down(mate_pending_t *pending, size_t at) {
    mate_held_t *entry = pending->heap[at];
    for (;;) {
        size_t child = 2 * at + 1;
        if (child >= pending->held_count) break;
        if (child + 1 < pending->held_count && pending->heap[child + 1]->mpos < pending->heap[child]->mpos) child++;
        if (entry->mpos <= pending->heap[child]->mpos) break;
        heap_place(pending, at, pending->heap[child]);
        at = child;
    }
    heap_place(pending, at, entry);
}

static void heap_remove(mate_pending_t *pending, mate_held_t *entry) {
    mate_held_t *last = pending->heap[--pending->held_count];
    size_t at = entry->heap_index;
    if (last == entry) return;
    heap_place(pending, at, last);
    heap_sift_up(pending, at);
    heap_sift_down(pending, last->heap_index);
}

static int pending_init(mate_pending_t *pending, size_t max_pending) {
    size_t entries = max_pending < MATE_TABLE_START ? max_pending : MATE_TABLE_START;
    size_t slots;
    memset(pending, 0, sizeof(*pending));
    pending->max_pending = max_pending;
    if (!duckhts_bam_site_overlap_slot_capacity(entries, &slots)) return 0;
    pending->table.slots = calloc(slots, sizeof(*pending->table.slots));
    if (!pending->table.slots) return 0;
    pending->table.max_entries = entries;
    pending->table.slot_capacity = slots;
    return duckhts_bam_site_overlap_begin(&pending->table);
}

static void pending_destroy(mate_pending_t *pending) {
    for (size_t block = 0; block < pending->block_count; block++) {
        mate_held_t *entries = pending->blocks[block];
        size_t count = pending->entry_count - block * MATE_ENTRY_BLOCK;
        if (count > MATE_ENTRY_BLOCK) count = MATE_ENTRY_BLOCK;
        for (size_t at = 0; at < count; at++) free(entries[at].data);
        free(entries);
    }
    free(pending->blocks);
    free(pending->free_entries);
    free(pending->heap);
    free(pending->table.slots);
    memset(pending, 0, sizeof(*pending));
}

/* Take a held entry out of the name table and the heap. The caller frees it. */
static void pending_unlink(mate_pending_t *pending, mate_held_t *entry) {
    size_t slot;
    if (duckhts_bam_site_overlap_find(&pending->table, entry->qname, entry->qname_len, entry->hash, &slot)) {
        duckhts_bam_site_overlap_remove(&pending->table, slot);
    }
    heap_remove(pending, entry);
}

static void pending_free(mate_pending_t *pending, mate_held_t *entry) {
    pending->free_entries[pending->free_count++] = entry;
}

/* Drop every held mate whose mate starts before pos0: in coordinate order
 * that mate has been passed without a match. */
static void pending_drop_before(mate_pending_t *pending, hts_pos_t pos0) {
    while (pending->held_count > 0 && pending->heap[0]->mpos < pos0) {
        mate_held_t *entry = pending->heap[0];
        pending_unlink(pending, entry);
        pending_free(pending, entry);
    }
}

/* Drop every held mate, at a change of contig. */
static int pending_clear(mate_pending_t *pending) {
    while (pending->held_count > 0) pending_free(pending, pending->heap[--pending->held_count]);
    return duckhts_bam_site_overlap_begin(&pending->table);
}

/* The held mate with this name, or NULL. It stays held. */
static mate_held_t *pending_find(const mate_pending_t *pending, const char *qname, size_t qname_len,
                                 uint64_t hash) {
    size_t slot;
    if (!duckhts_bam_site_overlap_find(&pending->table, qname, qname_len, hash, &slot)) return NULL;
    return held_from_qname(pending->table.slots[slot].qname);
}

/* Double the name table, up to max_pending entries, and insert the held names again. */
static int pending_grow_table(mate_pending_t *pending) {
    duckhts_bam_site_overlap_scratch_t grown;
    size_t entries = pending->table.max_entries * 2;
    if (entries > pending->max_pending) entries = pending->max_pending;
    memset(&grown, 0, sizeof(grown));
    if (!duckhts_bam_site_overlap_slot_capacity(entries, &grown.slot_capacity)) return 0;
    grown.slots = calloc(grown.slot_capacity, sizeof(*grown.slots));
    if (!grown.slots) return 0;
    grown.max_entries = entries;
    if (!duckhts_bam_site_overlap_begin(&grown)) {
        free(grown.slots);
        return 0;
    }
    for (size_t at = 0; at < pending->held_count; at++) {
        const mate_held_t *entry = pending->heap[at];
        if (!duckhts_bam_site_overlap_insert(&grown, entry->qname, entry->qname_len, entry->hash)) {
            free(grown.slots);
            return 0;
        }
    }
    free(pending->table.slots);
    pending->table = grown;
    return 1;
}

/* A free entry, from the free list or a new block. NULL when out of memory. */
static mate_held_t *pending_acquire(mate_pending_t *pending) {
    size_t count = pending->max_pending - pending->entry_count;
    size_t total;
    mate_held_t *block;
    mate_held_t **grown;
    if (pending->free_count > 0) return pending->free_entries[--pending->free_count];
    if (count > MATE_ENTRY_BLOCK) count = MATE_ENTRY_BLOCK;
    total = pending->entry_count + count;
    /* Grow the lists first: a list longer than the entries is harmless. */
    grown = realloc(pending->free_entries, total * sizeof(*pending->free_entries));
    if (!grown) return NULL;
    pending->free_entries = grown;
    grown = realloc(pending->heap, total * sizeof(*pending->heap));
    if (!grown) return NULL;
    pending->heap = grown;
    grown = realloc(pending->blocks, (pending->block_count + 1) * sizeof(*pending->blocks));
    if (!grown) return NULL;
    pending->blocks = grown;
    block = calloc(count, sizeof(*block));
    if (!block) return NULL;
    pending->blocks[pending->block_count++] = block;
    pending->entry_count = total;
    for (size_t at = count; at > 1; at--) pending->free_entries[pending->free_count++] = &block[at - 1];
    return &block[0];
}

/* Hold the alignment in hand until its mate arrives. The caller has checked
 * that it is not held already, that the pending count is below max_pending
 * and that its CIGAR matches its sequence. */
static int pending_hold(mate_pending_t *pending, const bam1_t *alignment, const char *qname, size_t qname_len,
                        uint64_t hash) {
    const size_t cigar_bytes = (size_t)alignment->core.n_cigar * sizeof(uint32_t);
    const size_t sequence_bytes = ((size_t)alignment->core.l_qseq + 1) / 2;
    const size_t quality_bytes = (size_t)alignment->core.l_qseq;
    const size_t bytes = cigar_bytes + sequence_bytes + quality_bytes;
    mate_held_t *entry;

    if (pending->held_count == pending->table.max_entries && !pending_grow_table(pending)) return 0;
    entry = pending_acquire(pending);
    if (!entry) return 0;
    if (bytes > entry->data_capacity) {
        uint8_t *data = realloc(entry->data, bytes);
        if (!data) {
            pending_free(pending, entry);
            return 0;
        }
        entry->data = data;
        entry->data_capacity = bytes;
    }
    memcpy(entry->data, bam_get_cigar(alignment), cigar_bytes);
    memcpy(entry->data + cigar_bytes, bam_get_seq(alignment), sequence_bytes);
    memcpy(entry->data + cigar_bytes + sequence_bytes, bam_get_qual(alignment), quality_bytes);
    memcpy(entry->qname, qname, qname_len);
    entry->qname[qname_len] = '\0';
    entry->qname_len = (uint32_t)qname_len;
    entry->hash = hash;
    entry->pos = alignment->core.pos;
    entry->mpos = alignment->core.mpos;
    entry->n_cigar = alignment->core.n_cigar;
    entry->l_qseq = alignment->core.l_qseq;
    entry->flag = alignment->core.flag;
    if (!duckhts_bam_site_overlap_insert(&pending->table, entry->qname, qname_len, hash)) {
        pending_free(pending, entry);
        return 0;
    }
    heap_place(pending, pending->held_count++, entry);
    heap_sift_up(pending, entry->heap_index);
    return 1;
}

/* ---- comparing two mates ---- */

static mate_view_t view_of_alignment(const bam1_t *alignment) {
    mate_view_t view;
    view.cigar = bam_get_cigar(alignment);
    view.sequence = bam_get_seq(alignment);
    view.qualities = bam_get_qual(alignment);
    view.pos = alignment->core.pos;
    view.n_cigar = alignment->core.n_cigar;
    view.l_qseq = alignment->core.l_qseq;
    view.has_qualities = view.l_qseq > 0 && view.qualities[0] != 0xff;
    return view;
}

static mate_view_t view_of_held(const mate_held_t *entry) {
    const size_t cigar_bytes = (size_t)entry->n_cigar * sizeof(uint32_t);
    const size_t sequence_bytes = ((size_t)entry->l_qseq + 1) / 2;
    mate_view_t view;
    view.cigar = (const uint32_t *)(const void *)entry->data;
    view.sequence = entry->data + cigar_bytes;
    view.qualities = entry->data + cigar_bytes + sequence_bytes;
    view.pos = entry->pos;
    view.n_cigar = entry->n_cigar;
    view.l_qseq = entry->l_qseq;
    view.has_qualities = view.l_qseq > 0 && view.qualities[0] != 0xff;
    return view;
}

/* A read whose CIGAR and sequence disagree in length has no trusted base. */
static int alignment_has_bases(const bam1_t *alignment) {
    const uint32_t *cigar = bam_get_cigar(alignment);
    int64_t query = 0;
    if (alignment->core.l_qseq <= 0 || alignment->core.n_cigar == 0) return 0;
    for (uint32_t op = 0; op < alignment->core.n_cigar; op++) {
        if (bam_cigar_type(bam_cigar_op(cigar[op])) & 1) query += bam_cigar_oplen(cigar[op]);
    }
    return query == alignment->core.l_qseq;
}

static int blocks_next(mate_blocks_t *blocks) {
    while (blocks->next_op < blocks->view->n_cigar) {
        const uint32_t item = blocks->view->cigar[blocks->next_op++];
        const int type = bam_cigar_type(bam_cigar_op(item));
        const int64_t length = bam_cigar_oplen(item);
        const hts_pos_t reference = blocks->reference;
        const int64_t query = blocks->query;
        if (type & 1) blocks->query += length;
        if (type & 2) blocks->reference += length;
        if (type == 3 && length > 0) {
            blocks->block_ref = reference;
            blocks->block_query = query;
            blocks->block_length = length;
            return 1;
        }
    }
    return 0;
}

static void blocks_start(mate_blocks_t *blocks, const mate_view_t *view) {
    memset(blocks, 0, sizeof(*blocks));
    blocks->view = view;
    blocks->reference = view->pos;
}

/* The reference base at pos0 as 0..3, or -1 when it is not A, C, G or T or
 * there is no reference. */
static int reference_base(mate_reference_t *reference, int tid, const char *contig, hts_pos_t pos0, int *code,
                          char *error, size_t error_size) {
    *code = -1;
    if (!reference->fai) return 1;
    if (!reference->bases || reference->tid != tid || pos0 < reference->beg0 || pos0 >= reference->end0) {
        hts_pos_t end0;
        hts_pos_t fetched = 0;
        if (!reference->bases || reference->tid != tid) {
            reference->contig_length = faidx_seq_len64(reference->fai, contig);
            if (reference->contig_length < 0) {
                snprintf(error, error_size,
                         "duckhts_bam_mate_overlap_counts: reference has no contig named %s; "
                         "contig names are compared byte for byte",
                         contig);
                return 0;
            }
        }
        if (pos0 >= reference->contig_length) {
            snprintf(error, error_size,
                     "duckhts_bam_mate_overlap_counts: an alignment on %s reaches %lld, past the reference contig "
                     "length %lld; the reference does not match the alignment header",
                     contig, (long long)pos0 + 1, (long long)reference->contig_length);
            return 0;
        }
        free(reference->bases);
        reference->bases = NULL;
        reference->tid = tid;
        end0 = pos0 + MATE_WINDOW_BASES;
        if (end0 > reference->contig_length) end0 = reference->contig_length;
        reference->bases = faidx_fetch_seq64(reference->fai, contig, pos0, end0 - 1, &fetched);
        if (!reference->bases || fetched != end0 - pos0) {
            snprintf(error, error_size, "duckhts_bam_mate_overlap_counts: cannot fetch %s:%lld-%lld from the reference",
                     contig, (long long)pos0 + 1, (long long)end0);
            return 0;
        }
        reference->beg0 = pos0;
        reference->end0 = end0;
    }
    *code = base_code(reference->bases[pos0 - reference->beg0]);
    return 1;
}

/* The read base at query offset `query`, aligned to pos0, as 0..3 or -1. A '='
 * base is the reference base. */
static int read_base(const mate_view_t *view, int64_t query, hts_pos_t pos0, mate_reference_t *reference, int tid,
                     const char *contig, int *code, char *error, size_t error_size) {
    const int nt16 = bam_seqi(view->sequence, query);
    if (nt16 == 0) return reference_base(reference, tid, contig, pos0, code, error, error_size);
    *code = nt16_code(nt16);
    return 1;
}

static int quality_index(const mate_view_t *view, int64_t query) {
    int quality;
    if (!view->has_qualities) return MATE_QUALITY_MISSING;
    quality = view->qualities[query];
    return quality > MATE_QUALITY_MAX ? MATE_QUALITY_MAX : quality;
}

/* Count each reference position where both mates have an aligned base. The
 * aligned blocks of each mate come in reference order, so the walk merges the
 * two lists and keeps one block of each. */
static int count_overlap(const mate_view_t *read1, const mate_view_t *read2, mate_reference_t *reference, int tid,
                         const char *contig, uint64_t *cells, char *error, size_t error_size) {
    mate_blocks_t first, second;
    int have_first, have_second;
    blocks_start(&first, read1);
    blocks_start(&second, read2);
    have_first = blocks_next(&first);
    have_second = blocks_next(&second);
    while (have_first && have_second) {
        const hts_pos_t first_end = first.block_ref + first.block_length;
        const hts_pos_t second_end = second.block_ref + second.block_length;
        const hts_pos_t from = first.block_ref > second.block_ref ? first.block_ref : second.block_ref;
        const hts_pos_t to = first_end < second_end ? first_end : second_end;
        for (hts_pos_t pos0 = from; pos0 < to; pos0++) {
            const int64_t query1 = first.block_query + (pos0 - first.block_ref);
            const int64_t query2 = second.block_query + (pos0 - second.block_ref);
            int base1, base2;
            if (!read_base(read1, query1, pos0, reference, tid, contig, &base1, error, error_size)) return 0;
            if (base1 < 0) continue;
            if (!read_base(read2, query2, pos0, reference, tid, contig, &base2, error, error_size)) return 0;
            if (base2 < 0) continue;
            cells[mate_cell(quality_index(read1, query1), quality_index(read2, query2), base1, base2)]++;
        }
        if (first_end <= second_end) {
            have_first = blocks_next(&first);
        } else {
            have_second = blocks_next(&second);
        }
    }
    return 1;
}

/* A held mate and the alignment in hand form a pair when each names the
 * other's start and one is read 1 (0x40) and the other read 2 (0x80). */
static int mates_mirror(const mate_held_t *held, const bam1_t *alignment) {
    const uint16_t mates = BAM_FREAD1 | BAM_FREAD2;
    const uint16_t held_mate = held->flag & mates;
    const uint16_t alignment_mate = alignment->core.flag & mates;
    return held->mpos == alignment->core.pos && alignment->core.mpos == held->pos &&
           (held_mate == BAM_FREAD1 || held_mate == BAM_FREAD2) && alignment_mate == (mates ^ held_mate);
}

/* ---- the scan ---- */

static int region_name2id(void *header, const char *name) {
    return sam_hdr_name2tid((sam_hdr_t *)header, name);
}

/* The whole count: open the inputs, walk the alignments, release everything. */
static int count_alignments(const mate_bind_t *bind, uint64_t *cells, char *error, size_t error_size) {
    const uint16_t non_primary = BAM_FSECONDARY | BAM_FSUPPLEMENTARY;
    samFile *file = NULL;
    sam_hdr_t *header = NULL;
    hts_idx_t *index = NULL;
    hts_itr_t *iterator = NULL;
    bam1_t *alignment = NULL;
    char **regions = NULL;
    char *reference_locator = NULL;
    unsigned int region_count = 0;
    mate_reference_t reference;
    mate_pending_t pending;
    int pending_ready = 0;
    int last_tid = -1;
    hts_pos_t last_pos = -1;
    int read_status = 0;
    int ok = 0;

    memset(&reference, 0, sizeof(reference));
    reference.tid = -1;

    file = sam_open(bind->path, "r");
    if (!file) {
        snprintf(error, error_size, "duckhts_bam_mate_overlap_counts: cannot open %s", bind->path);
        goto cleanup;
    }
    if (bind->reference) {
        /* The CRAM decoder reads the FASTA index from the locator, so an index
         * that is not the adjacent .fai travels with the reference path. */
        reference_locator = duckhts_reference_locator(bind->reference, bind->reference_index_path);
        if (!reference_locator) {
            snprintf(error, error_size, "duckhts_bam_mate_overlap_counts: out of memory");
            goto cleanup;
        }
        if (hts_set_fai_filename(file, reference_locator) != 0) {
            snprintf(error, error_size, "duckhts_bam_mate_overlap_counts: cannot use reference %s", bind->reference);
            goto cleanup;
        }
    } else if (hts_get_format(file)->format == cram) {
        snprintf(error, error_size, "duckhts_bam_mate_overlap_counts: CRAM input %s needs reference", bind->path);
        goto cleanup;
    }
    header = sam_hdr_read(file);
    alignment = bam_init1();
    if (!header || !alignment) {
        snprintf(error, error_size, "duckhts_bam_mate_overlap_counts: cannot read the header of %s", bind->path);
        goto cleanup;
    }
    if (bind->reference) {
        reference.fai = fai_load3(bind->reference, bind->reference_index_path, NULL, 0);
        if (!reference.fai) {
            snprintf(error, error_size,
                     "duckhts_bam_mate_overlap_counts: cannot load the FASTA index of %s; build it with fasta_index()",
                     bind->reference);
            goto cleanup;
        }
    }

    if (bind->region) {
        if (!duckhts_region_list_parse(bind->region, &regions, &region_count, error, error_size) ||
            !duckhts_region_list_validate(regions, region_count, region_name2id, header, error, error_size)) {
            goto cleanup;
        }
        for (unsigned int item = 0; item < region_count; item++) {
            if (strcmp(regions[item], ".") == 0 || strcmp(regions[item], "*") == 0) {
                snprintf(error, error_size,
                         "duckhts_bam_mate_overlap_counts: region items must name a contig or an interval");
                goto cleanup;
            }
        }
        index = sam_index_load3(file, bind->path, bind->index_path, HTS_IDX_SILENT_FAIL);
        if (!index) {
            snprintf(error, error_size, "duckhts_bam_mate_overlap_counts: region needs an index for %s", bind->path);
            goto cleanup;
        }
        iterator = sam_itr_regarray(index, header, regions, region_count);
        if (!iterator) {
            snprintf(error, error_size, "duckhts_bam_mate_overlap_counts: no region of '%s' is in %s",
                     bind->region, bind->path);
            goto cleanup;
        }
    }

    if (!pending_init(&pending, (size_t)bind->max_pending)) {
        pending_destroy(&pending);
        snprintf(error, error_size, "duckhts_bam_mate_overlap_counts: out of memory");
        goto cleanup;
    }
    pending_ready = 1;

    for (;;) {
        const char *qname;
        const char *qname_end;
        size_t qname_len;
        uint64_t hash;
        mate_held_t *held;

        read_status = iterator ? sam_itr_next(file, iterator, alignment) : sam_read1(file, header, alignment);
        if (read_status < 0) break;
        if (alignment->core.tid < 0) continue;

        /* Mate tracking relies on coordinate order; say so rather than miss pairs. */
        if (alignment->core.tid < last_tid || (alignment->core.tid == last_tid && alignment->core.pos < last_pos)) {
            snprintf(error, error_size,
                     "duckhts_bam_mate_overlap_counts: alignments of %s are not in coordinate order at %s:%lld",
                     bind->path, sam_hdr_tid2name(header, alignment->core.tid),
                     (long long)alignment->core.pos + 1);
            goto cleanup;
        }
        if (alignment->core.tid != last_tid && !pending_clear(&pending)) {
            snprintf(error, error_size, "duckhts_bam_mate_overlap_counts: invalid mate table");
            goto cleanup;
        }
        last_tid = alignment->core.tid;
        last_pos = alignment->core.pos;
        pending_drop_before(&pending, alignment->core.pos);

        if (alignment->core.flag & BAM_FUNMAP) continue;
        if (alignment->core.flag & bind->exclude_flags) continue;
        if ((alignment->core.flag & bind->require_flags) != bind->require_flags) continue;
        if (alignment->core.qual < bind->min_mapq) continue;
        /* Only primary alignments of a pair with both mates on this contig pair up. */
        if (!(alignment->core.flag & BAM_FPAIRED) || (alignment->core.flag & (BAM_FMUNMAP | non_primary)) ||
            alignment->core.mtid != alignment->core.tid) {
            continue;
        }

        qname = bam_get_qname(alignment);
        /* htslib stores the NUL and limits a name to 254 characters; check both
         * rather than trust them, since the held buffer has a fixed size. */
        qname_end = memchr(qname, '\0', alignment->core.l_qname);
        if (!qname_end) continue;
        qname_len = (size_t)(qname_end - qname);
        if (qname_len >= MATE_QNAME_BYTES) continue;
        hash = duckhts_bam_site_qname_hash(qname, qname_len);

        held = pending_find(&pending, qname, qname_len, hash);
        if (held) {
            /* A record that has the name of a held mate but does not pair with
             * it is left out. The held mate stays held for its own mate, so a
             * malformed record does not cost a valid pair. */
            if (!mates_mirror(held, alignment)) continue;
            pending_unlink(&pending, held);
            if (alignment_has_bases(alignment)) {
                const mate_view_t held_view = view_of_held(held);
                const mate_view_t alignment_view = view_of_alignment(alignment);
                const int held_is_read1 = (held->flag & BAM_FREAD1) != 0;
                const char *contig = sam_hdr_tid2name(header, alignment->core.tid);
                if (!count_overlap(held_is_read1 ? &held_view : &alignment_view,
                                   held_is_read1 ? &alignment_view : &held_view, &reference,
                                   alignment->core.tid, contig, cells, error, error_size)) {
                    pending_free(&pending, held);
                    goto cleanup;
                }
            }
            pending_free(&pending, held);
            continue;
        }

        /* The first mate: hold it when its mate starts inside it. */
        if (alignment->core.mpos < alignment->core.pos || alignment->core.mpos >= bam_endpos(alignment) ||
            !alignment_has_bases(alignment)) {
            continue;
        }
        if (pending.held_count == pending.max_pending) {
            snprintf(error, error_size,
                     "duckhts_bam_mate_overlap_counts: more than max_pending (%d) mates wait for their mate at "
                     "%s:%lld; raise max_pending",
                     bind->max_pending, sam_hdr_tid2name(header, alignment->core.tid),
                     (long long)alignment->core.pos + 1);
            goto cleanup;
        }
        if (!pending_hold(&pending, alignment, qname, qname_len, hash)) {
            snprintf(error, error_size, "duckhts_bam_mate_overlap_counts: out of memory");
            goto cleanup;
        }
    }
    if (read_status < -1) {
        snprintf(error, error_size, "duckhts_bam_mate_overlap_counts: cannot read an alignment of %s", bind->path);
        goto cleanup;
    }
    ok = 1;

cleanup:
    if (pending_ready) pending_destroy(&pending);
    free(reference.bases);
    if (reference.fai) fai_destroy(reference.fai);
    free(regions);
    free(reference_locator);
    if (iterator) hts_itr_destroy(iterator);
    if (index) hts_idx_destroy(index);
    if (alignment) bam_destroy1(alignment);
    if (header) sam_hdr_destroy(header);
    if (file) sam_close(file);
    return ok;
}

/* ---- DuckDB adapter ---- */

static void destroy_mate_bind(void *data) {
    mate_bind_t *bind = (mate_bind_t *)data;
    if (!bind) return;
    if (bind->path) duckdb_free(bind->path);
    if (bind->region) duckdb_free(bind->region);
    if (bind->index_path) duckdb_free(bind->index_path);
    if (bind->reference) duckdb_free(bind->reference);
    if (bind->reference_index_path) duckdb_free(bind->reference_index_path);
    duckdb_free(bind);
}

static void destroy_mate_state(void *data) {
    mate_state_t *state = (mate_state_t *)data;
    if (!state) return;
    free(state->cells);
    duckdb_free(state);
}

static char *named_text(duckdb_bind_info info, const char *name) {
    duckdb_value value = duckdb_bind_get_named_parameter(info, name);
    char *text = NULL;
    if (value && !duckdb_is_null_value(value)) text = duckdb_get_varchar(value);
    if (value) duckdb_destroy_value(&value);
    if (text && text[0] == '\0') {
        duckdb_free(text);
        text = NULL;
    }
    return text;
}

/* A named integer in [low, high], or the default when it is absent or NULL. */
static int named_integer(duckdb_bind_info info, const char *name, int fallback, int low, int high, int *out) {
    duckdb_value value = duckdb_bind_get_named_parameter(info, name);
    int64_t number = fallback;
    if (value && !duckdb_is_null_value(value)) number = duckdb_get_int64(value);
    if (value) duckdb_destroy_value(&value);
    if (number < low || number > high) return 0;
    *out = (int)number;
    return 1;
}

static void mate_overlap_counts_bind(duckdb_bind_info info) {
    duckdb_value path_value = duckdb_bind_get_parameter(info, 0);
    mate_bind_t *bind = (mate_bind_t *)duckdb_malloc(sizeof(*bind));
    duckdb_logical_type tinyint_type, bigint_type, varchar_type;

    if (bind) {
        memset(bind, 0, sizeof(*bind));
        if (!duckdb_is_null_value(path_value)) bind->path = duckdb_get_varchar(path_value);
    }
    duckdb_destroy_value(&path_value);
    if (!bind) {
        duckdb_bind_set_error(info, "duckhts_bam_mate_overlap_counts: out of memory");
        return;
    }
    if (!bind->path || bind->path[0] == '\0') {
        destroy_mate_bind(bind);
        duckdb_bind_set_error(info, "duckhts_bam_mate_overlap_counts requires an alignment path");
        return;
    }
    bind->region = named_text(info, "region");
    bind->index_path = named_text(info, "index_path");
    bind->reference = named_text(info, "reference");
    bind->reference_index_path = named_text(info, "reference_index_path");
    if (bind->reference_index_path && !bind->reference) {
        destroy_mate_bind(bind);
        duckdb_bind_set_error(info, "duckhts_bam_mate_overlap_counts: reference_index_path needs reference");
        return;
    }
    if (!named_integer(info, "min_mapq", 20, 0, 255, &bind->min_mapq) ||
        !named_integer(info, "require_flags", 0, 0, 65535, &bind->require_flags) ||
        !named_integer(info, "exclude_flags", 3844, 0, 65535, &bind->exclude_flags) ||
        !named_integer(info, "max_pending", 1000000, 1, MATE_PENDING_MAX, &bind->max_pending)) {
        destroy_mate_bind(bind);
        duckdb_bind_set_error(info,
            "duckhts_bam_mate_overlap_counts: min_mapq must be in [0, 255], the flag masks in [0, 65535] "
            "and max_pending in [1, 100000000]");
        return;
    }

    tinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_TINYINT);
    bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_bind_add_result_column(info, "read1_quality", tinyint_type);
    duckdb_bind_add_result_column(info, "read2_quality", tinyint_type);
    duckdb_bind_add_result_column(info, "read1_base", varchar_type);
    duckdb_bind_add_result_column(info, "read2_base", varchar_type);
    duckdb_bind_add_result_column(info, "positions", bigint_type);
    duckdb_destroy_logical_type(&tinyint_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_bind_set_bind_data(info, bind, destroy_mate_bind);
}

static void mate_overlap_counts_init(duckdb_init_info info) {
    mate_state_t *state = (mate_state_t *)duckdb_malloc(sizeof(*state));
    if (state) {
        memset(state, 0, sizeof(*state));
        state->cells = calloc(MATE_CELLS, sizeof(uint64_t));
    }
    if (!state || !state->cells) {
        destroy_mate_state(state);
        duckdb_init_set_error(info, "duckhts_bam_mate_overlap_counts: out of memory");
        return;
    }
    duckdb_init_set_init_data(info, state, destroy_mate_state);
    duckdb_init_set_max_threads(info, 1);
}

static void set_quality(duckdb_vector vector, int8_t *data, idx_t row, int quality) {
    if (quality == MATE_QUALITY_MISSING) {
        duckdb_vector_ensure_validity_writable(vector);
        duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vector), row);
    } else {
        data[row] = (int8_t)quality;
    }
}

static void mate_overlap_counts_scan(duckdb_function_info info, duckdb_data_chunk output) {
    static const char base_text[MATE_BASES] = {'A', 'C', 'G', 'T'};
    const mate_bind_t *bind = (const mate_bind_t *)duckdb_function_get_bind_data(info);
    mate_state_t *state = (mate_state_t *)duckdb_function_get_init_data(info);
    duckdb_vector vectors[MATE_VECTOR_COLUMNS];
    int8_t *quality1_data, *quality2_data;
    int64_t *positions_data;
    idx_t rows = 0;
    const idx_t capacity = duckdb_vector_size();

    if (!state->counted) {
        char error[512];
        error[0] = '\0';
        if (!count_alignments(bind, state->cells, error, sizeof(error))) {
            duckdb_function_set_error(info, error[0] ? error : "duckhts_bam_mate_overlap_counts failed");
            return;
        }
        state->counted = 1;
    }

    for (idx_t column = 0; column < MATE_VECTOR_COLUMNS; column++) {
        vectors[column] = duckdb_data_chunk_get_vector(output, column);
    }
    quality1_data = (int8_t *)duckdb_vector_get_data(vectors[0]);
    quality2_data = (int8_t *)duckdb_vector_get_data(vectors[1]);
    positions_data = (int64_t *)duckdb_vector_get_data(vectors[4]);

    while (state->next_cell < MATE_CELLS && rows < capacity) {
        size_t cell = state->next_cell++;
        const uint64_t count = state->cells[cell];
        int base1, base2, quality1, quality2;
        if (count == 0) continue;
        base2 = (int)(cell % MATE_BASES);
        cell /= MATE_BASES;
        base1 = (int)(cell % MATE_BASES);
        cell /= MATE_BASES;
        quality2 = (int)(cell % MATE_QUALITIES);
        quality1 = (int)(cell / MATE_QUALITIES);
        set_quality(vectors[0], quality1_data, rows, quality1);
        set_quality(vectors[1], quality2_data, rows, quality2);
        duckdb_vector_assign_string_element_len(vectors[2], rows, &base_text[base1], 1);
        duckdb_vector_assign_string_element_len(vectors[3], rows, &base_text[base2], 1);
        positions_data[rows] = count > (uint64_t)INT64_MAX ? INT64_MAX : (int64_t)count;
        rows++;
    }
    duckdb_data_chunk_set_size(output, rows);
}

void register_duckhts_bam_mate_overlap_counts_function(duckdb_connection connection) {
    duckdb_table_function function = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);

    duckdb_table_function_set_name(function, "duckhts_bam_mate_overlap_counts");
    duckdb_table_function_add_parameter(function, varchar_type);
    duckdb_table_function_add_named_parameter(function, "region", varchar_type);
    duckdb_table_function_add_named_parameter(function, "index_path", varchar_type);
    duckdb_table_function_add_named_parameter(function, "reference", varchar_type);
    duckdb_table_function_add_named_parameter(function, "reference_index_path", varchar_type);
    duckdb_table_function_add_named_parameter(function, "min_mapq", bigint_type);
    duckdb_table_function_add_named_parameter(function, "require_flags", bigint_type);
    duckdb_table_function_add_named_parameter(function, "exclude_flags", bigint_type);
    duckdb_table_function_add_named_parameter(function, "max_pending", bigint_type);
    duckdb_table_function_set_bind(function, mate_overlap_counts_bind);
    duckdb_table_function_set_init(function, mate_overlap_counts_init);
    duckdb_table_function_set_function(function, mate_overlap_counts_scan);
    duckdb_register_table_function(connection, function);
    duckdb_destroy_table_function(&function);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
}
