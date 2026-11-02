/*
 * Native test of the coordinate partitions in src/bam_partition.c.
 *
 * Writes a coordinate-sorted BAM (one 2 Mb reference with 200,000 reads, every
 * 50th read spliced across 200 kb so that reads cross piece boundaries, every
 * 97th read placed but unmapped, an empty reference, and reads without
 * coordinates), indexes it as BAI and as CSI, and checks that:
 *   - the pieces of each reference tile [0, HTS_POS_MAX) in order;
 *   - a reference with more indexed bytes than piece_bytes is split, and
 *     aligned boundaries are multiples of the alignment;
 *   - reading every piece with its own iterator and the ownership rule keeps
 *     each placed record exactly once;
 *   - piece_bytes 0 gives one piece per reference;
 *   - the ownership rule keeps a POS-0 record in a reference's first piece.
 */

#include "bam_partition.h"

#include <htslib/sam.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define N_PLACED 200000
#define N_UNPLACED 300
#define READ_LENGTH 100
#define REF_A_LENGTH 2000000
#define PIECE_BYTES (64 * 1024)

static int fail_check(const char *file, int line, const char *expression) {
    fprintf(stderr, "%s:%d: check failed: %s\n", file, line, expression);
    return 1;
}

#define CHECK(expression) \
    do { \
        if (!(expression)) return fail_check(__FILE__, __LINE__, #expression); \
    } while (0)

/* Deterministic bases so the BAM compresses like sequencing data, not like a repeat. */
static uint32_t next_random(uint32_t *state) {
    *state = *state * 1664525u + 1013904223u;
    return *state >> 8;
}

static int write_record(samFile *out, sam_hdr_t *hdr, bam1_t *rec, uint32_t *state,
                        long index, int tid, hts_pos_t pos, int spliced, int unmapped) {
    char qname[32];
    char sequence[READ_LENGTH];
    char quality[READ_LENGTH];
    uint32_t cigar[3];
    size_t n_cigar = 1;
    uint16_t flag = unmapped ? BAM_FUNMAP : 0;

    snprintf(qname, sizeof(qname), "r%07ld", index);
    for (int i = 0; i < READ_LENGTH; i++) {
        sequence[i] = "ACGT"[next_random(state) & 3u];
        quality[i] = (char)(20 + next_random(state) % 20u);
    }
    if (unmapped) {
        n_cigar = 0;
    } else if (spliced) {
        cigar[0] = bam_cigar_gen(READ_LENGTH / 2, BAM_CMATCH);
        cigar[1] = bam_cigar_gen(200000, BAM_CREF_SKIP);
        cigar[2] = bam_cigar_gen(READ_LENGTH / 2, BAM_CMATCH);
        n_cigar = 3;
    } else {
        cigar[0] = bam_cigar_gen(READ_LENGTH, BAM_CMATCH);
    }
    if (bam_set1(rec, strlen(qname), qname, flag, tid, pos, tid < 0 ? 0 : 60,
                 n_cigar, cigar, -1, -1, 0, READ_LENGTH, sequence, quality, 0) < 0) {
        return -1;
    }
    return sam_write1(out, hdr, rec) < 0 ? -1 : 0;
}

/* chrA holds the reads at positions 0 to 1.8 Mb; ref_a_length sets the length the
 * header declares, so a long one leaves the reads in a short prefix. */
static int write_bam(const char *path, long ref_a_length) {
    char header_text[160];
    int header_length = snprintf(header_text, sizeof(header_text),
                                 "@HD\tVN:1.6\tSO:coordinate\n"
                                 "@SQ\tSN:chrA\tLN:%ld\n"
                                 "@SQ\tSN:chrEmpty\tLN:1000\n", ref_a_length);
    samFile *out = sam_open(path, "wb");
    sam_hdr_t *hdr = header_length > 0 && (size_t)header_length < sizeof(header_text)
        ? sam_hdr_parse((size_t)header_length, header_text) : NULL;
    bam1_t *rec = bam_init1();
    uint32_t state = 12345u;
    int rc = -1;

    if (!out || !hdr || !rec || sam_hdr_write(out, hdr) < 0) goto done;
    for (long i = 0; i < N_PLACED; i++) {
        hts_pos_t pos = (hts_pos_t)(i * 9);
        if (write_record(out, hdr, rec, &state, i, 0, pos, i % 50 == 0, i % 97 == 0) < 0) goto done;
    }
    for (long i = 0; i < N_UNPLACED; i++) {
        if (write_record(out, hdr, rec, &state, N_PLACED + i, -1, -1, 0, 1) < 0) goto done;
    }
    rc = 0;
done:
    if (rec) bam_destroy1(rec);
    if (hdr) sam_hdr_destroy(hdr);
    if (out && sam_close(out) < 0) rc = -1;
    return rc;
}

/* Reads every coordinate piece with its own iterator and the ownership rule,
 * and checks that each placed record is kept exactly once. */
static int check_pieces_keep_each_record_once(const char *path, hts_idx_t *idx,
                                              const duckhts_bam_partitions *parts) {
    samFile *in = sam_open(path, "r");
    sam_hdr_t *hdr = in ? sam_hdr_read(in) : NULL;
    bam1_t *rec = bam_init1();
    unsigned char *seen = calloc(N_PLACED, 1);
    long kept = 0;
    long largest_piece = 0;
    int pieces_a = 0;

    CHECK(in && hdr && rec && seen);
    for (int p = 0; p < parts->n; p++) {
        hts_itr_t *itr = sam_itr_queryi(idx, parts->tid[p], parts->beg[p], parts->end[p]);
        long kept_before = kept;
        int ret;
        if (parts->tid[p] == 0) pieces_a++;
        if (!itr) continue; /* the empty reference has no chunk */
        while ((ret = sam_itr_next(in, itr, rec)) >= 0) {
            long index;
            if (!duckhts_bam_partition_owns(rec, parts->beg[p])) continue;
            CHECK(rec->core.pos < parts->end[p]);
            index = strtol(bam_get_qname(rec) + 1, NULL, 10);
            CHECK(index >= 0 && index < N_PLACED);
            CHECK(!seen[index]);
            seen[index] = 1;
            kept++;
        }
        CHECK(ret == -1);
        hts_itr_destroy(itr);
        if (kept - kept_before > largest_piece) largest_piece = kept - kept_before;
    }
    CHECK(kept == N_PLACED);
    /* Balance: once chrA is split into many pieces, none holds a quarter of its
     * reads, wherever on the reference the reads lie. A boundary shifts right by at
     * most the span of the reads overlapping it (200 kb for the spliced reads). */
    if (pieces_a >= 8) CHECK(largest_piece * 4 <= N_PLACED);
    free(seen);
    bam_destroy1(rec);
    sam_hdr_destroy(hdr);
    sam_close(in);
    return 0;
}

static void remove_indexes(const char *path) {
    char index_path[512];
    snprintf(index_path, sizeof(index_path), "%s.bai", path);
    unlink(index_path);
    snprintf(index_path, sizeof(index_path), "%s.csi", path);
    unlink(index_path);
}

static int check_index(const char *path, int min_shift, int align) {
    samFile *in = NULL;
    sam_hdr_t *hdr = NULL;
    hts_idx_t *idx = NULL;
    duckhts_bam_partitions parts = {0};
    int pieces_a = 0;

    CHECK(sam_index_build(path, min_shift) == 0);
    in = sam_open(path, "r");
    hdr = in ? sam_hdr_read(in) : NULL;
    idx = in ? sam_index_load(in, path) : NULL;
    CHECK(in && hdr && idx);

    /* Without a piece size, one piece per reference. */
    CHECK(duckhts_bam_partitions_build(idx, hdr, 0, 256, 0, &parts));
    CHECK(parts.n == 2);
    CHECK(parts.tid[0] == 0 && parts.beg[0] == 0 && parts.end[0] == HTS_POS_MAX);
    CHECK(parts.tid[1] == 1 && parts.beg[1] == 0 && parts.end[1] == HTS_POS_MAX);

    CHECK(duckhts_bam_partitions_build(idx, hdr, PIECE_BYTES, 256, align, &parts));
    for (int p = 0; p < parts.n; p++) {
        int first_of_reference = p == 0 || parts.tid[p] != parts.tid[p - 1];
        int last_of_reference = p == parts.n - 1 || parts.tid[p + 1] != parts.tid[p];
        CHECK(parts.beg[p] < parts.end[p]);
        CHECK(first_of_reference ? parts.beg[p] == 0 : parts.beg[p] == parts.end[p - 1]);
        if (last_of_reference) CHECK(parts.end[p] == HTS_POS_MAX);
        if (!first_of_reference && align > 1) CHECK(parts.beg[p] % align == 0);
        if (!first_of_reference) CHECK(parts.tid[p] >= parts.tid[p - 1]);
        if (parts.tid[p] == 0) pieces_a++;
    }
    CHECK(parts.tid[parts.n - 1] == 1); /* the empty reference stays one piece */
    CHECK(pieces_a > 4);
    if (check_pieces_keep_each_record_once(path, idx, &parts)) return 1;

    duckhts_bam_partitions_free(&parts);
    CHECK(parts.n == 0 && parts.beg == NULL);
    hts_idx_destroy(idx);
    sam_hdr_destroy(hdr);
    sam_close(in);
    remove_indexes(path); /* the next check must load the index it builds */
    return 0;
}

static int check_ownership_rule(void) {
    bam1_t *rec = bam_init1();
    CHECK(rec);
    rec->core.tid = 0;
    rec->core.pos = -1; /* a reference and POS 0 */
    CHECK(duckhts_bam_partition_owns(rec, 0));
    CHECK(!duckhts_bam_partition_owns(rec, 10));
    rec->core.pos = 9;
    CHECK(!duckhts_bam_partition_owns(rec, 10));
    rec->core.pos = 10;
    CHECK(duckhts_bam_partition_owns(rec, 10));
    rec->core.tid = -1; /* no reference */
    rec->core.pos = -1;
    CHECK(duckhts_bam_partition_owns(rec, 10));
    bam_destroy1(rec);
    return 0;
}

int main(void) {
    char path[] = "/tmp/duckhts-bam-partition-XXXXXX";
    char bam_path[sizeof(path) + 8];
    int fd = mkstemp(path);
    int rc = 1;

    if (fd < 0) return fail_check(__FILE__, __LINE__, "mkstemp");
    close(fd);
    snprintf(bam_path, sizeof(bam_path), "%s.bam", path);

    if (check_ownership_rule() == 0 && write_bam(bam_path, REF_A_LENGTH) == 0 &&
        check_index(bam_path, 0, 0) == 0 &&    /* BAI */
        check_index(bam_path, 14, 0) == 0 &&   /* CSI, default min_shift */
        check_index(bam_path, 0, 10000) == 0 && /* BAI, bin-aligned boundaries */
        /* The same reads in the first 0.4% of a 500 Mb reference: boundaries must
         * follow the data, not the declared length. */
        write_bam(bam_path, 500000000L) == 0 &&
        check_index(bam_path, 0, 0) == 0) {
        rc = 0;
    }

    unlink(path);
    unlink(bam_path);
    remove_indexes(bam_path);
    if (rc == 0) printf("bam_partition_test: OK\n");
    return rc;
}
