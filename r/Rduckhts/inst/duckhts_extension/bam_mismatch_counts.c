/* duckhts_bam_mismatch_counts: aligned read bases compared with the reference,
 * counted by mate, cycle, base quality and substitution.
 *
 * One pass over the alignments, one worker. Live state is fixed: the counter
 * cells, one reference window and its mask. Nothing grows with the reads. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include "include/hts_io_tuning.h"
#include "include/region_list.h"

enum {
    MISMATCH_MATES = 3,            /* 0 unpaired, 1 first mate, 2 second mate */
    MISMATCH_CYCLE_MAX = 1000,     /* a later cycle is counted at this cycle */
    MISMATCH_QUALITY_MAX = 93,     /* a higher base quality is counted at this quality */
    MISMATCH_QUALITY_MISSING = 94, /* the read stores no base qualities */
    MISMATCH_QUALITIES = 95,
    MISMATCH_BASES = 4,            /* A, C, G, T */
    MISMATCH_FLANK_MAX = 1000,
    MISMATCH_VECTOR_COLUMNS = 6
};

#define MISMATCH_CELLS \
    ((size_t)MISMATCH_MATES * (MISMATCH_CYCLE_MAX + 1) * MISMATCH_QUALITIES * MISMATCH_BASES * MISMATCH_BASES)
/* Reference bases held at once. The window moves along an alignment that is
 * longer, so this is also the most reference one alignment needs in memory. */
#define MISMATCH_WINDOW_BASES ((hts_pos_t)1 << 20)

_Static_assert(MISMATCH_CYCLE_MAX + 1 <= INT32_MAX && MISMATCH_QUALITY_MISSING == MISMATCH_QUALITY_MAX + 1 &&
               MISMATCH_QUALITIES == MISMATCH_QUALITY_MISSING + 1, "counter axes");

typedef struct {
    char *path;
    char *reference;
    char *region;
    char *mask;
    char *index_path;
    char *reference_index_path;
    char *mask_index_path;
    int min_mapq;
    int require_flags;
    int exclude_flags;
    int indel_flank;
} mismatch_bind_t;

typedef struct {
    uint64_t *cells; /* MISMATCH_CELLS counters, owned (malloc family) */
    size_t next_cell;
    int counted;
} mismatch_state_t;

/* Known variants: a BCF, or a bgzip-compressed VCF with a tabix index. */
typedef struct {
    htsFile *file;
    bcf_hdr_t *header;
    hts_idx_t *bcf_index; /* BCF only */
    tbx_t *tabix;         /* VCF only */
    bcf1_t *record;
    kstring_t line;
} mismatch_mask_t;

/* One reference window [beg0, end0) of contig tid, with one mask byte per base.
 * It holds at most MISMATCH_WINDOW_BASES bases. */
typedef struct {
    faidx_t *fai;
    int tid;
    hts_pos_t beg0;
    hts_pos_t end0;
    hts_pos_t contig_length; /* of contig tid in the reference */
    char *bases;     /* owned by htslib's allocator: free() */
    uint8_t *masked; /* owned (malloc family), MISMATCH_WINDOW_BASES bytes */
} mismatch_window_t;

static size_t mismatch_cell(int mate, int cycle, int quality, int from, int to) {
    return ((((size_t)mate * (MISMATCH_CYCLE_MAX + 1) + (size_t)cycle) * MISMATCH_QUALITIES +
             (size_t)quality) * MISMATCH_BASES + (size_t)from) * MISMATCH_BASES + (size_t)to;
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

static void mask_close(mismatch_mask_t *mask) {
    if (mask->record) bcf_destroy(mask->record);
    if (mask->tabix) tbx_destroy(mask->tabix);
    if (mask->bcf_index) hts_idx_destroy(mask->bcf_index);
    if (mask->header) bcf_hdr_destroy(mask->header);
    if (mask->file) hts_close(mask->file);
    ks_free(&mask->line);
    memset(mask, 0, sizeof(*mask));
}

static int mask_open(mismatch_mask_t *mask, const char *path, const char *index_path,
                     char *error, size_t error_size) {
    const htsFormat *format;
    memset(mask, 0, sizeof(*mask));
    mask->file = hts_open(path, "r");
    if (!mask->file) {
        snprintf(error, error_size, "duckhts_bam_mismatch_counts: cannot open mask %s", path);
        return 0;
    }
    format = hts_get_format(mask->file);
    if (format->format != bcf && !(format->format == vcf && format->compression == bgzf)) {
        snprintf(error, error_size,
                 "duckhts_bam_mismatch_counts: mask must be a BCF or a bgzip-compressed VCF: %s", path);
        return 0;
    }
    mask->header = bcf_hdr_read(mask->file);
    mask->record = bcf_init();
    if (!mask->header || !mask->record) {
        snprintf(error, error_size, "duckhts_bam_mismatch_counts: cannot read the header of mask %s", path);
        return 0;
    }
    if (format->format == bcf) {
        mask->bcf_index = bcf_index_load3(path, index_path, HTS_IDX_SILENT_FAIL);
    } else {
        mask->tabix = tbx_index_load3(path, index_path, HTS_IDX_SILENT_FAIL);
    }
    if (!mask->bcf_index && !mask->tabix) {
        snprintf(error, error_size, "duckhts_bam_mismatch_counts: mask needs an index: %s", path);
        return 0;
    }
    return 1;
}

/* Mark [beg0, end0) of the window, clipped to it. */
static void window_mark(mismatch_window_t *window, hts_pos_t beg0, hts_pos_t end0) {
    if (beg0 < window->beg0) beg0 = window->beg0;
    if (end0 > window->end0) end0 = window->end0;
    if (beg0 < end0) memset(window->masked + (size_t)(beg0 - window->beg0), 1, (size_t)(end0 - beg0));
}

/* A record whose alleles all have the length of REF masks its own positions.
 * A record with an indel or a symbolic allele also masks flank bases on each side. */
static void window_mark_variant(mismatch_window_t *window, hts_pos_t pos0, hts_pos_t reference_length,
                                int same_length, int flank) {
    hts_pos_t end0;
    if (reference_length < 1) reference_length = 1;
    /* A mark is clipped to the window, so an end past the window is the
     * window's end. A stated end (INFO/END) can be any number; stopping it
     * here keeps the sums below in range. */
    end0 = (pos0 >= window->end0 || reference_length > window->end0 - pos0) ? window->end0
                                                                            : pos0 + reference_length;
    if (same_length) {
        window_mark(window, pos0, end0);
    } else {
        window_mark(window, pos0 - flank, end0 + flank);
    }
}

/* REF, ALT and INFO/END of one VCF text line. Returns 0 when the line has no
 * POS, REF and ALT fields. */
static int mask_parse_vcf_line(const char *line, hts_pos_t *pos0, hts_pos_t *reference_length,
                               int *same_length) {
    const char *field[9] = {0};
    size_t length[9] = {0};
    int fields = 0;
    const char *cursor = line;
    char *number_end = NULL;
    long long pos1;
    while (fields < 9) {
        const char *tab = strchr(cursor, '\t');
        field[fields] = cursor;
        length[fields] = tab ? (size_t)(tab - cursor) : strlen(cursor);
        fields++;
        if (!tab) break;
        cursor = tab + 1;
    }
    if (fields < 5) return 0;
    pos1 = strtoll(field[1], &number_end, 10);
    if (number_end == field[1] || pos1 < 1) return 0;
    *pos0 = (hts_pos_t)pos1 - 1;
    *reference_length = (hts_pos_t)length[3];
    *same_length = 1;
    for (size_t at = 0, start = 0; at <= length[4]; at++) {
        if (at < length[4] && field[4][at] != ',') continue;
        const char *allele = field[4] + start;
        size_t allele_length = at - start;
        if (!(allele_length == 1 && allele[0] == '.') &&
            (allele_length != length[3] || allele[0] == '<' || allele[0] == '*' ||
             memchr(allele, '[', allele_length) || memchr(allele, ']', allele_length))) {
            *same_length = 0;
        }
        start = at + 1;
    }
    if (fields >= 8) {
        /* A symbolic allele states its end in INFO/END, one-based and inclusive. */
        const char *info = field[7];
        size_t at = 0;
        while (at < length[7]) {
            size_t stop = at;
            while (stop < length[7] && info[stop] != ';') stop++;
            if (stop - at > 4 && strncmp(info + at, "END=", 4) == 0) {
                long long end1 = strtoll(info + at + 4, &number_end, 10);
                if (number_end != info + at + 4 && end1 >= pos1) {
                    *reference_length = (hts_pos_t)(end1 - pos1 + 1);
                }
            }
            at = stop + 1;
        }
    }
    return 1;
}

/* Fill the mask of the loaded window. A contig that the mask does not know is
 * an error: names are compared byte for byte, and a silent miss would count
 * known variants as errors. */
static int window_fill_mask(mismatch_window_t *window, mismatch_mask_t *mask, const char *contig,
                            int flank, char *error, size_t error_size) {
    hts_pos_t beg0 = window->beg0 > flank ? window->beg0 - flank : 0;
    hts_pos_t end0 = window->end0 + flank;
    hts_itr_t *iterator = NULL;
    int status = 0;
    int header_id = bcf_hdr_name2id(mask->header, contig);
    if (mask->bcf_index) {
        if (header_id < 0) goto unknown_contig;
        iterator = bcf_itr_queryi(mask->bcf_index, header_id, beg0, end0);
        if (!iterator) goto failed;
        while ((status = bcf_itr_next(mask->file, iterator, mask->record)) >= 0) {
            int same_length = 1;
            if (bcf_unpack(mask->record, BCF_UN_STR) != 0) goto failed;
            for (int allele = 1; allele < mask->record->n_allele; allele++) {
                const char *text = mask->record->d.allele[allele];
                if (strlen(text) != strlen(mask->record->d.allele[0]) || text[0] == '<' || text[0] == '*' ||
                    strchr(text, '[') || strchr(text, ']')) {
                    same_length = 0;
                }
            }
            window_mark_variant(window, mask->record->pos, mask->record->rlen, same_length, flank);
        }
    } else {
        int tabix_id = tbx_name2id(mask->tabix, contig);
        if (tabix_id < 0) {
            if (header_id < 0) goto unknown_contig;
            return 1; /* a header contig without records */
        }
        iterator = tbx_itr_queryi(mask->tabix, tabix_id, beg0, end0);
        if (!iterator) goto failed;
        while ((status = tbx_itr_next(mask->file, mask->tabix, iterator, &mask->line)) >= 0) {
            hts_pos_t pos0 = 0, reference_length = 0;
            int same_length = 1;
            if (!mask_parse_vcf_line(mask->line.s, &pos0, &reference_length, &same_length)) goto failed;
            window_mark_variant(window, pos0, reference_length, same_length, flank);
        }
    }
    hts_itr_destroy(iterator);
    if (status < -1) {
        snprintf(error, error_size, "duckhts_bam_mismatch_counts: cannot read mask records of %s", contig);
        return 0;
    }
    return 1;

unknown_contig:
    snprintf(error, error_size,
             "duckhts_bam_mismatch_counts: mask has no contig named %s; contig names are compared byte for byte",
             contig);
    return 0;
failed:
    if (iterator) hts_itr_destroy(iterator);
    snprintf(error, error_size, "duckhts_bam_mismatch_counts: cannot read mask records of %s", contig);
    return 0;
}

/* Load the window that starts at beg0 on contig tid: MISMATCH_WINDOW_BASES
 * bases, or fewer at the end of the contig, with their mask. beg0 is on the
 * contig. */
static int window_fetch(mismatch_window_t *window, mismatch_mask_t *mask, int tid, const char *contig,
                        hts_pos_t beg0, int flank, char *error, size_t error_size) {
    hts_pos_t end0 = beg0 + MISMATCH_WINDOW_BASES;
    hts_pos_t fetched = 0;
    if (end0 > window->contig_length) end0 = window->contig_length;
    if (!window->masked) {
        window->masked = malloc((size_t)MISMATCH_WINDOW_BASES);
        if (!window->masked) {
            snprintf(error, error_size, "duckhts_bam_mismatch_counts: out of memory");
            return 0;
        }
    }
    free(window->bases);
    window->bases = NULL;
    window->tid = tid;
    window->beg0 = beg0;
    window->end0 = beg0; /* empty until the fetch succeeds */
    window->bases = faidx_fetch_seq64(window->fai, contig, beg0, end0 - 1, &fetched);
    if (!window->bases || fetched != end0 - beg0) {
        snprintf(error, error_size, "duckhts_bam_mismatch_counts: cannot fetch %s:%lld-%lld from the reference",
                 contig, (long long)beg0 + 1, (long long)end0);
        return 0;
    }
    window->end0 = end0;
    memset(window->masked, 0, (size_t)fetched);
    if (mask->file) return window_fill_mask(window, mask, contig, flank, error, error_size);
    return 1;
}

/* Check that the alignment [beg0, end0) lies on the reference contig and make
 * the window hold its start. Alignments arrive in coordinate order, so a new
 * window starts at the alignment and runs ahead. An alignment that reaches
 * past the window moves it during the walk. */
static int window_cover(mismatch_window_t *window, mismatch_mask_t *mask, int tid, const char *contig,
                        hts_pos_t beg0, hts_pos_t end0, int flank, char *error, size_t error_size) {
    if (window->bases && window->tid == tid && beg0 >= window->beg0 && end0 <= window->end0) return 1;
    if (window->tid != tid || !window->bases) {
        window->contig_length = faidx_seq_len64(window->fai, contig);
        if (window->contig_length < 0) {
            snprintf(error, error_size,
                     "duckhts_bam_mismatch_counts: reference has no contig named %s; "
                     "contig names are compared byte for byte",
                     contig);
            return 0;
        }
    }
    /* A contig the reference has but shorter than the alignment says is a
     * reference that does not match the alignment header; skipping the bases
     * past its end would give plausible but incomplete counts. */
    if (end0 > window->contig_length) {
        snprintf(error, error_size,
                 "duckhts_bam_mismatch_counts: an alignment on %s ends at %lld, past the reference contig length %lld; "
                 "the reference does not match the alignment header",
                 contig, (long long)end0, (long long)window->contig_length);
        return 0;
    }
    return window_fetch(window, mask, tid, contig, beg0, flank, error, error_size);
}

/* A CIGAR operation next to which bases are left out. */
static int cigar_is_gap(int code) {
    return code == BAM_CINS || code == BAM_CSOFT_CLIP || code == BAM_CDEL || code == BAM_CREF_SKIP;
}

/* Count the aligned bases of one alignment. Bases within flank query bases of
 * an insertion, a deletion, a reference skip or a soft clip are left out, with
 * masked positions and bases that are not A, C, G or T on either side.
 *
 * Gaps come in query order, so only two of them decide whether a base is left
 * out: the last gap behind it and the next gap ahead of it. The walk keeps
 * those two bounds and no list of gaps, so its state does not grow with the
 * number of CIGAR operations. */
static int count_alignment(const bam1_t *alignment, const char *contig, int flank,
                           mismatch_window_t *window, mismatch_mask_t *mask,
                           uint64_t *cells, char *error, size_t error_size) {
    const uint32_t *cigar = bam_get_cigar(alignment);
    const uint8_t *sequence = bam_get_seq(alignment);
    const uint8_t *qualities = bam_get_qual(alignment);
    const int query_length = alignment->core.l_qseq;
    const int reverse = (alignment->core.flag & BAM_FREVERSE) != 0;
    const int has_qualities = query_length > 0 && qualities[0] != 0xff;
    const hts_pos_t beg0 = alignment->core.pos;
    const hts_pos_t end0 = bam_endpos(alignment);
    int mate = 0;
    const uint32_t operations = alignment->core.n_cigar;
    int64_t hard_left = 0, hard_right = 0, query = 0;
    hts_pos_t reference = beg0;
    /* Bases are left out before left_out_until (the end of the last gap behind,
     * plus the flank) and from left_out_from on (the start of the next gap
     * ahead, minus the flank). next_gap is the operation that set
     * left_out_from, or `operations` when no gap is ahead. */
    int64_t left_out_until = 0;
    int64_t left_out_from = INT64_MAX;
    uint32_t next_gap = 0;

    if (query_length <= 0 || alignment->core.n_cigar == 0) return 1;
    if (alignment->core.flag & BAM_FPAIRED) {
        mate = (alignment->core.flag & BAM_FREAD1) ? 1 : (alignment->core.flag & BAM_FREAD2) ? 2 : 0;
    }
    if (!window_cover(window, mask, alignment->core.tid, contig, beg0, end0, flank, error, error_size)) return 0;

    for (uint32_t op = 0; op < operations; op++) {
        const int code = bam_cigar_op(cigar[op]);
        const int64_t length = bam_cigar_oplen(cigar[op]);
        if (code == BAM_CHARD_CLIP) {
            if (op == 0) hard_left = length; else hard_right = length;
        }
        if (bam_cigar_type(code) & 1) query += length;
    }
    if (query != query_length) return 1; /* CIGAR and sequence disagree: no base is trusted */

    query = 0;
    for (uint32_t op = 0; op < operations; op++) {
        const int code = bam_cigar_op(cigar[op]);
        const int type = bam_cigar_type(code);
        const int64_t length = bam_cigar_oplen(cigar[op]);
        if (flank > 0 && next_gap <= op) {
            /* The next gap after this operation. The search starts where the
             * last one ended, so every operation is looked at once. */
            int64_t ahead = query + ((type & 1) ? length : 0);
            next_gap = op + 1;
            while (next_gap < operations && !cigar_is_gap(bam_cigar_op(cigar[next_gap]))) {
                if (bam_cigar_type(bam_cigar_op(cigar[next_gap])) & 1) ahead += bam_cigar_oplen(cigar[next_gap]);
                next_gap++;
            }
            left_out_from = next_gap < operations ? ahead - flank : INT64_MAX;
        }
        if (flank > 0 && cigar_is_gap(code)) {
            /* An insertion or a soft clip holds query bases; a deletion or a
             * reference skip sits between two of them. */
            left_out_until = query + ((type & 1) ? length : 0) + flank;
        }
        if (type == 3) {
            for (int64_t step = 0; step < length; step++) {
                const int64_t at = query + step;
                const hts_pos_t pos0 = reference + step;
                int quality, from, to;
                int64_t cycle;
                if (at < left_out_until || at >= left_out_from) continue;
                if (pos0 < window->beg0 || pos0 >= window->end0) {
                    /* The alignment reaches past the window: move the window to this base. */
                    if (!window_fetch(window, mask, alignment->core.tid, contig, pos0, flank, error, error_size)) {
                        return 0;
                    }
                }
                if (window->masked[pos0 - window->beg0]) continue;
                from = base_code(window->bases[pos0 - window->beg0]);
                to = nt16_code(bam_seqi(sequence, at));
                if (from < 0 || to < 0) continue;
                if (reverse) {
                    from = 3 - from;
                    to = 3 - to;
                }
                quality = has_qualities ? qualities[at] : MISMATCH_QUALITY_MISSING;
                if (has_qualities && quality > MISMATCH_QUALITY_MAX) quality = MISMATCH_QUALITY_MAX;
                cycle = reverse ? hard_right + (query_length - at) : hard_left + at + 1;
                if (cycle > MISMATCH_CYCLE_MAX) cycle = MISMATCH_CYCLE_MAX;
                cells[mismatch_cell(mate, (int)cycle, quality, from, to)]++;
            }
        }
        if (type & 1) query += length;
        if (type & 2) reference += length;
    }
    return 1;
}

static int region_name2id(void *header, const char *name) {
    return sam_hdr_name2tid((sam_hdr_t *)header, name);
}

/* The whole count: open the inputs, walk the alignments, release everything. */
static int count_alignments(const mismatch_bind_t *bind, uint64_t *cells, char *error, size_t error_size) {
    samFile *file = NULL;
    sam_hdr_t *header = NULL;
    hts_idx_t *index = NULL;
    hts_itr_t *iterator = NULL;
    bam1_t *alignment = NULL;
    char **regions = NULL;
    char *reference_locator = NULL;
    unsigned int region_count = 0;
    mismatch_mask_t mask;
    mismatch_window_t window;
    int read_status = 0;
    int ok = 0;

    memset(&mask, 0, sizeof(mask));
    memset(&window, 0, sizeof(window));
    window.tid = -1;

    file = sam_open(bind->path, "r");
    if (!file) {
        snprintf(error, error_size, "duckhts_bam_mismatch_counts: cannot open %s", bind->path);
        goto cleanup;
    }
    /* The CRAM decoder reads the FASTA index from the locator, so an index
     * that is not the adjacent .fai travels with the reference path. */
    reference_locator = duckhts_reference_locator(bind->reference, bind->reference_index_path);
    if (!reference_locator) {
        snprintf(error, error_size, "duckhts_bam_mismatch_counts: out of memory");
        goto cleanup;
    }
    if (hts_set_fai_filename(file, reference_locator) != 0) {
        snprintf(error, error_size, "duckhts_bam_mismatch_counts: cannot use reference %s", bind->reference);
        goto cleanup;
    }
    header = sam_hdr_read(file);
    alignment = bam_init1();
    if (!header || !alignment) {
        snprintf(error, error_size, "duckhts_bam_mismatch_counts: cannot read the header of %s", bind->path);
        goto cleanup;
    }
    window.fai = fai_load3(bind->reference, bind->reference_index_path, NULL, 0);
    if (!window.fai) {
        snprintf(error, error_size,
                 "duckhts_bam_mismatch_counts: cannot load the FASTA index of %s; build it with fasta_index()",
                 bind->reference);
        goto cleanup;
    }
    if (bind->mask && !mask_open(&mask, bind->mask, bind->mask_index_path, error, error_size)) goto cleanup;

    if (bind->region) {
        if (!duckhts_region_list_parse(bind->region, &regions, &region_count, error, error_size) ||
            !duckhts_region_list_validate(regions, region_count, region_name2id, header, error, error_size)) {
            goto cleanup;
        }
        for (unsigned int item = 0; item < region_count; item++) {
            if (strcmp(regions[item], ".") == 0 || strcmp(regions[item], "*") == 0) {
                snprintf(error, error_size,
                         "duckhts_bam_mismatch_counts: region items must name a contig or an interval");
                goto cleanup;
            }
        }
        index = sam_index_load3(file, bind->path, bind->index_path, HTS_IDX_SILENT_FAIL);
        if (!index) {
            snprintf(error, error_size, "duckhts_bam_mismatch_counts: region needs an index for %s", bind->path);
            goto cleanup;
        }
        iterator = sam_itr_regarray(index, header, regions, region_count);
        if (!iterator) {
            snprintf(error, error_size, "duckhts_bam_mismatch_counts: no region of '%s' is in %s",
                     bind->region, bind->path);
            goto cleanup;
        }
    }

    for (;;) {
        const char *contig;
        read_status = iterator ? sam_itr_next(file, iterator, alignment) : sam_read1(file, header, alignment);
        if (read_status < 0) break;
        if (alignment->core.tid < 0 || (alignment->core.flag & BAM_FUNMAP)) continue;
        if (alignment->core.flag & bind->exclude_flags) continue;
        if ((alignment->core.flag & bind->require_flags) != bind->require_flags) continue;
        if (alignment->core.qual < bind->min_mapq) continue;
        contig = sam_hdr_tid2name(header, alignment->core.tid);
        if (!contig) continue;
        if (!count_alignment(alignment, contig, bind->indel_flank, &window, &mask, cells, error, error_size)) {
            goto cleanup;
        }
    }
    if (read_status < -1) {
        snprintf(error, error_size, "duckhts_bam_mismatch_counts: cannot read an alignment of %s", bind->path);
        goto cleanup;
    }
    ok = 1;

cleanup:
    free(window.masked);
    free(window.bases);
    if (window.fai) fai_destroy(window.fai);
    mask_close(&mask);
    free(regions);
    free(reference_locator);
    if (iterator) hts_itr_destroy(iterator);
    if (index) hts_idx_destroy(index);
    if (alignment) bam_destroy1(alignment);
    if (header) sam_hdr_destroy(header);
    if (file) sam_close(file);
    return ok;
}

static void destroy_mismatch_bind(void *data) {
    mismatch_bind_t *bind = (mismatch_bind_t *)data;
    if (!bind) return;
    if (bind->path) duckdb_free(bind->path);
    if (bind->reference) duckdb_free(bind->reference);
    if (bind->region) duckdb_free(bind->region);
    if (bind->mask) duckdb_free(bind->mask);
    if (bind->index_path) duckdb_free(bind->index_path);
    if (bind->reference_index_path) duckdb_free(bind->reference_index_path);
    if (bind->mask_index_path) duckdb_free(bind->mask_index_path);
    duckdb_free(bind);
}

static void destroy_mismatch_state(void *data) {
    mismatch_state_t *state = (mismatch_state_t *)data;
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
static int named_integer(duckdb_bind_info info, const char *name, int fallback, int low, int high,
                         int *out) {
    duckdb_value value = duckdb_bind_get_named_parameter(info, name);
    int64_t number = fallback;
    if (value && !duckdb_is_null_value(value)) number = duckdb_get_int64(value);
    if (value) duckdb_destroy_value(&value);
    if (number < low || number > high) return 0;
    *out = (int)number;
    return 1;
}

static void mismatch_counts_bind(duckdb_bind_info info) {
    duckdb_value path_value = duckdb_bind_get_parameter(info, 0);
    duckdb_value reference_value = duckdb_bind_get_parameter(info, 1);
    mismatch_bind_t *bind = (mismatch_bind_t *)duckdb_malloc(sizeof(*bind));
    duckdb_logical_type tinyint_type, integer_type, bigint_type, varchar_type;

    if (bind) {
        memset(bind, 0, sizeof(*bind));
        if (!duckdb_is_null_value(path_value)) bind->path = duckdb_get_varchar(path_value);
        if (!duckdb_is_null_value(reference_value)) bind->reference = duckdb_get_varchar(reference_value);
    }
    duckdb_destroy_value(&path_value);
    duckdb_destroy_value(&reference_value);
    if (!bind) {
        duckdb_bind_set_error(info, "duckhts_bam_mismatch_counts: out of memory");
        return;
    }
    if (!bind->path || bind->path[0] == '\0' || !bind->reference || bind->reference[0] == '\0') {
        destroy_mismatch_bind(bind);
        duckdb_bind_set_error(info,
            "duckhts_bam_mismatch_counts requires an alignment path and a reference FASTA path");
        return;
    }
    bind->region = named_text(info, "region");
    bind->mask = named_text(info, "mask");
    bind->index_path = named_text(info, "index_path");
    bind->reference_index_path = named_text(info, "reference_index_path");
    bind->mask_index_path = named_text(info, "mask_index_path");
    if (!named_integer(info, "min_mapq", 20, 0, 255, &bind->min_mapq) ||
        !named_integer(info, "require_flags", 0, 0, 65535, &bind->require_flags) ||
        !named_integer(info, "exclude_flags", 3844, 0, 65535, &bind->exclude_flags) ||
        !named_integer(info, "indel_flank", 5, 0, MISMATCH_FLANK_MAX, &bind->indel_flank)) {
        destroy_mismatch_bind(bind);
        duckdb_bind_set_error(info,
            "duckhts_bam_mismatch_counts: min_mapq must be in [0, 255], the flag masks in [0, 65535] "
            "and indel_flank in [0, 1000]");
        return;
    }

    tinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_TINYINT);
    integer_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_bind_add_result_column(info, "mate", tinyint_type);
    duckdb_bind_add_result_column(info, "cycle", integer_type);
    duckdb_bind_add_result_column(info, "base_quality", tinyint_type);
    duckdb_bind_add_result_column(info, "reference_base", varchar_type);
    duckdb_bind_add_result_column(info, "read_base", varchar_type);
    duckdb_bind_add_result_column(info, "bases", bigint_type);
    duckdb_destroy_logical_type(&tinyint_type);
    duckdb_destroy_logical_type(&integer_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_bind_set_bind_data(info, bind, destroy_mismatch_bind);
}

static void mismatch_counts_init(duckdb_init_info info) {
    mismatch_state_t *state = (mismatch_state_t *)duckdb_malloc(sizeof(*state));
    if (state) {
        memset(state, 0, sizeof(*state));
        state->cells = calloc(MISMATCH_CELLS, sizeof(uint64_t));
    }
    if (!state || !state->cells) {
        destroy_mismatch_state(state);
        duckdb_init_set_error(info, "duckhts_bam_mismatch_counts: out of memory");
        return;
    }
    duckdb_init_set_init_data(info, state, destroy_mismatch_state);
    duckdb_init_set_max_threads(info, 1);
}

static void mismatch_counts_scan(duckdb_function_info info, duckdb_data_chunk output) {
    static const char base_text[MISMATCH_BASES] = {'A', 'C', 'G', 'T'};
    const mismatch_bind_t *bind = (const mismatch_bind_t *)duckdb_function_get_bind_data(info);
    mismatch_state_t *state = (mismatch_state_t *)duckdb_function_get_init_data(info);
    duckdb_vector vectors[MISMATCH_VECTOR_COLUMNS];
    int8_t *mate_data, *quality_data;
    int32_t *cycle_data;
    int64_t *bases_data;
    idx_t rows = 0;
    const idx_t capacity = duckdb_vector_size();

    if (!state->counted) {
        char error[512];
        error[0] = '\0';
        if (!count_alignments(bind, state->cells, error, sizeof(error))) {
            duckdb_function_set_error(info, error[0] ? error : "duckhts_bam_mismatch_counts failed");
            return;
        }
        state->counted = 1;
    }

    for (idx_t column = 0; column < MISMATCH_VECTOR_COLUMNS; column++) {
        vectors[column] = duckdb_data_chunk_get_vector(output, column);
    }
    mate_data = (int8_t *)duckdb_vector_get_data(vectors[0]);
    cycle_data = (int32_t *)duckdb_vector_get_data(vectors[1]);
    quality_data = (int8_t *)duckdb_vector_get_data(vectors[2]);
    bases_data = (int64_t *)duckdb_vector_get_data(vectors[5]);

    while (state->next_cell < MISMATCH_CELLS && rows < capacity) {
        size_t cell = state->next_cell++;
        const uint64_t count = state->cells[cell];
        int to, from, quality, cycle, mate;
        if (count == 0) continue;
        to = (int)(cell % MISMATCH_BASES); cell /= MISMATCH_BASES;
        from = (int)(cell % MISMATCH_BASES); cell /= MISMATCH_BASES;
        quality = (int)(cell % MISMATCH_QUALITIES); cell /= MISMATCH_QUALITIES;
        cycle = (int)(cell % (MISMATCH_CYCLE_MAX + 1)); cell /= (MISMATCH_CYCLE_MAX + 1);
        mate = (int)cell;
        mate_data[rows] = (int8_t)mate;
        cycle_data[rows] = cycle;
        if (quality == MISMATCH_QUALITY_MISSING) {
            duckdb_vector_ensure_validity_writable(vectors[2]);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vectors[2]), rows);
        } else {
            quality_data[rows] = (int8_t)quality;
        }
        duckdb_vector_assign_string_element_len(vectors[3], rows, &base_text[from], 1);
        duckdb_vector_assign_string_element_len(vectors[4], rows, &base_text[to], 1);
        bases_data[rows] = count > (uint64_t)INT64_MAX ? INT64_MAX : (int64_t)count;
        rows++;
    }
    duckdb_data_chunk_set_size(output, rows);
}

void register_duckhts_bam_mismatch_counts_function(duckdb_connection connection) {
    duckdb_table_function function = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);

    duckdb_table_function_set_name(function, "duckhts_bam_mismatch_counts");
    duckdb_table_function_add_parameter(function, varchar_type);
    duckdb_table_function_add_parameter(function, varchar_type);
    duckdb_table_function_add_named_parameter(function, "region", varchar_type);
    duckdb_table_function_add_named_parameter(function, "mask", varchar_type);
    duckdb_table_function_add_named_parameter(function, "index_path", varchar_type);
    duckdb_table_function_add_named_parameter(function, "reference_index_path", varchar_type);
    duckdb_table_function_add_named_parameter(function, "mask_index_path", varchar_type);
    duckdb_table_function_add_named_parameter(function, "min_mapq", bigint_type);
    duckdb_table_function_add_named_parameter(function, "require_flags", bigint_type);
    duckdb_table_function_add_named_parameter(function, "exclude_flags", bigint_type);
    duckdb_table_function_add_named_parameter(function, "indel_flank", bigint_type);
    duckdb_table_function_set_bind(function, mismatch_counts_bind);
    duckdb_table_function_set_init(function, mismatch_counts_init);
    duckdb_table_function_set_function(function, mismatch_counts_scan);
    duckdb_register_table_function(connection, function);
    duckdb_destroy_table_function(&function);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
}
