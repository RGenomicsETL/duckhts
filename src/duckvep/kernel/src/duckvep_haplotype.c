/*
 * duckvep_haplotype.c — multi-edit CDS haplotype mutation helpers.
 * See duckvep_haplotype.h. No allocation; no DuckDB/htslib.
 */
#include "duckvep_haplotype.h"
#include "duckvep_dna.h"

#include <limits.h>
#include <string.h>

static char haplo_norm_base(char c) {
    return duckvep_dna_normalize(c, 0);
}

static char haplo_norm_cds_base(uint8_t b) {
    return duckvep_dna_normalize((char)b, 1);
}

static char haplo_complement(char b) {
    return duckvep_dna_complement(b);
}

static char haplo_oriented_base(const uint8_t *seq, uint32_t len, uint32_t idx,
                                int reverse_complement) {
    char b;
    if (seq == NULL || idx >= len) return '\0';
    b = haplo_norm_base((char)seq[reverse_complement ? (len - 1u - idx) : idx]);
    if (b == '\0') return '\0';
    return reverse_complement ? haplo_complement(b) : b;
}

static void haplo_result_init(duckvep_haplotype_result_t *result) {
    if (result != NULL) memset(result, 0, sizeof *result);
}

static duckvep_haplotype_status_t haplo_fail(duckvep_haplotype_result_t *result,
                                             size_t *cds_len_out,
                                             size_t *protein_len_out,
                                             duckvep_haplotype_status_t status) {
    haplo_result_init(result);
    if (cds_len_out != NULL) *cds_len_out = 0u;
    if (protein_len_out != NULL) *protein_len_out = 0u;
    return status;
}

static int64_t diff_i64(uint32_t alt_len, uint32_t ref_len, int *ok) {
    int64_t a = (int64_t)alt_len;
    int64_t r = (int64_t)ref_len;
    *ok = 1;
    return a - r;
}

static int haplo_add_i64(int64_t a, int64_t b, int64_t *out) {
    if ((b > 0 && a > INT64_MAX - b) ||
        (b < 0 && a < INT64_MIN - b)) return 0;
    *out = a + b;
    return 1;
}

static int haplo_shift_coordinate(uint64_t coordinate, int64_t shift,
                                  uint64_t *out) {
    if (shift >= 0) {
        uint64_t add = (uint64_t)shift;
        if (add > UINT64_MAX - coordinate) return 0;
        *out = coordinate + add;
        return 1;
    }
    {
        uint64_t subtract = (uint64_t)(-(shift + 1)) + UINT64_C(1);
        if (subtract > coordinate) return 0;
        *out = coordinate - subtract;
    }
    return 1;
}

static duckvep_haplotype_status_t haplo_partition_pass(
    const duckvep_haplotype_edit_t *edits,
    size_t edit_count,
    duckvep_haplotype_block_t *blocks,
    size_t *block_count) {

    size_t block_begin = 0u;
    size_t count = 0u;
    size_t i;
    uint32_t previous_start = 0u;
    uint32_t previous_end = 0u;
    uint32_t block_end = 0u;
    int have_previous = 0;
    int block_saw_frameshift = 0;
    int64_t block_difference = 0;
    int64_t shift_before = 0;
    uint32_t block_flags = 0u;

    *block_count = 0u;
    for (i = 0u; i < edit_count; i++) {
        const duckvep_haplotype_edit_t *edit = &edits[i];
        uint32_t edit_end;
        int64_t difference;
        int ok_difference;
        int flush;

        if (edit->cds_start == 0u ||
            (edit->ref_len != 0u && edit->ref == NULL) ||
            (edit->alt_len != 0u && edit->alt == NULL) ||
            (edit->variant_strand != (int8_t)1 &&
             edit->variant_strand != (int8_t)-1)) {
            return DUCKVEP_HAPLOTYPE_INVALID_ARG;
        }
        if (edit->ref_len == 0u) {
            edit_end = edit->cds_start;
        } else {
            if (edit->ref_len - 1u > UINT32_MAX - edit->cds_start) {
                return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
            }
            edit_end = edit->cds_start + edit->ref_len - 1u;
        }
        if (have_previous &&
            (edit->cds_start <= previous_start ||
             edit->cds_start <= previous_end)) {
            return DUCKVEP_HAPLOTYPE_EDIT_ORDER;
        }
        previous_start = edit->cds_start;
        previous_end = edit_end;
        have_previous = 1;
        if (edit_end > block_end) block_end = edit_end;

        difference = diff_i64(edit->alt_len, edit->ref_len, &ok_difference);
        if (!ok_difference ||
            !haplo_add_i64(block_difference, difference, &block_difference)) {
            return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        }
        if (difference != 0) {
            block_flags |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
            if (difference % 3 != 0) block_saw_frameshift = 1;
        }

        flush = i + 1u == edit_count;
        if (!flush && block_difference % 3 == 0) {
            uint64_t alternate_start;
            uint64_t alternate_end;
            uint64_t next_start;

            if (!haplo_shift_coordinate((uint64_t)edit->cds_start - 1u,
                                        shift_before, &alternate_start) ||
                !haplo_shift_coordinate((uint64_t)edits[i + 1u].cds_start - 1u,
                                        block_difference, &next_start)) {
                return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
            }
            alternate_end = alternate_start;
            if (edit->alt_len != 0u) {
                if ((uint64_t)edit->alt_len - 1u >
                    UINT64_MAX - alternate_end) {
                    return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
                }
                alternate_end += (uint64_t)edit->alt_len - 1u;
            }
            flush = alternate_end / 3u != next_start / 3u;
        }
        if (!haplo_add_i64(shift_before, difference, &shift_before)) {
            return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        }
        if (!flush) continue;

        if (block_saw_frameshift) {
            if (block_difference % 3 == 0) {
                block_flags |= DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT;
            } else {
                block_flags |= DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT;
            }
        }
        if (blocks != NULL) {
            blocks[count].edit_begin = block_begin;
            blocks[count].edit_count = i - block_begin + 1u;
            blocks[count].cds_start = edits[block_begin].cds_start;
            blocks[count].cds_end = block_end;
            blocks[count].length_diff = block_difference;
            blocks[count].flags = block_flags;
        }
        count++;
        block_begin = i + 1u;
        block_end = 0u;
        block_difference = 0;
        shift_before = 0;
        block_flags = 0u;
        block_saw_frameshift = 0;
    }
    *block_count = count;
    return DUCKVEP_HAPLOTYPE_OK;
}

duckvep_haplotype_status_t duckvep_haplotype_partition(
    const duckvep_haplotype_edit_t *edits,
    size_t edit_count,
    duckvep_haplotype_block_t *blocks,
    size_t block_capacity,
    size_t *required_blocks) {

    duckvep_haplotype_status_t status;
    size_t needed;

    if (required_blocks == NULL || (edit_count != 0u && edits == NULL) ||
        (block_capacity != 0u && blocks == NULL)) {
        return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    }
    *required_blocks = 0u;
    status = haplo_partition_pass(edits, edit_count, NULL, &needed);
    if (status != DUCKVEP_HAPLOTYPE_OK) return status;
    *required_blocks = needed;
    if (needed > block_capacity) return DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL;
    return haplo_partition_pass(edits, edit_count, blocks, required_blocks);
}

/* Exact input/output aliasing was accepted by the original implementation.
 * Keep that compatibility path, where memmove is required because cumulative
 * edit displacement can change sign within one haplotype. Normal callers use
 * distinct reference and scratch buffers and take the linear rebuild below. */
static duckvep_haplotype_status_t haplo_apply_exact_alias(
    uint8_t                         *cds,
    size_t                           ref_cds_len,
    const duckvep_haplotype_edit_t  *edits,
    size_t                           edit_count,
    int8_t                           transcript_strand,
    size_t                          *cds_len_out) {

    size_t cur_len = ref_cds_len;
    size_t i;

    for (i = 0u; i < ref_cds_len; i++) {
        char b = haplo_norm_cds_base(cds[i]);
        if (b == '\0') return DUCKVEP_HAPLOTYPE_INVALID_BASE;
        cds[i] = (uint8_t)b;
    }

    for (i = 0u; i < edit_count; i++) {
        const duckvep_haplotype_edit_t *e = &edits[i];
        size_t start0 = (size_t)e->cds_start - 1u;
        size_t tail_start = start0 + (size_t)e->ref_len;
        size_t tail_len = cur_len - tail_start;
        size_t new_len = cur_len - (size_t)e->ref_len + (size_t)e->alt_len;
        uint32_t j;
        int reverse = e->variant_strand != transcript_strand;

        memmove(cds + start0 + (size_t)e->alt_len,
                cds + tail_start, tail_len);
        for (j = 0u; j < e->alt_len; j++) {
            cds[start0 + (size_t)j] = (uint8_t)haplo_oriented_base(
                e->alt, e->alt_len, j, reverse);
        }
        cur_len = new_len;
    }

    *cds_len_out = cur_len;
    return DUCKVEP_HAPLOTYPE_OK;
}

duckvep_haplotype_status_t duckvep_haplotype_apply_cds_edits(
    const uint8_t                    *ref_cds,
    size_t                            ref_cds_len,
    const duckvep_haplotype_edit_t    *edits,
    size_t                            edit_count,
    int8_t                            transcript_strand,
    uint8_t                          *cds_out,
    size_t                            cds_cap,
    size_t                           *cds_len_out,
    duckvep_haplotype_result_t       *result) {

    size_t final_len = ref_cds_len;
    size_t peak_len = ref_cds_len;
    size_t i;
    uint32_t prev_start = UINT32_MAX;
    uint32_t prev_ref_len = 0u;
    int saw_frameshift_edit = 0;
    int64_t total_diff = 0;
    uint32_t flags = 0u;

    haplo_result_init(result);
    if (cds_len_out != NULL) *cds_len_out = 0u;

    if (ref_cds == NULL || cds_out == NULL || cds_len_out == NULL ||
        (edit_count > 0u && edits == NULL) ||
        (transcript_strand != (int8_t)1 && transcript_strand != (int8_t)-1)) {
        return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    }
    /* Validate the complete edit set and measure the final sequence before
     * publishing output. Coordinates stay on the original CDS axis. */
    for (i = 0u; i < edit_count; i++) {
        const duckvep_haplotype_edit_t *e = &edits[i];
        size_t start0;
        size_t new_len;
        uint32_t j;
        int reverse;
        int ok_diff = 0;
        int64_t d;
        uint32_t effective_end;

        if (e->cds_start == 0u ||
            (e->variant_strand != (int8_t)1 && e->variant_strand != (int8_t)-1) ||
            (e->ref_len > 0u && e->ref == NULL) ||
            (e->alt_len > 0u && e->alt == NULL)) {
            return haplo_fail(result, cds_len_out, NULL, DUCKVEP_HAPLOTYPE_INVALID_ARG);
        }

        if (e->ref_len != 0u &&
            e->ref_len - 1u > UINT32_MAX - e->cds_start) {
            return haplo_fail(result, cds_len_out, NULL, DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
        }
        effective_end = e->ref_len == 0u ? e->cds_start :
            e->cds_start + e->ref_len - 1u;
        if (i > 0u) {
            if (e->cds_start > prev_start) {
                return haplo_fail(result, cds_len_out, NULL, DUCKVEP_HAPLOTYPE_EDIT_ORDER);
            }
            if (e->ref_len > 0u && effective_end >= prev_start) {
                return haplo_fail(result, cds_len_out, NULL, DUCKVEP_HAPLOTYPE_EDIT_ORDER);
            }
            if (e->ref_len == 0u && e->cds_start == prev_start && prev_ref_len > 0u) {
                return haplo_fail(result, cds_len_out, NULL, DUCKVEP_HAPLOTYPE_EDIT_ORDER);
            }
        }
        prev_start = e->cds_start;
        prev_ref_len = e->ref_len;

        start0 = (size_t)e->cds_start - 1u;
        if (start0 > ref_cds_len ||
            (size_t)e->ref_len > ref_cds_len - start0) {
            return haplo_fail(result, cds_len_out, NULL, DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
        }

        reverse = (e->variant_strand != transcript_strand);
        for (j = 0u; j < e->ref_len; j++) {
            char expected = haplo_oriented_base(e->ref, e->ref_len, j, reverse);
            char observed = haplo_norm_cds_base(ref_cds[start0 + (size_t)j]);
            if (expected == '\0') return haplo_fail(result, cds_len_out, NULL, DUCKVEP_HAPLOTYPE_INVALID_BASE);
            if (observed == '\0') return haplo_fail(result, cds_len_out, NULL, DUCKVEP_HAPLOTYPE_INVALID_BASE);
            if (observed == 'N' || observed != expected) {
                return haplo_fail(result, cds_len_out, NULL, DUCKVEP_HAPLOTYPE_REF_MISMATCH);
            }
        }
        for (j = 0u; j < e->alt_len; j++) {
            if (haplo_oriented_base(e->alt, e->alt_len, j, reverse) == '\0') {
                return haplo_fail(result, cds_len_out, NULL, DUCKVEP_HAPLOTYPE_INVALID_BASE);
            }
        }

        d = diff_i64(e->alt_len, e->ref_len, &ok_diff);
        if (!ok_diff) return haplo_fail(result, cds_len_out, NULL, DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
        if (!haplo_add_i64(total_diff, d, &total_diff)) {
            return haplo_fail(result, cds_len_out, NULL,
                              DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
        }
        if (d != 0) flags |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
        if ((d % 3) != 0) saw_frameshift_edit = 1;

        if ((size_t)e->ref_len > final_len ||
            (size_t)e->alt_len >
                SIZE_MAX - (final_len - (size_t)e->ref_len)) {
            return haplo_fail(result, cds_len_out, NULL, DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
        }
        new_len = final_len - (size_t)e->ref_len + (size_t)e->alt_len;
        final_len = new_len;
        if (final_len > peak_len) peak_len = final_len;
    }

    if (final_len > cds_cap || (ref_cds == cds_out && peak_len > cds_cap)) {
        return haplo_fail(result, cds_len_out, NULL,
                          DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL);
    }

    if (ref_cds == cds_out) {
        duckvep_haplotype_status_t alias_status = haplo_apply_exact_alias(
            cds_out, ref_cds_len, edits, edit_count, transcript_strand,
            cds_len_out);
        if (alias_status != DUCKVEP_HAPLOTYPE_OK) {
            return haplo_fail(result, cds_len_out, NULL, alias_status);
        }
    } else {
        /* Rebuild once from right to left. Every reference byte and ALT byte is
         * written exactly once; edit count no longer multiplies CDS-tail moves. */
        size_t src_cursor = ref_cds_len;
        size_t dst_cursor = final_len;

        for (i = 0u; i < edit_count; i++) {
            const duckvep_haplotype_edit_t *e = &edits[i];
            size_t start0 = (size_t)e->cds_start - 1u;
            size_t ref_end = start0 + (size_t)e->ref_len;
            size_t suffix_len;
            uint32_t j;
            int reverse = e->variant_strand != transcript_strand;

            if (ref_end > src_cursor) {
                return haplo_fail(result, cds_len_out, NULL,
                                  DUCKVEP_HAPLOTYPE_EDIT_ORDER);
            }
            suffix_len = src_cursor - ref_end;
            if (suffix_len > dst_cursor ||
                (size_t)e->alt_len > dst_cursor - suffix_len) {
                return haplo_fail(result, cds_len_out, NULL,
                                  DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
            }
            while (src_cursor > ref_end) {
                char b = haplo_norm_cds_base(ref_cds[--src_cursor]);
                if (b == '\0') {
                    return haplo_fail(result, cds_len_out, NULL,
                                      DUCKVEP_HAPLOTYPE_INVALID_BASE);
                }
                cds_out[--dst_cursor] = (uint8_t)b;
            }
            src_cursor = start0;
            for (j = e->alt_len; j > 0u; j--) {
                cds_out[--dst_cursor] = (uint8_t)haplo_oriented_base(
                    e->alt, e->alt_len, j - 1u, reverse);
            }
        }
        if (src_cursor > dst_cursor) {
            return haplo_fail(result, cds_len_out, NULL,
                              DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
        }
        while (src_cursor > 0u) {
            char b = haplo_norm_cds_base(ref_cds[--src_cursor]);
            if (b == '\0') {
                return haplo_fail(result, cds_len_out, NULL,
                                  DUCKVEP_HAPLOTYPE_INVALID_BASE);
            }
            cds_out[--dst_cursor] = (uint8_t)b;
        }
        if (dst_cursor != 0u) {
            return haplo_fail(result, cds_len_out, NULL,
                              DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
        }
        *cds_len_out = final_len;
    }

    if (final_len < cds_cap) cds_out[final_len] = (uint8_t)'\0';
    if (saw_frameshift_edit) {
        if ((total_diff % 3) == 0) flags |= DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT;
        else flags |= DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT;
    }
    if (result != NULL) {
        result->cds_len = final_len;
        result->length_diff = total_diff;
        result->flags = flags;
        result->applied_edits = edit_count;
    }
    return DUCKVEP_HAPLOTYPE_OK;
}

duckvep_haplotype_status_t duckvep_haplotype_translate_cds(
    const uint8_t                    *cds,
    size_t                            cds_len,
    duckvep_codon_table_t             table,
    uint8_t                          *protein_out,
    size_t                            protein_cap,
    size_t                           *protein_len_out,
    duckvep_haplotype_result_t       *result) {

    size_t codon_count;
    size_t i;
    size_t out_len = 0u;
    uint32_t flags = 0u;

    haplo_result_init(result);
    if (protein_len_out != NULL) *protein_len_out = 0u;
    if (cds == NULL || protein_out == NULL || protein_len_out == NULL ||
        protein_cap == 0u || !duckvep_codon_table_supported(table)) {
        return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    }
    codon_count = cds_len / 3u;

    for (i = 0u; i < codon_count; i++) {
        char codon[4];
        char aa;
        if (out_len + 1u >= protein_cap) {
            return haplo_fail(result, NULL, protein_len_out, DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL);
        }
        codon[0] = haplo_norm_cds_base(cds[i * 3u]);
        codon[1] = haplo_norm_cds_base(cds[i * 3u + 1u]);
        codon[2] = haplo_norm_cds_base(cds[i * 3u + 2u]);
        codon[3] = '\0';
        if (codon[0] == '\0' || codon[1] == '\0' || codon[2] == '\0') {
            return haplo_fail(result, NULL, protein_len_out, DUCKVEP_HAPLOTYPE_INVALID_BASE);
        }
        aa = duckvep_translate_codon(codon, table);
        protein_out[out_len++] = (uint8_t)aa;
        if (aa == '*') {
            if (i + 1u < codon_count) flags |= DUCKVEP_HAPLOTYPE_FLAG_STOP_TRUNCATED;
            break;
        }
    }
    protein_out[out_len] = (uint8_t)'\0';
    *protein_len_out = out_len;
    if (result != NULL) {
        result->cds_len = cds_len;
        result->protein_len = out_len;
        result->flags = flags;
    }
    return DUCKVEP_HAPLOTYPE_OK;
}
