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
                                             duckvep_haplotype_status_t status) {
    haplo_result_init(result);
    if (cds_len_out != NULL) *cds_len_out = 0u;
    return status;
}

static int haplo_overlaps_output(const void *input, size_t input_len,
                                 const uint8_t *output, size_t output_cap) {
    uintptr_t src = (uintptr_t)input;
    uintptr_t dst = (uintptr_t)output;
    if (input_len == 0u || output_cap == 0u) return 0;
    /* Subtract addresses instead of forming a potentially wrapped end. */
    return src <= dst ? dst - src < input_len : src - dst < output_cap;
}

static int haplo_add_i64(int64_t a, int64_t b, int64_t *out) {
    if ((b > 0 && a > INT64_MAX - b) ||
        (b < 0 && a < INT64_MIN - b)) return 0;
    *out = a + b;
    return 1;
}

duckvep_haplotype_status_t duckvep_haplotype_compose_replacements(
    const uint8_t *reference, size_t reference_length,
    const duckvep_haplotype_edit_t *edits, size_t edit_count, int8_t transcript_strand,
    uint64_t *source_ids, uint8_t *cds, size_t cds_capacity,
    duckvep_haplotype_block_t *components, size_t component_capacity,
    size_t *component_count, duckvep_haplotype_result_t *result) {
    if (component_count) *component_count = 0u;
    haplo_result_init(result);
    if (!reference || !cds || !component_count || !result ||
        reference_length > UINT32_MAX ||
        (edit_count && (!edits || !source_ids || !components)) ||
        edit_count > SIZE_MAX / sizeof(*edits) ||
        edit_count > SIZE_MAX / sizeof(*source_ids) ||
        component_capacity > SIZE_MAX / sizeof(*components) ||
        (transcript_strand != 1 && transcript_strand != -1))
        return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    if (component_capacity < edit_count) return DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL;
    size_t edit_bytes = edit_count * sizeof(*edits);
    size_t id_bytes = edit_count * sizeof(*source_ids);
    size_t component_bytes = component_capacity * sizeof(*components);
#define OVERLAPS(p, n, q, m) haplo_overlaps_output((p), (n), (const uint8_t *)(q), (m))
    if (OVERLAPS(reference, reference_length, cds, cds_capacity) ||
        OVERLAPS(edits, edit_bytes, cds, cds_capacity) ||
        OVERLAPS(source_ids, id_bytes, cds, cds_capacity) ||
        OVERLAPS(components, component_bytes, cds, cds_capacity) ||
        OVERLAPS(reference, reference_length, source_ids, id_bytes) ||
        OVERLAPS(edits, edit_bytes, source_ids, id_bytes) ||
        OVERLAPS(reference, reference_length, components, component_bytes) ||
        OVERLAPS(edits, edit_bytes, components, component_bytes) ||
        OVERLAPS(source_ids, id_bytes, components, component_bytes))
        return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    size_t cursor = reference_length, suffix = 0u, peak = 0u;
    int64_t nominal_difference = 0;
    uint32_t flags = 0u;
    int frame_changed = 0;
    for (size_t i = 0u; i < reference_length; i++)
        if (!haplo_norm_cds_base(reference[i])) return DUCKVEP_HAPLOTYPE_INVALID_BASE;
    /* Length planning does not depend on sequence equality: replacing equal bytes
     * leaves the same representation and length, but creates no source contribution. */
    for (size_t i = 0u; i < edit_count; i++) {
        const duckvep_haplotype_edit_t *e = &edits[i];
        if (!e->cds_start || !e->ref_len || !e->ref || (e->alt_len && !e->alt) ||
            (e->variant_strand != 1 && e->variant_strand != -1) ||
            OVERLAPS(e->ref, e->ref_len, cds, cds_capacity) ||
            OVERLAPS(e->alt, e->alt_len, cds, cds_capacity) ||
            OVERLAPS(e->ref, e->ref_len, source_ids, id_bytes) ||
            OVERLAPS(e->alt, e->alt_len, source_ids, id_bytes) ||
            OVERLAPS(e->ref, e->ref_len, components, component_bytes) ||
            OVERLAPS(e->alt, e->alt_len, components, component_bytes))
            return DUCKVEP_HAPLOTYPE_INVALID_ARG;
        size_t start0 = (size_t)e->cds_start - 1u;
        if (start0 > reference_length || e->ref_len > reference_length - start0)
            return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        if (start0 > cursor) return DUCKVEP_HAPLOTYPE_EDIT_ORDER;
        int reverse = e->variant_strand != transcript_strand;
        for (uint32_t j = 0u; j < e->ref_len; j++) {
            char base = haplo_oriented_base(e->ref, e->ref_len, j, reverse);
            if (!base) return DUCKVEP_HAPLOTYPE_INVALID_BASE;
            if (base != haplo_norm_cds_base(reference[start0 + j]))
                return DUCKVEP_HAPLOTYPE_REF_MISMATCH;
        }
        for (uint32_t j = 0u; j < e->alt_len; j++)
            if (!haplo_oriented_base(e->alt, e->alt_len, j, reverse))
                return DUCKVEP_HAPLOTYPE_INVALID_BASE;
        size_t prefix = cursor - start0;
        if (e->ref_len <= prefix) {
            size_t gap = prefix - e->ref_len;
            if (gap > SIZE_MAX - suffix) return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
            suffix += gap;
        } else {
            size_t drop = e->ref_len - prefix;
            suffix -= drop < suffix ? drop : suffix;
        }
        if (e->alt_len > SIZE_MAX - suffix) return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        suffix += e->alt_len;
        if (suffix > INT64_MAX) return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        if (suffix > peak) peak = suffix;
        cursor = start0;
        int64_t delta = (int64_t)e->alt_len - e->ref_len;
        if (!haplo_add_i64(nominal_difference, delta, &nominal_difference))
            return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        if (delta) flags |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
        if (delta % 3) frame_changed = 1;
    }
#undef OVERLAPS
    if (cursor > SIZE_MAX - suffix || cursor + suffix > INT64_MAX)
        return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
    size_t final_length = cursor + suffix;
    if (final_length > peak) peak = final_length;
    if (peak > cds_capacity) return DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL;

    cursor = reference_length;
    size_t at = cds_capacity, groups = 0u, changed_count = 0u;
    for (size_t i = 0u; i < edit_count; i++) {
        const duckvep_haplotype_edit_t *e = &edits[i];
        uint64_t source_id = source_ids[i];
        size_t start0 = (size_t)e->cds_start - 1u, prefix = cursor - start0;
        /* min(ref_len, prefix + suffix), without forming an unbounded sum. */
        size_t removed = e->ref_len;
        if (removed > prefix && removed - prefix > cds_capacity - at)
            removed = prefix + (cds_capacity - at);
        int reverse = e->variant_strand != transcript_strand;
        int changed = removed != e->alt_len;
        for (size_t j = 0u; !changed && j < removed; j++) {
            char current = j < prefix ? haplo_norm_cds_base(reference[start0 + j])
                                      : (char)cds[at + j - prefix];
            changed = current != haplo_oriented_base(e->alt, e->alt_len, (uint32_t)j, reverse);
        }
        size_t ref_end = start0 + e->ref_len, output_end, keep = groups, begin = changed_count;
        if (e->ref_len <= prefix) {
            size_t gap = prefix - e->ref_len;
            at -= gap;
            for (size_t j = 0u; j < gap; j++)
                cds[at + j] = (uint8_t)haplo_norm_cds_base(reference[ref_end + j]);
            output_end = at;
        } else {
            size_t consumed = at + (removed - prefix);
            output_end = consumed;
            if (changed) {
                size_t mapped_ref = cursor, mapped_out = at;
                int inside = 0;
                while (keep) {
                    const duckvep_haplotype_block_t *g = &components[keep - 1u];
                    size_t end = g->alt_start0 + g->alt_len;
                    if (!(g->alt_start0 < consumed || (!g->alt_len && end == consumed))) break;
                    keep--;
                    if (consumed <= end) {
                        ref_end = (size_t)g->cds_start - 1u + g->ref_len;
                        output_end = end;
                        inside = 1;
                        break;
                    }
                    mapped_ref = (size_t)g->cds_start - 1u + g->ref_len;
                    mapped_out = end;
                }
                if (!inside) ref_end = mapped_ref + (consumed - mapped_out);
                if (keep < groups) begin = components[keep].edit_begin;
            }
            at = consumed;
        }
        at -= e->alt_len;
        for (uint32_t j = 0u; j < e->alt_len; j++)
            cds[at + j] = (uint8_t)haplo_oriented_base(e->alt, e->alt_len, j, reverse);
        cursor = start0;
        if (changed) {
            source_ids[changed_count++] = source_id;
            int64_t delta = (int64_t)(output_end - at) - (int64_t)(ref_end - start0);
            uint32_t component_flags = delta ? DUCKVEP_HAPLOTYPE_FLAG_INDEL : 0u;
            if (delta % 3) component_flags |= DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT;
            components[keep] = (duckvep_haplotype_block_t){begin, changed_count - begin,
                e->cds_start, (uint32_t)(ref_end - start0), at, output_end - at,
                delta, component_flags};
            groups = keep + 1u;
        }
    }
    at -= cursor;
    for (size_t j = 0u; j < cursor; j++) cds[at + j] = (uint8_t)haplo_norm_cds_base(reference[j]);
    memmove(cds, cds + at, final_length);
    for (size_t i = 0u; i < groups; i++) components[i].alt_start0 -= at;
    if (final_length < cds_capacity) cds[final_length] = 0u;
    if (frame_changed) flags |= nominal_difference % 3 ? DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT
                                                     : DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT;
    *component_count = groups;
    *result = (duckvep_haplotype_result_t){final_length, nominal_difference, flags, changed_count};
    return DUCKVEP_HAPLOTYPE_OK;
}

duckvep_haplotype_status_t duckvep_haplotype_reference_protein(
    const uint8_t *cds, size_t cds_length, duckvep_codon_table_t table,
    const uint32_t *edit_positions1, const uint8_t *edit_alternates, size_t edit_count,
    uint8_t *peptide, size_t capacity, size_t *length) {
    if (!length) return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    *length = 0u;
    if (!cds || !peptide || !duckvep_codon_table_supported(table) ||
        edit_count > SIZE_MAX / sizeof(*edit_positions1) ||
        (edit_count && (!edit_positions1 || !edit_alternates))) return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    if (cds_length < 3u) return DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE;
    size_t codons = cds_length / 3u;
    if (capacity < codons + 2u) return DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL;
    if (haplo_overlaps_output(cds, cds_length, peptide, capacity) ||
        haplo_overlaps_output(edit_positions1, edit_count * sizeof(*edit_positions1), peptide, capacity) ||
        haplo_overlaps_output(edit_alternates, edit_count, peptide, capacity))
        return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    for (size_t i = 0u; i < edit_count; i++) {
        uint32_t position = edit_positions1[i]; uint8_t aa = edit_alternates[i];
        if (!position || position > codons || (i && position <= edit_positions1[i - 1u]))
            return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        if (aa != '*' && (aa < 'A' || aa > 'Z')) return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    }
    duckvep_translation_t translated;
    duckvep_translation_status_t status = duckvep_translate_cds(cds, cds_length, table,
        DUCKVEP_TRANSLATION_N_CONSENSUS, peptide, capacity, &translated);
    if (status != DUCKVEP_TRANSLATION_OK) return status == DUCKVEP_TRANSLATION_INVALID_BASE
        ? DUCKVEP_HAPLOTYPE_INVALID_BASE : DUCKVEP_HAPLOTYPE_INVALID_ARG;
    size_t n = translated.length;
    if (peptide[n - 1u] == '*') n--;
    if (duckvep_codon_is_start(cds, table)) {
        peptide[0] = 'M';
        if (!n) n = 1u;
    }
    /* Descending single-residue edits match SeqEdit::apply_edit. Only an edit
     * at the stripped terminal position can extend this peptide. */
    for (size_t i = edit_count; i > 0u; i--) {
        size_t at = (size_t)edit_positions1[i - 1u] - 1u;
        if (at > n) return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        peptide[at] = edit_alternates[i - 1u];
        if (at == n) n++;
    }
    const uint8_t *last = cds + cds_length - 3u;
    if (!memcmp(last, "TAA", 3u) || !memcmp(last, "TAG", 3u) || !memcmp(last, "TGA", 3u))
        peptide[n++] = '*';
    peptide[n] = 0u; *length = n;
    return DUCKVEP_HAPLOTYPE_OK;
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
    int have_previous = 0;
    int block_saw_frameshift = 0;
    int64_t block_difference = 0;
    int64_t shift_before = 0;
    int64_t completed_shift = 0;
    uint32_t block_flags = 0u;

    *block_count = 0u;
    for (i = 0u; i < edit_count; i++) {
        const duckvep_haplotype_edit_t *edit = &edits[i];
        uint32_t edit_end;
        int64_t difference;
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

        difference = (int64_t)edit->alt_len - (int64_t)edit->ref_len;
        if (!haplo_add_i64(block_difference, difference, &block_difference)) {
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
        uint32_t cds_start = edits[block_begin].cds_start;
        uint64_t ref_len = (uint64_t)edit->cds_start - cds_start + edit->ref_len;
        uint64_t alt_start0, alt_len;
        if (ref_len > UINT32_MAX ||
            !haplo_shift_coordinate((uint64_t)cds_start - 1u, completed_shift, &alt_start0) ||
            !haplo_shift_coordinate(ref_len, block_difference, &alt_len) ||
            alt_start0 > SIZE_MAX || alt_len > SIZE_MAX - alt_start0 ||
            !haplo_add_i64(completed_shift, block_difference, &completed_shift)) {
            return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        }
        if (blocks != NULL) {
            blocks[count].edit_begin = block_begin;
            blocks[count].edit_count = i - block_begin + 1u;
            blocks[count].cds_start = cds_start;
            blocks[count].ref_len = (uint32_t)ref_len;
            blocks[count].alt_start0 = (size_t)alt_start0;
            blocks[count].alt_len = (size_t)alt_len;
            blocks[count].length_diff = block_difference;
            blocks[count].flags = block_flags;
        }
        count++;
        block_begin = i + 1u;
        block_difference = 0;
        shift_before = 0;
        block_flags = 0u;
        block_saw_frameshift = 0;
    }
    *block_count = count;
    return DUCKVEP_HAPLOTYPE_OK;
}

duckvep_haplotype_status_t duckvep_haplotype_block_frame_intersects(
    const duckvep_haplotype_edit_t  *edits,
    size_t                          edit_count,
    const duckvep_haplotype_block_t *block,
    size_t                          alt_start0,
    size_t                          alt_length,
    int                            *intersects) {

    if (intersects) *intersects = 0;
    if (!edits || !block || !intersects || !block->edit_count || !block->cds_start ||
        block->edit_begin > edit_count || block->edit_count > edit_count - block->edit_begin) {
        return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    }
    if (alt_length > SIZE_MAX - alt_start0 || block->alt_start0 > INT64_MAX) {
        return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
    }
    size_t end0 = alt_start0 + alt_length;
    int found = 0, open = 0, saw_displacement = 0;
    uint64_t open_start0 = 0u, previous_end1 = 0u, last_ref_end0 = 0u, last_alt_end0 = 0u;
    int64_t difference = 0;
    int64_t shift = (int64_t)block->alt_start0 - ((int64_t)block->cds_start - 1);
    uint32_t flags = 0u;
    if (shift % 3 != 0 || edits[block->edit_begin].cds_start != block->cds_start) {
        return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    }
    for (size_t i = 0u; i < block->edit_count; i++) {
        const duckvep_haplotype_edit_t *edit = edits + block->edit_begin + i;
        if (!edit->cds_start || (edit->ref_len && !edit->ref) ||
            (edit->alt_len && !edit->alt) ||
            (edit->variant_strand != 1 && edit->variant_strand != -1)) {
            return DUCKVEP_HAPLOTYPE_INVALID_ARG;
        }
        if (i && edit->cds_start <= previous_end1) return DUCKVEP_HAPLOTYPE_EDIT_ORDER;
        last_ref_end0 = (uint64_t)edit->cds_start - 1u + edit->ref_len;
        if (last_ref_end0 > UINT32_MAX) return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        previous_end1 = edit->ref_len ? last_ref_end0 : edit->cds_start;
        uint64_t edit_start0;
        if (!haplo_shift_coordinate((uint64_t)edit->cds_start - 1u, shift, &edit_start0) ||
            edit_start0 > SIZE_MAX || edit->alt_len > SIZE_MAX - edit_start0) {
            return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        }
        last_alt_end0 = edit_start0 + edit->alt_len;
        int64_t change = (int64_t)edit->alt_len - edit->ref_len;
        if (!haplo_add_i64(difference, change, &difference) ||
            !haplo_add_i64(shift, change, &shift)) {
            return DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        }
        if (change) flags |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
        if (change % 3) saw_displacement = 1;
        if (!open && difference % 3) {
            open_start0 = edit_start0;
            open = 1;
        } else if (open && difference % 3 == 0) {
            if (alt_length && last_alt_end0 > open_start0 &&
                (uint64_t)alt_start0 < last_alt_end0 &&
                (uint64_t)end0 > open_start0) found = 1;
            open = 0;
        }
    }
    if (saw_displacement) flags |= difference % 3
        ? DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT : DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT;
    if (difference != block->length_diff || flags != block->flags ||
        last_ref_end0 - ((uint64_t)block->cds_start - 1u) != block->ref_len ||
        last_alt_end0 < block->alt_start0 || last_alt_end0 - block->alt_start0 != block->alt_len) {
        return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    }
    if (open && alt_length && (uint64_t)end0 > open_start0) found = 1;
    *intersects = found;
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

    size_t final_len;
    uint64_t measured_len;
    size_t i;
    uint32_t prev_start = UINT32_MAX;
    int saw_frameshift_edit = 0;
    int64_t total_diff = 0;
    uint32_t flags = 0u;

    haplo_result_init(result);
    if (cds_len_out != NULL) *cds_len_out = 0u;

    if (ref_cds == NULL || cds_out == NULL || cds_len_out == NULL ||
        (edit_count > 0u && edits == NULL) ||
        edit_count > SIZE_MAX / sizeof *edits ||
        haplo_overlaps_output(ref_cds, ref_cds_len, cds_out, cds_cap) ||
        haplo_overlaps_output(edits, edit_count * sizeof *edits, cds_out, cds_cap) ||
        (transcript_strand != (int8_t)1 && transcript_strand != (int8_t)-1)) {
        return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    }
    /* Validate the complete edit set and measure the final sequence before
     * publishing output. Coordinates stay on the original CDS axis. */
    for (i = 0u; i < edit_count; i++) {
        const duckvep_haplotype_edit_t *e = &edits[i];
        size_t start0;
        uint32_t j;
        int reverse;
        int64_t d;
        uint32_t effective_end;

        if (e->cds_start == 0u ||
            (e->variant_strand != (int8_t)1 && e->variant_strand != (int8_t)-1) ||
            (e->ref_len > 0u && e->ref == NULL) ||
            (e->alt_len > 0u && e->alt == NULL) ||
            haplo_overlaps_output(e->ref, e->ref_len, cds_out, cds_cap) ||
            haplo_overlaps_output(e->alt, e->alt_len, cds_out, cds_cap)) {
            return haplo_fail(result, cds_len_out, DUCKVEP_HAPLOTYPE_INVALID_ARG);
        }

        if (e->ref_len != 0u &&
            e->ref_len - 1u > UINT32_MAX - e->cds_start) {
            return haplo_fail(result, cds_len_out, DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
        }
        effective_end = e->ref_len == 0u ? e->cds_start :
            e->cds_start + e->ref_len - 1u;
        /* The same occupied-site test as the ascending partitioner: an
         * insertion has no REF bases but still owns its interbase site. */
        if (i > 0u && effective_end >= prev_start) {
            return haplo_fail(result, cds_len_out, DUCKVEP_HAPLOTYPE_EDIT_ORDER);
        }
        prev_start = e->cds_start;

        start0 = (size_t)e->cds_start - 1u;
        if (start0 > ref_cds_len ||
            (size_t)e->ref_len > ref_cds_len - start0) {
            return haplo_fail(result, cds_len_out, DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
        }

        reverse = (e->variant_strand != transcript_strand);
        for (j = 0u; j < e->ref_len; j++) {
            char expected = haplo_oriented_base(e->ref, e->ref_len, j, reverse);
            char observed = haplo_norm_cds_base(ref_cds[start0 + (size_t)j]);
            if (expected == '\0') return haplo_fail(result, cds_len_out, DUCKVEP_HAPLOTYPE_INVALID_BASE);
            if (observed == '\0') return haplo_fail(result, cds_len_out, DUCKVEP_HAPLOTYPE_INVALID_BASE);
            if (observed == 'N' || observed != expected) {
                return haplo_fail(result, cds_len_out, DUCKVEP_HAPLOTYPE_REF_MISMATCH);
            }
        }
        for (j = 0u; j < e->alt_len; j++) {
            if (haplo_oriented_base(e->alt, e->alt_len, j, reverse) == '\0') {
                return haplo_fail(result, cds_len_out, DUCKVEP_HAPLOTYPE_INVALID_BASE);
            }
        }

        d = (int64_t)e->alt_len - (int64_t)e->ref_len;
        if (!haplo_add_i64(total_diff, d, &total_diff)) {
            return haplo_fail(result, cds_len_out,
                              DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
        }
        if (d != 0) flags |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
        if ((d % 3) != 0) saw_frameshift_edit = 1;

    }

    if (!haplo_shift_coordinate(ref_cds_len, total_diff, &measured_len) ||
        measured_len > SIZE_MAX) {
        return haplo_fail(result, cds_len_out, DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
    }
    final_len = (size_t)measured_len;
    if (final_len > cds_cap) {
        return haplo_fail(result, cds_len_out,
                          DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL);
    }

    {
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
                return haplo_fail(result, cds_len_out,
                                  DUCKVEP_HAPLOTYPE_EDIT_ORDER);
            }
            suffix_len = src_cursor - ref_end;
            if (suffix_len > dst_cursor ||
                (size_t)e->alt_len > dst_cursor - suffix_len) {
                return haplo_fail(result, cds_len_out,
                                  DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
            }
            while (src_cursor > ref_end) {
                char b = haplo_norm_cds_base(ref_cds[--src_cursor]);
                if (b == '\0') {
                    return haplo_fail(result, cds_len_out,
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
            return haplo_fail(result, cds_len_out,
                              DUCKVEP_HAPLOTYPE_OUT_OF_RANGE);
        }
        while (src_cursor > 0u) {
            char b = haplo_norm_cds_base(ref_cds[--src_cursor]);
            if (b == '\0') {
                return haplo_fail(result, cds_len_out,
                                  DUCKVEP_HAPLOTYPE_INVALID_BASE);
            }
            cds_out[--dst_cursor] = (uint8_t)b;
        }
        if (dst_cursor != 0u) {
            return haplo_fail(result, cds_len_out,
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
