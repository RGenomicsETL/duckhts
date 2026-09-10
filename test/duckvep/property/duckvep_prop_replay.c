#include "duckvep_property.h"

/* ===================================================================== *
 * Haplosaurus-shaped multi-edit CDS mutation + translation.
 *
 * Oracle: rebuild the observed CDS from left to right against original
 * coordinates, appending untouched reference segments and independently
 * oriented alternate alleles. The implementation applies edits in-place from
 * right to left; the oracle never memmoves the working sequence. These tests
 * pin the source-audited Haplosaurus behavior without claiming sample/phase-set
 * aggregation or VEP differential parity yet.
 * ===================================================================== */


char haplo_test_oriented_base(const uint8_t *seq, uint32_t len, uint32_t idx,
                                     int reverse_complement) {
    char b;
    if (seq == NULL || idx >= len) return '\0';
    b = (char)seq[reverse_complement ? (len - 1u - idx) : idx];
    if (b >= 'a' && b <= 'z') b = (char)(b - ('a' - 'A'));
    if (b == 'U') b = 'T';
    if (b != 'A' && b != 'C' && b != 'G' && b != 'T') return '\0';
    return reverse_complement ? coding_test_comp(b) : b;
}

uint8_t haplo_test_variant_from_tx_base(char b, int reverse_complement) {
    return (uint8_t)(reverse_complement ? coding_test_comp(b) : b);
}

int haplo_oracle_rebuild(const uint8_t *ref, size_t ref_len,
                                const duckvep_haplotype_edit_t *edits, size_t edit_count,
                                int8_t transcript_strand,
                                uint8_t *out, size_t out_cap, size_t *out_len,
                                int64_t *length_diff, uint32_t *flags) {
    size_t cursor = 0u;
    size_t n = 0u;
    size_t ei;
    int saw_frameshift = 0;
    int64_t diff_total = 0;
    uint32_t f = 0u;

    for (ei = edit_count; ei > 0u; ei--) {
        const duckvep_haplotype_edit_t *e = &edits[ei - 1u];
        size_t start0 = (size_t)e->cds_start - 1u;
        uint32_t j;
        int reverse = e->variant_strand != transcript_strand;
        int64_t d = (int64_t)e->alt_len - (int64_t)e->ref_len;
        if (e->cds_start == 0u || start0 < cursor || start0 > ref_len ||
            (size_t)e->ref_len > ref_len - start0) return 0;
        if (n + (start0 - cursor) > out_cap) return 0;
        memcpy(out + n, ref + cursor, start0 - cursor);
        n += start0 - cursor;
        for (j = 0u; j < e->ref_len; j++) {
            char expected = haplo_test_oriented_base(e->ref, e->ref_len, j, reverse);
            if (expected == '\0' || expected != (char)ref[start0 + (size_t)j]) return 0;
        }
        if (n + (size_t)e->alt_len > out_cap) return 0;
        for (j = 0u; j < e->alt_len; j++) {
            char b = haplo_test_oriented_base(e->alt, e->alt_len, j, reverse);
            if (b == '\0') return 0;
            out[n++] = (uint8_t)b;
        }
        cursor = start0 + (size_t)e->ref_len;
        diff_total += d;
        if (d != 0) f |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
        if ((d % 3) != 0) saw_frameshift = 1;
    }
    if (n + (ref_len - cursor) > out_cap) return 0;
    memcpy(out + n, ref + cursor, ref_len - cursor);
    n += ref_len - cursor;
    if (n < out_cap) out[n] = (uint8_t)'\0';
    if (saw_frameshift) {
        if ((diff_total % 3) == 0) f |= DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT;
        else f |= DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT;
    }
    *out_len = n;
    *length_diff = diff_total;
    *flags = f;
    return 1;
}

/* Deliberately full-matrix, signed-score, gapped-string oracle. Production
 * instead uses nonnegative costs, a proved band and ungapped output spans. */
static size_t difference_alignment_oracle(const uint8_t *a, size_t n,
    const uint8_t *b, size_t m, int align, uint8_t *oa, uint8_t *ob) {
    int score[8][8], direction[8][8];
    for (size_t i = 0u; i <= n; i++) { score[i][0] = -(int)i; direction[i][0] = 1; }
    for (size_t j = 0u; j <= m; j++) { score[0][j] = -(int)j; direction[0][j] = -1; }
    for (size_t i = 1u; i <= n; i++) for (size_t j = 1u; j <= m; j++) {
        int sub = score[i - 1u][j - 1u] + (a[i - 1u] == b[j - 1u] ? 1 : -1);
        int del = score[i][j - 1u] - 1, ins = score[i - 1u][j] - 1;
        if (sub > del && sub > ins) { score[i][j] = sub; direction[i][j] = 0; }
        else if (del > ins) { score[i][j] = del; direction[i][j] = -1; }
        else { score[i][j] = ins; direction[i][j] = 1; }
    }
    size_t length = 0u, i = n, j = m;
    while (i || j) {
        int step = align ? direction[i][j] : i == j ? 0 : i > j ? 1 : -1;
        oa[length] = step == -1 ? '-' : a[--i];
        ob[length++] = step == 1 ? '-' : b[--j];
    }
    for (size_t k = 0u; k < length / 2u; k++) {
        uint8_t swap = oa[k]; oa[k] = oa[length - 1u - k]; oa[length - 1u - k] = swap;
        swap = ob[k]; ob[k] = ob[length - 1u - k]; ob[length - 1u - k] = swap;
    }
    return length;
}

TEST haplotype_differences_match_full_matrix_exhaustively(void) {
    uint8_t a[7], b[7], oa[14], ob[14], trace[64];
    uint64_t scores[16];
    duckvep_sequence_difference_t actual[14], expected[14];
    duckvep_sequence_diff_scratch_t scratch = {scores, 16u, trace, sizeof(trace)};
    for (size_t n = 0u; n <= 6u; n++) for (size_t m = 0u; m <= 6u; m++) {
        for (unsigned x = 0u; x < (1u << n); x++) for (unsigned y = 0u; y < (1u << m); y++) {
            for (size_t i = 0u; i < n; i++) a[i] = (x >> i) & 1u ? 'A' : 'C';
            for (size_t j = 0u; j < m; j++) b[j] = (y >> j) & 1u ? 'A' : 'C';
            for (int align = 0; align <= 1; align++) {
                size_t length = difference_alignment_oracle(a, n, b, m, align, oa, ob);
                size_t count = 0u, ri = 0u, ai = 0u;
                for (size_t col = 0u; col < length;) {
                    if (oa[col] == ob[col]) { ri++; ai++; col++; continue; }
                    size_t begin = col, rn = 0u, an = 0u;
                    int rgap = oa[col] == '-', agap = ob[col] == '-';
                    do {
                        rn += oa[col] != '-'; an += ob[col] != '-'; col++;
                    } while (col < length && oa[col] != ob[col] &&
                        (oa[col] == '-') == rgap && (ob[col] == '-') == agap);
                    expected[count++] = (duckvep_sequence_difference_t){ri, ai, rn, an, begin};
                    ri += rn; ai += an;
                }
                duckvep_sequence_diff_result_t result;
                ASSERT_EQ(DUCKVEP_SEQUENCE_DIFF_OK, duckvep_sequence_differences(
                    a, n, b, m, align, &scratch, actual, 14u, &result));
                ASSERT_EQ(length, result.alignment_length); ASSERT_EQ(count, result.count);
                ASSERT_EQ(0, memcmp(actual, expected, count * sizeof(*actual)));
                if (count) {
                    memset(actual, 0xa5, sizeof(actual));
                    ASSERT_EQ(DUCKVEP_SEQUENCE_DIFF_OUTPUT_FULL, duckvep_sequence_differences(
                        a, n, b, m, align, &scratch, actual, count - 1u, &result));
                    ASSERT_EQ(count, result.count);
                    const uint8_t *bytes = (const uint8_t *)actual;
                    for (size_t k = 0u; k < sizeof(actual); k++) ASSERT_EQ(0xa5u, bytes[k]);
                }
            }
        }
    }
    PASS();
}

TEST haplotype_differences_bound_alignment_work_and_report_limits(void) {
    uint8_t ref[10000], alt[10001], trace[30003];
    uint64_t scores[20004];
    memset(ref, 'A', sizeof(ref)); memset(alt, 'A', sizeof(alt)); alt[5000] = 'C';
    duckvep_sequence_diff_scratch_t scratch = {scores, 20004u, trace, sizeof(trace)};
    duckvep_sequence_difference_t differences[2];
    duckvep_sequence_diff_result_t result;
    ASSERT_EQ(DUCKVEP_SEQUENCE_DIFF_OK, duckvep_sequence_differences(ref, sizeof(ref),
        alt, sizeof(alt), 1, &scratch, differences, 2u, &result));
    ASSERT_EQ(30003u, result.trace_cells); ASSERT_EQ(1u, result.count);
    ASSERT_EQ(5000u, differences[0].ref_start0); ASSERT_EQ(5000u, differences[0].alt_start0);
    ASSERT_EQ(0u, differences[0].ref_length); ASSERT_EQ(1u, differences[0].alt_length);
    scratch.trace_capacity--;
    ASSERT_EQ(DUCKVEP_SEQUENCE_DIFF_TRACE_FULL, duckvep_sequence_differences(ref, sizeof(ref),
        alt, sizeof(alt), 1, &scratch, differences, 2u, &result));
    ASSERT_EQ(30003u, result.trace_cells);
    scratch.trace_capacity++; scratch.score_capacity--;
    ASSERT_EQ(DUCKVEP_SEQUENCE_DIFF_SCORE_FULL, duckvep_sequence_differences(ref, sizeof(ref),
        alt, sizeof(alt), 1, &scratch, differences, 2u, &result));
    ASSERT_EQ(DUCKVEP_SEQUENCE_DIFF_INVALID_ARG, duckvep_sequence_differences(ref, SIZE_MAX,
        alt, sizeof(alt), 1, &scratch, differences, 2u, &result));
    ASSERT_EQ(DUCKVEP_SEQUENCE_DIFF_INVALID_ARG, duckvep_sequence_differences((const uint8_t *)"-", 1u,
        alt, sizeof(alt), 1, &scratch, differences, 2u, &result));
    PASS();
}

TEST haplotype_full_translation_matches_every_supported_codon_table(void) {
    const char alphabet[] = "ACGTUNacgtun";
    uint8_t cds[] = {'A', 'T', 'G', 'T', 'G', 'C', 'T', 'A', 'A', 'G', 'C', 'C', 'N'};
    uint8_t peptide[8];
    duckvep_translation_t result;
    for (unsigned table = 0u; table < 32u; table++) {
        if (!duckvep_codon_table_supported((duckvep_codon_table_t)table)) {
            ASSERT_EQ(DUCKVEP_TRANSLATION_INVALID_ARG,
                duckvep_translate_cds(cds, sizeof(cds), (duckvep_codon_table_t)table,
                    peptide, sizeof(peptide), &result));
            continue;
        }
        for (size_t a = 0u; a < sizeof(alphabet) - 1u; a++) {
            for (size_t b = 0u; b < sizeof(alphabet) - 1u; b++) {
                for (size_t c = 0u; c < sizeof(alphabet) - 1u; c++) {
                    cds[3] = (uint8_t)alphabet[a]; cds[4] = (uint8_t)alphabet[b];
                    cds[5] = (uint8_t)alphabet[c];
                    for (unsigned tail = 0u; tail < 2u; tail++) {
                        cds[12] = tail ? 'n' : 'u';
                        size_t first_stop = 0u;
                        uint8_t expected[4];
                        for (size_t i = 0u; i < 4u; i++) {
                            expected[i] = (uint8_t)duckvep_translate_codon((const char *)cds + i * 3u,
                                (duckvep_codon_table_t)table);
                            if (expected[i] == '*' && !first_stop) first_stop = i + 1u;
                        }
                        ASSERT_EQ(DUCKVEP_TRANSLATION_OK,
                            duckvep_translate_cds(cds, sizeof(cds), (duckvep_codon_table_t)table,
                                peptide, 5u, &result));
                        ASSERT_EQ(4u, result.length);
                        ASSERT_EQ(0, memcmp(peptide, expected, sizeof(expected)));
                        ASSERT_EQ(0u, peptide[4]);
                        ASSERT_EQ(first_stop, result.first_stop_position1);
                        uint32_t scanned_stop;
                        ASSERT(duckvep_cds_first_stop_position1(cds, sizeof(cds),
                            (duckvep_codon_table_t)table, &scanned_stop));
                        ASSERT_EQ(first_stop, scanned_stop);
                        ASSERT_EQ(!tail && !memchr(cds, 'N', 12u) && !memchr(cds, 'n', 12u),
                            result.unambiguous);
                    }
                }
            }
        }
    }
    /* A complete output allocation is checked before writing any peptide. */
    memset(peptide, 0xa5, sizeof(peptide));
    ASSERT_EQ(DUCKVEP_TRANSLATION_BUFFER_TOO_SMALL,
        duckvep_translate_cds(cds, sizeof(cds), STD, peptide, 4u, &result));
    for (size_t i = 0u; i < sizeof(peptide); i++) ASSERT_EQ(0xa5u, peptide[i]);
    ASSERT_EQ(0u, result.length); ASSERT_EQ(0u, result.first_stop_position1);
    ASSERT_EQ(0u, result.unambiguous);
    /* SIZE_MAX bytes exceed the address range before metadata may be cleared.
     * The ordinary output-capacity failure is checked separately above. */
    ASSERT_EQ(DUCKVEP_TRANSLATION_INVALID_ARG,
        duckvep_translate_cds(cds, SIZE_MAX, STD, peptide, sizeof(peptide), &result));
    ASSERT_EQ(DUCKVEP_TRANSLATION_INVALID_ARG,
        duckvep_translate_cds(cds, sizeof(cds), STD, cds, sizeof(cds), &result));
    ASSERT_EQ(DUCKVEP_TRANSLATION_INVALID_ARG,
        duckvep_translate_cds(cds, sizeof(cds), STD, peptide, sizeof(peptide), NULL));
    ASSERT_EQ(DUCKVEP_TRANSLATION_INVALID_ARG,
        duckvep_translate_cds(NULL, sizeof(cds), STD, peptide, sizeof(peptide), &result));
    ASSERT_EQ(DUCKVEP_TRANSLATION_INVALID_ARG,
        duckvep_translate_cds(cds, sizeof(cds), STD, NULL, sizeof(peptide), &result));
    ASSERT_EQ(DUCKVEP_TRANSLATION_INVALID_ARG,
        duckvep_translate_cds(cds, sizeof(cds), STD, peptide, SIZE_MAX, &result));
    /* Invalid bytes after an internal stop and in the trailing partial codon
     * must be rejected too. No valid result may describe a partial translation. */
    memcpy(cds, "ATGTAAGCCGCCT", sizeof(cds));
    for (size_t at = 0u; at < sizeof(cds); at++) {
        uint8_t before = cds[at]; cds[at] = '?';
        ASSERT_EQ(DUCKVEP_TRANSLATION_INVALID_BASE,
            duckvep_translate_cds(cds, sizeof(cds), STD, peptide, sizeof(peptide), &result));
        ASSERT_EQ(0u, result.length); ASSERT_EQ(0u, result.first_stop_position1);
        cds[at] = before;
    }
    ASSERT_EQ(DUCKVEP_TRANSLATION_OK,
        duckvep_translate_cds(cds, 0u, STD, peptide, 1u, &result));
    ASSERT_EQ(0u, result.length); ASSERT_EQ(0u, result.first_stop_position1);
    ASSERT_EQ(1u, result.unambiguous); ASSERT_EQ(0u, peptide[0]);
    PASS();
}

TEST haplotype_consensus_translation_matches_all_n_expansions(void) {
    const char alphabet[] = "ACGTNacgtn", literal[] = "ACGT";
    uint8_t cds[4], peptide[2];
    duckvep_translation_t result;
    uint8_t reference[3], coding_peptide[3];
    duckvep_translation_t coding;
    for (unsigned table = 1u; table < 32u; table++) {
        if (!duckvep_codon_table_supported((duckvep_codon_table_t)table)) continue;
        for (size_t a = 0u; a < 10u; a++) for (size_t b = 0u; b < 10u; b++)
            for (size_t c = 0u; c < 10u; c++) {
                cds[0] = (uint8_t)alphabet[a]; cds[1] = (uint8_t)alphabet[b];
                cds[2] = (uint8_t)alphabet[c]; cds[3] = 'N';
                unsigned char amino_acids[128] = {0};
                int any_start = 0;
                for (size_t x = 0u; x < 4u; x++) for (size_t y = 0u; y < 4u; y++)
                    for (size_t z = 0u; z < 4u; z++) {
                        char expansion[4] = {literal[x], literal[y], literal[z], 0};
                        int match = 1;
                        for (size_t i = 0u; i < 3u; i++) {
                            char base = (char)(cds[i] & 0xdfu);
                            if (base != 'N' && base != expansion[i]) match = 0;
                        }
                        if (match) {
                            amino_acids[(unsigned char)duckvep_translate_codon(
                                expansion, (duckvep_codon_table_t)table)] = 1u;
                            any_start |= duckvep_codon_is_start((const uint8_t *)expansion,
                                (duckvep_codon_table_t)table);
                        }
                    }
                unsigned count = 0u; uint8_t expected = 'X';
                for (unsigned i = 0u; i < 128u; i++) if (amino_acids[i]) {
                    count++; expected = (uint8_t)i;
                }
                if (count == 2u && amino_acids['D'] && amino_acids['N']) expected = 'B';
                else if (count == 2u && amino_acids['E'] && amino_acids['Q']) expected = 'Z';
                else if (count > 1u) expected = 'X';
                ASSERT_EQ(expected, duckvep_translate_codon((const char *)cds,
                    (duckvep_codon_table_t)table));
                char prepared[3];
                for (size_t i = 0u; i < 3u; i++) prepared[i] = (char)(cds[i] & 0xdfu);
                duckvep_codon_result_t change = duckvep_codon_change_prepared(
                    prepared, "ATG", (duckvep_codon_table_t)table);
                ASSERT_EQ(expected, change.aa_ref);
                ASSERT_EQ('M', change.aa_alt);
                uint32_t expected_change = expected == 'X' ? DUCKVEP_CODON_INVALID
                    : expected == 'M' ? DUCKVEP_CODON_SYNONYMOUS
                    : expected == '*' ? DUCKVEP_CODON_STOP_LOST : DUCKVEP_CODON_MISSENSE;
                ASSERT_EQ(expected_change, change.change);
                for (size_t length = 3u; length <= 4u; length++) {
                    ASSERT_EQ(DUCKVEP_TRANSLATION_OK, duckvep_translate_cds(cds, length,
                        (duckvep_codon_table_t)table,
                        peptide, sizeof(peptide), &result));
                    ASSERT_EQ(expected, peptide[0]); ASSERT_EQ(0u, peptide[1]);
                    ASSERT_EQ(1u, result.length);
                    ASSERT_EQ(expected == '*', result.first_stop_position1);
                    ASSERT_EQ(length == 3u && a % 5u != 4u && b % 5u != 4u && c % 5u != 4u,
                        result.unambiguous);
                    uint32_t scanned_stop;
                    ASSERT(duckvep_cds_first_stop_position1(cds, length,
                        (duckvep_codon_table_t)table, &scanned_stop));
                    ASSERT_EQ(expected == '*', scanned_stop);
                    size_t reference_length;
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, duckvep_haplotype_reference_proteins(
                        cds, length, (duckvep_codon_table_t)table, NULL, NULL, 0u,
                        reference, coding_peptide, sizeof(reference), &reference_length, &coding));
                    ASSERT_EQ(expected, coding_peptide[0]);
                    ASSERT_EQ(0u, coding_peptide[1]);
                    ASSERT_EQ(1u, coding.length);
                    ASSERT_EQ(expected == '*', coding.first_stop_position1);
                    ASSERT_EQ(result.unambiguous, coding.unambiguous);
                    size_t expected_length = expected == '*' ? 0u : 1u;
                    uint8_t curated[3] = {0u};
                    if (expected_length) curated[0] = any_start ? 'M' : expected;
                    if (memcmp(cds + length - 3u, "TAA", 3u) == 0 ||
                        memcmp(cds + length - 3u, "TAG", 3u) == 0 ||
                        memcmp(cds + length - 3u, "TGA", 3u) == 0) {
                        curated[expected_length++] = '*';
                    }
                    ASSERT_EQ(expected_length, reference_length);
                    ASSERT_MEM_EQ(curated, reference, expected_length + 1u);
                }
            }
    }
    ASSERT_EQ(DUCKVEP_TRANSLATION_INVALID_ARG, duckvep_translate_cds(cds, 3u,
        (duckvep_codon_table_t)8, peptide, sizeof(peptide), &result));
    ASSERT_EQ(0u, result.length);
    PASS();
}

TEST reference_translation_views_validate_distinct_bounded_storage(void) {
    uint8_t cds[] = "ATGGCNTGANN", reference[8], coding_peptide[8];
    duckvep_translation_t result;
    size_t reference_length;
    for (size_t length = 0u; length <= 10u; length++) {
        memset(reference, 0xa5, sizeof(reference));
        memset(coding_peptide, 0xa5, sizeof(coding_peptide));
        size_t codons = length / 3u;
        ASSERT_EQ(length < 3u ? DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE : DUCKVEP_HAPLOTYPE_OK,
            duckvep_haplotype_reference_proteins(cds, length, STD, NULL, NULL, 0u,
                reference, coding_peptide, codons + 2u, &reference_length, &result));
        size_t expected_reference_length = length == 10u ? 2u : codons;
        ASSERT_EQ(expected_reference_length, reference_length);
        ASSERT_MEM_EQ("MA*", reference, expected_reference_length);
        ASSERT_MEM_EQ("MA*", coding_peptide, codons);
        ASSERT_EQ(0u, reference[expected_reference_length]);
        ASSERT_EQ(0u, coding_peptide[codons]);
        ASSERT_EQ(0xa5u, reference[codons + 2u]);
        ASSERT_EQ(0xa5u, coding_peptide[codons + 1u]);
        ASSERT_EQ(codons, result.length);
        ASSERT_EQ(codons == 3u ? 3u : 0u, result.first_stop_position1);
        ASSERT_EQ(length < 6u, result.unambiguous);
    }
    for (unsigned invalid = 0u; invalid < 7u; invalid++) {
        memset(reference, 0xa5, sizeof(reference));
        memset(coding_peptide, 0xa5, sizeof(coding_peptide));
        uint8_t *raw_output = coding_peptide;
        if (invalid == 4u) raw_output = NULL;
        if (invalid == 5u) raw_output = cds;
        if (invalid == 6u) raw_output = reference + 1u;
        ASSERT_EQ(invalid == 0u ? DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL : DUCKVEP_HAPLOTYPE_INVALID_ARG,
            duckvep_haplotype_reference_proteins(invalid == 1u ? NULL : cds, 10u,
                invalid == 2u ? (duckvep_codon_table_t)8 : STD,
                NULL, NULL, 0u, invalid == 3u ? NULL : reference, raw_output,
                invalid == 0u ? 4u : sizeof(reference), &reference_length, &result));
        ASSERT_EQ(0u, result.length); ASSERT_EQ(0u, result.first_stop_position1);
        ASSERT_EQ(0u, result.unambiguous);
        ASSERT_EQ(0u, reference_length);
        for (size_t i = 0u; i < sizeof(reference); i++) {
            ASSERT_EQ(0xa5u, reference[i]);
            ASSERT_EQ(0xa5u, coding_peptide[i]);
        }
        ASSERT_STR_EQ("ATGGCNTGANN", (const char *)cds);
    }
    for (size_t at = 0u; at < 10u; at++) {
        uint8_t before = cds[at];
        cds[at] = '?';
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_BASE, duckvep_haplotype_reference_proteins(cds, 10u,
            STD, NULL, NULL, 0u, reference, coding_peptide, sizeof(reference),
            &reference_length, &result));
        ASSERT_EQ(0u, result.length); ASSERT_EQ(0u, result.first_stop_position1);
        ASSERT_EQ(0u, result.unambiguous);
        cds[at] = before;
    }
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG, duckvep_haplotype_reference_proteins(cds, 10u,
        STD, NULL, NULL, 0u, reference, coding_peptide, sizeof(reference), &reference_length, NULL));
    PASS();
}

TEST haplotype_reference_protein_applies_ensembl_rules_with_checked_storage(void) {
    static const struct { const char *cds; unsigned table; const char *protein; } cases[] = {
        {"CTGGCCTAA", 1u, "MA*"}, {"ATGTGAGCCTAA", 1u, "M*A*"},
        {"NNNGCCTAA", 1u, "MA*"}, {"ATGGCNTAA", 1u, "MA*"},
        {"GTGGCCTAA", 1u, "VA*"}, {"ATGGCCTGA", 2u, "MAW*"},
        {"ATGGCCAGA", 2u, "MA"}, {"ATGGCCTAAA", 1u, "MA"},
        {"atggcctaa", 1u, "MA"}, {"TAA", 1u, "*"}, {"TAACC", 1u, ""}
    };
    uint8_t peptide[32], coding_peptide[32]; size_t length;
    duckvep_translation_t coding_translation;
    for (size_t i = 0u; i < sizeof(cases) / sizeof(cases[0]); i++) {
        size_t n = strlen(cases[i].cds);
        memset(peptide, 0xa5, sizeof(peptide));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, duckvep_haplotype_reference_proteins(
            (const uint8_t *)cases[i].cds, n, (duckvep_codon_table_t)cases[i].table,
            NULL, NULL, 0u, peptide, coding_peptide, n / 3u + 2u, &length, &coding_translation));
        ASSERT_EQ(strlen(cases[i].protein), length);
        ASSERT_STR_EQ(cases[i].protein, (const char *)peptide);
        ASSERT_EQ(0xa5u, peptide[n / 3u + 2u]);
        ASSERT_EQ(n / 3u, coding_translation.length);
        size_t first_stop = 0u;
        for (size_t j = 0u; j < n / 3u; j++) {
            char aa = duckvep_translate_codon(cases[i].cds + j * 3u,
                (duckvep_codon_table_t)cases[i].table);
            ASSERT_EQ((uint8_t)aa, coding_peptide[j]);
            if (aa == '*' && !first_stop) first_stop = j + 1u;
        }
        ASSERT_EQ(first_stop, coding_translation.first_stop_position1);
        ASSERT_EQ(0u, coding_peptide[n / 3u]);
    }
    uint8_t cds[] = "GTGTGAGCCTAA";
    uint32_t positions[] = {1u, 2u, 3u, 4u};
    uint8_t alternates[] = "MUVW";
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, duckvep_haplotype_reference_proteins(cds, 12u, STD,
        positions, alternates, 4u, peptide, coding_peptide, 6u, &length, &coding_translation));
    ASSERT_EQ(5u, length); ASSERT_STR_EQ("MUVW*", (const char *)peptide);
    ASSERT_STR_EQ("V*A*", (const char *)coding_peptide);
    ASSERT_EQ(2u, coding_translation.first_stop_position1);
    memset(peptide, 0xa5, sizeof(peptide));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL, duckvep_haplotype_reference_proteins(cds, 12u, STD,
        positions, alternates, 4u, peptide, coding_peptide, 5u, &length, &coding_translation));
    ASSERT_EQ(0u, length);
    for (size_t i = 0u; i < sizeof(peptide); i++) ASSERT_EQ(0xa5u, peptide[i]);
    positions[1] = 1u;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OUT_OF_RANGE, duckvep_haplotype_reference_proteins(cds, 12u, STD,
        positions, alternates, 4u, peptide, coding_peptide, sizeof(peptide), &length, &coding_translation));
    positions[1] = 2u; positions[3] = 5u;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OUT_OF_RANGE, duckvep_haplotype_reference_proteins(cds, 12u, STD,
        positions, alternates, 4u, peptide, coding_peptide, sizeof(peptide), &length, &coding_translation));
    positions[3] = 4u; alternates[0] = '?';
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG, duckvep_haplotype_reference_proteins(cds, 12u, STD,
        positions, alternates, 4u, peptide, coding_peptide, sizeof(peptide), &length, &coding_translation));
    for (size_t i = 0u; i < sizeof(peptide); i++) ASSERT_EQ(0xa5u, peptide[i]);
    alternates[0] = 'M';
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG, duckvep_haplotype_reference_proteins(cds, 12u, STD,
        positions, alternates, 4u, cds, coding_peptide, sizeof(cds), &length, &coding_translation));
    ASSERT_STR_EQ("GTGTGAGCCTAA", (const char *)cds);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG, duckvep_haplotype_reference_proteins(cds, 12u, STD,
        positions, alternates, 4u, (uint8_t *)positions, coding_peptide, sizeof(positions), &length,
        &coding_translation));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG, duckvep_haplotype_reference_proteins(cds, 12u, STD,
        positions, peptide, 4u, peptide, coding_peptide, sizeof(peptide), &length, &coding_translation));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG, duckvep_haplotype_reference_proteins(cds, SIZE_MAX, STD,
        NULL, NULL, 0u, peptide, coding_peptide, sizeof(peptide), &length, &coding_translation));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG, duckvep_haplotype_reference_proteins(cds, 12u, STD,
        NULL, NULL, SIZE_MAX, peptide, coding_peptide, sizeof(peptide), &length, &coding_translation));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE, duckvep_haplotype_reference_proteins(cds, 2u, STD,
        NULL, NULL, 0u, peptide, coding_peptide, sizeof(peptide), &length, &coding_translation));
    ASSERT_EQ(0u, length);
    ASSERT_EQ(0u, coding_translation.length);
    ASSERT_EQ(1u, coding_translation.unambiguous);
    ASSERT_EQ(0u, coding_peptide[0]);
    cds[10] = '?';
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_BASE, duckvep_haplotype_reference_proteins(cds, 12u, STD,
        NULL, NULL, 0u, peptide, coding_peptide, sizeof(peptide), &length, &coding_translation));
    ASSERT_EQ(0u, length);
    ASSERT_EQ(0, duckvep_codon_is_start((const uint8_t *)"?TG", STD));
    ASSERT_EQ(0, duckvep_codon_is_start((const uint8_t *)"ATG", (duckvep_codon_table_t)8));
    PASS();
}

TEST reference_translation_result_aliases_preserve_all_storage(void) {
    union aligned_reference_storage {
        duckvep_translation_t alignment;
        uint8_t bytes[128];
    } storage[5], before[5];
    const size_t base = 2u * sizeof(duckvep_translation_t);
    const size_t alignment = _Alignof(duckvep_translation_t);
    ASSERT(base + sizeof(duckvep_translation_t) + alignment <= sizeof(storage[0].bytes));
    for (unsigned span = 0u; span < 5u; span++) {
        for (unsigned placement = 0u; placement < (span < 3u ? 3u : 2u); placement++) {
            memset(storage, 0xa5, sizeof(storage));
            uint8_t *cds = storage[0].bytes + base;
            uint8_t *peptide = storage[1].bytes + base;
            uint8_t *coding = storage[2].bytes + base;
            uint32_t *position = (uint32_t *)(void *)(storage[3].bytes + base);
            uint8_t *alternate = storage[4].bytes + base;
            memcpy(cds, "CTGGCNTAA", 9u);
            uint32_t position1 = 2u;
            memcpy(position, &position1, sizeof(position1));
            *alternate = 'W';
            memcpy(before, storage, sizeof(storage));
            size_t offset = placement == 0u ? base - alignment
                : placement == 1u ? base : base + alignment;
            duckvep_translation_t *aliased = (duckvep_translation_t *)(void *)(storage[span].bytes + offset);
            if (span < 3u) {
                size_t length = 99u;
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG, duckvep_haplotype_reference_proteins(
                    cds, 9u, STD, NULL, NULL, 0u, peptide, coding, 16u, &length, aliased));
                ASSERT_EQ(99u, length);
                ASSERT_MEM_EQ(before, storage, sizeof(storage));
            }
            if (span < 2u) {
                ASSERT_EQ(DUCKVEP_TRANSLATION_INVALID_ARG, duckvep_translate_cds(
                    cds, 9u, STD, peptide, 16u, aliased));
                ASSERT_MEM_EQ(before, storage, sizeof(storage));
            }
            size_t length = 99u;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG, duckvep_haplotype_reference_proteins(
                cds, 9u, STD, position, alternate, 1u, peptide, coding, 16u, &length, aliased));
            ASSERT_EQ(99u, length);
            ASSERT_MEM_EQ(before, storage, sizeof(storage));
            /* Length storage has a smaller span: test its exact overlap. */
            size_t *aliased_length = (size_t *)(void *)(storage[span].bytes + base);
            duckvep_translation_t result, result_before;
            memset(&result, 0xa5, sizeof(result));
            memcpy(&result_before, &result, sizeof(result));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG, duckvep_haplotype_reference_proteins(
                cds, 9u, STD, position, alternate, 1u, peptide, coding, 16u, aliased_length, &result));
            ASSERT_MEM_EQ(before, storage, sizeof(storage));
            ASSERT_MEM_EQ(&result_before, &result, sizeof(result));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG, duckvep_haplotype_reference_proteins(
                cds, 9u, STD, position, alternate, 1u, peptide, coding, 16u, &result.length, &result));
            ASSERT_MEM_EQ(before, storage, sizeof(storage));
            ASSERT_MEM_EQ(&result_before, &result, sizeof(result));
        }
    }
    PASS();
}

TEST haplotype_apply_and_translate_known_cases(void) {
    duckvep_haplotype_edit_t edits[2];
    duckvep_haplotype_result_t r;
    uint8_t cds[64];
    uint8_t protein[32];
    size_t cds_len;
    duckvep_translation_t translated;

    memset(edits, 0, sizeof edits);
    edits[0].cds_start = 2u; edits[0].ref_len = 1u; edits[0].ref = (const uint8_t *)"T";
    edits[0].alt_len = 1u; edits[0].alt = (const uint8_t *)"G"; edits[0].variant_strand = (int8_t)1;
    edits[1].cds_start = 1u; edits[1].ref_len = 1u; edits[1].ref = (const uint8_t *)"A";
    edits[1].alt_len = 1u; edits[1].alt = (const uint8_t *)"T"; edits[1].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATGAAATAA", 9u,
                                                edits, 2u, (int8_t)1, cds, sizeof cds, &cds_len, &r));
    ASSERT_EQ(9u, cds_len);
    ASSERT_STR_EQ("TGGAAATAA", (const char *)cds);
    ASSERT_EQ(DUCKVEP_TRANSLATION_OK,
              duckvep_translate_cds(cds, cds_len, STD, protein, sizeof protein, &translated));
    ASSERT_STR_EQ("WK*", (const char *)protein);
    ASSERT_EQ(3u, translated.first_stop_position1);

    /* Negative transcript strand: genomic C>A orients to transcript G>T. */
    memset(edits, 0, sizeof edits);
    edits[0].cds_start = 4u; edits[0].ref_len = 1u; edits[0].ref = (const uint8_t *)"C";
    edits[0].alt_len = 1u; edits[0].alt = (const uint8_t *)"A"; edits[0].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATGGCTTAA", 9u,
                                                edits, 1u, (int8_t)-1, cds, sizeof cds, &cds_len, &r));
    ASSERT_STR_EQ("ATGTCTTAA", (const char *)cds);
    ASSERT_EQ(DUCKVEP_TRANSLATION_OK,
              duckvep_translate_cds(cds, cds_len, STD, protein, sizeof protein, &translated));
    ASSERT_STR_EQ("MS*", (const char *)protein);

    /* Two non-triplet edits with net length diff 0 are a resolved frameshift. */
    memset(edits, 0, sizeof edits);
    edits[0].cds_start = 10u; edits[0].ref_len = 1u; edits[0].ref = (const uint8_t *)"G";
    edits[0].alt_len = 0u; edits[0].alt = NULL; edits[0].variant_strand = (int8_t)1;
    edits[1].cds_start = 4u; edits[1].ref_len = 0u; edits[1].ref = NULL;
    edits[1].alt_len = 1u; edits[1].alt = (const uint8_t *)"C"; edits[1].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATGAAACCCGGGTAA", 15u,
                                                edits, 2u, (int8_t)1, cds, sizeof cds, &cds_len, &r));
    ASSERT_STR_EQ("ATGCAAACCCGGTAA", (const char *)cds);
    ASSERT_EQ(0, r.length_diff);
    ASSERT((r.flags & DUCKVEP_HAPLOTYPE_FLAG_INDEL) != 0u);
    ASSERT((r.flags & DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT) != 0u);
    ASSERT((r.flags & DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT) == 0u);

    /* The visible protein still ends at the first stop; later residues remain
     * available to coding facts without translating the same sequence again. */
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATGAAATAACCC", 12u,
                                                NULL, 0u, (int8_t)1, cds, sizeof cds, &cds_len, &r));
    {
        duckvep_translation_t fresh_result;
        ASSERT_EQ(DUCKVEP_TRANSLATION_OK,
                  duckvep_translate_cds(cds, cds_len, STD, protein, sizeof protein, &fresh_result));
        ASSERT_EQ(3u, fresh_result.first_stop_position1);
        ASSERT_EQ(0, memcmp(protein, "MK*", fresh_result.first_stop_position1));
        ASSERT_EQ(4u, fresh_result.length);
        ASSERT_STR_EQ("MK*P", (const char *)protein);
        ASSERT(fresh_result.first_stop_position1 < fresh_result.length);
    }

    /* Contract checks: descending order and reference validation are explicit. */
    edits[0].cds_start = 1u; edits[0].ref_len = 1u; edits[0].ref = (const uint8_t *)"A";
    edits[0].alt_len = 1u; edits[0].alt = (const uint8_t *)"C"; edits[0].variant_strand = (int8_t)1;
    edits[1].cds_start = 2u; edits[1].ref_len = 1u; edits[1].ref = (const uint8_t *)"T";
    edits[1].alt_len = 1u; edits[1].alt = (const uint8_t *)"G"; edits[1].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_EDIT_ORDER,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATG", 3u, edits, 2u,
                                                (int8_t)1, cds, sizeof cds, &cds_len, &r));
    edits[0].cds_start = 1u; edits[0].ref = (const uint8_t *)"C";
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_REF_MISMATCH,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATG", 3u, edits, 1u,
                                                (int8_t)1, cds, sizeof cds, &cds_len, &r));
    edits[0].cds_start = 4u; edits[0].ref_len = 0u; edits[0].ref = NULL;
    edits[0].alt_len = 1u; edits[0].alt = (const uint8_t *)"A"; edits[0].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATG", 3u, edits, 1u,
                                                (int8_t)1, cds, sizeof cds, &cds_len, &r));
    ASSERT_STR_EQ("ATGA", (const char *)cds);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATG", 3u, edits, 1u,
                                                (int8_t)1, cds, 3u, &cds_len, &r));

    /* A later lower-coordinate deletion restores the final CDS length.
     * Only final-length scratch is needed; in-place mutation is rejected. */
    {
        uint8_t alias_small[17];
        uint8_t reference[17];
        uint8_t before[17];

        memset(edits, 0, sizeof edits);
        edits[0].cds_start = 11u;
        edits[0].alt_len = 4u;
        edits[0].alt = (const uint8_t *)"AAAA";
        edits[0].variant_strand = (int8_t)1;
        edits[1].cds_start = 2u;
        edits[1].ref_len = 4u;
        edits[1].ref = (const uint8_t *)"CGTA";
        edits[1].variant_strand = (int8_t)1;

        memset(alias_small, 0xa5, sizeof alias_small);
        memcpy(alias_small, "ACGTACGTACGT", 12u);
        memcpy(before, alias_small, sizeof before);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG,
                  duckvep_haplotype_apply_cds_edits(
                      alias_small, 12u, edits, 2u, (int8_t)1,
                      alias_small, 12u, &cds_len, &r));
        ASSERT_EQ(0, memcmp(alias_small, before, sizeof before));

        memset(reference, 0, sizeof reference);
        memcpy(reference, "ACGTACGTACGT", 12u);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                  duckvep_haplotype_apply_cds_edits(
                      reference, 12u, edits, 2u, (int8_t)1,
                      alias_small, 12u, &cds_len, &r));
        ASSERT_EQ(12u, cds_len);
        ASSERT_EQ(0, memcmp("ACGTACAAAAGT", alias_small, 12u));
        ASSERT_STR_EQ("ACGTACGTACGT", (const char *)reference);
    }

    edits[0].cds_start = 2u; edits[0].ref_len = 2u; edits[0].ref = (const uint8_t *)"TG";
    edits[0].alt_len = 1u; edits[0].alt = (const uint8_t *)"C"; edits[0].variant_strand = (int8_t)1;
    edits[1].cds_start = 3u; edits[1].ref_len = 1u; edits[1].ref = (const uint8_t *)"G";
    edits[1].alt_len = 1u; edits[1].alt = (const uint8_t *)"A"; edits[1].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_EDIT_ORDER,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATG", 3u, edits, 2u,
                                                (int8_t)1, cds, sizeof cds, &cds_len, &r));
    PASS();
}

TEST haplotype_edit_geometry_agrees_in_both_orders(void) {
    const uint8_t reference[] = "AAAAAAAAAAAA";
    uint8_t output[32], before[32];
    duckvep_haplotype_edit_t ascending[2], descending[2];
    duckvep_haplotype_block_t blocks[2];
    duckvep_haplotype_result_t result;

    /* Exhaust every pair of in-range replacement spans and interbase sites.
     * Alleles match by construction, so geometry is the only possible error.
     * Include both transcript orientations, distinct scratch and rejected alias,
     * CDS starts/ends, adjacent edits, and coincident zero-length insertions. */
    for (uint32_t left = 1u; left <= sizeof reference; left++) {
        for (uint32_t right = left; right <= sizeof reference; right++) {
            for (uint32_t left_len = 0u; left_len <= sizeof reference - left; left_len++) {
                for (uint32_t right_len = 0u; right_len <= sizeof reference - right; right_len++) {
                    for (int strand = -1; strand <= 1; strand += 2) {
                        ascending[0] = (duckvep_haplotype_edit_t){
                            left, left_len, reference, 1u, (const uint8_t *)"C", (int8_t)strand
                        };
                        ascending[1] = (duckvep_haplotype_edit_t){
                            right, right_len, reference, 1u, (const uint8_t *)"G", (int8_t)strand
                        };
                        descending[0] = ascending[1];
                        descending[1] = ascending[0];
                        /* A point insertion occupies its interbase site for
                         * conflict detection; two unordered edits cannot own it. */
                        uint32_t occupied = left_len ? left_len : 1u;
                        duckvep_haplotype_status_t expected = right - left >= occupied
                            ? DUCKVEP_HAPLOTYPE_OK : DUCKVEP_HAPLOTYPE_EDIT_ORDER;
                        size_t required = 999u;
                        ASSERT_EQ(expected == DUCKVEP_HAPLOTYPE_OK,
                                  duckvep_haplotype_partition(ascending, 2u, blocks, 2u,
                                                              &required) == DUCKVEP_HAPLOTYPE_OK);
                        for (int alias = 0; alias < 2; alias++) {
                            memset(output, 0xa5, sizeof output);
                            if (alias) memcpy(output, reference, sizeof reference);
                            memcpy(before, output, sizeof before);
                            size_t output_len = 999u;
                            memset(&result, 0xa5, sizeof result);
                            duckvep_haplotype_status_t apply_expected = alias
                                ? DUCKVEP_HAPLOTYPE_INVALID_ARG : expected;
                            ASSERT_EQ(apply_expected, duckvep_haplotype_apply_cds_edits(
                                alias ? output : reference, sizeof reference - 1u,
                                descending, 2u, (int8_t)strand, output, sizeof output,
                                &output_len, &result));
                            if (apply_expected != DUCKVEP_HAPLOTYPE_OK) {
                                if (expected == DUCKVEP_HAPLOTYPE_EDIT_ORDER) {
                                    ASSERT_EQ(0u, required);
                                }
                                ASSERT_EQ(0u, output_len);
                                ASSERT_EQ(0u, result.applied_edits);
                                ASSERT_EQ(0u, result.flags);
                                ASSERT_EQ(0, memcmp(before, output, sizeof output));
                            } else {
                                ASSERT_EQ(2u, result.applied_edits);
                                ASSERT_EQ(sizeof reference - 1u - left_len - right_len + 2u,
                                          output_len);
                            }
                        }
                    }
                }
            }
        }
    }
    PASS();
}

TEST haplotype_apply_rejects_overlapping_inputs(void) {
    const uint8_t reference[] = "AAAAAAAAAAAA";
    uint8_t storage[64], before[64];
    duckvep_haplotype_result_t result;
    size_t length;

    /* Independently enumerate exact overlap, either direction, and adjacent
     * spans for the CDS and each allele. No invalid call may mutate storage. */
    for (int kind = 0; kind < 3; kind++) {
        for (size_t src = 0u; src <= 20u; src++) {
            for (size_t dst = 0u; dst <= 20u; dst++) {
                size_t source_len = kind == 0 ? 12u : 1u;
                int overlap = src < dst + 13u && dst < src + source_len;
                duckvep_haplotype_edit_t edit = {
                    6u, 1u, kind == 1 ? storage + src : reference,
                    1u, kind == 2 ? storage + src : reference, (int8_t)1
                };
                memset(storage, 'A', sizeof storage);
                memcpy(before, storage, sizeof before);
                memset(&result, 0xa5, sizeof result);
                length = 999u;
                ASSERT_EQ(overlap ? DUCKVEP_HAPLOTYPE_INVALID_ARG : DUCKVEP_HAPLOTYPE_OK,
                    duckvep_haplotype_apply_cds_edits(
                        kind == 0 ? storage + src : reference, 12u, &edit, 1u,
                        (int8_t)1, storage + dst, 13u, &length, &result));
                if (overlap) {
                    ASSERT_EQ(0u, length);
                    ASSERT_EQ(0u, result.applied_edits);
                    ASSERT_EQ(0u, result.cds_len);
                    ASSERT_EQ(0u, result.flags);
                    ASSERT_EQ(0, memcmp(storage, before, sizeof storage));
                } else {
                    ASSERT_EQ(12u, length);
                    ASSERT_EQ(1u, result.applied_edits);
                    ASSERT_EQ(0, memcmp(storage + dst, reference, 12u));
                    ASSERT_EQ(0u, storage[dst + 12u]);
                    ASSERT_EQ(0, memcmp(storage + src, before + src, source_len));
                }
            }
        }
    }
    duckvep_haplotype_edit_t edit = {6u, 1u, reference, 1u, reference, (int8_t)1};
    unsigned char edit_before[sizeof edit];
    memcpy(edit_before, &edit, sizeof edit);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG,
        duckvep_haplotype_apply_cds_edits(reference, 12u, &edit, 1u, (int8_t)1,
            (uint8_t *)&edit, sizeof edit, &length, &result));
    ASSERT_EQ(0, memcmp(edit_before, &edit, sizeof edit));
    ASSERT_EQ(0u, length);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG,
        duckvep_haplotype_apply_cds_edits(reference, 12u, &edit,
            SIZE_MAX / sizeof edit + 1u, (int8_t)1, storage, sizeof storage, &length, &result));
    ASSERT_EQ(0u, length);
    PASS();
}

TEST haplotype_partition_known_cases(void) {
    duckvep_haplotype_edit_t edits[2];
    duckvep_haplotype_block_t blocks[2];
    size_t required = 0u;

    memset(edits, 0, sizeof edits);
    edits[0].cds_start = 1u; edits[0].ref_len = 1u; edits[0].ref = (const uint8_t *)"A";
    edits[0].alt_len = 1u; edits[0].alt = (const uint8_t *)"C"; edits[0].variant_strand = (int8_t)1;
    edits[1].cds_start = 3u; edits[1].ref_len = 1u; edits[1].ref = (const uint8_t *)"G";
    edits[1].alt_len = 1u; edits[1].alt = (const uint8_t *)"T"; edits[1].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_partition(edits, 2u, blocks, 2u, &required));
    ASSERT_EQ(1u, required);
    ASSERT_EQ(0u, blocks[0].edit_begin);
    ASSERT_EQ(2u, blocks[0].edit_count);
    ASSERT_EQ(1u, blocks[0].cds_start);
    ASSERT_EQ(3u, blocks[0].ref_len);
    ASSERT_EQ(0, blocks[0].length_diff);
    ASSERT_EQ(0u, blocks[0].flags);

    edits[1].cds_start = 7u;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_partition(edits, 2u, blocks, 2u, &required));
    ASSERT_EQ(2u, required);
    ASSERT_EQ(1u, blocks[0].edit_count);
    ASSERT_EQ(1u, blocks[1].edit_count);

    /* A downstream deletion cannot be flushed while the insertion has left the
     * frame displaced.  Together they restore it. */
    memset(edits, 0, sizeof edits);
    edits[0].cds_start = 4u; edits[0].alt_len = 1u;
    edits[0].alt = (const uint8_t *)"A"; edits[0].variant_strand = (int8_t)1;
    edits[1].cds_start = 10u; edits[1].ref_len = 1u;
    edits[1].ref = (const uint8_t *)"G"; edits[1].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_partition(edits, 2u, blocks, 2u, &required));
    ASSERT_EQ(1u, required);
    ASSERT_EQ(0, blocks[0].length_diff);
    ASSERT((blocks[0].flags & DUCKVEP_HAPLOTYPE_FLAG_INDEL) != 0u);
    ASSERT((blocks[0].flags & DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT) != 0u);
    ASSERT((blocks[0].flags & DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT) == 0u);

    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_partition(edits, 1u, blocks, 2u, &required));
    ASSERT_EQ(1u, required);
    ASSERT_EQ(1, blocks[0].length_diff);
    ASSERT((blocks[0].flags & DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT) != 0u);
    ASSERT((blocks[0].flags & DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT) == 0u);

    /* The sizing pass publishes the required count but never a partial array. */
    edits[0].cds_start = 1u; edits[0].ref_len = 1u;
    edits[0].ref = (const uint8_t *)"A"; edits[0].alt_len = 1u;
    edits[0].alt = (const uint8_t *)"C";
    edits[1].cds_start = 7u; edits[1].ref_len = 1u;
    edits[1].ref = (const uint8_t *)"G"; edits[1].alt_len = 1u;
    edits[1].alt = (const uint8_t *)"T";
    memset(blocks, 0xa5, sizeof blocks);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL,
              duckvep_haplotype_partition(edits, 2u, blocks, 1u, &required));
    ASSERT_EQ(2u, required);
    {
        const unsigned char *raw = (const unsigned char *)blocks;
        size_t i;
        for (i = 0u; i < sizeof blocks; i++) ASSERT_EQ(0xa5u, raw[i]);
    }

    edits[1].cds_start = 1u;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_EDIT_ORDER,
              duckvep_haplotype_partition(edits, 2u, blocks, 2u, &required));
    edits[1].cds_start = 2u;
    edits[0].ref_len = 2u; edits[0].ref = (const uint8_t *)"AT";
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_EDIT_ORDER,
              duckvep_haplotype_partition(edits, 2u, blocks, 2u, &required));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG,
              duckvep_haplotype_partition(NULL, 1u, blocks, 2u, &required));
    PASS();
}

TEST haplotype_partition_spans_track_both_cds_axes(void) {
    const uint8_t *a = (const uint8_t *)"A";
    duckvep_haplotype_edit_t edits[] = {
        {1u, 0u, NULL, 3u, (const uint8_t *)"AAA", 1},
        {4u, 1u, a, 1u, (const uint8_t *)"C", 1},
        {7u, 0u, NULL, 1u, (const uint8_t *)"T", 1},
        {10u, 1u, a, 0u, NULL, 1},
        {16u, 1u, a, 1u, (const uint8_t *)"G", 1}
    };
    duckvep_haplotype_block_t blocks[5];
    size_t count;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, duckvep_haplotype_partition(edits, 5u, blocks, 5u, &count));
    ASSERT_EQ(4u, count);
    const uint32_t starts[] = {1u, 4u, 7u, 16u}, ref_len[] = {0u, 1u, 4u, 1u};
    const size_t alt_start0[] = {0u, 6u, 9u, 18u}, alt_len[] = {3u, 1u, 4u, 1u};
    for (size_t i = 0u; i < count; i++) {
        ASSERT_EQ(starts[i], blocks[i].cds_start);
        ASSERT_EQ(ref_len[i], blocks[i].ref_len);
        ASSERT_EQ(alt_start0[i], blocks[i].alt_start0);
        ASSERT_EQ(alt_len[i], blocks[i].alt_len);
    }
    ASSERT_EQ(2u, blocks[2].edit_count);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_FLAG_INDEL | DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT,
        blocks[2].flags);
    /* Delete a complete codon, then substitute in the next: the alternate axis
     * moves backwards, but the reference-axis spans must not move. */
    edits[0] = (duckvep_haplotype_edit_t){1u, 3u, (const uint8_t *)"AAA", 0u, NULL, 1};
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, duckvep_haplotype_partition(edits, 2u, blocks, 5u, &count));
    /* Both edits touch the first alternate codon after deleting the first one. */
    ASSERT_EQ(1u, count);
    ASSERT_EQ(4u, blocks[0].ref_len);
    ASSERT_EQ(0u, blocks[0].alt_start0);
    ASSERT_EQ(1u, blocks[0].alt_len);
    PASS();
}

TEST haplotype_ordered_replacements_validate_before_mutating(void) {
    const uint8_t *ref = (const uint8_t *)"ACGTACGTACGT";
    duckvep_haplotype_edit_t edits[] = {
        {4u, 1u, ref + 3u, 1u, ref, 1},
        {1u, 4u, ref, 4u, ref, 1}
    };
    uint8_t cds[32]; uint64_t ids[2] = {11u, 12u};
    duckvep_haplotype_block_t blocks[2];
    size_t count;
    duckvep_haplotype_result_t result;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, duckvep_haplotype_compose_replacements(ref, 12u,
        edits, 2u, 1, ids, cds, 12u, blocks, 2u, &count, &result));
    ASSERT_MEM_EQ(ref, cds, 12u); ASSERT_EQ(2u, result.applied_edits); ASSERT_EQ(1u, count);
    ASSERT_EQ(0, blocks[0].length_diff); ASSERT_EQ(4u, blocks[0].ref_len);
    ASSERT_EQ(11u, ids[0]); ASSERT_EQ(12u, ids[1]);
    for (unsigned scenario = 0u; scenario < 8u; scenario++) {
        duckvep_haplotype_edit_t copy[2]; memcpy(copy, edits, sizeof(copy));
        memset(cds, 0xa5, sizeof(cds)); memset(blocks, 0xa5, sizeof(blocks));
        ids[0] = 11u; ids[1] = 12u;
        size_t cap = sizeof(cds), block_cap = 2u;
        duckvep_haplotype_status_t want = DUCKVEP_HAPLOTYPE_INVALID_ARG;
        if (!scenario) { cap = 11u; want = DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL; }
        if (scenario == 1u) { block_cap = 1u; want = DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL; }
        if (scenario == 2u) { copy[0].ref = ref; want = DUCKVEP_HAPLOTYPE_REF_MISMATCH; }
        if (scenario == 3u) { copy[0].alt = (const uint8_t *)"?"; want = DUCKVEP_HAPLOTYPE_INVALID_BASE; }
        if (scenario == 4u) { copy[1].cds_start = 5u; want = DUCKVEP_HAPLOTYPE_EDIT_ORDER; }
        if (scenario == 5u) { copy[1].ref_len = 13u; want = DUCKVEP_HAPLOTYPE_OUT_OF_RANGE; }
        if (scenario == 6u) copy[0].alt = cds;
        if (scenario == 7u) copy[1].ref_len = 0u;
        ASSERT_EQ(want, duckvep_haplotype_compose_replacements(ref, 12u, copy, 2u,
            1, ids, cds, cap, blocks, block_cap, &count, &result));
        ASSERT_EQ(0u, count); ASSERT_EQ(0u, result.cds_len); ASSERT_EQ(0u, result.applied_edits);
        ASSERT_EQ(11u, ids[0]); ASSERT_EQ(12u, ids[1]);
        for (size_t i = 0u; i < sizeof(cds); i++) ASSERT_EQ(0xa5u, cds[i]);
        for (size_t i = 0u; i < sizeof(blocks); i++) ASSERT_EQ(0xa5u, ((uint8_t *)blocks)[i]);
    }
    PASS();
}

enum { KPROP_REPLACEMENT_REF = 96, KPROP_REPLACEMENTS = 32, KPROP_REPLACEMENT_CAP = 4096 };
struct kprop_replacements {
    uint8_t reference[KPROP_REPLACEMENT_REF];
    uint8_t ref[KPROP_REPLACEMENTS][KPROP_REPLACEMENT_REF];
    uint8_t alt[KPROP_REPLACEMENTS][KPROP_REPLACEMENT_REF];
    duckvep_haplotype_edit_t edits[KPROP_REPLACEMENTS];
    size_t length, count;
    int8_t strand;
};
static struct {
    unsigned forward, reverse, ref_slot, clipped, noop, tied, empty, merged_sources;
} replacement_coverage;

static enum theft_alloc_res kprop_replacements_alloc(struct theft *t, void *env, void **instance) {
    (void)env;
    struct kprop_replacements *c = calloc(1u, sizeof(*c));
    if (!c) return THEFT_ALLOC_ERROR;
    c->length = 1u + kprop_bounded(t, KPROP_REPLACEMENT_REF);
    c->count = 1u + kprop_bounded(t, KPROP_REPLACEMENTS);
    c->strand = kprop_bounded(t, 2u) ? 1 : -1;
    for (size_t i = 0u; i < c->length; i++) c->reference[i] = (uint8_t)"ACGT"[kprop_bounded(t, 4u)];
    for (size_t i = 0u; i < c->count; i++) {
        duckvep_haplotype_edit_t *e = &c->edits[i];
        e->cds_start = 1u + (uint32_t)kprop_bounded(t, c->length);
        e->ref_len = 1u + (uint32_t)kprop_bounded(t, c->length - e->cds_start + 1u);
        e->alt_len = (uint32_t)kprop_bounded(t, 33u);
        e->variant_strand = kprop_bounded(t, 2u) ? 1 : -1;
        int reverse = e->variant_strand != c->strand;
        unsigned mode = (unsigned)kprop_bounded(t, 8u);
        if (!mode) e->alt_len = e->ref_len;
        for (uint32_t j = 0u; j < e->ref_len; j++)
            c->ref[i][reverse ? e->ref_len - 1u - j : j] =
                haplo_test_variant_from_tx_base((char)c->reference[e->cds_start - 1u + j], reverse);
        for (uint32_t j = 0u; j < e->alt_len; j++)
            c->alt[i][reverse ? e->alt_len - 1u - j : j] = haplo_test_variant_from_tx_base(
                !mode ? (char)c->reference[e->cds_start - 1u + j] : "ACGT"[kprop_bounded(t, 4u)], reverse);
        e->ref = c->ref[i]; e->alt = c->alt[i];
        if (i && mode == 1u) *e = c->edits[i - 1u];
    }
    /* Whole-CDS repeated deletion forces clipping and an empty final sequence. */
    if (!kprop_bounded(t, 16u)) {
        c->count = 2u;
        c->edits[0] = c->edits[1] = (duckvep_haplotype_edit_t){
            1u, (uint32_t)c->length, c->reference, 0u, NULL, c->strand};
    }
    for (size_t i = 1u; i < c->count; i++) {
        duckvep_haplotype_edit_t e = c->edits[i]; size_t j = i;
        while (j && c->edits[j - 1u].cds_start < e.cds_start) {
            c->edits[j] = c->edits[j - 1u]; j--;
        }
        c->edits[j] = e;
    }
    *instance = c;
    return THEFT_ALLOC_OK;
}

static enum theft_trial_res prop_replacements_match_literal_and_reconstruct(struct theft *t, void *arg) {
    (void)t;
    const struct kprop_replacements *c = arg;
    uint8_t expected[KPROP_REPLACEMENT_CAP], got[KPROP_REPLACEMENT_CAP], rebuilt[KPROP_REPLACEMENT_CAP];
    duckvep_haplotype_block_t blocks[KPROP_REPLACEMENTS];
    uint64_t ids[KPROP_REPLACEMENTS], changed = 0u, observed = 0u;
    size_t length = c->length, nblocks, nchanged = 0u;
    int64_t nominal = 0;
    uint32_t flags = 0u;
    int shifted = 0;
    memcpy(expected, c->reference, length);
    replacement_coverage.forward += c->strand == 1;
    replacement_coverage.reverse += c->strand == -1;
    for (size_t i = 0u; i < c->count; i++) {
        ids[i] = i;
        const duckvep_haplotype_edit_t *e = &c->edits[i];
        size_t start = e->cds_start - 1u;
        if (start > length) return THEFT_TRIAL_ERROR;
        size_t removed = e->ref_len < length - start ? e->ref_len : length - start;
        uint8_t alt[KPROP_REPLACEMENT_REF];
        int reverse = e->variant_strand != c->strand;
        for (uint32_t j = 0u; j < e->alt_len; j++) alt[j] =
            haplo_test_variant_from_tx_base((char)e->alt[reverse ? e->alt_len - 1u - j : j], reverse);
        replacement_coverage.ref_slot += e->alt_len == e->ref_len &&
            !memcmp(alt, c->reference + start, e->alt_len);
        replacement_coverage.clipped += removed != e->ref_len;
        replacement_coverage.tied += i && start == c->edits[i - 1u].cds_start - 1u;
        int64_t delta = (int64_t)e->alt_len - e->ref_len;
        nominal += delta;
        if (delta) flags |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
        shifted |= delta % 3 != 0;
        if (removed != e->alt_len || memcmp(expected + start, alt, removed)) {
            changed |= UINT64_C(1) << i; nchanged++;
            memmove(expected + start + e->alt_len, expected + start + removed, length - start - removed);
            memcpy(expected + start, alt, e->alt_len);
            length += e->alt_len; length -= removed;
        } else replacement_coverage.noop++;
    }
    replacement_coverage.empty += !length;
    if (shifted) flags |= nominal % 3 ? DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT
                                     : DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT;
    duckvep_haplotype_result_t result;
    if (duckvep_haplotype_compose_replacements(c->reference, c->length, c->edits, c->count,
            c->strand, ids, got, sizeof(got), blocks, KPROP_REPLACEMENTS, &nblocks, &result) !=
        DUCKVEP_HAPLOTYPE_OK) return THEFT_TRIAL_FAIL;
    if (result.cds_len != length || memcmp(got, expected, length) || result.flags != flags ||
        result.length_diff != nominal || result.applied_edits != nchanged) return THEFT_TRIAL_FAIL;
    size_t used = 0u, from = 0u, provenance = 0u;
    for (size_t i = nblocks; i > 0u; i--) {
        const duckvep_haplotype_block_t *b = &blocks[i - 1u];
        size_t start = b->cds_start - 1u;
        if (start < from || start > c->length || b->ref_len > c->length - start ||
            b->alt_start0 > length || b->alt_len > length - b->alt_start0 ||
            b->edit_begin > nchanged || b->edit_count > nchanged - b->edit_begin)
            return THEFT_TRIAL_FAIL;
        memcpy(rebuilt + used, c->reference + from, start - from); used += start - from;
        if (b->alt_start0 != used) return THEFT_TRIAL_FAIL;
        memcpy(rebuilt + used, got + b->alt_start0, b->alt_len); used += b->alt_len;
        from = start + b->ref_len;
        if (b->length_diff != (int64_t)b->alt_len - b->ref_len) return THEFT_TRIAL_FAIL;
        replacement_coverage.merged_sources += b->edit_count > 1u;
        for (size_t j = b->edit_begin; j < b->edit_begin + b->edit_count; j++) {
            if (ids[j] >= c->count || (observed & (UINT64_C(1) << ids[j]))) return THEFT_TRIAL_FAIL;
            observed |= UINT64_C(1) << ids[j]; provenance++;
        }
    }
    memcpy(rebuilt + used, c->reference + from, c->length - from); used += c->length - from;
    if (used != length || memcmp(rebuilt, got, length) || changed != observed || provenance != nchanged)
        return THEFT_TRIAL_FAIL;
    return THEFT_TRIAL_PASS;
}

static void kprop_replacements_free(void *instance, void *env) { (void)env; free(instance); }

TEST haplotype_ordered_replacements_match_literal_for_overlapping_records(void) {
    struct theft_type_info info = {.alloc = kprop_replacements_alloc, .free = kprop_replacements_free};
    struct theft_run_config cfg = {0};
    cfg.name = "ordered source replacements == literal replay, net spans and applied provenance";
    cfg.prop1 = prop_replacements_match_literal_and_reconstruct;
    cfg.type_info[0] = &info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&replacement_coverage, 0, sizeof(replacement_coverage));
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    fprintf(stderr, "[ordered-replacement coverage] forward=%u reverse=%u ref_slot=%u clipped=%u "
        "noop=%u tied=%u empty=%u merged_sources=%u\n", replacement_coverage.forward,
        replacement_coverage.reverse, replacement_coverage.ref_slot, replacement_coverage.clipped,
        replacement_coverage.noop, replacement_coverage.tied, replacement_coverage.empty,
        replacement_coverage.merged_sources);
    PASS();
}

static enum theft_alloc_res kprop_haplo_alloc(struct theft *t, void *env, void **instance) {
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    struct kprop_haplo_case *c = (struct kprop_haplo_case *)calloc(1u, sizeof *c);
    uint32_t starts_tmp[KPROP_HAPLO_MAX_EDITS];
    uint32_t ref_len_tmp[KPROP_HAPLO_MAX_EDITS];
    uint32_t alt_len_tmp[KPROP_HAPLO_MAX_EDITS];
    uint32_t cursor = 1u;
    size_t target;
    size_t n = 0u;
    size_t i;
    (void)env;
    if (c == NULL) return THEFT_ALLOC_ERROR;
    for (i = 0u; i < KPROP_HAPLO_CDS_LEN; i++) c->ref[i] = (uint8_t)bases[kprop_bounded(t, 4u)];
    c->ref[KPROP_HAPLO_CDS_LEN] = (uint8_t)'\0';
    target = (size_t)kprop_bounded(t, (uint64_t)KPROP_HAPLO_MAX_EDITS + 1u);
    c->transcript_strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;

    /* Build a valid non-overlapping edit set in ascending original CDS order,
     * allowing adjacency, insertions at the CDS end, and variable gaps; reverse
     * it below because the kernel contract is descending application order. */
    while (n < target && cursor <= KPROP_HAPLO_CDS_LEN + 1u) {
        uint32_t gap = (uint32_t)kprop_bounded(t, 3u); /* 0..2: adjacent is common. */
        uint32_t start = cursor + gap;
        uint32_t max_ref;
        uint32_t ref_len;
        if (start > KPROP_HAPLO_CDS_LEN + 1u) break;
        max_ref = start <= KPROP_HAPLO_CDS_LEN ? KPROP_HAPLO_CDS_LEN - start + 1u : 0u;
        if (max_ref > 3u) max_ref = 3u;
        ref_len = max_ref == 0u ? 0u : (uint32_t)kprop_bounded(t, (uint64_t)max_ref + 1u);
        starts_tmp[n] = start;
        ref_len_tmp[n] = ref_len;
        alt_len_tmp[n] = (uint32_t)kprop_bounded(t, KPROP_HAPLO_ALLELE_MAX + 1u); /* 0..4 */
        n++;
        cursor = start + (ref_len > 0u ? ref_len : 1u) + (uint32_t)kprop_bounded(t, 2u);
    }
    c->edit_count = n;

    for (i = 0u; i < c->edit_count; i++) {
        size_t src = c->edit_count - 1u - i;
        duckvep_haplotype_edit_t *e = &c->edits[i];
        uint32_t ref_len = ref_len_tmp[src];
        uint32_t alt_len = alt_len_tmp[src];
        uint32_t j;
        int reverse;
        e->cds_start = starts_tmp[src];
        e->ref_len = ref_len;
        e->alt_len = alt_len;
        e->variant_strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;
        reverse = e->variant_strand != c->transcript_strand;
        for (j = 0u; j < ref_len; j++) {
            char tx = (char)c->ref[(size_t)e->cds_start - 1u + (size_t)j];
            c->ref_alleles[i][reverse ? (ref_len - 1u - j) : j] = haplo_test_variant_from_tx_base(tx, reverse);
        }
        for (j = 0u; j < alt_len; j++) {
            char tx = bases[kprop_bounded(t, 4u)];
            c->alt_alleles[i][reverse ? (alt_len - 1u - j) : j] = haplo_test_variant_from_tx_base(tx, reverse);
        }
        e->ref = ref_len ? c->ref_alleles[i] : NULL;
        e->alt = alt_len ? c->alt_alleles[i] : NULL;
    }
    *instance = c;
    return THEFT_ALLOC_OK;
}
void kprop_haplo_free(void *instance, void *env) { (void)env; free(instance); }
static struct theft_type_info kprop_haplo_info = {
    .alloc = kprop_haplo_alloc,
    .free  = kprop_haplo_free,
};

static enum theft_trial_res prop_haplotype_apply_matches_rebuild_oracle(struct theft *t, void *arg1) {
    const struct kprop_haplo_case *c = (const struct kprop_haplo_case *)arg1;
    uint8_t got[KPROP_HAPLO_CAP];
    uint8_t in_place[KPROP_HAPLO_CAP];
    uint8_t want[KPROP_HAPLO_CAP];
    size_t got_len = 0u, in_place_len = 0u, want_len = 0u;
    int64_t want_diff = 0;
    uint32_t want_flags = 0u;
    duckvep_haplotype_result_t r;
    duckvep_haplotype_result_t in_place_result;
    (void)t;
    if (!haplo_oracle_rebuild(c->ref, KPROP_HAPLO_CDS_LEN, c->edits, c->edit_count,
                              c->transcript_strand, want, sizeof want, &want_len,
                              &want_diff, &want_flags)) {
        return THEFT_TRIAL_ERROR;
    }
    if (duckvep_haplotype_apply_cds_edits(c->ref, KPROP_HAPLO_CDS_LEN, c->edits, c->edit_count,
                                          c->transcript_strand, got, sizeof got, &got_len, &r) !=
        DUCKVEP_HAPLOTYPE_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (got_len != want_len || memcmp(got, want, got_len) != 0) return THEFT_TRIAL_FAIL;
    if (r.cds_len != want_len || r.length_diff != want_diff || r.flags != want_flags ||
        r.applied_edits != c->edit_count) {
        return THEFT_TRIAL_FAIL;
    }
    memcpy(in_place, c->ref, KPROP_HAPLO_CDS_LEN + 1u);
    if (duckvep_haplotype_apply_cds_edits(
            in_place, KPROP_HAPLO_CDS_LEN, c->edits, c->edit_count,
            c->transcript_strand, in_place, sizeof in_place, &in_place_len,
            &in_place_result) != DUCKVEP_HAPLOTYPE_INVALID_ARG) {
        return THEFT_TRIAL_FAIL;
    }
    if (in_place_len != 0u ||
        memcmp(in_place, c->ref, KPROP_HAPLO_CDS_LEN + 1u) != 0 ||
        in_place_result.cds_len != 0u ||
        in_place_result.length_diff != 0 ||
        in_place_result.flags != 0u ||
        in_place_result.applied_edits != 0u) {
        return THEFT_TRIAL_FAIL;
    }
    return THEFT_TRIAL_PASS;
}

TEST haplotype_apply_matches_rebuild_oracle_for_any_valid_edit_set(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "multi-edit CDS haplotype apply == left-to-right rebuild oracle";
    cfg.prop1 = prop_haplotype_apply_matches_rebuild_oracle;
    cfg.type_info[0] = &kprop_haplo_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

static enum theft_trial_res prop_haplotype_partition_preserves_interactions(
    struct theft *t, void *arg1) {
    const struct kprop_haplo_case *c = (const struct kprop_haplo_case *)arg1;
    duckvep_haplotype_edit_t ascending[KPROP_HAPLO_MAX_EDITS];
    duckvep_haplotype_block_t blocks[KPROP_HAPLO_MAX_EDITS];
    size_t block_count = 0u;
    size_t covered = 0u;
    size_t bi;
    (void)t;

    for (bi = 0u; bi < c->edit_count; bi++) {
        ascending[bi] = c->edits[c->edit_count - 1u - bi];
    }
    if (duckvep_haplotype_partition(ascending, c->edit_count, blocks,
                                    KPROP_HAPLO_MAX_EDITS, &block_count) !=
        DUCKVEP_HAPLOTYPE_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (block_count > c->edit_count ||
        (c->edit_count == 0u && block_count != 0u)) {
        return THEFT_TRIAL_FAIL;
    }

    for (bi = 0u; bi < block_count; bi++) {
        const duckvep_haplotype_block_t *block = &blocks[bi];
        int64_t difference = 0;
        int64_t shift_before = 0;
        uint32_t expected_flags = 0u;
        uint32_t expected_end = 0u;
        int saw_frameshift = 0;
        size_t j;

        if (block->edit_begin != covered || block->edit_count == 0u ||
            block->edit_begin + block->edit_count > c->edit_count) {
            return THEFT_TRIAL_FAIL;
        }
        for (j = block->edit_begin;
             j < block->edit_begin + block->edit_count; j++) {
            const duckvep_haplotype_edit_t *edit = &ascending[j];
            int64_t edit_difference = (int64_t)edit->alt_len - (int64_t)edit->ref_len;
            uint32_t edit_end = edit->cds_start +
                (edit->ref_len == 0u ? 0u : edit->ref_len - 1u);
            int should_stay_open = 0;

            if (edit_end > expected_end) expected_end = edit_end;
            difference += edit_difference;
            if (edit_difference != 0) {
                expected_flags |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
                if (edit_difference % 3 != 0) saw_frameshift = 1;
            }
            if (j + 1u < c->edit_count) {
                int64_t current_start = (int64_t)edit->cds_start - 1 + shift_before;
                int64_t current_end = current_start +
                    (edit->alt_len == 0u ? 0 : (int64_t)edit->alt_len - 1);
                int64_t next_start = (int64_t)ascending[j + 1u].cds_start - 1 + difference;
                should_stay_open = difference % 3 != 0 ||
                    current_end / 3 == next_start / 3;
            }
            shift_before += edit_difference;

            if (j + 1u < block->edit_begin + block->edit_count) {
                if (!should_stay_open) return THEFT_TRIAL_FAIL;
            } else if (j + 1u < c->edit_count && should_stay_open) {
                return THEFT_TRIAL_FAIL;
            }
        }
        if (saw_frameshift) {
            expected_flags |= difference % 3 == 0
                ? DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT
                : DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT;
        }
        if (block->cds_start != ascending[block->edit_begin].cds_start ||
            block->cds_start - 1u + block->ref_len +
                (ascending[block->edit_begin + block->edit_count - 1u].ref_len == 0u)
                != expected_end ||
            block->length_diff != difference ||
            block->flags != expected_flags) {
            return THEFT_TRIAL_FAIL;
        }
        covered += block->edit_count;
    }
    return covered == c->edit_count ? THEFT_TRIAL_PASS : THEFT_TRIAL_FAIL;
}

TEST haplotype_partition_preserves_interactions_for_any_valid_edit_set(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "haplotype blocks preserve every frame and same-codon interaction";
    cfg.prop1 = prop_haplotype_partition_preserves_interactions;
    cfg.type_info[0] = &kprop_haplo_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* Additional evidence lane: the original generator, partition property and
 * left-to-right oracle stay unchanged. Replacing each block's entire span
 * must reproduce the original physical edits on both transcript strands. */
static enum theft_trial_res prop_haplotype_block_spans_reconstruct(struct theft *t, void *arg) {
    const struct kprop_haplo_case *c = arg;
    duckvep_haplotype_edit_t ascending[KPROP_HAPLO_MAX_EDITS], composite[KPROP_HAPLO_MAX_EDITS];
    duckvep_haplotype_block_t blocks[KPROP_HAPLO_MAX_EDITS];
    uint8_t expected[KPROP_HAPLO_CAP], rebuilt[KPROP_HAPLO_CAP];
    size_t expected_len, rebuilt_len, count;
    int64_t difference;
    uint32_t flags;
    (void)t;
    if (!haplo_oracle_rebuild(c->ref, KPROP_HAPLO_CDS_LEN, c->edits, c->edit_count,
            c->transcript_strand, expected, sizeof(expected), &expected_len, &difference, &flags))
        return THEFT_TRIAL_ERROR;
    for (size_t i = 0u; i < c->edit_count; i++) ascending[i] = c->edits[c->edit_count - 1u - i];
    if (duckvep_haplotype_partition(ascending, c->edit_count, blocks,
            KPROP_HAPLO_MAX_EDITS, &count) != DUCKVEP_HAPLOTYPE_OK)
        return THEFT_TRIAL_FAIL;
    for (size_t i = 0u; i < count; i++) {
        const duckvep_haplotype_block_t *b = &blocks[i];
        size_t start0 = (size_t)b->cds_start - 1u;
        if (start0 > KPROP_HAPLO_CDS_LEN || b->ref_len > KPROP_HAPLO_CDS_LEN - start0 ||
            b->alt_start0 > expected_len || b->alt_len > expected_len - b->alt_start0 ||
            b->alt_len > UINT32_MAX ||
            b->length_diff != (int64_t)b->alt_len - (int64_t)b->ref_len)
            return THEFT_TRIAL_FAIL;
        composite[count - 1u - i] = (duckvep_haplotype_edit_t){b->cds_start, b->ref_len,
            c->ref + start0, (uint32_t)b->alt_len, expected + b->alt_start0, c->transcript_strand};
    }
    if (!haplo_oracle_rebuild(c->ref, KPROP_HAPLO_CDS_LEN, composite, count,
            c->transcript_strand, rebuilt, sizeof(rebuilt), &rebuilt_len, &difference, &flags))
        return THEFT_TRIAL_FAIL;
    return expected_len == rebuilt_len && !memcmp(expected, rebuilt, expected_len)
        ? THEFT_TRIAL_PASS : THEFT_TRIAL_FAIL;
}

TEST haplotype_block_spans_reconstruct_every_generated_edit_set(void) {
    struct theft_run_config cfg = {0};
    cfg.name = "haplotype block spans reconstruct the independently replayed CDS";
    cfg.prop1 = prop_haplotype_block_spans_reconstruct;
    cfg.type_info[0] = &kprop_haplo_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}
