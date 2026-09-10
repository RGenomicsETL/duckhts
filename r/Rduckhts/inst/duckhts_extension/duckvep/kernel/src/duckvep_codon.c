/*
 * duckvep_codon.c — codon translation + coding-change classification.
 * See duckvep_codon.h. No DuckDB/htslib; no allocation.
 *
 * Tables are indexed by codon = (b1<<4)|(b2<<2)|b3 with T=0,C=1,A=2,G=3.
 * The strings are the BioPerl 1.7.8 tables used by Ensembl VEP 116. Removed
 * NCBI ids remain NULL so an imported model cannot silently use table 1.
 */
#include "duckvep_codon.h"
#include "duckvep_dna.h"

#include <stddef.h>
#include <string.h>

static const char *const AA_TABLES[32] = {
    [1]  = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [2]  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    [3]  = "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [4]  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [5]  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    [6]  = "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [9]  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    [10] = "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [11] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [12] = "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [13] = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    [14] = "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    [16] = "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [21] = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    [22] = "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [23] = "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [24] = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    [25] = "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [26] = "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [27] = "FFLLSSSSYYQQCCWWLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [28] = "FFLLSSSSYYQQCCWWLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [29] = "FFLLSSSSYYYYCC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [30] = "FFLLSSSSYYEECC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    [31] = "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
};

/* BioPerl 1.7.8 STARTS entries marked M, with the same TCAG codon indexing.
 * Start identity is not implied by ordinary translation to methionine. */
static const uint64_t START_CODONS[32] = {
    [1] = UINT64_C(0x0000000800080008),
    [2] = UINT64_C(0x0008000f00000000),
    [3] = UINT64_C(0x0000000c00000000),
    [4] = UINT64_C(0x0008000f0008000c),
    [5] = UINT64_C(0x0008000f00000008),
    [6] = UINT64_C(0x0000000800000000),
    [9] = UINT64_C(0x0008000800000000),
    [10] = UINT64_C(0x0000000800000000),
    [11] = UINT64_C(0x0008000f00080008),
    [12] = UINT64_C(0x0000000800080000),
    [13] = UINT64_C(0x0008000c00000008),
    [14] = UINT64_C(0x0000000800000000),
    [16] = UINT64_C(0x0000000800000000),
    [21] = UINT64_C(0x0008000800000000),
    [22] = UINT64_C(0x0000000800000000),
    [23] = UINT64_C(0x0008000900000000),
    [24] = UINT64_C(0x0008000800080008),
    [25] = UINT64_C(0x0008000800000008),
    [26] = UINT64_C(0x0000000800080000),
    [27] = UINT64_C(0x0000000800000000),
    [28] = UINT64_C(0x0000000800000000),
    [29] = UINT64_C(0x0000000800000000),
    [30] = UINT64_C(0x0000000800000000),
    [31] = UINT64_C(0x0000000800000000)
};

static int codon_base_masks(const uint8_t *codon, uint8_t masks[3]) {
    if (!codon) return 0;
    for (size_t i = 0u; i < 3u; i++) {
        char base = duckvep_dna_normalize((char)codon[i], 1);
        if (!base) return 0;
        masks[i] = base == 'N' ? 15u : (uint8_t)(1u << duckvep_dna_codon_code(base));
    }
    return 1;
}

int duckvep_codon_is_start(const uint8_t *codon3, duckvep_codon_table_t table) {
    uint8_t masks[3];
    if (!duckvep_codon_table_supported(table) || !codon_base_masks(codon3, masks)) return 0;
    uint64_t starts = START_CODONS[(unsigned)table];
    for (unsigned a = 0u; a < 4u; a++) if (masks[0] & (1u << a))
        for (unsigned b = 0u; b < 4u; b++) if (masks[1] & (1u << b))
            for (unsigned c = 0u; c < 4u; c++) if (masks[2] & (1u << c))
                if (starts & (UINT64_C(1) << ((a << 4u) | (b << 2u) | c))) return 1;
    return 0;
}

/* Expand at most 64 codons. The shared literal table remains the only amino
 * acid authority; this is BioPerl 1.7.8 _translate_ambiguous_codon over ACGTN. */
static uint8_t codon_n_consensus(const uint8_t *codon, const char *amino_acids) {
    uint8_t masks[3];
    if (!codon_base_masks(codon, masks)) return 'X';
    uint32_t seen = 0u;
    for (unsigned a = 0u; a < 4u; a++) if (masks[0] & (1u << a))
        for (unsigned b = 0u; b < 4u; b++) if (masks[1] & (1u << b))
            for (unsigned c = 0u; c < 4u; c++) if (masks[2] & (1u << c)) {
                uint8_t aa = (uint8_t)amino_acids[(a << 4u) | (b << 2u) | c];
                seen |= UINT32_C(1) << (aa == '*' ? 26u : (unsigned)(aa - 'A'));
            }
    if (seen == ((UINT32_C(1) << ('D' - 'A')) | (UINT32_C(1) << ('N' - 'A')))) return 'B';
    if (seen == ((UINT32_C(1) << ('E' - 'A')) | (UINT32_C(1) << ('Q' - 'A')))) return 'Z';
    if (seen & (seen - 1u)) return 'X';
    for (unsigned i = 0u; i <= 26u; i++) if (seen == (UINT32_C(1) << i))
        return i == 26u ? (uint8_t)'*' : (uint8_t)('A' + i);
    return 'X';
}

static int translation_overlaps(const void *a, size_t a_size, const void *b, size_t b_size) {
    if (!a || !b || !a_size || !b_size) return 0;
    uintptr_t left = (uintptr_t)a, right = (uintptr_t)b;
    return left <= right ? right - left < a_size : left - right < b_size;
}

static duckvep_translation_status_t translate_cds(
    const uint8_t *cds, size_t cds_length, duckvep_codon_table_t table,
    duckvep_translation_ambiguity_t ambiguity,
    uint8_t *peptide, size_t peptide_capacity, duckvep_translation_t *result,
    uint8_t *conservative, duckvep_translation_t *conservative_result) {
    if (!result) return DUCKVEP_TRANSLATION_INVALID_ARG;
    if ((cds && cds_length > UINTPTR_MAX - (uintptr_t)cds) ||
        (peptide && peptide_capacity > UINTPTR_MAX - (uintptr_t)peptide) ||
        (conservative && peptide_capacity > UINTPTR_MAX - (uintptr_t)conservative))
        return DUCKVEP_TRANSLATION_INVALID_ARG;
    if (translation_overlaps(result, sizeof(*result), cds, cds_length) ||
        translation_overlaps(result, sizeof(*result), peptide, peptide_capacity) ||
        translation_overlaps(result, sizeof(*result), conservative, peptide_capacity) ||
        translation_overlaps(conservative_result, sizeof(*conservative_result), cds, cds_length) ||
        translation_overlaps(conservative_result, sizeof(*conservative_result), peptide, peptide_capacity) ||
        translation_overlaps(conservative_result, sizeof(*conservative_result), conservative, peptide_capacity) ||
        translation_overlaps(result, sizeof(*result), conservative_result, sizeof(*conservative_result)))
        return DUCKVEP_TRANSLATION_INVALID_ARG;
    memset(result, 0, sizeof(*result));
    if (conservative_result) memset(conservative_result, 0, sizeof(*conservative_result));
    if ((conservative != NULL) != (conservative_result != NULL))
        return DUCKVEP_TRANSLATION_INVALID_ARG;
    const char *amino_acids = duckvep_codon_table_amino_acids(table);
    if (!cds || !peptide || !amino_acids || (ambiguity != DUCKVEP_TRANSLATION_N_UNKNOWN &&
        ambiguity != DUCKVEP_TRANSLATION_N_CONSENSUS)) return DUCKVEP_TRANSLATION_INVALID_ARG;
    size_t codons = cds_length / 3u;
    if (peptide_capacity < codons + 1u) return DUCKVEP_TRANSLATION_BUFFER_TOO_SMALL;
    if (translation_overlaps(cds, cds_length, peptide, peptide_capacity) ||
        translation_overlaps(cds, cds_length, conservative, peptide_capacity) ||
        translation_overlaps(peptide, peptide_capacity, conservative, peptide_capacity))
        return DUCKVEP_TRANSLATION_INVALID_ARG;
    duckvep_translation_t translated = {codons, 0u, 1u};
    size_t conservative_stop = 0u;
    static const uint8_t normalized_code[8] = {0u, 2u, 0u, 1u, 0u, 0u, 0u, 3u};
    for (size_t i = 0u; i < codons; i++) {
        uint8_t code = 0u, has_n = 0u;
        for (size_t j = 0u; j < 3u; j++) {
            char base = duckvep_dna_normalize((char)cds[i * 3u + j], 1);
            if (!base) return DUCKVEP_TRANSLATION_INVALID_BASE;
            if (base == 'N') has_n = 1u;
            code = (uint8_t)((code << 2u) | normalized_code[(unsigned char)base & 7u]);
        }
        if (has_n) translated.unambiguous = 0u;
        uint8_t aa = has_n ? (uint8_t)'X' : (uint8_t)amino_acids[code];
        if (conservative) {
            conservative[i] = aa;
            if (aa == '*' && !conservative_stop) conservative_stop = i + 1u;
        }
        if (has_n && ambiguity == DUCKVEP_TRANSLATION_N_CONSENSUS)
            aa = codon_n_consensus(cds + i * 3u, amino_acids);
        peptide[i] = aa;
        if (aa == '*' && !translated.first_stop_position1) translated.first_stop_position1 = i + 1u;
    }
    for (size_t i = codons * 3u; i < cds_length; i++) {
        char base = duckvep_dna_normalize((char)cds[i], 1);
        if (!base) return DUCKVEP_TRANSLATION_INVALID_BASE;
        if (base == 'N') translated.unambiguous = 0u;
    }
    peptide[codons] = 0u;
    *result = translated;
    if (conservative) {
        conservative[codons] = 0u;
        *conservative_result = translated;
        conservative_result->first_stop_position1 = conservative_stop;
    }
    return DUCKVEP_TRANSLATION_OK;
}

duckvep_translation_status_t duckvep_translate_cds(
    const uint8_t *cds, size_t cds_length, duckvep_codon_table_t table,
    duckvep_translation_ambiguity_t ambiguity,
    uint8_t *peptide, size_t peptide_capacity, duckvep_translation_t *result) {
    return translate_cds(cds, cds_length, table, ambiguity,
        peptide, peptide_capacity, result, NULL, NULL);
}

duckvep_translation_status_t duckvep_translate_reference_cds(
    const uint8_t *cds, size_t cds_length, duckvep_codon_table_t table,
    uint8_t *consensus, uint8_t *conservative, size_t capacity,
    duckvep_translation_t *conservative_result) {
    if (!conservative_result) return DUCKVEP_TRANSLATION_INVALID_ARG;
    duckvep_translation_t consensus_result;
    return translate_cds(cds, cds_length, table, DUCKVEP_TRANSLATION_N_CONSENSUS,
        consensus, capacity, &consensus_result, conservative, conservative_result);
}

static int base2bit(char c) {
    return duckvep_dna_codon_code(c);
}

static char translate_codon_with_table(const char *codon3,
                                       const char *amino_acids) {
    int b1, b2, b3, idx;

    if (codon3 == NULL || amino_acids == NULL) return 'X';
    b1 = base2bit(codon3[0]);
    b2 = base2bit(codon3[1]);
    b3 = base2bit(codon3[2]);
    if (b1 < 0 || b2 < 0 || b3 < 0) return 'X';
    idx = (b1 << 4) | (b2 << 2) | b3;
    return amino_acids[idx];
}

static duckvep_codon_result_t codon_change_from_amino_acids(
    char aa_ref,
    char aa_alt) {

    duckvep_codon_result_t r;

    r.aa_ref = aa_ref;
    r.aa_alt = aa_alt;
    if (r.aa_ref == 'X' || r.aa_alt == 'X') {
        r.change = DUCKVEP_CODON_INVALID;
    } else if (r.aa_ref == r.aa_alt) {
        r.change = DUCKVEP_CODON_SYNONYMOUS;
    } else if (r.aa_alt == '*') {
        r.change = DUCKVEP_CODON_STOP_GAINED;
    } else if (r.aa_ref == '*') {
        r.change = DUCKVEP_CODON_STOP_LOST;
    } else {
        r.change = DUCKVEP_CODON_MISSENSE;
    }
    return r;
}

static int prepared_codon_index(const char *codon3, uint8_t *index_out) {
    static const uint8_t code[8] = {0u, 2u, 0u, 1u, 0u, 0u, 0u, 3u};
    unsigned char b0;
    unsigned char b1;
    unsigned char b2;

    if (codon3 == NULL || index_out == NULL) return 0;
    b0 = (unsigned char)codon3[0];
    b1 = (unsigned char)codon3[1];
    b2 = (unsigned char)codon3[2];
    if (b0 == (unsigned char)'N' || b1 == (unsigned char)'N' ||
        b2 == (unsigned char)'N') {
        return 0;
    }
    *index_out = (uint8_t)((code[b0 & 7u] << 4u) |
                           (code[b1 & 7u] << 2u) |
                           code[b2 & 7u]);
    return 1;
}

int duckvep_codon_table_supported(duckvep_codon_table_t table) {
    unsigned int id = (unsigned int)table;
    return id < sizeof AA_TABLES / sizeof AA_TABLES[0] &&
           AA_TABLES[id] != NULL;
}

const char *duckvep_codon_table_amino_acids(duckvep_codon_table_t table) {
    return duckvep_codon_table_supported(table)
        ? AA_TABLES[(unsigned int)table] : NULL;
}

int duckvep_cds_first_stop_position1(
    const uint8_t         *cds,
    size_t                 cds_length,
    duckvep_codon_table_t  table,
    uint32_t              *position1_out) {

    const char *amino_acids;
    size_t codon_count;
    size_t i;

    if (position1_out != NULL) *position1_out = 0u;
    if (position1_out == NULL || (cds_length != 0u && cds == NULL)) return 0;
    amino_acids = duckvep_codon_table_amino_acids(table);
    if (amino_acids == NULL) return 0;
    codon_count = cds_length / 3u;
    for (i = 0u; i < codon_count; i++) {
        uint8_t code = 0u;
        uint32_t j;
        int has_n = 0;

        for (j = 0u; j < 3u; j++) {
            char base = duckvep_dna_normalize(
                (char)cds[i * 3u + (size_t)j], 1);
            int base_code;

            if (base == '\0') return 0;
            if (base == 'N') {
                has_n = 1;
                base_code = 0;
            } else {
                base_code = duckvep_dna_codon_code(base);
                if (base_code < 0) return 0;
            }
            code = (uint8_t)((code << 2u) | (uint8_t)base_code);
        }
        if (!has_n && amino_acids[code] == '*') {
            if (i >= (size_t)UINT32_MAX) return 0;
            *position1_out = (uint32_t)i + 1u;
            return 1;
        }
    }
    for (i = codon_count * 3u; i < cds_length; i++) {
        if (duckvep_dna_normalize((char)cds[i], 1) == '\0') return 0;
    }
    return 1;
}

char duckvep_translate_codon(const char *codon3, duckvep_codon_table_t table) {
    return translate_codon_with_table(
        codon3, duckvep_codon_table_amino_acids(table));
}

duckvep_codon_result_t duckvep_codon_change(const char *ref3, const char *alt3,
                                            duckvep_codon_table_t table) {
    const char *amino_acids = duckvep_codon_table_amino_acids(table);

    return codon_change_from_amino_acids(
        translate_codon_with_table(ref3, amino_acids),
        translate_codon_with_table(alt3, amino_acids));
}

duckvep_codon_result_t duckvep_codon_change_prepared(
    const char *ref3,
    const char *alt3,
    duckvep_codon_table_t table) {

    const char *amino_acids = duckvep_codon_table_amino_acids(table);
    uint8_t ref_index;
    uint8_t alt_index;
    char aa_ref = 'X';
    char aa_alt = 'X';

    if (amino_acids != NULL) {
        if (prepared_codon_index(ref3, &ref_index))
            aa_ref = amino_acids[ref_index];
        if (prepared_codon_index(alt3, &alt_index))
            aa_alt = amino_acids[alt_index];
    }
    return codon_change_from_amino_acids(aa_ref, aa_alt);
}
