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
