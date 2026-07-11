/*
 * duckvep_codon.c — codon translation + coding-change classification.
 * See duckvep_codon.h. No DuckDB/htslib; no allocation.
 *
 * Tables are indexed by codon = (b1<<4)|(b2<<2)|b3 with T=0,C=1,A=2,G=3. The
 * vertebrate mitochondrial table is the standard table with four edits:
 *   TGA  (idx 14): '*' -> 'W'        ATA  (idx 34): 'I' -> 'M'
 *   AGA  (idx 46): 'R' -> '*'        AGG  (idx 47): 'R' -> '*'
 */
#include "duckvep_codon.h"
#include "duckvep_dna.h"

#include <stddef.h>

static const char AA_STANDARD[65] =
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
static const char AA_VERT_MITO[65] =
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";

static int base2bit(char c) {
    return duckvep_dna_codon_code(c);
}

char duckvep_translate_codon(const char *codon3, duckvep_codon_table_t table) {
    const char *tab;
    int b1, b2, b3, idx;
    if (codon3 == NULL) return 'X';
    b1 = base2bit(codon3[0]);
    b2 = base2bit(codon3[1]);
    b3 = base2bit(codon3[2]);
    if (b1 < 0 || b2 < 0 || b3 < 0) return 'X';
    idx = (b1 << 4) | (b2 << 2) | b3;
    tab = (table == DUCKVEP_CODON_TABLE_VERT_MITO) ? AA_VERT_MITO : AA_STANDARD;
    return tab[idx];
}

duckvep_codon_result_t duckvep_codon_change(const char *ref3, const char *alt3,
                                            duckvep_codon_table_t table) {
    duckvep_codon_result_t r;
    r.aa_ref = duckvep_translate_codon(ref3, table);
    r.aa_alt = duckvep_translate_codon(alt3, table);
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
