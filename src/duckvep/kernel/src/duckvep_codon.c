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

int duckvep_codon_table_supported(duckvep_codon_table_t table) {
    unsigned int id = (unsigned int)table;
    return id < sizeof AA_TABLES / sizeof AA_TABLES[0] &&
           AA_TABLES[id] != NULL;
}

char duckvep_translate_codon(const char *codon3, duckvep_codon_table_t table) {
    const char *tab;
    int b1, b2, b3, idx;
    if (codon3 == NULL || !duckvep_codon_table_supported(table)) return 'X';
    b1 = base2bit(codon3[0]);
    b2 = base2bit(codon3[1]);
    b3 = base2bit(codon3[2]);
    if (b1 < 0 || b2 < 0 || b3 < 0) return 'X';
    idx = (b1 << 4) | (b2 << 2) | b3;
    tab = AA_TABLES[(unsigned int)table];
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
